#!/usr/bin/env python3
"""
CASSL (Cell-type Annotation using Semi-Supervised Learning) Implementation
Based on the original paper methodology
"""

import numpy as np
import pandas as pd
import math
import random
from sklearn.decomposition import NMF
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, fowlkes_mallows_score
from kneed import KneeLocator
import argparse
import pickle
import sys

class CASSL:
    def __init__(self, verbose=True):
        self.verbose = verbose
        self.optimal_dim = None
        self.nmf_model = None
        self.final_clusters = []
        self.cluster_labels = []
        self.final_centroids = []
        
    def log(self, message):
        if self.verbose:
            print(message)
    
    def find_optimal_dimensions(self, data, max_dim=15):
        """Find optimal NMF dimensions using reconstruction error and knee detection"""
        self.log("Finding optimal NMF dimensions...")
        
        dimensions = range(2, min(max_dim + 1, min(data.shape)))
        errors = []
        
        for k in dimensions:
            self.log(f"Testing k={k}...")
            model = NMF(init="nndsvd", n_components=k, max_iter=100, random_state=42)
            model.fit(data)
            errors.append(model.reconstruction_err_)
        
        # Use knee detection to find optimal dimensions
        kn = KneeLocator(list(dimensions), errors, curve='convex', direction='decreasing')
        self.optimal_dim = kn.knee if kn.knee else dimensions[0]
        
        self.log(f"Optimal dimensions: {self.optimal_dim}")
        return self.optimal_dim
    
    def apply_nmf(self, data):
        """Apply NMF with optimal dimensions"""
        self.log("Applying NMF...")
        self.nmf_model = NMF(n_components=self.optimal_dim, init='nndsvd', max_iter=200, random_state=42)
        W = self.nmf_model.fit_transform(data)
        H = self.nmf_model.components_
        
        self.log(f"NMF shapes - W: {W.shape}, H: {H.shape}")
        return W, H
    
    def distance(self, a, b):
        """Calculate Euclidean distance between two points"""
        return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))
    
    def recursive_kmeans(self, data_points):
        """Recursive k-means clustering as described in CASSL paper"""
        if len(data_points) == 0:
            return
        
        # Find unique known labels in this cluster
        known_labels = set()
        for point in data_points:
            if point[-1] != -1:  # -1 indicates unknown label
                known_labels.add(point[-1])
        
        self.log(f"Cluster has {len(known_labels)} unique labels: {known_labels}")
        
        # Base cases
        if len(known_labels) <= 1:
            if len(known_labels) == 1:
                # All labeled data belong to one class
                self.final_clusters.append(data_points)
                self.cluster_labels.append(list(known_labels)[0])
            else:
                # No labeled data exists
                self.final_clusters.append(data_points)
                self.cluster_labels.append(-1)
            return
        
        # Need to split further - cluster into k=len(known_labels) clusters
        k = len(known_labels)
        self.log(f"Splitting into {k} clusters...")
        
        # Initialize centroids randomly
        random.seed(42)
        centroid_indices = random.sample(range(len(data_points)), min(k, len(data_points)))
        centroids = [data_points[i][:-2] for i in centroid_indices]  # Exclude label and index
        
        # K-means clustering
        max_iter = 10
        for iteration in range(max_iter):
            # Assign points to nearest centroid
            clusters = {i: [] for i in range(k)}
            
            for point in data_points:
                features = point[:-2]  # Exclude label and index
                distances = [self.distance(features, centroid) for centroid in centroids]
                closest_cluster = distances.index(min(distances))
                clusters[closest_cluster].append(point)
            
            # Update centroids
            new_centroids = []
            for i in range(k):
                if clusters[i]:
                    features_matrix = np.array([point[:-2] for point in clusters[i]])
                    new_centroid = np.mean(features_matrix, axis=0)
                    new_centroids.append(new_centroid.tolist())
                else:
                    new_centroids.append(centroids[i])
            
            centroids = new_centroids
        
        # Recursively process each cluster
        for cluster_points in clusters.values():
            if cluster_points:
                self.recursive_kmeans(cluster_points)
    
    def label_unknown_cells(self, data_with_labels):
        """Label unknown cells using nearest cluster approach"""
        self.log("Labeling unknown cells...")
        
        # Calculate final centroids
        self.final_centroids = []
        for cluster in self.final_clusters:
            if cluster:
                features_matrix = np.array([point[:-2] for point in cluster])
                centroid = np.mean(features_matrix, axis=0)
                self.final_centroids.append(centroid)
        
        # Calculate weighted radius threshold
        all_distances = []
        cluster_sizes = []
        
        for i, cluster in enumerate(self.final_clusters):
            distances = []
            for point in cluster:
                dist = self.distance(self.final_centroids[i], point[:-2])
                distances.append(dist)
                all_distances.append(dist)
            cluster_sizes.append(len(cluster))
        
        # Weighted radius
        weighted_radius = sum(max(all_distances[sum(cluster_sizes[:i]):sum(cluster_sizes[:i+1])]) * size 
                             for i, size in enumerate(cluster_sizes)) / sum(cluster_sizes)
        
        # Find unlabeled clusters and assign labels
        unlabeled_indices = [i for i, label in enumerate(self.cluster_labels) if label == -1]
        
        if unlabeled_indices:
            # Calculate distances between centroids
            centroid_distances = np.zeros((len(self.final_centroids), len(self.final_centroids)))
            for i in range(len(self.final_centroids)):
                for j in range(len(self.final_centroids)):
                    centroid_distances[i][j] = self.distance(self.final_centroids[i], self.final_centroids[j])
            
            # Assign labels to unlabeled clusters
            for unlabeled_idx in unlabeled_indices:
                distances = centroid_distances[unlabeled_idx].copy()
                
                # Find nearest labeled cluster
                nearest_labeled_idx = None
                min_distance = float('inf')
                
                for i, label in enumerate(self.cluster_labels):
                    if label != -1 and distances[i] < min_distance:
                        min_distance = distances[i]
                        nearest_labeled_idx = i
                
                if nearest_labeled_idx is not None:
                    # Assign label if within threshold
                    assigned_count = 0
                    for point in self.final_clusters[unlabeled_idx]:
                        point_distance = self.distance(self.final_centroids[nearest_labeled_idx], point[:-2])
                        if point_distance <= weighted_radius:
                            point[-1] = self.cluster_labels[nearest_labeled_idx]
                            assigned_count += 1
                    
                    # Update cluster label if all points were assigned
                    if assigned_count == len(self.final_clusters[unlabeled_idx]):
                        self.cluster_labels[unlabeled_idx] = self.cluster_labels[nearest_labeled_idx]
                        self.log(f"Cluster {unlabeled_idx} assigned label {self.cluster_labels[nearest_labeled_idx]}")
    
    def fit_predict(self, ref_data, ref_labels, test_data, p_missing=0.9):
        """Main CASSL pipeline"""
        self.log("Starting CASSL analysis...")
        
        # Step 1: Find optimal dimensions
        if self.optimal_dim is None:
            self.find_optimal_dimensions(ref_data)
        
        # Step 2: Apply NMF to reference data
        ref_W, ref_H = self.apply_nmf(ref_data)
        
        # Step 3: Project test data
        self.log("Projecting test data...")
        test_W = self.nmf_model.transform(test_data)
        
        # Step 4: Combine reference and test data for semi-supervised learning
        combined_data = np.vstack([ref_W, test_W])
        
        # Step 5: Prepare labels (encode and simulate missing labels)
        le = LabelEncoder()
        encoded_ref_labels = le.fit_transform(ref_labels)
        
        # Create combined labels (reference + unknown for test)
        combined_labels = np.concatenate([encoded_ref_labels, [-1] * len(test_data)])
        
        # Simulate missing labels in reference data
        n_missing = int(p_missing * len(ref_data))
        if n_missing > 0:
            missing_indices = random.sample(range(len(ref_data)), n_missing)
            for idx in missing_indices:
                combined_labels[idx] = -1
        
        # Step 6: Prepare data for recursive clustering
        # Add index column
        indices = np.arange(len(combined_data))
        data_with_labels = []
        
        for i in range(len(combined_data)):
            point = list(combined_data[i]) + [indices[i], combined_labels[i]]
            data_with_labels.append(point)
        
        # Step 7: Apply recursive k-means
        self.log("Applying recursive k-means clustering...")
        self.final_clusters = []
        self.cluster_labels = []
        
        # Determine number of unique classes
        unique_labels = len(set(encoded_ref_labels))
        
        # Initialize centroids for first clustering
        random.seed(42)
        initial_centroids = random.sample(range(len(data_with_labels)), unique_labels)
        centroids = [data_with_labels[i][:-2] for i in initial_centroids]
        
        # Initial k-means
        clusters = {i: [] for i in range(unique_labels)}
        for point in data_with_labels:
            features = point[:-2]
            distances = [self.distance(features, centroid) for centroid in centroids]
            closest_cluster = distances.index(min(distances))
            clusters[closest_cluster].append(point)
        
        # Apply recursive clustering to each initial cluster
        for cluster_points in clusters.values():
            if cluster_points:
                self.recursive_kmeans(cluster_points)
        
        # Step 8: Label unknown cells
        self.label_unknown_cells(data_with_labels)
        
        # Step 9: Extract predictions for test data
        test_predictions = []
        test_start_idx = len(ref_data)
        
        for cluster in self.final_clusters:
            for point in cluster:
                original_idx = int(point[-2])
                if original_idx >= test_start_idx:
                    test_idx = original_idx - test_start_idx
                    predicted_label = point[-1]
                    if predicted_label != -1:
                        # Decode label back to original
                        decoded_label = le.classes_[predicted_label] if predicted_label < len(le.classes_) else predicted_label
                        test_predictions.append((test_idx, decoded_label))
        
        # Sort by test index
        test_predictions.sort(key=lambda x: x[0])
        
        # Create final prediction array
        final_predictions = ['Unknown'] * len(test_data)
        for test_idx, pred_label in test_predictions:
            if test_idx < len(final_predictions):
                final_predictions[test_idx] = pred_label
        
        self.log(f"CASSL completed. Predicted {len(test_predictions)} out of {len(test_data)} test cells")
        
        return {
            'ref_embedding': ref_W,
            'test_embedding': test_W,
            'test_predictions': final_predictions,
            'optimal_dim': self.optimal_dim,
            'n_clusters': len(self.final_clusters)
        }

def main():
    parser = argparse.ArgumentParser(description='CASSL - Cell-type Annotation using Semi-Supervised Learning')
    parser.add_argument('--ref_data', required=True, help='Reference data file (CSV)')
    parser.add_argument('--ref_labels', required=True, help='Reference labels file (CSV)')
    parser.add_argument('--test_data', required=True, help='Test data file (CSV)')
    parser.add_argument('--output', required=True, help='Output file (pickle)')
    parser.add_argument('--p_missing', type=float, default=0.9, help='Proportion of missing labels')
    parser.add_argument('--max_dim', type=int, default=15, help='Maximum dimensions for NMF')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    try:
        # Load data
        ref_data = pd.read_csv(args.ref_data, index_col=0).values
        ref_labels = pd.read_csv(args.ref_labels, index_col=0).iloc[:, 0].values
        test_data = pd.read_csv(args.test_data, index_col=0).values
        
        # Run CASSL
        cassl = CASSL(verbose=args.verbose)
        results = cassl.fit_predict(ref_data, ref_labels, test_data, args.p_missing)
        
        # Save results
        with open(args.output, 'wb') as f:
            pickle.dump(results, f)
        
        print("CASSL analysis completed successfully!")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
