#!/usr/bin/env python3
"""
SIDA-based Method for Single-cell RNA-seq Integration and Cell Type Mapping
Adapted from the pancreas_5dataset_all_cls.py script for automated benchmarking
"""

import os
import sys
import argparse
import pickle
import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import normalize
from sklearn.metrics import accuracy_score
import random
from collections import Counter
from tensorflow.keras import utils
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Input, Lambda, Convolution1D, MaxPooling1D
from tensorflow.keras.layers import Activation, Dropout, Flatten, Dense
from tensorflow.keras.optimizers import Adam
from tensorflow.keras import backend as K
import warnings
warnings.filterwarnings('ignore')

# Set random seeds for reproducibility
np.random.seed(42)
tf.random.set_seed(42)
random.seed(42)

def normalize_and_pca(data, n_components=50):
    """Normalize data and apply PCA"""
    # Library size normalization and log transformation
    data_normalized = data / (data.sum(axis=1, keepdims=True) + 1e-6) * 10000
    data_log = np.log10(data_normalized + 1)
    
    # PCA
    pca = TruncatedSVD(n_components=n_components, random_state=42)
    data_pca = pca.fit_transform(data_log)
    
    return data_pca

def create_siamese_model(input_dim=50, nb_classes=None):
    """Create the siamese network model"""
    
    def euclidean_distance(vects):
        eps = 1e-5
        x, y = vects
        return K.sqrt(K.maximum(K.sum(K.square(x - y), axis=1, keepdims=True), eps))

    def eucl_dist_output_shape(shapes):
        shape1, shape2 = shapes
        return (shape1[0], 1)

    def contrastive_loss(y_true, y_pred):
        margin = 1
        return K.mean(y_true * K.square(y_pred) + (1 - y_true) * K.square(K.maximum(margin - y_pred, 0)))

    # Create base model
    input_shape = (input_dim, 1)
    model = Sequential()
    
    # Convolution layers
    model.add(Convolution1D(64, 5, padding='valid', input_shape=input_shape))
    model.add(Activation('relu'))
    model.add(Convolution1D(64, 5, padding='valid'))
    model.add(Activation('relu'))
    model.add(Convolution1D(64, 3, padding='valid'))
    model.add(Activation('relu'))
    model.add(Convolution1D(64, 3, padding='valid'))
    model.add(Activation('relu'))
    model.add(Dropout(0.2))
    
    model.add(Convolution1D(32, 3, padding='valid'))
    model.add(Activation('relu'))
    model.add(Convolution1D(32, 3, padding='valid'))
    model.add(Activation('relu'))
    model.add(Convolution1D(32, 3, padding='valid'))
    model.add(Activation('relu'))
    model.add(Convolution1D(32, 1, padding='valid'))
    model.add(Activation('relu'))
    model.add(Dropout(0.25))
    
    model.add(Convolution1D(32, 1, padding='valid'))
    model.add(Activation('relu'))
    model.add(Dense(256))
    model.add(Activation('relu'))
    model.add(Convolution1D(32, 1, padding='valid'))
    model.add(Activation('relu'))
    model.add(Dense(128))
    model.add(Activation('relu'))
    model.add(Dropout(0.25))
    
    model.add(Flatten())
    model.add(Dense(64))
    model.add(Activation('relu'))

    # Create siamese network
    input_a = Input(shape=input_shape)
    input_b = Input(shape=input_shape)
    
    processed_a = model(input_a)
    processed_b = model(input_b)
    
    # Classification head
    out1 = Dropout(0.5)(processed_a)
    out1 = Dense(64, activation='relu')(out1)
    out1 = Dense(32, activation='relu')(out1)
    out1 = Dense(nb_classes)(out1)
    out1 = Activation('softmax', name='classification')(out1)
    
    # Distance calculation
    distance = Lambda(euclidean_distance, output_shape=eucl_dist_output_shape, name='CSA')(
        [processed_a, processed_b])
    
    full_model = Model(inputs=[input_a, input_b], outputs=[out1, distance])
    
    # Compile model
    alpha = 0.6
    full_model.compile(
        loss={'classification': 'categorical_crossentropy', 'CSA': contrastive_loss},
        optimizer=Adam(0.00001),
        loss_weights={'classification': 1 - alpha, 'CSA': alpha}
    )
    
    return full_model, model, contrastive_loss

def encode_labels(labels, label_set):
    """Encode string labels to numeric"""
    label_to_code = {label: i for i, label in enumerate(label_set)}
    return [label_to_code[label] for label in labels]

def create_training_pairs(X1, y1, X2, y2, max_pairs_per_type=2000):
    """Create training pairs for siamese network"""
    positive_pairs = []
    negative_pairs = []
    
    # Create positive pairs (same cell type)
    for i, label1 in enumerate(y1):
        for j, label2 in enumerate(y2):
            if label1 == label2:
                positive_pairs.append([i, j])
            else:
                negative_pairs.append([i, j])
    
    # Balance pairs
    random.shuffle(negative_pairs)
    max_neg = min(len(negative_pairs), 2 * len(positive_pairs), max_pairs_per_type)
    training_pairs = positive_pairs + negative_pairs[:max_neg]
    random.shuffle(training_pairs)
    
    # Create training arrays
    X_a = np.zeros([len(training_pairs), X1.shape[1]])
    X_b = np.zeros([len(training_pairs), X2.shape[1]])
    y_a = np.zeros([len(training_pairs)])
    y_b = np.zeros([len(training_pairs)])
    y_contrastive = np.zeros([len(training_pairs)])
    
    for i, (idx1, idx2) in enumerate(training_pairs):
        X_a[i, :] = X1[idx1, :]
        X_b[i, :] = X2[idx2, :]
        y_a[i] = y1[idx1]
        y_b[i] = y2[idx2]
        y_contrastive[i] = 1 if y1[idx1] == y2[idx2] else 0
    
    return X_a, X_b, y_a, y_b, y_contrastive

def train_siamese_model(model, training_data_list, test_data_list, epochs=10, batch_size=512):
    """Train the siamese model"""
    best_acc = 0
    
    print(f'Training the model - Epochs: {epochs}')
    
    for epoch in range(epochs):
        print(f'Epoch {epoch + 1}/{epochs}')
        
        # Train on all batch pairs
        for train_data in training_data_list:
            X_a, X_b, y_a, y_b, y_contrastive = train_data
            
            # Reshape for CNN
            X_a_reshaped = X_a.reshape(X_a.shape[0], X_a.shape[1], 1)
            X_b_reshaped = X_b.reshape(X_b.shape[0], X_b.shape[1], 1)
            
            # Train in batches
            n_batches = len(y_a) // batch_size
            for i in range(n_batches):
                start_idx = i * batch_size
                end_idx = (i + 1) * batch_size
                
                model.train_on_batch(
                    [X_a_reshaped[start_idx:end_idx], X_b_reshaped[start_idx:end_idx]],
                    [y_a[start_idx:end_idx], y_contrastive[start_idx:end_idx]]
                )
        
        # Evaluate on test data
        total_acc = 0
        for i, (X_test, y_test) in enumerate(test_data_list):
            X_test_reshaped = X_test.reshape(X_test.shape[0], X_test.shape[1], 1)
            out = model.predict([X_test_reshaped, X_test_reshaped], verbose=0)
            pred = np.argmax(out[0], axis=1)
            true = np.argmax(y_test, axis=1)
            acc = np.mean(pred == true)
            total_acc += acc
            
        avg_acc = total_acc / len(test_data_list)
        print(f'Average accuracy: {avg_acc:.4f}')
        
        if avg_acc > best_acc:
            best_acc = avg_acc
    
    return best_acc

def get_embeddings(full_model, base_model, data):
    """Get embeddings from the trained model"""
    # Extract features using the base model (before the siamese network)
    data_reshaped = data.reshape(data.shape[0], data.shape[1], 1)
    embeddings = base_model.predict(data_reshaped, verbose=0)
    
    return embeddings

def main():
    parser = argparse.ArgumentParser(description='SIDA method for cell type mapping')
    parser.add_argument('--ref_data', required=True, help='Path to reference data CSV')
    parser.add_argument('--ref_labels', required=True, help='Path to reference labels CSV')
    parser.add_argument('--test_data', required=True, help='Path to test data CSV')
    parser.add_argument('--test_labels', required=False, help='Path to test labels CSV (for evaluation)')
    parser.add_argument('--output', required=True, help='Path to output pickle file')
    parser.add_argument('--epochs', type=int, default=10, help='Number of training epochs')
    parser.add_argument('--pca_dims', type=int, default=50, help='Number of PCA dimensions')
    
    args = parser.parse_args()
    
    try:
        # Load data
        print("Loading data...")
        ref_data = pd.read_csv(args.ref_data, index_col=0)
        ref_labels_df = pd.read_csv(args.ref_labels, index_col=0)
        test_data = pd.read_csv(args.test_data, index_col=0)
        
        ref_labels = ref_labels_df['labels'].tolist()
        
        # Load test labels if provided
        test_labels = None
        if args.test_labels:
            test_labels_df = pd.read_csv(args.test_labels, index_col=0)
            test_labels = test_labels_df['labels'].tolist()
        
        print(f"Reference data shape: {ref_data.shape}")
        print(f"Test data shape: {test_data.shape}")
        print(f"Reference cell types: {Counter(ref_labels)}")
        if test_labels:
            print(f"Test cell types: {Counter(test_labels)}")
        
        # Normalize and apply PCA
        print("Preprocessing data...")
        ref_data_pca = normalize_and_pca(ref_data.values, n_components=args.pca_dims)
        test_data_pca = normalize_and_pca(test_data.values, n_components=args.pca_dims)
        
        # Get unique cell types and encode labels
        all_cell_types = sorted(list(set(ref_labels)))
        nb_classes = len(all_cell_types)
        
        ref_labels_encoded = encode_labels(ref_labels, all_cell_types)
        if test_labels:
            test_labels_encoded = encode_labels(test_labels, all_cell_types)
        
        print(f"Number of cell types: {nb_classes}")
        
        # Create and train model
        print("Creating model...")
        model, base_model, contrastive_loss = create_siamese_model(
            input_dim=args.pca_dims, 
            nb_classes=nb_classes
        )
        
        # For simplicity, we'll treat ref data as batch 1 and test data as batch 2
        # In practice, you might split reference data into multiple batches
        print("Creating training pairs...")
        X_a, X_b, y_a, y_b, y_contrastive = create_training_pairs(
            ref_data_pca, ref_labels_encoded,
            ref_data_pca, ref_labels_encoded  # Self-pairs for now
        )
        
        # Convert labels to categorical
        y_a_cat = utils.to_categorical(y_a, nb_classes)
        y_b_cat = utils.to_categorical(y_b, nb_classes)
        
        # Prepare test data
        ref_test = ref_data_pca.reshape(ref_data_pca.shape[0], ref_data_pca.shape[1], 1)
        ref_labels_cat = utils.to_categorical(ref_labels_encoded, nb_classes)
        
        test_test = test_data_pca.reshape(test_data_pca.shape[0], test_data_pca.shape[1], 1)
        if test_labels:
            test_labels_cat = utils.to_categorical(test_labels_encoded, nb_classes)
            test_data_list = [(ref_test, ref_labels_cat), (test_test, test_labels_cat)]
        else:
            test_data_list = [(ref_test, ref_labels_cat)]
        
        # Train model
        print("Training model...")
        training_data_list = [(X_a, X_b, y_a_cat, y_b_cat, y_contrastive)]
        best_acc = train_siamese_model(
            model, training_data_list, test_data_list, 
            epochs=args.epochs
        )
        
        # Get embeddings
        print("Extracting embeddings...")
        ref_embeddings = get_embeddings(model, base_model, ref_data_pca)
        test_embeddings = get_embeddings(model, base_model, test_data_pca)
        
        # Predict test labels using the trained model
        test_predictions = None
        if test_labels:
            test_reshaped = test_data_pca.reshape(test_data_pca.shape[0], test_data_pca.shape[1], 1)
            pred_out = model.predict([test_reshaped, test_reshaped], verbose=0)
            test_predictions = np.argmax(pred_out[0], axis=1)
            
            # Calculate accuracy
            accuracy = np.mean(test_predictions == test_labels_encoded)
            print(f"Final test accuracy: {accuracy:.4f}")
        
        # Prepare results
        results = {
            'ref_embedding': ref_embeddings,
            'test_embedding': test_embeddings,
            'optimal_dim': ref_embeddings.shape[1],
            'cell_types': all_cell_types,
            'ref_labels': ref_labels_encoded,
            'test_predictions': test_predictions,
            'model_performance': best_acc
        }
        
        if test_labels:
            results['test_labels'] = test_labels_encoded
            results['accuracy'] = accuracy
        
        # Save results
        print(f"Saving results to {args.output}")
        with open(args.output, 'wb') as f:
            pickle.dump(results, f)
        
        print("SIDA method completed successfully!")
        
    except Exception as e:
        print(f"Error in SIDA method: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
