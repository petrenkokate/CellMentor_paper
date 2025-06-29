# Load all required libraries
print('LIBRARIES')
library(splatter)
library("scater")
library("VariantAnnotation")
library("ggplot2")
library(Seurat)
library(cowplot)
library(dplyr)
library(CellMentor)
library(qs)
library(dittoSeq)
library(profvis)
library(devtools)
library(NNLM)
library(pCMF)
library(glmpca)
library(mclust)
library(cluster)
library(aricode)
library(clevr)
library(SingleR)
library(reticulate)
library(sceasy)
library(harmony)
library(spam)
library(spam64)

SUFFIX <- '_batch_final_newNMF'

simulation_params <- list(
  # Simulation 1
  list(
    nGenes = 15000, 
    seed = 100, 
    batchCells = c(1000, 1000, 1000, 1000), 
    batch.facLoc = 0.05,                    
    batch.facScale = 0.1,
    out.prob = 0.01,                    
    group.prob = c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6), 
    de.prob = c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6), 
    de.downProb = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
    de.facLoc = c(1.2, 1.2, 1.2, 1.2, 1.2, 1.2), 
    de.facScale = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
    bcv.common = 0.2,                     
    bcv.df = 15,                           
    dropout.mid = c(0.5, 0.5, 0.5, 0.5),  
    dropout.shape = c(-1.0, -1.0, -1.0, -1.0),
    dropout.type = 'batch',
    lib.loc = 15, 
    lib.scale = 0.5                        
  ),
  # Simulation 2
  list(
    nGenes = 15000, 
    seed = 100, 
    batchCells = c(1000, 1000, 1000, 1000),
    batch.facLoc = 0.05,
    batch.facScale = 0.1,
    out.prob = 0.01,
    group.prob = c(0.30, 0.25, 0.20, 0.15, 0.05, 0.05), 
    de.prob = c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6),          
    de.downProb = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2),
    de.facLoc = c(1.2, 1.2, 1.2, 1.2, 1.2, 1.2),
    de.facScale = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
    bcv.common = 0.2,
    bcv.df = 15,
    dropout.mid = c(0.5, 0.5, 0.5, 0.5),
    dropout.shape = c(-1.0, -1.0, -1.0, -1.0),
    dropout.type = 'batch',
    lib.loc = 15, 
    lib.scale = 0.5
  ),
  # Simulation 3
  list(
    nGenes = 15000, 
    seed = 100, 
    batchCells = c(1000, 1000, 1000, 1000),
    batch.facLoc = 0.15,
    batch.facScale = 0.2,
    out.prob = 0.02,
    group.prob = c(0.30, 0.25, 0.20, 0.15, 0.05, 0.05),
    de.prob = c(0.7, 0.6, 0.5, 0.4, 0.3, 0.2),
    de.downProb = c(0.4, 0.3, 0.2, 0.3, 0.4, 0.5),
    de.facLoc = c(1.4, 1.2, 1.0, 0.8, 0.6, 0.4),
    de.facScale = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
    bcv.common = 0.3,
    bcv.df = 12,
    dropout.mid = c(0.3, 0.6, 0.9, 1.2),
    dropout.shape = c(-0.8, -1.0, -1.2, -1.4),
    dropout.type = 'batch',
    lib.loc = 15, 
    lib.scale = 1.5                     
  ),
  # Simulation 4
  list(
    nGenes = 15000, 
    seed = 100, 
    batchCells = c(1000, 1000, 1000, 1000),
    batch.facLoc = 0.2,
    batch.facScale = 0.3,
    out.prob = 0.05,                        
    group.prob = c(0.35, 0.30, 0.20, 0.10, 0.03, 0.02), 
    de.prob = c(0.4, 0.3, 0.5, 0.3, 0.2, 0.1),        
    de.downProb = c(0.4, 0.2, 0.5, 0.3, 0.2, 0.6),     
    de.facLoc = c(0.8, 0.6, 1.0, 0.7, 0.5, 0.4),        
    de.facScale = c(0.3, 0.2, 0.3, 0.2, 0.4, 0.3),    
    bcv.common = 0.5,
    bcv.df = 6,                                         
    dropout.mid = c(0.5, 0.8, 1.1, 1.4),                
    dropout.shape = c(-0.7, -0.9, -1.1, -1.3),
    dropout.type = 'batch',
    lib.loc = 14,                                     
    lib.scale = 1.8                                  
  ),
  # Simulation 5
  list(
    nGenes = 15000, 
    seed = 100, 
    batchCells = c(1000, 1000, 1000, 1000),
    batch.facLoc = 0.2,
    batch.facScale = 0.3,
    out.prob = 0.01,
    group.prob = c(0.20, 0.20, 0.20, 0.20, 0.10, 0.10),
    de.prob = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    de.downProb = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2),
    de.facLoc = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
    de.facScale = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
    bcv.common = 0.5,                     
    bcv.df = 8,                             
    dropout.mid = c(0.5, 0.5, 0.5, 0.5),
    dropout.shape = c(-1.0, -1.0, -1.0, -1.0),
    dropout.type = 'batch',
    lib.loc = 15, 
    lib.scale = 0.8
  ),
  # Simulation 6
  list(
    nGenes = 15000, 
    seed = 100, 
    batchCells = c(1000, 1000, 1000, 1000),
    batch.facLoc = 0.2,
    batch.facScale = 0.3,
    out.prob = 0.01,
    group.prob = c(0.20, 0.20, 0.20, 0.20, 0.10, 0.10),
    de.prob = c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4),      
    de.downProb = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2),
    de.facLoc = c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7),  
    de.facScale = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
    bcv.common = 0.5,
    bcv.df = 8,
    dropout.mid = c(0.5, 0.5, 0.5, 0.5),
    dropout.shape = c(-1.0, -1.0, -1.0, -1.0),
    dropout.type = 'batch',
    lib.loc = 15, 
    lib.scale = 0.8
  ),
  # Simulation 7
  list(
    nGenes = 15000, 
    seed = 100, 
    batchCells = c(1000, 1000, 1000, 1000),
    batch.facLoc = 0.2,
    batch.facScale = 0.4,
    out.prob = 0.1, 
    group.prob = c(0.18, 0.23, 0.12,0.07, 0.10, 0.3),
    de.prob = c(0.35, 0.15, 0.45, 0.55, 0.25, 0.15),
    de.downProb = c(0.1, 0.4, 0.2, 0.2, 0.1, 0.6),
    de.facLoc = c(0.6, 0.01, 0.5, 0.45, 0.4, 0.01),
    de.facScale = c(0.1, 0.4, 0.2, 0.25, 0.2, 0.5),
    bcv.common = 0.8,
    bcv.df = 5,
    dropout.mid = c(1, 1, 3, 0.5),
    dropout.shape = c(-1.5, -1.5, -1, -2),
    dropout.type = 'batch',
    lib.loc = 15, 
    lib.scale = 1.0
  ),
  # Simulation 8
  list(
    nGenes = 15000, 
    seed = 100, 
    batchCells = c(1000, 1000, 1000, 1000),
    batch.facLoc = 0.4,
    batch.facScale = 0.5,
    out.prob = 0.07,                           
    group.prob = c(0.35, 0.25, 0.15, 0.10, 0.10, 0.05),
    de.prob = c(0.3, 0.2, 0.4, 0.2, 0.3, 0.2), 
    de.downProb = c(0.3, 0.2, 0.4, 0.1, 0.3, 0.5),
    de.facLoc = c(0.5, 0.4, 0.6, 0.3, 0.5, 0.4),  
    de.facScale = c(0.3, 0.4, 0.3, 0.4, 0.3, 0.4), 
    bcv.common = 0.7,
    bcv.df = 6,
    dropout.mid = c(1.0, 1.5, 0.8, 1.2),
    dropout.shape = c(-1.2, -1.3, -1.1, -1.4),
    dropout.type = 'batch',
    lib.loc = 14,                              
    lib.scale = 1.2                            
  ),
  # Simulation 9
  list(
    nGenes = 10000, 
    seed = 100, 
    batchCells = c(1000, 1000, 1000, 1000),
    batch.facLoc = 0.2,
    batch.facScale = 0.4,
    out.prob = 0.1, 
    group.prob = c(0.18, 0.23, 0.12,0.07, 0.10, 0.3),
    de.prob = c(0.3, 0.1, 0.4, 0.5, 0.2, 0.1),
    de.downProb = c(0.1, 0.4, 0.2, 0.2, 0.1, 0.6),
    de.facLoc = c(0.6, 0.01, 0.5, 0.45, 0.4, 0.01),
    de.facScale = c(0.1, 0.4, 0.2, 0.25, 0.2, 0.5),
    bcv.common = 0.8,
    bcv.df = 5,
    dropout.mid = c(1, 1, 3, 0.5),
    dropout.shape = c(-1.5, -1.5, -1, -2),
    dropout.type = 'batch',
    lib.loc = 15, 
    lib.scale = 1.0
  ),
  # Simulation 10
  list(
    nGenes = 15000, 
    seed = 100, 
    batchCells = c(1000, 1000, 1000, 1000),
    batch.facLoc = 0.2,
    batch.facScale = 0.4,
    out.prob = 0.1, 
    group.prob = c(0.18, 0.23, 0.12,0.07, 0.10, 0.3),
    de.prob = c(0.3, 0.1, 0.4, 0.5, 0.2, 0.1),
    de.downProb = c(0.1, 0.4, 0.2, 0.2, 0.1, 0.6),
    de.facLoc = c(0.6, 0.01, 0.5, 0.45, 0.4, 0.01),
    de.facScale = c(0.1, 0.4, 0.2, 0.25, 0.2, 0.5),
    bcv.common = 0.8,
    bcv.df = 5,
    dropout.mid = c(1, 1, 3, 0.5),
    dropout.shape = c(-1.5, -1.5, -1, -2),
    dropout.type = 'batch',
    lib.loc = 15, 
    lib.scale = 1.0
  )
)

run_simulation <- function(sim_index, params) {
  SIM_NUMBER = paste0(SUFFIX, sim_index)
  set.seed(100)
  
  print(paste("START SIMULATION", sim_index))
  
  # Set simulation parameters
  base_params <- newSplatParams()
  sim_params <- params[[min(sim_index, length(params))]]  # Reuse last param set if index exceeds available sets
  
  # Apply parameters
  params <- setParams(base_params,
                      nGenes = sim_params$nGenes,
                      seed = 100,
                      batchCells = sim_params$batchCells,
                      batch.facLoc = sim_params$batch.facLoc,
                      batch.facScale = sim_params$batch.facScale,
                      out.prob = sim_params$out.prob,
                      group.prob = sim_params$group.prob,
                      de.prob = sim_params$de.prob,
                      de.downProb = sim_params$de.downProb,
                      de.facLoc = sim_params$de.facLoc,
                      de.facScale = sim_params$de.facScale,
                      bcv.common = sim_params$bcv.common,
                      bcv.df = sim_params$bcv.df,
                      dropout.mid = sim_params$dropout.mid,
                      dropout.shape = sim_params$dropout.shape,
                      dropout.type = 'batch',
                      lib.loc = sim_params$lib.loc,
                      lib.scale = sim_params$lib.scale
  )
  
  # Run simulation
  sim <- splatSimulate(params, method = 'groups')
  sim <- logNormCounts(sim)
  
  seu <- as.Seurat(sim)
  
  seu <- seu %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA(verbose=FALSE) %>% 
    RunUMAP(dims = 1:10)
  
  seu$celltype <- seu$Group
  print("SIMULATION DONE!")
  
  # Create reference and test datasets
  print("Create ref")
  seu_ref <- subset(seu, Batch %in% c('Batch1', 'Batch2'))
  
  seu_ref <- seu_ref %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA(verbose=FALSE) %>% 
    RunUMAP(dims = 1:10) %>% 
    FindNeighbors(dims = 1:10) %>% 
    FindClusters(resolution = 0.2)
  
  seu_ref$Seurat_clusters <- seu_ref$originalexp_snn_res.0.2
  
  print("Create test")
  seu_test <- subset(seu, Batch %in% c('Batch3', 'Batch4'))
  
  seu_test <- seu_test %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(verbose=FALSE) %>%
    RunUMAP(dims = 1:10) %>%
    FindNeighbors(dims = 1:10) %>%
    FindClusters(resolution = 0.2)
  
  seu_test$Seurat_clusters <- seu_test$originalexp_snn_res.0.2
  
  print('rssNMF')
  source("code/graph_rssNMF.R")
  source("code/initial_NMF.R")
  
  find_marker_genes <- function(seurat_obj, n_markers = 20) {
    Idents(seurat_obj) <- seurat_obj$Group
    markers <- FindAllMarkers(seurat_obj, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25,
                              verbose = FALSE)
    
    # Get top markers for each cluster
    top_markers <- markers %>%
      group_by(cluster) %>%
      top_n(n = n_markers, wt = avg_log2FC) %>%
      pull(gene) %>%
      unique()
    
    return(top_markers)
  }
  
  # Function to create weight matrix Q based on marker genes
  create_weight_matrix <- function(expr_matrix, marker_genes, sigma = 1) {
    # Filter expression matrix to only marker genes
    marker_expr <- expr_matrix[rownames(expr_matrix) %in% marker_genes, ]
    
    n_cells <- ncol(marker_expr)
    Q <- matrix(0, n_cells, n_cells)
    
    # Calculate pairwise distances using marker genes only
    for(i in 1:n_cells) {
      for(j in 1:n_cells) {
        if(i != j) {
          # Euclidean distance between cells based on marker genes
          dist_ij <- sqrt(sum((marker_expr[,i] - marker_expr[,j])^2))
          # Heat kernel weighting
          Q[i,j] <- exp(-dist_ij^2 / sigma)
        }
      }
    }
    
    return(Q)
  }
  
  # Get expression data
  ref_data <- as.matrix(GetAssayData(seu_ref, slot = "data"))
  test_data <- as.matrix(GetAssayData(seu_test, slot = "data"))
  
  # Find marker genes from reference data
  marker_genes <- find_marker_genes(seu_ref, n_markers = 20)
  print(paste("Found", length(marker_genes), "marker genes"))
  
  # Create weight matrix Q for reference data
  Q_ref <- create_weight_matrix(ref_data, marker_genes, sigma = 1)
  
  # Set parameters for rssNMF
  k_rss <- 20  # Number of factors
  alpha_rss <- 2  # Sparsity parameter
  beta_rss <- 2   # Graph regularization parameter
  iteration <- 100
  sigma <- 1e-6
  
  # Initialize W and H matrices for reference data
  initial_ref <- initial_NMF(ref_data, k_rss, "nndsvd")
  
  # Create eigens_set for the function (seems to be used for storing eigenvalues)
  eigens_set <- list()
  eigens_set[[1]] <- matrix(0, 1, ncol(Q_ref))
  
  print("Running rssNMF on reference data...")
  rss_result_ref <- graph_rssNMF(X = ref_data, 
                                 initial = initial_ref, 
                                 Q = Q_ref, 
                                 alpha = alpha_rss, 
                                 beta = beta_rss, 
                                 iteration = iteration, 
                                 sigma = sigma, 
                                 j = 1, 
                                 k = 1)
  
  W_ref <- rss_result_ref[[1]]
  H_ref <- rss_result_ref[[2]]
  
  # Add rssNMF results to reference Seurat object
  seu_ref[["rssNMF"]] <- CreateDimReducObject(
    embeddings = t(H_ref),
    key = "rssNMF_",
    assay = DefaultAssay(seu_ref)
  )
  
  # Run UMAP and clustering on reference
  seu_ref <- seu_ref %>%
    RunUMAP(dims = 1:k_rss, reduction = 'rssNMF', reduction.name = 'umap_rssNMF', verbose = F) %>%
    FindNeighbors(dims = 1:k_rss, reduction = 'rssNMF', verbose = F) %>%
    FindClusters(resolution = 0.2, verbose = F)
  
  seu_ref$rssNMF_clusters <- seu_ref$originalexp_snn_res.0.2
  
  # For test data, we need to project using the learned W matrix
  print("Projecting test data using learned W matrix...")
  
  # Function to project new data using learned W matrix
  project_rssNMF <- function(X_new, W_learned, max_iter = 100, tol = 1e-6) {
    # Initialize H for new data
    k <- ncol(W_learned)
    n_new <- ncol(X_new)
    H_new <- matrix(runif(k * n_new, 0, 1), nrow = k, ncol = n_new)
    
    # Project using non-negative least squares approach
    # Minimize ||X_new - W_learned * H_new||^2 subject to H_new >= 0
    
    for(iter in 1:max_iter) {
      H_old <- H_new
      
      # Update H using multiplicative updates
      numerator <- t(W_learned) %*% X_new
      denominator <- t(W_learned) %*% W_learned %*% H_new
      denominator[denominator < 1e-10] <- 1e-10  # Avoid division by zero
      
      H_new <- H_new * (numerator / denominator)
      
      # Check convergence
      if(norm(H_new - H_old, "F") < tol) {
        break
      }
    }
    
    return(H_new)
  }
  
  # Project test data
  H_test <- project_rssNMF(test_data, W_ref)
  
  # Add rssNMF results to test Seurat object
  seu_test[["rssNMF"]] <- CreateDimReducObject(
    embeddings = t(H_test),
    key = "rssNMF_",
    assay = DefaultAssay(seu_test)
  )
  
  # Run UMAP and clustering on test data
  seu_test <- seu_test %>%
    RunUMAP(dims = 1:k_rss, reduction = 'rssNMF', reduction.name = 'umap_rssNMF', verbose = F) %>%
    FindNeighbors(dims = 1:k_rss, reduction = 'rssNMF', verbose = F) %>%
    FindClusters(resolution = 0.2, verbose = F)
  
  seu_test$rssNMF_clusters <- seu_test$originalexp_snn_res.0.2
  
  # Use SingleR for annotation
  predictions <- SingleR::SingleR(
    test = as.SingleCellExperiment(seu_test),
    ref = as.SingleCellExperiment(seu_ref),
    labels = seu_ref$Group,
    de.n = 50,
    clusters = seu_test$rssNMF_clusters
  )
  
  seu_test$rssNMF_clusters_SingleR <- plyr::mapvalues(
    seu_test$rssNMF_clusters,
    from = rownames(predictions),
    to = predictions$labels
  )
  
  print('LIGER')
  
  library(rliger)
  
  # Prepare data for LIGER
  # LIGER expects raw count data, so we'll use the counts slot
  ref_counts <- as.matrix(GetAssayData(seu_ref, slot = "counts"))
  test_counts <- as.matrix(GetAssayData(seu_test, slot = "counts"))
  
  print("Creating LIGER object...")
  liger_obj <- createLiger(list(ref = ref_counts, test = test_counts))
  
  # Preprocessing steps
  print("Preprocessing data...")
  
  # Normalize data
  liger_obj <- normalize(liger_obj)
  
  # Select highly variable genes
  # This identifies genes that are variable across datasets
  liger_obj <- selectGenes(liger_obj, 
                           var.thresh = 0.1,  # Variance threshold
                           alpha.thresh = 0.99) # Alpha threshold for selection
  
  print(paste("Selected", length(liger_obj@varFeatures), "variable genes"))
  
  # Scale data - this is important for iNMF
  liger_obj <- scaleNotCenter(liger_obj)
  
  # Run integrative NMF (iNMF)
  # This is the core LIGER algorithm that finds shared (W) and dataset-specific (V) factors
  print("Running integrative NMF...")
  liger_obj <- runINMF(liger_obj, 
                       k = 20,           # Number of factors
                       lambda = 5.0,     # Regularization parameter for dataset-specific factors
                       nRandomStarts = 3, # Number of restarts
                       nIteration = 100)  # Maximum iterations
  
  # Quantile normalization to align datasets
  print("Performing quantile normalization...")
  liger_obj <- quantileNorm(liger_obj, 
                            reference = "ref",  # Use reference as the reference dataset
                            minCells = 10,      # Minimum cells per cluster
                            quantiles = 50,     # Number of quantiles
                            center = FALSE)
  
  # Run clustering on the aligned data
  print("Clustering aligned data...")
  liger_obj <- runCluster(liger_obj, 
                          method = "louvain",
                          resolution = 0.2,
                          seed = 1)
  
  # Extract results for reference data
  ref_indices <- which(liger_obj@cellMeta$dataset == "ref")
  ref_H_norm <- liger_obj@H.norm[ref_indices, ]  # Normalized H matrix for reference
  ref_clusters <- liger_obj@cellMeta$louvain_cluster[ref_indices]
  
  
  # Fix names if necessary
  rownames(ref_H_norm) <- colnames(seu_ref)
  # Add LIGER results to reference Seurat object
  seu_ref[["LIGER"]] <- CreateDimReducObject(
    embeddings = ref_H_norm,
    key = "LIGER_",
    assay = DefaultAssay(seu_ref)
  )
  
  # Run UMAP and update clustering for reference
  seu_ref <- seu_ref %>%
    RunUMAP(dims = 1:20, reduction = 'LIGER', reduction.name = 'umap_liger', verbose = F) %>%
    FindNeighbors(dims = 1:20, reduction = 'LIGER', verbose = F) %>%
    FindClusters(resolution = 0.2, verbose = F)
  
  seu_ref$LIGER_clusters <- seu_ref$originalexp_snn_res.0.2
  
  # Extract results for test data
  test_indices <- which(liger_obj@cellMeta$dataset =="test")
  test_H_norm <- liger_obj@H.norm[test_indices, ]  # Normalized H matrix for test
  test_clusters <- liger_obj@cellMeta$louvain_cluster[test_indices]
  
  rownames(test_H_norm) <- colnames(seu_test)
  # Add LIGER results to test Seurat object
  seu_test[["LIGER"]] <- CreateDimReducObject(
    embeddings = test_H_norm,
    key = "LIGER_",
    assay = DefaultAssay(seu_test)
  )
  
  # Run UMAP and clustering for test data
  seu_test <- seu_test %>%
    RunUMAP(dims = 1:20, reduction = 'LIGER', reduction.name = 'umap_liger', verbose = F) %>%
    FindNeighbors(dims = 1:20, reduction = 'LIGER', verbose = F) %>%
    FindClusters(resolution = 0.2, verbose = F)
  
  seu_test$LIGER_clusters <- seu_test$originalexp_snn_res.0.2
  
  # Use SingleR for annotation
  predictions <- SingleR::SingleR(
    test = as.SingleCellExperiment(seu_test),
    ref = as.SingleCellExperiment(seu_ref),
    labels = seu_ref$Group,
    de.n = 50,
    clusters = seu_test$LIGER_clusters
  )
  
  seu_test$LIGER_clusters_SingleR <- plyr::mapvalues(
    seu_test$LIGER_clusters,
    from = rownames(predictions),
    to = predictions$labels
  )
  
  print('CoGAPS')
  
  library(CoGAPS)
  
  # Extract raw count matrix from Seurat
  ref_counts <- as.matrix(GetAssayData(seu_ref, slot = "counts"))
  
  # Option 1: run CoGAPS directly on matrix
  gaps_result <- CoGAPS(ref_counts, nPatterns = 20, sparseOptimization = TRUE, nIterations = 500)
  
  # Option 2: if you want to use parameters more flexibly
  params <- new("CogapsParams")
  params <- setParam(params, "nPatterns", 20)
  params <- setParam(params, "nIterations", 500)
  params <- setParam(params, "sparseOptimization", TRUE)
  
  gaps_result <- CoGAPS(ref_counts, params = params)
  
  library(projectR)
  
  test_counts <- as.matrix(GetAssayData(seu_test, slot = "counts"))
  
  # Optional: ensure genes match
  common_genes <- intersect(rownames(ref_counts), rownames(test_counts))
  ref_amplitudes <- gaps_result@featureLoadings[common_genes, ]
  test_data <- test_counts[common_genes, ]
  
  # Project test data into CoGAPS-derived latent space
  projection_result <- projectR(
    data = test_data,
    loadings = ref_amplitudes,
    dataNames = rownames(test_data),
    loadingsNames = rownames(ref_amplitudes),
    full = TRUE
  )
  
  # Add projection to Seurat object
  seu_test[["scCoGAPS"]] <- CreateDimReducObject(
    embeddings = t(projection_result$projection),
    key = "scCoGAPS_",
    assay = DefaultAssay(seu_test)
  )
  
  # UMAP, clustering, etc.
  seu_test <- seu_test %>%
    RunUMAP(dims = 1:20, reduction = "scCoGAPS", reduction.name = "umap_scCoGAPS") %>%
    FindNeighbors(dims = 1:20, reduction = "scCoGAPS") %>%
    FindClusters(resolution = 0.2)
  
  seu_test$scCoGAPS_clusters <- seu_test$originalexp_snn_res.0.2
  
  predictions <- SingleR::SingleR(
    test = as.SingleCellExperiment(seu_test),
    ref = as.SingleCellExperiment(seu_ref),
    labels = seu_ref$Group,
    clusters = seu_test$scCoGAPS_clusters
  )
  
  seu_test$scCoGAPS_clusters_SingleR <- plyr::mapvalues(
    seu_test$scCoGAPS_clusters,
    from = rownames(predictions),
    to = predictions$labels
  )
  # Save final results
  qsave(seu_test, 
        paste0('save_data/simulation_6tools', SIM_NUMBER, '.qs'))
  
  return(list(seu_test = seu_test, 
              seu_ref = seu_ref, 
              optimal_params = optimal_params_meth,
              simulation_params = sim_params))
}

# Run 10 simulations
results_list <- list()
for(i in 1:10) {
  print(paste("Starting simulation", i, "of 10"))
  results_list[[i]] <- tryCatch({
    run_simulation(i, simulation_params)
  }, error = function(e) {
    print(paste("Error in simulation", i, ":", e))
    return(NULL)
  })
  
}

# Save final combined results
qsave(results_list, 
      paste0('save_data/final_all_simulations_results_', SUFFIX, format(Sys.time(), "%Y%m%d_%H%M"), '.qs'))

print("All simulations completed!")