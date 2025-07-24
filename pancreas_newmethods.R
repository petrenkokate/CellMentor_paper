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

SIM_NUMBER = '_baron'
set.seed(100)
my_colors <- c("#6495EDFF", "#FF69B4FF", "#BA55D3FF",  "#F08080FF", 
               "#9ACD32FF", "#4682B4FF",  "#DDA0DDFF", "#FFA07AFF", "#8FBC8BFF",
               "#40E0D0FF",  "#F0E68CFF", "#5F9EA0FF", "#D2B48CFF",  
               "#32CD32FF",  
               "#FFDAB9FF",  "#87CEEBFF")


print("Create ref")

reference =  h.baron_dataset()
test = muraro_dataset()
important_genes = c("GCG", "INS", "SST", "PPY", "GHRL", "KRT19", "CPA1", "PDGFRB", "PDGFRA", "RGS5",
                    "VWF", "PECAM1", "TPSAB1", "KIT", "CPA3", "TRAC", "SOX10", "SDS")

seu_ref <- CreateSeuratObject(counts = reference$data)
seu_ref$celltype <- reference$celltypes

# Standard preprocessing
seu_ref <- NormalizeData(seu_ref) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()
seu_ref <- RunUMAP(seu_ref, dims = 1:10) %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 0.1)

seu_ref$Seurat_clusters <- seu_ref$RNA_snn_res.0.1

print("Create test")
rownames(test$data) <- make.unique(rownames(test$data))
seu_test <- CreateSeuratObject(counts = test$data)
seu_test$celltype <- test$celltypes

# Standard preprocessing
seu_test <- NormalizeData(seu_test) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()

seu_test <- RunUMAP(seu_test, dims = 1:10) %>% 
  FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 0.2)

print('PCA')

seu_test$Seurat_clusters <- seu_test$RNA_snn_res.0.2

predictions <- SingleR::SingleR(
  test = as.SingleCellExperiment(seu_test),
  ref = as.SingleCellExperiment(seu_ref),
  labels = seu_ref$celltype,
  de.n = 50,
  clusters = seu_test$Seurat_clusters
)

seu_test$Seurat_clusters_SingleR <- plyr::mapvalues(
  seu_test$Seurat_clusters,
  from = rownames(predictions),
  to = predictions$labels)

# print('ICA')
# 
# library(fastICA)
# 
# # Get scaled data
# scaled_data <- GetAssayData(seu_test, slot = "scale.data")
# 
# # Run ICA
# ica_result <- fastICA(t(scaled_data), n.comp = 20)
# 
# # Create embedding matrix for Seurat
# ica_embeddings <- ica_result$S
# rownames(ica_embeddings) <- colnames(seu_test)
# colnames(ica_embeddings) <- paste0("IC_", 1:ncol(ica_embeddings))
# 
# # Add ICA as a reduction to Seurat object
# seu_test[["ica"]] <- CreateDimReducObject(
#   embeddings = ica_embeddings,
#   key = "IC_",
#   assay = DefaultAssay(seu_test)
# )
# 
# # Continue with your workflow
# seu_test <- seu_test %>%
#   RunUMAP(dims = 1:10, reduction = 'ica', reduction.name = 'umap_ica') %>%
#   FindNeighbors(dims = 1:10, reduction = 'ica') %>%
#   FindClusters(resolution = 0.2)
# 
# seu_test$ICA_clusters <- seu_test$RNA_snn_res.0.2
# 
# predictions <- SingleR::SingleR(
#   test = as.SingleCellExperiment(seu_test),
#   ref = as.SingleCellExperiment(seu_ref),
#   labels = seu_ref$celltype,
#   de.n = 50,
#   clusters = seu_test$ICA_clusters
# )
# 
# seu_test$ICA_clusters_SingleR <- plyr::mapvalues(
#   seu_test$ICA_clusters,
#   from = rownames(predictions),
#   to = predictions$labels)
# 
# print('NMF')
# 
# seu_test <- seu_test %>% 
#   NormalizeData()
# 
# # Get normalized data and convert to dense matrix
# data_matrix <- as.matrix(GetAssayData(seu_test, slot = "data"))
# 
# # Remove rows with all zeros or NAs
# keep_rows <- rowSums(data_matrix) > 0 & !apply(data_matrix, 1, function(x) any(is.na(x)))
# data_matrix <- data_matrix[keep_rows,]
# 
# # Run NNMF
# res_nmf <- nnmf(data_matrix, k = 20, method = "lee", loss = "mse", n.threads = 30)
# 
# # Add NMF as reduction (use H matrix for cell embeddings)
# seu_test[["nmf"]] <- CreateDimReducObject(
#   embeddings = t(res_nmf$H), 
#   key = "NMF_",
#   assay = DefaultAssay(seu_test)
# )
# 
# seu_test <- seu_test %>%
#   RunUMAP(dims = 1:20, reduction = 'nmf', reduction.name = 'umap_nmf') %>%
#   FindNeighbors(dims = 1:20, reduction = 'nmf') %>%
#   FindClusters(resolution = 0.1)
# 
# seu_test$NMF_clusters <- seu_test$RNA_snn_res.0.1
# 
# predictions <- SingleR::SingleR(
#   test = as.SingleCellExperiment(seu_test),
#   ref = as.SingleCellExperiment(seu_ref),
#   labels = seu_ref$celltype,
#   de.n = 50,
#   clusters = seu_test$NMF_clusters
# )
# 
# seu_test$NMF_clusters_SingleR <- plyr::mapvalues(
#   seu_test$NMF_clusters,
#   from = rownames(predictions),
#   to = predictions$labels)
# 
# print('pCMF')
# 
# # Get normalized data and transpose to cells x genes
# data_matrix <- t(as.matrix(GetAssayData(seu_test, slot = "counts")))
# 
# # Pre-filter genes
# kept_cols <- prefilter(data_matrix, prop = 0.05, quant_max = 0.95, 
#                        presel = TRUE, threshold = 0.2)
# data_matrix <- data_matrix[,kept_cols]
# 
# # Run pCMF
# res_pcmf <- pCMF(data_matrix, K = 20, verbose = FALSE, 
#                  zero_inflation = TRUE, sparsity = TRUE, 
#                  ncores = 30)
# 
# embeddings <- res_pcmf$factor$U
# rownames(embeddings) <- Cells(seu_test)
# colnames(embeddings) <- paste0("pCMF_", 1:ncol(embeddings))
# 
# # Create dimension reduction object
# seu_test[["pCMF"]] <- CreateDimReducObject(
#   embeddings = embeddings,
#   key = "pCMF_",
#   assay = DefaultAssay(seu_test)
# )
# 
# # Continue with workflow
# seu_test <- seu_test %>%
#   RunUMAP(dims = 1:10, reduction = 'pCMF', reduction.name = 'umap_pcmf', verbose = F) %>%
#   FindNeighbors(dims = 1:20, reduction = 'pCMF', verbose = F) %>%
#   FindClusters(resolution = 0.3, verbose = F)
# 
# seu_test$pCMF_clusters <- seu_test$RNA_snn_res.0.3
# 
# predictions <- SingleR::SingleR(
#   test = as.SingleCellExperiment(seu_test),
#   ref = as.SingleCellExperiment(seu_ref),
#   labels = seu_ref$celltype,
#   de.n = 50,
#   clusters = seu_test$pCMF_clusters
# )
# 
# seu_test$pCMF_clusters_SingleR <- plyr::mapvalues(
#   seu_test$pCMF_clusters,
#   from = rownames(predictions),
#   to = predictions$labels)
# 
# print('GLM-PCA')
# 
# # Get count data
# data_matrix <- as.matrix(GetAssayData(seu_test, slot = "counts"))
# 
# # Remove rows with all zeros
# keep_rows <- rowSums(data_matrix) > 0
# data_matrix <- data_matrix[keep_rows,]
# 
# # Run GLMPCA
# res <- glmpca(data_matrix, L = 20)
# 
# # Add rownames/colnames to factors matrix
# factors <- res$factors
# rownames(factors) <- Cells(seu_test)
# colnames(factors) <- paste0("GBM_", 1:ncol(factors))
# 
# # Add GLMPCA as reduction
# seu_test[["gbm"]] <- CreateDimReducObject(
#   embeddings = as.matrix(factors),
#   key = "GBM_",
#   assay = DefaultAssay(seu_test)
# )
# 
# # Continue with workflow
# seu_test <- seu_test %>%
#   RunUMAP(dims = 1:10, reduction = 'gbm', reduction.name = 'umap_gbm', verbose = F) %>%
#   FindNeighbors(dims = 1:20, reduction = 'gbm', verbose = F) %>%
#   FindClusters(resolution = 0.2, verbose = F)
# 
# seu_test$GBM_clusters <- seu_test$RNA_snn_res.0.2
# 
# predictions <- SingleR::SingleR(
#   test = as.SingleCellExperiment(seu_test),
#   ref = as.SingleCellExperiment(seu_ref),
#   labels = seu_ref$celltype,
#   de.n = 50,
#   clusters = seu_test$GBM_clusters
# )
# 
# seu_test$GBM_clusters_SingleR <- plyr::mapvalues(
#   seu_test$GBM_clusters,
#   from = rownames(predictions),
#   to = predictions$labels)
# 
# print('CellMentor')
# 
# object = CreateCSFNMFobject(reference$data, reference$celltypes, test$data, 
#                             gene_list = important_genes)
# 
# file_path <- paste0('save_data/simulation_method', SIM_NUMBER, '.qs')
# # Check if file exists
# if (file.exists(file_path)) {
#   print("BEST METHOD ALREADY EXISTS. LOADING FILE...")
#   optimal_params_meth <- qread(file_path)
# } else {
#   print("SELECT BEST METHOD")
#   optimal_params_meth <- CellMentor(object, num_cores = 20, init_methods = c('regulated'), gamma_range = c(1), delta_range = c(0.1))
#   qsave(optimal_params_meth, file_path)
# }
# K_VALUE = optimal_params_meth$best_params$k
# final_model <- optimal_params_meth$best_model
# 
# h_test <- project_data(
#   W = final_model@W,                    # Use learned W matrix
#   X = final_model@matrices@data,        # Query/test data matrix
#   seed = 1,
#   num_cores = 10,
#   chunk_size = NULL,
#   verbose = TRUE
# )
# 
# seu_test$CellMentor <- CreateDimReducObject(
#   embeddings = t(as.matrix(h_test)),
#   key = paste0('CellMentor', "_"),
#   assay = DefaultAssay(seu_test), 
#   loadings = as.matrix(final_model@W)
# ) 
# 
# seu_test <- seu_test %>% 
#   RunUMAP(reduction = 'CellMentor', dims= 1:K_VALUE, reduction.name = 'umap_cellmentor', verbose = F) %>% 
#   FindNeighbors(dims = 1:K_VALUE, reduction = 'CellMentor', verbose = F) %>% 
#   FindClusters(resolution = 0.1, verbose = F)
# 
# seu_test$CellMentor_clusters <- seu_test$RNA_snn_res.0.1
# 
# predictions <- SingleR::SingleR(
#   test = as.matrix(h_test),
#   ref = as.matrix(final_model@H),
#   labels = final_model@annotation$celltype,
#   de.n = 50,
#   clusters = seu_test$CellMentor_clusters
# )
# 
# seu_test$CellMentor_clusters_SingleR <- plyr::mapvalues(
#   seu_test$CellMentor_clusters,
#   from = rownames(predictions),
#   to = predictions$labels)
# 
# print('READ OBJECT')
# 
# seu_test2 <- seu_test
# seu_test2$celltype <- 'Unknown'
# seu_test2$Dataset <- 'test'
# seu_ref$Dataset <- 'ref'
# library(reticulate)
# library(sceasy)
# # Ensure scvi-tools is installed
# if (!py_module_available("scvi")) {
#   py_install("scvi-tools")
# }
# scvi <- import("scvi")
# 
# print('PREPROCESSING for SCANVI')
# # Prepare data
# seu_ref[["RNA"]] <- as(seu_ref[["RNA"]], "Assay")
# adata <- convertFormat(seu_ref, from="seurat", to="anndata", 
#                        main_layer="counts", drop_single_values=FALSE)
# adata$X = adata$X$tocsr()
# scvi$model$SCVI$setup_anndata(adata, batch_key = NULL, labels_key = "celltype")
# model = scvi$model$SCVI(adata)
# 
# print('SCANVI TRAIN')
# model$train()
# 
# scanvi_model = scvi$model$SCANVI$from_scvi_model(
#   model,
#   unlabeled_category="Unknown"
# )
# 
# scanvi_model$train()
# 
# print('SCANVI DONE')
# seu_test2[["RNA"]] <- as(seu_test2[["RNA"]], "Assay")
# adata_test <- convertFormat(seu_test2, from="seurat", to="anndata", 
#                             main_layer="counts", drop_single_values=FALSE)
# adata_test$X = adata_test$X$tocsr()
# scvi$model$SCVI$prepare_query_anndata(adata_test, scanvi_model)
# latent_scanvi = scanvi_model$get_latent_representation(adata_test)
# 
# latent_scanvi <- as.matrix(latent_scanvi)
# rownames(latent_scanvi) = colnames(seu_test)
# seu_test[["scanvi"]] <- CreateDimReducObject(embeddings = latent_scanvi, key = "scanvi_", assay = DefaultAssay(seu_test))
# 
# seu_test <- seu_test %>% 
#   RunUMAP(reduction = 'scanvi', dims= 1:10, reduction.name = 'umap_scanvi', verbose = F) %>%
#   FindNeighbors(dims = 1:10, reduction = 'scanvi', verbose = F) %>%
#   FindClusters(resolution = 0.2, verbose = F)
# 
# seu_test$scanvi_clusters <- seu_test$RNA_snn_res.0.2
# 
# predictions <- SingleR::SingleR(
#   test = as.SingleCellExperiment(seu_test),
#   ref = as.SingleCellExperiment(seu_ref),
#   labels = seu_ref$celltype,
#   clusters = seu_test$scanvi_clusters
# )
# 
# seu_test$scanvi_clusters_SingleR <- plyr::mapvalues(
#   seu_test$scanvi_clusters,
#   from = rownames(predictions),
#   to = predictions$labels)

source("~/nmf/paper_scripts/code/graph_rssNMF.R")
source("~/nmf/paper_scripts/code/initial_NMF.R")

seu_ref$Group <- seu_ref$celltype
seu_test$Group <- seu_test$celltype
library(qs)

print('rssNMF')

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

# Create eigens_set for the function (seems to be used for storing eigenvalues)
eigens_set <<- list()  # Note the <<- for global assignment
eigens_set[[1]] <<- matrix(0, 1, ncol(Q_ref))

# Set parameters for rssNMF
k_rss <- 20  # Number of factors
alpha_rss <- 2  # Sparsity parameter
beta_rss <- 2   # Graph regularization parameter
iteration <- 100
sigma <- 1e-6

# Initialize W and H matrices for reference data
initial_ref <- initial_NMF(ref_data, k_rss, "nndsvd")


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

print('LIGER - with proper matrix formatting')

# Ensure matrices are in the correct format for LIGER
ref_counts <- as.matrix(GetAssayData(seu_ref, slot = "counts"))
test_counts <- as.matrix(GetAssayData(seu_test, slot = "counts"))

# Remove genes that are all zeros
ref_counts <- ref_counts[rowSums(ref_counts) > 0, ]
test_counts <- test_counts[rowSums(test_counts) > 0, ]
common_genes <- intersect(rownames(ref_counts), rownames(test_counts))
ref_counts <- ref_counts[common_genes, ]
test_counts <- test_counts[common_genes, ]

# CRITICAL: Convert to dense matrix and ensure integer type
ref_counts <- as.matrix(ref_counts)
test_counts <- as.matrix(test_counts)

# Ensure data is in integer format (LIGER expects count data)
storage.mode(ref_counts) <- "integer"
storage.mode(test_counts) <- "integer"

# Additional filtering - remove genes with very low counts
min_total_count <- 10
ref_counts <- ref_counts[rowSums(ref_counts) >= min_total_count, ]
test_counts <- test_counts[rowSums(test_counts) >= min_total_count, ]

# Get common genes again after filtering
common_genes <- intersect(rownames(ref_counts), rownames(test_counts))
ref_counts <- ref_counts[common_genes, ]
test_counts <- test_counts[common_genes, ]

print(paste("Final dimensions - Ref:", nrow(ref_counts), "x", ncol(ref_counts)))
print(paste("Final dimensions - Test:", nrow(test_counts), "x", ncol(test_counts)))
print(paste("Data type - Ref:", typeof(ref_counts)))
print(paste("Data type - Test:", typeof(test_counts)))

liger_obj <- createLiger(list("ref" = ref_counts, "test" = test_counts))

print("Preprocessing data...")
liger_obj <- normalize(liger_obj)

# CRITICAL: Much more stringent gene selection
# The fact that you're selecting 14927/14935 genes is likely the problem
liger_obj <- selectGenes(liger_obj,
                         var.thresh = 2.0,    # Much higher threshold
                         alpha.thresh = 0.8)   # More stringent

print(paste("Selected", length(liger_obj@varFeatures), "variable genes"))

# Only proceed if we have reasonable number of genes
# Target: should be <2000 genes
if (length(liger_obj@varFeatures) > 2000) {
  print("Still too many genes, trying even higher threshold...")
  liger_obj <- selectGenes(liger_obj,
                           var.thresh = 2.5,    # Even higher
                           alpha.thresh = 0.7)
  print(paste("After second filter:", length(liger_obj@varFeatures), "variable genes"))
}

# Only proceed if reasonable number
if (length(liger_obj@varFeatures) > 3000) {
  stop(paste("Still too many variable genes:", length(liger_obj@varFeatures), "- LIGER will crash"))
}

if (length(liger_obj@varFeatures) < 100) {
  stop("Too few variable genes - lower the threshold")
}

liger_obj <- scaleNotCenter(liger_obj)

# Ultra-conservative iNMF parameters
print("Running integrative NMF with minimal parameters...")
liger_obj <- runINMF(liger_obj,
                     k = 20,                # Very small k
                     lambda = 100.0,       # Very high regularization
                     nRandomStarts = 1,    # Single start
                     nIteration = 10,      # Very few iterations
                     thresh = 1e-3,        # High convergence threshold
                     seed = 42)

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

# Get the actual number of LIGER dimensions
k_liger <- ncol(seu_ref[["LIGER"]]@cell.embeddings)
print(paste("LIGER produced", k_liger, "dimensions"))

# Use the correct number of dimensions
seu_ref <- seu_ref %>%
  RunUMAP(dims = 1:k_liger, reduction = 'LIGER', reduction.name = 'umap_liger', verbose = F) %>%
  FindNeighbors(dims = 1:k_liger, reduction = 'LIGER', verbose = F) %>%
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
  RunUMAP(dims = 1:k_liger, reduction = 'LIGER', reduction.name = 'umap_liger', verbose = F) %>%
  FindNeighbors(dims = 1:k_liger, reduction = 'LIGER', verbose = F) %>%
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

print('CASSL')

# Helper function to find Python executable (prioritizing reticulate environment)
find_python_for_cassl <- function() {
  if (requireNamespace("reticulate", quietly = TRUE)) {
    tryCatch({
      python_path <- reticulate::py_config()$python
      if (!is.null(python_path) && file.exists(python_path)) {
        test_cmd <- sprintf("%s -c 'import numpy, pandas, sklearn; print(\"OK\")'", python_path)
        result <- suppressWarnings(system(test_cmd, intern = TRUE, ignore.stderr = TRUE))
        if (length(result) > 0 && grepl("OK", result[1])) {
          return(python_path)
        }
      }
    }, error = function(e) NULL)
  }
  
  python_commands <- c("python3.10", "python3", "python", "python3.11", "python3.9", "python3.8")
  for (cmd in python_commands) {
    test_cmd <- sprintf("%s -c 'import numpy, pandas, sklearn; print(\"OK\")'", cmd)
    result <- suppressWarnings(system(test_cmd, intern = TRUE, ignore.stderr = TRUE))
    if (length(result) > 0 && grepl("OK", result[1])) {
      return(cmd)
    }
  }
  stop("No suitable Python executable found with required packages")
}

# Function to run CASSL method
run_cassl_method <- function(seu_ref, seu_test, python_script = "cassl_method.py", temp_dir = tempdir()) {
  
  # Find Python executable
  python_cmd <- find_python_for_cassl()
  
  # Create temporary files
  ref_data_file <- file.path(temp_dir, paste0("ref_data_", Sys.getpid(), ".csv"))
  ref_labels_file <- file.path(temp_dir, paste0("ref_labels_", Sys.getpid(), ".csv"))
  test_data_file <- file.path(temp_dir, paste0("test_data_", Sys.getpid(), ".csv"))
  output_file <- file.path(temp_dir, paste0("cassl_results_", Sys.getpid(), ".pkl"))
  
  tryCatch({
    # Extract and save reference data (transpose to get cells x genes)
    ref_data <- t(as.matrix(GetAssayData(seu_ref, slot = "data")))
    ref_labels <- data.frame(labels = seu_ref$Group)
    test_data <- t(as.matrix(GetAssayData(seu_test, slot = "data")))
    
    # Find common genes
    common_genes <- intersect(colnames(ref_data), colnames(test_data))
    ref_data <- ref_data[, common_genes]
    test_data <- test_data[, common_genes]
    
    # Save files
    write.csv(ref_data, ref_data_file)
    write.csv(ref_labels, ref_labels_file)
    write.csv(test_data, test_data_file)
    
    # Run Python CASSL script
    cmd <- sprintf("%s %s --ref_data %s --ref_labels %s --test_data %s --output %s --p_missing 0.5 --max_dim 15",
                   python_cmd, python_script, ref_data_file, ref_labels_file, test_data_file, output_file)
    
    system_result <- system(cmd, intern = FALSE, wait = TRUE)
    
    if (!file.exists(output_file)) {
      stop("CASSL failed to generate output file")
    }
    
    # Load results using reticulate
    library(reticulate)
    py_run_string("import pickle")
    py_run_string(sprintf("
with open(r'%s', 'rb') as f:
    cassl_results = pickle.load(f)
", output_file))
    
    results <- py$cassl_results
    return(results)
    
  }, finally = {
    # Clean up temporary files
    unlink(c(ref_data_file, ref_labels_file, test_data_file, output_file))
  })
}

# Check if CASSL script exists
python_script_path <- "cassl_method.py"
# Run CASSL
cassl_results <- run_cassl_method(seu_ref, seu_test, python_script_path)

# Extract embeddings
ref_embedding <- cassl_results$ref_embedding
test_embedding <- cassl_results$test_embedding
optimal_k <- cassl_results$optimal_dim

# Set proper row and column names
rownames(ref_embedding) <- Cells(seu_ref)
colnames(ref_embedding) <- paste0("CASSL_", 1:ncol(ref_embedding))

rownames(test_embedding) <- Cells(seu_test)
colnames(test_embedding) <- paste0("CASSL_", 1:ncol(test_embedding))

# Add to Seurat objects
seu_ref[["cassl"]] <- CreateDimReducObject(
  embeddings = ref_embedding,
  key = "CASSL_",
  assay = DefaultAssay(seu_ref)
)

seu_test[["cassl"]] <- CreateDimReducObject(
  embeddings = test_embedding,
  key = "CASSL_",
  assay = DefaultAssay(seu_test)
)

# Continue with standard pipeline
seu_test <- seu_test %>%
  RunUMAP(dims = 1:optimal_k, reduction = 'cassl', 
          reduction.name = 'umap_cassl', verbose = F) %>%
  FindNeighbors(dims = 1:optimal_k, reduction = 'cassl', verbose = F) %>%
  FindClusters(resolution = 0.2, verbose = F)

seu_test$cassl_clusters <- seu_test$originalexp_snn_res.0.2

qsave(seu_test, paste0('save_data/6tools', SIM_NUMBER, '_scanvi.qs'))
