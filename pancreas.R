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

print('ICA')

library(fastICA)

# Get scaled data
scaled_data <- GetAssayData(seu_test, slot = "scale.data")

# Run ICA
ica_result <- fastICA(t(scaled_data), n.comp = 20)

# Create embedding matrix for Seurat
ica_embeddings <- ica_result$S
rownames(ica_embeddings) <- colnames(seu_test)
colnames(ica_embeddings) <- paste0("IC_", 1:ncol(ica_embeddings))

# Add ICA as a reduction to Seurat object
seu_test[["ica"]] <- CreateDimReducObject(
  embeddings = ica_embeddings,
  key = "IC_",
  assay = DefaultAssay(seu_test)
)

# Continue with your workflow
seu_test <- seu_test %>%
  RunUMAP(dims = 1:10, reduction = 'ica', reduction.name = 'umap_ica') %>%
  FindNeighbors(dims = 1:10, reduction = 'ica') %>%
  FindClusters(resolution = 0.2)

seu_test$ICA_clusters <- seu_test$RNA_snn_res.0.2

predictions <- SingleR::SingleR(
  test = as.SingleCellExperiment(seu_test),
  ref = as.SingleCellExperiment(seu_ref),
  labels = seu_ref$celltype,
  de.n = 50,
  clusters = seu_test$ICA_clusters
)

seu_test$ICA_clusters_SingleR <- plyr::mapvalues(
  seu_test$ICA_clusters,
  from = rownames(predictions),
  to = predictions$labels)

print('NMF')

seu_test <- seu_test %>% 
  NormalizeData()

# Get normalized data and convert to dense matrix
data_matrix <- as.matrix(GetAssayData(seu_test, slot = "data"))

# Remove rows with all zeros or NAs
keep_rows <- rowSums(data_matrix) > 0 & !apply(data_matrix, 1, function(x) any(is.na(x)))
data_matrix <- data_matrix[keep_rows,]

# Run NNMF
res_nmf <- nnmf(data_matrix, k = 20, method = "lee", loss = "mse", n.threads = 30)

# Add NMF as reduction (use H matrix for cell embeddings)
seu_test[["nmf"]] <- CreateDimReducObject(
  embeddings = t(res_nmf$H), 
  key = "NMF_",
  assay = DefaultAssay(seu_test)
)

seu_test <- seu_test %>%
  RunUMAP(dims = 1:20, reduction = 'nmf', reduction.name = 'umap_nmf') %>%
  FindNeighbors(dims = 1:20, reduction = 'nmf') %>%
  FindClusters(resolution = 0.1)

seu_test$NMF_clusters <- seu_test$RNA_snn_res.0.1

predictions <- SingleR::SingleR(
  test = as.SingleCellExperiment(seu_test),
  ref = as.SingleCellExperiment(seu_ref),
  labels = seu_ref$celltype,
  de.n = 50,
  clusters = seu_test$NMF_clusters
)

seu_test$NMF_clusters_SingleR <- plyr::mapvalues(
  seu_test$NMF_clusters,
  from = rownames(predictions),
  to = predictions$labels)

print('pCMF')

# Get normalized data and transpose to cells x genes
data_matrix <- t(as.matrix(GetAssayData(seu_test, slot = "counts")))

# Pre-filter genes
kept_cols <- prefilter(data_matrix, prop = 0.05, quant_max = 0.95, 
                       presel = TRUE, threshold = 0.2)
data_matrix <- data_matrix[,kept_cols]

# Run pCMF
res_pcmf <- pCMF(data_matrix, K = 20, verbose = FALSE, 
                 zero_inflation = TRUE, sparsity = TRUE, 
                 ncores = 30)

embeddings <- res_pcmf$factor$U
rownames(embeddings) <- Cells(seu_test)
colnames(embeddings) <- paste0("pCMF_", 1:ncol(embeddings))

# Create dimension reduction object
seu_test[["pCMF"]] <- CreateDimReducObject(
  embeddings = embeddings,
  key = "pCMF_",
  assay = DefaultAssay(seu_test)
)

# Continue with workflow
seu_test <- seu_test %>%
  RunUMAP(dims = 1:10, reduction = 'pCMF', reduction.name = 'umap_pcmf', verbose = F) %>%
  FindNeighbors(dims = 1:20, reduction = 'pCMF', verbose = F) %>%
  FindClusters(resolution = 0.3, verbose = F)

seu_test$pCMF_clusters <- seu_test$RNA_snn_res.0.3

predictions <- SingleR::SingleR(
  test = as.SingleCellExperiment(seu_test),
  ref = as.SingleCellExperiment(seu_ref),
  labels = seu_ref$celltype,
  de.n = 50,
  clusters = seu_test$pCMF_clusters
)

seu_test$pCMF_clusters_SingleR <- plyr::mapvalues(
  seu_test$pCMF_clusters,
  from = rownames(predictions),
  to = predictions$labels)

print('GLM-PCA')

# Get count data
data_matrix <- as.matrix(GetAssayData(seu_test, slot = "counts"))

# Remove rows with all zeros
keep_rows <- rowSums(data_matrix) > 0
data_matrix <- data_matrix[keep_rows,]

# Run GLMPCA
res <- glmpca(data_matrix, L = 20)

# Add rownames/colnames to factors matrix
factors <- res$factors
rownames(factors) <- Cells(seu_test)
colnames(factors) <- paste0("GBM_", 1:ncol(factors))

# Add GLMPCA as reduction
seu_test[["gbm"]] <- CreateDimReducObject(
  embeddings = as.matrix(factors),
  key = "GBM_",
  assay = DefaultAssay(seu_test)
)

# Continue with workflow
seu_test <- seu_test %>%
  RunUMAP(dims = 1:10, reduction = 'gbm', reduction.name = 'umap_gbm', verbose = F) %>%
  FindNeighbors(dims = 1:20, reduction = 'gbm', verbose = F) %>%
  FindClusters(resolution = 0.2, verbose = F)

seu_test$GBM_clusters <- seu_test$RNA_snn_res.0.2

predictions <- SingleR::SingleR(
  test = as.SingleCellExperiment(seu_test),
  ref = as.SingleCellExperiment(seu_ref),
  labels = seu_ref$celltype,
  de.n = 50,
  clusters = seu_test$GBM_clusters
)

seu_test$GBM_clusters_SingleR <- plyr::mapvalues(
  seu_test$GBM_clusters,
  from = rownames(predictions),
  to = predictions$labels)

print('CellMentor')

object = CreateCSFNMFobject(reference$data, reference$celltypes, test$data, 
                            gene_list = important_genes)

file_path <- paste0('save_data/simulation_method', SIM_NUMBER, '.qs')
# Check if file exists
if (file.exists(file_path)) {
  print("BEST METHOD ALREADY EXISTS. LOADING FILE...")
  optimal_params_meth <- qread(file_path)
} else {
  print("SELECT BEST METHOD")
  optimal_params_meth <- CellMentor(object, num_cores = 20, init_methods = c('regulated'), gamma_range = c(1), delta_range = c(0.1))
  qsave(optimal_params_meth, file_path)
}
K_VALUE = optimal_params_meth$best_params$k
final_model <- optimal_params_meth$best_model

h_test <- project_data(
  W = final_model@W,                    # Use learned W matrix
  X = final_model@matrices@data,        # Query/test data matrix
  seed = 1,
  num_cores = 10,
  chunk_size = NULL,
  verbose = TRUE
)

seu_test$CellMentor <- CreateDimReducObject(
  embeddings = t(as.matrix(h_test)),
  key = paste0('CellMentor', "_"),
  assay = DefaultAssay(seu_test), 
  loadings = as.matrix(final_model@W)
) 

seu_test <- seu_test %>% 
  RunUMAP(reduction = 'CellMentor', dims= 1:K_VALUE, reduction.name = 'umap_cellmentor', verbose = F) %>% 
  FindNeighbors(dims = 1:K_VALUE, reduction = 'CellMentor', verbose = F) %>% 
  FindClusters(resolution = 0.1, verbose = F)

seu_test$CellMentor_clusters <- seu_test$RNA_snn_res.0.1

predictions <- SingleR::SingleR(
  test = as.matrix(h_test),
  ref = as.matrix(final_model@H),
  labels = final_model@annotation$celltype,
  de.n = 50,
  clusters = seu_test$CellMentor_clusters
)

seu_test$CellMentor_clusters_SingleR <- plyr::mapvalues(
  seu_test$CellMentor_clusters,
  from = rownames(predictions),
  to = predictions$labels)

print('READ OBJECT')

seu_test2 <- seu_test
seu_test2$celltype <- 'Unknown'
seu_test2$Dataset <- 'test'
seu_ref$Dataset <- 'ref'
library(reticulate)
library(sceasy)
# Ensure scvi-tools is installed
if (!py_module_available("scvi")) {
  py_install("scvi-tools")
}
scvi <- import("scvi")

print('PREPROCESSING for SCANVI')
# Prepare data
seu_ref[["RNA"]] <- as(seu_ref[["RNA"]], "Assay")
adata <- convertFormat(seu_ref, from="seurat", to="anndata", 
                       main_layer="counts", drop_single_values=FALSE)
adata$X = adata$X$tocsr()
scvi$model$SCVI$setup_anndata(adata, batch_key = NULL, labels_key = "celltype")
model = scvi$model$SCVI(adata)

print('SCANVI TRAIN')
model$train()

scanvi_model = scvi$model$SCANVI$from_scvi_model(
  model,
  unlabeled_category="Unknown"
)

scanvi_model$train()

print('SCANVI DONE')
seu_test2[["RNA"]] <- as(seu_test2[["RNA"]], "Assay")
adata_test <- convertFormat(seu_test2, from="seurat", to="anndata", 
                            main_layer="counts", drop_single_values=FALSE)
adata_test$X = adata_test$X$tocsr()
scvi$model$SCVI$prepare_query_anndata(adata_test, scanvi_model)
latent_scanvi = scanvi_model$get_latent_representation(adata_test)

latent_scanvi <- as.matrix(latent_scanvi)
rownames(latent_scanvi) = colnames(seu_test)
seu_test[["scanvi"]] <- CreateDimReducObject(embeddings = latent_scanvi, key = "scanvi_", assay = DefaultAssay(seu_test))

seu_test <- seu_test %>% 
  RunUMAP(reduction = 'scanvi', dims= 1:10, reduction.name = 'umap_scanvi', verbose = F) %>%
  FindNeighbors(dims = 1:10, reduction = 'scanvi', verbose = F) %>%
  FindClusters(resolution = 0.2, verbose = F)

seu_test$scanvi_clusters <- seu_test$RNA_snn_res.0.2

predictions <- SingleR::SingleR(
  test = as.SingleCellExperiment(seu_test),
  ref = as.SingleCellExperiment(seu_ref),
  labels = seu_ref$celltype,
  clusters = seu_test$scanvi_clusters
)

seu_test$scanvi_clusters_SingleR <- plyr::mapvalues(
  seu_test$scanvi_clusters,
  from = rownames(predictions),
  to = predictions$labels)

qsave(seu_test, paste0('save_data/6tools', SIM_NUMBER, '_scanvi.qs'))
