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
library(harmony)

SIM_NUMBER = '_newcelltype_merge2_6'
set.seed(100)
my_colors <- c("#6495EDFF", "#FF69B4FF", "#BA55D3FF",  "#F08080FF", 
               "#9ACD32FF", "#4682B4FF",  "#DDA0DDFF", "#FFA07AFF", "#8FBC8BFF",
               "#40E0D0FF",  "#F0E68CFF", "#5F9EA0FF", "#D2B48CFF",  
               "#32CD32FF",  
               "#FFDAB9FF",  "#87CEEBFF")

print("START SIMULATION")

params <- newSplatParams()
params <- setParams(params,
                    nGenes = 20000, 
                    seed = 100, 
                    batchCells = c(1000, 1000, 1000, 1000),
                    batch.facLoc = 0.1,
                    batch.facScale = 0.1,
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

sim <- splatSimulate(params, method = 'groups')
sim <- logNormCounts(sim)

seu <- as.Seurat(sim)

seu <- seu %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose=FALSE) %>% 
  RunUMAP(dims = 1:10)

print("SIMULATION DONE!")

print("Create ref")
seu$Group <- as.character(seu$Group)
seu$celltype <- seu$Group
seu_ref <- subset(seu, Batch %in% c('Batch1', 'Batch2'))

seu_ref$Group[seu_ref$Group == 'Group6'] <- 'Group2'

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


print("CreateCSFNMFobject !!!")
object = CreateCSFNMFobject(seu_ref@assays$originalexp$counts, seu_ref$Group, seu_test@assays$originalexp$counts)

file_path <- paste0('save_data/simulation_method', SIM_NUMBER, '_group_merge.qs')
if (file.exists(file_path)) {
  print("BEST METHOD ALREADY EXISTS. LOADING FILE...")
  optimal_params_meth <- qread(file_path)
} else {
  print("SELECT BEST METHOD")
  optimal_params_meth <- CellMentor(object, 
                                                   num_cores = 10, 
                                                   init_methods = c('regulated'),
                                                   gamma_range = c(1), 
                                                   delta_range = c(1))

  
  print("SAVE BEST METHOD")
  qsave(optimal_params_meth, file_path)
}

K_VALUE = optimal_params_meth$best_params$k

print("BEST METHOD !!!")
final_model <- optimal_params_meth$best_model

print("PREDICTION PROJECTION")
h_test <- project_data(
  W = final_model@W,                    # Use learned W matrix
  X = final_model@matrices@data,        # Query/test data matrix
  seed = 1,
  num_cores = 10,
  chunk_size = NULL,
  verbose = TRUE
)

print("TRAIN")
seu_ref$CellMentor <- CreateDimReducObject(
  embeddings = t(as.matrix(final_model@H)),
  key = paste0('CellMentor', "_"),
  assay = DefaultAssay(seu_ref)
)

seu_ref <- seu_ref %>%
  RunUMAP(reduction = 'CellMentor', dims= 1:K_VALUE, reduction.name = 'umap_cellmentor', verbose = F) %>%
  FindNeighbors(dims = 1:K_VALUE, reduction = 'CellMentor', verbose = F) %>%
  FindClusters(resolution = 0.2, verbose = F)

seu_ref$CellMentor_clusters <- seu_ref$originalexp_snn_res.0.2

print('PCA')
predictions <- SingleR::SingleR(
  test = as.SingleCellExperiment(seu_test),
  ref = as.SingleCellExperiment(seu_ref),
  labels = seu_ref$Group,
  de.n = 50,
  clusters = seu_test$Seurat_clusters
)

seu_test$Seurat_clusters_SingleR <- plyr::mapvalues(
  seu_test$Seurat_clusters,
  from = rownames(predictions),
  to = predictions$labels)

print('Harmony')

seu_test <- seu_test %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose=FALSE) %>%
  RunHarmony(., 'Batch',
             lambda = 1, verbose = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:20, verbose=F, reduction.name = 'umap_harmony') %>%
  FindNeighbors(dims = 1:10, reduction = 'harmony') %>%
  FindClusters(resolution = 0.2)

seu_test$harmony_clusters <- seu_test$originalexp_snn_res.0.2

predictions <- SingleR::SingleR(
  test = as.SingleCellExperiment(seu_test),
  ref = as.SingleCellExperiment(seu_ref),
  labels = seu_ref$Group,
  de.n = 50,
  clusters = seu_test$harmony_clusters
)

seu_test$harmony_clusters_SingleR <- plyr::mapvalues(
  seu_test$harmony_clusters,
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

seu_test$ICA_clusters <- seu_test$originalexp_snn_res.0.2

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
  RunUMAP(dims = 1:10, reduction = 'nmf', reduction.name = 'umap_nmf') %>%
  FindNeighbors(dims = 1:10, reduction = 'nmf') %>%
  FindClusters(resolution = 0.2)

seu_test$NMF_clusters <- seu_test$originalexp_snn_res.0.2

predictions <- SingleR::SingleR(
  test = as.SingleCellExperiment(seu_test),
  ref = as.SingleCellExperiment(seu_ref),
  labels = seu_ref$Group,
  de.n = 50,
  clusters = seu_test$NMF_clusters
)

seu_test$NMF_clusters_SingleR <- plyr::mapvalues(
  seu_test$NMF_clusters,
  from = rownames(predictions),
  to = predictions$labels)

print('pCMF')

data_matrix <- t(as.matrix(GetAssayData(seu_test, slot = "counts")))

# Pre-filter genes
kept_cols <- prefilter(data_matrix, prop = 0.05, quant_max = 0.95,
                       presel = TRUE, threshold = 0.2)
data_matrix <- data_matrix[,kept_cols]

# Run pCMF
res_pcmf <- pCMF(data_matrix, K = 20, verbose = FALSE,
                 zero_inflation = TRUE, sparsity = TRUE,
                 ncores = 10)

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
  FindNeighbors(dims = 1:10, reduction = 'pCMF', verbose = F) %>%
  FindClusters(resolution = 0.2, verbose = F)

seu_test$pCMF_clusters <- seu_test$originalexp_snn_res.0.2

predictions <- SingleR::SingleR(
  test = as.SingleCellExperiment(seu_test),
  ref = as.SingleCellExperiment(seu_ref),
  labels = seu_ref$Group,
  de.n = 50,
  clusters = seu_test$pCMF_clusters
)

seu_test$pCMF_clusters_SingleR <- plyr::mapvalues(
  seu_test$pCMF_clusters,
  from = rownames(predictions),
  to = predictions$labels)

print('GLM-PCA')

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
  FindNeighbors(dims = 1:10, reduction = 'gbm', verbose = F) %>%
  FindClusters(resolution = 0.2, verbose = F)

seu_test$GBM_clusters <- seu_test$originalexp_snn_res.0.2

predictions <- SingleR::SingleR(
  test = as.SingleCellExperiment(seu_test),
  ref = as.SingleCellExperiment(seu_ref),
  labels = seu_ref$Group,
  de.n = 50,
  clusters = seu_test$GBM_clusters
)

seu_test$GBM_clusters_SingleR <- plyr::mapvalues(
  seu_test$GBM_clusters,
  from = rownames(predictions),
  to = predictions$labels)

print('CELLMENTOR')

seu_test$CellMentor <- CreateDimReducObject(
  embeddings = t(as.matrix(h_test)),
  key = paste0('CellMentor', "_"),
  assay = DefaultAssay(seu_test)
)
# print('READ OBJECT')
seu_test <- seu_test %>% 
  RunUMAP(reduction = 'CellMentor', dims= 1:K_VALUE, reduction.name = 'umap_cellmentor', verbose = F) %>% 
  FindNeighbors(dims = 1:K_VALUE, reduction = 'CellMentor', verbose = F) %>% 
  FindClusters(resolution = 0.2, verbose = F)

seu_test$CellMentor_clusters <- seu_test$originalexp_snn_res.0.2

predictions <- SingleR::SingleR(
  test = as.matrix(h_test),
  ref = as.matrix(final_model@H),
  labels = final_model@annotation$celltype,
  # de.n = 50,
  clusters = seu_test$CellMentor_clusters
)

seu_test$CellMentor_clusters_SingleR <- plyr::mapvalues(
  seu_test$CellMentor_clusters,
  from = rownames(predictions),
  to = predictions$labels)

seu_test2 <- seu_test
seu_test2$Group <- 'Unknown'
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
# seu_test_counts <- GetAssayData(seu_test, slot = "counts")
seu_ref@assays$RNA <- seu_ref$originalexp

adata <- convertFormat(seu_ref, from="seurat", to="anndata",
                       main_layer="counts", drop_single_values=FALSE)
adata$X = adata$X$tocsr()
scvi$model$SCVI$setup_anndata(adata, batch_key = NULL, labels_key = "Group")
model = scvi$model$SCVI(adata)

print('SCANVI TRAIN')
model$train()

scanvi_model = scvi$model$SCANVI$from_scvi_model(
  model,
  unlabeled_category="Unknown"
)

scanvi_model$train()

print('SCANVI DONE')
seu_test2@assays$RNA <- seu_test2$originalexp
adata_test <- convertFormat(seu_test2, from="seurat", to="anndata",
                            main_layer="counts", drop_single_values=FALSE)
adata_test$X = adata_test$X$tocsr()
latent_scanvi = scanvi_model$get_latent_representation(adata_test)

latent_scanvi <- as.matrix(latent_scanvi)
rownames(latent_scanvi) = colnames(seu_test)
seu_test[["scanvi"]] <- CreateDimReducObject(embeddings = latent_scanvi, key = "scanvi_", assay = DefaultAssay(seu_test))

seu_test <- seu_test %>%
  RunUMAP(reduction = 'scanvi', dims= 1:10, reduction.name = 'umap_scanvi', verbose = F) %>%
  FindNeighbors(dims = 1:10, reduction = 'scanvi', verbose = F) %>%
  FindClusters(resolution = 0.2, verbose = F)

seu_test$scanvi_clusters <- seu_test$originalexp_snn_res.0.2

predictions <- SingleR::SingleR(
  test = as.SingleCellExperiment(seu_test),
  ref = as.SingleCellExperiment(seu_ref),
  labels = seu_ref$Group,
  clusters = seu_test$scanvi_clusters
)

seu_test$scanvi_clusters_SingleR <- plyr::mapvalues(
  seu_test$scanvi_clusters,
  from = rownames(predictions),
  to = predictions$labels)

qsave(seu_test, paste0('save_data/simulation_6tools', '_newcelltype', SIM_NUMBER, '_scanvi.qs'))
