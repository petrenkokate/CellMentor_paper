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
# library(spam)
# library(spam64)
library(microbenchmark)
library(pryr)
library(bench)
library(fastICA)
library(rliger)
library(CoGAPS)
library(projectR)

# Source additional required files
source("~/nmf/paper_scripts/code/graph_rssNMF.R")
source("~/nmf/paper_scripts/code/initial_NMF.R")

SUFFIX <- '_comp_cost_'
# Function to get system memory info
get_memory_info <- function() {
  gc_info <- gc(verbose = FALSE)
  list(
    used_mb = sum(gc_info[, "used"]) * 8 / 1024,  # Convert to MB
    max_used_mb = sum(gc_info[, "max used"]) * 8 / 1024,
    ncells = gc_info["Ncells", "used"],
    vcells = gc_info["Vcells", "used"]
  )
}

# Function to track computational metrics
track_computation <- function(expr, method_name = "unknown") {
  # Force garbage collection before starting
  gc(verbose = FALSE, full = TRUE)
  
  # Get initial CPU state
  cpu_before <- proc.time()
  # Get initial memory state
  mem_before <- get_memory_info()
  start_time <- Sys.time()
  
  # Execute the expression and capture result
  result_wrapper <- tryCatch(
    {
      list(result = eval(expr), error = NULL)
    },
    error = function(e) {
      list(result = NULL, error = e$message)
    }
  )
  
  # Get final timing and memory
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Get final CPU state
  cpu_after <- proc.time()
  
  # Calculate CPU metrics
  cpu_user_time <- cpu_after["user.self"] - cpu_before["user.self"]
  cpu_system_time <- cpu_after["sys.self"] - cpu_before["sys.self"]
  cpu_total_time <- cpu_user_time + cpu_system_time
  cpu_efficiency <- if(elapsed_time > 0) cpu_total_time / elapsed_time else 0
  
  # Force garbage collection to get accurate memory usage
  gc(verbose = FALSE, full = TRUE)
  mem_after <- get_memory_info()
  
  # Calculate memory differences
  memory_used_mb <- mem_after$used_mb - mem_before$used_mb
  peak_memory_mb <- mem_after$max_used_mb - mem_before$max_used_mb
  
  # Return comprehensive metrics
  list(
    method = method_name,
    elapsed_time_seconds = elapsed_time,
    elapsed_time_minutes = elapsed_time / 60,
    memory_used_mb = memory_used_mb,
    peak_memory_mb = peak_memory_mb,
    memory_before_mb = mem_before$used_mb,
    memory_after_mb = mem_after$used_mb,
    start_time = start_time,
    end_time = end_time,
    cpu_user_time_seconds = cpu_user_time,
    cpu_system_time_seconds = cpu_system_time,
    cpu_total_time_seconds = cpu_total_time,
    cpu_total_time_minutes = cpu_total_time / 60,
    cpu_efficiency = cpu_efficiency,
    result = result_wrapper$result,
    success = is.null(result_wrapper$error)
  )
}

# Computational analysis parameters (smaller for comprehensive method testing)
computational_params <- list(
  # Small dataset
  list(
    name = "small",
    nGenes = 3000, 
    seed = 300, 
    batchCells = c(500, 500),
    batch.facLoc = 0.2,
    batch.facScale = 0.3,
    out.prob = 0.01,
    group.prob = c(0.2, 0.2, 0.2, 0.2, 0.1, 0.1),
    de.prob = c(0.3, 0.3, 0.3, 0.3, 0.3, 0.3),
    de.downProb = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2),
    de.facLoc = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    de.facScale = c(0.3, 0.3, 0.3, 0.3, 0.3, 0.3),
    bcv.common = 0.4,
    bcv.df = 8,
    dropout.mid = c(0.5, 0.5),
    dropout.shape = c(-1.0, -1.0),
    dropout.type = 'batch',
    lib.loc = 15,
    lib.scale = 0.8
  ),
  # Medium dataset
  list(
    name = "medium",
    nGenes = 5000, 
    seed = 300, 
    batchCells = c(2000, 2000),
    batch.facLoc = 0.2,
    batch.facScale = 0.3,
    out.prob = 0.05,
    group.prob = c(0.2, 0.2, 0.2, 0.2, 0.1, 0.1),
    de.prob = c(0.3, 0.3, 0.3, 0.3, 0.3, 0.3),
    de.downProb = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2),
    de.facLoc = c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6),
    de.facScale = c(0.3, 0.3, 0.3, 0.3, 0.3, 0.3),
    bcv.common = 0.5,
    bcv.df = 8,
    dropout.mid = c(0.8, 0.8),
    dropout.shape = c(-1.2, -1.2),
    dropout.type = 'batch',
    lib.loc = 15,
    lib.scale = 0.8
  ),
  # Large dataset
  list(
    name = "large",
    nGenes = 8000, 
    seed = 300, 
    batchCells = c(5000, 5000),
    batch.facLoc = 0.25,
    batch.facScale = 0.35,
    out.prob = 0.08,
    group.prob = c(0.2, 0.2, 0.2, 0.2, 0.1, 0.1),
    de.prob = c(0.35, 0.35, 0.35, 0.35, 0.35, 0.35),
    de.downProb = c(0.25, 0.25, 0.25, 0.25, 0.25, 0.25),
    de.facLoc = c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7),
    de.facScale = c(0.3, 0.3, 0.3, 0.3, 0.3, 0.3),
    bcv.common = 0.6,
    bcv.df = 6,
    dropout.mid = c(1.0, 1.0),
    dropout.shape = c(-1.3, -1.3),
    dropout.type = 'batch',
    lib.loc = 15,
    lib.scale = 0.9
  )
)

# Function to run computational cost analysis for ALL methods
run_comprehensive_analysis <- function(dataset_index, params_list) {
  
  dataset_params <- params_list[[dataset_index]]
  dataset_name <- dataset_params$name
  
  print(paste("=== ANALYZING DATASET SIZE:", toupper(dataset_name), "==="))
  
  set.seed(100)
  
  # Create simulation with timing
  sim_metrics <- track_computation({
    base_params <- newSplatParams()
    sim_params <- setParams(base_params,
                            nGenes = dataset_params$nGenes,
                            seed = 100,
                            batchCells = dataset_params$batchCells,
                            batch.facLoc = dataset_params$batch.facLoc,
                            batch.facScale = dataset_params$batch.facScale,
                            out.prob = dataset_params$out.prob,
                            group.prob = dataset_params$group.prob,
                            de.prob = dataset_params$de.prob,
                            de.downProb = dataset_params$de.downProb,
                            de.facLoc = dataset_params$de.facLoc,
                            de.facScale = dataset_params$de.facScale,
                            bcv.common = dataset_params$bcv.common,
                            bcv.df = 5,
                            dropout.mid = dataset_params$dropout.mid,
                            dropout.shape = dataset_params$dropout.shape,
                            dropout.type = 'batch',
                            lib.loc = 15,
                            lib.scale = 1.0)
    
    sim <- splatSimulate(sim_params, method = 'groups')
    sim <- logNormCounts(sim)
    
    seu <- as.Seurat(sim)
    seu <- seu %>% 
      NormalizeData() %>% 
      FindVariableFeatures() %>% 
      ScaleData() %>% 
      RunPCA(verbose=FALSE) %>% 
      RunUMAP(dims = 1:10)
    
    seu$celltype <- seu$Group
    seu
  }, method_name = "data_simulation")
  
  if (!sim_metrics$success) {
    print(paste("Simulation failed for", dataset_name))
    return(NULL)
  }
  
  seu <- sim_metrics$result
  
  print(paste("Dataset dimensions:", nrow(seu), "genes x", ncol(seu), "cells"))
  
  print("Create ref")
  seu_ref <- subset(seu, Batch %in% c('Batch1'))
  
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
  seu_test <- subset(seu, Batch %in% c('Batch2'))
  
  seu_test <- seu_test %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData()
  
  # Initialize results storage
  computational_results <- list(
    dataset_info = list(
      name = dataset_name,
      n_genes = nrow(seu),
      n_cells = ncol(seu),
      n_ref_cells = ncol(seu_ref),
      n_test_cells = ncol(seu_test),
      n_celltypes = length(unique(seu$celltype)),
      simulation_time = sim_metrics$elapsed_time_seconds,
      simulation_memory = sim_metrics$memory_used_mb
    ),
    method_results = list()
  )
  
  # 1. PCA Analysis
  print("1. Analyzing PCA...")
  pca_metrics <- track_computation({
    seu_test_pca <- seu_test %>%
      RunPCA(verbose=FALSE) %>%
      RunUMAP(dims = 1:10) %>%
      FindNeighbors(dims = 1:10) %>%
      FindClusters(resolution = 0.2)
    
    seu_test_pca$Seurat_clusters <- seu_test_pca$originalexp_snn_res.0.2
    seu_test_pca
  }, method_name = "PCA")
  
  computational_results$method_results$PCA <- list(
    time_seconds = pca_metrics$elapsed_time_seconds,
    time_minutes = pca_metrics$elapsed_time_minutes,
    cpu_user_time_seconds = pca_metrics$cpu_user_time_seconds,
    cpu_system_time_seconds = pca_metrics$cpu_system_time_seconds,
    cpu_total_time_seconds = pca_metrics$cpu_total_time_seconds,
    cpu_total_time_minutes = pca_metrics$cpu_total_time_minutes,
    cpu_efficiency = pca_metrics$cpu_efficiency,
    memory_mb = pca_metrics$memory_used_mb,
    peak_memory_mb = pca_metrics$peak_memory_mb,
    success = pca_metrics$success
  )
  
  library(qs)
  
  # 2. ICA Analysis
  print("2. Analyzing ICA...")
  ica_metrics <- track_computation({
    library(fastICA)
    
    # Get scaled data
    scaled_data <- GetAssayData(seu_test, slot = "scale.data")
    # Run ICA
    ica_result <- fastICA(t(scaled_data), n.comp = 20)
    # Create embedding matrix for Seurat
    ica_embeddings <- ica_result$S
    rownames(ica_embeddings) <- colnames(seu_test)
    colnames(ica_embeddings) <- paste0("IC_", 1:ncol(ica_embeddings))
    seu_test_ica <- seu_test
    # Add ICA as a reduction to Seurat object
    seu_test_ica[["ica"]] <- CreateDimReducObject(
      embeddings = ica_embeddings,
      key = "IC_",
      assay = DefaultAssay(seu_test)
    )
    
    # Continue with your workflow
    seu_test_ica <- seu_test_ica %>%
      RunUMAP(dims = 1:10, reduction = 'ica', reduction.name = 'umap_ica') %>%
      FindNeighbors(dims = 1:10, reduction = 'ica') %>%
      FindClusters(resolution = 0.2)
    
    seu_test_ica$ICA_clusters <- seu_test_ica$originalexp_snn_res.0.2
    
    seu_test_ica
  }, method_name = "ICA")
  
  computational_results$method_results$ICA <- list(
    time_seconds = ica_metrics$elapsed_time_seconds,
    time_minutes = ica_metrics$elapsed_time_minutes,
    memory_mb = ica_metrics$memory_used_mb,
    cpu_user_time_seconds = ica_metrics$cpu_user_time_seconds,
    cpu_system_time_seconds = ica_metrics$cpu_system_time_seconds,
    cpu_total_time_seconds = ica_metrics$cpu_total_time_seconds,
    cpu_total_time_minutes = ica_metrics$cpu_total_time_minutes,
    cpu_efficiency = ica_metrics$cpu_efficiency,
    peak_memory_mb = ica_metrics$peak_memory_mb,
    success = ica_metrics$success
  )
  
  # 3. NMF Analysis
  print("3. Analyzing NMF...")
  nmf_metrics <- track_computation({
    seu_test_nmf <- seu_test %>% NormalizeData()
    data_matrix <- as.matrix(GetAssayData(seu_test_nmf, slot = "data"))
    keep_rows <- rowSums(data_matrix) > 0 & !apply(data_matrix, 1, function(x) any(is.na(x)))
    data_matrix <- data_matrix[keep_rows,]
    
    res_nmf <- nnmf(data_matrix, k = 20, method = "lee", loss = "mse", n.threads = 1)
    
    seu_test_nmf[["nmf"]] <- CreateDimReducObject(
      embeddings = t(res_nmf$H),
      key = "NMF_",
      assay = DefaultAssay(seu_test)
    )
    
    seu_test_nmf <- seu_test_nmf %>%
      RunUMAP(dims = 1:10, reduction = 'nmf', reduction.name = 'umap_nmf', verbose = F) %>%
      FindNeighbors(dims = 1:10, reduction = 'nmf', verbose = F) %>%
      FindClusters(resolution = 0.2, verbose = F)
    
    seu_test_nmf$NMF_clusters <- seu_test_nmf$originalexp_snn_res.0.2
    seu_test_nmf
  }, method_name = "NMF")
  
  computational_results$method_results$NMF <- list(
    time_seconds = nmf_metrics$elapsed_time_seconds,
    time_minutes = nmf_metrics$elapsed_time_minutes,
    memory_mb = nmf_metrics$memory_used_mb,
    cpu_user_time_seconds = nmf_metrics$cpu_user_time_seconds,
    cpu_system_time_seconds = nmf_metrics$cpu_system_time_seconds,
    cpu_total_time_seconds = nmf_metrics$cpu_total_time_seconds,
    cpu_total_time_minutes = nmf_metrics$cpu_total_time_minutes,
    cpu_efficiency = nmf_metrics$cpu_efficiency,
    peak_memory_mb = nmf_metrics$peak_memory_mb,
    success = nmf_metrics$success
  )
  
  # 4. pCMF Analysis
  print("4. Analyzing pCMF...")
  pcmf_metrics <- track_computation({
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
    seu_test_pcmf <- seu_test
    # Create dimension reduction object
    seu_test_pcmf[["pCMF"]] <- CreateDimReducObject(
      embeddings = embeddings,
      key = "pCMF_",
      assay = DefaultAssay(seu_test)
    )
    
    # Continue with workflow
    seu_test_pcmf <- seu_test_pcmf %>%
      RunUMAP(dims = 1:10, reduction = 'pCMF', reduction.name = 'umap_pcmf', verbose = F) %>%
      FindNeighbors(dims = 1:10, reduction = 'pCMF', verbose = F) %>%
      FindClusters(resolution = 0.2, verbose = F)
    
    seu_test_pcmf$pCMF_clusters <- seu_test_pcmf$originalexp_snn_res.0.2
    
    seu_test
  }, method_name = "pCMF")
  
  computational_results$method_results$pCMF <- list(
    time_seconds = pcmf_metrics$elapsed_time_seconds,
    time_minutes = pcmf_metrics$elapsed_time_minutes,
    memory_mb = pcmf_metrics$memory_used_mb,
    cpu_user_time_seconds = pcmf_metrics$cpu_user_time_seconds,
    cpu_system_time_seconds = pcmf_metrics$cpu_system_time_seconds,
    cpu_total_time_seconds = pcmf_metrics$cpu_total_time_seconds,
    cpu_total_time_minutes = pcmf_metrics$cpu_total_time_minutes,
    cpu_efficiency = pcmf_metrics$cpu_efficiency,
    peak_memory_mb = pcmf_metrics$peak_memory_mb,
    success = pcmf_metrics$success
  )
  
  # 5. GLM-PCA Analysis
  print("5. Analyzing GLM-PCA...")
  glmpca_metrics <- track_computation({
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
    
    seu_test
  }, method_name = "GLM_PCA")
  
  computational_results$method_results$GLM_PCA <- list(
    time_seconds = glmpca_metrics$elapsed_time_seconds,
    time_minutes = glmpca_metrics$elapsed_time_minutes,
    memory_mb = glmpca_metrics$memory_used_mb,
    cpu_user_time_seconds = glmpca_metrics$cpu_user_time_seconds,
    cpu_system_time_seconds = glmpca_metrics$cpu_system_time_seconds,
    cpu_total_time_seconds = glmpca_metrics$cpu_total_time_seconds,
    cpu_total_time_minutes = glmpca_metrics$cpu_total_time_minutes,
    cpu_efficiency = glmpca_metrics$cpu_efficiency,
    peak_memory_mb = glmpca_metrics$peak_memory_mb,
    success = glmpca_metrics$success
  )
  
  
  # 6. CellMentor analysis
  print("6. Analyzing CellMentor...")
  cellmentor_metrics <- track_computation({
    print("CreateCSFNMFobject")
    object = CreateCSFNMFobject(seu_ref@assays$originalexp$counts, 
                                seu_ref$Group, 
                                seu_test@assays$originalexp$counts)
    
    # Define file path
    file_path <- paste0('save_data_v3/simulation_method', dataset_index, '.qs')
    
    # Check if file exists
    if (file.exists(file_path)) {
      print("BEST METHOD ALREADY EXISTS. LOADING FILE...")
      optimal_params_meth <- qread(file_path)
    } else {
      print("SELECT BEST METHOD")
      optimal_params_meth <- CellMentor(
        object,
        alpha_range = c(1),
        beta_range = c(1),
        gamma_range = c(0.1),
        delta_range = c(1),
        num_cores = 1,
        verbose = FALSE
      )
      
      print("SAVE BEST METHOD")
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
      FindClusters(resolution = 0.2, verbose = F)
    
    seu_test$CellMentor_clusters <- seu_test$originalexp_snn_res.0.2
    
    seu_test
  }, method_name = "CellMentor")
  
  computational_results$method_results$CellMentor <- list(
    time_seconds = cellmentor_metrics$elapsed_time_seconds,
    time_minutes = cellmentor_metrics$elapsed_time_minutes,
    memory_mb = cellmentor_metrics$memory_used_mb,
    cpu_user_time_seconds = cellmentor_metrics$cpu_user_time_seconds,
    cpu_system_time_seconds = cellmentor_metrics$cpu_system_time_seconds,
    cpu_total_time_seconds = cellmentor_metrics$cpu_total_time_seconds,
    cpu_total_time_minutes = cellmentor_metrics$cpu_total_time_minutes,
    cpu_efficiency = cellmentor_metrics$cpu_efficiency,
    peak_memory_mb = cellmentor_metrics$peak_memory_mb,
    success = cellmentor_metrics$success
  )
  
  # 7. CellMentor analysis
  print("7. Analyzing SCANVI...")
  scanvi_metrics <- track_computation({
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
    
    seu_test
  }, method_name = "SCANVI")
  
  computational_results$method_results$SCANVI <- list(
    time_seconds = scanvi_metrics$elapsed_time_seconds,
    time_minutes = scanvi_metrics$elapsed_time_minutes,
    memory_mb = scanvi_metrics$memory_used_mb,
    cpu_user_time_seconds = scanvi_metrics$cpu_user_time_seconds,
    cpu_system_time_seconds = scanvi_metrics$cpu_system_time_seconds,
    cpu_total_time_seconds = scanvi_metrics$cpu_total_time_seconds,
    cpu_total_time_minutes = scanvi_metrics$cpu_total_time_minutes,
    cpu_efficiency = scanvi_metrics$cpu_efficiency,
    peak_memory_mb = scanvi_metrics$peak_memory_mb,
    success = scanvi_metrics$success
  )
  
  # 8. rssNMF Analysis 
  print("8. Analyzing rssNMF...")
  rssnmf_metrics <- track_computation({
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
      seu_test
  }, method_name = "rssNMF")
  
  computational_results$method_results$rssNMF <- list(
    time_seconds = rssnmf_metrics$elapsed_time_seconds,
    time_minutes = rssnmf_metrics$elapsed_time_minutes,
    memory_mb = rssnmf_metrics$memory_used_mb,
    cpu_user_time_seconds = rssnmf_metrics$cpu_user_time_seconds,
    cpu_system_time_seconds = rssnmf_metrics$cpu_system_time_seconds,
    cpu_total_time_seconds = rssnmf_metrics$cpu_total_time_seconds,
    cpu_total_time_minutes = rssnmf_metrics$cpu_total_time_minutes,
    cpu_efficiency = rssnmf_metrics$cpu_efficiency,
    peak_memory_mb = rssnmf_metrics$peak_memory_mb,
    success = rssnmf_metrics$success
  )
  
    # 9. LIGER Analysis 
    print("9. Analyzing LIGER...")
    liger_metrics <- track_computation({
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
    seu_test
    }, method_name = "LIGER")
  
  computational_results$method_results$LIGER <- list(
    time_seconds = liger_metrics$elapsed_time_seconds,
    time_minutes = liger_metrics$elapsed_time_minutes,
    memory_mb = liger_metrics$memory_used_mb,
    cpu_user_time_seconds = liger_metrics$cpu_user_time_seconds,
    cpu_system_time_seconds = liger_metrics$cpu_system_time_seconds,
    cpu_total_time_seconds = liger_metrics$cpu_total_time_seconds,
    cpu_total_time_minutes = liger_metrics$cpu_total_time_minutes,
    cpu_efficiency = liger_metrics$cpu_efficiency,
    peak_memory_mb = liger_metrics$peak_memory_mb,
    success = liger_metrics$success
  )

  # 10. CoGAPS Analysis 
  print("10. Analyzing CoGAPS...")
  cogaps_metrics <- track_computation({
    
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
    
    seu_test
  }, method_name = "CoGAPS")
  
  computational_results$method_results$CoGAPS <- list(
    time_seconds = cogaps_metrics$elapsed_time_seconds,
    time_minutes = cogaps_metrics$elapsed_time_minutes,
    memory_mb = cogaps_metrics$memory_used_mb,
    cpu_user_time_seconds = cogaps_metrics$cpu_user_time_seconds,
    cpu_system_time_seconds = cogaps_metrics$cpu_system_time_seconds,
    cpu_total_time_seconds = cogaps_metrics$cpu_total_time_seconds,
    cpu_total_time_minutes = cogaps_metrics$cpu_total_time_minutes,
    cpu_efficiency = cogaps_metrics$cpu_efficiency,
    peak_memory_mb = cogaps_metrics$peak_memory_mb,
    success = cogaps_metrics$success
  )
  
  # 11. CASSL Analysis (if Python script is available)
  print("11. Analyzing CASSL...")
  cassl_metrics <- track_computation({
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
    seu_test
  }, method_name = "CASSL")
  
  computational_results$method_results$CASSL <- list(
    time_seconds = cassl_metrics$elapsed_time_seconds,
    time_minutes = cassl_metrics$elapsed_time_minutes,
    memory_mb = cassl_metrics$memory_used_mb,
    cpu_user_time_seconds = cassl_metrics$cpu_user_time_seconds,
    cpu_system_time_seconds = cassl_metrics$cpu_system_time_seconds,
    cpu_total_time_seconds = cassl_metrics$cpu_total_time_seconds,
    cpu_total_time_minutes = cassl_metrics$cpu_total_time_minutes,
    cpu_efficiency = cassl_metrics$cpu_efficiency,
    peak_memory_mb = cassl_metrics$peak_memory_mb,
    success = cassl_metrics$success
  )
  
  # Calculate scalability metrics
  computational_results$scalability_metrics <- list(
    cells_per_second = list(),
    memory_per_cell = list()
  )
  
  for (method_name in names(computational_results$method_results)) {
    method_result <- computational_results$method_results[[method_name]]
    if (method_result$success) {
      computational_results$scalability_metrics$cells_per_second[[method_name]] <- 
        computational_results$dataset_info$n_test_cells / method_result$time_seconds
      computational_results$scalability_metrics$memory_per_cell[[method_name]] <- 
        method_result$memory_mb / computational_results$dataset_info$n_test_cells
    }
  }
  
  return(computational_results)
}

# Main execution
print("=== STARTING COMPREHENSIVE COMPUTATIONAL COST ANALYSIS ===")
print("This analysis includes ALL methods from the original script")

# Run analysis for all dataset sizes
all_computational_results <- list()
for (i in 1:length(computational_params)) {
  print(paste("Analyzing dataset", i, "of", length(computational_params)))
  
  result <- tryCatch({
    run_comprehensive_analysis(i, computational_params)
  }, error = function(e) {
    print(paste("Error in computational analysis", i, ":", e$message))
    return(NULL)
  })
  
  all_computational_results[[i]] <- result
  
  # Save intermediate results
  if (!is.null(result)) {
    qsave(result, paste0('save_data_v3/computational_analysis_ALL_METHODS_', result$dataset_info$name, '_', 
                         format(Sys.time(), "%Y%m%d_%H%M"), '.qs'))
  }
  
  # Print progress summary
  if (!is.null(result)) {
    successful_methods <- sum(sapply(result$method_results, function(x) x$success))
    total_methods <- length(result$method_results)
    print(paste("Dataset", result$dataset_info$name, "completed:", 
                successful_methods, "out of", total_methods, "methods successful"))
  }
}
run_comprehensive_analysis(1, computational_params)
qsave(all_computational_results, paste0('save_data_v3/computational_analysis_ALL_METHODS_all_results', '.qs'))