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

simulation_params <- list(
  # Simulation 6
  list(
    nGenes = 10000, 
    seed = 100, 
    batchCells = c(1000, 1000),
    batch.facLoc = 0.3,
    batch.facScale = 0.5,
    out.prob = 0.1, 
    group.prob = c(0.18, 0.23, 0.16,0.03, 0.10, 0.3),
    de.prob = c(0.3, 0.1, 0.4, 0.5, 0.2, 0.1),
    de.downProb = c(0.1, 0.4, 0.2, 0.2, 0.1, 0.6),
    de.facLoc = c(0.6, 0.01, 0.5, 0.45, 0.4, 0.01),
    de.facScale = c(0.2, 0.5, 0.3, 0.35, 0.3, 0.6),
    bcv.common = 0.8,
    bcv.df = 5,
    dropout.mid = c(3, 5),
    dropout.shape = c(-2.5, -2.5),
    dropout.type = 'batch',
    lib.loc = 15, 
    lib.scale = 1.0
  )
)

# Function to run a single simulation
run_simulation <- function(sim_index, params) {
  SIM_NUMBER = paste0('_sensitivity_', sim_index)
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
                      bcv.df = 5,
                      dropout.mid = sim_params$dropout.mid,
                      dropout.shape = sim_params$dropout.shape,
                      dropout.type = 'batch',
                      lib.loc = 15,
                      lib.scale = 1.0
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
    ScaleData() %>%
    RunPCA(verbose=FALSE) %>%
    RunUMAP(dims = 1:10) %>%
    FindNeighbors(dims = 1:10) %>%
    FindClusters(resolution = 0.2)
  
  seu_test$Seurat_clusters <- seu_test$originalexp_snn_res.0.2
  
  # Run CellMentor analysis
  print("CreateCSFNMFobject")
  object = CreateCSFNMFobject(seu_ref@assays$originalexp$counts, 
                              seu_ref$Group, 
                              seu_test@assays$originalexp$counts)
  
  library(qs)
  
  # Define file path
  file_path <- paste0('save_data_v2/simulation_method', SIM_NUMBER, '.qs')
  
  # Check if file exists
  if (file.exists(file_path)) {
    print("BEST METHOD ALREADY EXISTS. LOADING FILE...")
    optimal_params_meth <- qread(file_path)
  } else {
    print("SELECT BEST METHOD")
    optimal_params_meth <- CellMentor(object, 
                                                     num_cores = 10, 
                                                     init_methods = c("uniform", "regulated", "NNDSVD", "skmeanGenes", "skmeanCells"),
                                                    alpha_range = c(0.1, 0.5, 1, 2, 5, 10, 20),
                                                    beta_range = c(0.1, 0.5, 1, 2, 5, 10, 20),
                                                     gamma_range = c(0, 0.001, 0.01, 0.1, 0.5, 1, 2),, 
                                                     delta_range = c(0, 0.1, 0.5, 1, 2, 5, 10))
    
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
  
  print('CELLMENTOR')
  
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
  
  # Save final results
  qsave(seu_test, 
        paste0('save_data/simulation_6tools', SIM_NUMBER, '.qs'))
  
  qsave(seu_ref, 
        paste0('save_data/simulation_6tools_ref', SIM_NUMBER, '.qs'))
  
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

print("All simulations completed!")