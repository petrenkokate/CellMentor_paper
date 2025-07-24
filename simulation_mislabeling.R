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
library(NNLM)
library(pCMF)
library(glmpca)
library(mclust)
library(cluster)
library(aricode)
library(clevr)
library(SingleR)
library(harmony)

SIM_NUMBER = '_mislabeling_challenge'
set.seed(100)

# Mislabeling percentages to test
mislabeling_percentages <- c(0, 5, 10, 20, 30, 40, 50)

print("START SIMULATION")

# Create simulation parameters
params <- newSplatParams()
params <- setParams(params,
                    nGenes = 20000, 
                    seed = 100, 
                    batchCells = c(1000, 1000, 1000, 1000),
                    batch.facLoc = 0.1,
                    batch.facScale = 0.1,
                    out.prob = 0.1, 
                    group.prob = c(0.18, 0.23, 0.12, 0.07, 0.10, 0.3),
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

# Create reference and test datasets
print("Create ref and test datasets")
seu$Group <- as.character(seu$Group)
seu$celltype <- seu$Group
seu_ref <- subset(seu, Batch %in% c('Batch1', 'Batch2'))
seu_test <- subset(seu, Batch %in% c('Batch3', 'Batch4'))

# Function to introduce mislabeling
introduce_mislabeling <- function(labels, mislabel_percent) {
  if (mislabel_percent == 0) {
    return(list(mislabeled = labels, annotation = labels))  # Fix this line!
  }
  
  n_cells <- length(labels)
  n_mislabel <- round(n_cells * mislabel_percent / 100)
  all_labels <- unique(labels)
  # Randomly select which cells to mislabel
  cells_to_mislabel <- sample(seq_len(n_cells), n_mislabel)
  mislabeled <- labels
  mislabeled[cells_to_mislabel] <- "mislabeled group"
  
  mislabeled_groups <- labels
  for (i in cells_to_mislabel) {
    current_label <- labels[i]
    # Get all other labels (excluding current one)
    other_labels <- setdiff(all_labels, current_label)
    # Randomly pick one existing label
    mislabeled_groups[i] <- sample(other_labels, 1)
  }
  
  cat(paste("Changed", n_mislabel, "out of", n_cells, "cells to 'mislabeled group' (", mislabel_percent, "%)\n"))
  
  return(list(mislabeled = mislabeled,
              annotation = mislabeled_groups))
}

results_list <- list()

for (mislabel_pct in mislabeling_percentages) {
  print(paste("Testing mislabeling percentage:", mislabel_pct, "%"))
  
  # Create mislabeled reference
  seu_ref_mislabeled <- seu_ref
  new_labels <- introduce_mislabeling(
    seu_ref_mislabeled$Group, 
    mislabel_pct)
  seu_ref_mislabeled$mislabeled_cells <- new_labels$mislabeled
  seu_ref_mislabeled$Group_mislabeled <- new_labels$annotation
  
  seu_ref_mislabeled <- seu_ref_mislabeled %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA(verbose=FALSE) %>% 
    RunUMAP(dims = 1:10) %>% 
    FindNeighbors(dims = 1:10) %>% 
    FindClusters(resolution = 0.2)
  
  seu_test_processed <- seu_test %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(verbose=FALSE) %>%
    RunUMAP(dims = 1:10) %>%
    FindNeighbors(dims = 1:10) %>%
    FindClusters(resolution = 0.2)
  
  print("CreateCSFNMFobject")
  # Create CSFNMF object using mislabeled reference
  object <- CreateCSFNMFobject(
    seu_ref_mislabeled@assays$originalexp$counts, 
    seu_ref_mislabeled$Group_mislabeled,  # Use mislabeled annotations
    seu_test_processed@assays$originalexp$counts
  )
  
  print("SELECT BEST METHOD")
  # Run parameter optimization (simplified for speed)
  optimal_params <- CellMentor(
    object, 
    alpha_range = c(1, 5),
    beta_range = c(1, 5), 
    gamma_range = c(1), 
    delta_range = c(0.1),
    num_cores = 10
  )
  
  K_VALUE <- optimal_params$best_params$k
  final_model <- optimal_params$best_model
  
  print("PREDICTION PROJECTION")
  h_test <- project_data(
    W = final_model@W,
    X = final_model@matrices@data,
    seed = 1,
    num_cores = 10,
    chunk_size = NULL,
    verbose = TRUE
  )
  
  print("CELLMENTOR CLUSTERING")
  seu_test_processed$CellMentor <- CreateDimReducObject(
    embeddings = t(as.matrix(h_test)),
    key = "CellMentor_",
    assay = DefaultAssay(seu_test_processed)
  )
  
  seu_test_processed <- seu_test_processed %>% 
    RunUMAP(reduction = 'CellMentor', dims = 1:K_VALUE, 
            reduction.name = 'umap_cellmentor', verbose = F) %>% 
    FindNeighbors(dims = 1:K_VALUE, reduction = 'CellMentor', verbose = F) %>% 
    FindClusters(resolution = 0.2, verbose = F)
  
  seu_test_processed$CellMentor_clusters <- seu_test_processed$originalexp_snn_res.0.2
  
  seu_ref_mislabeled$CellMentor <- CreateDimReducObject(
    embeddings = t(as.matrix(optimal_params$best_model@H)),
    key = paste0('CellMentor', "_"),
    assay = DefaultAssay(seu_ref_mislabeled),
    loadings = as.matrix(final_model@W)
  )
  
  seu_ref_mislabeled <- seu_ref_mislabeled %>%
    RunUMAP(reduction = 'CellMentor', dims= 1:K_VALUE, reduction.name = 'umap_cellmentor', verbose = F)
  
  # Predict cell types using SingleR
  # predictions <- SingleR::SingleR(
  #   test = as.matrix(h_test),
  #   ref = as.matrix(final_model@H),
  #   labels = final_model@annotation$celltype,
  #   clusters = seu_test_processed$CellMentor_clusters
  # )
  # 
  # seu_test_processed$CellMentor_predicted <- plyr::mapvalues(
  #   seu_test_processed$CellMentor_clusters,
  #   from = rownames(predictions),
  #   to = predictions$labels
  # )
  
  # Calculate performance metrics
  true_labels <- seu_test_processed$Group  # True labels
  predicted_labels <- seu_test_processed$CellMentor_clusters
  
  # Calculate NMI (Normalized Mutual Information)
  nmi_score <- aricode::NMI(true_labels, predicted_labels)
  
  # Calculate ARI (Adjusted Rand Index)
  ari_score <- aricode::ARI(true_labels, predicted_labels)
  
  # Calculate accuracy (if labels match exactly)
  accuracy <- mean(true_labels == predicted_labels)
  
  # Store results
  results_list[[paste0("mislabel_", mislabel_pct)]] <- list(
    mislabeling_percentage = mislabel_pct,
    nmi = nmi_score,
    ari = ari_score,
    accuracy = accuracy,
    n_cells_mislabeled = sum(seu_ref_mislabeled$Group != seu_ref_mislabeled$Group_mislabeled),
    seu_test = seu_test_processed,
    seu_ref = seu_ref_mislabeled,
    optimal_params = optimal_params
  )
  
  print(paste("Mislabeling", mislabel_pct, "% - NMI:", round(nmi_score, 3), 
              "ARI:", round(ari_score, 3), "Accuracy:", round(accuracy, 3)))
}

# Save results
print("SAVING RESULTS")
qsave(results_list, paste0('save_data_v3/mislabeling_challenge_results', SIM_NUMBER, '.qs'))
print("SAVED")
# Create summary results dataframe
summary_results <- data.frame(
  mislabeling_percentage = sapply(results_list, function(x) x$mislabeling_percentage),
  nmi = sapply(results_list, function(x) x$nmi),
  ari = sapply(results_list, function(x) x$ari),
  accuracy = sapply(results_list, function(x) x$accuracy),
  n_cells_mislabeled = sapply(results_list, function(x) x$n_cells_mislabeled)
)

print("SUMMARY RESULTS:")
print(summary_results)


# Plot results
library(ggplot2)
library(reshape2)

# Melt data for plotting
plot_data <- melt(summary_results[, c("mislabeling_percentage", "nmi", "ari", "accuracy")], 
                  id.vars = "mislabeling_percentage")

# Create plot
p <- ggplot(plot_data, aes(x = mislabeling_percentage, y = value, color = variable)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "CellMentor Performance vs Reference Annotation Quality",
    x = "Mislabeling Percentage (%)",
    y = "Performance Metric",
    color = "Metric"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.position = "bottom"
  ) +
  scale_x_continuous(breaks = mislabeling_percentages) +
  scale_y_continuous(limits = c(0, 1))

print(p)

# Save plot
ggsave(paste0('save_data_v3/mislabeling_challenge_plot', SIM_NUMBER, '.pdf'), 
       plot = p, width = 10, height = 6)

print("MISLABELING CHALLENGE COMPLETE!")