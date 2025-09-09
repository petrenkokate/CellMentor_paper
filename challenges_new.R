# Corrected CellMentor Simulation Script
# Addresses reviewer concerns with coherent batch effects and correlated mislabeling

# Load required libraries
library(Matrix)
library(Seurat)
library(SingleR)
library(CellMentor)
library(qs)
library(dplyr)
library(parallel)
library(aricode)
library(ggplot2)
library(harmony)
library(plyr)
library(digest)

# Set global parameters
set.seed(100)
N_REPETITIONS <- 3
BASE_SEEDS <- c(100, 200, 300)
setwd('~/nmf/paper_scripts/')
# ============================================================================
# PART 1: IMPROVED SIMULATION FUNCTIONS
# ============================================================================

#' Generate base expression data with more realistic parameters
generate_base_expression <- function(n_genes, n_cells, cell_type_id, 
                                     base_expression = 0.5, noise_level = 0.8) {
  
  # Create base expression matrix with lower baseline
  base_matrix <- matrix(
    rpois(n_genes * n_cells, lambda = base_expression),
    nrow = n_genes, ncol = n_cells
  )
  
  # Add cell-type specific marker genes (top 15% of genes)
  n_markers <- ceiling(n_genes * 0.15)
  marker_indices <- sample(n_genes, n_markers)
  
  # More modest enhancement to avoid perfect separation
  enhancement_factor <- runif(n_markers, 1.5, 3.0)
  base_matrix[marker_indices, ] <- base_matrix[marker_indices, ] * enhancement_factor
  
  # Add substantial noise
  noise_matrix <- matrix(
    rnorm(n_genes * n_cells, 0, noise_level),
    nrow = n_genes, ncol = n_cells
  )
  
  # Ensure non-negative counts
  final_matrix <- pmax(base_matrix + noise_matrix, 0)
  
  # Set row and column names
  rownames(final_matrix) <- paste0("Gene_", 1:n_genes)
  colnames(final_matrix) <- paste0(cell_type_id, "_Cell_", 1:n_cells)
  
  return(as(final_matrix, "CsparseMatrix"))
}

#' Identify discriminatory genes for each cell type
find_discriminatory_genes <- function(expression_matrix, cell_labels, top_n = 150) {
  
  unique_types <- unique(cell_labels)
  discriminatory_genes <- list()
  
  for (cell_type in unique_types) {
    type_cells <- which(cell_labels == cell_type)
    other_cells <- which(cell_labels != cell_type)
    
    if (length(type_cells) < 3 || length(other_cells) < 3) next
    
    # Calculate mean expression in this type vs others
    type_mean <- Matrix::rowMeans(expression_matrix[, type_cells, drop = FALSE])
    other_mean <- Matrix::rowMeans(expression_matrix[, other_cells, drop = FALSE])
    
    # Calculate fold change and significance
    fold_change <- log2((type_mean + 0.1) / (other_mean + 0.1))
    
    # Select top discriminatory genes based on absolute fold change
    top_genes_idx <- order(abs(fold_change), decreasing = TRUE)[1:min(top_n, length(fold_change))]
    top_genes <- rownames(expression_matrix)[top_genes_idx]
    
    discriminatory_genes[[cell_type]] <- top_genes
  }
  
  return(discriminatory_genes)
}

apply_coherent_batch_effects <- function(expression_matrix, cell_labels, batch_labels,
                                         fraction_celltypes_affected = 0.5,
                                         fraction_features_affected = 0.3,
                                         batch_strength = 3.0) {
  
  discriminatory_genes <- find_discriminatory_genes(expression_matrix, cell_labels)
  unique_types <- unique(cell_labels)
  affected_types <- sample(unique_types, ceiling(length(unique_types) * fraction_celltypes_affected))
  
  modified_matrix <- expression_matrix
  
  for (cell_type in affected_types) {
    if (!cell_type %in% names(discriminatory_genes)) next
    
    type_genes <- discriminatory_genes[[cell_type]]
    n_affected_genes <- ceiling(length(type_genes) * fraction_features_affected)
    affected_genes <- sample(type_genes, n_affected_genes)
    
    # Get test batch cells of this type
    type_cells <- which(cell_labels == cell_type)
    test_cells <- intersect(type_cells, which(batch_labels %in% c("Batch3", "Batch4")))
    
    if (length(test_cells) > 0) {
      
      # CORRUPT discriminatory features with batch noise
      # This makes cell type discrimination harder
      for (gene in affected_genes) {
        gene_idx <- which(rownames(modified_matrix) == gene)
        if (length(gene_idx) > 0) {
          
          # Add systematic batch noise that corrupts discrimination
          batch_noise <- rnorm(length(test_cells), mean = 0, sd = batch_strength)
          modified_matrix[gene_idx, test_cells] <- 
            pmax(0, modified_matrix[gene_idx, test_cells] + batch_noise)
        }
      }
    }
  }
  
  return(modified_matrix)
}

# --- REPLACE the old apply_coherent_batch_effects() with this ---
apply_targeted_batch_shifts <- function(expression_matrix, cell_labels, batch_labels,
                                        fraction_celltypes_affected = 0.5,
                                        fraction_features_affected = 0.3,
                                        batch_logfc = log2(3),              # ~3x median shift
                                        mode = c("opposite", "one_sided"), # "opposite": B3 up, B4 down (or vice versa)
                                        batches_test = c("Batch3","Batch4"),
                                        jitter_sdlog = 0.25,                # per-cell jitter around the shift
                                        de_top_n = 150) {
  
  mode <- match.arg(mode)
  
  # 1) identify discriminatory genes per cell type (like you already do)
  discr <- find_discriminatory_genes(expression_matrix, cell_labels, top_n = de_top_n)
  
  unique_types <- unique(cell_labels)
  n_affect_types <- ceiling(length(unique_types) * fraction_celltypes_affected)
  affected_types <- sample(unique_types, n_affect_types)
  
  # Bookkeeping for sanity checks
  perturbed <- list()
  
  # 2) only perturb TEST batches
  test_idx <- which(batch_labels %in% batches_test)
  
  X <- expression_matrix # Csparse ok; we'll assign into it
  
  for (ct in affected_types) {
    if (!ct %in% names(discr)) next
    ct_genes <- discr[[ct]]
    if (length(ct_genes) == 0) next
    
    n_affect_genes <- ceiling(length(ct_genes) * fraction_features_affected)
    genes_hit <- sample(ct_genes, n_affect_genes)
    
    ct_cells_all <- which(cell_labels == ct)
    ct_cells_test <- intersect(ct_cells_all, test_idx)
    if (length(ct_cells_test) == 0) next
    
    # split by test batches for opposite/one_sided behavior
    b3 <- intersect(ct_cells_test, which(batch_labels == batches_test[1]))
    b4 <- intersect(ct_cells_test, which(batch_labels == batches_test[2]))
    
    # Define multiplicative factors per batch
    # Use log-normal jitter around the target logFC so it's not perfectly coherent.
    if (mode == "opposite") {
      # Batch3 up, Batch4 down (or the other way around randomly per cell type)
      flip <- sample(c(-1, 1), 1)
      lfc_b3 <-  flip * batch_logfc
      lfc_b4 <- -flip * batch_logfc
    } else { # one_sided
      lfc_b3 <- batch_logfc
      lfc_b4 <- 0
    }
    
    fac_b3 <- if (length(b3)) rlnorm(length(b3), meanlog = lfc_b3 * log(2), sdlog = jitter_sdlog) else numeric(0)
    fac_b4 <- if (length(b4)) rlnorm(length(b4), meanlog = lfc_b4 * log(2), sdlog = jitter_sdlog) else numeric(0)
    
    # 3) multiplicatively perturb ONLY the discriminatory genes for this cell type
    #    This creates a real centroid shift that PCA can capture.
    for (g in genes_hit) {
      gi <- which(rownames(X) == g)
      if (!length(gi)) next
      
      if (length(b3)) X[gi, b3] <- X[gi, b3] * fac_b3
      if (length(b4)) X[gi, b4] <- X[gi, b4] * fac_b4
      # no pmax(0, ...) needed; multiplicative keeps non-negativity
    }
    
    perturbed[[ct]] <- list(genes = genes_hit,
                            b3_cells = b3, b4_cells = b4,
                            lfc_b3 = lfc_b3, lfc_b4 = lfc_b4)
  }
  
  attr(X, "perturbed") <- perturbed
  return(X)
}
# Fast centroids once, reused
calculate_cell_type_distances <- function(expression_matrix, cell_labels) {
  types <- unique(cell_labels)
  # turn to dense once (7000x~600 is fine in RAM)
  X <- as.matrix(expression_matrix)
  # centroids as columns
  centroids <- vapply(types, function(tp) {
    cols <- cell_labels == tp
    rowMeans(X[, cols, drop = FALSE])
  }, numeric(nrow(X)))
  colnames(centroids) <- types
  
  # squared distances via matrix ops
  # ||x - c||^2 = ||x||^2 + ||c||^2 - 2 x·c
  x_sq <- colSums(X * X)
  c_sq <- colSums(centroids * centroids)
  dot  <- t(X) %*% centroids              # cells x types
  D2   <- sweep(sweep(-2 * dot, 2, -c_sq, `+`), 1, x_sq, `+`)
  
  # for each cell, own vs closest other
  own_idx  <- match(cell_labels, types)
  own_d2   <- D2[cbind(seq_len(nrow(D2)), own_idx)]
  D2_oth   <- D2
  D2_oth[cbind(seq_len(nrow(D2_oth)), own_idx)] <- Inf
  min_oth_d2 <- apply(D2_oth, 1, min)
  min_oth_id <- apply(D2_oth, 1, which.min)
  
  distances <- lapply(seq_len(ncol(X)), function(i) list(
    cell_idx = i,
    current_type = cell_labels[i],
    own_distance = sqrt(own_d2[i]),
    other_distances = setNames(sqrt(D2[i, types != cell_labels[i]]), types[types != cell_labels[i]]),
    closest_other_type = types[min_oth_id[i]],
    closest_other_distance = sqrt(min_oth_d2[i])
  ))
  attr(distances, "centroids") <- centroids
  return(distances)
}

create_correlated_mislabeling <- function(expression_matrix, cell_labels, fraction_mislabeled = 0.15,
                                          blend_strength = 0.3, max_markers = 50) {
  # dense once
  X <- as.matrix(expression_matrix)
  
  # distances + reuse centroids
  cell_distances <- calculate_cell_type_distances(X, cell_labels)
  centroids <- attr(cell_distances, "centroids")
  types <- colnames(centroids)
  rn <- rownames(X)
  
  # boundary cells (vectorized)
  own <- sapply(cell_distances, `[[`, "own_distance")
  clo <- sapply(cell_distances, `[[`, "closest_other_distance")
  boundary_cells <- which(own / clo > (1/1.3))
  
  if (!length(boundary_cells)) {
    # fallback: pick most ambiguous cells by ratio
    ratios <- sapply(cell_distances, function(z) z$own_distance / z$closest_other_distance)
    n_mis <- ceiling(length(cell_labels) * fraction_mislabeled)
    boundary_cells <- order(ratios, decreasing = TRUE)[seq_len(min(n_mis, length(ratios)))]
  }
  
  # precompute discriminatory genes once
  discr <- find_discriminatory_genes(X, cell_labels)
  
  n_mis <- ceiling(length(boundary_cells) * fraction_mislabeled)
  cells_to_mislabel <- sample(boundary_cells, n_mis)
  
  mislabeled_labels <- cell_labels
  for (ci in cells_to_mislabel) {
    tgt <- cell_distances[[ci]]$closest_other_type
    if (!length(discr[[tgt]])) next
    
    # indices of target markers (vectorized)
    genes <- discr[[tgt]][seq_len(min(max_markers, length(discr[[tgt]])))]
    gi <- match(genes, rn)
    gi <- gi[!is.na(gi)]
    if (!length(gi)) next
    
    # blend in one shot
    X[gi, ci] <- X[gi, ci] * (1 - blend_strength) + centroids[gi, tgt] * blend_strength
    mislabeled_labels[ci] <- tgt
  }
  
  return(list(
    mislabeled_labels = mislabeled_labels,
    modified_expression = as(Matrix::Matrix(X, sparse = TRUE), "dgCMatrix")
  ))
}

#' Identify boundary cells that are candidates for realistic mislabeling
identify_boundary_cells <- function(cell_distances, ambiguity_threshold = 1.3) {
  
  boundary_cells <- c()
  
  for (i in 1:length(cell_distances)) {
    cell_info <- cell_distances[[i]]
    
    # Cell is "boundary" if it's closer to another cell type than expected
    distance_ratio <- cell_info$own_distance / cell_info$closest_other_distance
    
    if (distance_ratio > (1 / ambiguity_threshold)) {
      boundary_cells <- c(boundary_cells, i)
    }
  }
  
  return(boundary_cells)
}



# ============================================================================
# PART 2: METHOD EVALUATION FUNCTIONS
# ============================================================================

#' Run CellMentor with UMAP generation
run_cellmentor <- function(seu_ref, seu_test) {
  
  print("Running CellMentor...")
  
  tryCatch({
    # Create CSFNMF object
    object <- CreateCSFNMFobject(
      seu_ref@assays$RNA$counts, 
      seu_ref$celltype, 
      seu_test@assays$RNA$counts,
      verbose = FALSE
    )
    
    # Run CellMentor with streamlined parameters
    optimal_params <- CellMentor(
      object, 
      num_cores = 4,
      init_methods = c('regulated'),
      alpha_range = c(1, 5),
      beta_range = c(1, 5),
      gamma_range = c(0.1),
      delta_range = c(1),
      verbose = FALSE
    )
    
    final_model <- optimal_params$best_model
    K_VALUE <- optimal_params$best_params$k
    
    # Project test data
    h_test <- project_data(
      W = final_model@W,
      X = final_model@matrices@data,
      seed = 1,
      num_cores = 4,
      verbose = FALSE
    )
    
    # Add to Seurat object with UMAP
    seu_test[["CellMentor"]] <- CreateDimReducObject(
      embeddings = t(as.matrix(h_test)),
      key = "CellMentor_",
      assay = DefaultAssay(seu_test)
    )
    
    # Run UMAP, neighbors, and clustering
    seu_test <- seu_test %>%
      RunUMAP(reduction = 'CellMentor', dims = 1:K_VALUE, 
              reduction.name = 'umap_CellMentor', verbose = FALSE) %>%
      FindNeighbors(dims = 1:K_VALUE, reduction = 'CellMentor', 
                    graph.name = "CellMentor_snn", verbose = FALSE) %>%
      FindClusters(resolution = 0.2, graph.name = "CellMentor_snn", verbose = FALSE)
    
    seu_test$CellMentor_clusters <- seu_test$CellMentor_snn_res.0.2
    
    # SingleR annotation
    predictions <- SingleR::SingleR(
      test = as.matrix(h_test),
      ref = as.matrix(final_model@H),
      labels = final_model@annotation$celltype,
      clusters = seu_test$CellMentor_clusters,
      de.n = 50
    )
    
    seu_test$CellMentor_SingleR <- plyr::mapvalues(
      seu_test$CellMentor_clusters,
      from = rownames(predictions),
      to = predictions$labels,
      warn_missing = FALSE
    )
    
  }, error = function(e) {
    print(paste("CellMentor failed:", e))
    seu_test$CellMentor_clusters <- rep("1", ncol(seu_test))
    seu_test$CellMentor_SingleR <- rep("Unknown", ncol(seu_test))
  })
  
  return(seu_test)
}

#' Run PCA with UMAP generation
run_pca <- function(seu_ref, seu_test) {
  
  print("Running PCA...")
  
  tryCatch({
    seu_test <- seu_test %>%
      NormalizeData(verbose = FALSE) %>%
      FindVariableFeatures(verbose = FALSE) %>%
      ScaleData(verbose = FALSE) %>%
      RunPCA(verbose = FALSE) %>%
      RunUMAP(dims = 1:20, reduction = 'pca', 
              reduction.name = 'umap_PCA', verbose = FALSE) %>%
      FindNeighbors(dims = 1:20, reduction = 'pca', 
                    graph.name = "PCA_snn", verbose = FALSE) %>%
      FindClusters(resolution = 0.2, graph.name = "PCA_snn", verbose = FALSE)
    
    seu_test$PCA_clusters <- seu_test$PCA_snn_res.0.2
    
    # SingleR annotation
    predictions <- SingleR::SingleR(
      test = as.SingleCellExperiment(seu_test),
      ref = as.SingleCellExperiment(seu_ref),
      labels = seu_ref$celltype,
      clusters = seu_test$PCA_clusters,
      de.n = 50
    )
    
    seu_test$PCA_SingleR <- plyr::mapvalues(
      seu_test$PCA_clusters,
      from = rownames(predictions),
      to = predictions$labels,
      warn_missing = FALSE
    )
    
  }, error = function(e) {
    print(paste("PCA failed:", e))
    seu_test$PCA_clusters <- rep("1", ncol(seu_test))
    seu_test$PCA_SingleR <- rep("Unknown", ncol(seu_test))
  })
  
  return(seu_test)
}

#' Run Harmony with UMAP generation
run_harmony <- function(seu_ref, seu_test) {
  
  print("Running Harmony...")
  
  tryCatch({
    # Temporarily combine datasets for harmony
    seu_ref$dataset <- "ref"
    seu_test$dataset <- "test"
    seu_combined <- merge(seu_ref, seu_test)
    
    seu_combined <- seu_combined %>%
      NormalizeData(verbose = FALSE) %>%
      FindVariableFeatures(verbose = FALSE) %>%
      ScaleData(verbose = FALSE) %>%
      RunPCA(verbose = FALSE)
    
    # Run Harmony
    seu_combined <- RunHarmony(seu_combined, group.by.vars = "dataset", verbose = FALSE)
    
    # Extract test results
    test_cells <- colnames(seu_test)
    harmony_embeddings <- seu_combined@reductions$harmony@cell.embeddings[test_cells, ]
    
    seu_test[["harmony"]] <- CreateDimReducObject(
      embeddings = harmony_embeddings,
      key = "Harmony_",
      assay = DefaultAssay(seu_test)
    )
    
    seu_test <- seu_test %>%
      RunUMAP(dims = 1:20, reduction = 'harmony', 
              reduction.name = 'umap_Harmony', verbose = FALSE) %>%
      FindNeighbors(dims = 1:20, reduction = 'harmony', 
                    graph.name = "Harmony_snn", verbose = FALSE) %>%
      FindClusters(resolution = 0.2, graph.name = "Harmony_snn", verbose = FALSE)
    
    seu_test$Harmony_clusters <- seu_test$Harmony_snn_res.0.2
    
    # SingleR annotation
    predictions <- SingleR::SingleR(
      test = as.SingleCellExperiment(seu_test),
      ref = as.SingleCellExperiment(seu_ref),
      labels = seu_ref$celltype,
      clusters = seu_test$Harmony_clusters,
      de.n = 50
    )
    
    seu_test$Harmony_SingleR <- plyr::mapvalues(
      seu_test$Harmony_clusters,
      from = rownames(predictions),
      to = predictions$labels,
      warn_missing = FALSE
    )
    
  }, error = function(e) {
    print(paste("Harmony failed:", e))
    seu_test$Harmony_clusters <- rep("1", ncol(seu_test))
    seu_test$Harmony_SingleR <- rep("Unknown", ncol(seu_test))
  })
  
  return(seu_test)
}

#' Calculate performance metrics with proper column checking
calculate_performance_metrics <- function(seu_test, true_labels) {
  
  methods <- c("CellMentor", "PCA", "Harmony")
  results <- data.frame(
    method = methods,
    cluster_ari = NA,
    singler_ari = NA,
    cluster_nmi = NA,
    singler_nmi = NA
  )
  
  for (i in seq_along(methods)) {
    method <- methods[i]
    
    # Check clustering performance
    cluster_col <- paste0(method, "_clusters")
    if (cluster_col %in% colnames(seu_test@meta.data)) {
      pred_clusters <- seu_test@meta.data[[cluster_col]]
      if (length(unique(pred_clusters)) > 1) {
        results$cluster_ari[i] <- aricode::ARI(true_labels, pred_clusters)
        results$cluster_nmi[i] <- aricode::NMI(true_labels, pred_clusters)
      }
    }
    
    # Check SingleR performance
    singler_col <- paste0(method, "_SingleR")
    if (singler_col %in% colnames(seu_test@meta.data)) {
      pred_singler <- seu_test@meta.data[[singler_col]]
      if (length(unique(pred_singler)) > 1) {
        results$singler_ari[i] <- aricode::ARI(true_labels, pred_singler)
        results$singler_nmi[i] <- aricode::NMI(true_labels, pred_singler)
      }
    }
  }
  
  return(results)
}

# ============================================================================
# PART 3: PARAMETER GRID
# ============================================================================

create_parameter_grid <- function() {
  
  # Batch effect evaluation (5x5 grid as reviewer requested)
  batch_grid <- expand.grid(
    frac_types_batch = seq(0, 1, 0.25),
    frac_features_batch = seq(0, 1, 0.25),
    batch_strength = 3.0,
    fraction_mislabeled = 0.0,
    evaluation_type = "batch_effects"
  )
  
  # Mislabeling evaluation (3 levels)
  mislabel_grid <- expand.grid(
    frac_types_batch = 0.0,
    frac_features_batch = 0.0,
    batch_strength = 1.0,
    fraction_mislabeled = seq(0, 1, 0.25),
    evaluation_type = "mislabeling"
  )
  
  return(rbind(batch_grid, mislabel_grid))
}

generate_complete_simulation <- function(
    n_genes = 5000, n_cells_per_type = 200, n_cell_types = 6,
    fraction_celltypes_affected_batch = 0.5,
    fraction_features_affected_batch  = 0.3,
    batch_strength = 3.0, fraction_mislabeled = 0.15,
    seed = 100, batch_where = "test",
    mode_ref = "opposite", mode_test = "opposite"
) {
  set.seed(seed_salt(
    seed,
    n_genes, n_cells_per_type, n_cell_types,
    fraction_celltypes_affected_batch,
    fraction_features_affected_batch,
    batch_strength, fraction_mislabeled,
    batch_where, mode_ref, mode_test
  ))
  
  base_types <- paste0("CellType_", 1:n_cell_types)
  all_mats <- list(); all_true <- c(); all_batch <- c()
  
  for (ct in base_types) {
    M <- generate_base_expression(n_genes, n_cells_per_type, ct)
    all_mats[[ct]] <- M
    all_true <- c(all_true, rep(ct, ncol(M)))
    all_batch <- c(all_batch, sample(rep(paste0("Batch", 1:4), length.out = ncol(M))))
  }
  
  X <- do.call(cbind, all_mats)
  
  if (fraction_celltypes_affected_batch > 0 || fraction_features_affected_batch > 0) {
    if (batch_where %in% c("ref","both")) {
      message("Applying targeted batch shifts to REFERENCE (Batch1/Batch2)...")
      X <- apply_targeted_batch_shifts(
        X, cell_labels = all_true, batch_labels = all_batch,
        fraction_celltypes_affected = fraction_celltypes_affected_batch,
        fraction_features_affected  = fraction_features_affected_batch,
        batch_logfc = log2(batch_strength),
        mode = mode_ref, batches_test = c("Batch1","Batch2")
      )
    }
    if (batch_where %in% c("test","both")) {
      message("Applying targeted batch shifts to TEST (Batch3/Batch4)...")
      X <- apply_targeted_batch_shifts(
        X, cell_labels = all_true, batch_labels = all_batch,
        fraction_celltypes_affected = fraction_celltypes_affected_batch,
        fraction_features_affected  = fraction_features_affected_batch,
        batch_logfc = log2(batch_strength),
        mode = mode_test, batches_test = c("Batch3","Batch4")
      )
    }
  }
  
  ref_idx <- which(all_batch %in% c("Batch1","Batch2"))
  ref_labels <- all_true[ref_idx]
  
  if (fraction_mislabeled > 0) {
    message("Creating correlated mislabeling in REFERENCE...")
    mis <- create_correlated_mislabeling(X[, ref_idx], ref_labels, fraction_mislabeled)
    ref_labels_mislabeled <- mis$mislabeled_labels
    X[, ref_idx] <- mis$modified_expression
  } else {
    ref_labels_mislabeled <- ref_labels
  }
  
  list(
    expression_matrix = X,
    true_labels = all_true,
    reference_labels_mislabeled = ref_labels_mislabeled,
    batch_labels = all_batch,
    ref_cell_indices = ref_idx
  )
}

# ============================================================================
# PART 4: MAIN SIMULATION FUNCTION
# ============================================================================

run_single_simulation <- function(param_set, seed = 100) {
  
  set.seed(seed)
  
  print(paste(
    "Running simulation with:",
    "batch_types =",     param_set$frac_types_batch,
    "batch_features =",  param_set$frac_features_batch,
    "batch_strength =",  param_set$batch_strength,
    "where =",           param_set$batch_where,
    "mode_ref =",        param_set$mode_ref,
    "mode_test =",       param_set$mode_test,
    "mislabeled =",      param_set$fraction_mislabeled,
    "type =",            param_set$evaluation_type
  ))
  
  # Generate simulation data (ALL from param_set now)
  sim_data <- generate_complete_simulation(
    n_genes                          = 7000,
    n_cells_per_type                 = 200,
    n_cell_types                     = 6,
    fraction_celltypes_affected_batch= param_set$frac_types_batch,
    fraction_features_affected_batch = param_set$frac_features_batch,
    batch_strength                   = param_set$batch_strength,
    fraction_mislabeled              = param_set$fraction_mislabeled,
    seed                             = seed,
    batch_where                      = param_set$batch_where,  # << from grid
    mode_ref                         = param_set$mode_ref,     # << from grid
    mode_test                        = param_set$mode_test     # << from grid
  )
  
  
  # Create Seurat objects
  ref_cells <- sim_data$ref_cell_indices
  test_cells <- setdiff(1:ncol(sim_data$expression_matrix), ref_cells)
  
  # Reference dataset
  seu_ref <- CreateSeuratObject(
    counts = sim_data$expression_matrix[, ref_cells]
  )
  seu_ref$celltype <- sim_data$reference_labels_mislabeled
  seu_ref$true_celltype <- sim_data$true_labels[ref_cells]
  
  # Test dataset
  seu_test <- CreateSeuratObject(
    counts = sim_data$expression_matrix[, test_cells]
  )
  seu_test$true_celltype <- sim_data$true_labels[test_cells]
  
  # Prepare reference data
  seu_ref <- seu_ref %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
  
  # Run methods
  seu_test <- run_cellmentor(seu_ref, seu_test)
  seu_test <- run_pca(seu_ref, seu_test)
  seu_test <- run_harmony(seu_ref, seu_test)
  
  # Calculate performance
  performance <- calculate_performance_metrics(seu_test, seu_test$true_celltype)
  
  # Add parameter information
  performance$frac_types_batch <- param_set$frac_types_batch
  performance$frac_features_batch <- param_set$frac_features_batch
  performance$batch_strength <- param_set$batch_strength
  performance$fraction_mislabeled <- param_set$fraction_mislabeled
  performance$evaluation_type <- param_set$evaluation_type
  performance$seed <- seed
  performance$batch_where         <- param_set$batch_where
  performance$mode_ref            <- param_set$mode_ref
  performance$mode_test           <- param_set$mode_test
  
  return(list(
    performance = performance,
    seu_test = seu_test,
    seu_ref = seu_ref,
    sim_data = sim_data
  ))
}

# ============================================================================
# PART 5: MAIN EXECUTION
# ============================================================================

seed_salt <- function(seed, ...) {
  h <- digest(list(...), algo = "xxhash32")
  salt <- strtoi(substr(h, 1, 7), 16L)  # 7 hex ~ 2^28 < .Machine$integer.max
  as.integer((seed + salt) %% .Machine$integer.max)
}

# Create parameter grid
create_parameter_grid <- function() {
  batch_grid <- expand.grid(
    frac_types_batch   = seq(0, 1, 0.25),
    frac_features_batch= seq(0, 1, 0.25),
    batch_strength     = c(2),        # ~2x, 4x, 8x multiplicative on affected markers
    fraction_mislabeled= 0.0,
    batch_where        = c("ref","test","both"),
    mode_ref           = "opposite",        # reviewer’s confounding case
    mode_test          = "opposite",
    evaluation_type    = "batch_effects",
    stringsAsFactors   = FALSE
  )
  
  mislabel_grid <- expand.grid(
    frac_types_batch   = 0.0,
    frac_features_batch= 0.0,
    batch_strength     = 1.0,
    fraction_mislabeled= seq(0, 1, 0.1),
    batch_where        = "none",            # no targeted batch in mislabel runs
    mode_ref           = "opposite",
    mode_test          = "opposite",
    evaluation_type    = "mislabeling",
    stringsAsFactors   = FALSE
  )
  
  rbind(batch_grid, mislabel_grid)
}
param_grid <- create_parameter_grid()
print(paste("Total parameter combinations:", nrow(param_grid)))
param_grid <- subset(param_grid,
                     !(evaluation_type == "batch_effects" &
                         ( (frac_types_batch == 0 & frac_features_batch > 0) |
                             (frac_types_batch > 0 & frac_features_batch == 0) ))
)

# Initialize results storage
all_results <- list()
result_counter <- 1

# Create output directory
dir.create("corrected_simulation_results", showWarnings = FALSE)

# Run simulations
for (i in 1:nrow(param_grid)) {
  
  print(paste("=== Parameter combination", i, "of", nrow(param_grid), "==="))
  param_set <- param_grid[i, ]
  
  # Run multiple seeds for this parameter combination
  for (seed_idx in 1:length(BASE_SEEDS)) {
    
    current_seed <- BASE_SEEDS[seed_idx]
    print(paste("Running seed", seed_idx, "of", length(BASE_SEEDS)))
    
    # Check if results already exist
    result_file <- file.path("corrected_simulation_results", 
                             paste0("result_", i, "_seed_", current_seed, ".qs"))
    
    result <- tryCatch({
      run_single_simulation(param_set, seed = current_seed)
    }, error = function(e) {
      print(paste("Simulation failed:", e))
      list(
        performance = data.frame(
          method = c("CellMentor", "PCA", "Harmony"),
          cluster_ari = NA, singler_ari = NA,
          cluster_nmi = NA, singler_nmi = NA,
          frac_types_batch = param_set$frac_types_batch,
          frac_features_batch = param_set$frac_features_batch,
          batch_strength = param_set$batch_strength,
          fraction_mislabeled = param_set$fraction_mislabeled,
          evaluation_type = param_set$evaluation_type,
          seed = current_seed
        ),
        error = as.character(e)
      )
    })
    qsave(result, result_file)
    
    all_results[[result_counter]] <- result
    result_counter <- result_counter + 1
  }
}

# ============================================================================
# PART 6: RESULTS COMPILATION
# ============================================================================

print("Compiling results...")

# Extract performance data
performance_list <- lapply(all_results, function(x) x$performance)
performance_list <- performance_list[!sapply(performance_list, is.null)]
final_results <- do.call(rbind, performance_list)

# Save final results
final_file <- file.path("corrected_simulation_results", 
                        paste0("final_corrected_results_", 
                               format(Sys.time(), "%Y%m%d_%H%M"), ".qs"))

qsave(list(
  performance_data = final_results,
  all_results = all_results,
  parameter_grid = param_grid
), final_file)