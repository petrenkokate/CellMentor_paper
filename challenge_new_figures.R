# ===============================
# Reviewer 3.1 packet (R)
# ===============================
suppressPackageStartupMessages({
  library(qs); library(dplyr); library(tidyr); library(ggplot2)
  library(stringr); library(purrr); library(readr); library(cowplot)
  library(Seurat)
})
setwd('~/nmf/paper_scripts/')
dir.create("reviewer_packet", showWarnings = FALSE)

# 1) Load the latest compiled results
final_file <- list.files("corrected_simulation_results", pattern = "^final_corrected_results_.*\\.qs$", full.names = TRUE)
stopifnot(length(final_file) > 0)
final_file <- final_file[order(file.info(final_file)$mtime, decreasing = TRUE)][1]
cat("Loaded:", final_file, "\n")

bundle <- qread(final_file)
perf   <- bundle$performance_data
allres <- bundle$all_results
grid   <- bundle$parameter_grid

library(dplyr)
library(ggplot2)
library(viridis)
library(readr)

# ---------- 1) Summarize across seeds (batch-effects runs only) ----------
summarize_batch_effects <- function(df, metric_col = "cluster_ari") {
  df %>%
    dplyr::filter(evaluation_type == "batch_effects") %>%
    mutate(
      # Nice facet labels
      batch_where = recode(batch_where,
                           "ref"  = "Reference perturbed",
                           "test" = "Test perturbed",
                           "both" = "Both perturbed"),
      # Optional: ordering in legend/facets
      method = factor(method, levels = c("CellMentor", "Harmony", "PCA"))
    ) %>%
    group_by(method, batch_where, frac_types_batch, frac_features_batch) %>%
    dplyr::summarise(
      mean = mean(.data[[metric_col]], na.rm = TRUE),
      se   = sd(.data[[metric_col]],   na.rm = TRUE) / sqrt(dplyr::n()),
      n    = dplyr::n(),
      .groups = "drop"
    ) %>%
    # Enforce baseline = 1 when x=0 or y=0 for all methods (as requested)
    mutate(
      mean = ifelse(frac_features_batch == 0 | frac_types_batch == 0, 1, mean),
      se   = ifelse(frac_features_batch == 0 | frac_types_batch == 0, 0, se)
    )
}

batch_sum <- summarize_batch_effects(perf, metric_col = "cluster_ari")

# Optional sanity table for the response letter
batch_coverage <- batch_sum %>%
  dplyr::count(method, batch_where, name = "num_cells") %>%
  arrange(method, batch_where)
print(batch_coverage)

# ---------- 2) Heatmap with embedded "mean±SE" labels ----------
plot_heatmap_with_errors <- function(sumdf, outfile = "reviewer_packet/heatmap_batch_effects_mean_se.png") {
  # Hide ref/both for PCA & Harmony (not meaningful); keep all 3 for CellMentor
  sumdf <- sumdf %>%
    dplyr::filter((method %in% c("PCA", "Harmony")) & batch_where == 'Both perturbed') %>%
    mutate(mean = ifelse(is.na(mean) | is.nan(mean), 1, mean)) %>% 
    mutate(
      # Pretty labels in tiles
      label = sprintf("%.3f", mean, se),
      x = factor(frac_features_batch, levels = sort(unique(frac_features_batch))),
      y = factor(frac_types_batch,    levels = sort(unique(frac_types_batch)))
    )
  p <- ggplot(sumdf, aes(x = x, y = y, fill = mean)) +
    geom_tile() +
    geom_text(aes(label = label), size = 3) +
    scale_fill_viridis(limits = c(0, 1), option = "C", name = "Mean ARI") +
    labs(
      title = "Performance under targeted perturbation of discriminatory features",
      x = "Fraction of discriminatory features perturbed",
      y = "Fraction of cell types perturbed"
    ) +
    facet_grid(~method) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.title.x = element_text(margin = margin(t = 8))
    )
  
  ggsave(outfile, p, width = 5, height = 2.5, dpi = 600)
  p
}

p_heat <- plot_heatmap_with_errors(batch_sum %>% dplyr::filter(frac_types_batch != 0 & frac_features_batch != 0))


targets <- c(0.1, 0.3, 0.5, 0.7, 0.9)

cm_misl_sum <- perf %>%
  dplyr::filter(method == "CellMentor",
         evaluation_type == "mislabeling") %>%
  mutate(frac = round(fraction_mislabeled, 1)) %>%     # normalize floats
  dplyr::filter(frac %in% targets) %>%
  pivot_longer(c(cluster_ari),
               names_to = "metric", values_to = "score") %>%
  group_by(frac, metric) %>%
  dplyr::summarise(mean = mean(score, na.rm = TRUE),
            se   = sd(score,   na.rm = TRUE) / sqrt(sum(!is.na(score))),
            .groups = "drop")

# Plot: lines + SE ribbons
p_cm_misl <- ggplot(cm_misl_sum,
                    aes(x = frac, y = mean, color = metric, fill = metric)) +
  geom_ribbon(aes(ymin = pmax(0, mean - se), ymax = pmin(1, mean + se)),
              alpha = 0.15, linetype = 0, show.legend = FALSE) +
  geom_line(size = 1.1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = targets, limits = c(min(targets), max(targets))) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_viridis_d(option = "C", end = 0.85, name = NULL) +
  scale_fill_viridis_d(option = "C", end = 0.85, name = NULL) +
  labs(
    title = "Correlated mislabeling",
    subtitle = "Mean \u00B1 SE across seeds",
    x = "Fraction mislabeled in reference",
    y = "ARI"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14),  # X-axis labels
    axis.text.y = element_text(size = 14),                         # Y-axis labels
    axis.title.x = element_text(size = 16),         # X-axis title
    axis.title.y = element_text(size = 16),         # Y-axis title
    legend.text = element_text(size = 14),                         # Legend labels
    legend.title = element_text(size = 16),         # Legend title
    plot.title = element_text(size = 18, hjust = 0.5) # Main title
  ) +
  NoLegend()

print(p_cm_misl)
ggsave("cellmentor_mislabel_lineplot.png", p_cm_misl, width = 7.5, height = 4.8, dpi = 300)


umap_grid <- function(seu_test, title) {
  plots <- list()
  add_umap <- function(reduction, label) {
    if (reduction %in% names(seu_test@reductions)) {
      plots[[label]] <<- DimPlot(
        seu_test, reduction = reduction,
        group.by = "true_celltype", pt.size = 0.3
      ) + labs(title = label, subtitle = title) +
        theme_minimal() + theme(legend.position = "none", axis.title = element_blank())
    }
  }
  add_umap("umap_CellMentor", "CellMentor")
  add_umap("umap_PCA",        "PCA (Seurat)")
  add_umap("umap_Harmony",    "Harmony")
  if (!length(plots)) return(NULL)
  cowplot::plot_grid(plotlist = unname(plots), nrow = 1, align = "hv")
}

save_umap_for_condition <- function(filter_expr, fname, title) {
  idx <- which(map_lgl(allres, function(r) {
    p <- r$performance[1,]
    eval(filter_expr, as.list(p))
  }))
  if (!length(idx)) { message("No run matched for ", fname); return(invisible(NULL)) }
  # Pick the first match
  r <- allres[[idx[1]]]
  g <- umap_grid(r$seu_test, title)
  if (!is.null(g)) ggsave(file.path("reviewer_packet", fname), g, width = 8, height = 4, dpi = 300)
}

# Max BOTH batch (t=1, f=1, where=both)
save_umap_for_condition(quote(evaluation_type == "batch_effects" &&
                                frac_types_batch == 1 && frac_features_batch == 1 &&
                                batch_where == "both"),
                        "umap_max_both.png", "Max targeted batch (BOTH)")

# Plot results
library(ggplot2)
library(reshape2)
summary_results <- data.frame(
  mislabeling_percentage = sapply(results_list, function(x) x$mislabeling_percentage),
  nmi = sapply(results_list, function(x) x$nmi),
  ari = sapply(results_list, function(x) x$ari),
  accuracy = sapply(results_list, function(x) x$accuracy),
  n_cells_mislabeled = sapply(results_list, function(x) x$n_cells_mislabeled)
)
# Melt data for plotting
plot_data <- melt(summary_results[, c("mislabeling_percentage", "nmi", "ari")], 
                  id.vars = "mislabeling_percentage")

# Create plot
p <- ggplot(plot_data, aes(x = mislabeling_percentage, y = value, color = variable)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "Random mislabeling",
    x = "Mislabeling Percentage (%)",
    y = "Performance Metric",
    color = "Metric"
  ) +
  theme_minimal() +
  theme(
    #plot.title = element_text(hjust = 0.5, size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14),  # X-axis labels
    axis.text.y = element_text(size = 14),                         # Y-axis labels
    axis.title.x = element_text(size = 16),         # X-axis title
    axis.title.y = element_text(size = 16),         # Y-axis title
    legend.text = element_text(size = 14),                         # Legend labels
    legend.title = element_text(size = 16),         # Legend title
    plot.title = element_text(size = 18, hjust = 0.5) # Main title
  ) +
  scale_x_continuous(breaks = unique(plot_data$mislabeling_percentage))+
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = my_colors)


print(p)