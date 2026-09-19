###Beta diversity now::
# Load Packages:
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)
library(compositions)
library(zCompositions)
library(ggpubr)

# Load ps object
setwd("C:/Users/mateu/OneDrive/Area_de_Trabalho/2025-12_Gustavo_UFV")
ps_orig <- readRDS("./ASV_tables/For_use/ps_final.rds")

## Rizosphere and Fruits are too different, we must plot as separate subsets.
ps_melip <- subset_samples(ps_orig, Genus %in% "Melipona")
ps_melip <- prune_taxa(taxa_sums(ps_melip) > 0, ps_melip)
ps_melip <- filter_taxa(ps_melip, function(x) sum(x) > 2, TRUE)
ps_nomelip <- subset_samples(ps_orig, !(Genus %in% "Melipona"))
ps_nomelip <- prune_taxa(taxa_sums(ps_nomelip) > 0, ps_nomelip)
ps_nomelip <- filter_taxa(ps_nomelip, function(x) sum(x) > 2, TRUE)

## Set parameters:
ps <- ps_orig                 # Change based on which dataset you want
g_var <- "Local_Name"        # Change to desired variable in your metadata
g_var_e <- "Local_Name"      # Change if you want the shapes as a different variable
out_dir <- "Beta_diversity/Local_name" # Output directory

# Colors for plotting -> optional but recommended
colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#78e",
  "#e31a1c", "#fdbf6f", "#f00", "#cab2d6", "#6a3d9a",
  "#ffff99", "#b15928", "#ffd92f", "#b3b3b3", "#ccebc5",
  "#decbe4", "#fed9a6", "#ffffcc", "#e5d8bd", "#32e", "yellow", "#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231",
  "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe",
  "#008080", "#e6beff", "#aa6e28", "#800000", "#aaffc3",
  "#808000", "#ffd8b1", "#000080", "#808080", "#032bec",
  "#000000", "#7cfc00", "#40e0d0", "#6495ed", "#ff4500",
  "#32cd32", "#8a2be2", "#ff1493", "#00ff7f", "#dc143c",
  "#ffd700", "#ff69b4", "#00ced1", "#9932cc", "#ff6347"
  
)
shapes <- c(15, 16, 17, 18, 0, 1, 2, 5, 7, 10, 9, 8, 6, 3)
#------- Functions----
# Helper Functions
get_metadata_df <- function(ps) {
  meta <- data.frame(sample_data(ps), check.names = FALSE) %>%
    rownames_to_column("Sample") %>%
    mutate(
      LibrarySize = sample_sums(ps),
      Log_LibrarySize = log10(LibrarySize)
    )
  return(meta)
}
# Aitchison
compute_robust_aitchison <- function(ps) {
  # 1. Extract OTU matrix
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)
  # 2. Zero replacement 
  # Set z.delete = FALSE to prevent dropping samples/taxa
  otu_zr <- zCompositions::cmultRepl(
    otu, 
    method = "CZM",        # More robust for sparse data
    output = "p-counts", 
    z.delete = FALSE,      # CRITICAL: Stops the function from deleting data
    z.warning = FALSE
  )
  # 3. CLR transform and Euclidean distance
  otu_clr <- compositions::clr(otu_zr)
  dist(otu_clr, method = "euclidean")}
# Score Extraction
extract_ordination_scores <- function(ord, method) {
  if (method == "PCoA") {
    # Calculate percentages for PCoA
    pc_var <- round(100 * (ord$eig / sum(ord$eig)), 1)
    df <- as.data.frame(ord$points[, 1:2])
    attr(df, "axis_labels") <- paste0("Axis ", 1:2, " [", pc_var[1:2], "%]")
  } else {
    # Extract NMDS scores and Stress
    df <- as.data.frame(vegan::scores(ord, display = "sites"))
    attr(df, "stress") <- round(ord$stress, 3)
  }
  colnames(df) <- c("Dim1", "Dim2")
  return(df %>% rownames_to_column("Sample"))
}
#Plotting
plot_ordination <- function(ord, ps, Variable, ShapeVar, method, distance_name, outdir) {
  meta <- get_metadata_df(ps)
  scores_df_raw <- extract_ordination_scores(ord, method)
  
  scores_df <- scores_df_raw %>%
    left_join(meta, by = "Sample")
  
  # Ensure factors
  scores_df[[Variable]] <- factor(scores_df[[Variable]])
  scores_df[[ShapeVar]] <- factor(scores_df[[ShapeVar]])
  
  # Set Axis Labels
  x_lab <- "Dimension 1"
  y_lab <- "Dimension 2"
  subtitle_text <- paste("Color:", Variable, "| Shape:", ShapeVar)
  
  if (method == "PCoA") {
    labs_attr <- attr(scores_df_raw, "axis_labels")
    x_lab <- labs_attr[1]
    y_lab <- labs_attr[2]
  } else if (method == "NMDS") {
    stress_val <- attr(scores_df_raw, "stress")
    subtitle_text <- paste(subtitle_text, "| Stress:", stress_val)
  }
  
  p <- ggplot(scores_df, aes(x = Dim1, y = Dim2)) +
    geom_point(aes(color = !!sym(Variable), shape = !!sym(ShapeVar)), size = 4, alpha = 0.8) +
    scale_color_manual(values = colors) + # Uses your custom colors
    scale_shape_manual(values = shapes) +
    stat_ellipse(aes(color = !!sym(Variable)), linetype = 2, size = 0.5) + ### Toggle this for elipses
    geom_point(aes(color = !!sym(Variable), shape = !!sym(ShapeVar)), size = 4, alpha = 0.8) + ### Toggle this for elipses
    theme_bw(base_size = 14) +
    labs(
      title = paste(method, "-", distance_name),
      subtitle = subtitle_text,
      x = x_lab,
      y = y_lab,
      color = Variable,
      shape = ShapeVar
    )
  
  # 1. Extract the legend
  leg <- get_legend(p)
  p_legend <- as_ggplot(leg)
  
  # 2. Remove legend from the main plot
  p_no_legend <- p + theme(legend.position = "none")
  
  # 3. Save the Plot
  file_base <- file.path(outdir, paste0(distance_name, "_", Variable, "_", method))
  ggsave(paste0(file_base, "_plot.pdf"), p_no_legend, width = 8, height = 6)
  ggsave(paste0(file_base, ".png"), p_no_legend, width = 8, height = 6, dpi = 300)
  
  # 4. Save the Legend (adjust width/height based on legend size)
  ggsave(paste0(file_base, "_legend.pdf"), p_legend, width = 4, height = 8)
  
  return(p_no_legend)
}
# Wrapper Function
run_beta_pipeline <- function(ps, Variable, ShapeVar, outdir) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  metrics <- c("bray", "wunifrac", "aitchison")
  results_list <- list()
  
  for (m in metrics) {
    message("Processing Distance: ", m)
    
    # Compute Distance
    dist_mat <- if (m == "bray") phyloseq::distance(ps, method = "bray") else if (m == "wunifrac") phyloseq::distance(ps, method = "wunifrac") else compute_robust_aitchison(ps)
    
    # Run PERMANOVA
    meta <- get_metadata_df(ps)
    ids <- attr(dist_mat, "Labels")
    meta_sub <- meta %>% filter(Sample %in% ids) %>% arrange(match(Sample, ids))
    
    # Formula uses both variables
    form <- as.formula(paste("dist_mat ~", Variable, "+", ShapeVar, "+ Log_LibrarySize"))
    perm <- vegan::adonis2(form, data = meta_sub, permutations = 999, by = "margin")
    write.csv(as.data.frame(perm), file.path(outdir, paste0(m, "_", Variable, "_PERMANOVA.csv")))
    
    # Ordinations
    pcoa <- cmdscale(dist_mat, k = 2, eig = TRUE)
    nmds <- vegan::metaMDS(dist_mat, k = 2, trymax = 50, autotransform = FALSE, trace = FALSE)
    
    # Plot & Save
    plot_ordination(pcoa, ps, Variable, ShapeVar, "PCoA", m, outdir)
    plot_ordination(nmds, ps, Variable, ShapeVar, "NMDS", m, outdir)
    
    results_list[[m]] <- list(dist = dist_mat, permanova = perm)
  }
  return(results_list)
}
#------- Call function-------
results <- run_beta_pipeline(ps = ps, Variable = g_var, ShapeVar = g_var_e, outdir = out_dir)