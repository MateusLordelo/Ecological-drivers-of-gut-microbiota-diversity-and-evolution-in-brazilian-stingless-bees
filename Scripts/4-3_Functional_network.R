library(phyloseq)
library(igraph)
library(Hmisc)
library(dplyr)
library(tibble)
library(compositions)

setwd("C:/Users/mateu/OneDrive/Area_de_Trabalho/shared_lis/29-spp")

# ------------ Load the files ----------
ps_KO <- readRDS("./ASV_tables/For_use/ps_KO.rds")
ps_Metacyc <- readRDS("./ASV_tables/For_use/ps_Metacyc.rds")

generate_functional_network_robust <- function(ps_obj, rank_name, output_prefix, top_n) {
  
  message(paste("--- Generating Robust Network for:", output_prefix, "---"))
  # 0. Filter "Unclassified" Taxonomy
  # We extract the tax table, find rows that AREN'T unclassified, and prune.
  tax_mat <- as(tax_table(ps_obj), "matrix")
  keep_taxa <- !(tax_mat[, rank_name] %in% c("Unclassified", "unclassified", "Unknown", "unknown", ""))
  # Update the ps_obj to only include these taxa
  ps_obj <- prune_taxa(keep_taxa, ps_obj)
  
  # 1. High Prevalence Filtering (Found in 30% of bees)
  ps_filt <- filter_taxa(ps_obj, function(x) sum(x > 0) > (0.30 * nsamples(ps_obj)), TRUE)
  
  # 2. Feature Selection: Keep top N most variable KOs/Pathways
  # Functions that don't change across samples just add noise.
  cv_values <- apply(otu_table(ps_filt), 1, function(x) sd(x) / mean(x))
  keep_top <- names(sort(cv_values, decreasing = TRUE))[1:min(top_n, length(cv_values))]
  ps_filt <- prune_taxa(keep_top, ps_filt)
  
  # 3. CLR Transformation (The "Biomass Bias" Killer)
  # We add a small pseudocount (1) then apply CLR
  otu <- as.matrix(otu_table(ps_filt))
  otu_clr <- t(apply(otu + 1, 2, function(x) log(x / exp(mean(log(x))))))
  
  tax <- as.data.frame(tax_table(ps_filt))
  
  # 4. Correlation (Spearman on CLR data)
  cor_res <- rcorr(otu_clr, type = "spearman")
  cor_rho <- cor_res$r
  cor_p <- cor_res$P
  
  # Diagnostic
  filt_neg <- sum(cor_rho < 0)
  filt_pos <- sum(cor_rho > 0)
  message(paste("Features analyzed:", ncol(otu_clr)))
  message(paste("Raw Negative Corrs:", sum(cor_rho < 0, na.rm = TRUE)))
  message(paste("Filtered Correlations (>0.6) - Positive:", filt_pos, "| Negative:", filt_neg))
  
  # 5. Strict Thresholding
  keep_mask <- (cor_p < 0.01 & abs(cor_rho) > 0.6)
  cor_rho[!keep_mask] <- 0
  diag(cor_rho) <- 0 
  
  # 6. Create igraph
  net <- graph_from_adjacency_matrix(cor_rho, mode = "undirected", weighted = TRUE)
  net <- delete_vertices(net, V(net)[degree(net) == 0])
  
  # 7. Edge List
  edges <- igraph::as_data_frame(net, what = "edges") %>%
    rename(Source = from, Target = to, Weight = weight) %>%
    mutate(Direction = ifelse(Weight > 0, "Positive", "Negative"))
  
  # 8. Node Metadata
  nodes <- data.frame(ID = V(net)$name) %>%
    left_join(rownames_to_column(tax, "ID"), by = "ID")
  
  # 9. Stats
  net_abs <- net
  E(net_abs)$weight <- abs(E(net_abs)$weight)
  clusters <- cluster_fast_greedy(net_abs)
  
  stats <- data.frame(
    Metric = c("Nodes", "Edges", "Positive_Edges", "Negative_Edges", 
               "Density", "Avg_Degree", "Modularity", "Transitivity"),
    Value = c(
      vcount(net), ecount(net),
      sum(edges$Direction == "Positive"), sum(edges$Direction == "Negative"),
      edge_density(net), mean(degree(net)),
      modularity(net_abs, membership(clusters)), transitivity(net, type = "global")
    )
  )
  
  # 10. Rank-Level Edge Counts
  rank_edge_counts <- edges %>%
    left_join(select(nodes, ID, !!sym(rank_name)), by = c("Source" = "ID")) %>%
    group_by(!!sym(rank_name)) %>%
    summarise(Edges = n(), .groups = "drop") %>%
    rename(Pathway = !!sym(rank_name))
  
  # 11. Export
  out_dir <- "./Functional_prediction/Functional_network/"
    # Create the directory if it doesn't exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  write.csv(edges, file.path(out_dir, paste0(output_prefix, "_edges.csv")), row.names = FALSE)
  write.csv(nodes, file.path(out_dir, paste0(output_prefix, "_nodes_metadata.csv")), row.names = FALSE)
  write.csv(stats, file.path(out_dir, paste0(output_prefix, "_network_stats.csv")), row.names = FALSE)
  write.csv(rank_edge_counts, file.path(out_dir, paste0(output_prefix, "_rank_edge_counts.csv")), row.names = FALSE)
  
  message(paste("Success! Files for", output_prefix, " saved to", out_dir))
}
  

generate_functional_network_robust(ps_KO, "Pathway_L3", "Network_KO", 1000)
generate_functional_network_robust(ps_Metacyc, "Second_Level", "Network_MetaCyc", 1000)

