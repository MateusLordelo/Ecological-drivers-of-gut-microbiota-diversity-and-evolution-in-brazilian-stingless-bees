
#----- Old ------
#####analysis with the bee tree:
library(ggtree)
library(ggplot2)
library(phyloseq)
library(ape)
setwd("C:/Users/mateu/OneDrive/Area_de_Trabalho/2025-12_Gustavo_UFV")
ps <- readRDS("./ASV_tables/For_use/ps_final.rds")
tree_file <- readLines("./ASV_tables/For_use/Bee_tree.nwk")
writeLines(tree_file[1], "./ASV_tables/For_use/bee_tree_best.nwk")  # save only the first tree
bee_tree <- read.tree("./ASV_tables/For_use/bee_tree_best.nwk")

# Plot with ggtree and uniform branch lengths
p <- ggtree(bee_tree, branch.length = "none") +
  geom_tiplab(size = 2.5) +  # Tip labels (species names)
  
  # Internal node numbers
  geom_nodelab(aes(label = node), 
               hjust = -1.3, 
               size = 2.5, 
               color = "blue") +
  
  # Branch length labels (on internal branches)
  geom_label2(aes(subset = !isTip, 
                  label = round(branch.length, 3)), 
              size = 2, 
              fill = "white", 
              label.padding = unit(0.15, "lines")) +
  
  ggtitle("Bee Tree with Uniform Branches and Labeled Nodes") +
  theme_tree2()

# Automatically add right margin padding so tip labels don’t get cropped
max_depth <- max(node.depth(bee_tree))
plot(p + xlim(NA, max_depth + 5))




##compute bray-curtis::


library(phyloseq)
library(vegan)
library(dplyr)

# Get distance matrix between samples
bray_dist <- phyloseq::distance(ps, method = "bray")
bray_mat <- as.matrix(bray_dist)

# Extract sample metadata
meta <- data.frame(sample_data(ps))

# Get average dissimilarity between species (centroid-based)
species_list <- unique(meta$Sample_species)
species_dist <- matrix(NA, nrow = length(species_list), ncol = length(species_list),
                       dimnames = list(species_list, species_list))

for (i in species_list) {
  for (j in species_list) {
    samples_i <- rownames(meta[meta$Sample_species == i, ])
    samples_j <- rownames(meta[meta$Sample_species == j, ])
    dists <- bray_mat[samples_i, samples_j, drop = FALSE]
    if (all(is.na(dists))) {
      species_dist[i, j] <- NA
    } else {
      species_dist[i, j] <- mean(dists, na.rm = TRUE)
    }
  }
}

species_dist <- as.dist(species_dist)

# Extract the current species labels
species_labels <- attr(species_dist, "Labels")

# Substitute spaces for underscores
species_labels <- gsub(" ", "_", species_labels)

# Assign them back
attr(species_dist, "Labels") <- species_labels

# Convert to matrix
species_dist_mat <- as.matrix(species_dist)

# Extract host tree distance matrix
host_dist <- cophenetic(bee_tree)
host_dist_mat <- as.matrix(host_dist)

# Find common species
common_species <- intersect(rownames(host_dist_mat), rownames(species_dist_mat))

# Subset both matrices
species_dist_mat <- species_dist_mat[common_species, common_species]
host_dist_mat <- host_dist_mat[common_species, common_species]





############tentar extrair a arvore de dissimilaridade::
library(ggtree)
library(ape)
library(ape)

# Perform hierarchical clustering on the Bray-Curtis dissimilarity matrix
hc_microbiome <- hclust(species_dist, method = "average")  # or method = "ward.D2", "complete", etc.

# Convert to phylo object (tree format)
microbiome_tree <- as.phylo(hc_microbiome)

p <- ggtree(microbiome_tree, branch.length = "none") +
  geom_tiplab(size = 2.5) +
  geom_nodelab(aes(label = node), hjust = -1.3, size = 2.5, color = "blue") +
  geom_label2(aes(subset = !isTip, label = round(branch.length, 3)), 
              size = 2, fill = "white", label.padding = unit(0.15, "lines")) +
  ggtitle("Microbiome Tree with Uniform Branch Lengths and Labeled Distances") +
  theme_tree2()
plot(p)
max_depth <- max(node.depth(microbiome_tree))
p + xlim(NA, max_depth + 5)  # Add 5 units of padding to the right





####mantels test for correlation
library(vegan)
mantel_result <- mantel(host_dist_mat, species_dist_mat, permutations = 999)
print(mantel_result)



#####procrustes graph::

procrustes_result <- protest(host_dist_mat, species_dist_mat, permutations = 999)
plot(procrustes_result)




####dendograms???


library(ape)
library(dendextend)

# Convert both matrices to dendrograms
host_dend <- as.dendrogram(hclust(as.dist(host_dist_mat)))
microbiome_dend <- as.dendrogram(hclust(as.dist(species_dist_mat)))

# Combine into dendlist for untangling
dend_list <- dendlist(host_dend, microbiome_dend)

# Find the optimal layout (minimize entanglement)
dend_list <- untangle(dend_list, method = "step2side")  # or "step1side", "step2side", "ladderize", try all

# Assuming dend_list[[1]] and dend_list[[2]] are your dendrograms
dend_left_ranked <- rank_branches(dend_list[[1]])
dend_right_ranked <- rank_branches(dend_list[[2]])

pdf("tanglegram_plot_ranked.pdf", width = 12, height = 8)
tanglegram(dend_left_ranked, dend_right_ranked,
           main_left = "Bee Phylogeny",
           main_right = "Microbiome Dissimilarity",
           common_subtrees_color_lines = TRUE,
           margin_inner = 12,
           highlight_distinct_edges = FALSE,
           lwd = 2,
           branches_lwd = 1,
           lab.cex = 0.7,
           lab.gap = 0.02,
           lab.pos = 1
)
dev.off()


# 1. Prepare the phylo object
micro_phylo <- as.phylo(hclust(as.dist(species_dist_mat)))

# 2. Extract the REAL distances BEFORE we plot
# edge.length contains the actual Bray-Curtis values
real_distances <- round(micro_phylo$edge.length, 3)

# 3. Plot as a CLADOGRAM
# 'use.edge.length = FALSE' forces all branches to be the same length
plot(micro_phylo, 
     use.edge.length = FALSE, 
     main = "Microbiome Dissimilarity (Uniform Branches)",
     cex = 0.8)

# 4. Add the REAL values to the edges
# We use edgelabels instead of nodelabels to place them on the branches
edgelabels(real_distances, 
           adj = c(0.5, -0.5), 
           frame = "none", 
           cex = 0.7, 
           col = "darkred")

####heatmaps???


library(pheatmap)

# Heatmaps using same species order
species_order <- rownames(host_dist_mat)

# Normalize both distance matrices
scale_matrix <- function(mat) {
  (mat - min(mat)) / (max(mat) - min(mat))
}

host_scaled <- scale_matrix(host_dist_mat)
microbiome_scaled <- scale_matrix(species_dist_mat)

# Compute similarity score: 1 - abs difference
pairwise_agreement <- 1 - abs(host_scaled - microbiome_scaled)

# Plot heatmap of correlation between phylogeny and microbiome
# Plot heatmap of correlation between phylogeny and microbiome
pheatmap(pairwise_agreement[species_order, species_order],
         main = "Phylogeny vs Microbiome Correlation (per pair)",
         color = colorRampPalette(c("yellow", "white", "darkgreen"))(100),
         clustering_method = "average",
         display_numbers = TRUE,  # This will display the numbers in the cells
         number_color = "black",  # Set the color of the numbers
         fontsize_number = 8,     # Adjust the font size as needed
         number_format = "%.2f"   # Format the numbers to two decimal places
)



####scatterplot??

library(ggplot2)
df_mantel <- data.frame(
  Host = as.vector(host_dist_mat),
  Microbiome = as.vector(species_dist_mat)
)

ggplot(df_mantel, aes(x = Host, y = Microbiome)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  labs(title = "Mantel Test Scatterplot",
       x = "Host Phylogenetic Distance",
       y = "Microbiome Bray-Curtis Dissimilarity") +
  theme_minimal()




####save the custom scale so we can use to colour other stuff:

library(tidyverse)

# Step 1: Convert matrix to data frame
pairwise_df <- as.data.frame(as.table(pairwise_agreement))

# Step 2: Rename columns for clarity
colnames(pairwise_df) <- c("Species1", "Species2", "CorrelationScore")

# Step 3: Remove duplicate pairs (keep upper triangle only)
pairwise_df <- pairwise_df %>%
  filter(as.character(Species1) < as.character(Species2)) %>%
  arrange(desc(CorrelationScore))

# Step 4: Save to CSV
write.csv(pairwise_df, "pairwise_phylogeny_microbiome_agreement.csv", row.names = FALSE)


#####this can be used to colour other stuff, like maybe the tabglegram's lines. 






#do the same thing, but using unifrac weighted



# Calculate weighted UniFrac distance between samples
unifrac_dist <- phyloseq::distance(ps, method = "unifrac")  # weighted = TRUE by default

unifrac_mat <- as.matrix(unifrac_dist)


# Extract sample metadata
meta <- data.frame(sample_data(ps))

# Get average dissimilarity between species (centroid-based)
species_list <- unique(meta$Sample_species)
species_dist <- matrix(NA, nrow = length(species_list), ncol = length(species_list),
                       dimnames = list(species_list, species_list))

for (i in species_list) {
  for (j in species_list) {
    samples_i <- rownames(meta[meta$Sample_species == i, ])
    samples_j <- rownames(meta[meta$Sample_species == j, ])
    dists <- unifrac_mat[samples_i, samples_j, drop = FALSE]
    if (all(is.na(dists))) {
      species_dist[i, j] <- NA
    } else {
      species_dist[i, j] <- mean(dists, na.rm = TRUE)
    }
  }
}

species_dist <- as.dist(species_dist)

# Extract the current species labels
species_labels <- attr(species_dist, "Labels")

# Substitute spaces for underscores
species_labels <- gsub(" ", "_", species_labels)

# Assign them back
attr(species_dist, "Labels") <- species_labels

# Convert to matrix
species_dist_mat <- as.matrix(species_dist)

# Extract host tree distance matrix
host_dist <- cophenetic(bee_tree)
host_dist_mat <- as.matrix(host_dist)

# Find common species
common_species <- intersect(rownames(host_dist_mat), rownames(species_dist_mat))

# Subset both matrices
species_dist_mat <- species_dist_mat[common_species, common_species]
host_dist_mat <- host_dist_mat[common_species, common_species]



####mantels test for correlation
library(vegan)
mantel_result <- mantel(host_dist_mat, species_dist_mat, permutations = 999)
print(mantel_result)



#####procrustes graph::

procrustes_result <- protest(host_dist_mat, species_dist_mat, permutations = 999)
plot(procrustes_result)




####dendograms???


library(ape)
library(dendextend)

# Convert both matrices to dendrograms
host_dend <- as.dendrogram(hclust(as.dist(host_dist_mat)))
microbiome_dend <- as.dendrogram(hclust(as.dist(species_dist_mat)))

# Combine into dendlist for untangling
dend_list <- dendlist(host_dend, microbiome_dend)

# Find the optimal layout (minimize entanglement)
dend_list <- untangle(dend_list, method = "step2side")  # or "step1side", "step2side", "ladderize", try all

# Plot the optimized tanglegram
pdf("tanglegram_plot-unifrac.pdf", width = 12, height = 8)

tanglegram(dend_list[[1]], dend_list[[2]],
           main_left = "Bee Phylogeny",
           main_right = "Microbiome Dissimilarity",
           common_subtrees_color_lines = TRUE,
           margin_inner = 12,
           highlight_distinct_edges = FALSE,
           lwd = 2,
           lab.cex = 0.7,
           lab.gap = 0.02,
           lab.pos = 1
)

dev.off()



####heatmaps???


library(pheatmap)

# Heatmaps using same species order
species_order <- rownames(host_dist_mat)

pheatmap(host_dist_mat[species_order, species_order],
         main = "Host Phylogenetic Distance", clustering_method = "average")

pheatmap(species_dist_mat[species_order, species_order],
         main = "Microbiome Unifrac Dissimilarity", clustering_method = "average")


#######trying a combined approach:::

####my idea was this::For each species pair: compare how similar their phylogenetic distance and microbiome dissimilarity are.
#→ If both are high or both are low → correlation
#→ If one is high and the other is low → mismatch

# Normalize both distance matrices
scale_matrix <- function(mat) {
  (mat - min(mat)) / (max(mat) - min(mat))
}

host_scaled <- scale_matrix(host_dist_mat)
microbiome_scaled <- scale_matrix(species_dist_mat)

# Compute similarity score: 1 - abs difference
pairwise_agreement <- 1 - abs(host_scaled - microbiome_scaled)

# Plot heatmap of correlation between phylogeny and microbiome
pheatmap(pairwise_agreement[species_order, species_order],
         main = "Phylogeny vs Microbiome Correlation (per pair)",
         color = colorRampPalette(c("yellow", "white", "darkgreen"))(100),
         clustering_method = "average")


####scatterplot??

library(ggplot2)
df_mantel <- data.frame(
  Host = as.vector(host_dist_mat),
  Microbiome = as.vector(species_dist_mat)
)

ggplot(df_mantel, aes(x = Host, y = Microbiome)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  labs(title = "Mantel Test Scatterplot",
       x = "Host Phylogenetic Distance",
       y = "Microbiome Unifrac Dissimilarity") +
  theme_minimal()

