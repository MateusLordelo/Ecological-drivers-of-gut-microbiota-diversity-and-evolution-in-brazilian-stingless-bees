# Downstream analysis of Picrust2 output
# Load Packages:
library(dplyr)
library(phyloseq)
library(DESeq2)
library(clusterProfiler)
library(KEGGREST)
library(ggplot2)
library(enrichplot)
library(tibble)
library(ggplot2)

setwd("C:/Users/mateu/OneDrive/Area_de_Trabalho/shared_lis/29-spp")
# Path for normalised ASV counts (seqtab_norm.tzv)
seqtab_perASV <- read.delim("./Functional_prediction/picrust2_out_my_data/KO_metagenome_out/seqtab_norm.tsv", check.names = FALSE, row.names = 1)
# Path for ASV NSTI values
nsti_perASV <- read.delim("./Functional_prediction/picrust2_out_my_data/marker_predicted_and_nsti.tsv", check.names = FALSE, row.names = 1)
# Path to KO counts:
KO_perASV <- read.delim("./Functional_prediction/picrust2_out_my_data/KO_predicted.tsv", check.names = FALSE, row.names = 1)
# Path to sample metadata
metadata <- read.csv("./ASV_tables/For_use/metadata.csv", row.names = 1)

# Filtering parameters
abundance_filter <- 0.001 # Set a threshold for low abundance ASVs
nsti_filter <- 0.35 # Set a threshold for high NSTI ASVs -> recommended 0.35 for environmental samples

# Analysis parameters ##### IMPORTANT, ONLY RUN AFTER LOADING ALL HELPER FUNCTIONS####
## First subset the ps object 
# Subset for Melipona
ps_melip <- subset_samples(ps_KO, Genus == "Melipona")
ps_melip <- prune_taxa(taxa_sums(ps_melip) > 0, ps_melip)
# Subset for Non-melipona
ps_nomel <- subset_samples(ps_KO, Genus != "Melipona")
ps_nomel <- prune_taxa(taxa_sums(ps_nomel) > 0, ps_nomel)

# Set parameters
ps_obj <- ps_KO # set to desired dataset object
variable <- "Genus" # Set to desired metadata variable
analysis_name <- "FullDataset" # Set to match dataset prefix


#----------------ASV Filtering------------
rel_abun_table <- sweep(seqtab_perASV, 2, colSums(seqtab_perASV), FUN = "/")
asvs_to_keep_abun <- rownames(rel_abun_table)[rowSums(rel_abun_table >= abundance_filter) >= 1]
asvs_to_keep_nsti <- rownames(nsti_perASV)[nsti_perASV$metadata_NSTI <= nsti_filter]
final_asv_list <- intersect(asvs_to_keep_abun, asvs_to_keep_nsti)
final_count_table <- seqtab_perASV[final_asv_list, ]

#----------------Check filtering results and adjust accordingly:--------------
# --- Initial ASV Counts ---
print("Dimensions of the unfiltered table:")
print(dim(seqtab_perASV))
asv_counts_initial <- colSums(seqtab_perASV > 0)
print("Number of unique ASVs per sample before filtering:")
print(asv_counts_initial)
# --- Final Results ---
print("Dimensions of the final filtered table:")
print(dim(final_count_table))
asv_counts_final <- colSums(final_count_table > 0)
print("Number of unique ASVs per sample in the final table:")
print(asv_counts_final)

#----------------Matrix multiplication for a final table count------
# Ensure ASV IDs match between count table and function tables
common_asvs <- intersect(rownames(final_count_table), rownames(KO_perASV))
final_counts <- final_count_table[common_asvs, ]
KO_filtered  <- KO_perASV[common_asvs, ]

KO_perSample <- t(final_counts) %*% as.matrix(KO_filtered)
# Transpose back so rows = functions, cols = samples
KO_perSample <- t(KO_perSample)
# --- Save results ---
write.table(KO_perSample, "KO_perSample.tsv", sep = "\t", quote = FALSE, col.names = NA)

#----------------Phyloseq objects of Picrust2 outputs-------------------------
### create phyloseq objects
library(phyloseq)
# --- Step 1: Prepare metadata ---
# Ensure sample names match between metadata and count tables
rownames(metadata) <- rownames(metadata)  # already sample IDs
metadata_ps <- sample_data(metadata)
# --- Step 2: Build KO phyloseq object ---
KO_otu <- otu_table(KO_perSample, taxa_are_rows = TRUE)
ps_KO <- phyloseq(KO_otu, metadata_ps)
# --- Step 3: Check objects ---
ps_KO
# Check sample names align perfectly
all(sample_names(ps_KO) %in% rownames(metadata))


#----------------Helper Functions DESeq2 of Picrust2---------------------------------------

# Function to run DESeq2
run_deseq <- function(ps_obj, design_formula, contrast_vector, shrink=TRUE) {
  dds <- phyloseq_to_deseq2(ps_obj, design_formula)
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = contrast_vector)
  
  if (shrink) {
    res <- lfcShrink(dds, contrast = contrast_vector, res = res, type = "ashr")
  }
  
  res <- as.data.frame(res[order(res$padj), ])
  return(res)
}


#----------------Helper Functions Enrichment of DESeq2 output--------------

### 1. Get KEGG pathways and build KO → pathway mapping
pathways <- keggList("pathway", "ko")
pathway_ids <- names(pathways)
# build TERM2GENE table
KO_to_pathway <- do.call(rbind, lapply(pathway_ids, function(pid) {
  k <- keggLink("ko", pid)
  if (length(k) > 0) {
    data.frame(pathway = pid, KO = sub("ko:", "", k))
  }
}))

# Helper function for enrichment
run_and_plot_kegg <- function(deseq_res, 
                              ko_to_pathway_df, 
                              p_cutoff = 0.05, 
                              output_prefix = "Contrast",
                              output_dir = "./kegg_enrichment_results/") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # ---- 1. Run enrichment ----
  sig <- rownames(deseq_res[which(deseq_res$padj < p_cutoff), ])
  universe <- rownames(deseq_res)
  
  if (length(sig) == 0) {
    message(paste("No significant KOs found for:", output_prefix))
    return(NULL)
  }
  
  enr <- enricher(sig,
                  universe = universe,
                  TERM2GENE = ko_to_pathway_df,
                  pvalueCutoff = p_cutoff)
  if (is.null(enr)) return(NULL)
  
  # ---- 2. Prepare DESeq2 results ----
  deseq_res_df <- deseq_res %>%
    as.data.frame() %>%
    rownames_to_column("KO") %>%
    filter(!is.na(log2FoldChange))
  
  sig_with_pathways <- deseq_res_df %>%
    filter(padj < p_cutoff) %>%
    inner_join(ko_to_pathway_df, by = "KO")
  
  # ---- 3. Calculate pathway-level LFC ----
  pathway_summary <- sig_with_pathways %>%
    group_by(pathway) %>%
    summarise(
      mean_LFC = mean(log2FoldChange, na.rm = TRUE),
      median_LFC = median(log2FoldChange, na.rm = TRUE),
      n_sig_KOs = n(),
      .groups = "drop"
    )
  
  # ---- 4. Add descriptions & join ----
  enr@result <- enr@result %>%
    left_join(pathway_summary, by = c("ID" = "pathway"))
  
  enr_df <- as.data.frame(enr)
  
  # Add pathway names from KEGG
  pathways <- keggList("pathway", "ko")
  pathway_map <- data.frame(
    ID = names(pathways),
    PathwayName = pathways,
    stringsAsFactors = FALSE
  )
  
  enr_df <- enr_df %>%
    left_join(pathway_map, by = "ID")
  
  # ---- 5. Save CSV ----
  file_name <- paste0(output_dir, "KEGG_Enrichment_LFC_", output_prefix, ".csv")
  write.csv(enr_df, file_name, row.names = FALSE)
  message(paste("Saved enrichment results to:", file_name))
  
  # ---- 6. Plot ----
  # Take top 20 categories by p.adjust
  df_plot <- enr_df
  
  p <- ggplot(df_plot, aes(x = GeneRatio, y = reorder(PathwayName, GeneRatio))) +
    geom_point(aes(size = Count, color = p.adjust)) +
    geom_text(aes(label = format(round(mean_LFC, 2), nsmall = 2)),
              vjust = 1.5, size = 3.2, color = "black") +
    scale_color_gradient(low = "red", high = "blue", name = "Adj. p") +
    labs(
      x = "Gene Ratio (Count / Background Genes)",
      y = NULL,
      size = "Gene Count",
      title = paste("KO Pathway Enrichment:", output_prefix)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.y = element_text(size = 9),
      legend.title = element_text(size = 10)
    )
  
  # Save plot
  plot_file <- paste0(output_dir, "KEGG_Enrichment_Plot_", output_prefix, ".pdf")
  ggsave(plot_file, plot = p, width = 9, height = 7)
  message(paste("Saved plot to:", plot_file))
  
  return(list(table = enr_df, plot = p))
}


#----------------Helper Master Call the functions-------
#Master function
run_full_analysis <- function(ps_obj, 
                              variable, 
                              numerator, 
                              denominator, 
                              subset_label = "All", 
                              ko_map = KO_to_pathway,
                              output_dir = "./Picrust_downstream/Kegg_enrichment/") {
  
  # 1. Create a distinct prefix for file naming
  # Format: numerator_vs_denominator_subsetLabel
  prefix <- paste(numerator, "vs", denominator, subset_label, sep = "_")
  deseq_folder <- "Picrust_downstream/Deseq_of_Kegg_Picrust"
  # Create the DESeq folder if it doesn't exist
  if (!dir.exists(deseq_folder)) {
    dir.create(deseq_folder, recursive = TRUE)
  }
  message(paste("--- Starting Analysis for:", prefix, "---"))
  
  # 2. Define the contrast vector for DESeq2
  contrast <- c(variable, numerator, denominator)
  
  # 3. Run DESeq2 (using your existing function)
  # design_formula is created dynamically using as.formula
  deseq_results <- run_deseq(ps_obj, 
                             design_formula = as.formula(paste0("~", variable)), 
                             contrast_vector = contrast)
  # --- Save DESeq results to CSV ---
  deseq_file_name <- paste0("DESeq2_", variable, "_", numerator, "_vs_", denominator, ".csv")
  write.csv(deseq_results, file = file.path(deseq_folder, deseq_file_name))
  message(paste("Successfully saved DESeq table to:", deseq_folder))
  # 4. Run KEGG Enrichment and Plotting (using your existing function)
  # This returns a list containing the table and the plot
  kegg_output <- run_and_plot_kegg(deseq_res = deseq_results, 
                                   ko_to_pathway_df = ko_map, 
                                   output_prefix = prefix,
                                   output_dir = output_dir)
  
  # 5. Return everything in case you want to inspect in R
  return(list(
    deseq_table = deseq_results,
    kegg_results = kegg_output$table,
    kegg_plot = kegg_output$plot
  ))
}

#----------------Final call:----
run_automated_analysis <- function(ps_obj, variable, analysis_name) {
  
  # 1. Extract unique levels from the specified metadata variable
  levels <- as.character(unique(sample_data(ps_obj)[[variable]]))
  
  # 2. Generate all possible pairwise combinations (contrasts)
  # combn creates a matrix where each column is a pair
  comparisons <- combn(levels, 2)
  
  # 3. Iterate through the comparisons and run the analysis
  for (i in 1:ncol(comparisons)) {
    cont_1 <- comparisons[1, i]
    cont_2 <- comparisons[2, i]
    
    message(paste("Running:", analysis_name, "| Var:", variable, "|", cont_1, "vs", cont_2))
    
    run_full_analysis(ps_obj, variable, cont_1, cont_2, analysis_name)
  }
}
run_automated_analysis(ps_obj, variable, analysis_name)



#----------Additional code to save ps_KO object----------------

# --- Step 1: Prepare the Pathway Map (as before) ---
pathways <- keggList("pathway", "ko")
pathway_map <- data.frame(
  pathway = names(pathways),
  PathwayName = pathways,
  stringsAsFactors = FALSE
)

# --- Step 2: Create a 3-Slot Taxonomy Map ---
ko_tax_data <- KO_to_pathway %>%
  left_join(pathway_map, by = "pathway") %>%
  group_by(KO) %>%
  summarize(
    # Get all unique pathways for this KO
    all_paths = list(unique(PathwayName)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    # Logic to fill 3 slots and "double up" if missing
    P1 = all_paths[1],
    P2 = ifelse(length(all_paths) >= 2, all_paths[2], P1),
    P3 = ifelse(length(all_paths) >= 3, all_paths[3], P2)
  ) %>%
  ungroup() %>%
  select(KO, P1, P2, P3)

# --- Step 3: Align with Phyloseq Object ---
kos_in_ps <- taxa_names(ps_KO)

tax_matrix <- data.frame(KO = kos_in_ps) %>%
  left_join(ko_tax_data, by = "KO") %>%
  column_to_rownames("KO") %>%
  as.matrix()

# Clean up NAs
tax_matrix[is.na(tax_matrix)] <- "Unclassified"

# --- Step 4: Build and Merge ---
colnames(tax_matrix) <- c("Pathway_L1", "Pathway_L2", "Pathway_L3")
KO_taxonomy <- tax_table(tax_matrix)
ps_KO <- merge_phyloseq(ps_KO, KO_taxonomy)

# Verify the doubling logic
head(tax_table(ps_KO))
# --- Step 5: Save and Check ---
print(ps_KO)
saveRDS(ps_KO, "./ASV_tables/For_use/ps_KO.rds")



