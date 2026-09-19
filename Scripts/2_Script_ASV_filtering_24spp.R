library(matrixStats)
library(dplyr)

# 1) Load Data
setwd("C:/Users/mateu/OneDrive/Area_de_Trabalho/2025-12_Gustavo_UFV")
otu <- readRDS("./ASV_tables/For_use/seqtab_all_nochim.rds")    # samples x ASVs (columns = DNA sequences)
taxa <- read.csv("./ASV_tables/Unfiltered_dada_tables/ASV_table_taxonomy_SEQ.csv", row.names = 1, check.names = FALSE)
asv_map <- read.csv("./ASV_tables/Unfiltered_dada_tables/ASV_sequence_to_ID_map.csv", stringsAsFactors = FALSE)
cat("Loaded seqtab:", nrow(otu), "samples x", ncol(otu), "ASVs\n")
cat("Loaded taxa:", nrow(taxa), "ASVs x", ncol(taxa), "tax ranks\n")
cat("Loaded ASV map:", nrow(asv_map), "ASVs\n")

# 2) Filter Step 1 — Taxonomy-based filtering
# Remove non-Bacteria, chloroplast, mitochondria, archaea, eukaryotes
is_bad <- apply(taxa, 1, function(x){
  k <- x["Kingdom"]
  (is.na(k) || k == "" || grepl("Eukaryota|Mitochondria|Chloroplast|Archaea", paste(x, collapse=" "), ignore.case = TRUE))})
keep_asvs <- rownames(taxa)[!is_bad]

# seqtab columns are sequences → match these DNA strings
keep_cols <- colnames(otu) %in% keep_asvs
otu2 <- otu[, keep_cols, drop=FALSE]
taxa2 <- taxa[rownames(taxa) %in% colnames(otu2), , drop=FALSE]
cat("After taxonomic filter:", ncol(otu2), "ASVs retained\n")

# 3) Filter Step 2 — Remove ASVs with <10 total reads
min_total <- 10
asv_totals <- colSums(otu2)
keep_tot <- asv_totals >= min_total
otu3 <- otu2[, keep_tot, drop=FALSE]
taxa3 <- taxa2[rownames(taxa2) %in% colnames(otu3), , drop=FALSE]
cat("After abundance filter (≥", min_total, "reads):", ncol(otu3), "ASVs retained\n")

# 4) Filter Step 3 — Prevalence filter (≥2 samples)
prev_thresh <- 2
prev_counts <- colSums(otu3 > 0)
keep_prev <- prev_counts >= prev_thresh
otu4 <- otu3[, keep_prev, drop=FALSE]
taxa4 <- taxa3[rownames(taxa3) %in% colnames(otu4), , drop=FALSE]
cat("After prevalence filter (≥", prev_thresh, "samples):", ncol(otu4), "ASVs retained\n")

# 5) CREATE UPDATED ASV MAP
filtered_sequences <- colnames(otu4)
filtered_ids <- asv_map$ASV_ID[match(filtered_sequences, asv_map$Sequence)]
asv_map_filtered <- data.frame(
  Sequence = filtered_sequences,
  ASV_ID = filtered_ids,
  stringsAsFactors = FALSE
)

# 6) Prepare count table for CSV (rows = ASV IDs, cols = samples)
asv_counts_filtered <- t(otu4)          # now ASVs as rows
rownames(asv_counts_filtered) <- asv_map_filtered$ASV_ID
# 6.5) Update taxonomy table row names to ASV IDs
# Create a copy to keep taxa4 intact for other analyses if needed
taxa_filtered_ids <- taxa4
# Match the DNA sequence (rownames) to the ID map to get the correct ASV_IDs
taxa_filtered_ids$ASV_ID <- asv_map_filtered$ASV_ID[match(rownames(taxa4), asv_map_filtered$Sequence)]

# Move ASV_ID to the first column and set as rownames
taxa_filtered_ids <- taxa_filtered_ids %>%
  relocate(ASV_ID) %>%
  `rownames<-`(.$ASV_ID)

# 7) Save ALL OUTPUTS as CSV
write.csv(asv_counts_filtered, "ASV_counts_filtered.csv", quote = FALSE)
write.csv(taxa4, "ASV_taxonomy_filtered.csv", quote = FALSE)
write.csv(taxa_filtered_ids, "ASV_taxonomy_IDs_filtered.csv", row.names = FALSE, quote = FALSE)
write.csv(asv_map_filtered, "ASV_sequence_to_ID_map_filtered.csv", row.names = FALSE, quote = FALSE)

# 8) Filtering diagnostics summary CSV
summary_df <- data.frame(
  Step = c("Initial", "After taxonomic filter", "After abundance filter", "After prevalence filter"),
  ASVs = c(ncol(otu), ncol(otu2), ncol(otu3), ncol(otu4)),
  Reads = c(sum(otu), sum(otu2), sum(otu3), sum(otu4)),
  Percent_ASVs_Removed = c(
    0,
    100 * (1 - ncol(otu2)/ncol(otu)),
    100 * (1 - ncol(otu3)/ncol(otu2)),
    100 * (1 - ncol(otu4)/ncol(otu3))
  ),
  Total_Samples = nrow(otu),
  Min_Total_Count_Filter = min_total,
  Prevalence_Threshold = prev_thresh,
  stringsAsFactors = FALSE
)

write.csv(summary_df, "ASV_filtering_diagnostics.csv", row.names = FALSE)
cat("All filtered outputs saved as CSV files.\n")

# 9) Rank-abundance and prevalence diagnostics (visual checks)
pdf(file="Cumulative_abundance_ASVs_per_reads.pdf", width=7, height=7)
# Rank abundance plot (cumulative)
asv_totals_full <- sort(colSums(otu), decreasing=TRUE)
plot(cumsum(asv_totals_full)/sum(asv_totals_full), type='l',
     xlab="ASV rank", ylab="Cumulative fraction of reads",
     main="Cumulative abundance (ASV rank)")
abline(h=0.9, col="red") # where top ASVs give 90% of reads

dev.off()

# Prevalence vs abundance scatter
pdf(file="Prevalence_abundance_ASVs_per_samples.pdf", width=7, height=7)
mean_abund <- colMeans(otu)
prev <- colSums(otu > 0)
plot(log10(mean_abund+1), prev, pch=20, cex=0.6,
     xlab="log10(mean abundance)", ylab="prevalence (# samples)",
     main="Prevalence vs abundance")

dev.off()

# 10) Contribution / assignment analysis
# taxa: rows = ASV (sequence or ID), columns = Kingdom, Phylum, Class, Order, Family, Genus, Species
ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
res <- data.frame(rank = ranks, asv_unassigned = NA, reads_unassigned = NA, stringsAsFactors=FALSE)

for(i in seq_along(ranks)){
  r <- ranks[i]
  # which ASVs have NA or empty at that rank
  unassigned <- is.na(taxa[, r]) | taxa[, r] == "" | taxa[, r] == "unidentified"
  # ASV fraction:
  res$asv_unassigned[i] <- sum(unassigned) / nrow(taxa)
  # reads fraction:
  # ensure OTU matrix columns match taxa rows; if not, match names
  reads_unassigned <- sum(colSums(otu[, colnames(otu) %in% rownames(taxa)[unassigned], drop=FALSE]))
  res$reads_unassigned[i] <- reads_unassigned / sum(otu)
}
print(res)

write.csv(res, file="Taxonomic_Assignment_Summary.csv", row.names=FALSE)
