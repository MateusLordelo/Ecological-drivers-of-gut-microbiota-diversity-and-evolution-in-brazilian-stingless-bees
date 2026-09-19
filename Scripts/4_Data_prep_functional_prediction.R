# Prepare files for Picrust2 and Faprotax
# Load libraries
library(phyloseq)
library(biomformat)
library(Biostrings)
library(dplyr)
library(tidyr)
library(tibble)

# Read ps object
setwd("C:/Users/mateu/OneDrive/Area_de_Trabalho/2025-12_Gustavo_UFV")
ps <- readRDS("./ASV_tables/For_use/ps_final.rds")

# For Picrust2
# Extract the OTU table
count_table <- otu_table(ps)
# Extract the reference sequences
ref_sequences <- refseq(ps)
# Convert the count table to a BIOM object
biom_object <- make_biom(data = count_table)
# Write the BIOM object to a file
write_biom(biom_object, "count_table_picrust2.biom")
# Write the reference sequences to a FASTA file
writeXStringSet(ref_sequences, "sequences.fna", format = "fasta")

# For Faprotax
# FAPROTAX needs a table with samples as columns and the LAST column as the taxonomy string
tax_df <- as.data.frame(tax_table(ps)) %>% rownames_to_column("ASV_ID") %>%
  # Create the semicolon-separated taxonomy path
  unite("Taxonomy", Kingdom:Species, sep = ";", na.rm = TRUE, remove = FALSE) %>% dplyr::select(ASV_ID, Taxonomy)
faprotax_table <- as.data.frame(count_table) %>% rownames_to_column("ASV_ID") %>% left_join(tax_df, by = "ASV_ID")
# Export for FAPROTAX
write.table(faprotax_table, "faprotax_asv_taxa_table.tsv", sep = "\t", row.names = FALSE, quote = FALSE)