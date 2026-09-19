library(phyloseq)
library(dplyr)
library(Biostrings)

# Read data from CSV files
setwd("C:/Users/mateu/OneDrive/Area_de_Trabalho/2025-12_Gustavo_UFV")
asv_table <- read.csv("./ASV_tables/Filtered_dada_tables/ASV_counts_filtered.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("./ASV_tables/For_use/metadata.csv", row.names = 2)
taxa_table <- read.csv("./ASV_tables/Filtered_dada_tables/ASV_taxonomy_IDs_filtered.csv", row.names = 1, check.names = FALSE)
sequence_map <- read.csv("./ASV_tables/Filtered_dada_tables/ASV_sequence_to_ID_map_filtered.csv")

# Format components for phyloseq
# Convert all columns to character
metadata <- metadata %>%
  mutate(across(everything(), as.character))
asv_otu_table <- otu_table(asv_table, taxa_are_rows = TRUE)
sample_data_df <- sample_data(metadata)
taxa_table_matrix <- tax_table(as.matrix(taxa_table))
# Create the DNAStringSet object
# Ensure 'Sequence' and 'ASV_ID' match your CSV column names
dna_seqs <- DNAStringSet(sequence_map$Sequence)
names(dna_seqs) <- sequence_map$ASV_ID

# Ensure all component names match (a critical step for phyloseq)
# Sort to make sure they are in the same order
# Check sample names
if (!identical(sort(sample_names(asv_otu_table)), sort(sample_names(sample_data_df)))) {
  stop("Sample names do not match.")
}

# Check taxa names
if (!identical(sort(taxa_names(asv_otu_table)), sort(taxa_names(taxa_table_matrix)))) {
  stop("Taxa names do not match.")
}

# Create the phyloseq object
ps <- phyloseq(asv_otu_table, sample_data_df, taxa_table_matrix)
# Add reference sequences
ps <- merge_phyloseq(ps, dna_seqs)

print(ps)

###save the ps object
saveRDS(ps, file = "my_phyloseq_object.rds")
#then load using:
#ps <- readRDS("my_phyloseq_object.rds")