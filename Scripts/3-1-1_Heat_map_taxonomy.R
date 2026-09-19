# Load libraries (assuming they are already installed)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

setwd("C:/Users/mateu/OneDrive/Area_de_Trabalho/2025-12_Gustavo_UFV")

# Define the list of bacterial genera you are interested in (exact order for plotting)
selected_genera <- c("Snodgrassella", "Lactobacillus", "Bifidobacterium", "Floricoccus", "Bombella",
                     "Apilactobacillus", "Bombilactobacillus", "Commensalibacter", "Fructobacillus", "Weissella",
                     "Acinetobacter", "Pectinatus", "Leuconostoc", "Apibacter", "Asaia",  "Enterococcus",
                     "Ligilactobacillus", "Swingsia", "Fructilactobacillus","Lonsdalea", "Gilliamella", "Orbaceae_NA")

# --- 1. Load and Prepare Data ---
df_s <- read.csv("ASV_tables/For_use/abundance_table_long_full_data_orb.csv")

# Filter data and ensure all selected genera are present for each sample_Genus
df_filtered_s <- df_s %>%
  filter(Genus %in% selected_genera)

complete_grid_s <- expand.grid(
  Sample_species = unique(df_filtered_s$Sample_species),
  Genus = selected_genera
)

df_complete_s <- complete_grid_s %>%
  left_join(df_filtered_s, by = c("Sample_species", "Genus")) %>%
  mutate(RelAbundance = replace_na(RelAbundance, 0)) # Ensure 0 for missing

# --- 2. Heatmap of Relative Abundance with Adjustments ---

ggplot_heatmap_data_s <- df_complete_s

# Reorder the Genus factor according to 'selected_genera' list
ggplot_heatmap_data_s$Genus <- factor(ggplot_heatmap_data_s$Genus, levels = selected_genera)

# --- REORDERING Sample_species (Corrected Logic) ---
# Get unique bee species names
bee_species_names <- unique(ggplot_heatmap_data_s$Sample_species)

# Separate Melipona species
melipona_species <- sort(grep("^Melipona", bee_species_names, value = TRUE))

# Separate remaining species and sort them alphabetically
other_species <- sort(setdiff(bee_species_names, melipona_species))

# Combine to create the desired order: Melipona first, then others
desired_species_order <- c(melipona_species, other_species)

# Set the factor levels for Sample_species
ggplot_heatmap_data_s$Sample_species <- factor(ggplot_heatmap_data_s$Sample_species, levels = desired_species_order)
# --- END REORDERING ---


# Create a new column for coloring: NA for 0 abundance, otherwise the actual abundance
ggplot_heatmap_data_s <- ggplot_heatmap_data_s %>%
  mutate(RelAbundance_for_color = ifelse(RelAbundance == 0, NA, RelAbundance))

# Set factor levels in reverse order for top-down plotting (Melipona on top)
ggplot_heatmap_data_s$Sample_species <- factor(
  ggplot_heatmap_data_s$Sample_species,
  levels = rev(desired_species_order)
)

p1_s <- ggplot(ggplot_heatmap_data_s, aes(x = Genus, y = Sample_species, fill = RelAbundance_for_color)) +
  geom_tile(color = "black", linewidth = 0.2) +
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlGnBu"),
    name = "Relative\nAbundance",
    limits = c(0, max(ggplot_heatmap_data_s$RelAbundance, na.rm = TRUE)),
    values = scales::rescale(c(0, 0.1, max(ggplot_heatmap_data_s$RelAbundance, na.rm = TRUE))),
    na.value = "white"
  ) +
  # Add text labels to the cells
  geom_text(
    aes(
      label = ifelse(RelAbundance == 0, "0", # Display "0" if abundance is exactly 0
                     ifelse(RelAbundance < 0.01, "<0.01",
                            formatC(RelAbundance, format = "f", digits = 2)))
    ),
    color = "black",
    size = 2.5
  ) +
  labs(
    title = "Relative Abundance of Bacterial Genera in Bee Guts",
    x = "Bacterial Genus",
    y = "Bee Species"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )


print(p1_s)





#####subset pra melipona apenas:::


# Load libraries (assuming they are already installed)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# Define the list of bacterial genera for MELIPONA bees
selected_genera_m <- c("Snodgrassella", "Lactobacillus", "Bifidobacterium", "Floricoccus", "Bombella",
                       "Apilactobacillus", "Commensalibacter", "Weissella", "Fructobacillus", 
                       "Acinetobacter", "Pectinatus", "Leuconostoc",  "Rosenbergiella", "Ligilactobacillus",
                        "Dellaglioa", "Asaia", "Swingsia","Candidatus Kinetoplastibacterium", "Gilliamella", "Orbaceae_NA")

# Define the list of bacterial genera for NON-MELIPONA bees
selected_genera_nm <- c("Snodgrassella", "Lactobacillus", "Bifidobacterium", "Floricoccus", "Bombella",
                        "Apilactobacillus", "Bombilactobacillus", "Commensalibacter", "Fructobacillus", "Weissella",
                        "Acinetobacter", "Pectinatus", "Leuconostoc", "Apibacter", "Zymobacter",  "Enterococcus",
                        "Ligilactobacillus", "Lonsdalea", "Fructilactobacillus","Swingsia", "Gilliamella", "Orbaceae_NA")


# --- 1. Load and Prepare Data (Common Steps) ---
df_original <- read.csv("ASV_tables/For_use/abundance_table_long_full_data_orb.csv")

# Get all unique bee species names for later separation
all_bee_species_names <- unique(df_original$Sample_species) # Use original df for full list

# Separate Melipona species from others
melipona_species_list <- sort(grep("^Melipona", all_bee_species_names, value = TRUE))
non_melipona_species_list <- sort(setdiff(all_bee_species_names, melipona_species_list))

# --- Heatmap for MELIPONA Bees ONLY (_m) ---
message("Generating heatmap for Melipona species...")

# Filter data for Melipona species AND relevant genera
df_m_filtered <- df_original %>%
  filter(Sample_species %in% melipona_species_list,
         Genus %in% selected_genera_m) # Use selected_genera_m here

# Create complete grid for Melipona data
complete_grid_m <- expand.grid(
  Sample_species = unique(df_m_filtered$Sample_species),
  Genus = selected_genera_m # Use selected_genera_m here
)

df_complete_m <- complete_grid_m %>%
  left_join(df_m_filtered, by = c("Sample_species", "Genus")) %>%
  mutate(RelAbundance = replace_na(RelAbundance, 0))

ggplot_heatmap_data_m <- df_complete_m

# Reorder Genus factor using selected_genera_m
ggplot_heatmap_data_m$Genus <- factor(ggplot_heatmap_data_m$Genus, levels = selected_genera_m)

# Set factor levels for Melipona species (alphabetical order, reversed for top-down plot)
ggplot_heatmap_data_m$Sample_species <- factor(
  ggplot_heatmap_data_m$Sample_species,
  levels = rev(sort(unique(ggplot_heatmap_data_m$Sample_species)))
)

# Create color column (NA for 0)
ggplot_heatmap_data_m <- ggplot_heatmap_data_m %>%
  mutate(RelAbundance_for_color = ifelse(RelAbundance == 0, NA, RelAbundance))

p_m <- ggplot(ggplot_heatmap_data_m, aes(x = Genus, y = Sample_species, fill = RelAbundance_for_color)) +
  geom_tile(color = "black", linewidth = 0.2) +
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlGnBu"),
    name = "Relative\nAbundance",
    limits = c(0, max(ggplot_heatmap_data_m$RelAbundance, na.rm = TRUE)),
    values = scales::rescale(c(0, 0.1, max(ggplot_heatmap_data_m$RelAbundance, na.rm = TRUE))),
    na.value = "white"
  ) +
  geom_text(
    aes(
      label = ifelse(RelAbundance == 0, "0",
                     ifelse(RelAbundance < 0.01, "<0.01",
                            formatC(RelAbundance, format = "f", digits = 2)))
    ),
    color = "black",
    size = 2.5
  ) +
  labs(
    title = "Relative Abundance of Bacterial Genera in Melipona Bee Guts",
    x = "Bacterial Genus",
    y = "Melipona Bee Species"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

print(p_m)


# --- Heatmap for NON-MELIPONA Bees ONLY (_nm) ---
message("Generating heatmap for Non-Melipona species...")

# Filter data for Non-Melipona species AND relevant genera
df_nm_filtered <- df_original %>%
  filter(Sample_species %in% non_melipona_species_list,
         Genus %in% selected_genera_nm) # Use selected_genera_nm here

# Create complete grid for Non-Melipona data
complete_grid_nm <- expand.grid(
  Sample_species = unique(df_nm_filtered$Sample_species),
  Genus = selected_genera_nm # Use selected_genera_nm here
)

df_complete_nm <- complete_grid_nm %>%
  left_join(df_nm_filtered, by = c("Sample_species", "Genus")) %>%
  mutate(RelAbundance = replace_na(RelAbundance, 0))

ggplot_heatmap_data_nm <- df_complete_nm

# Reorder Genus factor using selected_genera_nm
ggplot_heatmap_data_nm$Genus <- factor(ggplot_heatmap_data_nm$Genus, levels = selected_genera_nm)

# Set factor levels for Non-Melipona species (alphabetical order, reversed for top-down plot)
ggplot_heatmap_data_nm$Sample_species <- factor(
  ggplot_heatmap_data_nm$Sample_species,
  levels = rev(sort(unique(ggplot_heatmap_data_nm$Sample_species)))
)

# Create color column (NA for 0)
ggplot_heatmap_data_nm <- ggplot_heatmap_data_nm %>%
  mutate(RelAbundance_for_color = ifelse(RelAbundance == 0, NA, RelAbundance))

p_nm <- ggplot(ggplot_heatmap_data_nm, aes(x = Genus, y = Sample_species, fill = RelAbundance_for_color)) +
  geom_tile(color = "black", linewidth = 0.2) +
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlGnBu"),
    name = "Relative\nAbundance",
    limits = c(0, max(ggplot_heatmap_data_nm$RelAbundance, na.rm = TRUE)),
    values = scales::rescale(c(0, 0.1, max(ggplot_heatmap_data_nm$RelAbundance, na.rm = TRUE))),
    na.value = "white"
  ) +
  geom_text(
    aes(
      label = ifelse(RelAbundance == 0, "0",
                     ifelse(RelAbundance < 0.01, "<0.01",
                            formatC(RelAbundance, format = "f", digits = 2)))
    ),
    color = "black",
    size = 2.5
  ) +
  labs(
    title = "Relative Abundance of Bacterial Genera in Non-Melipona Bee Guts",
    x = "Bacterial Genus",
    y = "Non-Melipona Bee Species"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

print(p_nm)