library(phyloseq)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("C:/Users/mateu/OneDrive/Area_de_Trabalho/shared_lis/29-spp")

# ------------ Load the files ----------
ps <- readRDS("./ASV_tables/For_use/ps_final.rds")
### For orbaceae_NA
ps_orb <- ps
# 1. Extrair a tabela taxonômica
tax <- as.data.frame(tax_table(ps_orb))
# 2. Definir os IDs dos ASVs que queremos renomear
asvs_to_edit <- c("ASV_67", "ASV_127", "ASV_135", "ASV_143", "ASV_260", 
                  "ASV_272", "ASV_327", "ASV_329", "ASV_341", "ASV_410", 
                  "ASV_597", "ASV_659", "ASV_900", "ASV_970", "ASV_1017", "ASV_1217")
# 3. Alterar o gênero para "Orbaceae_NA" apenas para esses IDs
tax[rownames(tax) %in% asvs_to_edit, "Genus"] <- "Orbaceae_NA"
# 4. Inserir a tabela modificada de volta no objeto phyloseq
tax_table(ps_orb) <- tax_table(as.matrix(tax))
# Verificação
head(tax_table(ps_orb)[asvs_to_edit, ])




#---------------
# 1. Your color palette assigned to each unique species
species_colors <- c(
  "Trigona hyalinata"           = "#1f77b4",
  "Tetragonisca angustula"      = "#ff7f0e",
  "Tetragona clavipes"          = "#2ca02c",
  "Schwarziana quadripunctata"  = "#d62728",
  "Scaptotrigona xanthotricha"  = "#9467bd",
  "Scaptotrigona polysticta"    = "#8c564b",
  "Scaptotrigona depilis"       = "#e377c2",
  "Scaptotrigona bipunctata"    = "#7f7f7f",
  "Plebeia phrynostoma"         = "#bcbd22",
  "Plebeia lucii"               = "#17becf",
  "Plebeia droryana"            = "#a6cee3",
  "Partamona helleri"           = "#1f78b4",
  "Oxytrigona tataira"          = "#b2df8a",
  "Nannotrigona testaceicornis" = "#33a02c",
  "Lestrimelitta limao"          = "#78e",
  "Frieseomelitta varia"        = "#e31a1c",
  "Friesella schrottkyi"      = "#fdbf6f",
  "Cephalotrigona capitata"     = "#f00",
  "Melipona scutellaris"        = "#cab2d6",
  "Melipona quadrifasciata"     = "#6a3d9a",
  "Melipona mondury"            = "#ffff99",
  "Melipona marginata"          = "#b15928",
  "Melipona mandacaia"          = "#ffd92f",
  "Melipona fulva"              = "#b3b3b3",
  "Melipona captiosa"           = "#ccebc5",
  "Melipona capixaba"           = "#decbe4",
  "Melipona bicolor"            = "#fed9a6",
  "Melipona asilvai"            = "#ffffcc",
  "Trigona spinipes"            = "#e5d8bd"
)

# 2. Distinct R plotting symbols (0-14) assigned to each unique genus
genus_shapes <- c(
  "Trigona"         = 15, # Filled square
  "Tetragonisca"    = 16, # Filled circle
  "Tetragona"       = 17, # Filled triangle up
  "Schwarziana"     = 18, # Filled diamond
  "Scaptotrigona"   = 21, # Circle with border
  "Plebeia"         = 22, # Square with border
  "Partamona"       = 23, # Diamond with border
  "Oxytrigona"      = 24, # Triangle up with border
  "Nannotrigona"    = 25, # Triangle down with border
  "Lestrimelitta"    = 0,  # Open square
  "Frieseomelitta"  = 1,  # Open circle
  "Friesella"       = 2,  # Open triangle
  "Cephalotrigona"  = 5,  # Open diamond
  "Melipona"        = 8   # Star/Plus cross
)
#------------
get_subset_tree <- function(physeq, target_strains, rank = "Genus") {
  # Extract taxonomy table
  tax_df <- as.data.frame(tax_table(physeq))
  
  # Find ASV IDs matching your target criteria
  matched_asvs <- rownames(tax_df[tax_df[[rank]] %in% target_strains, ])
  
  if (length(matched_asvs) == 0) {
    stop("No matching strains found for the specified targets.")
  }
  
  # Prune phyloseq object to keep only target ASVs and relevant tree branches
  sub_ps <- prune_taxa(matched_asvs, physeq)
  return(sub_ps)
}


plot_bee_tree <- function(sub_ps, host_col = "Sample_species") {
  # 1. Extract OTU matrix
  otu_mat <- as.matrix(otu_table(sub_ps))
  if (!taxa_are_rows(sub_ps)) { otu_mat <- t(otu_mat) }
  otu_df <- as.data.frame(otu_mat, check.names = FALSE)
  
  # 2. Extract metadata and clean short codes
  meta <- data.frame(sample_data(sub_ps), check.names = FALSE)
  meta$SampleID <- sample_names(sub_ps) 
  meta$Host_Code <- stringr::word(meta$SampleID, 1, sep = "-")
  
  # 3. Map ASVs to sharing hosts
  asv_host_map <- otu_df %>%
    mutate(ASV = rownames(.)) %>%
    pivot_longer(-ASV, names_to = "SampleID", values_to = "Abundance") %>%
    filter(Abundance > 0) %>%
    left_join(meta, by = "SampleID") %>%
    group_by(ASV) %>%
    summarise(
      Composite_Abbr = paste(sort(unique(Host_Code)), collapse = "_"),
      .groups = "drop"
    )
  
  # 4. Get pristine tip order for sequential numbering
  tree_tips <- phy_tree(sub_ps)$tip.label
  
  tip_metadata <- data.frame(
    label = tree_tips, 
    Num_ID = seq_along(tree_tips), 
    stringsAsFactors = FALSE
  ) %>%
    left_join(asv_host_map, by = c("label" = "ASV")) %>%
    mutate(Composite_Label = paste0(Num_ID, "_", Composite_Abbr))
  
  # 5. Plot monochrome tree matching classic publication styles
  p <- ggtree(phy_tree(sub_ps)) %<+% tip_metadata +
    # Branch lengths
    geom_text2(aes(label = ifelse(!is.na(branch.length) & branch.length > 0, 
                                  sprintf("%.3f", branch.length), "")), 
               vjust = -0.4, hjust = 0.5, size = 2.2, color = "black") +
    # Tip circle highlight
    geom_tippoint(size = 1.5, color = "black", shape = 16) +
    # Standard black text label
    geom_tiplab(aes(label = Composite_Label), color = "black", 
                hjust = -0.15, size = 3.2) +
    hexpand(0.4) + 
    theme_tree2()
  
  return(p)
}
# Step 1: Subset your tree for target genera
bac_list <- c("Bifidobacterium", "Lactobacillus", "Snodgrassella", "Orbaceae_NA",
              "Floricoccus", "Bombilactobacillus", "Bombella", "Acinetobacter", "Commensalibacter") 
ps_subset <- get_subset_tree(ps, target_strains = bac_list, rank = "Genus")

my_plots <- list()

for (bac in bac_list) {
  message("Generating plot for: ", bac)
  
  tryCatch({
    # 1. Subset for the single genus
    ps_subset <- get_subset_tree(ps_orb, target_strains = bac, rank = "Genus")
    
    # 2. Guard rail: Check if tree has enough tips to be drawn
    if (ntaxa(ps_subset) < 2) {
      message("   -> Skipped ", bac, " because it only contains 1 ASV (ggtree requires >= 2 tips).")
      next
    }
    
    # 3. Generate base plot
    p <- plot_bee_tree(ps_subset, host_col = "Sample_species")
    
    # 4. Add custom styles
    p_custom <- p + 
      labs(title = paste("Phylogenetic Tree of Genus:", bac)) +
      theme(
        legend.position = "none" # Completely removes any accidental legend elements
      )
    
    # 5. Save and display
    my_plots[[bac]] <- p_custom
    print(p_custom)
    
  }, error = function(e) {
    # conditionMessage(e) successfully extracts modern Tidyverse/ggtree error details
    message("   -> Skipped ", bac, " due to error: ", conditionMessage(e))
  })
}