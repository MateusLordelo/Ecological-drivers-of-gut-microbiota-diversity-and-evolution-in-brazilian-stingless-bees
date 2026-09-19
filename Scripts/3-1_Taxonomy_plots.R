library(phyloseq)
library(dplyr)
library(ggplot2)
library(rlang)

setwd("C:/Users/mateu/OneDrive/Area_de_Trabalho/2025-12_Gustavo_UFV")
ps <- readRDS("./ASV_tables/For_use/my_phyloseq_object.rds")

# Subset for Melipona
ps_melipona <- subset_samples(ps, Genus == "Melipona")
ps_melipona <- prune_taxa(taxa_sums(ps_melipona) > 0, ps_melipona)

# Subset for Fruit
ps_nomel <- subset_samples(ps, Genus != "Melipona")
ps_nomel <- prune_taxa(taxa_sums(ps_nomel) > 0, ps_nomel)


### set the variable and level of the plot:
ps_obj <- ps_orb   # Change if you subset, use 'ps' for whole dataset
variable_name <- "Sample_species"   # must be a variable in you metadata
level_name <- "Genus"      # the desired taxonomic rank
cutoff <- 0.99               # set to 90% of total abundance

### Function to plot abundance wrapped by sample type
plot_tax_abundance_pdf <- function(ps, variable, level, abundance_cutoff, outdir = "./Taxonomy_plots") {
  colors <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#78e",
    "#e31a1c", "#fdbf6f", "#f00", "#cab2d6", "#6a3d9a",
    "#ffff99", "#b15928", "#ffd92f", "#b3b3b3", "#ccebc5",
    "#decbe4", "#fed9a6", "#ffffcc", "#e5d8bd", "#32e",
    "yellow", "#e6194b", "#3cb44b", "#ffe119", "#0082c8",
    "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c",
    "#fabebe", "#008080", "#e6beff", "#aa6e28", "#800000",
    "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080",
    "#032bec", "#000000", "#7cfc00", "#40e0d0", "#6495ed",
    "#ff4500", "#32cd32", "#8a2be2", "#ff1493", "#00ff7f",
    "#dc143c", "#ffd700", "#ff69b4", "#00ced1", "#9932cc",
    "#ff6347"
  )
  # 1. Normalize to Relative Abundance (0 to 1) first
  ps_rel <- phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
  ps_melted <- phyloseq::psmelt(ps_rel)
  
  var_sym   <- sym(variable)
  level_sym <- sym(level)
  
  # 2. Identify top taxa globally (Fixed Logic)
  # Using lag() ensures we keep the taxon that actually pushes us past the cutoff
  top_taxa <- ps_melted %>% 
    group_by(!!level_sym) %>% 
    summarize(GlobalAbundance = sum(Abundance), .groups = "drop") %>%
    arrange(desc(GlobalAbundance)) %>% 
    mutate(RelativeAbundance = GlobalAbundance / sum(GlobalAbundance),
           CumulativeAbundance = cumsum(RelativeAbundance)) %>% 
    filter(lag(CumulativeAbundance, default = 0) < abundance_cutoff) %>% 
    pull(!!level_sym)
  
  # 3. Collapse low-abundance taxa into "Other"
  grouped_abundance <- ps_melted %>% 
    mutate(!!level_sym := if_else(!!level_sym %in% top_taxa, as.character(!!level_sym), "Other")) %>% 
    group_by(!!var_sym, !!level_sym) %>% 
    summarize(MeanRelAbundance = mean(Abundance), .groups = "drop")
  
  # 4. Renormalize per group (ensures bars hit exactly 1.0)
  grouped_abundance <- grouped_abundance %>% 
    group_by(!!var_sym) %>% 
    mutate(NormalizedRelAbundance = MeanRelAbundance / sum(MeanRelAbundance)) %>% 
    ungroup()
  
  # Output filename
  if(!dir.exists(outdir)) dir.create(outdir)
  outfile <- file.path(outdir, paste0(tolower(level), "_abundance_grouped_by_", variable, "_of_top_", cutoff*100, "_percent.pdf"))
  
  # Plot
  p <- ggplot(grouped_abundance, aes(x = !!var_sym, y = NormalizedRelAbundance, fill = !!level_sym)) +
    geom_bar(stat = "identity", position = "stack") + 
    labs(x = variable, y = "Relative Abundance", fill = level,
         title = paste("Taxonomic Composition by", level, "(Top", cutoff*100, "%)")) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_fill_manual(values = colors)
  
  # Save as PDF
  ggsave(filename = outfile, plot = p, width = 10, height = 7, units = "in")
  return(p)
}

##Call the function
plot_tax_abundance_pdf(
     ps               = ps_obj,
     variable         = variable_name,
     level            = level_name,
     abundance_cutoff = cutoff
   )



# Code to extract long table format data 

p_obj <- plot_tax_abundance_pdf(
  ps               = ps_obj,
  variable         = variable_name,
  level            = level_name,
  abundance_cutoff = cutoff
)

# 2. Extract the data used for the plot
abundance_table_long <- p_obj$data %>%
  select(!!sym(variable_name), !!sym(level_name), NormalizedRelAbundance) %>%
  rename(RelAbundance = NormalizedRelAbundance)

# 3. View the first few rows
head(abundance_table_long)

# 4. (Optional) Save to CSV
write.csv(abundance_table_long, "./Taxonomy_plots/abundance_table_long_full_data_orb.csv", row.names = FALSE)


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





