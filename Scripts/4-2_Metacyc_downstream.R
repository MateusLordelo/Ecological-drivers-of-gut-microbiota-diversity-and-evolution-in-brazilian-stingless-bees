#metacyc analysis,  using the pathway output from picrust

# Load Packages
library(phyloseq)
library(ggplot2)
library(vegan)
library(ANCOMBC)
library(dplyr)

#---------Set file Paths-------
# Read the unstratified pathway abundance table
setwd("C:/Users/mateu/OneDrive/Area_de_Trabalho/shared_lis/29-spp")
pathway_perSample <- read.delim("./Functional_prediction/picrust2_out_my_data/pathways_out/path_abun_unstrat.tsv", check.names = FALSE, row.names = 1)
# Pathway info file
pathway_meta_df <- read.delim("./Functional_prediction/picrust2_out_my_data/pathways_out/consolidated_pathway_metadata.tsv", row.names = 1, stringsAsFactors = FALSE)
# Path to sample metadata
metadata <- read.csv("ASV_tables/For_use/metadata.csv", row.names = 1)


#---------Build phyloseq objects using path info as a taxa_table ------
Pathway_otu <- otu_table(pathway_perSample, taxa_are_rows = TRUE)
pathway_tax_table <- phyloseq::tax_table(as.matrix(pathway_meta_df))
metadata_ps <- sample_data(metadata)
ps_Pathway <- phyloseq(Pathway_otu, metadata_ps, pathway_tax_table)

# --- Check objects ---
print(ps_Pathway)
# Check for sample alignment
all(sample_names(ps_Pathway) %in% rownames(metadata))

## First subset the ps object 
# Subset for Melipona
ps_melip <- subset_samples(ps_Pathway, Genus == "Melipona")
ps_melip <- prune_taxa(taxa_sums(ps_melip) > 0, ps_melip)
# Subset for Non-melipona
ps_nomel <- subset_samples(ps_Pathway, Genus != "Melipona")
ps_nomel <- prune_taxa(taxa_sums(ps_nomel) > 0, ps_nomel)

ps_Pathway2 <- subset_samples(ps_Pathway, !Genus %in% c("Oxytrigona", "Lestrimelitta"))
#---------Set parameters:------
ps <- ps_Pathway2         # Change based on which dataset you want (Pathway, fr, rz)
g_var <- "Genus"        # Change to desired variable in your metadata -> coloured samples
g_var_e <- "Sample_species"      # Change to desired variable in your metadata -> shaped samples
out_dir <- "Metacyc/Beta_diversity" # Output directory
tsuffix <- "Full_Total" # Change to match your analysis
form_var <- "Genus" # Change to desired variable in metadata -> ANCOM-BC formula

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

#------ Functions------

# Helper funtion for Plotting
plot_pcoa_split <- function(ps_obj, color_var, shape_var, out_path, title_suffix = "") {
  
  # 1. Calculate PCoA
  dist_mat <- phyloseq::distance(ps_obj, method = "bray")
  ord_pcoa <- ordinate(ps_obj, method = "PCoA", distance = dist_mat)
  
  # ---- ADD: Extract % variance explained ----
  eig_vals <- ord_pcoa$values$Eigenvalues
  eig_perc <- round(100 * eig_vals / sum(eig_vals), 2)
  
  xlab_txt <- paste0("PCoA1 (", eig_perc[1], "%)")
  ylab_txt <- paste0("PCoA2 (", eig_perc[2], "%)")
  
  # 2. Extract DataFrames manually
  df_samp <- plot_ordination(ps_obj, ord_pcoa, type = "samples", justDF = TRUE)
  df_taxa <- plot_ordination(ps_obj, ord_pcoa, type = "taxa", justDF = TRUE)
  
  # DIAGNOSTIC
  cat("Number of Samples found:", nrow(df_samp), "\n")
  cat("Number of Pathways found:", nrow(df_taxa), "\n")
  
  if (nrow(df_taxa) == 0) {
    stop("No pathways found in the phyloseq object! Check your prune_taxa steps.")
  }
  
  # 3. Panel labeling
  df_samp$panel_type <- "Samples"
  df_taxa$panel_type <- "Pathways"
  
  # Ensure metadata columns exist in taxa DF
  if (!(color_var %in% colnames(df_taxa)))  df_taxa[[color_var]]  <- NA
  if (!(shape_var %in% colnames(df_taxa)))  df_taxa[[shape_var]]  <- NA
  
  # 4. Merge
  plot_data <- dplyr::bind_rows(df_samp, df_taxa)
  
  # 5. Unified color variable
  plot_data <- plot_data %>%
    dplyr::mutate(
      Color_By = dplyr::case_when(
        panel_type == "Samples"  ~ as.character(.data[[color_var]]),
        panel_type == "Pathways" ~ as.character(Second_Level),
        TRUE ~ "Unknown"
      ),
      Color_By = ifelse(is.na(Color_By), "Unclassified", Color_By)
    )
  
  # 6. Plot
  p <- ggplot(plot_data, aes(x = Axis.1, y = Axis.2)) +
    geom_point(
      data = dplyr::filter(plot_data, panel_type == "Samples"),
      aes(color = Color_By, shape = .data[[shape_var]]),
      size = 4
    ) +
    geom_point(
      data = dplyr::filter(plot_data, panel_type == "Pathways"),
      aes(color = Color_By),
      size = 2,
      alpha = 0.6
    ) +
    facet_wrap(~ panel_type, scales = "free") +
    scale_color_manual(values = colors) +
    theme_bw() +
    labs(
      title = paste("MetaCyc PCoA:", title_suffix),
      color = "Group / Pathway Class",
      shape = shape_var,
      x = xlab_txt,
      y = ylab_txt
    ) +
    theme(strip.text = element_text(size = 14, face = "bold"))
  
  # 7. Save
  ggsave(out_path, p, width = 26, height = 6)
  return(p)
}

# Plot and PERMANOVA function
run_beta_pipeline <- function(ps_obj, variable1, variable2, outdir, label) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # Calculate distance
  dist_mat <- phyloseq::distance(ps_obj, method = "bray")
  
  # PERMANOVA (adonis2)
  metadata_df <- data.frame(sample_data(ps_obj))
  form <- as.formula(paste("dist_mat ~", variable1))
  perm <- vegan::adonis2(form, data = metadata_df, permutations = 999)
  write.csv(as.data.frame(perm), file.path(outdir, paste0(label, "_PERMANOVA.csv")))
  
  # Generate Plot
  plot_path <- file.path(outdir, paste0(label, "_PCoA_Split.pdf"))
  plot_pcoa_split(ps_obj, color_var = variable1, shape_var = variable2, 
                  out_path = plot_path, title_suffix = label)
  
  return(list(dist = dist_mat, permanova = perm))
}

# run ANCOMBC2
run_ancombc2_and_export <- function(ps_obj,
                                    fix_formula,
                                    group_var = NULL,
                                    pairwise = FALSE,
                                    out_prefix = "ANCOMBC_result",
                                    p_adj_method = "fdr",
                                    prv_cut = 0.10,
                                    lib_cut = 0,
                                    s0_perc = 0.05,
                                    n_cl = 1,
                                    alpha = 0.05,
                                    tax_level) {
  
  # run ANCOM-BC2
  ancom_res <- ancombc2(
    data = ps_obj,
    assay.type = "counts",
    tax_level = tax_level,    # since ps_obj is genus-collapsed this is safe; if using ASV-level you can still set tax_level="Genus" if taxonomy present
    fix_formula = fix_formula,
    rand_formula = NULL,
    p_adj_method = p_adj_method,
    prv_cut = prv_cut,
    lib_cut = lib_cut,
    s0_perc = s0_perc,
    group = group_var,
    pairwise = pairwise,
    global = FALSE,
    dunnet = FALSE,
    trend = FALSE,
    n_cl = n_cl,
    alpha = alpha,
    verbose = TRUE
  )
  
  # Primary per-taxon results
  res_df <- ancom_res$res
  # If pairwise was requested and produced results, extract res_pair
  res_pair_df <- NULL
  if (pairwise && !is.null(ancom_res$res_pair)) {
    res_pair_df <- ancom_res$res_pair
  }
  
  # Tidy the main results: add genus column if not present
  if (!("taxon" %in% colnames(res_df)) && !is.null(tax_table(ps_obj))) {
    # ensure rownames correspond to genus/taxon names
    res_df$Genus <- rownames(res_df)
  }
  
  # Convert natural log fold-changes to log2 (optional)
  # find columns that start with 'lfc' and convert a copy to log2
  lfc_cols <- grep("^lfc", colnames(res_df), value = TRUE)
  if (length(lfc_cols) > 0) {
    for (cname in lfc_cols) {
      newname <- paste0(cname, "_log2")
      res_df[[newname]] <- res_df[[cname]] / log(2)
    }
  }
  
  # Export main table
  csv_main <- paste0(out_prefix, "_main.csv")
  write.csv(res_df, csv_main, row.names = TRUE)
  message("Wrote: ", csv_main)
  
  # Export pairwise if present
  if (!is.null(res_pair_df)) {
    csv_pair <- paste0(out_prefix, "_pairwise.csv")
    write.csv(res_pair_df, csv_pair, row.names = TRUE)
    message("Wrote: ", csv_pair)
  }
  
  # Return the raw ANCOM object and exported tables
  return(list(ancom_raw = ancom_res, main = res_df, pairwise = res_pair_df))
}

# Final wrapper function
run_full_pipeline <- function(
    ps,
    g_var,
    g_var_e,
    form_var,
    tsuffix,
    out_dir,
    pairwise = TRUE
) {
  
  run_beta_pipeline(
    ps        = ps,
    variable1 = g_var,
    variable2 = g_var_e,
    outdir    = out_dir,
    label     = paste0("Initial_", tsuffix)
  )
  res_ancom <- run_ancombc2_and_export(
    ps_obj      = ps,
    fix_formula = form_var,
    group_var   = g_var,
    pairwise    = pairwise,
    out_prefix  = file.path(out_dir, paste0("ANCOM_", tsuffix)),
    tax_level   = NULL
  )
  # Extract metadata
  meta <- data.frame(sample_data(ps))
  # Extract and order levels of the variable used in ANCOM-BC
  form_levels <- meta[[form_var]] |> 
    as.character() |> 
    unique() |> 
    sort()
  if (length(form_levels) < 2) {
    stop("ANCOM-BC requires at least two levels in ", form_var)
  }
  # Reference is the first level
  ref_level <- form_levels[1]
  comp_levels <- form_levels[-1]
  # Build expected diff_* column names
  diff_cols <- paste0("diff_", form_var, comp_levels)
  # Check which actually exist (defensive programming)
  diff_cols <- diff_cols[diff_cols %in% colnames(res_ancom$main)]
  if (length(diff_cols) == 0) {
    stop("No diff_* columns found in ANCOM-BC output. Check model specification.")
  }
  message("Using ANCOM-BC diff columns: ", paste(diff_cols, collapse = ", "))
  # Identify taxa significant in at least one comparison
  sig_pathways <- res_ancom$main %>%
    filter(apply(across(all_of(diff_cols)), 1, any)) %>%
    pull(taxon)  # This gets the actual pathway names (e.g., PWY-xxx)
  message("Significant pathways retained: ", length(sig_pathways))
  # Add a check for zero results before the Match Error check
  if (length(sig_pathways) == 0) {
    message("!!! WARNING: No significant pathways found at alpha = 0.05.")
    return(list(ancom = res_ancom, ps_pruned = NULL, sig_taxa = NULL))
  }
  
  # Now check for matches
  if (!any(sig_pathways %in% taxa_names(ps))) {
    # Print a diagnostic to see why they don't match
    message("Taxon IDs in result: ", head(sig_pathways, 1))
    message("Taxa names in PS: ", head(taxa_names(ps), 1))
    stop("Match Error: IDs do not match.")
  }
  message("=== 4) Pruning phyloseq object ===")
  ps_sig <- prune_taxa(sig_pathways, ps)
  # Also prune any samples that might now have 0 total reads
  ps_sig <- prune_samples(sample_sums(ps_sig) > 0, ps_sig)
  message("=== 5) Beta-diversity split plot after pruning ===")
  run_beta_pipeline(
    ps        = ps_sig,
    variable1 = g_var,
    variable2 = g_var_e,
    outdir    = out_dir,
    label     = paste0("After_pruning_", tsuffix)
  )
  return(
    list(
      ancom      = res_ancom,
      ps_pruned  = ps_sig,
      sig_taxa   = sig_pathways,
      diff_cols  = diff_cols
    )
  )
}


# ------ Call function ------
run_full_pipeline(ps, g_var, g_var_e, form_var, tsuffix, out_dir)


#--------------Visualise composition--------
library(phyloseq)
library(dplyr)
library(ggplot2)
library(rlang)

ps <- ps_Pathway

# Subset for Melipona
ps_melipona <- subset_samples(ps, Genus == "Melipona")
ps_melipona <- prune_taxa(taxa_sums(ps_melipona) > 0, ps_melipona)

# Subset for Fruit
ps_nomel <- subset_samples(ps, Genus != "Melipona")
ps_nomel <- prune_taxa(taxa_sums(ps_nomel) > 0, ps_nomel)


### set the variable and level of the plot:
ps_obj <- ps_nomel   # Change if you subset, use 'ps' for whole dataset
variable_name <- "Genus"   # must be a variable in you metadata
level_name <- "Second_Level"      # the desired taxonomic rank
cutoff <- 0.99               # set to 90% of total abundance

### Function to plot abundance wrapped by sample type
plot_tax_abundance_pdf <- function(ps, variable, level, abundance_cutoff, outdir = "./Metacyc/Taxonomy_plots") {
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



#---------- Additional code to save ps_Metacyc --------
saveRDS(ps_Pathway, "./ASV_tables/For_use/ps_Metacyc.rds")
