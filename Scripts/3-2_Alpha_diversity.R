
####Alpha diversity:
# Load Packages
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(multcompView)
library(emmeans)
library(multcomp)
library(ggpubr)

# Load ps object
setwd("C:/Users/mateu/OneDrive/Area_de_Trabalho/shared_lis/29-spp")
ps <- readRDS("./ASV_tables/For_use/ps_final.rds")
####Alpha diversity:
# Prepare data
prepare_alpha_diversity_df <- function(ps, measures = c("Shannon", "Simpson", "Observed", "Chao1")) {
  
  library(phyloseq)
  library(dplyr)
  library(tibble)
  
  # Alpha diversity
  alpha_df <- estimate_richness(ps, measures = measures) %>%
    rownames_to_column("Sample")
  
  # Metadata (robust extraction)
  metadata_df <- sample_data(ps) %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("Sample")
  
  # Library size
  libsize_df <- data.frame(
    Sample = sample_names(ps),
    LibrarySize = sample_sums(ps)
  ) %>%
    mutate(Log_LibrarySize = log10(LibrarySize))
  
  # Ensure consistent Sample IDs
  metadata_df$Sample <- gsub("-", ".", metadata_df$Sample)
  libsize_df$Sample  <- gsub("-", ".", libsize_df$Sample)
  
  # Merge all
  alpha_df %>%
    left_join(metadata_df, by = "Sample") %>%
    left_join(libsize_df, by = "Sample")
}
# Call function to prepare data
alpha_df <- prepare_alpha_diversity_df(ps)

### In case you need to subset by groups:
subset_melip <- alpha_df %>%
  filter(Genus %in% "Melipona")
subset_nomelip <- alpha_df %>%
  filter(!(Genus %in% "Melipona"))

## Set parameters:
response <- "Simpson"    # Change to desired response (e.g. Simpson, Chao1...)
g_var <- "Genus"        # Change to desired variable in your metadata
out_pref <- paste(response, "_by_", g_var)   # Prefix for output files
dataframe <- subset_nomelip    # Change to specific dataframe in case you subset by groups
out_dir <- "Alpha_diversity" # Output directory

# Diagnostics
run_alpha_diagnostics <- function(df, Metric, Variable, prefix, outdir) {
  
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  library(car)
  # Ensure output directory exists
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  metric_sym   <- rlang::sym(Metric)
  variable_sym <- rlang::sym(Variable)
  # Histogram
  ggsave(
    paste0(outdir, prefix, "_histogram.pdf"),
    ggplot(df, aes(x = !!metric_sym)) +
      geom_histogram(bins = 10, fill = "steelblue", color = "black") +
      facet_wrap(vars(!!variable_sym), scales = "free") +
      labs(title = paste("Histogram of", Metric)),
    width = 8, height = 6
  )
  # QQ plot
  ggsave(
    paste0(outdir, prefix, "_qqplot.pdf"),
    ggqqplot(df[[Metric]], title = paste("QQ plot of", Metric)),
    width = 6, height = 6
  )
  # Shapiro per group
  shapiro_df <- df %>%
    group_by(!!variable_sym) %>%
    filter(n() >= 3) %>%
    summarise(
      p_value = shapiro.test(.data[[Metric]])$p.value,
      .groups = "drop"
    )
  write.csv(
    shapiro_df,
    file.path(outdir, paste0(prefix, "_shapiro.csv")),
    row.names = FALSE
  )
  # Levene test
  levene <- leveneTest(
    df[[Metric]] ~ as.factor(df[[Variable]])
  )
  capture.output(
    levene,
    file = file.path(outdir, paste0(prefix, "_levene.txt"))
  )
  invisible(list(shapiro = shapiro_df, levene = levene))
}
# Call function to run diagnostics
run_alpha_diagnostics(df = dataframe, Metric = response, Variable = g_var, prefix = out_pref, outdir = out_dir)

# Gamma GLM
run_alpha_gamma_glm <- function(df, Metric, Variable, prefix, outdir) {
  
  library(dplyr)
  library(ggplot2)
  library(emmeans)
  library(multcompView)
  library(car) # For Anova()
  
  message("Running Gamma GLM: ", prefix)
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # Data Cleaning
  df <- df %>%
    filter(!is.na(.data[[Metric]]), 
           !is.na(.data[[Variable]]), 
           !is.na(Log_LibrarySize)) %>%
    mutate(!!Variable := factor(.data[[Variable]]))
  
  # Fit Gamma GLM (Additive model for covariate adjustment)
  # Formula: $Metric \sim Variable + Log\_LibrarySize$
  form <- as.formula(paste0(Metric, " ~ ", Variable, " + Log_LibrarySize"))
  
  glm_model <- glm(form, data = df, family = Gamma(link = "log"))
  
  # Significance Testing (Type III Wald Chisq)
  glm_anova <- car::Anova(glm_model, type = "III")
  capture.output(
    summary(glm_model),
    glm_anova,
    file = file.path(outdir, paste0(prefix, "_Gamma_GLM_summary.txt"))
  )
  
  # Post-hoc Comparisons (using response scale)
  emm <- emmeans(glm_model, as.formula(paste("~", Variable)), type = "response")
  contrasts <- contrast(emm, method = "pairwise")
  
  write.csv(as.data.frame(contrasts), 
            file.path(outdir, paste0(prefix, "_pairwise_contrasts.csv")), 
            row.names = FALSE)
  # Compact Letter Display (CLD)
  cld_results <- cld(emm, alpha = 0.05, Letters = letters, adjust = "tukey") %>%
    as.data.frame() %>% dplyr::mutate(Letters = trimws(.group))
  # Dynamically identify the mean column (usually 'response' for Gamma GLMs)
  res_col <- if ("response" %in% names(cld_results)) "response" else "emmean"
  # Finalize cld_df for plotting
  cld_df <- cld_results %>%
    dplyr::select(dplyr::all_of(Variable), Letters, !!sym(res_col)) %>%
    dplyr::rename(Plot_Mean = !!sym(res_col))
  # Merge with max values for label positioning
  plot_df <- df %>%
    group_by(!!sym(Variable)) %>%
    summarise(max_val = max(.data[[Metric]], na.rm = TRUE)) %>%
    left_join(cld_df, by = Variable)
  # Plotting
  p <- ggplot(df, aes_string(x = Variable, y = Metric, fill = Variable)) +
    geom_boxplot(alpha = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3) +
    # Plot the Estimated Marginal Means
    geom_point(data = plot_df, aes_string(x = Variable, y = "Plot_Mean"), 
               color = "red", size = 3) +
    # Add CLD letters above the boxes
    geom_text(data = plot_df, 
              aes_string(x = Variable, y = "max_val * 1.05", label = "Letters"),
              size = 5, fontface = "bold") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Gamma GLM Adjusted:", Metric),
         subtitle = "Red dots = Model-estimated means (adjusted for Library Size)",
         y = Metric)
  
  ggsave(file.path(outdir, paste0(prefix, "_plot.png")), p, width = 7, height = 5)
  ggsave(file.path(outdir, paste0(prefix, "_plot.pdf")), p, width = 7, height = 5)
  message("Finished: ", prefix)
  invisible(list(model = glm_model, anova = glm_anova, emm = emm, plot = p))
}
# Call function to run Gamma_GLM
run_alpha_gamma_glm( df = dataframe, Metric = response, Variable = g_var, prefix = out_pref, outdir = out_dir)


## Save the alpha diversity values
export_alpha_species_tables <- function(df_all, df_melip, df_nomelip, 
                                        species_col = "Sample_species",
                                        measures = c("Shannon", "Simpson", "Observed", "Chao1"), 
                                        outdir = "Alpha_diversity") {
  library(dplyr)
  library(tidyr)
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # Internal helper function to group by species and calculate stats
  get_species_summary <- function(df, dataset_name) {
    # Check if the species column exists
    if (!species_col %in% colnames(df)) {
      stop(paste0("Column '", species_col, "' not found in dataset: ", dataset_name))
    }
    
    df %>%
      dplyr::select(dplyr::all_of(c(species_col, measures))) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(measures), 
        names_to = "Metric", 
        values_to = "Value"
      ) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(species_col)), Metric) %>%
      dplyr::summarise(
        N      = sum(!is.na(Value)),
        Mean   = round(mean(Value, na.rm = TRUE), 3),
        SD     = round(sd(Value, na.rm = TRUE), 3),
        Median = round(median(Value, na.rm = TRUE), 3),
        Min    = round(min(Value, na.rm = TRUE), 3),
        Max    = round(max(Value, na.rm = TRUE), 3),
        .groups = "drop"
      )
  }
  
  # List of datasets to process
  datasets <- list(
    "All_Samples"  = df_all,
    "Melipona"     = df_melip,
    "Non_Melipona" = df_nomelip
  )
  
  results <- list()
  
  # Loop through datasets, generate summary, and save individual CSVs
  for (name in names(datasets)) {
    summary_df <- get_species_summary(datasets[[name]], name)
    filepath <- file.path(outdir, paste0("Alpha_Diversity_Species_", name, ".csv"))
    
    write.csv(summary_df, filepath, row.names = FALSE)
    message("Saved: ", filepath)
    
    results[[name]] <- summary_df
  }
  
  return(results)
}

# Run the function:
species_summaries <- export_alpha_species_tables(
  df_all      = alpha_df, 
  df_melip    = subset_melip, 
  df_nomelip  = subset_nomelip,
  species_col = "Sample_species", # Adjust if your species column name differs
  outdir      = out_dir
)