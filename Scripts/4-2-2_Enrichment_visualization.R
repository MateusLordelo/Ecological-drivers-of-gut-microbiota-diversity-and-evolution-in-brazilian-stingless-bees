# ------------------------------------------------------------
# LOAD LIBRARIES
# ------------------------------------------------------------
library(tidyverse)

# ------------------------------------------------------------
# IMPORT MATRIX
# Rows = groups, Columns = pathways
# ------------------------------------------------------------
setwd("C:/Users/mateu/OneDrive/Area_de_Trabalho/shared_lis/29-spp")
df <- read.csv("./ASV_tables/For_use/Metacyc_matrix.csv", row.names = 1, check.names = FALSE)

# Convert to long format (tidy format for ggplot)
df_long <- df %>%
  rownames_to_column("Group") %>%
  pivot_longer(
    cols = -Group,
    names_to = "Pathway",
    values_to = "MeanLFC"
  )

# Inspect
head(df_long)

# Define pathways to exclude manually
exclude_pathways <- c()

# Your predefined color palette
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
  "#ff6347", "yellow", "#e6194b", "#3cb44b", "#ffe119", "#0082c8",
  "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c",
  "#fabebe", "#008080", "#e6beff", "#aa6e28", "#800000",
  "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080",
  "#032bec", "#000000", "#7cfc00", "#40e0d0", "#6495ed",
  "#ff4500", "#32cd32", "#8a2be2", "#ff1493", "#00ff7f",
  "#dc143c", "#ffd700", "#ff69b4", "#00ced1", "#9932cc",
  "#ff6347"
)

# Map pathways → colors
pathways <- unique(df_long$Pathway)

# Ensure enough colors
if(length(colors) < length(pathways)){
  stop("Not enough colors for all pathways!")
}

pathway_colors <- setNames(colors[1:length(pathways)], pathways)


# Exlude pathways

df_filtered <- df_long %>%
  filter(!Pathway %in% exclude_pathways)

# Optional: remove near-zero values (noise filtering)
df_filtered <- df_filtered %>%
  filter(MeanLFC > 0.1)



#Bubble plot

ggplot(df_filtered, aes(x = Pathway, y = Group)) +
  geom_point(
    aes(size = MeanLFC, color = Pathway),
    alpha = 0.8
  ) +
  scale_color_manual(values = pathway_colors) +
  scale_size(range = c(1, 10)) +
  guides(color = "none") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "MetaCyc pathways enriched across genera",
    x = "Pathway",
    y = "Genus",
    size = "Mean LFC"
  )