library(voronoiTreemap)
library(dplyr)
library(purrr) # For advanced color mapping


# 1. Load data
df <- read.csv("edge_counts.csv")

# 2. Updated Color Mapping (Ensure these match your Area names EXACTLY)
# If a name is missing, it will default to Grey.
base_colors <- c(
  "Cellular Maintenance and Information Processing" = "#E63946",
  "Central Carbon and Energy Metabolism" = "#457B9D",
  "Cofactor and Secondary Metabolism" = "#2A9D8F",
  "Core Biosynthesis" = "#F4A261", # Fixed name
  "Metabolic Flexibility and Environmental Adaptation" = "#9B59B6",
  "Nutrient Acquisition and Catabolism" = "#2ECC71",
  "Stress Resistance and Detoxification" = "#BDC3C7"
)

# 3. Create Shaded Palette
df_shaded <- df %>%
  group_by(Area) %>%
  mutate(
    n_pathways = n(),
    # Use unname() to prevent naming conflicts in the palette
    base_col = ifelse(Area %in% names(base_colors), base_colors[Area], "#BDC3C7"),
    # Generate shades. Using pmax(..., 1) to avoid errors on single-item groups
    color = colorRampPalette(c(unname(base_col[1]), "white"))(n_pathways + 2)[row_number()]
  ) %>%
  ungroup()

# 4. Reorder by Edge Value and Update Labels
area_order <- df_shaded %>%
  group_by(Area) %>%
  summarise(TotalEdges = sum(Edges)) %>%
  arrange(desc(TotalEdges)) %>%
  pull(Area)

df_tree <- df_shaded %>%
  mutate(
    h1 = "Microbiome Functions",
    h2 = factor(Area, levels = area_order),
    h3 = Pathway,
    weight = Edges,
    # CHANGE THIS: Set codes to the numeric value of edges
    codes = as.character(Edges) 
  ) %>%
  select(h1, h2, h3, color, weight, codes) %>%
  arrange(h2, desc(weight)) %>%
  as.data.frame()

# 5. Build and Plot
vt_data <- vt_input_from_df(df_tree)
vt_json <- vt_export_json(vt_data)

vt_d3(
  vt_json,
  title = "Functional Pathway Connectivity",
  legend = TRUE,
  label = TRUE,           # Ensure labels are turned on
  color_border = "#ffffff", 
  size_border = "1.5px"
)