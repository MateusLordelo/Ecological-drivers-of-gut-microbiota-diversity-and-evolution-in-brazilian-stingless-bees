#####compute geo distanc pairwise1

setwd("C:/Users/mateu/OneDrive/Area_de_Trabalho/2025-12_Gustavo_UFV")
ps <- readRDS("./ASV_tables/For_use/ps_final.rds")

library(vegan)
library(geosphere)#For geographic distances
library(dplyr)
library(tidyr)
library(ggplot2)
library(phyloseq)
# First, add rownames as a column
metadata <- data.frame(sample_data(ps))
metadata$Sample_id <- rownames(metadata)
# Select the SampleID, lon, and lat (and remove NAs)
coords <- metadata %>%
  dplyr::select(Sample_id, Longitude, Latitude) %>%
  drop_na()

# 1. Ensure coords are numeric
coords$Longitude <- as.numeric(as.character(coords$Longitude))
coords$Latitude <- as.numeric(as.character(coords$Latitude))

# 2. Calculate Geographic Distance Matrix (in meters)
# distm expects Longitude in the first column, Latitude in the second
geo_dist_mat <- distm(coords[, c("Longitude", "Latitude")], 
                      fun = distHaversine)

# 3. Add names back
rownames(geo_dist_mat) <- coords$Sample_id
colnames(geo_dist_mat) <- coords$Sample_id

# 4. Optional: Convert to Kilometers
geo_dist_mat_km <- geo_dist_mat / 1000

# Convert to dist object
geo_dist <- as.dist(geo_dist_mat)

####beta diversity:

bray_dist <- phyloseq::distance(ps, method = "wunifrac")

####mantel test:::

# Convert bray_dist (dist object) to a full matrix
bray_mat <- as.matrix(bray_dist)

# Find the shared sample IDs between both matrices
common_ids <- intersect(rownames(geo_dist_mat), rownames(bray_mat))

# Subset both distance matrices
bray_common <- as.dist(bray_mat[common_ids, common_ids])
geo_common <- as.dist(geo_dist_mat[common_ids, common_ids])

# Run the Mantel test
mantel_res <- mantel(bray_common, geo_common, permutations = 999)
print(mantel_res)

#scatterplot for this correlation

library(ggplot2)
df_mantel <- data.frame(
  Host = as.vector(geo_dist_mat_km),
  Microbiome = as.vector(bray_mat)
)

ggplot(df_mantel, aes(x = Host, y = Microbiome)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  labs(title = "Mantel Test Scatterplot",
       x = "Sample Geographic Distance",
       y = "Microbiome Weighted Unifrac Dissimilarity") +
  theme_minimal()
