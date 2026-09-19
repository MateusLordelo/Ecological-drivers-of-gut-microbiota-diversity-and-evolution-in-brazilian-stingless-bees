# Stacked bar plot of ASV novel sequences
# Load Packages
library(tidyverse)

# Read Files
setwd("C:/Users/mateu/OneDrive/Area_de_Trabalho/2025-12_Gustavo_UFV")
df <- read.csv("novelty_vsearch_out/novelty_summary_bins.csv",
               stringsAsFactors = FALSE)

# Ready the data
df$Bin <- factor(df$Bin,
                 levels = c(">=99", ">=97", ">=95", ">=90", "<90", "no_hit"))
df_long <- df %>%
  select(Bin, ASV_Fraction, Read_Fraction) %>%
  pivot_longer(cols = c(ASV_Fraction, Read_Fraction),
               names_to = "Metric",
               values_to = "Fraction")
df_long$Metric <- recode(df_long$Metric,
                         ASV_Fraction = "ASVs",
                         Read_Fraction = "Reads")

# Plot
p <- ggplot(df_long, aes(x = Metric, y = Fraction, fill = Bin)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  labs(
    x = NULL,
    y = "Fraction (%)",
    fill = "Best-hit identity",
    title = "Sequence novelty relative to SILVA reference database"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p

# Save
ggsave("novelty_vsearch_out/novelty_stacked_barplot.pdf",
       p, width = 8, height = 6, dpi = 300)
