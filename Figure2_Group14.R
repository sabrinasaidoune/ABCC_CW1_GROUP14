
# MA Plot to create Figure 2 using Table S1

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# Set working directory
setwd("~/Desktop")

# Load Table S1

# Skip the first row (not include the title)
table1 <- read.csv("Copy of Table 1.csv", skip = 1)

# Check column names
names(table1)

# Rename columns to standard names for convenience
df <- table1 %>%
  rename(
    gene = USA300.numbers,
    log2FoldChange = M.value,
    baseMean = Base.mean,
    padj = Adjusted.P.value
  )

# Categorise genes based on thresholds in the paper
df <- df %>%
  mutate(
    sig = case_when(
      padj <= 0.05 & log2FoldChange >= 0.6 ~ "Upregulated",
      padj <= 0.05 & log2FoldChange <= -0.6 ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    plot_group = case_when(
      sig == "Upregulated" ~ Regulon,          # colour by regulon
      sig == "Downregulated" ~ "Repressed",
      TRUE ~ "Not significant"
    )
  )

# Define colours for each regulon

regulon_colours <- c(
  "TetR" = "#E4C22A",      # mustard yellow
  "HypR" = "#8B0000",      # dark red/brown
  "MhqR" = "#CC5500",      # burnt orange
  "CidR" = "#97D700",      # green-yellow
  "QsrR" = "#FFFF00",      # bright yellow
  "CtsR" = "#7B2CBF",      # purple
  "HrcA" = "#C77DFF",      # lavender
  "CymR" = "#009688",      # teal
  "PerR" = "#00B8D9",      # turquoise
  "Fur"  = "#005BBB",      # dark blue
  "CsoR" = "#5BC0EB",      # cyan
  "CstR" = "#89CFF0",      # light blue
  "Zur"  = "#00FF7F",      # mint green
  
  # fallback categories
  "Repressed" = "darkgrey",
  "Not significant" = "lightgrey"
)

# Pick top genes to label (highest fold-change)

top_genes <- df %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(12)

# Create MA plot

p <- ggplot(df, aes(x = A.value, y = log2FoldChange)) +
  geom_point(aes(colour = plot_group), alpha = 0.9, size = 2.5) +
  scale_colour_manual(values = regulon_colours, na.value = "lightgrey") +
  theme_bw(base_size = 14) +
  labs(
    title = "MA plot of S. aureus USA300 under AGXX® stress",
    x = "A-value (log2 base mean)",
    y = "M-value / log2(fold change)",
    colour = "Regulon"
  ) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_hline(yintercept = c(-0.6, 0.6), linetype = "dashed") +
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3.5,
    max.overlaps = 30
  )


# Show plot
print(p)

# Save to file
ggsave("MA_plot_from_TableS1.png", plot = p, width = 10, height = 7, dpi = 300)




