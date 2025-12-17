library(ggplot2)
library(dplyr)
library(ggrepel)

#Loading Table S1

setwd("~/Desktop")   # Change if needed

table1 <- read.csv("Copy of Table 1.csv", skip = 1)

# Inspect to confirm structure
print(names(table1))

#renaming standard variable names

df <- table1 %>% 
  rename(
    gene           = USA300.numbers,     # gene ID
    A.value        = A.value,            # average log2 expression
    log2FoldChange = M.value,            # M-value (log2 fold change)
    padj           = Adjusted.P.value,   # adjusted p-value
    Regulon        = Regulon             # regulon label
  )

#categorise genes according to thresholds

df <- df %>%
  mutate(
    sig = case_when(
      padj <= 0.05 & log2FoldChange >=  0.6 ~ "Upregulated",
      padj <= 0.05 & log2FoldChange <= -0.6 ~ "Downregulated",
      TRUE                                   ~ "Not significant"
    ),
    
    # Groups used for colour:
    plot_group = case_when(
      sig == "Upregulated"   ~ Regulon,      # coloured by regulon
      sig == "Downregulated" ~ "Repressed",  # grey
      TRUE                   ~ "Not significant"
    )
  )

#regulon colours used in the paper

regulon_colours <- c(
  "TetR" = "#E4C22A",
  "HypR" = "#8B0000",
  "MhqR" = "#CC5500",
  "CidR" = "#97D700",
  "QsrR" = "#FFFF00",
  "CtsR" = "#7B2CBF",
  "HrcA" = "#C77DFF",
  "CymR" = "#009688",
  "PerR" = "#00B8D9",
  "Fur"  = "#005BBB",
  "CsoR" = "#5BC0EB",
  "CstR" = "#89CFF0",
  "Zur"  = "#00FF7F",
  
  # fallback groups:
  "Repressed"       = "darkgrey",
  "Not significant" = "lightgrey"
)


# selection of strongly responding genes 

top_genes <- df %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(20)

#To make the plot

p <- ggplot(df, aes(x = A.value, y = log2FoldChange)) +
  
  # Main coloured points
  geom_point(aes(colour = plot_group), size = 2.6, alpha = 0.9) +
  
  # Colour mapping
  scale_colour_manual(values = regulon_colours, na.value = "lightgrey") +
  
  # Figure styling
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 10)
  ) +
  
  # Axis labels and title
  labs(
    title = "Replicated MA Plot of S. aureus USA300 Under AGXX® Stress",
    x = "A-value (log2 mean expression)",
    y = "M-value (log2 fold change)",
    colour = "Regulon"
  ) +
  
  # Horizontal reference lines
  geom_hline(yintercept = 0, colour = "black", size = 0.4) +
  geom_hline(yintercept = c(-0.6, 0.6), linetype = "dashed", size = 0.4) +
  
  # Text labels for top genes
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3.2,
    box.padding = 0.4,
    max.overlaps = 200,
    segment.size = 0.2
  )

print(p)

ggsave(
  "Replicated_Figure2_From_TableS1.png",
  plot = p,
  width = 10,
  height = 7,
  dpi = 300
)




