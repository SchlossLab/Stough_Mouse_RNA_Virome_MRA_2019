library(tidyverse)
library(ggtree)
library(ape)
library(phangorn)

input_mito_tree <- "data/process/trees/mitovirus/trees/mito_tree.treefile"

mito_tree <- read.tree(file = input_mito_tree)
midpoint_mito_tree <- midpoint(mito_tree)
mito_tree_data <- ggtree(midpoint_mito_tree)
mito_tree_table <- mito_tree_data$data 

mito_tip_labels <- mito_tree_data %>% filter(isTip == TRUE) %>%
  mutate(label = str_replace_all(label, pattern = "_", replacement = " "))

mito_node_labels <- mito_tree_data %>% filter(isTip == FALSE) %>%
  mutate(label = as.numeric(label)) %>%
  filter(label >= 50)

mito_tree_figure <- ggtree(midpoint_mito_tree) +
  geom_hilight(node = 48, fill = "blue", extend = 0.01) +
  geom_hilight(node = 58, fill = "purple", extend = 0.8) +
  geom_hilight(node = 62, fill = "green", extend = 0.35) +
  geom_hilight(node = 83, fill = "pink", alpha = 0.8, extend = 0.01) +
  geom_hilight(node = 72, fill = "chocolate", alpha = 0.7, extend = 0.12) +
  geom_hilight(node = 66, fill = "darkgray", alpha = 0.7, extend = 0.52) +
  geom_text(data = mito_node_labels, aes(label = node), hjust = 1.3, vjust = -0.6, size = 1) +
  geom_tiplab(data = mito_tip_labels, aes(label = label), align = TRUE, size = 2) +
  geom_treescale(x = 0.4, y = 35, color = "red") +
  xlim_tree(12)

mito_tree_figure

ggsave(filename = "data/process/trees/mitovirus/trees/mito_tree.png", plot = mito_tree_figure, width = 8, height = 8)
