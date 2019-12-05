library(tidyverse)
library(ggtree)
library(ape)
library(phangorn)

input_tree <- "data/process/trees/mitovirus/trees/mito_tree.treefile"

tree <- read.tree(file = input_tree)
midpoint_tree <- midpoint(tree)
tree_data <- ggtree(midpoint_tree)
tree_table <- tree_data$data 

tip_labels <- tree_data %>% filter(isTip == TRUE) %>%
  mutate(label = str_replace_all(label, pattern = "_", replacement = " "))

node_labels <- tree_data %>% filter(isTip == FALSE) %>%
  mutate(label = as.numeric(label)) %>%
  filter(label >= 50)

tree_figure <- ggtree(midpoint_tree) +
  geom_hilight(node = 48, fill = "blue", extend = 0.01) +
  geom_hilight(node = 58, fill = "purple", extend = 0.88) +
  geom_hilight(node = 62, fill = "green", extend = 0.3) +
  geom_hilight(node = 83, fill = "pink", alpha = 0.8, extend = 0.01) +
  geom_hilight(node = 72, fill = "chocolate", alpha = 0.7, extend = 0.08) +
  geom_hilight(node = 66, fill = "darkgray", alpha = 0.7, extend = 0.52) +
  geom_text(data = node_labels, aes(label = node), hjust = 1.3, vjust = -0.6, size = 2) +
  geom_tiplab(data = tip_labels, aes(label = label), align = TRUE, size = 3) +
  geom_treescale(x = 0.4, y = 35, color = "red") +
  xlim_tree(12)

tree_figure

ggsave(filename = "data/process/trees/mitovirus/trees/mito_tree.png", plot = tree_figure, width = 8, height = 8)
