library(tidyverse)
library(ggtree)
library(ape)
library(phangorn)

input_tree <- "data/process/trees/astrovirus/trees/astro_tree.treefile"

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
  geom_text(data = node_labels, aes(label = node), hjust = 1.3, vjust = -0.6, size = 2) +
  geom_tiplab(data = tip_labels, aes(label = label), align = TRUE, size = 3) +
  geom_treescale(x = 0.05, y = 16, color = "red") +
  xlim_tree(3)

tree_figure

ggsave(filename = "data/process/trees/mitovirus/trees/astro_tree.png", plot = tree_figure, width = 8, height = 8)
