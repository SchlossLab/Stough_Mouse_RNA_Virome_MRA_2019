library(tidyverse)
library(ggtree)
library(ape)
library(phangorn)

get_tip_labels <- function(mito_tree){
mito_tip_labels <- ggtree(midpoint(mito_tree)) %>% filter(isTip == TRUE) %>%
  mutate(label = str_replace_all(label, pattern = "_", replacement = " "))
}

get_node_labels <- function(mito_tree){
mito_node_labels <- ggtree(midpoint(mito_tree)) %>% filter(isTip == FALSE) %>%
  mutate(label = as.numeric(label)) %>%
  filter(label >= 50)
}
