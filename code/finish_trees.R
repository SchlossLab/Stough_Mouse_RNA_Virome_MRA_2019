library(cowplot)
source("code/edit_astro_tree.R")
source("code/edit_mito_tree.R")

figure1 <- plot_grid(astro_tree_figure, mito_tree_figure, 
                     labels = c("A", "B"), label_size = 30)

ggsave("results/figures/Figure1.tiff", figure1, device = "tiff", width = 12, height = 8)
