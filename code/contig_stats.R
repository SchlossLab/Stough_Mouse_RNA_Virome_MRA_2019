library(tidyverse)

num_hits <- read_tsv(file = "results/tables/rdrp_hits.tsv",
                     col_names = FALSE) %>%
  summarise(n())

scaffold_stats <- read_tsv(file = "results/tables/contig_stats_raw.tsv",
                           col_names = FALSE) %>%
  separate(X1, into = c("a", "b", "c", "d", "e", "length", "g", "coverage",
                        "j", "k"), sep = '_') %>%
  unite("id", a:d, sep = '_') %>%
	mutate(length = as.numeric(length), coverage = as.numeric(coverage)) %>%
  select(id, length, coverage)

write.table(scaffold_stats, file = "results/tables/scaffold_stats.tsv", sep = '_')

num_scaffolds <- scaffold_stats %>%
  summarise(n())

astro_length <- scaffold_stats %>%
  filter(id == "cefoperazone_630_NODE_270") %>%
  select(length)

astro_coverage <- scaffold_stats %>%
  filter(id == "cefoperazone_630_NODE_270") %>%
  select(coverage)

mito_length <- scaffold_stats %>%
  filter(id %in% c("cefoperazone_mock_NODE_3040",
                   "clindamycin_mock_NODE_4406",
                   "germ_free_NODE_1298",
                   "streptomycin_630_NODE_11363",
                   "streptomycin_mock_NODE_4960")) %>%
  summarise("low" = min(as.numeric(length)), "high" = max(as.numeric(length)))

mito_cov <- scaffold_stats %>%
  filter(id %in% c("cefoperazone_mock_NODE_3040",
                   "clindamycin_mock_NODE_4406",
                   "germ_free_NODE_1298",
                   "streptomycin_630_NODE_11363",
                   "streptomycin_mock_NODE_4960")) %>%
  summarise("low" = min(as.numeric(coverage)), "high" = max(as.numeric(coverage)))
