args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)

orf_file = args[1]
ann_file = args[2]
orf_name = str_remove(basename(orf_file), "_orfs.sco")
ann_path = dirname(ann_file)

pred_orfs <- read_delim(file = orf_file, delim = "_", 
                        col_names = c("orf_num", "start", "stop", "strand"),
                        col_types = "ciic") %>%
  filter(str_detect(orf_num, "^>")) %>%
  mutate(orf_num = str_remove_all(orf_num, ">")) %>%
  mutate(orfname = orf_name) %>%
  unite(id, orfname, orf_num, sep = "_") %>%
  select(-strand)

annotations <- read_tsv(file = ann_file, col_names = c("id", "ip_code", "length", "method", "accession", 
                                                       "description", "start_base", "stop_base", "evalue", 
                                                       "status", "date", "ip_accession", "ip_description")) %>%
                 select(id, method, description, evalue) %>%
                 filter(!is.na(description)) %>%
                 filter(!is.na(as.numeric(evalue))) %>%
                 select(-evalue)

annotation_table <- right_join(pred_orfs, annotations, by = "id") %>%
  arrange(start) %>%
  mutate(feature_key = "CDS") %>%
  mutate(qualifier_key = "product") %>%
  select(start, stop, feature_key, qualifier_key, description)

write.table(annotation_table, col.names = FALSE, sep = '\t', 
            file = paste0(ann_path, "/", orf_name, "_annotation_table.tsv"),
            row.names = FALSE, quote = FALSE)

