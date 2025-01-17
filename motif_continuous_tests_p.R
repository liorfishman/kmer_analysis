rm(list = ls())

library(dplyr)
library(data.table)
library(tidyr)

run_ks_tests <- function(all_kmers, params) {
  kmer_test <- all_kmers %>%
    left_join(params) %>%
    filter(!is.na(param_val)) %>%
    group_by(kmer, param_name) %>%
    filter(max(ingroup) == 1) %>%
    filter(min(ingroup) == 0) %>%
    summarise(
      ks_pval_gr = ks.test(param_val[ingroup], param_val[!ingroup], alternative = 'greater')$p.value,
      ks_pval_ls = ks.test(param_val[ingroup], param_val[!ingroup], alternative = 'less')$p.value,
      kmer_mean = mean(param_val[ingroup], na.rm = T),
      bg_mean = mean(param_val[!ingroup], na.rm = T),
      tot_sd = sd(param_val, na.rm = T),
      kmer_size = sum(ingroup, na.rm = T),
      bg_size = sum(!ingroup, na.rm = T)
    )
  
  return(kmer_test)
  
}

get_ks_res <- function(path, params, gene_list) {
  gc()
  data <- load_data(path, gene_list)
  gc()
  return(run_ks_tests(data, params))
}

load_data <- function(path, gene_list) {
  all_kmers <- readRDS(path) %>%
    filter(ensembl_gene_id %in% gene_list) %>%
    mutate(ingroup = counts > 0) %>%
    select(-counts, -ensembl_transcript_id) %>%
    complete(nesting(ensembl_gene_id, gene_name), kmer, fill = list(ingroup = 0))
  
  k = nchar(all_kmers$kmer[[1]])
  print(paste0(Sys.time(), ", k=", k, ", ", all_kmers %>% select(kmer) %>% distinct() %>% nrow()))
  return(all_kmers)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) <= 3) {
  stop("Please supply path to folder with kmer tables, parameter table and output path", call. = FALSE)
}

if (file.exists(args[1])) {
  kmer.file <- args[1]
} else {
  stop(paste("file", args[1], "doesn't exist"))
}

if (file.exists(args[2])) {
  params.path <- args[2]
} else {
  stop(paste("file", args[2], "doesn't exist"))
}

if (file.exists(args[3])) {
  out.base <- args[3]
} else {
  stop(paste("folder", args[3], "doesn't exist"))
}

# kmer.file <- file.path("kebabs_pipeline", "files", "cluster_analysed_filtered_longest_folded_4mer_p1.rds")
# kmer.file <- file.path("kebabs_pipeline", "files", "cluster_all_filtered_longest_10mer_p20.rds")
# params.path <- file.path("kebabs_pipeline", "paramaters_for_kmer_hld.tsv")
# params.path <- file.path("dina_params_avg2.tsv")
# filter.term <- "all"
# out.base <- file.path("kebabs_pipeline", "files2")

params_df <- fread(params.path)

if (length(args) == 4) {
  filter.term <- args[4]
  
  params_df <- params_df %>%
    mutate(param_name = paste(category, param_name)) %>%
    filter(category == filter.term) %>%
    select(-category)
}

genes_to_analyze <- (params_df %>% select(-param_name, -param_val) %>% distinct())$ensembl_gene_id #TODO enable gene_name

print(paste0(length(genes_to_analyze), " genes for analysis"))
print(paste0(unique(params_df$param_name), " parameters for analysis"))

all_ks <- get_ks_res(kmer.file, params_df, genes_to_analyze)

print(paste0(Sys.time(), basename(kmer.file), " completed analysis"))

out.path <- file.path(out.base, paste0(filter.term, '_', gsub(".rds$", ".tmpkmers", basename(kmer.file))))
fwrite(all_ks, out.path, sep = '\t')

print(paste0("saved to ", out.path))