rm(list = ls())

library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)

load_data <- function(path, gene_list) {
  all_kmers <- readRDS(path) %>%
    filter(ensembl_gene_id %in% gene_list) %>%
    mutate(ingroup = counts > 0) %>%
    select(-counts, -ensembl_transcript_id) %>%
    complete(nesting(ensembl_gene_id, gene_name), kmer, fill = list(ingroup = 0))

  if (!'ensembl_gene_id' %in% colnames(all_kmers)) stop("Missing ensembl_gene_id column")
  if (nrow(all_kmers) == 0) stop("No data for gene list")

  k = nchar(all_kmers$kmer[[1]])
  print(paste0(Sys.time(), ", k=", k, ", ", all_kmers %>% select(kmer) %>% distinct() %>% nrow()))
  return(all_kmers)
}

run_hg_tests <- function(all_kmers, params, cls_to_test) {
    hg_res <- all_kmers %>%
        left_join(params) %>%
        group_by(kmer) %>%
        summarise(
          sample.success = sum((cluster == cls_to_test)[ingroup]),
          population.success = sum(ingroup),
          population.failure = sum(!ingroup),
          sample.size = sum(cluster == cls_to_test),
          population.size = n()
        ) %>%
        mutate(
          fold.enrichment = (sample.success / sample.size) / (population.success / population.size),
          p.value.more = phyper( # observed more than sample.success
            q = sample.success,
            m = population.success,
            n = population.failure,
            k = sample.size,
            lower.tail = FALSE
          ),
          p.value.less = phyper( # observed less than sample.success
            q = sample.success - 1,
            m = population.success,
            n = population.failure,
            k = sample.size,
            lower.tail = TRUE
          ),
          p.value.eq = dhyper( # observed sample.success
            x = sample.success,
            m = population.success,
            n = population.failure,
            k = sample.size
          ),
          cluster = cls_to_test
        )
  return(hg_res)
}

get_hg_res <- function(cls_to_test, path, params, gene_list) {
  gc()
  data <- load_data(path, gene_list)
  gc()
  return(run_hg_tests(data, params, cls_to_test))
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Please supply path kmer table file, parameter table, output path and run id", call. = FALSE)
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

if (file.exists(dirname(args[3]))) {
  out.base <- args[3]
} else {
  stop(paste("folder", args[3], "doesn't exist"))
}

run.id <- args[4]
out.path <- file.path(out.base, paste0(run.id, '_', gsub(".rds$", ".tmpkmers", basename(kmer.file))))

# kmer.file <- '..\\kmer_analysis\\kebabs_pipeline\\files\\cluster_all_filtered_longest_5mer_p1.rds'
# params.path <- "data_files/clustering/cluster_assignments_rowmeans.csv"

params_df <- fread(params.path)
if (!'ensembl_gene_id' %in% colnames(params_df)) stop("Missing ensembl_gene_id column")
genes_to_analyze <- (params_df %>% select(-cluster) %>% distinct())$ensembl_gene_id #TODO enable gene_name

t0 <- Sys.time()
hg_list <- lapply(unique(params_df$cluster), get_hg_res, path = kmer.file, params = params_df, gene_list = genes_to_analyze)
t1 <- Sys.time()

all_hg <- rbindlist(hg_list)

all_hg_plot <- all_hg %>%
  mutate(p.value = ifelse(fold.enrichment < 1, p.value.less + p.value.eq, p.value.more + p.value.eq))

fwrite(all_hg_plot, out.path, sep = '\t')
