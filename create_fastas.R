#### imports ####

library(seqinr)
library(dplyr)
library(stringr)
library(data.table)
library(tidyverse)

MISSING <- "sequence unavailable"
MIN_LENGTH <- 9
SHORT_SEQ_LENGTH <- 7500

write_fasta_from_df <- function (df, path) {
  for_list <- df %>%
    mutate(name = paste(ensembl_id, ensembl_id_stable, ensembl_transcript_id, ensembl_transcript_id_stable, gene_name, sep = "|")) %>%
    select(name, sequence) %>%
    column_to_rownames("name")

  list_for_fasta <- as.list(as.data.frame(t(for_list)))

  write.fasta(
    list_for_fasta,
    names = names(list_for_fasta),
    as.string = TRUE,
    file.out = path
  )
}

write_cluster_fastas <- function(df, cluster_path) {
  clusters <- fread(cluster_path)

  for (cls in unique(clusters$cluster)) {
    cur_cls <- clusters %>%
      filter(cluster == cls)

    cluster_sequences <- cur_cls %>%
      left_join(df) %>%
      select(-cluster)

    write_fasta_from_df(cluster_sequences, file.path('fastas', paste0('cluster_', cls, '_filtered_longest.fa')))
  }
}

# currently function returns longest, can change table if needed
filter_fasta <- function(fasta.path) {
  fastafile <-
    read.fasta(
      file = fasta.path,
      seqtype = "DNA",
      as.string = TRUE,
      set.attributes = FALSE
    )

  df <- data.frame(sequence = unlist(fastafile)) %>%
    rownames_to_column("name") %>%
    separate(name, sep = "\\|", into = c("ensembl_id", "ensembl_id_stable", "ensembl_transcript_id", "ensembl_transcript_id_stable", "gene_name"))

  names <- df %>%
    select(ensembl_id, ensembl_id_stable, ensembl_transcript_id, ensembl_transcript_id_stable, gene_name)

  df_no_missing <- df %>%
    filter(sequence != MISSING) %>%
    filter(!grepl('^n+$', sequence))

  df_long_enough <- df_no_missing %>%
    filter(nchar(sequence) > MIN_LENGTH)

  df_longest <- df_long_enough %>%
    mutate(seq_len = nchar(sequence)) %>%
    group_by(ensembl_id) %>%
    summarise(
     longest = max(seq_len),
     sequence = sequence[seq_len == longest],
     ensembl_transcript_id_stable = ensembl_transcript_id_stable[seq_len == longest]
    ) %>%
    distinct(ensembl_id, .keep_all = TRUE) %>% # this is very important, otherwise different transcripts of the same gene with the same sequence are included
    ungroup()

  df_longest_named <- df_longest %>%
    left_join(names)

  df_short <- df_long_enough %>%
    filter(nchar(sequence) < SHORT_SEQ_LENGTH)

  df_short_longest <- df_short %>%
    mutate(seq_len = nchar(sequence)) %>%
    group_by(ensembl_id) %>%
    summarise(
      longest = max(seq_len),
      sequence = sequence[seq_len == longest],
      ensembl_transcript_id_stable = ensembl_transcript_id_stable[seq_len == longest]
    ) %>%
    distinct(ensembl_id, .keep_all = TRUE) %>% # this is very important, otherwise different transcripts of the same gene with the same sequence are included
    ungroup()

  df_short_longest_named <- df_short_longest %>%
    left_join(names)

  return(df_longest_named)
}

#### get and validate command line arguments ####

args <- commandArgs(trailingOnly = TRUE)

# test if there are two arguments: if not, return an error
if (length(args) < 1) {
  stop("Please supply fasta file path", call. = FALSE)
}

if (file.exists(args[1])) {
  fasta.path <- args[1]
} else {
  stop(paste("file", args[1], "doesn't exist"))
}

sequence_df <- filter_fasta(fasta.path)
write_fasta_from_df(sequence_df, file.path('fastas', 'cluster_all_filtered_longest.fa'))

if (length(args) == 2) {
  if (file.exists(args[2])) {
    clusters.path <- args[2]
  } else {
    stop(paste("file", args[2], "doesn't exist"))
  }
  write_cluster_fastas(sequence_df, clusters.path)
}

# fasta.path <- file.path('motif_search_pipeline', '3utr', 'biomart_GRCz11_3UTR_ensembl103.fasta')
# clusters.path <- file.path('motif_search_pipeline', 'top_bottom_10_perc_clusters.tsv')
