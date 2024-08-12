library(kebabs)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

# test if there are two arguments: if not, return an error
if (length(args) != 2) {
  stop("Please supply fasta file and kmer length", call. = FALSE)
}

if (file.exists(args[1])) {
  fasta.path <- args[1]
} else {
  stop(paste("file", args[1], "doesn't exist"))
}

if (!is.na(as.numeric(args[2]))) {
  k.len <- as.numeric(args[2])
} else {
  stop(paste(args[2], "isn't a valid kmer length"))
}

# fasta.path = file.path("fasta_examples", "cluster_all_filtered_longest.fa")
# k.len = 8

all.obj <- readDNAStringSet(fasta.path)
skernel <- spectrumKernel(k=k.len, normalized=FALSE)
kmerFreq <- drop(getExRep(all.obj, skernel))

kmerFreqCounts <- as(kmerFreq, "TsparseMatrix")
long_format <- data.frame(ids = rownames(kmerFreqCounts)[kmerFreqCounts@i + 1], 
                          kmer = colnames(kmerFreqCounts)[kmerFreqCounts@j + 1], 
                          counts = kmerFreqCounts@x)
counts_df <- long_format %>%
  separate_wider_delim(ids, 
                       delim = "|", 
                       names = c('ensembl_gene_id', 'ensembl_gene_id_uniq', 'ensembl_transcript_id', 'ensembl_transcript_id_uniq', 'gene_name')
                       ) %>%
  select(-ensembl_gene_id_uniq, -ensembl_transcript_id_uniq)

all_genes <- counts_df %>%
  select(-kmer, -counts) %>%
  distinct() %>%
  mutate(kmer = strrep("N", k.len),
         counts = 1)

with_group <- counts_df %>%
  mutate(
    grouping = substr(kmer, 1, nchar(kmer) - 6))
df_list <- group_split(with_group %>% group_by(grouping))

for (ele in 1:length(df_list)) {
  outpath <- gsub('\\.\\w+$', paste0('_', k.len, "mer_p", ele, ".rds"), fasta.path)
  print(outpath)
  saveRDS(df_list[[ele]] %>% select(-grouping) %>% ungroup() %>% rbind(all_genes), outpath)
}


# outpath <- gsub('\\.\\w+$', paste0('_', k.len, "_kmer.rds"), fasta.path)
# saveRDS(counts_df, outpath)

