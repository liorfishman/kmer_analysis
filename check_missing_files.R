args <- commandArgs(trailingOnly = TRUE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Please supply path to folder with kmer tables, parameter table and output path", call. = FALSE)
}

if (file.exists(args[1])) {
  out.base <- args[1]
} else {
  stop(paste("folder", args[1], "doesn't exist"))
}

if (length(args) == 2) {
  filter.term <- args[2]
}

#### code ####

output_files <- list.files(out.base)

all_ks <- sort(as.numeric(unique(gsub(paste0(filter.term, "_cluster_all_filtered_longest_(\\d+)mer_p(\\d+).tmpkmers"), "\\1", output_files))))

print("k lengths in analysis:")
print(all_ks)

all_output_files <- c()
for (k in all_ks) {
  all_output_files <- c(all_output_files, paste0(filter.term, "_cluster_all_filtered_longest_", k, "mer_p", seq(1, (4^k)/4096),".tmpkmers"))
}

missing_files <- setdiff(all_output_files, output_files)

if (length(missing_files) == 0) {
  print("ran succesfully, continue analysis")
} else {
  print("missing the following files:")
  print(missing_files)
}
