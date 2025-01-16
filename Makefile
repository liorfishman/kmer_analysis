
# ----------------------------------------------------------
# make sure to update all following parameters
# datatypes tested are - 3utr 5utr coding
# make sure to note in the file name what ensembl version was downloaded
# ----------------------------------------------------------


directory = /path/to/files
ORGANISM = drerio_gene_ensembl
dtype = 5utr

FASTA_FILE = biomart_$(ORGANISM)_$(dtype)_ensembl103.fasta


# ----------------------------------------------------------
# Check R requirements
# (1) make check_requirements
# ----------------------------------------------------------

check_requirements:
	source $$MODULESHOME/init/bash; module load R4/4.1.3; \
	Rscript check_requirements.R

# ----------------------------------------------------------
# Prepare fastas
# (1 - 1) make download_fasta
# (1 - 2) if using fasta - make with_existing_fasta
# (2) make create_fastas_all_job dtype=3utr
# ----------------------------------------------------------

download_fasta:
	mkdir $(dtype); \
	wget -O $(dtype)/$(FASTA_FILE) 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "$(ORGANISM)" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "ensembl_gene_id_version" /><Attribute name = "ensembl_transcript_id" /><Attribute name = "ensembl_transcript_id_version" /><Attribute name = "$(dtype)" /><Attribute name = "external_gene_name" /></Dataset></Query>'

download_fasta_job:
	sbatch.pl -R 8 -t 1:0:0 "make download_fasta dtype=$(dtype);";

with_existing_fasta:
	mkdir $(dtype); \
	echo "for analysis to work, add fasta to $(dtype) directory"

create_fastas_all:
	cd $(dtype); \
	mkdir fastas; \
	Rscript ../create_fastas.R $(FASTA_FILE); \
	cd ..

create_fastas_all_job:
	sbatch.pl -R 8 -t 1:0:0 "module load hurcs; module load R4/4.1.3; make create_fastas_all dtype=$(dtype);";

# ----------------------------------------------------------
# kmer analysis
# (1) make create_kmer_folders dtype=3utr filter=_longest
# (2) make create_kmer_data_all_kebabs
# (3) make all_ks_tests_kebabs category=poly
# ----------------------------------------------------------

# four columns gene_name/ensembl_id, param_name, param_val and category
PARAMETERS_FILE = parameters.txt

#changes according to file, could be multiple
CATEGORY = poly

create_kmer_folders:
	mkdir $(dtype)/kmer_matrices_kebabs;
	mkdir $(dtype)/kmer_out_kebabs;
	ln -s $(directory)/$(dtype)/fastas/cluster_all_filtered$(filter).fa $(directory)/$(dtype)/kmer_matrices_kebabs/cluster_all_filtered$(filter).fa;

create_kmer_data_kebabs:
	Rscript create_kmer_matrices_k.R $(dtype)/kmer_matrices_kebabs/cluster_all_filtered$(filter).fa $(klen)

create_kmer_data_all_kebabs:
	rm -rf kmer_table.txt; \
	$(foreach k, $(shell seq 4 10), \
		echo "make create_kmer_data_kebabs klen=$(k) dtype=$(dtype) filter=$(filter)" >> kmer_table.txt; \
	) \
	sbatch.pl -a kmer_table.txt -R 16 -t 1:00:0 "module load hurcs; module load R4/4.1.3";

run_ks_test_kebabs:
	rm -rf test_table_$(term)_$(n)$(p).txt; \
	$(foreach k, $(shell ls $(dir)/*$(n)mer_p$(p)*.rds), \
		echo "Rscript motif_continuous_tests_p.R $(k) $(PARAMETERS_FILE) $(dtype)/kmer_out_kebabs $(term)" >> test_table_$(term)_$(n)$(p).txt; \
	) \
	sbatch.pl -a test_table_$(term)_$(n)$(p).txt -R 32 -t 1:00:0 "module load hurcs; module load R4/4.1.3";

all_ks_tests_kebabs:
	make run_ks_test_kebabs dir=$(dtype)/kmer_matrices_kebabs n=4 term=$(category)
	make run_ks_test_kebabs dir=$(dtype)/kmer_matrices_kebabs n=5 term=$(category)
	make run_ks_test_kebabs dir=$(dtype)/kmer_matrices_kebabs n=6 term=$(category)
	make run_ks_test_kebabs dir=$(dtype)/kmer_matrices_kebabs n=7 term=$(category)
	make run_ks_test_kebabs dir=$(dtype)/kmer_matrices_kebabs n=8 term=$(category)
	make run_ks_test_kebabs dir=$(dtype)/kmer_matrices_kebabs n=9 term=$(category)
	make run_ks_test_kebabs dir=$(dtype)/kmer_matrices_kebabs n=10 term=$(category)

run_ks_test_kebabs_comp:
	mkdir $(dtype)/kmer_out_kebabs/WT1KO; \
	mkdir $(dtype)/kmer_out_kebabs/WT2KO; \
	mkdir $(dtype)/kmer_out_kebabs/1KO2KO; \
	rm -rf comp_test_table_$(term)_$(n)$(p).txt; \
	$(foreach k, $(shell ls $(dir)/*$(n)mer_p$(p)*.rds), \
		echo "Rscript motif_continuous_tests_comp.R $(k) MS_Riboseq_TE_params_WT.txt MS_Riboseq_TE_params_1KO.txt $(dtype)/kmer_out_kebabs/WT1KO $(term)" >> comp_test_table_$(term)_$(n)$(p).txt; \
		echo "Rscript motif_continuous_tests_comp.R $(k) MS_Riboseq_TE_params_WT.txt MS_Riboseq_TE_params_2KO.txt $(dtype)/kmer_out_kebabs/WT2KO $(term)" >> comp_test_table_$(term)_$(n)$(p).txt; \
		echo "Rscript motif_continuous_tests_comp.R $(k) MS_Riboseq_TE_params_1KO.txt MS_Riboseq_TE_params_2KO.txt $(dtype)/kmer_out_kebabs/1KO2KO $(term)" >> comp_test_table_$(term)_$(n)$(p).txt; \
	) \
	sbatch.pl -a comp_test_table_$(term)_$(n)$(p).txt -R 16 -t 1:00:0 "module load hurcs; module load R4/4.1.3";

all_ks_tests_kebabs_comp:
	make run_ks_test_kebabs_comp dir=$(dtype)/kmer_matrices_kebabs n=4 term=$(category)
	make run_ks_test_kebabs_comp dir=$(dtype)/kmer_matrices_kebabs n=5 term=$(category)
	make run_ks_test_kebabs_comp dir=$(dtype)/kmer_matrices_kebabs n=6 term=$(category)
	make run_ks_test_kebabs_comp dir=$(dtype)/kmer_matrices_kebabs n=7 term=$(category)
	make run_ks_test_kebabs_comp dir=$(dtype)/kmer_matrices_kebabs n=8 term=$(category)
	make run_ks_test_kebabs_comp dir=$(dtype)/kmer_matrices_kebabs n=9 term=$(category)
	make run_ks_test_kebabs_comp dir=$(dtype)/kmer_matrices_kebabs n=10 term=$(category)

run_hg_test_kebabs:
	rm -rf test_table_$(id)_$(n)$(p).txt; \
	$(foreach k, $(shell ls $(dir)/*$(n)mer_p$(p)*.rds), \
		echo "Rscript motif_hg_test.R $(k) $(param_file) $(dtype)/kmer_out_kebabs_hg $(id)" >> test_table_$(id)_$(n)$(p).txt; \
	) \
	sbatch.pl -a test_table_$(id)_$(n)$(p).txt -R 8 -t 1:00:0 "module load R4";

all_hg_tests_kebabs:
	make run_hg_test_kebabs dir=$(dtype)/kmer_matrices_kebabs n=4 param_file=$(param_file) id=$(run_id)
	make run_hg_test_kebabs dir=$(dtype)/kmer_matrices_kebabs n=5 param_file=$(param_file) id=$(run_id)
	make run_hg_test_kebabs dir=$(dtype)/kmer_matrices_kebabs n=6 param_file=$(param_file) id=$(run_id)
	make run_hg_test_kebabs dir=$(dtype)/kmer_matrices_kebabs n=7 param_file=$(param_file) id=$(run_id)
	make run_hg_test_kebabs dir=$(dtype)/kmer_matrices_kebabs n=8 param_file=$(param_file) id=$(run_id)
	make run_hg_test_kebabs dir=$(dtype)/kmer_matrices_kebabs n=9 param_file=$(param_file) id=$(run_id)
	make run_hg_test_kebabs dir=$(dtype)/kmer_matrices_kebabs n=10 param_file=$(param_file) id=$(run_id)

# ----------------------------------------------------------
# kmer analysis - quality
# (1) make check_all_files category=$(CATEGORY) - and follow instructions
# (2) make combine_all_outputs
# ----------------------------------------------------------


check_all_files:
	Rscript check_missing_files.R $(dtype)/kmer_out_kebabs $(category)

combine_all_outputs:
	sbatch.pl "awk '(NR == 1) || (FNR > 1)' $(dtype)/kmer_out_kebabs/$(category)*.tmpkmers > $(dtype)/kmer_out_kebabs/$(category)_ks_raw_with_stats.tsv; rm -rf $(dtype)/kmer_out_kebabs/$(category)*.tmpkmers"

check_all_files_comp:
	Rscript check_missing_files.R $(dtype)/kmer_out_kebabs/$(set) $(category)

combine_all_outputs_comp:
	sbatch.pl "awk '(NR == 1) || (FNR > 1)' $(dtype)/kmer_out_kebabs/$(set)/$(category)*.tmpkmers > $(dtype)/kmer_out_kebabs/$(set)_$(category)_ks_raw_with_stats.tsv; rm -rf $(dtype)/kmer_out_kebabs/$(set)/$(category)*.tmpkmers"

check_all_files_hg:
	Rscript check_missing_files.R $(dtype)/kmer_out_kebabs_hg $(run_id)

combine_all_outputs_hg:
	sbatch.pl "awk '(NR == 1) || (FNR > 1)' $(dtype)/kmer_out_kebabs_hg/$(run_id)*.tmpkmers > $(dtype)/kmer_out_kebabs_hg/$(run_id)_ks_raw_with_stats.tsv; rm -rf $(dtype)/kmer_out_kebabs_hg/$(run_id)*.tmpkmers"


