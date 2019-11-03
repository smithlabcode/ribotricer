#/bin/bash
set -eox pipefail
#FASTA=/home/cmb-panasas2/skchoudh/genomes/R64-1-1/fasta/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
#GTF=/home/cmb-panasas2/skchoudh/genomes/R64-1-1/annotation/Saccharomyces_cerevisiae.R64-1-1.96.gtf
#START_CODONS=ATG,AAG,ACG,ATC,GTG,AGG,ATA,ATT,CTG,TTG 
wget -c https://www.dropbox.com/s/vztvkq45k8gdv5k/ribotricer_test_data_R64.zip
unzip ribotricer_test_data_R64.zip
#ribotricer prepare-orfs --gtf $GTF --prefix ribotricer_test_data_R64/R64_index_out --fasta $FASTA --start_codons $START_CODONS --longest
ribotricer detect-orfs --bam ribotricer_test_data_R64/bams_unique/SRP028552/ribo_SRX332185.bam --ribotricer_index ribotricer_test_data_R64/index/ribotricer_v96_annotation_longest_candidate_orfs.tsv --prefix ribotricer_test_data_R64/SRP028552_out/ribo
#ribotricer detect-orfs --bam ribotricer_test_data_R64/bams_unique/SRP028552/rna_SRX332188.bam --ribotricer_index ribotricer_test_data_R64/index/ribotricer_v96_annotation_longest_candidate_orfs.tsv --prefix ribotricer_test_data_R64/SRP028552_out/rna
ribotricer learn-cutoff --ribo_bams ribotricer_test_data_R64/bams_unique/SRP028552/ribo_SRX332185.bam,ribotricer_test_data_R64/bams_unique/SRP075766/ribo_SRX1801603.bam --rna_bams ribotricer_test_data_R64/bams_unique/SRP028552/rna_SRX332188.bam,ribotricer_test_data_R64/bams_unique/SRP075766/rna_SRX1801650.bam --prefix ribotricer_test_data_R64/learn_cutoff_2_datasets --ribotricer_index ribotricer_test_data_R64/index/ribotricer_v96_annotation_longest_candidate_orfs.tsv
ribotricer learn-cutoff --ribo_tsvs ribotricer_test_data_R64/learn_cutoff_2_datasets__ribo_bam_1_translating_ORFs.tsv,ribotricer_test_data_R64/learn_cutoff_2_datasets__ribo_bam_2_translating_ORFs.tsv --rna_tsvs ribotricer_test_data_R64/learn_cutoff_2_datasets__rna_bam_1_translating_ORFs.tsv,ribotricer_test_data_R64/learn_cutoff_2_datasets__rna_bam_2_translating_ORFs.tsv 

