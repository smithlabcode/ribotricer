#/bin/bash
set -eox pipefail
wget -c https://www.dropbox.com/s/8iynb0k5dbhv0b1/ribotricer_test_data_tair10.zip
unzip ribotricer_test_data_tair10.zip
ribotricer detect-orfs --bam ribotricer_test_data_tair10/bams_unique/SRX219170.bam --ribotricer_index ribotricer_test_data_tair10/index/ribotricer_v44_annotation_longest_candidate_orfs.tsv --prefix ribotricer_test_data_tair10/SRX219170_generated
MD5_expected=$(md5sum ribotricer_test_data_tair10/translating_ORFs/SRX219170_translating_ORFs.tsv | awk '{ print $1 }')
MD5_observed=$(md5sum ribotricer_test_data_tair10/SRX219170_generated_translating_ORFs.tsv | awk '{ print $1 }')
if [ $MD5_expected != $MD5_observed ]; then
echo "Detect ORFs MD5 mismatch"
fi
ribotricer count-orfs --detected_orfs ribotricer_test_data_tair10/SRX219170_generated_translating_ORFs.tsv --ribotricer_index ribotricer_test_data_tair10/index/ribotricer_v44_annotation_longest_candidate_orfs.tsv --features annotated --out ribotricer_test_data_tair10/SRX219170_generated_annotated_counts.tsv
MD5_expected=$(md5sum ribotricer_test_data_tair10/orfs_count/SRX219170_annotated_counts_cnt.txt | awk '{ print $1 }')
MD5_observed=$(md5sum ribotricer_test_data_tair10/SRX219170_generated_annotated_counts.tsv | awk '{ print $1 }')
if [ $MD5_expected != $MD5_observed ]; then
echo "Detect ORFs MD5 mismatch"
fi
