# V1.3.3 (2022-02-09)

- Print start codons when preparing orfs
- Fix for custom start codons

# v1.3.2 (2020-05-03)

- Better support for extracting sequences from non-conventional GTFs
- Support `--protein` flag for extracting protein sequences

# v1.3.1 (2019-12-12)

- Detecting uniquely mapping reads now raises a warning instead of aborting. We expect the user to ensure all reads in the bam are uniquely mapping unless it was produced for aligners that support either the NH or standard uniquely mapping flags (See common.py for a list of flags).
- Clean setup.py to improve classifiers


# v1.3.0 (2019-11-01)

- A small bug fix in metagene.py where the first annotated CDS was being skipped
- Silenced progress bars after they are complete
- More informative progress messages

# v1.2.0 (2019-10-09)

- Added additional flags for filtering:
	- `--min_valid_codons_ratio` (default=0): Minimum ratio of codons with non-zero reads to total codons for determining active translation
	- `--min_reads_per_codon` (default=0): Minimum number of reads per codon for determining active translation
	- `--min_read_density` (default=0.0): Minimum read density (total_reads/length) over an ORF total codons for determining active translation
- Added a new command-line argument to learn dataset-specific parameters from Ribo-seq and RNA-seq bams learn-cutoff
- Travis based CI testing

# v1.1.0 (2019-09-17)

- Allow specifying custom phase score cutoffs
- Added documentation for species-specific cutoffs
- Helper functions for determining best cutoff given golden datasets

# v0.1.0 (2019-05-05)

- First release of ribotricer on pypi
