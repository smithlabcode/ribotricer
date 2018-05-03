# TODO

## riboraptor

- [ ] Add Snakemake and config template
- [ ] Add tests with demo sra files
- [ ] Make conda requirement strict; Add more detailed installation docs
- [ ] Add annotation BED files for each genome's GTF and associated notebook
- [ ] Improve documentation:
    - [ ] QC - Fragment length, Periodicity, Summary of STAR mapping results
    - [ ] CDS/(UTR5+CDS+UTR3) distribution
    - [ ] Instructions for parsing gene_coverages.tsv file in both Python and R
- [ ] Add differential translation analysis with Riborex
- [ ] Add ribotapor for uORF or dORF identification

## ribopod


- [ ] Genome browser track links organized by each project (or publication bibtex key)
- [ ] Visualization of fragment length distribution
- [ ] Visualization of metagene plots
- [ ] Visualization of genewise coverage
- [ ] Visualization of CDS/UTR enrichment across genes
- [ ] MultiQC report 

## cleanup

- Commands to keep;
  - [ ] export-gene-coverage (Rename?)
  - [ ] metagene-coverage 
  - [ ] read-length-dist (Also change pickle => tsv)
  - [ ] periodicity 
  - [ ] plot-read-dist (Rename?)
  - [ ] plot-read-counts (Rename to plot-metagene)
  - [ ] export-bed-fasta 
  - [ ] create-bed-region
 
## comments
- [ ] remove /tmp
- [ ] make snakemake work locally
- [ ] make better doc for usage mode


## From group meeting (04/27/2018)

- [ ] Emulate databases with gene name support
- [ ] For template use one empty line and one commented version of the line
- [ ] Change /tmp paths


# Features missing

- [ ] Metagene plots should be separated for fragment length
- [ ] Add a module for guessing strand protocol of the library
- [ ] Examples where strand information is essential
- [ ] Fragment length wise distribution of 3'UTR/CDS/5'UTR => Which fragments are truly periodic?

# Install instructions

- [ ] Add instructions to output logs in OUT_DIR
- [ ] Add insttuctions about slurm/qsub change
- [ ] riboraptor_output.html should be reanmed?
- [ ] Pipeline requires storage for sqlite
