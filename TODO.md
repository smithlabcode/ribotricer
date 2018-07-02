# TODO

## riboraptor

- [x] Add Snakemake and config template
- [x] Add tests - S
- [x] Make conda requirement strict; Add more detailed installation docs
- [x] Add annotation BED files for each genome's GTF and associated notebook
- [ ] Improve documentation:
    - [ ] QC - Fragment length, Periodicity, Summary of STAR mapping results
    - [ ] CDS/(UTR5+CDS+UTR3) distribution
    - [ ] Instructions for parsing gene_coverages.tsv file in both Python and R
- [ ] Add differential translation analysis with Riborex
- [ ] Add translating ORF detection
- [ ] Add RSEM to Snakemake, make use of mapping weights to use multimapping reads? - S

## ribopod


- [ ] Genome browser track links organized by each project (or publication bibtex key) - S 
- [ ] Visualization of fragment length distribution - S
- [ ] Visualization of metagene plots - S
- [ ] Visualization of genewise coverage - S
- [ ] Visualization of CDS/UTR enrichment across genes - S
- [ ] MultiQC report - S

## cleanup

- Commands to keep;
  - [x] export-gene-coverage (Rename?)
  - [x] metagene-coverage 
  - [x] read-length-dist (Also change pickle => tsv)
  - [x] periodicity 
  - [x] plot-read-dist (Rename?)
  - [x] plot-read-counts (Rename to plot-metagene)
  - [x] export-bed-fasta 
  - [x] create-bed-region
 
## comments
- [x] remove /tmp 
- [x] make snakemake work locally - S
- [ ] make better doc for usage mode - S


## From group meeting (04/27/2018)

- [ ] Emulate databases with gene name support
- [x] For template use one empty line and one commented version of the line
- [x] Change /tmp paths


# Features missing
- [ ] Compare with [plastid](https://plastid.readthedocs.io/en/latest)
- [ ] Metagene plots should be separated for fragment length - S
- [ ] Add a module for guessing strand protocol of the library - S
- [ ] Examples where strand information is essential - S
- [ ] Fragment length wise distribution of 3'UTR/CDS/5'UTR => Which fragments are truly periodic? - S

# Install instructions

- [ ] Add instructions to output logs in OUT_DIR
- [ ] Add instructions about slurm/qsub change
- [ ] riboraptor_output.html should be renamed?
- [ ] Pipeline requires storage for sqlite
