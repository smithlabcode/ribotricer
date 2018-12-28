## Robustly detecting actively translating ORFs from Ribo-seq data

### Installation
We highly recommend that you install RiboCop via conda
```bash
conda install -c bioconda ribocop
```

### Workflow of RiboCop
In order to run RiboCop, you need to have three input files prepared including:
* **genome annotation file** in GTF format, supporting both GENCODE and Ensembl annotation
* **reference genome file** in FASTA format
* **alignment file** in BAM format
