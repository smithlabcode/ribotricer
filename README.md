# Robustly detecting actively translating ORFs from Ribo-seq data

## Installation
We highly recommend that you install RiboCop via conda
```bash
conda install -c bioconda ribocop
```

## Workflow of RiboCop
In order to run RiboCop, you need to have the following three files prepared including:
* **genome annotation file** in GTF format, supporting both GENCODE and Ensembl annotation
* **reference genome file** in FASTA format
* **alignment file** in BAM format

### Prepare candidate ORFs
The first step of RiboCop is to take the GTF file and the FASTA file to find all candidate ORFs.
In order to generate all candidate ORFs, please run
```bash
RiboCop prepare-orfs --gtf {GTF} --fasta {FASTA} --prefix {PREFIX}
```
The command above by default only includes ORFs with length longer than 60 nts,
and only uses 'ATG' as start codon. You can change the setting by including
options ```--min_orf_length``` and ```--start_codons```.

## LICENSE
RiboCop for detecting actively translating ORFs from Ribo-seq data
Copyright (C) 2018 Andrew D Smith, Wenzheng Li, Saket Choudhary and
the University of Southern California

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
