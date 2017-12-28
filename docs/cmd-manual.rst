============================================================================
The Smithlab Ribo-seq Data Analysis Pipeline (RiboPipe)
============================================================================


Contents
========
The RiboPipe software package is a comprehensive pipeline and set of tools for
analyzing Ribo-seq sequencing data. This manual explains
the stages in our pipeline, how to use the analysis tools, and how to modify
the pipeline for your specific context.


Assumptions
===========
Our pipeline was designed to run in a cluster computing context, with many
processing nodes available, and a job submission system like PBS or SGE.
Much of this analysis is computationally intensive. We assume that individual
nodes will have several GB of memory available for processing. Typically the
data we deal with amounts to a minimum of 100GB for a mammalian methylome
at 10x coverage. Intermediate files may cause this amount to more than double
during execution of the pipeline, and likely at the end of the pipeline
the total size of files will amount to almost double the size of the raw data.
Users are assumed to be quite familiar with UNIX/Linux and related concepts
(e.g. building software from source,
using the command line, shell environment variables, etc.).


Translation profile construction
================================
Ribo-seq experiments are always single-end sequenced. Ribosome protected fragments range
from 28-31 nucleotides and hence most experiments involve 50 bp single end reads. Before mapping,
we need to get rid of the adapters which are ligated at the 3' end of the fragments as part of the
library protocol.

Trimming Reads
--------------

We use trim_galore_ for trimming. It automates adapter trimming:


.. code-block:: sh

   trim_galore -o <out_dir> -q <min_quality> <input.fq.gz>

-o out_dir        Output directory
-q min_quality    Trim low-quality ends from reads in addition to adapter removal

Mapping Reads
-------------

We use STAR_ to map reads. The first step is to create an index, preferably
using a GTF file. If the index step is run without a GTF file (which is optional), 
STAR_ will not be splice-aware.


Creating Index
~~~~~~~~~~~~~~
.. code-block:: sh

   STAR --runThreadN <threads>\
        --runMode genomeGenerate\
        --genomeDir <index_out_dir>\
        --genomeSAindexNbases <SA_INDEX_Nbases>\
        --genomeFastaFiles <input.fasta>\
        --sjdbGTFfile <input.gtf>


--runThreadN threads                     Number of threads to use
--runMode genomeGenerate                 Flag to set for index mode
--genomeDir index_out_dir                Directory to write index files to
--genomeSAindexNbases SA_INDEX_Nbases    min(14, log2(GenomeLength)/2 - 1),
                                         this **must** be scaled down for
                                         small genomes
--genomeFastaFiles input_fasta           Path to reference fasta
--sjdbGTFfile input_gtf                  Path to GTf file


Mapping
~~~~~~~

.. code-block:: sh

    STAR --runThreadN <threads>\
         --genomeDir <input.index>\
         --outFilterMismatchNmax 2\
         --alignIntronMin <ALIGN_INTRON_Nmin>\
         --alignIntronMax <ALIGN_INTRON_Nmax>\
         --outFileNamePrefix <params.prefix> --readFilesIn <input.R1>\
         --outSAMtype BAM Unsorted\
         --readFilesCommand zcat\
         --quantMode TranscriptomeSAM\
         --outTmpDir /tmp/<params.name>_tmp\
         --outReadsUnmapped Fastx\


--runThreadN threads                  Number of threads to use
--genomeDir index_out_dir             Path to index directory
--outFilterMismatchNmax mismatches    Allow a maximum of mismatches=2
--alignIntronMin ALIGN_INTRON_Nmin    Minimum intron size. Any genomic gap
                                      is considered intron if its
                                      length >= alignIntronMin.
--alignIntronMax ALIGN_INTRON_Nmax    Maximum intron size
--outFileNamePrefix prefix            Prefix for output files
--readFilesIn input_fq_gz             Path to input fastq.gz
--outSAMtype outtype                  Output an unsorted BAM file (outtype=BAM Unsorted)
--readFilesCommand zcat               Since input is gzipped use zcat to
                                      decompress it on the fly
--quantMode TranscriptomeSAM          Also output BAM aligned to the
                                      transcriptome
--outTmpDir tpmdir                    Directory to use for writing 
                                      temporary files
--outReadsUnmapped Fastx              Write unmapped reads to separate 
                                      fastq file


Sorting and Indexing
~~~~~~~~~~~~~~~~~~~~

STAR outputted BAM files are not sorted. We need a BAM file sorted
by coordinates.

.. code-block:: sh

   samtools sort <prefix>Aligned.out.bam -o <output.bam> -T <tmpdir>_sort &&\
   samtools index <prefix>Aligned.out.bam

Additionaly, we also need BAM file sorted by name, since htseq-counts_
(and featureCounts_) prefer a BAM sorted by name in their default mode.

.. code-block:: sh

    samtools sort -on <input.bam> -T <tmpdir> -o <output.bam> &&\
    samtools index <output.bam>


Translation profile analysis
============================

Once we have the bams, we are ready for downstream analysis. We will use our riboraptor_ tool 
for all downstream analysis.

The first step is to simply caculate number of uniquely mapping reads.
We recommend a minimum of 5 million reads for any downstream analysis.

.. code-block:: sh

    riboraptor uniq-mapping-count --bam <input.bam>

--bam input.bam    Path to bam file


Example
-------
We will use two samples from GSE94454_ as examples for examples that follow.

.. code-block:: console

   $ riboraptor uniq-mapping-count --bam data/SRR5227310.bam
   28637667
   $


This is a pretty deep library. Next, we check out what the fragment distribution looks like:

.. code-block:: console

   $ riboraptor read-length-dist --bam data/SRR5227310.bam | riboraptor plot-read-dist --saveto SRR5227310.png

.. figure:: images/SRR5227310.png
    :align: center
    :alt: Fragment length distribution SRR5227310
    :figclass: align center

    Fragment length distribution for SRR5227310

An ideal Ribo-seq library is expected to have 28-31 nt long fragments most enriched.
We can calculate enrichment and plot the fragment size distribution using riboraptor:

.. code-block:: console

   $ riboraptor read-length-dist --bam data/SRR5227310.bam | riboraptor read-enrichment


So the fragment length distribution doesn't seem to be enriched. We next perform metagene 
analysis. Ribo-seq data is expected to have an inherent periodicity of 3, since ribosomes move
one codon at a time during active translation.

.. code-block:: console

   $ riboraptor metagene 

This is not likely a Ribo-seq sample.

Let's try another sample: SRR5227306.

.. code-block:: console

   $ riboraptor uniq-mapping-count --bam data/SRR5227306.bam
   10658208

.. code-block:: console

   $ riboraptor read-length-dist --bam data/SRR5227306.bam | riboraptor plot-read-dist --saveto SRR5227306.png

.. figure:: images/SRR5227306.png
    :align: center
    :alt: Fragment length distribution SRR5227306
    :figclass: align center

    Fragment length distribution for SRR5227306

.. code-block:: console

   $ riboraptor read-length-dist --bam data/SRR5227306.bam | riboraptor read-enrichment


Metagene counts : Calculate Periodicity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Distributio of 5'UTR/CDS/3'UTR counts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





.. _trim_galore: https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
.. _STAR: https://github.com/alexdobin/STAR
.. _riboraptor: https://github.com/saketkc/riboraptor
.. _GSE94454: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94454
.. _htseq-counts: https://htseq.readthedocs.io/
.. _featureCounts: http://bioinf.wehi.edu.au/featureCounts/
