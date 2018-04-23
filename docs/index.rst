Welcome to riboraptor's documentation!
======================================

``riboraptor`` library is a Python library for analysis of Ribo-seq data.
It assumes the reads have already been aligned to a reference and are available as BAM/SAM for input.

Current capabilities of ``riboraptor`` include:

:Visualization:
    - Read length distribution
    - Metagene coverage
    - Possible pausing sites

:Utilities:
    - Mapping summary statistics
    - 5'UTR/CDS/3'UTR coverage
    - P-site offset calculation

.. code-block:: console

    $ Usage: riboraptor [OPTIONS] COMMAND [ARGS]...                              

      riboraptor: Tool for ribosome profiling analysis                         

     Options:                             
        --version  Show the version and exit.                                    
        --help     Show this message and exit.                                   

     Commands:                            
        bam-to-bedgraph              Convert bam to bedgraph                     
        bedgraph-to-bigwig           Convert bedgraph to bigwig                  
        collapse-gene-coverage       Collapse gene coverage to metagene of target...                                                                        
        count-all-features           Count reads in 5'UTR/CDS/3'UTR regions      
        count-in-feature             Count reads in given feature bed file       
        diff-region-enrichment       calculate enrichment of cds over utr3/utr5...                                                                          
        export-bed-fasta             Export gene level fasta from specified bed...                                                                          
        export-complete-fasta        Export gene level fasta from specified bed...                                                                          
        export-gene-coverages        Export gene coverages                       
        export-single-gene-coverage  Export coverage for a gene                  
        extract-star-logs            collpase star logs to a dataframe           
        gene-coverage                Calculate coverage across a gene            
        htseq-to-tpm                 Convert HTSeq counts file to TPMs sorted... 
        mapping-summary              Mapping summary                             
        metagene-coverage            Calculate metagene coverage                 
        periodicity                  Calculate periodicity                       
        plot-framewise-counts        Plot read counts highlighting frames        
        plot-read-counts             Plot read counts distribution across a gene 
        plot-read-dist               Plot read length distribution               
        read-enrichment              Calculate read length enrichment            
        read-length-dist             Calculate read length distribution          
        uniq-mapping-count           Count number of unique mapping reads     


Contents:

.. toctree::
   :maxdepth: 2

   installation
   cmd-manual
   api-usage
   scores
   modules
   history
   contributing
   authors



===============================================================================

Site map
========
 - :ref:`genindex`
 - :ref:`modindex`
 - :ref:`search`

