"""Generate sequences for ribotricer annotation"""

import pandas as pd
import pyfaidx

from tqdm import tqdm


from .interval import Interval
from .fasta import FastaReader


def orf_seq(ribotricer_index, genome_fasta, saveto):
    """Generate sequence for ribotricer annotation.
  
  Parameters
  -----------

  ribotricer_index: string
                         Path to ribotricer generate annotation 
  genome_Fasta: string
                Path to genome fasta

  saveto: string
          Path to output
  """
    fasta = FastaReader(genome_fasta)
    annotation_df = pd.read_csv(ribotricer_index, sep="\t")
    with open(saveto, "w") as fh:
        fh.write("ORF_ID\tsequence\n")
        for idx, row in tqdm(annotation_df.iterrows(), total=annotation_df.shape[0]):
            chrom = str(row.chrom)
            orf_id = row.ORF_ID
            coordinates = row.coordinate.split(",")
            strand = row.strand
            intervals = []
            seq = ""
            for coordinate in coordinates:
                start, stop = coordinate.split("-")
                start = int(start)
                stop = int(stop)
                interval = Interval(chrom, start, stop, strand)
                intervals.append(interval)

            seq = ("").join(fasta.query(intervals))
            if strand == "-":
                seq = fasta.reverse_complement(seq)
            fh.write("{}\t{}\n".format(orf_id, seq))
