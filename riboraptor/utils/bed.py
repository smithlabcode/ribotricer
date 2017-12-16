from collections import OrderedDict
import pybedtools
import pandas as pd

def region_sizes(bed):
    """Get collapsed sizes of regions in bed indexed by name

    Parameters
    ----------
    bed : str
          Path to bed file

    region_sizes : dict
                   OrderedDict sorted in descending order by values

    """
    region_bed_grouped = pybedtools.BedTool(bed).to_dataframe().groupby('name')
    region_sizes = {}
    for gene_name, gene_group in region_bed_grouped:
        ## Get rid of trailing dots
        gene_name = re.sub(r'\.[0-9]+', '', gene_name)
        # Collect all intervals at once
        intervals = zip(gene_group['chrom'], gene_group['start'],
                        gene_group['end'], gene_group['strand'])
        for interval in intervals:
            if gene_name not in region_sizes:
                # End is always 1-based so does not require +1
                region_sizes[gene_name] = interval[2] - interval[1]
            else:
                region_sizes[gene_name] += interval[2] - interval[1]
    sizes_df = pd.DataFrame(region_sizes.items(), columns=['name', 'length'])
    sizes_df.sort_values(by='length', ascending=False, inplace=True)
    region_sizes = OrderedDict([tuple(x) for x in sizes_df[['name', 'length']].values])
    return region_sizes

def get_fasta_sequence(fasta, intervals):
    """Extract fasta sequence given a list of intervals

    Parameters
    ----------
    fasta : str
            Path to fasta file

    intervals : list(tuple)
                A list of tuple in the form [(chrom, start, stop, strand)]


    Returns
    -------
    seq : list
          List of sequences at intervals
    """
    if isinstance(fasta, str):
        fasta = Fasta(fasta)
    sequence = []
    for interval in intervals:
        chrom, start, stop, strand = interval
        if strand == '+':
            seq = fasta[chrom][int(start):int(stop)].seq
        elif strand == '-':
            seq = fasta[chrom][int(start):int(stop)].reverse.seq
        sequence.append(seq)
    return sequence

