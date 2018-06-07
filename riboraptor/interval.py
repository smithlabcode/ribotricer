class Interval:
    """Class for bed interval"""

    def __init__(self, chrom=None, start=0, end=0, strand=None):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand

    def __eq__(self, other):
        """Override the default Equals behavior"""
        return (self.chrom == other.chrom and self.start == other.start
                and self.end == other.end and self.strand == other.strand)

    def __repr__(self):
        return '{}\t{}\t{}\t{}'.format(self.chrom, self.start, self.end,
                                       self.strand)
