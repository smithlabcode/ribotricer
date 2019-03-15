class Interval:
    """Class for interval
       All the intervals used in this project is 1-based and closed
    """

    def __init__(self, chrom=None, start=1, end=1, strand=None):
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
