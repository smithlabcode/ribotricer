import os
import warnings

from pyfaidx import Fasta


class FastaReader(object):
    """Class for reading and querying fasta file."""

    complement_letters = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}

    def __init__(self, fasta_location):
        """
        Parameters
        ---------
        fasta_location : string
                         Path to fasta file

        """
        self.fasta_location = fasta_location
        try:
            self.fasta = Fasta(
                fasta_location, as_raw=True, sequence_always_upper=True)
        except Exception as e:
            raise Exception('Error reading fasta file {} : {}'.format(
                os.path.abspath(self.fasta_location), e))

    def query(self, intervals):
        """ Query regions for sequence.

        Parameters
        ----------
        intervals: list of Interval 
                   The intervals for fasta is one-based and full-closed

        Returns
        -------
        sequences: list(str)
                   An array containing scores for each Interval
                   This function is agnostic of the strand information,
                   the position in the scores is corresponding to the interval

        .. currentmodule:: .FastaReader
        .. autosummary::
            .FastaReader

        """
        sequences = []
        chrom_lengths = self.chromosomes
        for i in intervals:
            if i.chrom not in list(chrom_lengths.keys()):
                warnings.warn(
                    'Chromosome {} does not appear in the fasta'.format(
                        i.chrom), UserWarning)
                continue

            chrom_length = chrom_lengths[i.chrom]
            if i.start > chrom_length:
                raise Exception(
                    'Chromsome start point exceeds chromosome length: {}>{}'.
                    format(i.start, chrom_length))
            elif i.end > chrom_length:
                raise Exception(
                    'Chromsome end point exceeds chromosome length: {}>{}'.
                    format(i.end, chrom_length))
            seq = self.fasta.get_seq(i.chrom, i.start, i.end)
            sequences.append(seq)
        return sequences

    def complement(self, seq):
        seq = seq.uppper()
        comp = []
        for c in seq:
            if c in complement_letters:
                comp.append(complement_letters[c])
            else:
                comp.append(c)
        return ''.join(comp)

    def reverse_complement(self, seq):
        seq = seq.upper()
        return self.complement(seq)[::-1]

    @property
    def chromosomes(self):
        """Return list of chromsome and their sizes
        as in the fasta file.

        Returns
        -------
        chroms : dict
                 Dictionary with {"chr": "Length"} format


        .. currentmodule:: .FastaReader
        .. autosummary::
            .FastaReader
        """
        chroms = {}
        for chrom in self.fasta.keys():
            chroms[chrom] = len(self.fasta[chrom])
        return chroms
