import os
import warnings

import numpy as np
import pyBigWig


class WigReader(object):
    """Class for reading and querying wigfiles."""

    def __init__(self, wig_location):
        """
        Parameters
        ---------
        wig_location : string
                       Path to wig file

        """
        self.wig_location = wig_location
        try:
            self.wig = pyBigWig.open(self.wig_location)
        except Exception as e:
            raise Exception('Error reading wig file {} : {}'.format(
                os.path.abspath(self.wig_location), e))

    def query(self, intervals):
        """ Query regions for scores.

        Parameters
        ----------
        intervals : list of Interval 

        Returns
        -------
        scores : array of array
                 A numpy array containing scores for each Interval
                 This function is agnostic of the strand information,
                 the position in the scores is corresponding to the interval

        .. currentmodule:: .WigReader
        .. autosummary::
            .WigReader

        """
        scores = []
        chrom_lengths = self.chromosomes
        for i in intervals:
            if i.chrom not in list(chrom_lengths.keys()):
                warnings.warn(
                    'Chromosome {} does not appear in the bigwig'.format(
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
            score = self.wig.values(i.chrom, i.start, i.end)
            scores.append(score)
        return np.array(scores)

    @property
    def chromosomes(self):
        """Return list of chromsome and their sizes
        as in the wig file.

        Returns
        -------
        chroms : dict
                 Dictionary with {"chr": "Length"} format


        .. currentmodule:: .WigReader
        .. autosummary::
            .WigReader
        """
        return self.wig.chroms()
