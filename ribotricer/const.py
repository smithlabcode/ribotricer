"""Constants used in ribotricer"""
# Part of ribotricer software
#
# Copyright (C) 2020 Saket Choudhary, Wenzheng Li, and Andrew D Smith
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# ribotricer default cutoff for lavbeling ORFs 'translating'
CUTOFF = 0.428571428571
# p-site offset
TYPICAL_OFFSET = 12
# minimum number of valid codons required in an ORF to label
# it 'translating'
MINIMUM_VALID_CODONS = 5
# minimum number of reads required per codon in an ORF to label
# it 'translating'
# default: 0 (decided by CUTOFF and MINIMUM_VALID_CODONS)
MINIMUM_READS_PER_CODON = 0
# fraction of codons with non zero reads
MINIMUM_VALID_CODONS_RATIO = 0
# Minimum read density over ORF
# defined as the number of reads per unit length of the ORF
MINIMUM_DENSITY_OVER_ORF = 0.0
