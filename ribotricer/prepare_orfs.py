"""Functions for finding all candidate ORFs"""

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

from __future__ import annotations

import datetime
import re
from collections import defaultdict

from tqdm.autonotebook import tqdm

from .common import merge_intervals
from .fasta import FastaReader
from .gtf import GTFReader, GTFTrack
from .interval import Interval
from .orf import ORF

tqdm.pandas()

# Type alias for CDS ORF dictionary
CDSOrfs = defaultdict[str, defaultdict[str, ORF]]


def tracks_to_ivs(tracks: list[GTFTrack]) -> list[Interval]:
    """Convert GTF tracks to intervals.

    Parameters
    ----------
    tracks : list[GTFTrack]
        List of GTF tracks.

    Returns
    -------
    list[Interval]
        List of Interval objects.
    """
    chrom_set = {track.chrom for track in tracks}
    strand_set = {track.strand for track in tracks}
    if len(chrom_set) != 1 or len(strand_set) != 1:
        print("fail to fetch seq: inconsistent chrom or strand")
        intervals: list[Interval] = []
    else:
        chrom = list(chrom_set)[0]
        strand = list(strand_set)[0]
        intervals = [
            Interval(chrom, track.start, track.end, strand) for track in tracks
        ]
        intervals = merge_intervals(intervals)
    return intervals


def transcript_to_genome_iv(
    start: int,
    end: int,
    intervals: list[Interval],
    reverse: bool = False,
) -> list[Interval]:
    """Convert transcript coordinates to genome coordinates.

    Parameters
    ----------
    start : int
        Start position in transcript (0-based closed).
    end : int
        End position in transcript (0-based closed).
    intervals : list[Interval]
        Coordinate in genome (1-based closed).
    reverse : bool, optional
        Whether on the reverse strand, by default False.

    Returns
    -------
    list[Interval]
        The coordinate for start, end in genome.
    """
    total_len = sum(i.end - i.start + 1 for i in intervals)
    if reverse:
        start, end = total_len - end - 1, total_len - start - 1
    ivs: list[Interval] = []
    start_genome: int | None = None
    end_genome: int | None = None

    # find start in genome
    cur = 0
    for i in intervals:
        i_len = i.end - i.start + 1
        if cur + i_len > start:
            start_genome = i.start + start - cur
            break
        cur += i_len

    # find end in genome
    cur = 0
    for i in intervals:
        i_len = i.end - i.start + 1
        if cur + i_len > end:
            end_genome = i.start + end - cur
            break
        cur += i_len

    # find overlap with (start_genome, end_genome)
    for i in intervals:
        s = max(i.start, start_genome)
        e = min(i.end, end_genome)
        if s <= e:
            ivs.append(Interval(i.chrom, s, e, i.strand))
    return ivs


def fetch_seq(fasta: FastaReader | str, tracks: list[GTFTrack]) -> str:
    """Fetch sequence for GTF tracks.

    Parameters
    ----------
    fasta : FastaReader | str
        Instance of FastaReader or path to FASTA file.
    tracks : list[GTFTrack]
        List of GTF tracks.

    Returns
    -------
    str
        Combined sequence for the region.
    """
    intervals = tracks_to_ivs(tracks)
    if not isinstance(fasta, FastaReader):
        fasta = FastaReader(fasta)
    sequences = fasta.query(intervals)
    merged_seq = "".join(sequences)
    strand = tracks[0].strand
    if strand == "-":
        merged_seq = fasta.reverse_complement(merged_seq)
    return merged_seq


def search_orfs(
    fasta: FastaReader | str,
    intervals: list[Interval],
    min_orf_length: int,
    start_codons: set[str],
    stop_codons: set[str],
    longest: bool,
) -> list[tuple[list[Interval], str, str, str]]:
    """Search for ORFs in intervals.

    Parameters
    ----------
    fasta : FastaReader | str
        Instance of FastaReader or path to FASTA file.
    intervals : list[Interval]
        List of intervals.
    min_orf_length : int
        Minimum length (nts) of ORF to include.
    start_codons : set[str]
        Set of start codons.
    stop_codons : set[str]
        Set of stop codons.
    longest : bool
        Whether to choose the most upstream start codon when multiple
        in-frame ones exist.

    Returns
    -------
    list[tuple[list[Interval], str, str, str]]
        List of (intervals, seq, leader, trailer) tuples.
    """
    if not intervals:
        return []

    orfs = []
    if not isinstance(fasta, FastaReader):
        fasta = FastaReader(fasta)
    intervals = merge_intervals(intervals)
    sequences = fasta.query(intervals)
    merged_seq = "".join(sequences)
    reverse = False
    strand = intervals[0].strand
    if strand == "-":
        merged_seq = fasta.reverse_complement(merged_seq)
        reverse = True

    start_stop_idx = []
    if "ATG" in start_codons:
        start_stop_idx += [(m.start(0), "ATG") for m in re.finditer("ATG", merged_seq)]
    if start_codons - {"ATG"}:
        alternative_start_regx = re.compile(
            r"(?=({}))".format("|".join(start_codons - {"ATG"}))
        )
        start_stop_idx += [
            (m.start(0), "start")
            for m in re.finditer(alternative_start_regx, merged_seq)
        ]
    stop_regx = re.compile(r"(?=({}))".format("|".join(stop_codons)))
    start_stop_idx += [(m.start(0), "stop") for m in re.finditer(stop_regx, merged_seq)]
    start_stop_idx.sort(key=lambda x: x[0])
    for frame in [0, 1, 2]:
        inframe_codons = [x for x in start_stop_idx if x[0] % 3 == frame]
        starts = []
        for idx, label in inframe_codons:
            if label != "stop":
                starts.append(idx)
            elif starts:
                for start in starts:
                    if idx - start >= min_orf_length:
                        ivs = transcript_to_genome_iv(
                            start, idx - 1, intervals, reverse
                        )
                        seq = merged_seq[start:idx]
                        leader = merged_seq[:start]
                        trailer = merged_seq[idx + 3 :]
                        if ivs:
                            orfs.append((ivs, seq, leader, trailer))
                    if longest:
                        break
                starts = []
    return orfs


def check_orf_type(orf: ORF, cds_orfs: CDSOrfs) -> str:
    """Determine the ORF type relative to annotated CDS.

    Parameters
    ----------
    orf : ORF
        Instance of ORF.
    cds_orfs : CDSOrfs
        Dictionary of annotated CDS ORFs.

    Returns
    -------
    str
        Type of the candidate ORF.

    Notes
    -----
    This method uses a fail-fast strategy and hence multiple returns.
    """
    if orf.gid not in cds_orfs:
        return "novel"
    if orf.tid not in cds_orfs[orf.gid]:
        return "novel"
    matched_cds = cds_orfs[orf.gid][orf.tid]
    if orf.intervals == matched_cds.intervals:
        return "annotated"
    gene_cds = [cds_orfs[orf.gid][tid] for tid in cds_orfs[orf.gid]]
    gene_start = min([gc.intervals[0].start for gc in gene_cds])
    gene_end = max([gc.intervals[-1].end for gc in gene_cds])
    if orf.intervals[-1].end < gene_start:
        return "super_uORF" if orf.strand == "+" else "super_dORF"
    if orf.intervals[0].start > gene_end:
        return "super_dORF" if orf.strand == "+" else "super_uORF"
    if orf.intervals[0].start < matched_cds.intervals[0].start:
        if orf.intervals[-1].end < matched_cds.intervals[0].start:
            return "uORF" if orf.strand == "+" else "dORF"
        if orf.intervals[-1].end < matched_cds.intervals[-1].end:
            return "overlap_uORF" if orf.strand == "+" else "overlap_dORF"
    if orf.intervals[-1].end > matched_cds.intervals[-1].end:
        if orf.intervals[0].start > matched_cds.intervals[-1].end:
            return "dORF" if orf.strand == "+" else "uORF"
        if orf.intervals[0].start > matched_cds.intervals[0].start:
            return "overlap_dORF" if orf.strand == "+" else "overlap_uORF"
    return "internal"


def prepare_orfs(
    gtf: GTFReader | str,
    fasta: FastaReader | str,
    prefix: str,
    min_orf_length: int,
    start_codons: set[str],
    stop_codons: set[str],
    longest: bool,
) -> None:
    """Prepare candidate ORFs from GTF and FASTA files.

    Parameters
    ----------
    gtf : GTFReader | str
        Instance of GTFReader or path to GTF file.
    fasta : FastaReader | str
        Instance of FastaReader or path to FASTA file.
    prefix : str
        Prefix for output file.
    min_orf_length : int
        Minimum length (nts) of ORF to include.
    start_codons : set[str]
        Set of start codons.
    stop_codons : set[str]
        Set of stop codons.
    longest : bool
        Whether to choose the most upstream start codon when multiple
        in-frame ones exist.
    """

    now = datetime.datetime.now()
    print(now.strftime("%b %d %H:%M:%S ..... started ribotricer prepare-orfs"))
    now = datetime.datetime.now()
    print(now.strftime("%b %d %H:%M:%S ... starting to parse GTF file"))
    candidate_orfs = []
    if not isinstance(gtf, GTFReader):
        gtf = GTFReader(gtf)
    if not isinstance(fasta, FastaReader):
        fasta = FastaReader(fasta)

    # process CDS gtf
    now = datetime.datetime.now()
    print(now.strftime("%b %d %H:%M:%S ... starting extracting annotated ORFs"))
    cds_orfs = defaultdict(lambda: defaultdict(ORF))
    for gid in tqdm(gtf.cds, unit="lines", leave=False):
        for tid in gtf.cds[gid]:
            tracks = gtf.cds[gid][tid]
            seq = fetch_seq(fasta, tracks)
            orf = ORF.from_tracks(tracks, "annotated", seq=seq[:3])
            if orf:
                cds_orfs[gid][tid] = orf
                candidate_orfs.append(orf)

    now = datetime.datetime.now()
    print(
        "{} ... {}".format(
            now.strftime("%b %d %H:%M:%S"),
            "starting searching transcriptome-wide ORFs. This may take a long time...",
        )
    )
    for tid in tqdm(gtf.transcript, unit="transcripts", leave=False):
        tracks = gtf.transcript[tid]
        ttype = tracks[0].transcript_type
        gid = tracks[0].gene_id
        gname = tracks[0].gene_name
        gtype = tracks[0].gene_type
        chrom = tracks[0].chrom
        strand = tracks[0].strand
        ivs = tracks_to_ivs(tracks)
        orfs = search_orfs(
            fasta, ivs, min_orf_length, start_codons, stop_codons, longest
        )
        for ivs, seq, leader, trailer in orfs:
            orf = ORF(
                "unknown",
                tid,
                ttype,
                gid,
                gname,
                gtype,
                chrom,
                strand,
                ivs,
                seq=seq[:3],
            )
            orf.category = check_orf_type(orf, cds_orfs)
            if orf.category != "annotated" and orf.category != "internal":
                candidate_orfs.append(orf)

    # save to file
    now = datetime.datetime.now()
    print(now.strftime("%b %d %H:%M:%S ... saving candidate ORFs into disk"))
    columns = [
        "ORF_ID",
        "ORF_type",
        "transcript_id",
        "transcript_type",
        "gene_id",
        "gene_name",
        "gene_type",
        "chrom",
        "strand",
        "start_codon",
        "coordinate\n",
    ]

    to_write = "\t".join(columns)
    formatter = "{}\t" * (len(columns) - 1) + "{}\n"
    for orf in tqdm(candidate_orfs, unit="ORFs"):
        coordinate = ",".join([f"{iv.start}-{iv.end}" for iv in orf.intervals])
        if orf.start_codon in start_codons:
            to_write += formatter.format(
                orf.oid,
                orf.category,
                orf.tid,
                orf.ttype,
                orf.gid,
                orf.gname,
                orf.gtype,
                orf.chrom,
                orf.strand,
                orf.start_codon,
                coordinate,
            )

    with open(f"{prefix}_candidate_orfs.tsv", "w") as output:
        output.write(to_write)
    now = datetime.datetime.now()
    print(now.strftime("%b %d %H:%M:%S ... finished ribotricer prepare-orfs"))
