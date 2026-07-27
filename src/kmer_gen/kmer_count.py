"""
Simple kmer calculation function.

Authors: Aran Garrod <aran.garrod@gmail.com>
Date last updated: July 2026
"""
from collections import Counter

def calc_kmer(
    seq: str,
    k: int,
) -> Counter:
    """
    Kmer count calculation.

    Parameters
    ----------
    seq: str
        Sequence.
    k: int
        Desired kmer length.

    Returns
    -------
    Counter
        Dict subclass containing all present kmers in sequence.
    """
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    counts = Counter(kmers)
    return counts