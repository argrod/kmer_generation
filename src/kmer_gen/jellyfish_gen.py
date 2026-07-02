"""
Module of functions to utilise Jellyfish to generate standardised kmer analysis.

Generated kmer tables are written to the passed output directory.

Usage:
```
jellyfish_gen.py --
"""

import gzip
import os
import subprocess
import tempfile
import argparse
import re
import numpy as np
from glob import glob
from Bio import SeqIO
from itertools import product
from pathlib import Path

JELLYFISH_BINARY = "jellyfish"
# MODULES_LOADING = "module load use.dev intel/2022.2 mkl gcc13/13.1.0"


def create_complete_kmer_template(
    k: int,
) -> dict[str, int]:
    """Create template with all possible kmers

    Params
    ------
    k : int
        Length of k-mer.

    Returns
    -------
    dict[str, int]
        List of all possible kmers.
    """
    kmers = ["".join(p) for p in product("ACGT", repeat=k)]
    template = {}
    for kmer in sorted(kmers):
        template[kmer] = 0
    return template


# Instead of default temp directory
temp_dir = os.path.expanduser("~/tmp")
os.makedirs(temp_dir, exist_ok=True)


def process_sequence(
    seq_id: str,
    sequence: str,
    kmer_size: int,
    outdir: Path,
    template: dict[str, int] | None = None,
    identifier: str | None = None,
) -> None:
    """
    Process single sequence and return complete kmer counts.

    Parameters
    ----------
    seq_id : str
        Specific sequence ID.
    sequence : str
        RNA sequence.
    kmer_size : int
        Length if kmer.
    template : list[str] | None
        Zero kmer template (all kmers).
    outdir : str
        Desired directory for output.
    identifier : str | None
        Optional filename. If None, passes seq_id.
    """
    # define temporary directory for intermediate (fasta) files
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)

        # define the absolute Path to temp file
        fasta_path = temp_path / f"{seq_id}.fasta"
        jf_file_path = temp_path / f"{seq_id}.jf"

        # create temporary fasta file
        with open(fasta_path, "w") as temp_fasta:
            temp_fasta.write(f">{seq_id}\n{sequence}\n")

        # run jellyfish count with absolute paths
        count_cmd = (
            f"{JELLYFISH_BINARY} count -m {kmer_size} -s 100M -C "
            f"-o {jf_file_path.absolute()} {fasta_path.absolute()}"
        )
        subprocess.run(
            count_cmd,
            shell=True,
            check=True,
        )

        # Run jellyfish dump to capture the k-mer table to result
        dump_cmd = f"{JELLYFISH_BINARY} dump -c -t {jf_file_path.absolute()}"
        result = subprocess.run(
            dump_cmd,
            shell=True,
            capture_output=True,
            text=True,
            check=True,
        )

    # Parse observed counts
    if template is None:
        observed_counts = {}
    else:
        observed_counts = template.copy()

    for line in result.stdout.strip().split("\n"):
        if line:
            kmer, count = line.split("\t")
            observed_counts[kmer] = int(count)

    if sum(list(observed_counts.values())) == 0:
        raise ValueError(f"No kmers counted for {seq_id}")

    if identifier is None:
        identifier = seq_id

    # Write complete table
    output_file = f"{identifier}_{kmer_size}mers.txt"
    with open(outdir / output_file, "w") as f:
        for kmer in sorted(observed_counts.keys()):
            f.write(f"{kmer}\t{observed_counts[kmer]}\n")

    print(
        f"{kmer_size}mer counts calculated for {seq_id}\n"
        f"\tData saved at: {output_file}"
    )


def parse_fasta_gz(
    file_path: Path | str,
):
    """
    Parse gzipped fasta file.

    Parameters
    ----------
    file_path: Path | str
        Path to gzipped fasta file.
    """
    if str(file_path).split(".")[-1] == "gz":
        with gzip.open(file_path, "rt") as f:
            for read in SeqIO.parse(f, "fasta"):
                seq_id = read.id
                sequence = read.seq
                yield seq_id, sequence
    elif str(file_path).split(".")[-1] == "fa":
        for read in SeqIO.parse(file_path, "fasta"):
            seq_id = read.id
            sequence = read.seq
            yield seq_id, sequence


def gen_kmer_files(
    fasta_filepath: Path | str,
    kmer_length: int,
    outdir: Path,
    seq_ids: list[str] | None = None,
    accession: str | None = None,
    seq_subsets: dict[str, int | str] | None = None,
) -> None:
    """Generate standardised kmer files for a fasta file and given sequence IDs.

    Params
    ------
    species : str
        Species name.
    fasta_filepath : Path | str
        Path to fasta file.
    kmer_length : int
        Desired length of kmers.
    outdir : Path
        Output directory for kmer files.
    seq_ids : list[str], optional
        List of sequence IDs. If None, all kmers are calculated for all
        sequences.
    accession : str | None
        Optional accession name. If None given, output named the same as
        sequence ID. Defaults to None.
    seq_subsets : list[int | str] | None
        Optional argument to generate kmer tables for subsets of sequences.
        Should be of same length as sequence IDs. Format should be start and end
        indices, or string of index with colon prefix or suffix to indicate all
        prior or all following, respectively.
    """
    # generate a template
    template = create_complete_kmer_template(kmer_length)

    for seq_id, sequence in parse_fasta_gz(fasta_filepath):
        # create an identifier if details given
        if accession is None:
            identifier = seq_id
        if seq_ids is None:
            run_sequence = True
        elif seq_id in seq_ids:
            run_sequence = True
        else:
            run_sequence = False
        if run_sequence:
            if seq_subsets is not None:
                seq_subset = seq_subsets[seq_id]
                for seq_sub_indices, seq_category in zip(seq_subset, ["upstream", "CDS", "downstream"]):
                    if isinstance(seq_sub_indices[0], str):
                        if seq_sub_indices[0][0] == ":":
                            seq_sub_indices = [0, int(seq_sub_indices[0][1:])]
                        elif seq_sub_indices[0][-1] == ":":
                            seq_sub_indices = [
                                int(seq_sub_indices[0][:-1]),
                                len(sequence),
                            ]
                        if seq_sub_indices[1] < seq_sub_indices[0]:
                            raise ValueError(
                                f"Sequence index error: \
                                    end index smaller than start index for {seq_id}"
                            )
                    # pass if sequence length shorter than kmer length
                    if seq_sub_indices[1] - seq_sub_indices[0] < kmer_length:
                        print(f"No kmers generated for {seq_id}_{seq_category}")
                        continue
                    else:
                        subsequence = sequence[seq_sub_indices[0] : seq_sub_indices[1]]
                        process_sequence(
                            seq_id=seq_id,
                            sequence=subsequence,
                            kmer_size=kmer_length,
                            outdir=outdir,
                            template=template,
                            identifier=f"{identifier}_{seq_category}",
                        )
            else:
                process_sequence(
                    seq_id=seq_id,
                    sequence=sequence,
                    kmer_size=kmer_length,
                    outdir=outdir,
                    template=template,
                    identifier=f"{identifier}_complete",
                )


def build_kmer_count_table(
    counts_directory: Path,
) -> tuple[list[str], np.ndarray]:
    """
    Convert output of `gen_kmer_files` (txt files) to a single numpy array of
    counts and list of (file-based) identifiers.

    Parameters
    ----------
    counts_directory : Path
        Path to directory containing txt files.

    Returns
    -------
    list[str]
        List of identifiers corresponding to count array.
    np.ndarray
        Array of kmer counts.
    """
    count_files = glob(str(counts_directory) + "/*.txt")
    identifiers, gene_segment, nmers = map(np.array, zip(*(
        re.search(r".*/([A-Za-z0-9.]+)_([A-Za-z]+)_([0-9]+)mers.txt", x).groups()
        for x in count_files
    )))
    nmers = nmers.astype(int)
    if len(np.unique(nmers)) != 1:
        raise ValueError(
            f"Multiple kmer lengths found for txt files in {counts_directory}\n"
            f"Expected 1 kind of kmer, got '{np.unique(nmers).tolist()}"
        )
    if len(np.unique(gene_segment)) != 1:
        raise ValueError(
            f"Kmer files for multiple gene segments found in {counts_directory}\n"
            f"Expected 1 type of segment, got '{np.unique(gene_segment).tolist()}"
        )
    kmer_counts = np.zeros((len(count_files), 4**np.unique(nmers).item()))
    for idx, seq_f in enumerate(count_files):
        kmer_count_row = np.loadtxt(
            seq_f,
            usecols=[1],
            delimiter="\t",
        )
        kmer_counts[idx] = kmer_count_row
    return identifiers.tolist(), kmer_counts