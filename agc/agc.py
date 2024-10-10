#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import gzip
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "BEN RADHIA Takwa"
__copyright__ = "Universite Paris Cite"
__credits__ = ["BEN RADHIA Takwa"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "BEN RADHIA Takwa"
__email__ = "takwa.ben-radhia@etu.u-paris.fr"
__status__ = "Developpement"



def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__,
                usage=f"{sys.argv[0]} -h")
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file',
            type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    with gzip.open(amplicon_file, "rt") as f:
        sequence = ""
        for line in f:
            if line.startswith(">"):
                if len(sequence) >= minseqlen:
                    yield sequence
                sequence = ""  # Reset sequence for the next one
            else:
                sequence += line.strip()  # Accumulate the sequence
        if len(sequence) >= minseqlen:
            yield sequence


def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence with a 
            count >= mincount and a length >= minseqlen.
    """
    seq_counts = Counter(read_fasta(amplicon_file, minseqlen))

    sorted_seq_counts = sorted(seq_counts.items(), key=lambda x: x[1], reverse=True)

    for sequence, count in sorted_seq_counts:
        if count >= mincount:
            yield [sequence, count]

def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences 
            in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    seq1, seq2 = alignment_list
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a!= '-' and b != '-')
    identity = matches / len(seq1) * 100
    return identity

def abundance_greedy_clustering(
    amplicon_file: Path,
    minseqlen: int,
    mincount: int,
    chunk_size: int = 0,
    kmer_size: int = 0) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    otu_list = []

    for sequence, count in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        is_otu = True

        for otu in otu_list:
            alignment = nw.global_align(sequence, otu[0], gap_open=-1, gap_extend=-1,
                        matrix=str(Path(__file__).parent / "MATCH"))
            identity = get_identity(alignment)

            if identity > 97:
                is_otu = False
                break

        if is_otu:
            otu_list.append([sequence, count])

    return otu_list


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, 'w', encoding="utf-8") as f:
        for i, (sequence, count) in enumerate(OTU_list, 1):
            f.write(f">OTU_{i} occurrence:{count}\n")
            f.write(textwrap.fill(sequence, width=80) + "\n")


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    otu_list = abundance_greedy_clustering(
        amplicon_file=args.amplicon_file,
        minseqlen=args.minseqlen,
        mincount=args.mincount,
        chunk_size=50,  # Chunk size (par défaut 50)
        kmer_size=8     # Kmer size (non utilisé cette année)
    )

    # Écrire les résultats des OTUs dans le fichier de sortie
    write_OTU(otu_list, args.output_file)


if __name__ == '__main__':
    main()
