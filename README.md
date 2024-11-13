# OTU Clustering

## Introduction

This Python script performs OTU (Operational Taxonomic Units) clustering from DNA sequence data, using deduplication and alignment methods. It identifies representative OTU sequences from an amplicon dataset based on their abundance and similarity.

## Author
Takwa Ben Radhia

## Description
The script performs the following tasks:

1. Reading and Filtering Sequences: Reads a compressed FASTA file containing amplicon DNA sequences and filters those that meet length and abundance criteria.
2. Dereplication and Counting: Groups identical sequences and counts their frequency.
3. Greedy Abundance Clustering: Uses a greedy clustering method based on sequence abundance and identity to identify OTUs.
4. Writing OTU Sequences: Saves the identified OTU sequences to an output file in FASTA format.

## Dependencies
- Python 3
- nwalign3: A library for global sequence alignment
- Alignment matrix file for nwalign3 (e.g., BLOSUM or MATCH)

## Usage
Running the Program
To run the script, use the following command while specifying the required arguments:

```bash
python3 script.py -i <amplicon_file> -s <minseqlen> -m <mincount> -o <output_file>
```

## Arguments
- -i, --amplicon_file: The compressed amplicon file in .fasta.gz format (required).
- -s, --minseqlen: The minimum sequence length for dereplication (default: 400).
- -m, --mincount: The minimum occurrence count for dereplication (default: 10).
- -o, --output_file: The output file where the OTU sequences will be saved (default: OTU.fasta).

## Example Command
```bash
python3 script.py -i example_amplicon.fasta.gz -s 400 -m 10 -o OTU_output.fasta
```
## Main Functions
1. isfile(path): Checks if the provided path is an existing file.
2. get_arguments(): Parses command-line arguments and returns an argument object.
3. read_fasta(amplicon_file, minseqlen): Reads a compressed FASTA file and extracts sequences that meet the minimum length requirement.
4. dereplication_fulllength(amplicon_file, minseqlen, mincount): Dereplicates sequences and counts their occurrences.
5. get_identity(alignment_list): Computes the identity rate between two aligned sequences.
6. abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size): Performs greedy clustering based on sequence abundance and identity to identify OTUs.
7. write_OTU(OTU_list, output_file): Writes the identified OTU sequences to an output FASTA file.

## Notes
- The amplicon file must be in .fasta.gz compressed FASTA format.
- The parameters chunk_size and kmer_size are present but not used in this version of the script.
- The script uses the MATCH matrix for alignment, which should be available in the same directory as the script.
## Support
For any questions, contact: takwa.ben-radhia@etu.u-paris.fr

