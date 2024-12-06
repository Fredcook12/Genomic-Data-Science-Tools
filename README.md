# Genomic-Data-Science-Tools
This repository provides a collection of tools and scripts designed to assist in analysing genomic data.

# DNA Tools for Sequence Analysis

This script provides a set of tools for DNA sequence analysis. It is designed to work with DNA sequences in FASTA format and offers functionality to analyse sequence lengths, identify open reading frames (ORFs), and detect repeated sequences.

## Features

1. **Count Records**: Count the number of DNA sequences in a multi-FASTA file.
2. **Check Sequence Lengths**: Determine the lengths of sequences, including the longest and shortest.
3. **ORF Identifier**: Identify open reading frames (ORFs) in DNA sequences.
4. **Reverse Complement Generator**: Compute the reverse complement of a DNA sequence.
5. **Repeat Finder**: Detect and count repeats of a specified length in a DNA sequence.

## Example usage

1. **Count Records**
  print("Count Records")
  dna_tools = dna_tool_sets("../data/example.fasta")  # Replace with your file path
  dna_tools.count_records()  # Outputs the number of records in the file

2. **Check Sequence Lengths**
  print("Check Sequence Lengths")
  dna_tools.check_length()  # Outputs the length of the longest and shortest sequences

3. **Identify ORFs**
  print("Identify ORFs")
  dna_tools.orf_identifier()  # Outputs ORFs' start positions and lengths for each sequence

4.**Find Repeats**
  print("\Find Repeats")
  dna_sequence = "ATGATGATG"  # Replace with your DNA sequence
  repeats = dna_tools.find_repeats(dna_sequence, 3)  # Finds repeats of length 3
  print(repeats)  # Outputs a dictionary of repeated sequences and their counts
