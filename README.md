# 🔬 Semper Recode

## Introduction

SemperRecode is a Python package designed for computational biologists and genetic engineers working with gene expression modulation. This package facilitates the detection of internal start codons (TIS) within a gene sequence and allows for tuning gene expression by replacing the TIS with a sequence that minimizes expression efficiency.

****

## Features

- Detection of internal start codons within gene sequences.
- Identification of the TIS position and associated efficiency score.
- Substitution of the TIS sequence with the lowest possible efficiency.
- Optimized gene expression modulation for experimental and synthetic biology applications.

****

## Installation

1. Clone the repository to your local machine. You can do this by running the following command in your terminal or command prompt:

```shell
git clone https://github.com/ishaanjdev/semper-recode.git
```

2. SemperRecode can be installed via pip using the following command:

```shell
pip install SemperRecode
```

3. After the installation is complete, you can import the library in your Python script or interactive session by using the following statement:

```shell
import semper_recode
```

Alternatively

```shell
git clone https://github.com/ishaanjdev/semper-recode.git
cd SemperRecode
pip install .
```

****

## Prerequisites

Before using this package, make sure you have the following dependencies installed:

- [pandas](https://pandas.pydata.org/): A powerful data manipulation library.
- [Biopython](https://biopython.org/): A set of freely available tools for biological computation.
- [pytest](https://docs.pytest.org/): A Python testing framework

You can install these dependencies using pip:

```shell
pip install pandas biopython pytest
```

or

```shell
pip install -r requirements.txt
```

***

## To get started
**Parse .fasta file sequence by sequence**

```shell
with open({input_file_path}, 'r') as file:
    for seq in SeqIO.parse(file, 'fasta'): # Parsing sequence line by line
        # Proceed with calling desired functions
```

For example:
```shell
with open("tests/sample_file/sample_file_inputs.fasta", 'r') as file:
    for seq in SeqIO.parse(file, 'fasta'):
        obj = SemperRecode(seq)
        modified_seq = obj.process_sequence()
```
****

## Sample workflow

Example 
```shell
with open("tests/sample_file/sample_inputs.fasta", 'r') as file:
    for line in SeqIO.parse(file, 'fasta'):
        input = str(line.seq)
        obj = SemperRecode(input)
        new_seq = obj.process_sequence()
        modified_seq = SeqRecord(Seq(new_seq), id=f"{line.id}_semper_recode", description='')
        output.append(modified_seq)

output_file = "tests/sample_file/sample_outputs.fasta"

with open(output_file, 'w') as file:
    SeqIO.write(output, file, 'fasta')
```

Sample output
```shell
>EmGFP_mARG2.0_semper_recode
ATGGTGTCCAAGGGCGAGGAACTGTTCACCGGCGTGGTGCCCATCCTGGTGGAACTGGAT
GGCGACGTGAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAAGGCGACGCCACATAC
GGAAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCTTGGCCTACC
CTCGTGACCACACTGACCTACGGCGTGCAGTGCTTCGCCAGATACCCCGATCATATGAAA
CAGCACGATTTCTTCAAGAGTGCTATGCCTGAGGGCTACGTGCAGGAACGGACCATCTTC
TTCAAGGACGACGGCAACTACAAGACAAGAGCCGAAGTGAAGTTCGAGGGCGACACCCTC
GTGAACCGGATCGAGCTGAAGGGCATCGACTTCAAAGAGGATGGCAACATCCTGGGCCAC
AAGCTGGAGTACAACTACAACAGCCACAAGGTGTACATCACCGCCGACAAGCAGAAAAAC
GGCATCAAAGTGAACTTCAAGACCCGGCACAACATCGAGGACGGCAGCGTGCAGCTGGCC
GACCACTACCAGCAGAACACCCCCATCGGAGATGGCCCCGTGCTGCTGCCCGACAACCAC
TACCTGAGCACACAAAGCGCCCTGAGCAAGGACCCCAACGAGAAGCGGGATCATATGGTT
CTGCTGGAATTTGTGACCGCCGCTGGCATCACCCTTGGTATGGATGAGCTGTACAAGTGA

>emiRFP670_semper_recode
ATGGCGGAAGGCTCCGTCGCCAGGCAGCCTGACCTCTTGACCTGCGAACATGAAGAGATC
CACCTCGCCGGCTCGATCCAGCCGCATGGCGCGCTTCTGGTCGTCAGCGAACATGATCAT
CGCGTCATCCAGGCCAGCGCCAACGCCGCGGAATTTCTGAATCTCGGAAGCGTACTCGGC
GTTCCGCTCGCCGAGATCGACGGCGATCTGTTGATCAAGATCCTGCCGCATCTCGATCCC
ACCGCCGAAGGGATGCCGGTCGCGGTGCGCTGCCGGATCGGCAATCCCTCTACGGAGTAC
TGCGGGTTGATGCATCGGCCTCCGGAAGGCGGGCTGATCATCGAACTCGAACGTGCCGGC
CCGTCGATCGATCTGTCAGGCACGCTGGCGCCGGCGCTGGAGCGGATCCGCACGGCGGGT
TCACTGCGCGCGCTGTGCGATGACACCGTGCTGCTGTTTCAGCAGTGCACCGGCTACGAC
CGTGTTATGGTGTATCGTTTCGATGAGCAAGGCCACGGCCTGGTATTCTCCGAGTGCCAT
GTGCCTGGGCTCGAATCCTATTTCGGCAACCGCTATCCGTCGTCGACTGTCCCACAGATG
GCGCGGCAGCTGTACGTGCGGCAGCGCGTCCGCGTGCTGGTCGACGTCACCTATCAGCCG
GTGCCGCTGGAGCCGCGGCTGTCGCCGCTGACCGGGCGCGATCTTGATATGAGTGGCTGC
TTCCTGCGGTCTATGAGTCCGTGCCATCTGCAGTTCCTGAAGGATATGGGCGTGCGCGCC
ACCCTGGCGGTGTCGCTGGTGGTCGGCGGCAAGCTGTGGGGCCTGGTTGTCTGTCACCAT
TATCTGCCGCGCTTCATCCGTTTCGAGCTGCGGGCGATCTGCAAACGGCTCGCCGAAAGG
ATCGCGACGCGGATCACCGCGCTTGAGAGCTAA

>GvpA_mARG2.0_semper_recode
ATGGCCGTGGAAAAGACCAACAGCAGCAGCTCCCTGGCCGAAGTGATCGACAGAATCCTG
GACAAGGGCATCGTGATCGACGCCTGGGTGCGCGTGTCCCTCGTGGGAATTGAGCTGCTG
GCCATCGAGGCCCGGATCGTGATTGCCAGCGTGGAAACCTACCTGAAGTACGCCGAGGCC
GTGGGCCTGACACAGAGTGCTGCTGTGCCTGCTTGA
```