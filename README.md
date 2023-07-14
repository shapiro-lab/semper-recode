# 🔬 Semper Recode

## Introduction

SemperRecode is a Python package designed for computational biologists and genetic engineers working with gene expression modulation. This package facilitates the detection of internal start codons (TIS) within a gene sequence and allows for tuning gene expression by replacing the TIS with a sequence that minimizes expression efficiency.

## Features

- Detection of internal start codons within gene sequences.
- Identification of the TIS position and associated efficiency score.
- Substitution of the TIS sequence with the lowest possible efficiency.
- Optimized gene expression modulation for experimental and synthetic biology applications.

## Installation

SemperRecode can be installed via pip using the following command:

```shell
pip install SemperRecode
```

Alternatively

```shell
git clone https://github.com/ishaanjdev/semper-recode.git
cd SemperRecode
pip install .
```

## To get started
How to parse .fasta file sequence by sequence

```shell
with open([input_file_path], 'r') as file:
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
