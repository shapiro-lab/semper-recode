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

****

## Prerequisites

Before using this package, make sure you have the following dependencies installed:

- [pandas](https://pandas.pydata.org/): A powerful data manipulation library.
- [Biopython](https://biopython.org/): A set of freely available tools for biological computation.
- [pytest](https://docs.pytest.org/): A Python testing framework

You can install these dependencies using pip:

```shell
pip install pandas numpy biopython pytest
```

After installing the package and its dependencies, you can import the required modules in your Python script as follows:

```python
import os
import pickle
import warnings
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
```

***

## To get started
**To parse .fasta file sequence by sequence**

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

        # Assert the generated sequence match the correctly modified sequence
        assert new_seq in groud_truth 

        modified_seq = SeqRecord(Seq(new_seq), id=f"{line.id}_semper_recode", description='')
        output.append(modified_seq)

output_file = "tests/sample_file/sample_outputs.fasta"

with open(output_file, 'w') as file:
    SeqIO.write(output, file, 'fasta')
```