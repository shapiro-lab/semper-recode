# 🔬 Semper Recode

## Introduction

SemperRecode is a Python package designed for computational biologists and genetic engineers working with gene expression modulation. This package facilitates the detection of internal start codons in the translation initiation site (TIS) within a gene sequence and allows for tuning gene expression by replacing the TIS with a sequence that minimizes expression efficiency.

### What's SEMPER?
`also known as Stoichiometric Expression of Messenger Polycistrons by Eukaryotic Ribosomes`

SEMPER is an approach in mammalian synthetic biology that enables the expression of multiple proteins from a single transcript by leveraging the phenomenon knows as leaky ribosomal scanning. It involves placing short upstream open reading frames (uORFs) before genes of interest (GOIs) to divert ribosome flux and control the translation of downstream proteins. By varying the initiation strength of each uORF, SEMPER allows tunable and user-determined expression levels of multiple proteins. Compared to other methods like IRES or self-cleaving peptides (such as P2A), SEMPER offers compact genetic constructs without fusion proteins or reduced downstream expression. This versatile tool opens up new possibilities for engineering functional multimeric protein complexes and complex genetic circuits in mammalian cells.


## Features

- Detection of internal start codons within gene sequences.
- Detection of out-of-frame and in-frame translation initiation sites within gene sequences.
- Identification of the TIS position and associated efficiency score.
- Substitution of the TIS sequence with the lowest possible efficiency.


## Installation

1. Clone the repository to your local machine. You can do this by running the following command in your terminal or command prompt:

```shell
git clone https://github.com/ishaanjdev/semper-recode.git
```

2. SemperRecode can be installed via pip using the following command (to be implemented):

```shell
pip install SemperRecode
```

3. After the installation is complete, you can import the library in your Python script or interactive session by using the following statement (to be implemented):

```shell
import semper_recode
```

Alternatively

```shell
git clone https://github.com/ishaanjdev/semper-recode.git
pip install -r requirements.txt
pip install .
```

## Prerequisites

Before using this package, make sure you have the following dependencies installed:

- [Pandas](https://pandas.pydata.org/): A powerful data manipulation library.
- [Biopython](https://biopython.org/): A set of freely available tools for biological computation.
- [Pytest](https://docs.pytest.org/): A Python testing framework

You can install these dependencies using pip:

```shell
pip install -r requirements.txt
```

## To get started
**Parse .fasta file sequence by sequence**

```shell
with open({input_file_path}, 'r') as file:
    for seq in SeqIO.parse(file, 'fasta'): # Parsing fasta file sequence by sequence
        # Proceed with calling desired functions
```

Example:
```shell
with open("tests/sample_file/sample_file_inputs.fasta", 'r') as file:
    for seq in SeqIO.parse(file, 'fasta'):
        obj = SemperRecode(seq)
        modified_seq = obj.process_sequence()
```

**Convert the output into data frame (modified sequence and error list)**
Example:
```
error_list, data_list  = [], []
    
    # Read input sequences from the FASTA file and process them using the SemperRecode class
    with open("tests/sample_file/sample_file_inputs.fasta", 'r') as file:
        for line in SeqIO.parse(file, 'fasta'):
            input = str(line.seq)
            obj = SemperRecode(input)
            new_seq, error_list = obj.process_sequence()
            
            # Accumulate the processed data as dictionaries in a list
            data = {'ID': line.id, 'sequence': new_seq, 'error': error_list}
            data_list.append(data)

    # Create a DataFrame from the list of processed data dictionaries
    df = pd.DataFrame(data_list)

    # Export the DataFrame to a CSV file
    output_file = "tests/sample_file/sample_file_df_outputs.csv"
    df.to_csv(output_file, index=False)
```


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