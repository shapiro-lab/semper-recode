=============
Semper Recode
=============

Introduction
------------
SemperRecode is a Python package designed for computational biologists and genetic engineers working with gene expression modulation. 
This package facilitates the detection of internal translation initiation sites (TIS), composed of a start codon and adjacent nucleotides, within a gene sequence and allows for tuning gene expression by replacing the TIS with a sequence that minimizes the probability of translation initiation.

.. image:: https://www.biorxiv.org/content/biorxiv/early/2023/05/26/2023.05.26.541240/F1.large.jpg
   :alt: A) Architecture of an mRNA transcript produced by transfected 2-ORF SEMPER plasmid DNA. Cap-dependent ribosomes translate ORF 1, msfGFP[r5M], or ORF 2, mEBFP2, with frequencies dependent on the trinucleotide (NNN) upstream of each ORF. The relative strengths of the trinucleotides and TISs they represent are depicted. *** (TTTCCAT) does not contain a start codon, preventing translation of the ORF. IRES-mCherry is included to normalize fluorescence measurements. B) mEBFP2 distribution plots for constructs containing different versions of monomeric Superfolder GFP encoded in ORF 1. Removal of internal methionines from the msfGFP leads to much stronger expression of the downstream mEBFP2. Both ORFs used the strong ACC TIS (N=1, representative of four replicates). 
   :align: center
   :width: 370

------------

Features
--------
- Detection of internal start codons within gene sequences.
- Identification of the TIS position and associated efficiency score.
- Substitution of the TIS sequence with the lowest possible efficiency.
- Optimized gene expression modulation for experimental and synthetic biology applications.

------------

Installation
------------
1. Clone the repository to your local machine. You can do this by running the following command in your terminal or command prompt:

.. code-block:: shell

   git clone https://github.com/ishaanjdev/semper-recode.git

2. SemperRecode can be installed via pip using the following command:

.. code-block:: shell

   pip install SemperRecode

3. After the installation is complete, you can import the library in your Python script or interactive session by using the following statement:

.. code-block:: shell

   import semper_recode

Alternatively

.. code-block:: shell

   git clone https://github.com/ishaanjdev/semper-recode.git
   cd SemperRecode
   pip install .

------------

Prerequisites
-------------
Before using this package, make sure you have the following dependencies installed:

- `pandas <https://pandas.pydata.org/>`_: A powerful data manipulation library.
- `Biopython <https://biopython.org/>`_: A set of freely available tools for biological computation.
- `pytest <https://docs.pytest.org/>`_: A Python testing framework

You can install these dependencies using pip:

.. code-block:: shell

   pip install pandas biopython pytest

or

.. code-block:: shell

   pip install -r requirements.txt

------------

To get started
--------------
**Parse .fasta file sequence by sequence**

.. code-block:: shell

   with open({input_file_path}, 'r') as file:
       for seq in SeqIO.parse(file, 'fasta'): # Parsing sequence line by line
           # Proceed with calling desired functions

For example:

.. code-block:: shell

   with open("tests/sample_file/sample_file_inputs.fasta", 'r') as file:
       for seq in SeqIO.parse(file, 'fasta'):
           obj = SemperRecode(seq)
           modified_seq = obj.process_sequence()

------------

Sample workflow
---------------

.. code-block:: shell

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

Sample output

.. code-block:: shell

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

------------

Credits
-------
Ishaan Dev and Gayvalin Sujaritchai

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
