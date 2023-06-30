'''Install biopython and pytest using pip if not already installed

  1. pip install biopython
  2. pip install pytest

'''

import pandas as pd
import numpy as np
import re
from collections import Counter
import pickle
import itertools
from itertools import product
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# import ishaan_utils

class SemperRecode:
    def __init__(self, user_seq = None, current_path = "../data/", input_file_path = None):
        """
        Initializes the SemperRecode object.

        Parameters
        ----------
        user_seq : str, optional
            User-defined sequence, by default None. In case user want to tune single 
            sequence instead of importing a whole FASTA file
        current_path : str
            File path to user's data which will be use to append to file names when file path is required
        input_file_path : str
            File path to user's input fasta file
        """
        self.start_codon = ['AUG', 'ATG', 'GUG', 'UUG'] # Use in find_in_frame()
        self.four_letter_codes = None # Use in filtered_sequence()
        self.master_df = None
        self.load_data(current_path)
        self.input_file = input_file_path # Path to the input fasta file
        self.seq = Seq(user_seq) if user_seq is not None else 'ATGCTGACGGTAUGGACTTACCTGTATGCGTGCTAAATGCTAAGGCTGGTGCCGACCGGACCGTTGGGAGCGCTGTTGACCGGATGCTAAAGGGCCCGAGTCTTGTAGTACCGGACTTAAATGCGTTGTTTGACACCTGTT'
    
    def load_data(self, path):
        """
        Loads data required for the analysis. Reads various data files and initializes necessary attributes.
        
        Parameters
        ----------
        path : str
            User-defined sequence, by default None. In case user want to tune single sequence 
            instead of importing a whole FASTA file

        Returns
        -------
        None

        Raises
        ------
        FileNotFoundError
            If any of the required data files are not found.

        """
        # Master dataframe
        self.master_df = pd.read_csv(path + 'master_df_os_2023.csv')
        
        # Four letters codes
        with open(path + 'four_letter_codes.pkl', 'rb') as file:
            self.four_letter_codes = pickle.load(file)

    def process_sequence(self, file_path):
        """
        Takes in a fasta file path and returns the modified sequence with lower efficiency (if possible).

        Parameters
        ----------
        file_path : str
            Path to the input fasta file.

        Returns
        -------
        list
            List of modified sequences of the input file.

        Raises
        ------
        FileNotFoundError
            If the file_path does not exist.

        """
        modified_sequences = []

        with open(file_path, 'r') as file:
            for record in SeqIO.parse(file, 'fasta'):
                sequence = str(record.seq)
                modified_sequences.append(self.modify_TIS(sequence))

        export_file_name = input("Enter the name of export file")
        self.to_fasta(self, modified_sequences, export_file_name)   

        return modified_sequences

    def modify_TIS(self, sequence):
        '''
        Takes in sequence (str) and return modified sequence with lower TIS efficiency

        Arguments:
            sequence (str)
        Return:
            modified_sequence (str)
        '''

    def find_in_frame(self, sequence):
        """
        Takes in a sequence (str) and returns a list of indices where the AUG codon is found.

        Parameters
        ----------
        sequence : str
            The input sequence to search for AUG codons.

        Returns
        -------
        list
            A list of integers representing the indices of in-frame AUG codons in the sequence.

        """
        pos = []  # position(s) of in-frame AUG in the string (sequence)

        # Loop through codon by codon (start from 0, end at last codon, increase by 3)
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i+3]

            if codon in self.start_codon:
                pos.append(i)

        return pos

    def to_fasta(self, sequence, output_file_name):
        '''
        Takes in a list of modified sequences and converts them back to FASTA format for exporting to users.

        Parameters
        ----------
        sequence : list
            The modified list of sequences to be converted to a FASTA file.
        output_file_name : str
            The desired name of the output FASTA file.

        Returns
        -------
        None

        '''
        output = []

        for i, seq in enumerate(sequence):
            temp = SeqRecord(Seq(seq), id=f'Seq{i+1}', description='')
            output.append(temp)
        '''
        i = index as wel iterated through sequence 

        In each iteration, a new SeqRecord object (Seq object but with additional metadata)
        is being created 

        Metadata:
            id: index+1
            description: None
1
        '''
        with open(output_file_name, 'w') as file:
            SeqIO.write(output, file, 'fasta')

    def filtered_sequence(self, efficiency):
        """
        Takes in efficiency level input and finds the sequences in the dataframe by filtering,
        then displays all the possible sequences that have an efficiency lower than the given value
        and the percentage of found and missing amino acid sequences in the form of a dataframe.

        Parameters
        ----------
        efficiency : float
            The efficiency level to filter the dataframe.

        Returns
        -------
        pd.DataFrame
            A dataframe containing the filtered sequences.
        """
        df = self.master_df

        filtered = df[df['efficiency'] < efficiency]
        count = filtered.shape[0]
        percent = (count / df.shape[0]) * 100

        groups = filtered.groupby('4-codons')
        seq_names = groups['4-letters'].apply(lambda x: ', '.join(x.astype(str).unique())).tolist()

        set1 = set(self.four_letter_codes)
        set2 = set(seq_names)

        missing_elements = list(set1 - set2)

        found_percent = len(set2) / len(set1) * 100

        print(f"Total: {count} sequence(s) ({percent:.2f}%)")
        print(f"Found ({found_percent:.2f}%): {', '.join(set2)}")
        print(f"Missing ({(100-found_percent):.2f}%): {', '.join(missing_elements)}\n")

        for each in groups:
            print(each, "\n")

        return filtered.DataFrame
    
    def find_lower_eff_sequence(self, efficiency, seq):
        """
        Takes in a nucleotide sequence and checks if there is another sequence that produces
        the same protein but with a lower efficiency level than the input level.

        Parameters
        ----------
        efficiency : float
            The efficiency level to filter the dataframe.
        seq : Bio.Seq
            nucleotide sequence

        Raises
        ------
        ValueError
            If the given sequence doesn't have replacing sequence with efficiency lower
            than the input efficiency level.


        Returns
        -------
        Bio.Seq
            Nucleotide sequence with lower efficiency level
        """

        possible_seq_df = self.filtered_sequence(efficiency, self.master_df)


