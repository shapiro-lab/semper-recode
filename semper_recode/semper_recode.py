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
import os

# import ishaan_utils

class SemperRecode:
    def __init__(self, user_seq = None, current_path = "data/", input_file_path = None):
        """
        Initializes the SemperRecode object.

        Parameters
        ----------
        user_seq : str or Seq object, optional
            User-defined sequence, by default None. In case user want to tune single 
            sequence instead of importing a whole FASTA file
        current_path : str
            File path to user's data which will be use to append to file names when file path is required
        input_file_path : str
            File path to user's input fasta file

        """
        self.start_codon = ['AUG', 'ATG', 'GUG', 'UUG'] # Use in find_in_frame()
        self.four_letter_codes = [] # Use in filtered_sequence()
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
        master_df_path = os.path.join(path, 'master_df_os_2023.csv')
        if not os.path.exists(master_df_path):
            raise FileNotFoundError(f"Master dataframe file not found: {master_df_path}")
        
        self.master_df = pd.read_csv(path + 'master_df_os_2023.csv')
        self.master_df = self.master_df.drop(columns="sequence") 

        # Four letters codes
        self.four_letter_codes = list(self.master_df['4-letters'].unique())  # Get unique four-letter codes from master_df

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
            for record in SeqIO.parse(file, 'fasta'): # Parsing line by line
                sequence = str(record.seq)
                replace_sequence = self.modify_TIS(sequence)
                modified_sequences.append(replace_sequence)

        return modified_sequences

    def modify_TIS(self, sequence):
        '''
        Takes in sequence (str or Seq object) and return modified sequence
        with lower TIS efficiency (str ot Seq object)

        Parameters
        ----------
            sequence : str or Seq object

        Returns
        -------
            modified_sequence (str)
        '''
        new_seq = list(sequence)
        index = self.find_in_frame(sequence) # Get the indices of AUG(s)
        df = self.master_df

        '''
        Iterate through the length of the sequence (len(sequence)) and if the index in found in
        index (list of positions from find_in_frame function), find a sequence with lower efficiency,
        else, concatnate the old sequence to the output list (new_seq)
        '''

        for pos in index:
            #Ignore the first AUG
            # Iterate through
            if pos != 0:
                internal_TIS_seq = Seq(sequence[pos-6:pos+5])
                aa4 = Seq(sequence[pos-6:pos+6]).translate()

                # Get the current efficiency level 
                current_eff = self.efficiency_level(internal_TIS_seq)

                filtered = df[df['4-letters'] == str(aa4)]
                new_eff = filtered['efficiency'].iloc[0]

                '''
                If a new sequence with a lower efficiency is found,
                replace the sequence in the string
                '''
                if(int(new_eff) < current_eff):
                    new_seq[pos-6:pos+6] = filtered["4-codons"][0]

        return ''.join(new_seq)

    def efficiency_level(self, sequence):
        '''
        Takes in sequence and return efficiency level

        Parameters
        ----------
        sequence : str or Seq obj
            The input sequence to find efficiency level.
        
        Returns
        -------
        eff : int
            An interger represent the efficiency level of the sequence when
            being compared with the master dataframe

        '''
        df = self.master_df

        '''
        Find the matching row in master_df with the same value in tis-sequence
        as the input sequence 
        '''
        match = df[df["tis-sequence"] == str(sequence)]
        if match.empty:
            # No matching is found
            raise ValueError(f"The sequence {sequence} is not found in the dataframe")

        return match['efficiency'].iloc[0]

    def find_in_frame(self, sequence):
        """
        Takes in a sequence (str) and returns a list of indices where the AUG codon is found.

        Parameters
        ----------
        sequence : str
            The input sequence to search for AUG codons.

        Returns
        -------
        pos : list
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
        
        '''
        Output: 
        >Sequence1
        ATCGATCGATCG

        >Sequence2
        GCTAGCTAGCTA

        >Sequence3
        TGCAATGCAATG
        '''
        

    def filtered_sequence_eff(self, efficiency):
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


