import pandas as pd
import numpy as np
import re
from collections import Counter
import pickle
from itertools import product
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os

'''
Global variables declaration

Loads data required for the analysis. Reads various data files and initializes necessary attributes.

Raises
------
FileNotFoundError
    If any of the required data files are not found.

'''
start_codon = ['ATG']
path = "data/" 

# Master dataframe
master_df_path = os.path.join(path, 'master_df_os_2023.csv')
if not os.path.exists(master_df_path):
    raise FileNotFoundError(f"Master dataframe file not found: {master_df_path}")
        
master_df = pd.read_csv(path + 'master_df_os_2023.csv')

# Load codon list
codon_list_path = os.path.join(path, 'codon_list.pkl')
if not os.path.exists(codon_list_path):
    raise FileNotFoundError(f"Codon list (.pkl) file not found: {codon_list_path}")

with open(codon_list_path, 'rb') as file:
    codon_list = pickle.load(file)

if codon_list == {}:
    raise ValueError("codon_list is empty")

# Load codon dictionary
codon_dict_path = os.path.join(path, 'codon_dict.pkl')

with open(codon_dict_path, 'rb') as file:
    codon_dict = pickle.load(file)
        
# Four letters codes
four_letter_codes = list(master_df['4-letters'].unique())  # Get unique four-letter codes from master_df


'''
Functions:
    > data_prep(self, path)
    > process_sequence(self)
        > modify_TIS_in_frame(self, sequence)
            > find_in_frame(self, sequence)
            > efficiency_level(self, sequence)
        > modify_TIS_out_of_frame(self, sequence)
            > find_out_of_frame(self, sequence)
            > get_aa_key(self, sequence)
    > to_fasta(self, sequence, output_file_name)
    > filtered_sequence_eff(self, efficiency)
    > find_lower_eff_sequence(self, efficiency, seq)
'''

class SemperRecode:
    def __init__(self, user_seq):
        """
        Initializes the SemperRecode object.

        Parameters
        ----------
        user_seq : str
            User-input sequence as the result of .fasta file being parsed by users

        Raises
        ------
        ValueError
            1. If input sequence is blank.
            2. If input sequence contains "u" or "U"

        """
        global master_df
        self.master_df = master_df

        global start_codon
        self.start_codon = start_codon

        # Convert user input sequence into str (string)
        self.seq = str(user_seq)

        # Check if user input sequence is empty
        if len(self.seq) == 0 or self.seq.isspace():
            raise ValueError("No sequence input")
        
        # Raise ValueError if user input sequence contains "u" or "U"
        if "U" in self.seq or "u" in self.seq:
            raise ValueError(f"U or u found in the input sequence {self.seq}")
    
    def process_sequence(self):
        """
        Uses the nucleotide sequence input through the constructor and returns the modified sequence with lower efficiency (if any).
        Loop through modify_TIS_in_frame() and modify_TIS_out_of_frame() to ensure there's no out-of-frame internal AUG and that the
        expressional level is the lowest it could possibly be

        Parameters
        ----------
        None

        Returns
        -------
        replace_sequence : str
            String of modified sequences of the input sequence.

        Raises
        ------
        None

        """

        # Find modified sequence which is returned by modify_TIS_in_frame()
        replace_sequence = self.modify_TIS_in_frame(self.seq)

        return replace_sequence

    def modify_TIS_in_frame(self, sequence):
        '''
        Takes in sequence (str) and return modified sequence with lower TIS efficiency (str)
        by getting all the index location of internal AUG from find_in_frame() and modify the sequence accordingly

        Parameters
        ----------
            sequence : str

        Returns
        -------
            modified_sequence (str)
        '''
        new_seq = list(sequence)
        index = self.find_in_frame(sequence) # Get the indices of in-frame AUG(s)
        df = self.master_df

        '''
        Iterate through the length of the sequence (len(sequence)) and if the index in found in
        index (list of positions from find_in_frame function), find a sequence with lower efficiency,
        else, concatnate the old sequence to the output list (new_seq)
        '''
        for pos in index:
            # Ignore the first and last AUG
            if pos != 0 and pos != len(sequence) - 3:
                internal_TIS_seq = Seq(sequence[pos-6:pos+5])
                aa4 = Seq(sequence[pos-6:pos+6]).translate()

                # Get the current efficiency level 
                current_eff = self.efficiency_level(internal_TIS_seq)

                filtered = df[df['4-letters'] == str(aa4)]
                new_eff = filtered['efficiency'].iloc[0]

                '''
                If a new sequence with a lower efficiency is found,
                replace the sequence in the string. 
                Else, the sequence remain unchanged
                '''
                if(int(new_eff) < current_eff):
                    new_seq[pos-6:pos+6] = filtered["4-codons"].iloc[0]
                else:
                    print(f"No sequence with lower efficiency is found for {internal_TIS_seq}\nConsider mutate/remove the sequence")

        return ''.join(new_seq)

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
        for i in range(0, len(sequence)-2, 3):
            codon = sequence[i:i+3]

            if codon in self.start_codon:
                pos.append(i)

        return pos
    
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

    
    def modify_TIS_out_of_frame(self, sequence):
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
        index = self.find_out_of_frame(sequence) # Get the indices fo out-of-frame AUG(s)
        df = self.master_df

        '''
        Iterate through sequene and modify sequence to get rid of any out-of-frame AUG(s)
        using the index from find_out_of_frame()
        '''

        for pos in index:
            # Break the it down into 2 codons according to its proper codon frame
            first_aa = sequence[pos%3: pos%3 + 2]
            second_aa = sequence[pos%3 + 2 : pos%3 + 4]


            '''
            Find the other codon sequence (if any) which produce the same 
            amino acid as the key using the self.codon_dict which contains
            amino acid's name, all possible codon sequence, and its fraction

            {'A': {'GCC': 0.3975938884126084, 'GCT': 0.2626759965707805,
            'GCA': 0.2301414614526459, 'GCG': 0.1095886535639651}, .......

            Ex: 
                Given: 'N': {'AAC': 0.5184461245612413, 'AAT': 0.4815538754387587},
                Say our first_aa is 'AAC', then we get 'AAT'
            '''
            # first_aa_new_codon, first_aa_new_codon_value = self.get_aa_alternative(first_aa)
            # second_aa_new_codon, second_aa_new_codon_value = self.get_aa_alternative(second_aa)
            
            '''
            Compare which pair of codon (the old and the new one) has the least difference in fraction value
            When the pair is found, replace the old codon with the new codon
            '''

            
            

        return ''.join(new_seq)
    
    def find_out_of_frame(self, sequence):
        """
        Takes in a sequence (str) and returns a list of indices where the out-of-frame AUG codon is found.

        Parameters
        ----------
        sequence : str
            The input sequence to search for AUG codons.

        Returns
        -------
        pos : list
            A list of integers representing the indices of out-of-frame AUG codons in the sequence.

        """
        pos = []
        index = -1

        for codon in self.start_codon:
            index = -1

            while True:
                index = sequence.find(codon, index + 1)
                if index == -1:
                    break

                if index%3 != 0:
                    pos.append(index)
    
        return pos

    def get_aa_alternative(self, original_codon):
        '''
        Takes in original_codon (ex: 'GCC', 'GAG', 'TCA') then return 1. the codon which produce the key amino acid 
        with highest fraction (other than the original codon itself) and 2. The fraction of that codon

        Sample outout
        Input: get_aa_key('GCC')
        Output: 'GCT', 0.2301414614526459

        Parameters
        ----------
            original_codon : str

        Returns
        -------
            codon (str)
            value (int)

        '''
        # Find the corresponding amino acid (ex: 'A', 'R', 'N')
        aa = ""

        for key, value in codon_list.items():
            if original_codon in value:
                aa = key
                break

        # codon = list(codon_dict[aa].keys())
        # index = 0

        # # If the codon is the same as the original one, get to the next one
        # index += 1 if codon == original_codon else 0

        # value = codon_dict[aa][codon[index]]

        # return codon, int(value)

        return aa
    def to_fasta(self, sequence, output_file_name):
        # self.to_fasta(modified_sequences, "sample_file_outputs")
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

        Raises
        ------
        ValueError
            If the sequence is empty.

        '''
        if not sequence:
            raise ValueError("Sequence list is empty, unable to export")

        output = []

        # for seq in enumerate(sequence):
        #     temp = SeqRecord(Seq(seq), id=f"{seq.id}_semper_recode", description='')
        #     output.append(temp)
        '''
        i = index as wel iterated through sequence 

        In each iteration, a new SeqRecord object (Seq object but with additional metadata)
        is being created 

        Metadata:
            id: The sequence original sequence name stored in self.seq_id
            description: None
1
        '''
        output_file = "tests/sample_file/" + output_file_name + ".fasta"

        with open(output_file, 'w') as file:
            SeqIO.write(sequence, file, 'fasta')
        
        '''
        Output: 
        >Sequence1
        ATCGATCGATCG

        >Sequence2
        GCTAGCTAGCTA

        >Sequence3
        TGCAATGCAATG
        '''
        
    ''' 
    =======================================================================================
                                        
                                        Exploratory analysis

    =======================================================================================
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

    ''' 
    =======================================================================================
                                        
                                        Exploratory analysis
                                        
    =======================================================================================
    '''