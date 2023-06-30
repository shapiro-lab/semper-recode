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
    def __init__(self, user_seq = None):
        """
        Initializes the SemperRecode object.

        Parameters
        ----------
        user_seq : str, optional
            User-defined sequence, by default None. In case user want to tune single sequence instead of importing a whole FASTA file

        Returns
        -------
        None

        Raises
        ------
        None
        """
        self.pwm_df = None
        self.aa_chart = None
        self.codon_usage = None
        self.start_codon = None
        self.codon_dict = None
        self.four_letter_codes = None
        self.master_df = None
        self.start_codon = None
        self.load_data()
        
        self.seq = Seq(user_seq) if user_seq is not None else 'ATGCTGACGGTAUGGACTTACCTGTATGCGTGCTAAATGCTAAGGCTGGTGCCGACCGGACCGTTGGGAGCGCTGTTGACCGGATGCTAAAGGGCCCGAGTCTTGTAGTACCGGACTTAAATGCGTTGTTTGACACCTGTT'

    def load_data(self):
        """
        Loads data required for the analysis. Reads various data files and initializes necessary attributes.

        Returns
        -------
        None

        Raises
        ------
        FileNotFoundError
            If any of the required data files are not found.

        """
        # Position weight matrix
        self.pwm_df = pd.read_csv('../data/pwm_df.csv') # Position Weight Matrix

        # Amino acids (full-name, 3-letter, 1-letter)
        self.aa_chart = pd.read_csv('../data/aa_chart.csv')

        # Codon chart
        self.codon_usage = pd.read_csv('../data/codon_usage.csv')

        # List of start codons
        self.start_codon = ['AUG', 'ATG', 'GUG', 'UUG']

        # Master dataframe
        self.master_df = pd.read_csv('../data/master_df_os_2023.csv')
        
        # Four letters codes
        with open('../data/four_letter_codes.pkl', 'rb') as file:
            self.four_letter_codes = pickle.load(file)

        # Dictionary of all codon that produce the same amino acid ex: {'A' : {'GCC', 'GCT', 'GCA', 'GCG'}}
        with open('../data/codon_dict.pkl', 'rb') as file:
            self.codon_dict = pickle.load(file)

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

    def get_nuc_TIS_combos(self, aa_combo, pwm_samples_df):
        """
        Takes in a four-character string representative of amino acids and returns a list of the possible
        nucleotide sequences from -6 to +5 with all possible degenerate codons used along with the predicted
        efficiency from the PWM.

        Parameters
        ----------
        aa_combo : str
            A four-character string representing amino acids.
        pwm_samples_df : pd.DataFrame
            A DataFrame containing the PWM samples.

        Returns
        -------
        pd.DataFrame
            A DataFrame merging the combinations of nucleotide sequences with the PWM samples.
        
        """

        aa_list_0 = list(self.codon_dict[aa_combo[0]].keys())
        aa_list_1 = list(self.codon_dict[aa_combo[1]].keys())
        aa_list_2 = list(self.codon_dict[aa_combo[2]].keys())
        aa_list_3 = list(self.codon_dict[aa_combo[3]].keys())

        aa_lists = [aa_list_0, aa_list_1, aa_list_2, aa_list_3]

        # removes one nucleotide at +6 position
        combinations = [''.join(p) for p in itertools.product(*aa_lists)]
        tis_combinations = [''.join(p)[0:-1] for p in itertools.product(*aa_lists)]

        combo_df = pd.DataFrame({'4-codons':combinations, '4-letters':aa_combo, 'tis-sequence': tis_combinations})

        return combo_df.merge(pwm_samples_df, how='left', left_on = 'tis-sequence', right_on = 'sequence')

    def generate_aa_combinations(self):
        """
        Generate all possible amino acid combinations and generate a list of unique four-letter codes.
        The combinations are formed by appending the amino acid sequence with 'M' as the 3rd character.
        After generating, the method retrieves nucleotide TIS combinations for each four-letter code
        using the `get_nuc_TIS_combos` method.

        Parameters
        ----------
        None

        Raises
        ------
        ValueError
            If duplicates are found in the generated list of amino acid combinations.

        Returns
        -------
        List of four-letter codes from master_df

        Examples
        --------
        >>> obj = MyClass()
        >>> obj.explore_aa_combinations()

        See Also
        --------
        get_nuc_TIS_combos : Get nucleotide TIS (Translation Initiation Site) combinations.
        pd.DataFrame : pandas DataFrame object for data manipulation.

        """
        aa_list = list(self.aa_chart['1-letter'].unique())  # Get unique 1-letter amino acid codes
        aa_TIS_list = []  # Initialize list for amino acid combinations with 'M'

        for combo in product(aa_list, repeat=3):  # Generate all combinations of 3 amino acids
            aa_TIS_list.append(combo[0] + combo[1] + 'M' + combo[2])  # Append combination with 'M'

        dup = pd.Series(aa_TIS_list).duplicated()  # Check for duplicates in the combination list

        if dup.any():  # If duplicates are found
            raise ValueError("Duplicates found in the list!")

        master_df = pd.DataFrame(columns=['4-codons', '4-letters'])  # Initialize a master DataFrame

        for item in aa_TIS_list:  # Iterate over each amino acid combination
            df = self.get_nuc_TIS_combos(item, self.pwm_df)  # Get nucleotide TIS combos for the combination
            master_df = pd.concat([master_df, df])  # Concatenate the resulting DataFrame with master_df

        master_df = master_df.sort_values(by=['efficiency'], ascending=True).drop_duplicates(keep='first')  # Sort and drop duplicates in master_df
        master_df = master_df.reset_index(drop=True)  # Reset index of master_df

        return list(master_df['4-letters'].unique())  # Get unique four-letter codes from master_df

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


        Returns
        -------
        Bio.Seq
            Nucleotide sequence with lower efficiency level
        """

        possible_seq_df = self.filtered_sequence(efficiency, self.master_df)


