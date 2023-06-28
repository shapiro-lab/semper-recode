'''Install biopython and pytest using pip if not already installed

  1. pip install biopython
  2. pip install pytest

'''

import pandas as pd
import numpy as np
import re
from collections import Counter
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
from itertools import combinations_with_replacement
from itertools import product
from Bio.Seq import Seq
import ishaan_utils

class BioSeq:
    def __init__(self):
        self.four_letter_codes = []
        self.pwm_df, self.aa_chart, self.codon_usage, self.codon_dict = self.load_data()

    def load_data(self):
        # Position weight matrix
        pwm_df = pd.read_csv('/content/drive/MyDrive/Caltech Connection/package_materials/noderer-dinuc-pwm.txt', delimiter='\t', skiprows=1)
        pwm_df['sequence'] = pwm_df['sequence'].str.replace('U', 'T')
        pwm_df = pwm_df.sort_values(by=['efficiency'], ascending=True)
        pwm_df = pwm_df.reset_index(drop=True)

        # Amino acids (full-name, 3-letter, 1-letter)
        aa_chart = pd.read_csv('/content/drive/MyDrive/Caltech Connection/package_materials/amino-acid-codes.csv')
        aa_chart['1-letter'] = np.where(aa_chart['1-letter'] == 'STOP', '*', aa_chart['1-letter'])
        one_letter_codes = list(aa_chart['1-letter'].unique())


        # Codon chart (merge codon chart with aa_chart)
        codon_usage = pd.read_csv('/content/drive/MyDrive/Caltech Connection/package_materials/nuclear_codon_statistics.tsv', delimiter='\t', skiprows=8)
        codon_usage = codon_usage.merge(aa_chart, how='right', left_on='Amino acid', right_on='3-letter')
        codon_usage = codon_usage.sort_values(by=['1-letter', 'Fraction'], ascending=[True, False])
        codon_usage = codon_usage.reset_index(drop=True)

        # Create dictionary of all codon that produce the same amino acid ex: {'A' : {'GCC', 'GCT', 'GCA', 'GCG'}}
        codon_dict = {}
        for aa in one_letter_codes:
        codon_dict[aa] = dict(codon_usage.loc[codon_usage['1-letter'] == aa, ['CODON', 'Fraction']].values)

        # List of start codons
        start_codon = ['AUG', 'ATG', 'GUG', 'UUG']

        return pwm_df, aa_chart, codon_usage, codon_dict


    def get_nuc_TIS_combos(self, aa_combo, pwm_samples_df):
        """
        Takes in a four-character string representative of amino acids and returns a list of the possible
        nucleotide sequences from -6 to +5 with all possible degenerate codons used along with the predicted
        efficiency from the PWM.

        Args:
            aa_combo (str): A four-character string representing amino acids.
            pwm_samples_df (pd.DataFrame): A DataFrame containing the PWM samples.

        Returns:
            pd.DataFrame: A DataFrame merging the combinations of nucleotide sequences with the PWM samples.

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
        using the `get_nuc_TIS_combos` method

        Args:
            None

        Raises:
            ValueError: If duplicates are found in the generated list of amino acids combinations.

        Returns:
            None

        Examples:
            >>> obj = MyClass()
            >>> obj.explore_aa_combinations(efficiency_threshold=40)

        See Also:
            - :func:`get_nuc_TIS_combos`: Get nucleotide TIS (Translation Initiation Site) combinations.
            - :class:`pd.DataFrame`: pandas DataFrame object for data manipulation.
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

        self.four_letter_codes = list(master_df['4-letters'].unique())  # Get unique four-letter codes from master_df

    def troublesome_sequence(self, efficiency, df):
        """
        Takes in efficiency level input and finds the sequences in df (dataframe) by filtering, then displays all the possible sequences that have an efficiency lower than the given value and the percentage of found and missing amino acid sequences in the form of a dataframe.

        Args:
            efficiency (float): The efficiency level to filter the dataframe.
            df (pd.DataFrame): The dataframe containing the sequences and efficiencies.

        Returns:
            pd.DataFrame: A dataframe containing the filtered sequences.

        """
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

# Example usage
explorer = BioSeq()
