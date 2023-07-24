import os
import re
import pickle
import warnings
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

'''
Global variables declaration

Loads data required for the analysis. Reads various data files and initializes necessary attributes.

Raises
------
FileNotFoundError
    If any of the required data files are not found.

'''
START_CODON = ['ATG']
PATH = "data/" 

# Master dataframe
MASTER_DF_PATH = os.path.join(PATH, 'master_df_os_2023.csv')
if not os.path.exists(MASTER_DF_PATH):
    raise FileNotFoundError(f"Master dataframe file not found: {MASTER_DF_PATH}")
        
MASTER_DF = pd.read_csv(PATH + 'master_df_os_2023.csv')

# Load codon list
CODON_LIST_PATH = os.path.join(PATH, 'codon_list.pkl')
if not os.path.exists(CODON_LIST_PATH):
    raise FileNotFoundError(f"Codon list (.pkl) file not found: {CODON_LIST_PATH}")

with open(CODON_LIST_PATH, 'rb') as file:
    CODON_LIST = pickle.load(file)

if CODON_LIST == {}:
    raise ValueError("codon_list is empty")

# Load codon dictionary
CODON_DICT_PATH = os.path.join(PATH, 'codon_dict.pkl')

with open(CODON_DICT_PATH, 'rb') as file:
    CODON_DICT = pickle.load(file)
        
# Four letters codes
FOUR_LETTERS_CODE = list(MASTER_DF['4-letters'].unique())  # Get unique four-letter codes from master_df

'''
Functions:
    > data_prep(self, path)
    > process_sequence(self)
        > modify_TIS_in_frame(self, sequence)
            > find_in_frame(self, sequence)
            > efficiency_level(self, sequence)
        > modify_TIS_out_of_frame(self, sequence)
            > find_out_of_frame_list_list(self, sequence)
            > find_key(self, codon)
            > get_aa_key(self, original_codon)
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

        Raises
        ------
        ValueError
            1. If input sequence is blank.
            2. If input sequence contains "u" or "U"

        """
        self.master_df = MASTER_DF
        self.start_codon = START_CODON

        # Convert user input sequence into str, remove blank space, and capitalize the string
        self.seq = str(user_seq).replace(" ", "").upper()

        # Check if user input sequence is empty
        if len(self.seq) == 0 or self.seq.isspace():
            raise ValueError("No sequence input")
        
        # Raise ValueError if user input sequence contains "u" or "U"
        if "U" in self.seq or "u" in self.seq:
            raise ValueError(f"U or u found in the input sequence {self.seq}")
        
        # Raise ValueError if the sequence contains any character other than A, T, C, G
        pattern = r'[^ATCG]'
        if re.search(pattern, self.seq):
            raise ValueError("Invalid character found in the string.")
    
        self.error_list = []
    
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

        replace_sequence = self.seq
        round = 0

        # Loop through modify_TIS_in_frame() and modify_TIS_out_of_frame() while there is out-of-frame AUGs or the loop hasn't reached 10 times
        while self.find_out_of_frame_list_list(replace_sequence) != [] and round < 10:
            # Find modified sequence which is returned by modify_TIS_in_frame()
            temp = self.modify_TIS_in_frame(self.seq)
            replace_sequence = self.modify_TIS_out_of_frame(temp)
            round += 1

        # Return modified sequence along with error list
        return replace_sequence, self.error_list

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
        Convert the sequence to a list (new_seq) so that's it's mutable and when we iterate through the sequence, 
        we can replace the TIS sequence with the modified sequence with lower efficiency according to the index
        of AUG return by find_in_frame()

        This way, we'll iterate through the updated list of sequence (new_seq) with modified sequence
        '''
        for pos in index:
            # Ignore the first and last AUG
            if pos != 0 and pos != len(sequence) - 3:

                if pos == len(sequence) - 3:
                    self.raiseWarning(f"There's an AUG at the end of the sequence which cannot be modified")


                sub_1 = ''.join(new_seq[pos-6:pos+5])
                internal_TIS_seq = Seq(sub_1)

                sub_2 = ''.join(new_seq[pos-6:pos+6])
                aa4 = Seq(sub_2).translate()

                # Get the current efficiency level 
                current_eff = self.efficiency_level(internal_TIS_seq)

                filtered = df[df['4-letters'] == str(aa4)]
                new_eff = filtered['efficiency'].iloc[0]

                '''
                Replace the designated TIS sequence with the modified sequence
                If a new sequence with a lower efficiency is not found,
                print a message telling the user to consider mutate/remove the sequence
                '''
                new_seq[pos-6:pos+6] = filtered["4-codons"].iloc[0]
                
                if(int(new_eff) == current_eff):
                    self.raiseWarning(f"No sequence with lower efficiency is found for {internal_TIS_seq} at [{pos}], consider mutate/remove the sequence")

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
        sequence : str
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
        Takes in sequence (str) and return modified sequence
        with lower TIS efficiency (str)

        Parameters
        ----------
            sequence : str

        Returns
        -------
            modified_sequence (str)
        '''
        new_aa = ''
        new_seq = list(sequence)
        index = self.find_out_of_frame_list_list(sequence) # Get the indices of out-of-frame AUG(s)

        '''
        Iterate through sequene and modify sequence to get rid of any out-of-frame AUG(s)
        using the index from find_out_of_frame_list_list()
        '''

        for pos in index:
            # Break the it down into 2 codons according to its proper codon frame
            start = pos - pos%3
            first_aa = sequence[start : start + 3]
            second_aa = sequence[start + 3 : start + 6]

            # Find the amino acid of the codon
            first_aa_key = self.return_key(first_aa)
            second_aa_key = self.return_key(second_aa)

            # Case 1: xAT Gxx - only the second codon can be modified
            if first_aa[1:3] == 'AT' and second_aa[0] == 'G':
                # Since there's no other codons that starts with G besides 'D','E','V','A','G', we'll always modify the first codon
                while new_aa[1:3] != 'AT':
                    new_aa = self.get_alternative_codon(first_aa_key)

                if new_aa != '':
                    new_seq[start : start + 3] = new_aa # Replace the first codon

            # Case 2: xxA TGx
            elif first_aa[2] == 'A' and second_aa[0:2] == 'TG': # There's no codon that always ends with A
                '''
                Compare which codon (the old and the new one) has the least difference in fraction value
                When the pair is found, replace the old codon with the new codon
                '''
                # Get the fraction value of the first and second old codon
                first_old_val = CODON_DICT[first_aa_key][first_aa]
                second_old_val = CODON_DICT[second_aa_key][second_aa]
                
                # Get the codons and fraction values of the two new codons
                first_new_codon, first_new_val = self.get_alternative_codon(first_aa)
                second_new_codon, second_new_val = self.get_alternative_codon(second_aa)

                # Find the difference between each pair of codons (old and new)
                first_codon_diff = abs(first_new_val - first_old_val)
                second_codon_diff = abs(second_new_val - second_old_val)

                '''
                Replace the codon with the least fraction difference
                i.e. if the fraction difference between the old and the new amino acids of first aa is greater,
                replace the second codon with the new codon
                '''
                
                if second_codon_diff > first_codon_diff:
                    new_seq[start : start + 3] = first_new_codon # Replace the first codon
                else: 
                    new_seq[start + 3 : start + 6] = second_new_codon # Replace the second codon

        return ''.join(new_seq)      

    def return_key(self, codon):
        '''
        Takes in codon sequence and return the key (str) which is the abbreviation of the codon

        Parameters
        ----------
            codon : str

        Returns
        -------
            key (str)

        '''
        key = ""

        for char, value in CODON_LIST.items():
            if codon in value:
                key = char
                break

        return key
    
    def find_out_of_frame_list(self, sequence):
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

    def find_out_of_frame_index(self, sequence):
        """
        Takes in a sequence (str) and returns an integer of index where the out-of-frame AUG codon is found.

        Parameters
        ----------
        sequence : str
            The input sequence to search for AUG codons.

        Returns
        -------
        pos : int
            An integer representing the index of out-of-frame AUG codons in the sequence.

        """
        index = -1

        for codon in self.start_codon:
            index = -1

            while True:
                index = sequence.find(codon, index + 1)
                if index == -1:
                    break

                if index%3 != 0:
                    return index
    
        return None

    def get_alternative_codon(self, original_codon):
        '''
        Takes in original_codon (ex: 'GCC', 'GAG', 'TCA') then return:
            1. the codon which produce the key amino acid with highest fraction (other than the original codon itself)
            2. The fraction of that codon

        Sample code: get_aa_key('GCC')
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
        aa = self.return_key(original_codon)

        # Get all the possible codons which product the key (ex: ['CAC': 0.5742015169781323, 'CAT': 0.4257984830218677])
        codons = list(CODON_DICT[aa].keys())
        index = 0

        # If the codon is the same as the original one and there's an alternative, get to the next one
        if codons[0] == original_codon and len(codons) > 1:
            index += 1

        # The fraction of the new codon
        value = CODON_DICT[aa][codons[index]]

        return codons[index], value
    
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

        Raises
        ------
        ValueError
            If the sequence is empty.

        '''
        if not sequence:
            raise ValueError("Sequence list is empty, unable to export")

        output = []

        for i, seq in enumerate(sequence):
            temp = SeqRecord(Seq(seq), id=f'Seq{i+1}', description='')
            output.append(temp)
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
        
    

    def raiseWarning(self, error):
        '''
        Raise error and append errors in error_list using 
        
        Format "error" (i - start index, j -  last index) so user can parse the list and get the index if needed
        Ex: AUG at [2] cannot be modified. Consider mutate/remove the sequence. (0, 1)

        Parameters
        ----------
            error : f string

        Returns
        -------
            warning

        '''
        error_message = error
        self.error_list.append(error_message) # Append error message to error_list
        warnings.warn(error_message)
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