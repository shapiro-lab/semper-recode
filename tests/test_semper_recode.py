"""Tests for `semper_recode` package."""

from semper_recode.semper_recode import SemperRecode # Import SemperRecode class
import pytest
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

@pytest.fixture
def obj():
    return SemperRecode("ATGCTGATGCTAATGCGTACGTAGCTAA")

# ================== TEST CONSTRUCTOR ==================

def test_constructor_valid_sequence():
    """
    Purpose: Verify the functionality of the `SemperRecode` class constructor with a valid gene sequence.
    Goal: Ensure that the constructor initializes a `SemperRecode` object properly with the user-provided sequence.
    """
    user_seq = "ATGCATGC"
    recode = SemperRecode(user_seq)
    assert isinstance(recode, SemperRecode)
    assert recode.seq == user_seq

def test_constructor_empty_sequence():
    """
    Purpose: Verify the error handling of the `SemperRecode` class constructor with an empty gene sequence.
    Goal: Ensure that the constructor raises a `TypeError` when an empty sequence is provided as input.
    """
    with pytest.raises(TypeError):
        SemperRecode()

def test_constructor_whitespace_sequence():
    """
    Purpose: Verify the error handling of the `SemperRecode` class constructor with a whitespace-containing gene sequence.
    Goal: Ensure that the constructor raises a `ValueError` when a sequence containing only whitespace is provided as input.
    """
    with pytest.raises(ValueError):
        SemperRecode("    ")

def test_constructor_blank_sequence():
    """
    Purpose: Verify the error handling of the `SemperRecode` class constructor with a blank gene sequence.
    Goal: Ensure that the constructor raises a `ValueError` when a blank sequence (empty string) is provided as input.
    """
    with pytest.raises(ValueError):
        SemperRecode("")

# ================== TEST FIND_IN_FRAME ==================

def test_find_in_frame():
    """
    Purpose: Verify the functionality of the `find_in_frame` method in the SemperRecode class.
    Goal: Ensure that the method correctly identifies in-frame AUGs and returns their positions in the given sequence as list.
    """

    # In-frame AUGs
    sample1 = "ATGCTGATGCTAATGCGTACGTAGCTAA" # Have one extra nucleotide
    obj1 = SemperRecode(sample1)
    assert obj1.find_in_frame(sample1) == [0,6,12] 
    
    # Out-of-frame AUGs
    sample2 = "ATTGTTTATGGCCATTTGGAT"
    obj2 = SemperRecode(sample2)
    assert obj2.find_in_frame(sample2) == []

    sample3 = "AAAATG"
    obj3 = SemperRecode(sample3)
    assert obj3.find_in_frame(sample3) == [3]

    sample4 = "ATGGTGATCAAGAACATCCAGGTGTTCTTTATGAAAACCATCAGCAACCGGTCCATCAGCCGGGCCAAGATCAGCACCATGCCCAGACCTATG"
    print(len(sample4))
    obj4 = SemperRecode(sample4)
    assert obj4.find_in_frame(sample4) == [0, 30, 78, 90]

# ================== TEST FIND_OUT_OF_FRAME ==================

def test_find_out_of_frame_list(obj):
    """
    Purpose: Verify the functionality of the `find_out_of_frame_list` method in the SemperRecode class.
    Goal: Ensure that the method correctly identifies out-of-frame AUGs and even when there's extra nucleotide 
          in the sequence & in-frame AUG present and returns their positions in the given sequence as list.
    """
    assert obj.find_out_of_frame_list("AAAATG") == [] # No out-of-frame AUG
    assert obj.find_out_of_frame_list("ATGAATGAT") == [4]
    assert obj.find_out_of_frame_list("ATTGTTTATGGCCATTTGGAT") == [7]
    assert obj.find_out_of_frame_list("AATGTTTATGGCATGTTGGAT") == [1,7]
    assert obj.find_out_of_frame_list("AATGTATATGGCATGTTGGAT") == [1, 7]
    assert obj.find_out_of_frame_list("AATGTTTATGGCATGTTGGATATG") == [1, 7]
    assert obj.find_out_of_frame_list("ATTGTTTAGGGCCATTTGGAT") == [] # No AUG
    assert obj.find_out_of_frame_list("AATGTTTATGGCATGTTGAT") == [1, 7] # 2 extra nucleotide (partial nucleotide)
    assert obj.find_out_of_frame_list("AATGTTTATGGCCTGTTGGATATG") == [1, 7] # In-frame and out-of-frame sequence

# ================== TEST FIND_OUT_OF_FRAME ==================

def test_find_out_of_frame_index(obj):
    """
    Purpose: Verify the functionality of the `find_out_of_frame_list` method in the SemperRecode class.
    Goal: Ensure that the method correctly identifies out-of-frame AUGs and even when there's extra nucleotide 
          in the sequence & in-frame AUG present and returns their positions in the given sequence as list.
    """
    assert obj.find_out_of_frame_index("AAAATG") == None # No out-of-frame AUG
    assert obj.find_out_of_frame_index("ATGAATGAT") == 4
    assert obj.find_out_of_frame_index("ATTGTTTATGGCCATTTGGAT") == 7
    assert obj.find_out_of_frame_index("AATGTTTATGGCATGTTGGAT") == 1
    assert obj.find_out_of_frame_index("AATGTATATGGCATGTTGGAT") == 1
    assert obj.find_out_of_frame_index("AATGTTTATGGCATGTTGGATATG") == 1
    assert obj.find_out_of_frame_index("ATTGTTTAGGGCCATTTGGAT") == None # No AUG
    assert obj.find_out_of_frame_index("AATGTTTATGGCATGTTGAT") == 1 # 2 extra nucleotide
    assert obj.find_out_of_frame_index("AATGTTTATGGCCTGTTGGATATG") == 1 # In-frame and out-of-frame sequence


# ============= TEST EFFICIENCY_LEVEL =============

def test_efficiency_level_with_exist_sequence(obj):
    """
    Purpose: Test the `efficiency_level` method in the SemperRecode class with existing sequences.
    Goal: Verify that the efficiency level is calculated correctly for the given sequences.
    """

    assert(obj.efficiency_level("CATTGTATGCT") == 12)
    assert(obj.efficiency_level("AAACGTATGGG") == 62)
    assert(obj.efficiency_level("ATTAGAATGAC") == 95)
    assert(obj.efficiency_level("CCTATGATGTA") == 102)
    assert(obj.efficiency_level("TTCATCATGCA") == 150)
    
def test_efficiency_level_with_non_exist_sequence(obj):
    """
    Purpose: Test the `efficiency_level` method in the SemperRecode class with non-existing sequences.
    Goal: Ensure that a ValueError is raised when a sequence doesn't exist in the master_df.
    """

    # ValueError is expected to raise as the sequence doesn't exist in master_df
    with pytest.raises(ValueError):
        obj.efficiency_level("CCCTGTACGCT")
        obj.efficiency_level("CCCTTTAAACT")

# ============= TEST MODIFY_TIS_IN_FRAME =============

def test_modify_TIS_with_exist_sequence(obj):
    """
    Purpose: Test the `modify_TIS_in_frame` method in the SemperRecode class with existing sequences.
    Goal: Ensure that the TIS modification is performed correctly for the given sequences.
    """

    assert obj.modify_TIS_in_frame("CACTGCATGTTA") == "CATTGTATGCTG"  # Changing efficiency from 34 -> 12
    assert obj.modify_TIS_in_frame("ATGCACTGCATGTTA") == "ATGCATTGTATGCTG"  # Changing efficiency from 34 -> 12 + the first AUG remain the same as expected
    assert obj.modify_TIS_in_frame("CACTGCATGTTAATG") == "CATTGTATGCTGATG"  # Changing efficiency from 34 -> 12 + the last AUG remain the same as expected
    assert obj.modify_TIS_in_frame("AATGAAATGCTG") == "AATGAGATGCTG"  # Changing efficiency from 82 -> 80

def test_modify_TIS_with_non_exist_sequence(obj):
    """
    Purpose: Test the `modify_TIS_in_frame` method in the SemperRecode class with non-existing sequences.
    Goal: Ensure that no modification is performed when a sequence doesn't contain any AUGs.
    """

    # Tammy: Use this test cases when testing process_sequence() to see if both modify_TIS_in_frame and modify_TIS_out_of_frame works

    assert obj.modify_TIS_in_frame("ATGCATGTTA") == "ATGCATGTTA"  # No AUG in the sequence, so no modification expected
    assert obj.modify_TIS_in_frame("ATGCACTGCATG") == "ATGCACTGCATG"  # No AUG in the internal region, so no modification expected   

def test_modify_TIS_with_more_than_one_AUG(obj):
    """
    Purpose: Test the `modify_TIS_in_frame` method in the SemperRecode class with sequences that has more than 1 TIS sequences.
    Goal: Ensure that the modification is performed correctly even when there's 2 or more AUG being present in the same TIS sequence
    """

    assert obj.modify_TIS_in_frame("ATGCATGTTATGCATATGCAC") == "ATGCATGTTATGCATATGCAC" # There's 2 AUG in the second TIS sequence
    '''
    ATG | CAT | GTT | ATG | CAT | ATG | CAC
    in-frame AUG = [9, 15]

    pos = 9:
        orig_TIS = CATGTTATGCAT (HVMH | eff: 64)
        new_TIS  = CATGTTATGCAC (HVMH) | eff: 64)
        no change
        -> seq = orig_TIS = ATGCATGTTATGCATATGCAC

    pos = 15:
        orig_TIS = ATGCATATGCAC (MHMH | eff: 65)
        new_TIS  = ATGCATATGCAT (MHMH | eff: 65)
        no change
        -> seq = ATGCATGTTATGCATATGCAC

    => modified_seq = ATGCATGTTATGCATATGCAC
    '''
    assert obj.modify_TIS_in_frame("ATGGTCTCGATGAATTATTGCATGATG") == "ATGGTGTCTATGAATTATTGTATGATG"
    '''
    ATG | GTC | TCG | ATG | AAT | TAT | TGC | ATG | ATG
    in-frame AUG = [9, 21]

    pos = 9:
        orig_TIS =  GTCTCGATGAAT (VSMN | eff: 86)
        new_TIS  =  GTGTCTATGAAT (VSMN | eff: 53)
        changed GTCTCGATGAAT -> GTGTCTATGAAT
        -> seq = ATGGTGTCTATGAATTATTGCATGATG

    pos = 21:
        orig_TIS =  TATTGCATGATG (YCMM | eff: 17)
        new_TIS  =  TATTGTATGATG (MHMH | eff: 15)
        changed TATTGCATGATG -> TATTGTATGATG
        -> seq = ATGGTGTCTATGAATTATTGTATGATG

    => modified_seq = ATGGTGTCTATGAATTATTGTATGATG
    '''

# ============= TEST GET_AA_ALTERNATIVE =============

def test_get_alternative_codon(obj):
    """
    Purpose: Test the `get_alternative_codon` method in the SemperRecode class.
    Goal: Ensure that the alternative codon and its fraction are returned correctly when given the codon.
    """

    aa, val = obj.get_alternative_codon("AGA")
    assert aa == "AGG" and val == 0.2134612754706659

    aa, val = obj.get_alternative_codon("TGA")
    assert aa == "TAG" and val == 0.2235323759133283

    aa, val = obj.get_alternative_codon("AAT")
    assert aa == "AAC" and val == 0.5184461245612413

    aa, val = obj.get_alternative_codon("CCA")
    assert aa == "CCC" and val == 0.3214351373711047

    aa, val = obj.get_alternative_codon("CTG")
    assert aa == "CTC" and val == 0.1917680920552213

    aa, val = obj.get_alternative_codon("TGG") # There's only 1 codon that produce W (Trp)
    assert aa == "TGG" and val == 1.0

    aa, val = obj.get_alternative_codon("ATG") # There's only 1 codon that produce M (Met)
    assert aa == "ATG" and val == 1.0

    
# ============= TEST RETURN_KEY =============

def test_return_key(obj):
    """
    Purpose: Test the `return_key` method in the SemperRecode class.
    Goal: Ensure that the key (ex: A, R, N, P) is returned correctly when given the codon.
    """
    assert obj.return_key("CCA") == "P"
    assert obj.return_key("CTG") == "L"
    assert obj.return_key("GAG") == "E"
    assert obj.return_key("GCA") == "A"
    assert obj.return_key("TCG") == "S"
    assert obj.return_key("GAA") == "E"

    assert obj.return_key("TGG") == "W"
    assert obj.return_key("TAG") == "*"
    assert obj.return_key("CGT") == "R"
    assert obj.return_key("AAA") == "K"

# ============= TEST MODIFY_TIS_OUT_OF_FRAME =============

def test_modify_TIS_out_of_frame(obj):
    assert obj.modify_TIS_out_of_frame("TTGAATGATTTG") == "TTGAACGATTTG" # Out-of-frame AUG
    assert obj.modify_TIS_out_of_frame("GGGATGAATGAT") == "GGGATGAACGAT" # In-frame AUG at the start 
    assert obj.modify_TIS_out_of_frame("AGCAATGATATG") == "AGCAACGATATG" # In-frame AUG at the end 
    assert obj.modify_TIS_out_of_frame("ATGAATGATATG") == "ATGAACGATATG" # Out-of-frame AUG
    assert obj.modify_TIS_out_of_frame("GCCCAATGCGCC") == "GCCCAGTGCGCC" # Out-of-frame AUG

    # Consecutive out-of-frame AUG
    assert obj.modify_TIS_out_of_frame("CATGATGATGATGCC") == "CACGACGACGACGCC" # 4 consecutive type 1 out-of-frame AUG
    assert obj.modify_TIS_out_of_frame("GCATGATGATGATGG") == "GCCTGATAGTAGTGG" # 4 consecutive type 2 out-of-frame AUG
    assert obj.modify_TIS_out_of_frame("CTTCCATGTTATGGG") == "CTTCCCTGTTACGGG" # 2 consecutive type 2 + type 1 out-of-frame AUG

# ============= TEST PROCESS_SEQUENCE =============

def test_single_seq():
    """
    Purpose: Test the `process_sequence` method in the SemperRecode class.
    Goal: Verify the functionality of the SemperRecode class by processing a sample sequence and comparing the result with the ground truth (manually modified sequence).
    """

    seq = "ATGGTGTCCAAGGGCGAGGAACTGTTCACCGGCGTGGTGCCCATCCTGGTGGAACTGGATGGCGACGTGAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAAGGCGACGCCACATACGGAAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCTTGGCCTACCCTCGTGACCACACTGACCTACGGCGTGCAGTGCTTCGCCAGATACCCCGACCACATGAAGCAGCACGATTTCTTCAAGAGCGCCATGCCCGAGGGCTAC"
    my_obj = SemperRecode(seq)
    expected = "ATGGTGTCCAAGGGCGAGGAACTGTTCACCGGCGTGGTGCCCATCCTGGTGGAACTGGACGGCGACGTGAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAAGGCGACGCCACATACGGAAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCTTGGCCTACCCTCGTGACCACACTGACCTACGGCGTGCAGTGCTTCGCCAGATACCCCGATCATATGAAACAGCACGATTTCTTCAAGAGTGCTATGCCTGAGGGCTAC"
    assert my_obj.process_sequence()[0] == expected

def test_fasta_file():
    """
    Purpose: Test the overall implementation in the SemperRecode class (from file parsing to exporting).
    Goal: Verify the functionality of the SemperRecode class when implementing in real world scenerio.
    """
    modified_seq = "" # Used to store modified sequence with metadata (Seq object)
    output = [] # Used to store all the Seq object of modified sequence for file exporting
    groud_truth = ["ATGGTGTCCAAGGGCGAGGAACTGTTCACCGGCGTGGTGCCCATCCTGGTGGAACTGGACGGCGACGTGAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAAGGCGACGCCACATACGGAAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCTTGGCCTACCCTCGTGACCACACTGACCTACGGCGTGCAGTGCTTCGCCAGATACCCCGATCATATGAAACAGCACGATTTCTTCAAGAGTGCTATGCCTGAGGGCTACGTGCAGGAACGGACCATCTTCTTCAAGGACGACGGCAACTACAAGACAAGAGCCGAAGTGAAGTTCGAGGGCGACACCCTCGTGAACCGGATCGAGCTGAAGGGCATCGACTTCAAAGAGGACGGCAACATCCTGGGCCACAAGCTGGAGTACAACTACAACAGCCACAAGGTGTACATCACCGCCGACAAGCAGAAAAACGGCATCAAAGTGAACTTCAAGACCCGGCACAACATCGAGGACGGCAGCGTGCAGCTGGCCGACCACTACCAGCAGAACACCCCCATCGGAGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACACAAAGCGCCCTGAGCAAGGACCCCAACGAGAAGCGGGATCATATGGTTCTGCTGGAATTTGTGACCGCCGCTGGCATCACCCTTGGTATGGACGAGCTGTACAAGTGA"
                   , "ATGGCGGAAGGCTCCGTCGCCAGGCAGCCTGACCTCTTGACCTGCGAACACGAAGAGATCCACCTCGCCGGCTCGATCCAGCCGCACGGCGCGCTTCTGGTCGTCAGCGAACACGATCATCGCGTCATCCAGGCCAGCGCCAACGCCGCGGAATTTCTGAATCTCGGAAGCGTACTCGGCGTTCCGCTCGCCGAGATCGACGGCGATCTGTTGATCAAGATCCTGCCGCATCTCGATCCCACCGCCGAAGGGATGCCGGTCGCGGTGCGCTGCCGGATCGGCAATCCCTCTACGGAGTACTGCGGGTTGATGCATCGGCCTCCGGAAGGCGGGCTGATCATCGAACTCGAACGTGCCGGCCCGTCGATCGATCTGTCAGGCACGCTGGCGCCGGCGCTGGAGCGGATCCGCACGGCGGGTTCACTGCGCGCGCTGTGCGACGACACCGTGCTGCTGTTTCAGCAGTGCACCGGCTACGACCGTGTTATGGTGTATCGTTTCGACGAGCAAGGCCACGGCCTGGTATTCTCCGAGTGCCACGTGCCTGGGCTCGAATCCTATTTCGGCAACCGCTATCCGTCGTCGACTGTCCCACAGATGGCGCGGCAGCTGTACGTGCGGCAGCGCGTCCGCGTGCTGGTCGACGTCACCTATCAGCCGGTGCCGCTGGAGCCGCGGCTGTCGCCGCTGACCGGGCGCGATCTTGATATGAGTGGCTGCTTCCTGCGGTCTATGAGTCCGTGCCATCTGCAGTTCCTGAAGGATATGGGCGTGCGCGCCACCCTGGCGGTGTCGCTGGTGGTCGGCGGCAAGCTGTGGGGCCTGGTTGTCTGTCACCATTATCTGCCGCGCTTCATCCGTTTCGAGCTGCGGGCGATCTGCAAACGGCTCGCCGAAAGGATCGCGACGCGGATCACCGCGCTTGAGAGCTAA"
                   , "ATGGCCGTGGAAAAGACCAACAGCAGCAGCTCCCTGGCCGAAGTGATCGACAGAATCCTGGACAAGGGCATCGTGATCGACGCCTGGGTGCGCGTGTCCCTCGTGGGAATTGAGCTGCTGGCCATCGAGGCCCGGATCGTGATTGCCAGCGTGGAAACCTACCTGAAGTACGCCGAGGCCGTGGGCCTGACACAGAGTGCTGCTGTGCCTGCTTGA"
                   ]
    # Used to store all the ground truth sequence which we expected the SemperRecode() to return to compare with the generated result.

    '''
    Flow:
    1. Parse fasta file line by line
    2. Create an SemperRecode instance called "obj" for each line of sequence ("seq")
    3. Assign a str of modified sequence with lower efficiency (which is returned from process_sequence function) to "new_seq"
    4. Convert "new_seq" into Seq object with id and store it in "modified_seq" (Seq object)
    5. Append "modified_seq" to the list "output" of Seq object which store all the generated modified sequence
    6. Use SeqIO.write to convert every sequence in "output" into a single fasta file
       and store the file in the designated location according to "output_file" (file path)
    '''

    with open("tests/sample_file/3_sample_inputs.fasta", 'r') as file:
        for line in SeqIO.parse(file, 'fasta'):
            input = str(line.seq)
            obj = SemperRecode(input)
            new_seq = obj.process_sequence()[0]

            # Assert the generated sequence match the correctly modified sequence
            assert new_seq in groud_truth 

            modified_seq = SeqRecord(Seq(new_seq), id=f"{line.id}_semper_recode", description='')
            output.append(modified_seq)

    output_file = "tests/sample_file/3_sample_outputs.fasta"

    with open(output_file, 'w') as file:
        SeqIO.write(output, file, 'fasta')

def test_dataframe_export():
    """
    Purpose: Test the `process_sequence` method in the SemperRecode class and ensure that the dataframe is exported successfully 
             with all the required information.
    Goal: Verify the functionality of the SemperRecode class by processing a sample sequence and the store the return data 
          from process_sequence() as a csv.file
    """
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
