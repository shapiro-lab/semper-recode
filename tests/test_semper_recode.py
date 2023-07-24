"""Tests for `semper_recode` package."""

from semper_recode.semper_recode import SemperRecode # Import SemperRecode class
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

@pytest.fixture
def obj():
    return SemperRecode("ATGCTGATGCTAATGCGTACGTAGCTAA")

# ================== TEST CONSTRUCTOR ==================

def test_constructor_valid_sequence():
        user_seq = "ATGCATGC"
        recode = SemperRecode(user_seq)
        assert isinstance(recode, SemperRecode)
        assert recode.seq == user_seq

def test_constructor_empty_sequence():
    with pytest.raises(TypeError):
        SemperRecode()

def test_constructor_whitespace_sequence():
    with pytest.raises(ValueError):
        SemperRecode("    ")

def test_constructor_blank_sequence():
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

def test_find_out_of_frame(obj):
    """
    Purpose: Verify the functionality of the `find_out_of_frame` method in the SemperRecode class.
    Goal: Ensure that the method correctly identifies out-of-frame AUGs and even when there's extra nucleotide 
          in the sequence & in-frame AUG present and returns their positions in the given sequence as list.
    """
    assert obj.find_out_of_frame("AAAATG") == [] # No out-of-frame AUG
    assert obj.find_out_of_frame("ATGAATGAT") == [4]
    assert obj.find_out_of_frame("ATTGTTTATGGCCATTTGGAT") == [7]
    assert obj.find_out_of_frame("AATGTTTATGGCATGTTGGAT") == [1,7]
    assert obj.find_out_of_frame("AATGTATATGGCATGTTGGAT") == [1, 7]
    assert obj.find_out_of_frame("AATGTTTATGGCATGTTGGATATG") == [1, 7]
    assert obj.find_out_of_frame("ATTGTTTAGGGCCATTTGGAT") == [] # No AUG
    assert obj.find_out_of_frame("AATGTTTATGGCATGTTGAT") == [1, 7] # 2 extra nucleotide
    assert obj.find_out_of_frame("AATGTTTATGGCCTGTTGGATATG") == [1, 7] # In-frame and out-of-frame sequence


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

    assert obj.modify_TIS_in_frame("CACTGCATGTTA") == "CATTGTATGCTG"  # 34 -> 12
    assert obj.modify_TIS_in_frame("ATGCACTGCATGTTA") == "ATGCATTGTATGCTG"  # The first AUG is expected to remain the same 
    assert obj.modify_TIS_in_frame("CACTGCATGTTAATG") == "CATTGTATGCTGATG"  # The last AUG is expected to remain the same
    assert obj.modify_TIS_in_frame("AATGAAATGCTG") == "AATGAGATGCTG"  # 82 -> 80

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

    assert obj.modify_TIS_in_frame("ATGCATGTTATGCATATGCAC") == "ATGCATGTTATGCATATGCAT" # There's 2 AUG in the second TIS sequence
    '''
    ATG | CAT | GTT | ATG | CAT | ATG | CAC
    TIS sequences:
        1. CATGTTATGCAT -> CATGTTATGCAC
        --> seq = "ATGCATGTTATGCACATGCAC"
        2. ATGCATATGCAC -> ATGCACATGCAC -> ATGCATATGCAT
        --> seq = ATGCATGTTATGCATATGCAT
    '''


# ============= TEST GET_AA_ALTERNATIVE =============

def test_get_aa_alternative(obj):
    """
    Purpose: Test the `get_aa_alternative` method in the SemperRecode class.
    Goal: Ensure that the alternative codon and its fraction are returned correctly when given the old codon and .
    """

    aa, val = obj.get_aa_alternative("AGA")
    assert aa == "AGG" and val == 0.2134612754706659

    aa, val = obj.get_aa_alternative("AAT")
    assert aa == "AAC" and val == 0.5184461245612413

    aa, val = obj.get_aa_alternative("CCA")
    assert aa == "CCC" and val == 0.3214351373711047

    aa, val = obj.get_aa_alternative("TGG") # There's only 1 codon that produce W (Trp)
    assert aa == "TGG" and val == 1.0

    aa, val = obj.get_aa_alternative("ATG") # There's only 1 codon that produce M (Met)
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
    assert obj.modify_TIS_out_of_frame("AATGATATG") == "AACGATATG" # In-frame AUG at the end 
    assert obj.modify_TIS_out_of_frame("ATGAATGATATG") == "ATGAACGATATG" # Out-of-frame AUG
    assert obj.modify_TIS_out_of_frame("GCCCAATGCGCC") == "GCCCAATGTGCC" # Out-of-frame AUG
    # assert obj.modify_TIS_out_of_frame("CATGATGATGATGCC") == "####" # Consecutives Out-of-frame AUG
    # Test the function on seq with consecutive out-of-frame AUG

# ============= TEST RULE_BASED_MODIFICATION =============
def test_rule_based_modification(obj):
    assert obj.rule_based_modification("GCCTTGAATGATTTG") == "GCCTTGAACGATTTG" # Out-of-frame AUG [7]
    # Test edge acases for cases 2, case 2


# ============= TEST PROCESS_SEQUENCE =============

def test_single_seq():
    """
    Purpose: Test the `process_sequence` method in the SemperRecode class.
    Goal: Verify the functionality of the SemperRecode class by processing a sample sequence and comparing the result with the ground truth (manually modified sequence).
    """

    seq = "ATGGTGTCCAAGGGCGAGGAACTGTTCACCGGCGTGGTGCCCATCCTGGTGGAACTGGATGGCGACGTGAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAAGGCGACGCCACATACGGAAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCTTGGCCTACCCTCGTGACCACACTGACCTACGGCGTGCAGTGCTTCGCCAGATACCCCGACCACATGAAGCAGCACGATTTCTTCAAGAGCGCCATGCCCGAGGGCTACGTGCAGGAACGGACCATCTTCTTCAAGGACGACGGCAACTACAAGACAAGAGCCGAAGTGAAGTTCGAGGGCGACACCCTCGTGAACCGGATCGAGCTGAAGGGCATCGACTTCAAAGAGGATGGCAACATCCTGGGCCACAAGCTGGAGTACAACTACAACAGCCACAAGGTGTACATCACCGCCGACAAGCAGAAAAACGGCATCAAAGTGAACTTCAAGACCCGGCACAACATCGAGGACGGCAGCGTGCAGCTGGCCGACCACTACCAGCAGAACACCCCCATCGGAGATGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACACAAAGCGCCCTGAGCAAGGACCCCAACGAGAAGCGGGACCACATGGTGCTGCTGGAATTTGTGACCGCCGCTGGCATCACCCTGGGCATGGACGAGCTGTACAAGTGA"
    my_obj = SemperRecode(seq)
    expected = "ATGGTGTCCAAGGGCGAGGAACTGTTCACCGGCGTGGTGCCCATCCTGGTGGAACTGGACGGCGACGTGAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAAGGCGACGCCACATACGGAAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCTTGGCCTACCCTCGTGACCACACTGACCTACGGCGTGCAGTGCTTCGCCAGATACCCCGATCATATGAAACAGCACGATTTCTTCAAGAGTGCTATGCCTGAGGGCTACGTGCAGGAACGGACCATCTTCTTCAAGGACGACGGCAACTACAAGACAAGAGCCGAAGTGAAGTTCGAGGGCGACACCCTCGTGAACCGGATCGAGCTGAAGGGCATCGACTTCAAAGAGGACGGCAACATCCTGGGCCACAAGCTGGAGTACAACTACAACAGCCACAAGGTGTACATCACCGCCGACAAGCAGAAAAACGGCATCAAAGTGAACTTCAAGACCCGGCACAACATCGAGGACGGCAGCGTGCAGCTGGCCGACCACTACCAGCAGAACACCCCCATCGGAGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACACAAAGCGCCCTGAGCAAGGACCCCAACGAGAAGCGGGATCATATGGTTCTGCTGGAATTTGTGACCGCCGCTGGCATCACCCTTGGTATGGACGAGCTGTACAAGTGA"
    assert my_obj.process_sequence() == expected

def test_fasta_file():
    """
    Purpose: Test the overall implementation in the SemperRecode class (from file parsing to exporting).
    Goal: Verify the functionality of the SemperRecode class when implementing in real world scenerio.
    """
    modified_seq = "" # Used to store modified sequence with metadata (Seq object)
    output = [] # Used to store all the Seq object of modified sequence for file exporting
    groud_truth = ["ATGGTGTCCAAGGGCGAGGAACTGTTCACCGGCGTGGTGCCCATCCTGGTGGAACTGGACGGCGACGTGAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAAGGCGACGCCACATACGGAAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCTTGGCCTACCCTCGTGACCACACTGACCTACGGCGTGCAGTGCTTCGCCAGATACCCCGATCATATGAAACAGCACGATTTCTTCAAGAGTGCTATGCCTGAGGGCTACGTGCAGGAACGGACCATCTTCTTCAAGGACGACGGCAACTACAAGACAAGAGCCGAAGTGAAGTTCGAGGGCGACACCCTCGTGAACCGGATCGAGCTGAAGGGCATCGACTTCAAAGAGGACGGCAACATCCTGGGCCACAAGCTGGAGTACAACTACAACAGCCACAAGGTGTACATCACCGCCGACAAGCAGAAAAACGGCATCAAAGTGAACTTCAAGACCCGGCACAACATCGAGGACGGCAGCGTGCAGCTGGCCGACCACTACCAGCAGAACACCCCCATCGGAGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACACAAAGCGCCCTGAGCAAGGACCCCAACGAGAAGCGGGATCATATGGTTCTGCTGGAATTTGTGACCGCCGCTGGCATCACCCTTGGTATGGACGAGCTGTACAAGTGA"
                   , "ATGGCGGAAGGCTCCGTCGCCAGGCAGCCTGACCTCTTGACCTGCGAACATGAGGAGATCCACCTCGCCGGCTCGATCCAGCCGCATGGAGCGCTTCTGGTCGTCAGCGAACATGACCATCGCGTCATCCAGGCCAGCGCCAACGCCGCGGAATTTCTGAATCTCGGAAGCGTACTCGGCGTTCCGCTCGCCGAGATCGACGGCGATCTGTTGATCAAGATCCTGCCGCATCTCGATCCCACCGCCGAAGGGATGCCGGTCGCGGTGCGCTGCCGGATCGGCAATCCCTCTACGGAGTACTGCGGGTTGATGCATCGGCCTCCGGAAGGCGGGCTGATCATCGAACTCGAACGTGCCGGCCCGTCGATCGATCTGTCAGGCACGCTGGCGCCGGCGCTGGAGCGGATCCGCACGGCGGGTTCACTGCGCGCGCTGTGCGATGATACCGTGCTGCTGTTTCAGCAGTGCACCGGCTACGACCGTGTTATGGTGTATCGTTTCGACGAGCAAGGCCACGGCCTGGTATTCTCCGAGTGCCACGTGCCTGGGCTCGAATCCTATTTCGGCAACCGCTATCCGTCGTCGACTGTCCCACAGATGGCGCGGCAGCTGTACGTGCGGCAGCGCGTCCGCGTGCTGGTCGACGTCACCTATCAGCCGGTGCCGCTGGAGCCGCGGCTGTCGCCGCTGACCGGGCGCGATCTTGATATGAGTGGCTGCTTCCTGCGGTCTATGAGTCCGTGCCATCTGCAGTTCCTGAAGGATATGGGCGTGCGCGCCACCCTGGCGGTGTCGCTGGTGGTCGGCGGCAAGCTGTGGGGCCTGGTTGTCTGTCACCATTATCTGCCGCGCTTCATCCGTTTCGAGCTGCGGGCGATCTGCAAACGGCTCGCCGAAAGGATCGCGACGCGGATCACCGCGCTTGAGAGCTAA"
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
            new_seq = obj.process_sequence()

            # Assert the generated sequence match the correctly modified sequence
            # assert new_seq in groud_truth 

            modified_seq = SeqRecord(Seq(new_seq), id=f"{line.id}_semper_recode", description='')
            output.append(modified_seq)

    output_file = "tests/sample_file/3_sample_outputs.fasta"

    with open(output_file, 'w') as file:
        SeqIO.write(output, file, 'fasta')




