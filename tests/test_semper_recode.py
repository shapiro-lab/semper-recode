"""Tests for `semper_recode` package."""

from semper_recode.semper_recode import SemperRecode # Import SemperRecode class
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

@pytest.fixture
def obj():
    return SemperRecode("ATGCTGATGCTAATGCGTACGTAGCTAA")

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
    sample2 = "ATTGTTTATGGCCATTTGSATG"
    obj2 = SemperRecode(sample2)
    assert obj2.find_in_frame(sample2) == []

    sample3 = "AAAATG"
    obj3 = SemperRecode(sample3)
    assert obj3.find_in_frame(sample3) == [3]

# ================== TEST CONSTRUCTOR ==================

def test_constructor_valid_sequence():
        user_seq = "ATGCATGC"
        recode = SemperRecode(user_seq)
        assert isinstance(recode, SemperRecode)
        assert recode.seq == user_seq

def test_constructor_empty_sequence():
    with pytest.raises(ValueError):
        SemperRecode("")

def test_constructor_whitespace_sequence():
    with pytest.raises(ValueError):
        SemperRecode("    ")

def test_constructor_blank_sequence():
    with pytest.raises(ValueError):
        SemperRecode("")

def test_constructor_seq_id_initialization():
    user_seq = "ATGCATGC"
    recode = SemperRecode(user_seq)

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

    assert obj.modify_TIS_in_frame("ATGCATGTTA") == "ATGCATGTTA"  # No AUG in the sequence, so no modification expected
    assert obj.modify_TIS_in_frame("ATGCACTGCATG") == "ATGCACTGCATG"  # No AUG in the internal region, so no modification expected   

# ============= TEST GET_AA_ALTERNATIVE =============

def test_get_aa_alternative(obj):
    """
    Purpose: Test the `get_aa_alternative` method in the SemperRecode class.
    Goal: Ensure that the alternative codon and its fraction are returned correctly when given the old codon and .
    """

    assert obj.get_aa_alternative("GCC") == "A"

# ============= TEST PROCESS_SEQUENCE =============

def test_single_seq():
    """
    Purpose: Test the `process_sequence` method in the SemperRecode class.
    Goal: Verify the functionality of the SemperRecode class by processing a sample sequence and comparing the result with the ground truth (manually modified sequence).
    """

    seq = "ATGGTGTCCAAGGGCGAGGAACTGTTCACCGGCGTGGTGCCCATCCTGGTGGAACTGGATGGCGACGTGAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAAGGCGACGCCACATACGGAAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCTTGGCCTACCCTCGTGACCACACTGACCTACGGCGTGCAGTGCTTCGCCAGATACCCCGACCACATGAAGCAGCACGATTTCTTCAAGAGCGCCATGCCCGAGGGCTACGTGCAGGAACGGACCATCTTCTTCAAGGACGACGGCAACTACAAGACAAGAGCCGAAGTGAAGTTCGAGGGCGACACCCTCGTGAACCGGATCGAGCTGAAGGGCATCGACTTCAAAGAGGATGGCAACATCCTGGGCCACAAGCTGGAGTACAACTACAACAGCCACAAGGTGTACATCACCGCCGACAAGCAGAAAAACGGCATCAAAGTGAACTTCAAGACCCGGCACAACATCGAGGACGGCAGCGTGCAGCTGGCCGACCACTACCAGCAGAACACCCCCATCGGAGATGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACACAAAGCGCCCTGAGCAAGGACCCCAACGAGAAGCGGGACCACATGGTGCTGCTGGAATTTGTGACCGCCGCTGGCATCACCCTGGGCATGGACGAGCTGTACAAGTGA"
    my_obj = SemperRecode(seq)
    expected_output = "ATGGTGTCCAAGGGCGAGGAACTGTTCACCGGCGTGGTGCCCATCCTGGTGGAACTGGATGGCGACGTGAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAAGGCGACGCCACATACGGAAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCTTGGCCTACCCTCGTGACCACACTGACCTACGGCGTGCAGTGCTTCGCCAGATACCCCGATCATATGAAACAGCACGATTTCTTCAAGAGTGCTATGCCTGAGGGCTACGTGCAGGAACGGACCATCTTCTTCAAGGACGACGGCAACTACAAGACAAGAGCCGAAGTGAAGTTCGAGGGCGACACCCTCGTGAACCGGATCGAGCTGAAGGGCATCGACTTCAAAGAGGATGGCAACATCCTGGGCCACAAGCTGGAGTACAACTACAACAGCCACAAGGTGTACATCACCGCCGACAAGCAGAAAAACGGCATCAAAGTGAACTTCAAGACCCGGCACAACATCGAGGACGGCAGCGTGCAGCTGGCCGACCACTACCAGCAGAACACCCCCATCGGAGATGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACACAAAGCGCCCTGAGCAAGGACCCCAACGAGAAGCGGGATCATATGGTTCTGCTGGAATTTGTGACCGCCGCTGGCATCACCCTTGGTATGGATGAGCTGTACAAGTGA"
    assert my_obj.process_sequence() == expected_output

def test_fasta_file():
    """
    Purpose: Test the overall implementation in the SemperRecode class (from file parsing to exporting).
    Goal: Verify the functionality of the SemperRecode class when implementing in real world scenerio.
    """
    modified_seq = "" # Used to store modified sequence with metadata (Seq object)
    output = [] # Used to store all the Seq object of modified sequence for file exporting
    groud_truth = ["ATGGTGTCCAAGGGCGAGGAACTGTTCACCGGCGTGGTGCCCATCCTGGTGGAACTGGATGGCGACGTGAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAAGGCGACGCCACATACGGAAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCTTGGCCTACCCTCGTGACCACACTGACCTACGGCGTGCAGTGCTTCGCCAGATACCCCGATCATATGAAACAGCACGATTTCTTCAAGAGTGCTATGCCTGAGGGCTACGTGCAGGAACGGACCATCTTCTTCAAGGACGACGGCAACTACAAGACAAGAGCCGAAGTGAAGTTCGAGGGCGACACCCTCGTGAACCGGATCGAGCTGAAGGGCATCGACTTCAAAGAGGATGGCAACATCCTGGGCCACAAGCTGGAGTACAACTACAACAGCCACAAGGTGTACATCACCGCCGACAAGCAGAAAAACGGCATCAAAGTGAACTTCAAGACCCGGCACAACATCGAGGACGGCAGCGTGCAGCTGGCCGACCACTACCAGCAGAACACCCCCATCGGAGATGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACACAAAGCGCCCTGAGCAAGGACCCCAACGAGAAGCGGGATCATATGGTTCTGCTGGAATTTGTGACCGCCGCTGGCATCACCCTTGGTATGGATGAGCTGTACAAGTGA", "ATGGCGGAAGGCTCCGTCGCCAGGCAGCCTGACCTCTTGACCTGCGAACATGAAGAGATC"
    "CACCTCGCCGGCTCGATCCAGCCGCATGGCGCGCTTCTGGTCGTCAGCGAACATGATCAT"
    "CGCGTCATCCAGGCCAGCGCCAACGCCGCGGAATTTCTGAATCTCGGAAGCGTACTCGGC"
    "GTTCCGCTCGCCGAGATCGACGGCGATCTGTTGATCAAGATCCTGCCGCATCTCGATCCC"
    "ACCGCCGAAGGGATGCCGGTCGCGGTGCGCTGCCGGATCGGCAATCCCTCTACGGAGTAC"
    "TGCGGGTTGATGCATCGGCCTCCGGAAGGCGGGCTGATCATCGAACTCGAACGTGCCGGC"
    "CCGTCGATCGATCTGTCAGGCACGCTGGCGCCGGCGCTGGAGCGGATCCGCACGGCGGGT"
    "TCACTGCGCGCGCTGTGCGATGACACCGTGCTGCTGTTTCAGCAGTGCACCGGCTACGAC"
    "CGTGTTATGGTGTATCGTTTCGATGAGCAAGGCCACGGCCTGGTATTCTCCGAGTGCCAT"
    "GTGCCTGGGCTCGAATCCTATTTCGGCAACCGCTATCCGTCGTCGACTGTCCCACAGATG"
    "GCGCGGCAGCTGTACGTGCGGCAGCGCGTCCGCGTGCTGGTCGACGTCACCTATCAGCCG"
    "GTGCCGCTGGAGCCGCGGCTGTCGCCGCTGACCGGGCGCGATCTTGATATGAGTGGCTGC"
    "TTCCTGCGGTCTATGAGTCCGTGCCATCTGCAGTTCCTGAAGGATATGGGCGTGCGCGCC"
    "ACCCTGGCGGTGTCGCTGGTGGTCGGCGGCAAGCTGTGGGGCCTGGTTGTCTGTCACCAT"
    "TATCTGCCGCGCTTCATCCGTTTCGAGCTGCGGGCGATCTGCAAACGGCTCGCCGAAAGG"
    "ATCGCGACGCGGATCACCGCGCTTGAGAGCTAA", 
    "ATGGCCGTGGAAAAGACCAACAGCAGCAGCTCCCTGGCCGAAGTGATCGACAGAATCCTG"
    "GACAAGGGCATCGTGATCGACGCCTGGGTGCGCGTGTCCCTCGTGGGAATTGAGCTGCTG"
    "GCCATCGAGGCCCGGATCGTGATTGCCAGCGTGGAAACCTACCTGAAGTACGCCGAGGCC"
    "GTGGGCCTGACACAGAGTGCTGCTGTGCCTGCTTGA"]
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
            assert new_seq in groud_truth 

            modified_seq = SeqRecord(Seq(new_seq), id=f"{line.id}_semper_recode", description='')
            output.append(modified_seq)

    output_file = "tests/sample_file/3_sample_outputs.fasta"

    with open(output_file, 'w') as file:
        SeqIO.write(output, file, 'fasta')




