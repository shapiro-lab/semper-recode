"""Tests for `semper_recode` package."""

from semper_recode.semper_recode import SemperRecode # Import SemperRecode class
import pytest

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
    assert isinstance(recode.seq_id, list)
    assert recode.seq_id == []

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

def test_all():
    """
    Purpose: Test the `process_sequence` method in the SemperRecode class.
    Goal: Verify the functionality of the method by processing a sample file and comparing the result with the ground truth.
    """

    # Use Ishaan's sequence (gmail)
    # Create ground truth seq to compare with the result from package
    seq = "ATGGTGTCCAAGGGCGAGGAACTGTTCACCGGCGTGGTGCCCATCCTGGTGGAACTGGATGGCGACGTGAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAAGGCGACGCCACATACGGAAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCTTGGCCTACCCTCGTGACCACACTGACCTACGGCGTGCAGTGCTTCGCCAGATACCCCGACCACATGAAGCAGCACGATTTCTTCAAGAGCGCCATGCCCGAGGGCTACGTGCAGGAACGGACCATCTTCTTCAAGGACGACGGCAACTACAAGACAAGAGCCGAAGTGAAGTTCGAGGGCGACACCCTCGTGAACCGGATCGAGCTGAAGGGCATCGACTTCAAAGAGGATGGCAACATCCTGGGCCACAAGCTGGAGTACAACTACAACAGCCACAAGGTGTACATCACCGCCGACAAGCAGAAAAACGGCATCAAAGTGAACTTCAAGACCCGGCACAACATCGAGGACGGCAGCGTGCAGCTGGCCGACCACTACCAGCAGAACACCCCCATCGGAGATGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACACAAAGCGCCCTGAGCAAGGACCCCAACGAGAAGCGGGACCACATGGTGCTGCTGGAATTTGTGACCGCCGCTGGCATCACCCTGGGCATGGACGAGCTGTACAAGTGA"
    my_obj = SemperRecode(seq)
    expected_output = "ATGGTGTCCAAGGGCGAGGAACTGTTCACCGGCGTGGTGCCCATCCTGGTGGAACTGGATGGCGACGTGAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAAGGCGACGCCACATACGGAAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCTTGGCCTACCCTCGTGACCACACTGACCTACGGCGTGCAGTGCTTCGCCAGATACCCCGATCATATGAAACAGCACGATTTCTTCAAGAGTGCTATGCCTGAGGGCTACGTGCAGGAACGGACCATCTTCTTCAAGGACGACGGCAACTACAAGACAAGAGCCGAAGTGAAGTTCGAGGGCGACACCCTCGTGAACCGGATCGAGCTGAAGGGCATCGACTTCAAAGAGGATGGCAACATCCTGGGCCACAAGCTGGAGTACAACTACAACAGCCACAAGGTGTACATCACCGCCGACAAGCAGAAAAACGGCATCAAAGTGAACTTCAAGACCCGGCACAACATCGAGGACGGCAGCGTGCAGCTGGCCGACCACTACCAGCAGAACACCCCCATCGGAGATGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACACAAAGCGCCCTGAGCAAGGACCCCAACGAGAAGCGGGATCATATGGTTCTGCTGGAATTTGTGACCGCCGCTGGCATCACCCTTGGTATGGATGAGCTGTACAAGTGA"
    my_obj.process_sequence() 

      




