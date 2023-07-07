"""Tests for `semper_recode` package."""

from semper_recode.semper_recode import SemperRecode # Import SemperRecode class
import pytest

@pytest.fixture
def obj():
    return SemperRecode()

# ================== TEST FIND_IN_FRAME ==================

def test_find_in_frame():
    # In-frame AUGs
    sample1 = "ATGCTGATGCTAATGCGTACGTAGCTAA"
    obj1 = SemperRecode(sample1)
    assert obj1.find_in_frame(sample1) == [0,6,12] 
    
    # Out-of-frame AUGs
    sample2 = "AUTGTTTAUGGCCATTUGSAUG"
    obj2 = SemperRecode(sample2)
    assert obj2.find_in_frame(sample2) == []

    sample3 = "AAAATG"
    obj3 = SemperRecode(sample3)
    assert obj3.find_in_frame(sample3) == [3]

# ================== TEST FIND_OUT_OF_FRAME ==================

def test_find_out_of_frame():
    sample1 = "AUGCCCAAATGUCCCUUU"
    obj1 = SemperRecode(sample1)
    assert obj1.find_out_of_frame(sample1) == [8]

    sample1 = "AUGAAATGATGUCCCNNN"
    obj1 = SemperRecode(sample1)
    assert obj1.find_out_of_frame(sample1) == [5,8]


# ================== TEST CONSTRUCTOR ==================

def test_constructor_with_default_values():
    obj = SemperRecode()
    assert obj.start_codon == ['AUG', 'ATG']
    assert obj.master_df is not None
    assert obj.seq == 'ATGCTGACGGTAUGGACTTACCTGTATGCGTGCTAAATGCTAAGGCTGGTGCCGACCGGACCGTTGGGAGCGCTGTTGACCGGATGCTAAAGGGCCCGAGTCTTGTAGTACCGGACTTAAATGCGTTGTTTGACACCTGTT'

def test_constructor_with_user_seq():
    user_seq = 'ATCGTA'
    obj = SemperRecode(user_seq=user_seq)
    assert obj.master_df is not None
    assert obj.seq == user_seq

def test_constructor_with_input_file_path():
    input_file_path = 'path/to/input.fasta'
    obj = SemperRecode(input_file_path=input_file_path)
    assert obj.input_file == input_file_path

def test_constructor_with_all_arguments():
    user_seq = 'ATCGTA'
    current_path = 'data/'
    input_file_path = 'path/to/input.fasta'
    obj = SemperRecode(user_seq=user_seq, current_path=current_path, input_file_path=input_file_path)
    assert obj.seq == user_seq
    assert obj.input_file == input_file_path

def test_load_data_with_existing_files():
    current_path = 'data/'
    obj = SemperRecode(current_path=current_path)
    assert obj.master_df is not None

# ============= TEST EFFICIENCY_LEVEL =============

def test_efficiency_level_with_exist_sequence(obj):
    assert(obj.efficiency_level("CATTGTATGCT") == 12)
    assert(obj.efficiency_level("AAACGTATGGG") == 62)
    assert(obj.efficiency_level("ATTAGAATGAC") == 95)
    assert(obj.efficiency_level("CCTATGATGTA") == 102)
    assert(obj.efficiency_level("TTCATCATGCA") == 150)
    
def test_efficiency_level_with_non_exist_sequence(obj):

    # ValueError is expected to raise as the sequence doesn't exist in master_df
    with pytest.raises(ValueError):
        obj.efficiency_level("CCCTGTACGCT")
        obj.efficiency_level("CCCTTTAAACT")

# ============= TEST MODIFY_TIS_IN_FRAME =============

def test_modify_TIS_with_exist_sequence(obj):
    assert obj.modify_TIS_in_frame("CACTGCATGTTA") == "CATTGTATGCTG"  # 34 -> 12
    assert obj.modify_TIS_in_frame("AUGCACTGCATGTTA") == "AUGCATTGTATGCTG"  # The first AUG is expected to remain the same 
    assert obj.modify_TIS_in_frame("CACTGCATGTTAATG") == "CATTGTATGCTGATG"  # The last AUG is expected to remain the same
    assert obj.modify_TIS_in_frame("AATGAAATGCTG") == "AATGAGATGCTG"  # 82 -> 80


def test_modify_TIS_with_non_exist_sequence(obj):
    assert obj.modify_TIS_in_frame("ATGCATGTTA") == "ATGCATGTTA"  # No AUG in the sequence, so no modification expected
    assert obj.modify_TIS_in_frame("AUGCACTGCATG") == "AUGCACTGCATG"  # No AUG in the internal region, so no modification expected   

# ============= TEST PROCESS_SEQUENCE =============

def test_all():
    # Use Ishaan's sequence (gmail)
    # Create ground truth seq to compare with the result from package
    my_obj = SemperRecode(input_file_path = "tests/sample_file/sample_file.fasta")
    my_obj.process_sequence()


