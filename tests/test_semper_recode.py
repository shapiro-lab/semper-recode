# !/usr/bin/env python

"""Tests for `semper_recode` package."""
from semper_recode.semper_recode import SemperRecode 

# !/usr/bin/env python

"""Tests for `semper_recode` package."""
from semper_recode.semper_recode import SemperRecode 

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

# ================== TEST CONSTRUCTOR ==================

def test_constructor_with_default_values():
    obj = SemperRecode()
    assert obj.start_codon == ['AUG', 'ATG', 'GUG', 'UUG']
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

# =================== CONSTRUCTOR =======================

