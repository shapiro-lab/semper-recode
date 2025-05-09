from semper_recode.semper_recode import SemperRecode
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def handle_out_of_frame_atg(sequence):
    """
    Process a sequence to remove out-of-frame ATGs while preserving protein sequence.
    
    Args:
        sequence (str): Input nucleotide sequence starting with ATG and ending with stop codon
        
    Returns:
        str: Modified sequence with out-of-frame ATGs removed
    """
    # Create SemperRecode object
    recoder = SemperRecode(sequence)
    
    # Process sequence to remove out-of-frame ATGs
    modified_seq, errors = recoder.process_sequence()
    
    return modified_seq

def save_to_fasta(sequence, output_file="output.fasta"):
    """
    Save sequence to FASTA file.
    
    Args:
        sequence (str): Sequence to save
        output_file (str): Output FASTA file path
    """
    # Create SeqRecord object
    record = SeqRecord(Seq(sequence), id="modified_sequence", description="Sequence with out-of-frame ATGs removed")
    
    # Write to FASTA file
    with open(output_file, 'w') as handle:
        SeqIO.write(record, handle, 'fasta')

# Test with provided sequence
if __name__ == "__main__":
    test_seq = "ATGGTTAGTAAAGGAGAAGAACTCTTCACCGGCGTAGTTCCGATTCTCGTTGAGCTGGATGGAGATGTTAACGGGCACAAATTCTCTGTAAGGGGTGAAGGGGAAGGCGATGCTACCAATGGAAAACTGACCTTAAAATTTATTTGTACAACAGGTAAGCTCCCAGTTCCCTGGCCCACCCTTGTTACCACTCTGACTCACGGCGTTCAGTGTTTCTCTCGTTATCCAGATCATATAAAGCGGCATGACTTCTTTAAGTCCGCTCTCCCCGAAGGCTACGTCCAAGAAAGAACCATTAGTTTTAAGGACGATGGAACTTACAAAACCAGAGCAGAGGTCAAGTTCGAGGGAGATACTCTCGTCAATCGGATCGAACTGAAAGGAATAGACTTTAAAGAGGACGGCAATATCCTGGGTCACAAGCTGGAGTATAATTTCAATAGCCATAATGTGTATATCACGGCCGATAAGCAGAAGAACGGTATTAAGGCAAATTTTAAGATTCGCCATAATGTAGAAGATGGTTCCGTCCAATTGGCGGACCATTATCAACAAAACACTCCAATTGGAGATGGCCCGGTATTACTGCCTGATAATCACTACCTGTCCACCCAATCAAAACTCTCCAAAGATCCCAACGAGAAACGGGACCACGCGGTTCTGTTAGAATTTGTAACGGCCGCTGGTATTACGCATGGTAAGGATGAACTGTATAAGTGA"
    
    # Process sequence
    modified_seq = handle_out_of_frame_atg(test_seq)
    
    # Save to FASTA
    save_to_fasta(modified_seq, "sfbfp_modified.fasta") 