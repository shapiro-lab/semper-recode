from Bio import SeqIO
from Bio.Seq import Seq
import os
import argparse

def count_out_of_frame_atgs(sequence):
    """Count the number of out-of-frame ATGs in a sequence."""
    count = 0
    for i in range(len(sequence) - 2):
        if i % 3 != 0 and sequence[i:i+3] == "ATG":
            count += 1
    return count

def analyze_sequence(original, modified=None, seq_id="Sequence"):
    """
    Analyze ATG codons in a sequence and optionally compare with a modified version.
    
    Args:
        original (str): Original sequence
        modified (str, optional): Modified sequence for comparison
        seq_id (str): Sequence identifier
    """
    print(f"\n=== Analysis for {seq_id} ===")
    print(f"Sequence length: {len(original)}")
    
    # Find in-frame ATGs
    in_frame_atgs = []
    for i in range(0, len(original), 3):
        if original[i:i+3] == "ATG":
            in_frame_atgs.append(i//3 + 1)
    
    print(f"In-frame ATGs: {len(in_frame_atgs)} at positions {in_frame_atgs}")
    
    # Find out-of-frame ATGs
    out_frame_atgs = []
    for i in range(len(original) - 2):
        if i % 3 != 0 and original[i:i+3] == "ATG":
            out_frame_atgs.append(i)
    
    print(f"Out-of-frame ATGs: {len(out_frame_atgs)} at positions {out_frame_atgs}")
    
    # If modified sequence is provided, compare them
    if modified:
        print("\n--- Comparison with modified sequence ---")
        print(f"Modified sequence length: {len(modified)}")
        
        # Check if protein translation is preserved
        original_protein = str(Seq(original).translate())
        modified_protein = str(Seq(modified).translate())
        
        if original_protein == modified_protein:
            print("✓ Protein sequence is preserved")
        else:
            print("✗ Protein sequence changed")
            print(f"Original protein: {original_protein[:50]}..." if len(original_protein) > 50 else original_protein)
            print(f"Modified protein: {modified_protein[:50]}..." if len(modified_protein) > 50 else modified_protein)
        
        # Find in-frame ATGs in modified sequence
        mod_in_frame_atgs = []
        for i in range(0, len(modified), 3):
            if modified[i:i+3] == "ATG":
                mod_in_frame_atgs.append(i//3 + 1)
        
        # Find out-of-frame ATGs in modified sequence
        mod_out_frame_atgs = []
        for i in range(len(modified) - 2):
            if i % 3 != 0 and modified[i:i+3] == "ATG":
                mod_out_frame_atgs.append(i)
        
        print(f"Modified in-frame ATGs: {len(mod_in_frame_atgs)} at positions {mod_in_frame_atgs}")
        print(f"Modified out-of-frame ATGs: {len(mod_out_frame_atgs)} at positions {mod_out_frame_atgs}")
        
        # Check if all out-of-frame ATGs were removed
        if len(mod_out_frame_atgs) == 0:
            print("✓ All out-of-frame ATGs were successfully removed")
        else:
            print(f"❌ {len(mod_out_frame_atgs)} out-of-frame ATGs remain in the modified sequence")
        
        # Find codon differences
        differences = []
        for i in range(0, min(len(original), len(modified)), 3):
            if i+3 <= len(original) and i+3 <= len(modified):
                orig_codon = original[i:i+3]
                mod_codon = modified[i:i+3]
                
                if orig_codon != mod_codon:
                    differences.append((i//3 + 1, orig_codon, mod_codon))
        
        if differences:
            print(f"\nFound {len(differences)} codon changes:")
            for pos, orig, mod in differences[:10]:  # Show only first 10 changes
                print(f"Position {pos}: {orig} → {mod}")
            
            if len(differences) > 10:
                print(f"... and {len(differences) - 10} more changes")
        else:
            print("\nNo codon differences found")

def analyze_fasta_file(original_file, modified_file=None):
    """
    Analyze sequences in a FASTA file and optionally compare with modified versions.
    
    Args:
        original_file (str): Path to original FASTA file
        modified_file (str, optional): Path to modified FASTA file for comparison
    """
    if not os.path.exists(original_file):
        print(f"Error: Original file '{original_file}' not found")
        return
    
    # Load original sequences
    original_records = list(SeqIO.parse(original_file, "fasta"))
    print(f"\nLoaded {len(original_records)} sequences from {original_file}")
    
    # Load modified sequences if available
    modified_records = None
    if modified_file:
        if not os.path.exists(modified_file):
            print(f"Warning: Modified file '{modified_file}' not found")
        else:
            modified_records = list(SeqIO.parse(modified_file, "fasta"))
            print(f"Loaded {len(modified_records)} sequences from {modified_file}")
    
    # Analyze each sequence
    total_original_out_of_frame = 0
    total_modified_out_of_frame = 0
    
    for i, orig_record in enumerate(original_records):
        orig_seq = str(orig_record.seq)
        mod_seq = None
        
        # Get corresponding modified sequence if available
        if modified_records and i < len(modified_records):
            mod_seq = str(modified_records[i].seq)
        
        # Count out-of-frame ATGs
        orig_out_of_frame = count_out_of_frame_atgs(orig_seq)
        total_original_out_of_frame += orig_out_of_frame
        
        if mod_seq:
            mod_out_of_frame = count_out_of_frame_atgs(mod_seq)
            total_modified_out_of_frame += mod_out_of_frame
        
        # Analyze the sequence
        analyze_sequence(orig_seq, mod_seq, orig_record.id)
    
    # Print summary
    print("\n=== Summary ===")
    print(f"Total sequences: {len(original_records)}")
    print(f"Total out-of-frame ATGs in original sequences: {total_original_out_of_frame}")
    
    if modified_records:
        print(f"Total out-of-frame ATGs in modified sequences: {total_modified_out_of_frame}")
        print(f"Total out-of-frame ATGs removed: {total_original_out_of_frame - total_modified_out_of_frame}")
        
        # Calculate success rate
        if total_original_out_of_frame > 0:
            success_rate = ((total_original_out_of_frame - total_modified_out_of_frame) / total_original_out_of_frame) * 100
            print(f"Success rate: {success_rate:.1f}%")

def main():
    """
    Main function that parses command line arguments and analyzes FASTA files.
    """
    parser = argparse.ArgumentParser(description='Analyze out-of-frame ATGs in FASTA files.')
    parser.add_argument('original', help='Path to original FASTA file')
    parser.add_argument('-m', '--modified', help='Path to modified FASTA file (optional)')
    args = parser.parse_args()
    
    analyze_fasta_file(args.original, args.modified)

if __name__ == "__main__":
    # If run directly, parse command line arguments
    main() 