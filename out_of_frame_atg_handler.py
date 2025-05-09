from semper_recode.semper_recode import SemperRecode
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import argparse

def process_out_of_frame_atgs(sequence):
    """
    Process a nucleotide sequence to remove out-of-frame ATGs.
    
    Args:
        sequence (str): Nucleotide sequence starting with ATG and ending with a stop codon
        
    Returns:
        str: Modified sequence without out-of-frame ATGs
        list: Any errors encountered
    """
    # Create SemperRecode object with the input sequence
    recoder = SemperRecode(sequence)
    
    # Process the sequence - this will handle both in-frame and out-of-frame ATGs
    modified_seq, errors = recoder.process_sequence()
    
    return modified_seq, errors

def save_to_fasta(sequences, output_file="modified_sequence.fasta"):
    """
    Save sequences to a FASTA file.
    
    Args:
        sequences (list): List of SeqRecord objects to save
        output_file (str): Output file path
    """
    # Write to FASTA file
    with open(output_file, "w") as handle:
        SeqIO.write(sequences, handle, "fasta")

def process_fasta_file(input_file, output_file=None):
    """
    Process all sequences in a FASTA file to remove out-of-frame ATGs.
    
    Args:
        input_file (str): Path to input FASTA file
        output_file (str, optional): Path to output FASTA file. 
                                    If None, will use input_file_modified.fasta
    
    Returns:
        list: List of processed SeqRecord objects
        dict: Statistics about the processing
    """
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file {input_file} not found")
    
    if output_file is None:
        # Generate output filename based on input filename
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}_modified.fasta"
    
    # Read sequences from input file
    records = list(SeqIO.parse(input_file, "fasta"))
    processed_records = []
    stats = {
        "total_sequences": len(records),
        "processed_sequences": 0,
        "out_of_frame_atgs_removed": 0,
        "errors": []
    }
    
    # Process each sequence
    for record in records:
        try:
            # Get the sequence string
            seq_str = str(record.seq)
            
            # Process the sequence to remove out-of-frame ATGs
            modified_seq, errors = process_out_of_frame_atgs(seq_str)
            
            # Count out-of-frame ATGs in original sequence
            out_of_frame_atgs = 0
            for i in range(len(seq_str) - 2):
                if i % 3 != 0 and seq_str[i:i+3] == "ATG":
                    out_of_frame_atgs += 1
            
            stats["out_of_frame_atgs_removed"] += out_of_frame_atgs
            
            # Create a new SeqRecord with the modified sequence
            new_record = SeqRecord(
                Seq(modified_seq),
                id=record.id,
                description=f"{record.description} - Out-of-frame ATGs removed"
            )
            
            # Add the new record to the list
            processed_records.append(new_record)
            stats["processed_sequences"] += 1
            
            # Add any errors to the stats
            if errors:
                stats["errors"].append({
                    "sequence_id": record.id,
                    "errors": errors
                })
        
        except Exception as e:
            stats["errors"].append({
                "sequence_id": record.id,
                "errors": [str(e)]
            })
    
    # Save processed sequences to output file
    save_to_fasta(processed_records, output_file)
    
    # Add output file path to stats
    stats["output_file"] = output_file
    
    # Print processing summary
    print(f"✓ Processed {stats['processed_sequences']} sequences")
    print(f"✓ Removed {stats['out_of_frame_atgs_removed']} out-of-frame ATGs")
    print(f"✓ Results saved to {stats['output_file']}")
    
    if stats["errors"]:
        print(f"Encountered {len(stats['errors'])} errors during processing")
    
    return processed_records, stats

def main():
    """
    Main function that asks for input file path and processes the file.
    """
    parser = argparse.ArgumentParser(description='Process a FASTA file to remove out-of-frame ATGs.')
    parser.add_argument('-i', '--input', help='Input FASTA file path')
    parser.add_argument('-o', '--output', help='Output FASTA file path (optional)')
    args = parser.parse_args()
    
    # If input file path is not provided as command line argument, ask for it
    input_file = args.input
    if not input_file:
        input_file = input("Enter the path to your FASTA file: ")
    
    # Process the file
    try:
        process_fasta_file(input_file, args.output)
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main() 