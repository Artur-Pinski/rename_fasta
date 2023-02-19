import glob
import os

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def get_id_new_id(file_name: str, new_name: str):
    """
    Read the FASTA format file, extract the sequence names,
    add new column "new" with the names "new_name_n" and
    return a dictionary of old and new names.
    
    Args:
        file_name (str): The name of the input file.
        new_name (str): The prefix to be used for the new names.

    Returns:
        dict: A dictionary containing the old and new sequence names.
    """
    # Read the fasta file and extract the sequence names
    sequence_ids_df = pd.DataFrame([record.id for record in SeqIO.parse(file_name, 'fasta')], columns=["old_id"])

    # Check for duplicate names
    if len(sequence_ids_df['old_id']) != len(sequence_ids_df['old_id'].unique()):
        raise ValueError("Duplicate sequence names found in input file.")

    # Add new column "new" with the names "new_name_n"
    sequence_ids_df["new"] = [f"{new_name}_{i}" for i in range(1, sequence_ids_df.shape[0] + 1)]    
   
    old_new_names = sequence_ids_df.set_index('old_id')['new'].to_dict()
    
    return old_new_names


def rename_fasta(file_name: str, new_name: str, output_file_name: str = None):
    """
    Rename the sequence names in a given FASTA file to 'new_name_n' format.

    Args:
        file_name (str): Path to the input FASTA file to be renamed.
        new_name (str): The new prefix to be used for renaming the sequences.
        output_file_name (str): Optional. Path to the output file. If not provided, the output file will be named
                               '{file_name}_renamed.fasta' and saved in the same directory as the input file.

    Returns:
        Renamed fasta file.

    Raises:
        ValueError: If duplicate sequence names are found in the input file.
    """
    if output_file_name is None:
        base = os.path.basename(file_name)
        output_file_name = f"{os.path.splitext(base)[0]}_renamed.fasta"

    new_seq_record = []
    old_new_names = get_id_new_id(file_name, new_name)

    with open(output_file_name, 'w') as output_file:
        for record in SeqIO.parse(file_name, 'fasta'):
            # Replace the old id with the new one
            new_id = old_new_names[record.id]
            new_seq_record.append(SeqRecord(record.seq, id=new_id, description=''))
        SeqIO.write(new_seq_record, output_file, 'fasta')

def batch_rename_fasta_files(input_dir: str, new_name: str) -> None:
    """
    Batch rename the sequence names in all FASTA files in a given directory.

    Args:
        input_dir (str): Path to the input directory containing the FASTA files.
        new_name (str): The new prefix to be used for renaming.

    Returns:
        Renamed fasta files.
    """
    # Get a list of all the FASTA files in the input directory
    fasta_files = glob.glob(os.path.join(input_dir, '*.fasta'))
    
    # Check if any FASTA files were found
    if not fasta_files:
        print("No FASTA files found in the directory.")
        return

    # Rename the sequences in each file
    for file_name in fasta_files:
        fasta_file_name = os.path.splitext(os.path.basename(file_name))[0]
        output_file_name = os.path.join(input_dir, f"{fasta_file_name}_renamed.fasta")
        rename_fasta(file_name, new_name, output_file_name)
        print(f"{fasta_file_name}.fasta have been processed.") 
    print("\nAll files have been processed.")    
