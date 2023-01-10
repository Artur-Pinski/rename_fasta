import os
import pandas as pd
from Bio import SeqIO
import glob
from Bio.SeqRecord import SeqRecord

def get_id_new_id(file_name: str, new_name: str):
    """
    Read the FASTA format file, extract the sequence names,
    add new column "new" with the names "new_name_n" and 
    return a dictionary of old and new names
    """
    # Read the fasta file and extract the sequence names
    data = pd.DataFrame([record.id for record in SeqIO.parse(file_name, 'fasta')], columns=["old_id"])
        
    # Add new column "new" with the names "new_name_n"
    data["new"] = [f"{new_name}_{i}" for i in range(1, data.shape[0] + 1)]
                        
    data = data.set_index('old_id').T.to_dict('index')['new']     
    return data

def rename_fasta(file_name: str, new_name: str):
    """
    Rename the fasta file's sequence names to 'new_name_n' format
    """
    new_seqs = []
    old_new_names = get_id_new_id(file_name, new_name)
    base, _ = os.path.splitext(file_name)
    for record in SeqIO.parse(file_name, 'fasta'):
        # Replace the old id with the new one
        new_id = old_new_names[record.id]
        new_seqs.append(SeqRecord(record.seq, id=new_id, description=record.description))
    # write the new sequences in new fasta file 
    SeqIO.write(new_seqs, f"{base}_renamed.fasta", 'fasta')
    
for files in glob.glob('*.fasta'):
    rename_fasta(files, 'Contig')
