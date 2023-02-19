# Rename fasta

Rename fasta includes two functions, `rename_fasta()` and `batch_rename_fasta_files()`, for renaming the sequence names in one or multiple FASTA files.

## Dependencies:


* pandas
* Biopython

## Input:

`rename_fasta()` function requires three arguments:
* **file_name**: A string specifying the path to the input FASTA file to be renamed.
* **new_name**: A string specifying the new prefix to be used for renaming the sequences.
* **output_file_name** (optional): A string specifying the path to the output file. If not provided, the output file will be named 'file_name_renamed.fasta' and saved in the same directory as the input file.

`batch_rename_fasta_files()` function requires two arguments:
* **input_dir**: A string specifying the path to the directory containing the FASTA files to be renamed.
* **new_name**: A string specifying the new prefix to be used for renaming.

## Output:

`rename_fasta()` function returns a new FASTA file with renamed sequence names.
`batch_rename_fasta_files()` function renames the sequence names in all the FASTA files in the specified directory.

### Example 1:
To rename the sequence names in the FASTA file named **K-12.fasta** with the prefix *Contig*, use the following code:

`rename_fasta('example\K-12.fasta', 'Contig')`

### Example 2:

To rename the sequence names in all the FASTA files in the directory **\example** with the prefix *Contig*, use the following code:

`batch_rename_fasta_files('example', 'Contig')`
