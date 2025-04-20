from Bio import SeqIO
import os
import argparse

def modify_headers(input_fasta, output_dir=None):
    # Read in the FASTA file
    records = list(SeqIO.parse(input_fasta, "fasta"))

    # Get the base name of the file (excluding the directory and extension)
    file_name = os.path.splitext(os.path.basename(input_fasta))[0]

    # Modify the headers with the file name
    for record in records:
        record.id = file_name
        record.description = ""

    # Specify the output directory
    if output_dir:
        output_path = os.path.join(output_dir, file_name + "_newheaders.fasta")
    else:
        output_path = file_name + "_newheaders.fasta"

    # Write the modified records to the output file
    SeqIO.write(records, output_path, "fasta")

    print(f"Headers in {input_fasta} have been replaced with '{output_path}' and saved to {output_path}")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Modify headers in a FASTA file.')
    parser.add_argument('fasta_file', metavar='FASTA_FILE', type=str, help='Path to the input FASTA file')
    parser.add_argument('--output-dir', '-o', type=str, help='Specify the output directory')

    args = parser.parse_args()

    # Call the function with the provided arguments
    modify_headers(args.fasta_file, args.output_dir)