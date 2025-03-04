"""
This script demonstrates sequence searching using the iv2py library.
It performs the following steps:
1. Reads and normalizes sequences from a reference FASTA file (e.g., a partial hg38 genome).
2. Constructs an FM-index from the reference with a specified sampling rate.
3. Saves the FM-index to disk.
4. Performs an initial search using a hardcoded query sequence ("CGTACGTA").
5. Iterates over query sequences from a query FASTA file, rebuilds the FM-index for each query (for demonstration purposes),
   searches for the query in the reference, and prints the found position.
6. Stops processing after 11 queries.

Example usage:
    Ensure that 'hg38_partial.fasta.gz' and 'illumina_reads_100.fasta.gz' are available in the working directory,
    then run:
    
        python script.py
"""

import iv2py as iv
import numpy as np
from pathlib import Path

# Initialize a dummy variable (not used later).
x = 0

# Define the file paths for the reference genome and query sequences.
# You can uncomment the input lines below to allow user input.
# file1 = input("Enter the path (absolute works best) to your reference file: ")
file1 = 'hg38_partial.fasta.gz'
# file2 = input("Enter the path (absolute works best) to your queries file: ")
file2 = 'illumina_reads_100.fasta.gz'

# Read and normalize the reference sequences from the FASTA file.
# iv.fasta.reader iterates over records in the file, and iv.normalize processes each sequence.
reference = [iv.normalize(rec.seq) for rec in iv.fasta.reader(file1)]

# Create an FM-index from the reference sequences with a sampling rate of 16.
index = iv.fmindex(reference=reference, samplingRate=16)

# Create a Path object for the FM-index file to be saved.
path = iv.Path("file.fmindex")
print(type(path))  # Debug: print the type of the Path object.

# Save the constructed FM-index to disk.
index.save(str(iv.Path("file.fmindex")))

# Perform an initial search on the FM-index using a hardcoded query sequence.
result = index.search("CGTACGTA")

# Initialize a counter for tracking the number of queries processed.
counter = 0

# Iterate over query sequences from the query FASTA file.
for query in iv.fasta.reader(file=file2):
    # For each query, re-read and normalize the reference sequences.
    # Then rebuild the FM-index (note: rebuilding for each query is inefficient and is for demonstration only).
    reference = [iv.normalize(rec.seq) for rec in iv.fasta.reader(file1)]
    index = iv.fmindex(reference=reference, samplingRate=16)
    
    # Search the FM-index for the current query sequence.
    result = index.search(query.seq)

    # Print the query counter and the position where the query was found.
    print("The ", counter, "th query was found at position ", result)

    # Stop processing after more than 10 queries have been processed.
    if counter > 10: 
        print(counter, "th query reached!")
        break 
    
    # Increment the query counter.
    counter = counter + 1