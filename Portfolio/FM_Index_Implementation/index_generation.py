"""
This script creates an FM-index from a reference FASTA file using the iv2py library.
It performs the following steps:
1. Parses the command line argument to obtain the reference FASTA file path.
2. Starts memory tracking using tracemalloc.
3. Loads and normalizes the sequences from the FASTA file.
4. Creates an FM-index from the normalized sequences while measuring runtime.
5. Takes memory snapshots before and after index creation and reports the top memory differences.
6. Reports current and peak memory usage.
7. Saves the FM-index to disk.
8. Stops memory tracking.

Example usage:
    python script.py reference.fasta
"""

import iv2py as iv
import difflib
import numpy as np
import math
import time
import argparse
from memory_profiler import memory_usage
import tracemalloc

# Set up command line argument parsing.
parser = argparse.ArgumentParser(description='Create an FM-index from a reference FASTA file.')
parser.add_argument("file1", type=str, help="Reference FASTA file")
args = parser.parse_args()

# Extract the reference file path from command line arguments.
file1 = args.file1

# Start tracking memory usage with tracemalloc.
tracemalloc.start()

# Take a snapshot of current memory usage before index creation.
snapshot1 = tracemalloc.take_snapshot()

# Load the reference FASTA file and normalize sequences.
# iv.fasta.reader iterates over FASTA records, and iv.normalize cleans each sequence.
reference1 = [iv.normalize(rec.seq) for rec in iv.fasta.reader(file1)]

# Record the start time for FM-index creation.
index_start_time = time.time()

# Create the FM-index from the normalized reference sequences with a sampling rate of 16.
index = iv.fmindex(reference=reference1, samplingRate=16)

# Record the end time after FM-index creation.
index_end_time = time.time()

# Take another memory snapshot after the FM-index has been created.
snapshot2 = tracemalloc.take_snapshot()

# Compare memory snapshots to get differences, grouped by code line.
stats = snapshot2.compare_to(snapshot1, 'lineno')

# Print the top 3 memory usage differences.
print("Memory usage:")
for stat in stats[:3]:
    print(f"{stat.size_diff / 1024 / 1024:.2f} MB: {stat.count_diff} blocks: {stat}")

# Get overall current and peak memory usage from tracemalloc.
current, peak = tracemalloc.get_traced_memory()
print(f"\nCurrent memory usage: {current / 1024 / 1024:.2f} MB")
print(f"Peak memory usage: {peak / 1024 / 1024:.2f} MB")

# Calculate and print the runtime for FM-index creation.
index_runtime = index_end_time - index_start_time
print(f"\nFM Index Runtime: {index_runtime:.2f} seconds")

# Save the constructed FM-index to disk at the specified path.
index.save("./data/file.idx")

# Stop tracking memory.
tracemalloc.stop()