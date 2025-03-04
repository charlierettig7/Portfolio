"""
This script demonstrates approximate sequence matching using the iv2py library,
combined with the pigeonhole principle to allow for a specified number of mismatches.
It performs the following steps:

1. Parses command line arguments to obtain:
   - A reference FASTA file.
   - A query FASTA file.
   - The number of queries to run.
   - The allowed number of mismatches (errors).

2. Loads and normalizes the sequences from the reference file.

3. Loads a pre-built FM-index from disk (from "./data/file.idx").

4. Implements two search functions:
   - fm_search: Performs a simple FM-index search on a given query.
   - pigeonhole: Uses the pigeonhole principle to split the query into segments,
     searches each segment, and then refines candidate matches by checking the
     entire query against the reference.

5. Iterates over query sequences (up to the specified repeats) and applies the pigeonhole search,
   printing the matching positions.

6. Monitors runtime and memory usage using time, tracemalloc, and memory_profiler.

Example usage:
    python script.py reference.fasta query.fasta 100 2
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
parser = argparse.ArgumentParser(description='Approximate sequence matching using pigeonhole principle.')
parser.add_argument("file1", type=str, help="Reference FASTA file")
parser.add_argument("file2", type=str, help="Query FASTA file")
parser.add_argument("repeats", type=int, help="Number of queries to run")
parser.add_argument("errors", type=int, help="Number of mismatches allowed")
args = parser.parse_args()

# Extract command line arguments.
file1 = args.file1       # Reference file path.
file2 = args.file2       # Query file path.
repeats = args.repeats   # Maximum number of queries to process.
errors = args.errors     # Allowed mismatches per query.

# Load the reference file and normalize sequences.
# iv.fasta.reader iterates over FASTA records and iv.normalize cleans each sequence.
reference1 = [iv.normalize(rec.seq) for rec in iv.fasta.reader(file1)]

# Load a pre-built FM-index from disk (using a fixed path here).
index = iv.fmindex(path="./data/file.idx")

def find_contained_positions(s1, s2, k, spot):
    """
    Slide over string s2 to find positions where s1 matches with at most k mismatches.
    Each matching position is adjusted by the 'spot' value.
    
    Parameters:
        s1 (str): The query sequence segment.
        s2 (str): The portion of the reference sequence to search.
        k (int): Maximum allowed mismatches.
        spot (int): Offset to add to found positions.
    
    Returns:
        List of positions (integers) where s1 is found within s2 with <= k mismatches.
    """
    positions = []
    # Iterate over every possible alignment of s1 in s2.
    for i in range(len(s2) - len(s1) + 1):
        # Count mismatches between s1 and the current window in s2.
        differences = sum(c1 != c2 for c1, c2 in zip(s1, s2[i:i+len(s1)]))
        if differences <= k:
            positions.append(i + spot)
    return positions

def fm_search(index, querry):
    """
    Perform a search using the FM-index for the given query.
    
    Parameters:
        index: The pre-built FM-index.
        querry (str): The query sequence.
    
    Returns:
        List of match positions (extracted from search results).
    """
    result1 = index.search(querry)
    result2 = []
    # Extract the match positions from the search results.
    for res1 in result1:
        result2.append(res1[1])
    return result2

def pigeonhole(index, querry, errors, reference):
    """
    Use the pigeonhole principle to search for the query sequence in the reference.
    The query is split into segments based on allowed errors, and then each segment is
    searched via the FM-index. Candidate matches are refined by checking the entire query.
    
    Parameters:
        index: The pre-built FM-index.
        querry (str): The query sequence.
        errors (int): Allowed number of mismatches.
        reference (str): The reference sequence.
    
    Returns:
        List of positions in the reference where the query is found with <= errors mismatches.
    """
    used_matches = []
    # If no mismatches are allowed, perform a direct FM-index search.
    if errors == 0:
        return fm_search(index, querry)
    
    results = []
    # Calculate the segment length using ceiling division.
    seg_len = math.ceil(len(querry) / (errors + 1))
    # Split the query into segments of approximately equal length.
    segments = [querry[i: i + math.ceil(seg_len)] for i in range(0, len(querry), math.ceil(seg_len))]
    
    # For each segment, perform an FM-index search.
    for segment in segments:
        matches = fm_search(index, segment)
        current_errors = 0
        
        # Create new segments by concatenating adjacent segments (for refinement).
        new_segments = []
        for ii in range(math.floor(len(segments) / 2)):
            if ii == math.floor(len(segments) / 2):
                new_segments.append(segments[ii] + segments[ii+1] + segments[ii+2])
            else:
                new_segments.append(segments[ii] + segments[ii+1])
        
        # Iteratively refine candidate matches.
        while matches:
            # Remove duplicate matches.
            super_temp_matches = []
            for super_temp in matches:
                if super_temp in super_temp_matches:
                    continue
                else:
                    super_temp_matches.append(super_temp)
            matches = super_temp_matches
            
            temp_matches = []
            # Update current allowed error threshold.
            current_errors = current_errors * 2 + 1
            if current_errors >= errors:
                # If current errors exceed allowed errors, verify candidate matches.
                temp_results = []
                for match in matches:
                    if match in used_matches:
                        continue
                    else:
                        used_matches.append(match)
                        begin = max(0, match - len(querry) - errors)
                        end = min(len(reference), match + len(querry) + errors)
                        for ii in find_contained_positions(querry, reference[begin: end], errors, begin):
                            if ii in results or ii in temp_results:
                                continue
                            else:
                                temp_results.append(ii)
                for ii in temp_results:
                    results.append(ii)
                matches = []
            else:
                # Refine matches by checking new segments.
                for match in matches:
                    for new_segment in new_segments:
                        begin = max(0, match - len(new_segment) - errors)
                        end = min(len(reference[0]), match + len(new_segment) + errors)
                        for ii in find_contained_positions(reference[begin: end], querry, errors, begin):
                            if ii in temp_matches:
                                continue
                            else:
                                temp_matches.append(ii)
                temp_segments = []
                for ii in range(math.floor(len(segments) / 2)):
                    if ii == math.floor(len(segments) / 2):
                        temp_segments.append(segments[ii] + segments[ii+1] + segments[ii+2])
                    else:
                        temp_segments.append(segments[ii] + segments[ii+1])
                new_segments = temp_segments
                matches = temp_matches
    return results

# Set a timeout duration in seconds.
timeout_duration = 600
pex_start_time = time.time()
# Start tracing memory allocations.
tracemalloc.start()
pex_snapshot1 = tracemalloc.take_snapshot()

try:
    counter = 0
    # Process each query from the query FASTA file after normalizing the sequences.
    for file in [iv.normalize(rec.seq) for rec in iv.fasta.reader(file2)]:
        elapsed_time = time.time() - pex_start_time
        # If the processing time exceeds the timeout, raise an error.
        if elapsed_time > timeout_duration:
            raise TimeoutError("Timeout exceeded")
        # Process only the specified number of queries.
        if counter >= repeats:
            break
        else:
            counter += 1
            # Run the pigeonhole search on the query and print the results.
            result_positions = pigeonhole(index, file, errors, reference1[0])
            print("sequence:" + str(counter) + " " + str(result_positions))
    
    pex_end_time = time.time()
    pex_snapshot2 = tracemalloc.take_snapshot()
    
    # Compare memory snapshots to obtain memory differences.
    pex_stats = pex_snapshot2.compare_to(pex_snapshot1, 'lineno')
    
    # Print top 3 memory usage differences.
    print("\nDetailed memory usage:")
    for stat in pex_stats[:3]:
        print(f"{stat.size_diff / 1024 / 1024:.2f} MB: {stat.count_diff} blocks: {stat}")
    
    # Report current and peak memory usage.
    current, peak = tracemalloc.get_traced_memory()
    print(f"\nCurrent memory usage: {current / 1024 / 1024:.2f} MB")
    print(f"Peak memory usage: {peak / 1024 / 1024:.2f} MB")
    
    # Calculate the total memory difference.
    total_memory = sum(stat.size_diff for stat in pex_stats) / (1024 * 1024)  # in MB
    
    pex_runtime = pex_end_time - pex_start_time
    print(f"\nPEX Runtime: {pex_runtime:.2f} seconds")
    print(f"Total Memory Difference: {total_memory:.2f} MB")

except TimeoutError as e:
    # Handle timeout errors if the script runs for too long.
    print(str(e))

finally:
    # Stop memory tracing and indicate program completion.
    tracemalloc.stop()
    print('Program completed or terminated')