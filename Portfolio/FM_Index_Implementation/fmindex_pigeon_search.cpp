/*
 * This program performs approximate pattern matching using a combination of
 * an FM-index (loaded from a file) and the pigeonhole principle for verification.
 *
 * Main Steps:
 * 1. Parse command line arguments to obtain paths for the FM-index, reference, and query files,
 *    along with the number of queries to execute and the number of allowed errors.
 * 2. Load the reference and query sequences from input files using SeqAn3.
 * 3. Load a pre-built FM-index from a binary file using cereal for deserialization.
 * 4. Duplicate query sequences if there are not enough to reach the required count.
 * 5. For each query:
 *    - Split the query into parts (using the pigeonhole principle) with the function sliceQuery.
 *    - For each part, search the FM-index.
 *    - Use the verify function to check if the match in the reference (at the given position) is valid.
 *    - If a valid match is found, write the result to an output file and move to the next query.
 *
 * Example usage:
 *   ./fmindex_pigeon_search --index index_file.bin --reference ref_file.fasta --query query_file.fasta --query_ct 100 --errors 1
 */

#include <sstream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/search/search.hpp>

// Function: sliceQuery
// Splits a given vector (query) into parts based on the parameter k (number of allowed errors).
// The query is divided into (k+1) parts to use the pigeonhole principle during approximate search.
std::vector<std::vector<seqan3::dna5>> sliceQuery(const std::vector<seqan3::dna5>& vec, size_t k)
{
    std::vector<std::vector<seqan3::dna5>> slices;
    size_t n = vec.size();

    // If the query is empty or k is zero, return the whole vector as a single slice.
    if (n == 0 || k == 0) {
        slices.emplace_back(vec.begin(), vec.end());
        return slices;
    }

    // Calculate the number of parts. We want k+1 parts, but not more than n parts.
    size_t parts_count = k + 1;
    if (parts_count > n) {
        parts_count = n; // Ensure at most n parts if query is very short.
    }

    // Calculate the chunk size using ceiling division.
    size_t chunkSize = (n + parts_count - 1) / parts_count;

    size_t start = 0;
    // Slice the vector into parts.
    for (size_t i = 0; i < parts_count; ++i)
    {
        size_t end = std::min(start + chunkSize, n);
        slices.emplace_back(vec.begin() + start, vec.begin() + end);
        start = end;
        if (start >= n) break; // Prevent creating empty parts.
    }

    return slices;
}

// Function: verify
// Checks if the candidate match at position 'pos' in the reference (r) is valid.
// It verifies the mismatches for each sliced part of the query based on its expected offset.
bool verify(const std::vector<seqan3::dna5>& r, const std::vector<std::vector<seqan3::dna5>>& parts, size_t i, int pos, int errors)
{
    int found_errors = 0; 
       
    // Loop over each part of the sliced query.
    for (size_t j = 0; j < parts.size(); ++j)
    {
        const auto& part = parts[j];

        // Calculate the expected starting position for the current part in the reference.
        const int start_pos = pos + (j - i) * parts[0].size();
        // Check if the calculated region is within the bounds of the reference.
        if (start_pos < 0 || start_pos + part.size() > static_cast<int>(r.size()))
        {
            return false;
        }
        // Compare each character in the part to the corresponding reference region and count mismatches.
        for (size_t k = 0; k < part.size(); ++k)
            found_errors += (r[start_pos + k] != part[k]);
    }
        
    // Return true if the total number of mismatches is within the allowed errors.
    return (found_errors <= errors);
}

int main(int argc, char const* const* argv)
{
    // Set up argument parser with program name and disable update notifications.
    seqan3::argument_parser parser{"fmindex_pigeon_search", argc, argv, seqan3::update_notifications::off};

    // Set author and version metadata.
    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    // Define command line options.
    auto index_path = std::filesystem::path{};
    parser.add_option(index_path, '\0', "index", "path to the query file");

    auto reference_file = std::filesystem::path{};
    parser.add_option(reference_file, '\0', "reference", "path to the reference file");

    auto query_file = std::filesystem::path{};
    parser.add_option(query_file, '\0', "query", "path to the query file");

    auto number_of_queries = size_t{100};
    parser.add_option(number_of_queries, '\0', "query_ct", "number of queries; if not enough queries, these will be duplicated");

    auto number_of_errors = uint8_t{0};
    parser.add_option(number_of_errors, '\0', "errors", "number of allowed Hamming distance errors");

    // Parse the command line arguments and handle any errors.
    try {
         parser.parse();
    } catch (seqan3::argument_parser_error const& ext) {
        seqan3::debug_stream << "Parsing error. " << ext.what() << "\n";
        return EXIT_FAILURE;
    }

    // Load the reference sequences from the reference file.
    auto reference_stream = seqan3::sequence_file_input{reference_file};
    // Load the query sequences from the query file.
    auto query_stream     = seqan3::sequence_file_input{query_file};

    // Store reference sequences in memory.
    std::vector<std::vector<seqan3::dna5>> reference;
    for (auto& record : reference_stream)
    {
        reference.push_back(record.sequence());
    }

    // Store query sequences in memory.
    std::vector<std::vector<seqan3::dna5>> queries;
    for (auto& record : query_stream)
    {
        queries.push_back(record.sequence());
    }

    // Load the FM-index from a binary file using cereal for deserialization.
    using Index = decltype(seqan3::fm_index{std::vector<std::vector<seqan3::dna5>>{}}); // Type deduction hack.
    Index index; // Construct the FM-index.
    {
        seqan3::debug_stream << "Loading 2FM-Index ... " << std::flush;
        std::ifstream is{index_path, std::ios::binary};  // Open the index file in binary mode.
        cereal::BinaryInputArchive iarchive{is};         // Create binary archive for deserialization.
        iarchive(index);                                 // Load the index from the file.
        seqan3::debug_stream << "done\n";
    }

    // Duplicate queries until we have at least the required number.
    while (queries.size() < number_of_queries)
    {
        auto old_count = queries.size();
        queries.resize(2 * old_count);
        std::copy_n(queries.begin(), old_count, queries.begin() + old_count);
    }
    // Trim the queries vector if it contains more than needed.
    queries.resize(number_of_queries);

    // Configure search with zero errors; actual allowed errors are handled by pigeonhole verification.
    seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{0}};

    // Open an output file to write the search results.
    std::ofstream outFile("Pidgeonhole_FM_output_Full.txt");
            
    // Verify that the output file is open.
    if (!outFile.is_open())
    {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }
    
    // Write a header line to the output file.
    outFile << "File containing the positions of the discovered markers " << std::endl;

    int query_count = 0; 
    bool dummy_check = false; 
    // Iterate over each query.
    for (auto& query : queries)
    {
        dummy_check = false; 
        // Slice the query into parts using the allowed number of errors.
        auto parts = sliceQuery(query, number_of_errors);
        // For each part, perform a search using the FM-index.
        for (size_t i = 0; i < parts.size(); ++i)
        {
            // Search for the current part in the FM-index.
            for (auto && result : seqan3::search(parts[i], index, cfg))
            {
                // Verify the candidate match using the verify function.
                if (verify(reference[0], parts, i, result.reference_begin_position(), number_of_errors))
                {
                     outFile << "Query " << query_count << " found at position: " << result.reference_begin_position() << std::endl; 
                     dummy_check = true; 
                }
            }
            // If a valid match was found for any part, skip to the next query.
            if (dummy_check) break; 
        }
        query_count++; 
    }
    return 0;
}