/*
 * This program demonstrates how to perform approximate pattern matching using
 * a FM-index loaded from a file. It uses the SeqAn3 library for command line 
 * argument parsing, file I/O, and searching.
 *
 * Main Steps:
 * 1. Parse command line arguments to get the paths for the FM-index file and 
 *    the query file, along with the number of queries to run and allowed errors.
 * 2. Load the query sequences from the query file into memory.
 * 3. Deserialize and load the FM-index from a binary file.
 * 4. Duplicate query sequences if the loaded queries are fewer than required.
 * 5. Configure the search with the allowed Hamming distance errors.
 * 6. Perform the search using the FM-index and write the results (query id and 
 *    match positions) to an output file ("FM_output.txt").
 *
 * Example usage:
 *   ./fmindex_search --index path/to/index_file --query path/to/query_file --query_ct 100 --errors 1
 */

#include <sstream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

int main(int argc, char const* const* argv)
{
    // Set up the argument parser with the program name and disable update notifications.
    seqan3::argument_parser parser{"fmindex_search", argc, argv, seqan3::update_notifications::off};

    // Set parser metadata.
    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    // Define variables for command line options.
    auto index_path = std::filesystem::path{};
    parser.add_option(index_path, '\0', "index", "path to the index file");

    auto query_file = std::filesystem::path{};
    parser.add_option(query_file, '\0', "query", "path to the query file");

    auto number_of_queries = size_t{100};
    parser.add_option(number_of_queries, '\0', "query_ct", "number of query, if not enough queries, these will be duplicated");

    auto number_of_errors = uint8_t{0};
    parser.add_option(number_of_errors, '\0', "errors", "number of allowed hamming distance errors");

    // Parse the command line arguments; catch and report any parsing errors.
    try
    {
         parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "Parsing error. " << ext.what() << "\n";
        return EXIT_FAILURE;
    }

    // Load query sequences from the query file.
    auto query_stream = seqan3::sequence_file_input{query_file};

    // Store the sequences in a vector.
    std::vector<std::vector<seqan3::dna5>> queries;
    for (auto & record : query_stream)
    {
        queries.push_back(record.sequence());
    }

    // Define a type alias for the FM-index (this is a hack to deduce the type).
    using Index = decltype(seqan3::fm_index{std::vector<std::vector<seqan3::dna5>>{}});
    Index index; // Construct an empty FM-index.

    {
        // Load the FM-index from a binary file using cereal for deserialization.
        seqan3::debug_stream << "Loading 2FM-Index ... " << std::flush;
        std::ifstream is{index_path, std::ios::binary};  // Open index file in binary mode.
        cereal::BinaryInputArchive iarchive{is};         // Create a binary archive for deserialization.
        iarchive(index);                                 // Deserialize the FM-index.
        seqan3::debug_stream << "done\n";
    }

    // Duplicate query sequences until we have at least the required number.
    while (queries.size() < number_of_queries)
    {
        auto old_count = queries.size();
        queries.resize(2 * old_count);
        std::copy_n(queries.begin(), old_count, queries.begin() + old_count);
    }
    // If there are more queries than needed, trim the vector.
    queries.resize(number_of_queries);

    // Create a search configuration with the allowed total error count.
    seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{number_of_errors}};

    // Open an output file to write the search results.
    std::ofstream outFile("FM_output.txt");
            
    // Check if the file was successfully opened.
    if (!outFile.is_open())
    {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }
    
    // Write a header line to the output file.
    outFile << "File containing the positions of the discovered markers " << std::endl;

    // Perform the search using the FM-index on the query sequences.
    for (auto && result : seqan3::search(queries, index, cfg))
    {
        // Write the result: query id and the beginning position in the reference sequence.
        outFile << "Query " << result.query_id() << " found at position: " << result.reference_begin_position() << std::endl;
    }

    return 0;
}