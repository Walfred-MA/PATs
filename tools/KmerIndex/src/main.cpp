#include "kmer_index.hpp"

#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>

namespace {

void print_usage(std::ostream& output, const char* program) {
    output
        << "Usage: " << program << " --fasta assembly.fa --output index.bin [options]\n"
        << "\n"
        << "Required:\n"
        << "  -f, --fasta PATH       Uncompressed FASTA file\n"
        << "  -o, --output PATH      Binary index to create\n"
        << "\n"
        << "Options:\n"
        << "      --fai PATH         FASTA index (default: FASTA.fai)\n"
        << "  -t, --threads N        FASTA worker threads (default: CPU count)\n"
        << "      --maxcover N       Omit k-mers occurring more than N times\n"
        << "                           (default: 255, maximum: 65535)\n"
        << "      --maxspan N        Omit k-mers whose merged +/-10 kb intervals\n"
        << "                           exceed N bases (default: 1000000)\n"
        << "  -h, --help             Show this help\n";
}

std::uint64_t parse_number(const std::string& value, const std::string& option) {
    if (value.empty() || value[0] == '-') {
        throw std::runtime_error("invalid value for " + option + ": " + value);
    }
    std::size_t used = 0;
    unsigned long long parsed = 0;
    try {
        parsed = std::stoull(value, &used, 10);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid value for " + option + ": " + value);
    }
    if (used != value.size()) {
        throw std::runtime_error("invalid value for " + option + ": " + value);
    }
    return static_cast<std::uint64_t>(parsed);
}

std::string take_value(int argc, char* argv[], int& index,
                       const std::string& option) {
    if (++index >= argc) {
        throw std::runtime_error("missing value after " + option);
    }
    return argv[index];
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        kmer_index::Options options;
        const unsigned hardware_threads = std::thread::hardware_concurrency();
        options.threads = hardware_threads == 0 ? 1 : hardware_threads;

        for (int i = 1; i < argc; ++i) {
            const std::string argument = argv[i];
            if (argument == "-h" || argument == "--help") {
                print_usage(std::cout, argv[0]);
                return EXIT_SUCCESS;
            }
            if (argument == "-f" || argument == "--fasta") {
                options.fasta_path = take_value(argc, argv, i, argument);
            } else if (argument == "--fai") {
                options.fai_path = take_value(argc, argv, i, argument);
            } else if (argument == "-o" || argument == "--output") {
                options.output_path = take_value(argc, argv, i, argument);
            } else if (argument == "-t" || argument == "--threads") {
                const std::uint64_t value = parse_number(
                    take_value(argc, argv, i, argument), argument);
                if (value == 0 || value > std::numeric_limits<std::uint32_t>::max()) {
                    throw std::runtime_error("--threads must be between 1 and " +
                                             std::to_string(
                                                 std::numeric_limits<std::uint32_t>::max()));
                }
                options.threads = static_cast<std::uint32_t>(value);
            } else if (argument == "--maxcover") {
                const std::uint64_t value = parse_number(
                    take_value(argc, argv, i, argument), argument);
                if (value == 0 ||
                    value > std::numeric_limits<std::uint16_t>::max()) {
                    throw std::runtime_error(
                        "--maxcover must be between 1 and 65535");
                }
                options.max_cover = static_cast<std::uint16_t>(value);
            } else if (argument == "--maxspan") {
                options.max_span = parse_number(
                    take_value(argc, argv, i, argument), argument);
            } else {
                throw std::runtime_error("unknown option: " + argument);
            }
        }

        if (options.fasta_path.empty() || options.output_path.empty()) {
            print_usage(std::cerr, argv[0]);
            throw std::runtime_error("--fasta and --output are required");
        }
        if (options.fai_path.empty()) {
            options.fai_path = options.fasta_path + ".fai";
        }

        const kmer_index::IndexStats stats = kmer_index::build_index(options);
        std::cerr << "Indexed " << stats.contigs << " contigs (" << stats.bases
                  << " bases)\n"
                  << "Generated occurrences: " << stats.occurrences << '\n'
                  << "Distinct 31-mers: " << stats.distinct_kmers << '\n'
                  << "Written 31-mers: " << stats.written_kmers << " ("
                  << stats.written_locations << " locations)\n"
                  << "Filtered by max cover: " << stats.filtered_by_cover << '\n'
                  << "Filtered by max span: " << stats.filtered_by_span << '\n'
                  << "Output bytes: " << stats.output_bytes << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "kmerindex: error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
