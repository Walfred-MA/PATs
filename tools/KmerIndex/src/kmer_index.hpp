#ifndef KMER_INDEX_HPP
#define KMER_INDEX_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace kmer_index {

constexpr std::uint32_t kKmerLength = 31;
constexpr std::uint32_t kSpanFlank = 10'000;

#pragma pack(push, 1)
struct KmerOccurrence {
    std::uint64_t kmer;
    std::uint32_t position;
};
#pragma pack(pop)

static_assert(sizeof(KmerOccurrence) == 12,
              "KmerOccurrence must occupy exactly 12 bytes");

struct FaiEntry {
    std::string name;
    std::uint64_t length = 0;
    std::uint64_t file_offset = 0;
    std::uint32_t line_bases = 0;
    std::uint32_t line_width = 0;
    std::uint64_t global_offset = 0;
};

struct Options {
    std::string fasta_path;
    std::string fai_path;
    std::string output_path;
    std::uint32_t threads = 1;
    std::uint16_t max_cover = 255;
    std::uint64_t max_span = 1'000'000;
};

struct IndexStats {
    std::uint64_t contigs = 0;
    std::uint64_t bases = 0;
    std::uint64_t occurrences = 0;
    std::uint64_t distinct_kmers = 0;
    std::uint64_t written_kmers = 0;
    std::uint64_t written_locations = 0;
    std::uint64_t filtered_by_cover = 0;
    std::uint64_t filtered_by_span = 0;
    std::uint64_t output_bytes = 0;
};

std::vector<FaiEntry> read_fai(const std::string& path);
IndexStats build_index(const Options& options);
std::uint64_t encode_kmer(const std::string& sequence);

}  // namespace kmer_index

#endif
