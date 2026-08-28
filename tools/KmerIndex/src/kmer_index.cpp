#include "kmer_index.hpp"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstring>
#include <exception>
#include <fstream>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <utility>

namespace kmer_index {
namespace {

constexpr std::uint64_t kKmerMask = (std::uint64_t{1} << (2 * kKmerLength)) - 1;
std::uint64_t parse_u64(const std::string& text, const std::string& field,
                        std::size_t line_number) {
    if (text.empty()) {
        throw std::runtime_error("empty " + field + " on .fai line " +
                                 std::to_string(line_number));
    }
    std::size_t used = 0;
    unsigned long long value = 0;
    try {
        value = std::stoull(text, &used, 10);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid " + field + " on .fai line " +
                                 std::to_string(line_number) + ": " + text);
    }
    if (used != text.size()) {
        throw std::runtime_error("invalid " + field + " on .fai line " +
                                 std::to_string(line_number) + ": " + text);
    }
    return static_cast<std::uint64_t>(value);
}

int base_code(char base) {
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

void append_contig_occurrences(std::ifstream& fasta, const FaiEntry& contig,
                               std::vector<KmerOccurrence>& output) {
    output.clear();
    if (contig.length < kKmerLength) {
        return;
    }
    output.reserve(static_cast<std::size_t>(contig.length - kKmerLength + 1));

    fasta.clear();
    fasta.seekg(static_cast<std::streamoff>(contig.file_offset), std::ios::beg);
    if (!fasta) {
        throw std::runtime_error("cannot seek to contig " + contig.name +
                                 " at FASTA offset " +
                                 std::to_string(contig.file_offset));
    }

    std::vector<char> line(contig.line_bases);
    std::uint64_t remaining = contig.length;
    std::uint64_t local_position = 0;
    std::uint32_t valid_bases = 0;
    std::uint64_t forward = 0;
    std::uint64_t reverse = 0;

    while (remaining != 0) {
        const auto bases_this_line = static_cast<std::uint32_t>(
            std::min<std::uint64_t>(remaining, contig.line_bases));
        fasta.read(line.data(), static_cast<std::streamsize>(bases_this_line));
        if (fasta.gcount() != static_cast<std::streamsize>(bases_this_line)) {
            throw std::runtime_error("FASTA ended while reading contig " +
                                     contig.name + "; its .fai may be stale");
        }

        for (std::uint32_t i = 0; i < bases_this_line; ++i, ++local_position) {
            const int code = base_code(line[i]);
            if (code < 0) {
                valid_bases = 0;
                forward = 0;
                reverse = 0;
                continue;
            }

            forward = ((forward << 2) | static_cast<std::uint64_t>(code)) &
                      kKmerMask;
            reverse = (reverse >> 2) |
                      (static_cast<std::uint64_t>(3 - code) <<
                       (2 * (kKmerLength - 1)));
            if (valid_bases < kKmerLength) {
                ++valid_bases;
            }
            if (valid_bases == kKmerLength) {
                const std::uint64_t canonical = std::max(forward, reverse);
                const std::uint64_t global_position =
                    contig.global_offset + local_position - kKmerLength + 1;
                output.push_back({canonical,
                                  static_cast<std::uint32_t>(global_position)});
            }
        }

        remaining -= bases_this_line;
        if (remaining != 0) {
            const std::uint32_t newline_bytes =
                contig.line_width - contig.line_bases;
            fasta.seekg(static_cast<std::streamoff>(newline_bytes), std::ios::cur);
            if (!fasta) {
                throw std::runtime_error("cannot advance between FASTA lines for " +
                                         contig.name + "; its .fai may be stale");
            }
        }
    }
}

std::size_t contig_for_position(const std::vector<FaiEntry>& contigs,
                                std::uint32_t position) {
    const auto it = std::upper_bound(
        contigs.begin(), contigs.end(), static_cast<std::uint64_t>(position),
        [](std::uint64_t value, const FaiEntry& entry) {
            return value < entry.global_offset;
        });
    if (it == contigs.begin()) {
        throw std::runtime_error("internal error: position precedes first contig");
    }
    const std::size_t index = static_cast<std::size_t>(it - contigs.begin() - 1);
    if (position >= contigs[index].global_offset + contigs[index].length) {
        throw std::runtime_error("internal error: position is outside a contig");
    }
    return index;
}

std::uint64_t union_span(const KmerOccurrence* first,
                         const KmerOccurrence* last,
                         const std::vector<FaiEntry>& contigs,
                         std::uint64_t stop_after) {
    std::uint64_t total = 0;
    std::uint64_t interval_start = 0;
    std::uint64_t interval_end = 0;
    std::size_t current_contig = std::numeric_limits<std::size_t>::max();

    for (const KmerOccurrence* occurrence = first; occurrence != last;
         ++occurrence) {
        const std::size_t contig_index =
            contig_for_position(contigs, occurrence->position);
        const FaiEntry& contig = contigs[contig_index];
        const std::uint64_t contig_start = contig.global_offset;
        const std::uint64_t contig_end = contig_start + contig.length;
        const std::uint64_t position = occurrence->position;
        const std::uint64_t start =
            position > contig_start + kSpanFlank
                ? position - kSpanFlank
                : contig_start;
        const std::uint64_t end = std::min(contig_end, position + kSpanFlank);

        if (contig_index != current_contig || start > interval_end) {
            if (current_contig != std::numeric_limits<std::size_t>::max()) {
                total += interval_end - interval_start;
                if (total > stop_after) {
                    return total;
                }
            }
            current_contig = contig_index;
            interval_start = start;
            interval_end = end;
        } else {
            interval_end = std::max(interval_end, end);
        }
    }
    if (current_contig != std::numeric_limits<std::size_t>::max()) {
        total += interval_end - interval_start;
    }
    return total;
}

class BinaryWriter {
public:
    explicit BinaryWriter(const std::string& path) : stream_(path, std::ios::binary) {
        if (!stream_) {
            throw std::runtime_error("cannot open output file " + path + ": " +
                                     std::strerror(errno));
        }
        buffer_.reserve(1 << 20);
    }

    ~BinaryWriter() {
        try {
            flush();
        } catch (...) {
        }
    }

    void put_u32(std::uint32_t value) {
        for (int shift = 0; shift < 32; shift += 8) {
            buffer_.push_back(static_cast<char>((value >> shift) & 0xff));
        }
        maybe_flush();
    }

    void put_u64(std::uint64_t value) {
        for (int shift = 0; shift < 64; shift += 8) {
            buffer_.push_back(static_cast<char>((value >> shift) & 0xff));
        }
        maybe_flush();
    }

    void finish() {
        flush();
        stream_.flush();
        if (!stream_) {
            throw std::runtime_error("failed while finalizing binary output");
        }
    }

private:
    void maybe_flush() {
        if (buffer_.size() >= (1 << 20)) {
            flush();
        }
    }

    void flush() {
        if (buffer_.empty()) {
            return;
        }
        stream_.write(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
        if (!stream_) {
            throw std::runtime_error("failed while writing binary output");
        }
        buffer_.clear();
    }

    std::ofstream stream_;
    std::vector<char> buffer_;
};

}  // namespace

std::vector<FaiEntry> read_fai(const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open FASTA index " + path);
    }

    std::vector<FaiEntry> entries;
    std::string line;
    std::uint64_t global_offset = 0;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }

        std::vector<std::string> fields;
        std::size_t start = 0;
        while (true) {
            const std::size_t tab = line.find('\t', start);
            fields.push_back(line.substr(start, tab - start));
            if (tab == std::string::npos) {
                break;
            }
            start = tab + 1;
        }
        if (fields.size() < 5 || fields[0].empty()) {
            throw std::runtime_error("expected at least five tab-separated fields "
                                     "on .fai line " +
                                     std::to_string(line_number));
        }

        FaiEntry entry;
        entry.name = fields[0];
        entry.length = parse_u64(fields[1], "length", line_number);
        entry.file_offset = parse_u64(fields[2], "file offset", line_number);
        const std::uint64_t line_bases =
            parse_u64(fields[3], "line bases", line_number);
        const std::uint64_t line_width =
            parse_u64(fields[4], "line width", line_number);
        if (line_bases > std::numeric_limits<std::uint32_t>::max() ||
            line_width > std::numeric_limits<std::uint32_t>::max()) {
            throw std::runtime_error("FASTA line width exceeds uint32 on .fai line " +
                                     std::to_string(line_number));
        }
        entry.line_bases = static_cast<std::uint32_t>(line_bases);
        entry.line_width = static_cast<std::uint32_t>(line_width);
        if (entry.length != 0 &&
            (entry.line_bases == 0 || entry.line_width < entry.line_bases)) {
            throw std::runtime_error("invalid FASTA line geometry on .fai line " +
                                     std::to_string(line_number));
        }
        if (entry.length > std::numeric_limits<std::uint32_t>::max() -
                               global_offset) {
            throw std::runtime_error(
                "total assembly length must be less than 2^32 bases");
        }
        entry.global_offset = global_offset;
        global_offset += entry.length;
        entries.push_back(std::move(entry));
    }
    if (entries.empty()) {
        throw std::runtime_error("FASTA index is empty: " + path);
    }
    return entries;
}

std::uint64_t encode_kmer(const std::string& sequence) {
    if (sequence.size() != kKmerLength) {
        throw std::invalid_argument("encode_kmer requires exactly 31 bases");
    }
    std::uint64_t forward = 0;
    std::uint64_t reverse = 0;
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        const int code = base_code(sequence[i]);
        if (code < 0) {
            throw std::invalid_argument("encode_kmer accepts only A, C, G, and T");
        }
        forward = (forward << 2) | static_cast<std::uint64_t>(code);
        reverse |= static_cast<std::uint64_t>(3 - code)
                   << (2 * i);
    }
    return std::max(forward, reverse);
}

IndexStats build_index(const Options& options) {
    const std::vector<FaiEntry> contigs = read_fai(options.fai_path);
    IndexStats stats;
    stats.contigs = contigs.size();

    std::uint64_t occurrence_upper_bound = 0;
    std::uint64_t largest_possible_position = 0;
    bool has_possible_kmer = false;
    for (const FaiEntry& contig : contigs) {
        stats.bases += contig.length;
        if (contig.length >= kKmerLength) {
            occurrence_upper_bound += contig.length - kKmerLength + 1;
            largest_possible_position = std::max(
                largest_possible_position,
                contig.global_offset + contig.length - kKmerLength);
            has_possible_kmer = true;
        }
    }
    if (occurrence_upper_bound > std::numeric_limits<std::size_t>::max()) {
        throw std::runtime_error("number of k-mers exceeds this system's size_t");
    }

    if (options.max_cover >= 2 && has_possible_kmer) {
        const std::uint32_t lowest_marker =
            std::numeric_limits<std::uint32_t>::max() - options.max_cover;
        if (largest_possible_position >= lowest_marker) {
            throw std::runtime_error(
                "assembly coordinates overlap the binary count-marker range; "
                "reduce --maxcover or use an assembly shorter than " +
                std::to_string(lowest_marker) + " bases");
        }
    }

    std::vector<KmerOccurrence> occurrences;
    occurrences.reserve(static_cast<std::size_t>(occurrence_upper_bound));
    std::mutex occurrence_mutex;
    std::mutex error_mutex;
    std::exception_ptr worker_error;
    std::atomic<std::size_t> next_contig{0};
    std::atomic<bool> stop{false};

    const std::size_t thread_count = std::max<std::size_t>(
        1, std::min<std::size_t>(options.threads, contigs.size()));
    std::vector<std::thread> workers;
    workers.reserve(thread_count);
    for (std::size_t worker_index = 0; worker_index < thread_count;
         ++worker_index) {
        workers.emplace_back([&] {
            try {
                std::ifstream fasta(options.fasta_path, std::ios::binary);
                if (!fasta) {
                    throw std::runtime_error("cannot open FASTA " +
                                             options.fasta_path);
                }
                std::vector<KmerOccurrence> local;
                while (!stop.load(std::memory_order_relaxed)) {
                    const std::size_t index =
                        next_contig.fetch_add(1, std::memory_order_relaxed);
                    if (index >= contigs.size()) {
                        break;
                    }
                    append_contig_occurrences(fasta, contigs[index], local);
                    {
                        std::lock_guard<std::mutex> lock(occurrence_mutex);
                        occurrences.insert(occurrences.end(), local.begin(),
                                           local.end());
                    }
                    std::vector<KmerOccurrence>().swap(local);
                }
            } catch (...) {
                stop.store(true, std::memory_order_relaxed);
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!worker_error) {
                    worker_error = std::current_exception();
                }
            }
        });
    }
    for (std::thread& worker : workers) {
        worker.join();
    }
    if (worker_error) {
        std::rethrow_exception(worker_error);
    }

    stats.occurrences = occurrences.size();
    std::sort(occurrences.begin(), occurrences.end(),
              [](const KmerOccurrence& left, const KmerOccurrence& right) {
                  if (left.kmer != right.kmer) {
                      return left.kmer < right.kmer;
                  }
                  return left.position < right.position;
              });

    BinaryWriter writer(options.output_path);
    std::size_t group_start = 0;
    while (group_start < occurrences.size()) {
        std::size_t group_end = group_start + 1;
        while (group_end < occurrences.size() &&
               occurrences[group_end].kmer == occurrences[group_start].kmer) {
            ++group_end;
        }
        ++stats.distinct_kmers;
        const std::size_t count = group_end - group_start;
        if (count > options.max_cover) {
            ++stats.filtered_by_cover;
            group_start = group_end;
            continue;
        }
        const std::uint64_t span = union_span(
            occurrences.data() + group_start, occurrences.data() + group_end,
            contigs, options.max_span);
        if (span > options.max_span) {
            ++stats.filtered_by_span;
            group_start = group_end;
            continue;
        }

        writer.put_u64(occurrences[group_start].kmer);
        if (count == 1) {
            writer.put_u32(occurrences[group_start].position);
            stats.output_bytes += sizeof(std::uint64_t) + sizeof(std::uint32_t);
        } else {
            const std::uint32_t marker =
                std::numeric_limits<std::uint32_t>::max() -
                static_cast<std::uint32_t>(count);
            writer.put_u32(marker);
            for (std::size_t i = group_start; i < group_end; ++i) {
                writer.put_u32(occurrences[i].position);
            }
            stats.output_bytes += sizeof(std::uint64_t) + sizeof(std::uint32_t) +
                                  count * sizeof(std::uint32_t);
        }
        ++stats.written_kmers;
        stats.written_locations += count;
        group_start = group_end;
    }
    writer.finish();
    return stats;
}

}  // namespace kmer_index
