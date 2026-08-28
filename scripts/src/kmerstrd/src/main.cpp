#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

constexpr int KMER_SIZE = 31;

// Keep last (KMER_SIZE-1) bases when shifting in a new base (2 bits each).
constexpr int  DROP_BASES = KMER_SIZE - 1;
constexpr int  TOP_SHIFT  = 2 * (KMER_SIZE - 1);
constexpr uint64_t FORWARD_MOD  = (1ULL << (2 * DROP_BASES));      // 2^(2*(K-1))
constexpr uint64_t FORWARD_MASK = (FORWARD_MOD - 1ULL);            // low bits mask

static inline char complement(char b) {
    switch (b) {
        case 'A': return 'T'; case 'T': return 'A';
        case 'C': return 'G'; case 'G': return 'C';
        case 'a': return 't'; case 't': return 'a';
        case 'c': return 'g'; case 'g': return 'c';
        default:  return 'N';
    }
}

static inline std::string reverse_complement(const std::string& seq) {
    std::string rc(seq.rbegin(), seq.rend());
    std::transform(rc.begin(), rc.end(), rc.begin(), complement);
    return rc;
}

static inline int base_to_int(char base) {
    switch (std::toupper(static_cast<unsigned char>(base))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return -1;
    }
}

struct FastaEntry {
    std::string name;   // first token after '>'
    std::string header; // full header without leading '>'
    std::string seq;    // concatenated sequence
};

static std::vector<FastaEntry> read_fasta(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin) throw std::runtime_error("Failed to open input FASTA: " + filename);

    std::vector<FastaEntry> entries;
    std::string line, current_name, current_header, current_seq;

    while (std::getline(fin, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back(); // CRLF safety
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!current_name.empty()) {
                entries.push_back({current_name, current_header, current_seq});
                current_seq.clear();
            }
            std::istringstream iss(line.substr(1));
            iss >> current_name;
            current_header = line.substr(1);
        } else {
            current_seq += line;
        }
    }

    if (!current_name.empty()) {
        entries.push_back({current_name, current_header, current_seq});
    }
    return entries;
}

static std::string trim_copy(const std::string& text) {
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

static std::vector<std::string> split_priority_patterns(const std::string& text) {
    std::vector<std::string> patterns;
    std::stringstream stream(text);
    std::string pattern;
    while (std::getline(stream, pattern, ',')) {
        pattern = trim_copy(pattern);
        if (!pattern.empty()) patterns.push_back(std::move(pattern));
    }
    return patterns;
}

static size_t sequence_priority(
    const std::string& name,
    const std::vector<std::string>& patterns
) {
    for (size_t i = 0; i < patterns.size(); ++i) {
        if (name.find(patterns[i]) != std::string::npos) return i;
    }
    return patterns.size();
}

static void write_fasta(const std::string& filename, const std::vector<FastaEntry>& entries) {
    std::ofstream fout(filename);
    if (!fout) throw std::runtime_error("Failed to open output FASTA: " + filename);

    for (const auto& e : entries) {
        fout << ">" << e.header << "\n";
        fout << e.seq << "\n";
    }
}

// Parse contig from header: take 2nd whitespace token, strip trailing +/- if present, then take prefix before ':'
static std::string parse_contig_from_header(const std::string& header) {
    std::istringstream iss(header);
    std::string tok1, tok2;
    if (!(iss >> tok1)) return "NA";
    if (!(iss >> tok2)) return "NA";

    if (!tok2.empty() && (tok2.back() == '+' || tok2.back() == '-')) tok2.pop_back();
    auto p = tok2.find(':');
    if (p != std::string::npos) tok2.resize(p);
    if (tok2.empty()) return "NA";
    return tok2;
}

class KmerAnnotation {
public:
    // Signed 8-bit: deterministic across platforms (unlike plain char).
    std::unordered_map<uint64_t, int8_t> kmersign;

    // scratch for "unknown" kmers in annotate(): store signed canonical (sign encodes strand)
    std::vector<int64_t> scratch;

    int annotate(const std::string& seq, int totalsign_hint) {
        int totalsign = totalsign_hint;

        uint64_t forward_k = 0, reverse_k = 0;
        int have = 0;

        scratch.clear();
        scratch.reserve(seq.size()); // upper bound; okay if large

        for (char base : seq) {
            int val = base_to_int(base);
            if (val < 0) {
                have = 0;
                forward_k = reverse_k = 0;
                continue;
            }

            if (have >= KMER_SIZE) {
                // keep last 2*(K-1) bits then shift/add new base
                forward_k = ((forward_k & FORWARD_MASK) << 2) | (uint64_t)val;
                reverse_k = (reverse_k >> 2) | ((uint64_t)(3 - val) << TOP_SHIFT);
            } else {
                forward_k = (forward_k << 2) | (uint64_t)val;
                reverse_k |= (uint64_t)(3 - val) << (2 * have);
                ++have;
            }

            if (have >= KMER_SIZE) {
                uint64_t canonical = (forward_k >= reverse_k) ? forward_k : reverse_k;
                int strand = (canonical == forward_k) ? 1 : -1;

                auto it = kmersign.find(canonical);
                if (it != kmersign.end()) {
                    totalsign += strand * (int)it->second;
                } else {
                    scratch.push_back((int64_t)canonical * (int64_t)strand);
                }
            }
        }

        // same behavior as your original code
        totalsign = (totalsign >= 0) ? std::max(1, totalsign) : totalsign;

        const int8_t assign_sign = (totalsign >= 0) ? (int8_t)1 : (int8_t)-1;
        for (int64_t v : scratch) {
            int strand = (v >= 0) ? 1 : -1;
            uint64_t key = (v >= 0) ? (uint64_t)v : (uint64_t)(-v);
            kmersign[key] = (int8_t)(strand * assign_sign);
        }

        return totalsign;
    }

    int annotate_passive(const std::string& seq, int totalsign_hint) const {
        int totalsign = totalsign_hint;

        uint64_t forward_k = 0, reverse_k = 0;
        int have = 0;

        for (char base : seq) {
            int val = base_to_int(base);
            if (val < 0) {
                have = 0;
                forward_k = reverse_k = 0;
                continue;
            }

            if (have >= KMER_SIZE) {
                forward_k = ((forward_k & FORWARD_MASK) << 2) | (uint64_t)val;
                reverse_k = (reverse_k >> 2) | ((uint64_t)(3 - val) << TOP_SHIFT);
            } else {
                forward_k = (forward_k << 2) | (uint64_t)val;
                reverse_k |= (uint64_t)(3 - val) << (2 * have);
                ++have;
            }

            if (have >= KMER_SIZE) {
                uint64_t canonical = (forward_k >= reverse_k) ? forward_k : reverse_k;
                int strand = (canonical == forward_k) ? 1 : -1;

                auto it = kmersign.find(canonical);
                if (it != kmersign.end()) {
                    totalsign += strand * (int)it->second;
                }
            }
        }

        totalsign = (totalsign >= 0) ? std::max(1, totalsign) : totalsign;
        return totalsign;
    }

    void kmerassign(const std::string& seq, int inc = 1) {
        uint64_t forward_k = 0, reverse_k = 0;
        int have = 0;

        for (char base : seq) {
            int val = base_to_int(base);
            if (val < 0) {
                have = 0;
                forward_k = reverse_k = 0;
                continue;
            }

            if (have >= KMER_SIZE) {
                forward_k = ((forward_k & FORWARD_MASK) << 2) | (uint64_t)val;
                reverse_k = (reverse_k >> 2) | ((uint64_t)(3 - val) << TOP_SHIFT);
            } else {
                forward_k = (forward_k << 2) | (uint64_t)val;
                reverse_k |= (uint64_t)(3 - val) << (2 * have);
                ++have;
            }

            if (have >= KMER_SIZE) {
                uint64_t canonical = (forward_k >= reverse_k) ? forward_k : reverse_k;
                int strand = (canonical == forward_k) ? 1 : -1;

                int delta = inc * strand;
                int cur = 0;
                auto it = kmersign.find(canonical);
                if (it != kmersign.end()) cur = (int)it->second;

                int next = cur + delta;
                if (next > 127) next = 127;
                if (next < -127) next = -127;

                kmersign[canonical] = (int8_t)next;
            }
        }
    }
};

static void usage() {
    std::cerr
        << "Usage: kmerstrd -i input.fasta -o output.fasta "
        << "[-r ref.fasta] [-p pattern1,pattern2,...]\n"
        << "  -p, --priority  Prioritize input records by case-sensitive substring matches\n"
        << "                  against the first FASTA header token. Earlier patterns win.\n";
}

int main(int argc, char** argv) {
    std::string input_file, output_file, ref_file, priority_text;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-i" || arg == "--in") && i + 1 < argc) input_file = argv[++i];
        else if ((arg == "-o" || arg == "--out") && i + 1 < argc) output_file = argv[++i];
        else if ((arg == "-r" || arg == "--ref") && i + 1 < argc) ref_file = argv[++i];
        else if ((arg == "-p" || arg == "--priority") && i + 1 < argc) priority_text = argv[++i];
        else if (arg == "-h" || arg == "--help") { usage(); return 0; }
        else { std::cerr << "Unknown/invalid arg: " << arg << "\n"; usage(); return 1; }
    }

    if (input_file.empty() || output_file.empty()) {
        usage();
        return 1;
    }

    KmerAnnotation annotator;

    // Optional reference: learn kmer directionality
    if (!ref_file.empty()) {
        auto refs = read_fasta(ref_file);
        for (const auto& e : refs) annotator.kmerassign(e.seq, 1);

        // normalize to {-1,0,+1}
        for (auto& kv : annotator.kmersign) {
            if (kv.second > 0) kv.second = (int8_t)1;
            else if (kv.second < 0) kv.second = (int8_t)-1;
            else kv.second = (int8_t)0;
        }
    }

    auto entries = read_fasta(input_file);
    const auto priority_patterns = split_priority_patterns(priority_text);

    // contig bookkeeping
    std::map<std::string, int> contig_sign;
    std::unordered_map<std::string, std::string> contig_from_name;
    std::map<std::string, int64_t> contig_size;

    contig_from_name.reserve(entries.size());

    for (const auto& e : entries) {
        std::string contig = parse_contig_from_header(e.header);
        contig_from_name[e.name] = contig;
        contig_size[contig] += (int64_t)e.seq.size();
    }

    // Priority patterns are primary; retain the previous size ordering within each rank.
    std::vector<FastaEntry> sorted = entries;
    std::sort(sorted.begin(), sorted.end(), [&](const FastaEntry& a, const FastaEntry& b) {
        const size_t pa = sequence_priority(a.name, priority_patterns);
        const size_t pb = sequence_priority(b.name, priority_patterns);
        if (pa != pb) return pa < pb;

        const std::string& ca = contig_from_name[a.name];
        const std::string& cb = contig_from_name[b.name];

        int64_t sa = contig_size[ca];
        int64_t sb = contig_size[cb];

        if (sa != sb) return sa > sb;
        return a.seq.size() > b.seq.size();
    });

    std::unordered_map<std::string, int> name_sign;
    name_sign.reserve(entries.size());

    // annotate all but last with learning, last with passive
    for (size_t i = 0; i + 1 < sorted.size(); ++i) {
        const auto& e = sorted[i];
        const std::string& contig = contig_from_name[e.name];
        int hint = std::min(500, contig_sign[contig] / 100);
        int s = annotator.annotate(e.seq, hint);
        contig_sign[contig] += s;
        name_sign[e.name] = s;
    }
    if (!sorted.empty()) {
        const auto& e = sorted.back();
        const std::string& contig = contig_from_name[e.name];
        int hint = std::min(500, contig_sign[contig] / 100);
        int s = annotator.annotate_passive(e.seq, hint);
        contig_sign[contig] += s;
        name_sign[e.name] = s;
    }

    // Apply flips + fix strand marker in header deterministically.
    for (auto& e : entries) {
        int sign = 0;
        auto it = name_sign.find(e.name);
        if (it != name_sign.end()) sign = it->second;

        // Tokenize header by whitespace
        std::istringstream iss(e.header);
        std::vector<std::string> fields;
        std::string tok;
        while (iss >> tok) {
            if (!tok.empty() && tok.back() == '\r') tok.pop_back();
            fields.push_back(std::move(tok));
        }
        if (fields.empty()) continue;

        bool had_strand = false;
        char old_strand = '+';

        if (fields.size() >= 2 && !fields[1].empty()) {
            char last = fields[1].back();
            if (last == '+' || last == '-') {
                had_strand = true;
                old_strand = last;
                fields[1].pop_back();
            }
        }

        if (sign < 0) {
            // flip sequence
            e.seq = reverse_complement(e.seq);

            // toggle strand marker if we will write one
            char new_strand = (old_strand == '+') ? '-' : '+';

            // only force-write strand if it existed OR we flipped
            if (fields.size() >= 2) fields[1].push_back(new_strand);
            else fields[0].push_back(new_strand);

        } else {
            // no flip: preserve strand if it existed; otherwise do not modify header by adding '+'
            if (had_strand) {
                if (fields.size() >= 2) fields[1].push_back(old_strand);
                else fields[0].push_back(old_strand);
            }
        }

        // Reconstruct header with single spaces (deterministic)
        std::ostringstream oss;
        for (size_t i = 0; i < fields.size(); ++i) {
            if (i) oss << ' ';
            oss << fields[i];
        }
        e.header = oss.str();
    }

    write_fasta(output_file, entries);
    return 0;
}

