#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>

constexpr int KMER_SIZE = 31;
constexpr uint64_t KMER_MASK = (1ULL << (2 * (KMER_SIZE - 1)));

char complement(char base) {
    switch (base) {
        case 'A': return 'T'; case 'T': return 'A';
        case 'C': return 'G'; case 'G': return 'C';
        case 'a': return 't'; case 't': return 'a';
        case 'c': return 'g'; case 'g': return 'c';
        default: return 'N';
    }
}

std::string reverse_complement(const std::string& seq) {
    std::string rc(seq.rbegin(), seq.rend());
    std::transform(rc.begin(), rc.end(), rc.begin(), complement);
    return rc;
}

int base_to_int(char base) {
    switch (std::toupper(base)) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

class KmerAnnotation {
public:
    std::unordered_map<uint64_t, char> kmersign;
    int annotate(const std::string& seq, int totalsign_hint) {
        int totalsign = totalsign_hint;
        uint64_t forward_k = 0, reverse_k = 0;
        int size = 0;
        
        static std::vector<long long> unsigned_kmers(20000000);
        uint32_t unsigned_number = 0;
        
        for (char base : seq) {
            int val = base_to_int(base);
            if (val == -1) {
                size = 0;
                forward_k = reverse_k = 0;
                continue;
            }
            
            if (size >= KMER_SIZE) {
                forward_k = ((forward_k % KMER_MASK) << 2) + val;
                reverse_k = (reverse_k >> 2) + ((3 - val) * KMER_MASK);
            } else {
                forward_k = (forward_k << 2) + val;
                reverse_k += (3 - val) << (2 * size);
                size++;
            }
            
            if (size >= KMER_SIZE) {
                uint64_t canonical = std::max(forward_k, reverse_k);
                int strand = (canonical == forward_k) ? 1 : -1;
                if (kmersign.count(canonical))
                {
                    totalsign += strand * kmersign[canonical];
                }
                else
                {
                    unsigned_kmers[unsigned_number++] = ((long long) canonical) * strand;
                }
                   
            }
        }
        
        totalsign = (totalsign >= 0) ? std::max(1, totalsign) : totalsign;
        for (int index = 0; index < unsigned_number; index++)
        {
            int strand = (unsigned_kmers[index] >= 0) ? 1 : -1;
            kmersign[abs(unsigned_kmers[index])] = strand * ((totalsign >= 0) ? 1 : -1);
        }
        
        return totalsign;
    }
};

struct FastaEntry {
    std::string name;
    std::string header;
    std::string seq;
};

std::vector<FastaEntry> read_fasta(const std::string& filename) {
    std::ifstream fin(filename);
    std::vector<FastaEntry> entries;
    std::string line, current_name, current_header, current_seq;
    
    while (std::getline(fin, line)) {
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

void write_fasta(const std::string& filename, const std::vector<FastaEntry>& entries) {
    std::ofstream fout(filename);
    for (const auto& entry : entries) {
        fout << ">" << entry.header << "\n";
        fout << entry.seq << "\n";
    }
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: ./kmer_annotator input.fasta output.fasta\n";
        return 1;
    }
    
    std::string input_file, output_file;

    for (int i = 1; i < argc - 1; ++i) {
        std::string arg = argv[i];
        if (arg == "-i") {
            input_file = argv[++i];
        } else if (arg == "-o") {
            output_file = argv[++i];
        }
    }
    
    auto entries = read_fasta(input_file);
    KmerAnnotation annotator;
    
    std::map<std::string, int> contig_sign;
    std::unordered_map<std::string, std::string> contig_from_name;
    std::map<std::string, int> contig_size;
    
    for (const auto& entry : entries) {
        std::string contig = entry.header.substr(entry.header.find(" ") + 1);
        contig = contig.substr(0, contig.find(":"));
        contig_from_name[entry.name] = contig;
        contig_size[contig] += entry.seq.size();
    }
    
    std::vector<FastaEntry> sorted_entries = entries;
    std::sort(sorted_entries.begin(), sorted_entries.end(), [&](const FastaEntry& a, const FastaEntry& b) {
        const std::string& contig_a = contig_from_name[a.name];
        const std::string& contig_b = contig_from_name[b.name];

        int contig_size_a = contig_size[contig_a];
        int contig_size_b = contig_size[contig_b];

        if (contig_size_a != contig_size_b) {
            return contig_size_a > contig_size_b;
        }
        return a.seq.size() > b.seq.size();
    });
    
    std::map<std::string, int> name_sign;
    for (auto& entry : sorted_entries) {
        std::string contig = contig_from_name[entry.name];
        int sign = annotator.annotate(entry.seq, std::min(500, contig_sign[contig] / 100));
        contig_sign[contig] += sign;
        name_sign[entry.name] = sign;
    }
    
    for (auto& entry : entries) {
        int sign = name_sign[entry.name];
        std::string& header = entry.header;
        std::string strand = "+";

        // Parse header into fields
        std::istringstream iss(header);
        std::vector<std::string> fields;
        std::string token;
        while (iss >> token) {
            fields.push_back(token);
        }

        if (fields.size() >= 2) {
            // Remove old strand from second field if exists
            if (!fields[1].empty() && (fields[1].back() == '+' || fields[1].back() == '-')) {
                strand = std::string(1, fields[1].back());  // remember old strand
                fields[1].pop_back();
            }

            // Adjust strand based on sign
            if (sign < 0) {
                entry.seq = reverse_complement(entry.seq);
                strand = (strand == "-" ? "+" : "-");
            } else {
                strand = (strand == "-" ? "-" : "+");
            }

            // Append new strand to second field
            fields[1] += strand;
        }

        // Reconstruct header
        std::ostringstream oss;
        for (size_t i = 0; i < fields.size(); ++i) {
            if (i > 0) oss << '\t';
            oss << fields[i];
        }
        entry.header = oss.str();
    }
    
    write_fasta(output_file, entries);
    return 0;
}

