//
//  KmerAssigner.hpp
//  KmerPartition
//
//  Created by walfred on 6/1/26.
//

#ifndef KmerAssigner_hpp
#define KmerAssigner_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <thread>
#include <atomic>
#include <mutex>
#include <sstream>
#include <math.h>
#include <numeric>
#include <cstdint>
#include <cctype>
#include <filesystem>
#include <limits>
#include <cstdlib>
#include <future>
#include "fasta.hpp"
#include "KmerStruct.hpp"

extern int cutoff;
extern int nthreads;
extern bool ifmultipleKmer;

using kmer64_dict = std::unordered_map<u128, uint, hash_128> ;
using kmer32_dict = std::unordered_map<ull, uint > ;

using kmer64_dict_nt = std::unordered_map<u128, uint8, hash_128> ;
using kmer32_dict_nt = std::unordered_map<ull, uint8>;

using kmer64_dict_mul = std::unordered_map<u128, uint*, hash_128> ;
using kmer32_dict_mul = std::unordered_map<ull, uint* > ;

template <int dictsize>
class kmer_assigner
{
    kmer32_dict_nt kmer_hash;

    std::unordered_map<std::string, uint8> group_index;
    std::unordered_map<std::string, vector<uint8>> sequence_groups;
    vector<std::string> group_names = {""};
    vector<vector<uint8>> sequence_groups_by_order;
public:
    
    kmer_assigner()
    {};
    ~kmer_assigner()
    {};

    template <class typefile>
    void read_kmertarget(typefile &fastafile);
    
    template <class typefile>
    void read_fasta(typefile &fastafile);
    
    void assigngroup(string &inputfile, string &kmerfile, string &outputfolder);
    void assigngroup_mul(string &inputfile, vector<std::string> &kmerfiles, string &outputfile);
};

inline static std::string assigner_trim(const std::string &text)
{
    size_t start = 0;
    while (start < text.size() && std::isspace(static_cast<unsigned char>(text[start]))) start++;

    size_t end = text.size();
    while (end > start && std::isspace(static_cast<unsigned char>(text[end - 1]))) end--;

    return text.substr(start, end - start);
}

template <typename T>
static void assigner_kmer_read_c(char base, int &current_size, T &current_kmer, T &reverse_kmer)
{
    int converted = 0;
    T reverse_converted;

    if (base == '\n' || base == ' ') return;

    if (base_to_int(base, converted))
    {
        current_kmer <<= (8 * sizeof(current_kmer) - 2 * klen + 2);
        current_kmer >>= (8 * sizeof(current_kmer) - 2 * klen);
        current_kmer += converted;

        reverse_kmer >>= 2;
        reverse_converted = 0b11 - converted;
        reverse_converted <<= (2 * klen - 2);
        reverse_kmer += reverse_converted;

        current_size++;
    }
    else
    {
        current_size = 0;
        current_kmer = 0;
        reverse_kmer = 0;
    }
}

inline static bool assigner_has_lowercase_base(const std::string &line)
{
    for (char base : line)
    {
        switch (base)
        {
            case 'a': case 't': case 'c': case 'g':
                return true;
            default:
                break;
        }
    }

    return false;
}

inline static std::string assigner_group_to_bed_score(const std::string &group_name, uint8 group_index)
{
    if (group_name.size() > 1 && (group_name[0] == 'p' || group_name[0] == 'P') &&
        std::isdigit(static_cast<unsigned char>(group_name[1])))
    {
        return group_name.substr(1);
    }
    
    return std::to_string(group_index);
}

inline static std::string assigner_true_filename(const std::string &path)
{
    return std::filesystem::path(path).filename().string();
}

inline static std::string assigner_prefix_from_kmer_path(const std::string &path)
{
    std::string filename = assigner_true_filename(path);
    size_t underscore = filename.find('_');
    if (underscore != std::string::npos) return filename.substr(0, underscore);
    return std::filesystem::path(filename).stem().string();
}

inline static std::string assigner_matrix_candidate_from_text(const std::string &text)
{
    std::string filename = assigner_true_filename(text);
    size_t end = filename.find_first_of("_:|,;");
    if (end == std::string::npos) return filename;
    return filename.substr(0, end);
}

inline static std::string assigner_find_bed_matrix(
    const std::vector<std::string> &fields,
    const std::unordered_map<std::string, std::string> &kmer_by_prefix)
{
    auto try_field = [&](const std::string &field) -> std::string
    {
        auto exact = kmer_by_prefix.find(field);
        if (exact != kmer_by_prefix.end()) return exact->first;
        
        std::string candidate = assigner_matrix_candidate_from_text(field);
        auto find = kmer_by_prefix.find(candidate);
        if (find != kmer_by_prefix.end()) return find->first;
        
        std::string parent = std::filesystem::path(field).parent_path().filename().string();
        find = kmer_by_prefix.find(parent);
        if (find != kmer_by_prefix.end()) return find->first;
        
        return "";
    };
    
    if (fields.size() > 3)
    {
        std::string matched = try_field(fields[3]);
        if (!matched.empty()) return matched;
    }
    
    for (const auto &field : fields)
    {
        std::string matched = try_field(field);
        if (!matched.empty()) return matched;
    }
    
    return "";
}

class bed_subset_fasta
{
    const std::vector<std::string> &bed_lines;
    const std::vector<size_t> &bed_indexes;
    size_t row_index = 0;
    std::vector<std::string> pending_lines;
    size_t pending_index = 0;
    
public:
    bed_subset_fasta(const std::vector<std::string> &lines, const std::vector<size_t> &indexes):
        bed_lines(lines), bed_indexes(indexes)
    {}
    
    bool nextLine(std::string &StrLine)
    {
        if (pending_index < pending_lines.size())
        {
            StrLine = pending_lines[pending_index++];
            return true;
        }
        
        pending_lines.clear();
        pending_index = 0;
        
        while (row_index < bed_indexes.size())
        {
            size_t bed_index = bed_indexes[row_index++];
            std::vector<std::string> fields;
            if (bed_index >= bed_lines.size() || !parse_bed_record_fields(bed_lines[bed_index], fields)) continue;
            
            unsigned long long start = std::stoull(fields[1]);
            unsigned long long end = std::stoull(fields[2]);
            std::string seq = fetch_query_sequence(fields[0], start, end);
            
            pending_lines.push_back(">" + bed_lines[bed_index]);
            for (size_t i = 0; i < seq.size(); i += 80)
            {
                pending_lines.push_back(seq.substr(i, 80));
            }
            
            StrLine = pending_lines[pending_index++];
            return true;
        }
        
        return false;
    }
    
    void Close()
    {
        pending_lines.clear();
        pending_index = 0;
    }
};

template <int dictsize>
template <class typefile>
void kmer_assigner<dictsize>::read_kmertarget(typefile &fastafile)
{
    std::string StrLine;
    uint8 current_group = 0;
    size_t line_number = 0;

    while (fastafile.nextLine(StrLine))
    {
        line_number++;
        if (StrLine.empty()) continue;
        if (line_number == 1 && StrLine[0] == '#') continue;

        switch (StrLine[0])
        {
            case '>':
            {
                std::string partition_name = assigner_trim(StrLine.substr(1));

                if (partition_name == "p0")
                {
                    current_group = 0;
                    continue;
                }

                auto find = group_index.find(partition_name);
                if (find == group_index.end())
                {
                    if (group_names.size() > std::numeric_limits<uint8>::max())
                    {
                        std::cerr << "ERROR: Too many predefined k-mer groups; uint8 supports at most 255 saved groups.\n";
                        std::_Exit(EXIT_FAILURE);
                    }

                    current_group = static_cast<uint8>(group_names.size());
                    group_index[partition_name] = current_group;
                    group_names.push_back(partition_name);
                }
                else
                {
                    current_group = find->second;
                }

                continue;
            }
            case '+': case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }

        if (current_group == 0) continue;
        if (ifmask && assigner_has_lowercase_base(StrLine)) continue;

        int current_size = 0;
        ull current_kmer = 0;
        ull reverse_kmer = 0;

        for (char base : StrLine)
        {
            if (base == '\0') break;
            if (base == '\n' || base == ' ') continue;

            assigner_kmer_read_c(base, current_size, current_kmer, reverse_kmer);

            if (current_size < klen) continue;

            ull larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
            kmer_hash[larger_kmer] = current_group;
        }
    }

    fastafile.Close();
}

template <int dictsize>
template <class typefile>
void kmer_assigner<dictsize>::read_fasta(typefile &fastafile)
{
    std::string StrLine;
    std::string current_name;
    vector<int> group_counts(std::numeric_limits<uint8>::max() + 1, 0);
    int current_size = 0;
    ull current_kmer = 0;
    ull reverse_kmer = 0;
    bool in_sequence = false;

    auto save_current_sequence = [&]()
    {
        if (!in_sequence) return;

        vector<uint8> assigned_groups;
        for (size_t group = 1; group < group_names.size(); ++group)
        {
            if (group_counts[group] >= cutoff)
            {
                assigned_groups.push_back(static_cast<uint8>(group));
            }
        }

        sequence_groups[current_name] = assigned_groups;
        sequence_groups_by_order.push_back(assigned_groups);
        std::fill(group_counts.begin(), group_counts.end(), 0);
    };

    while (fastafile.nextLine(StrLine))
    {
        if (StrLine.empty()) continue;

        switch (StrLine[0])
        {
            case '>':
                save_current_sequence();
                current_name = StrLine.substr(1);
                current_size = 0;
                current_kmer = 0;
                reverse_kmer = 0;
                in_sequence = true;
                continue;
            case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }

        if (!in_sequence) continue;

        for (char base : StrLine)
        {
            if (base == '\0') break;
            if (base == '\n' || base == ' ') continue;

            assigner_kmer_read_c(base, current_size, current_kmer, reverse_kmer);

            if (current_size < klen) continue;

            ull larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
            auto find = kmer_hash.find(larger_kmer);

            if (find != kmer_hash.end())
            {
                group_counts[find->second]++;
            }
        }
    }

    save_current_sequence();
    fastafile.Close();
}

template <int dictsize>
void kmer_assigner<dictsize>::assigngroup(string &inputfile, string &kmerfile, string &outputfolder)
{
    kmer_hash.clear();
    group_index.clear();
    sequence_groups.clear();
    sequence_groups_by_order.clear();
    group_names.clear();
    group_names.push_back("");

    fasta targetfile(kmerfile.c_str());
    read_kmertarget(targetfile);

    fasta readsfile(inputfile.c_str());
    read_fasta(readsfile);

    if (ifbed)
    {
        vector<std::string> bed_scores(sequence_groups_by_order.size(), "0");
        for (size_t i = 0; i < sequence_groups_by_order.size(); ++i)
        {
            vector<std::string> one_scores;
            for (uint8 group : sequence_groups_by_order[i])
            {
                if (group < group_names.size())
                {
                    one_scores.push_back(assigner_group_to_bed_score(group_names[group], group));
                }
            }
            
            if (!one_scores.empty())
            {
                std::ostringstream joined;
                for (size_t j = 0; j < one_scores.size(); ++j)
                {
                    if (j) joined << ",";
                    joined << one_scores[j];
                }
                bed_scores[i] = joined.str();
            }
        }
        
        write_bed_assignments(inputfile, outputfolder, bed_scores);
        return;
    }

    std::filesystem::path output_root(outputfolder.empty() ? "." : outputfolder);
    std::filesystem::create_directories(output_root);

    std::vector<std::ofstream> fasta_files(group_names.size());
    std::vector<std::ofstream> kmer_files(group_names.size());

    for (size_t group = 1; group < group_names.size(); ++group)
    {
        std::filesystem::path group_folder = output_root / group_names[group];
        std::filesystem::create_directories(group_folder);

        std::string fasta_name = group_names[group] + ".fasta";
        std::filesystem::path fasta_path = group_folder / fasta_name;
        std::filesystem::path kmer_path = group_folder / (fasta_name + "_kmer.list");

        fasta_files[group].open(fasta_path);
        if (!fasta_files[group].is_open())
        {
            std::cerr << "Failed to open file: " << fasta_path << "\n";
            return;
        }

        kmer_files[group].open(kmer_path);
        if (!kmer_files[group].is_open())
        {
            std::cerr << "Failed to open file: " << kmer_path << "\n";
            return;
        }
    }

    fasta input_fasta(inputfile.c_str());
    std::string StrLine;
    int sequence_index = -1;
    vector<uint8> *active_groups = nullptr;

    while (input_fasta.nextLine(StrLine))
    {
        if (!StrLine.empty() && StrLine[0] == '>')
        {
            sequence_index++;
            if (sequence_index >= 0 && static_cast<size_t>(sequence_index) < sequence_groups_by_order.size())
            {
                active_groups = &sequence_groups_by_order[sequence_index];
            }
            else
            {
                active_groups = nullptr;
            }
        }

        if (active_groups == nullptr) continue;

        for (uint8 group : *active_groups)
        {
            if (group < fasta_files.size() && fasta_files[group].is_open())
            {
                fasta_files[group] << StrLine << "\n";
            }
        }
    }
    input_fasta.Close();

    fasta target_kmers(kmerfile.c_str());
    uint8 current_group = 0;

    while (target_kmers.nextLine(StrLine))
    {
        if (StrLine.empty()) continue;

        switch (StrLine[0])
        {
            case '>':
            {
                std::string partition_name = assigner_trim(StrLine.substr(1));
                auto find = group_index.find(partition_name);
                current_group = (find == group_index.end()) ? 0 : find->second;
                continue;
            }
            case '+': case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }

        if (current_group == 0) continue;

        if (current_group < kmer_files.size() && kmer_files[current_group].is_open())
        {
            kmer_files[current_group] << ">\n" << StrLine << "\n";
        }
    }
    target_kmers.Close();
}

template <int dictsize>
void kmer_assigner<dictsize>::assigngroup_mul(string &inputfile, vector<std::string> &kmerfiles, string &outputfile)
{
    if (!ifbed)
    {
        std::cerr << "ERROR: -K/--kmer-list multi-matrix mode is only supported for BED input.\n";
        std::_Exit(EXIT_FAILURE);
    }
    
    std::unordered_map<std::string, std::string> kmer_by_prefix;
    std::vector<std::string> matrix_order;
    for (const auto &kmerfile : kmerfiles)
    {
        if (kmerfile.empty()) continue;
        
        std::string prefix = assigner_prefix_from_kmer_path(kmerfile);
        if (prefix.empty()) continue;
        
        auto inserted = kmer_by_prefix.emplace(prefix, kmerfile);
        if (!inserted.second && inserted.first->second != kmerfile)
        {
            std::cerr << "ERROR: Multiple kmer files share prefix " << prefix
                      << ": " << inserted.first->second << " and " << kmerfile << "\n";
            std::_Exit(EXIT_FAILURE);
        }
        
        if (inserted.second) matrix_order.push_back(prefix);
    }
    
    if (kmer_by_prefix.empty())
    {
        std::cerr << "ERROR: No kmer files were loaded for -K/--kmer-list mode.\n";
        std::_Exit(EXIT_FAILURE);
    }
    
    std::ifstream bed(inputfile);
    if (!bed)
    {
        std::cerr << "ERROR: Could not open BED input " << inputfile << " for reading.\n";
        std::_Exit(EXIT_FAILURE);
    }
    
    std::vector<std::string> bed_lines;
    std::vector<std::string> bed_scores;
    std::unordered_map<std::string, std::vector<size_t>> bed_rows_by_matrix;
    
    std::string line;
    while (std::getline(bed, line))
    {
        std::vector<std::string> fields;
        if (!parse_bed_record_fields(line, fields)) continue;
        
        size_t bed_index = bed_lines.size();
        bed_lines.push_back(line);
        bed_scores.push_back("0");
        
        std::string matrix = assigner_find_bed_matrix(fields, kmer_by_prefix);
        if (!matrix.empty())
        {
            bed_rows_by_matrix[matrix].push_back(bed_index);
        }
    }
    
    auto run_one_matrix = [&](const std::string &matrix) -> std::vector<std::pair<size_t, std::string>>
    {
        std::vector<std::pair<size_t, std::string>> result;
        auto rows_find = bed_rows_by_matrix.find(matrix);
        auto kmer_find = kmer_by_prefix.find(matrix);
        
        if (rows_find == bed_rows_by_matrix.end() || rows_find->second.empty() || kmer_find == kmer_by_prefix.end())
        {
            return result;
        }
        
        kmer_assigner<dictsize> assigner;
        fasta targetfile(kmer_find->second.c_str());
        assigner.read_kmertarget(targetfile);
        
        bed_subset_fasta subset_fasta(bed_lines, rows_find->second);
        assigner.read_fasta(subset_fasta);
        
        result.reserve(assigner.sequence_groups_by_order.size());
        for (size_t i = 0; i < assigner.sequence_groups_by_order.size() && i < rows_find->second.size(); ++i)
        {
            std::vector<std::string> one_scores;
            for (uint8 group : assigner.sequence_groups_by_order[i])
            {
                if (group < assigner.group_names.size())
                {
                    one_scores.push_back(assigner_group_to_bed_score(assigner.group_names[group], group));
                }
            }
            
            std::string score = "0";
            if (!one_scores.empty())
            {
                std::ostringstream joined;
                for (size_t j = 0; j < one_scores.size(); ++j)
                {
                    if (j) joined << ",";
                    joined << one_scores[j];
                }
                score = joined.str();
            }
            
            result.push_back({rows_find->second[i], score});
        }
        
        return result;
    };
    
    const size_t max_workers = std::max<size_t>(1, static_cast<size_t>(std::max(1, nthreads)));
    std::vector<std::future<std::vector<std::pair<size_t, std::string>>>> active;
    
    auto collect_front = [&]()
    {
        auto matrix_scores = active.front().get();
        active.erase(active.begin());
        
        for (const auto &one : matrix_scores)
        {
            if (one.first < bed_scores.size())
            {
                bed_scores[one.first] = one.second;
            }
        }
    };
    
    for (const auto &matrix : matrix_order)
    {
        if (bed_rows_by_matrix.find(matrix) == bed_rows_by_matrix.end()) continue;
        
        active.push_back(std::async(std::launch::async, run_one_matrix, matrix));
        if (active.size() >= max_workers)
        {
            collect_front();
        }
    }
    
    while (!active.empty())
    {
        collect_front();
    }
    
    write_bed_assignments(inputfile, outputfile, bed_scores);
}

#endif /* KmerAssigner_hpp */
