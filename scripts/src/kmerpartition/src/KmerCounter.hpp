//
//  KmerCounter.hpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 2/2/23.
//  Copyright © 2023 USC_Mark. All rights reserved.
//

#ifndef KmerCounter_hpp
#define KmerCounter_hpp

#define MAX_GROUP_FILES 99

#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstring>
#include <unordered_set>
#include <algorithm>
#include <thread>
#include <atomic>
#include <mutex>
#include <sstream>
#include <math.h>
#include <numeric>
#include <cstdint>
#include <omp.h>
#include "fasta.hpp"
#include "KmerStruct.hpp"

using namespace std;

extern bool ifweight;
extern bool ifmask;
extern bool ifhighres;
extern int cutoff;
extern int nthreads;
extern float similar;

using kmer64_dict = std::unordered_map<u128, uint, hash_128> ;
using kmer32_dict = std::unordered_map<ull, uint > ;

using kmer64_dict_nt = std::unordered_map<u128, uint8, hash_128> ;
using kmer32_dict_nt = std::unordered_map<ull, uint8>;

using kmer64_dict_mul = std::unordered_map<u128, uint*, hash_128> ;
using kmer32_dict_mul = std::unordered_map<ull, uint* > ;

#define sington_weight 0.05
#define sington_weight_r 20
#define MAX_UINT16 65535
template <int dictsize>
class kmer_counter
{
    kmer32_dict kmer_hash, masked_hash;
    
    vector<vector<bool>> norms, masknorms;
    vector<std::string> sample_names;
    vector<std::string> kmer_texts, masked_kmer_texts;
    
    ull totalkmers = 1, totalmasked = 1;
    ull totalsamples = 0;
public:
    
    kmer_counter()
    {};
    ~kmer_counter()
    {};

    template <class typefile>
    void read_target(typefile &fastafile);
    
    template <class typefile>
    void read_file(typefile &fastafile);
    
    void calculatenorm(vector<std::atomic<uint16_t>> & thenorm);
    void determinegroup(vector<uint> &groups, vector<std::atomic<uint16_t>> & thenorm);
    
    void write(string &inputfile, string &outputfile, vector<uint> & groups, vector<vector<bool>> &norms);
    void write_target_partition(string &inputfile, string &outputfile, const std::unordered_set<std::string> &target_names);
    
    void getnorm(string &inputfile, string &kmerfile, string &outputfile);
    void gettargetpartition(string &inputfile, string &kmerfile, string &outputfile, const std::unordered_set<std::string> &target_names);
};

inline static void update_counter(kmer32_dict &target_map, ull &larger_kmer, uint8* samplevecs)
{
    auto map_find = target_map.find(larger_kmer);

    if (map_find != target_map.end())
    {
        if (samplevecs[map_find->second] < 255) samplevecs[map_find->second] ++;
    }
}


template <typename T>
static void kmer_read_c(char base, int &current_size, T &current_kmer, T &reverse_kmer)
{
    int converted = 0;
    T reverse_converted;
    
    if (base == '\n' || base == ' ') return;
    
    if (base_to_int(base, converted))
    {
        
        current_kmer <<= ( 8*sizeof(current_kmer) - 2*klen + 2 );
        current_kmer >>=  ( 8*sizeof(current_kmer) - 2*klen );
        current_kmer += converted;
        
        reverse_kmer >>= 2;
        reverse_converted = 0b11-converted;
        reverse_converted <<= (2*klen-2);
        reverse_kmer += reverse_converted;
        
        current_size ++;
        
    }
    
    else
    {
        current_size = 0;
        current_kmer = 0;
        reverse_kmer = 0;
    }
}

inline static std::string counter_sequence_name_from_header(const std::string &header)
{
    std::string body = header;
    if (!body.empty() && body[0] == '>') body = body.substr(1);
    
    if (ifbed)
    {
        std::vector<std::string> fields = split_text_fields(body);
        if (fields.size() >= 4) return fields[3];
        if (!fields.empty()) return fields[0];
        return "";
    }
    
    size_t tab = body.find('\t');
    if (tab != std::string::npos) body = body.substr(0, tab);
    
    std::vector<std::string> fields = split_text_fields(body);
    if (!fields.empty()) return fields[0];
    
    return body;
}


template <int dictsize>
template <class typefile>
void kmer_counter<dictsize>::read_target(typefile &fastafile)
{
        
    int current_size = 0;
    
    ull current_kmer = 0;
    ull reverse_kmer = 0;
    std::string current_text_kmer;

    std::string StrLine;
    
    vector<vector<uint8_t>> norms;
    
    kmer32_dict* usehash = &kmer_hash;
    ull* usecounter = &totalkmers;
    vector<std::string>* usetexts = &kmer_texts;
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '>': case '+':
                current_size = 0;
                current_text_kmer.clear();
                continue;
            case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }
        
        for (auto base: StrLine)
        {
            if (base == '\0') break;
            
            if (base == '\n' || base == ' ') continue;

            int converted = 0;
            if (base_to_int(base, converted))
            {
                current_text_kmer.push_back(base);
                if (current_text_kmer.size() > klen) current_text_kmer.erase(0, current_text_kmer.size() - klen);
            }
            else
            {
                current_text_kmer.clear();
            }

            kmer_read_c(base, current_size, current_kmer, reverse_kmer);
           
            if (current_size < klen ) continue;
                               
            if ( ! ( ifmask && base >= 'a'))
            {
                usehash = &kmer_hash;
                usecounter = &totalkmers;
                usetexts = &kmer_texts;
            }
            else
            {
                usehash = &masked_hash;
                usecounter = &totalmasked;
                usetexts = &masked_kmer_texts;
            }

            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
            
            auto find = usehash->find(larger_kmer);
            
            if (find == usehash->end())
            {
                uint index = static_cast<uint>((*usecounter)++);
                (*usehash)[larger_kmer] = index;
                if (usetexts->size() <= index) usetexts->resize(index + 1);
                (*usetexts)[index] = current_text_kmer;
            }
            
        }
    }
    
    fastafile.Close();

};


template <int dictsize>
template <class typefile>
void kmer_counter<dictsize>::read_file(typefile &fastafile)
{
        
    int current_size = 0;
    
    ull current_kmer = 0;
    ull reverse_kmer = 0;
    std::string StrLine;
    totalsamples = 0;
    sample_names.clear();
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '>':
                sample_names.push_back(counter_sequence_name_from_header(StrLine));
                totalsamples++;
                continue;
            case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }
        
    }
    fastafile.Reset();
    
    norms.resize(totalkmers);
    for (size_t i =0 ; i < norms.size(); ++i)
    {
        norms[i].resize(totalsamples);
    }
    
    masknorms.resize(totalmasked);
    for (size_t i =0 ; i < masknorms.size(); ++i)
    {
        masknorms[i].resize(totalsamples);
    }
    
    int sample_index = -1;
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '>':
                current_size = 0;
                sample_index ++;
                continue;
            case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }
        
        
        for (auto base: StrLine)
        {
            if (base == '\0') break;
            
            if (base == '\n' || base == ' ') continue;

            kmer_read_c(base, current_size, current_kmer, reverse_kmer);
           
            if (current_size < klen ) continue;
                                
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
            
            auto find = kmer_hash.find(larger_kmer);
            
            if (find != kmer_hash.end())
            {
                norms[find->second][sample_index] = 1;
            }
            
            if (totalmasked > 1)
            {
                auto find_m = masked_hash.find(larger_kmer);
                
                if (find_m != masked_hash.end())
                {
                    masknorms[find_m->second][sample_index] = 1;
                }
            }
            
        }
    }
    
    fastafile.Close();
    
};

inline static void capped_atomic_increment(std::atomic<uint16_t>& atom) {
    uint16_t old = atom.load(std::memory_order_relaxed);
    while (old < MAX_UINT16 &&
           !atom.compare_exchange_weak(old, old + 1, std::memory_order_relaxed))
    {
        // compare_exchange_weak refreshes `old` with the current value on
        // failure; retry until the increment lands or the cap is reached.
    }
}

inline static void offsite_operation(std::atomic<ull> &offsite, const vector<bool> norm, vector<std::atomic<int>>& vec_offsite, vector<std::atomic<uint16_t>>&normmatrix, size_t colsize)
{
    vector<int> indexes(colsize, 0);
    size_t total = std::accumulate(norm.begin(), norm.end(), 0u);
    size_t index_ct = 0;
    if (2*total > norm.size())
    {
        offsite.fetch_add(1, std::memory_order_relaxed);
        
        for (size_t i = 0; i < norm.size(); ++i)
        {
            if (norm[i] == 0)
            {
                vec_offsite[i].fetch_sub(1, std::memory_order_relaxed);
                indexes[index_ct ++] = i;
            }
        }
        
        for (size_t i = 0; i < index_ct; ++i)
        {
            size_t index_i = indexes[i];
            size_t index_corr = ((index_i + 1)*index_i)/2;
            for (size_t j = i; j < index_ct; ++j)
            {
                auto& value = normmatrix[index_i*colsize + indexes[j]  - index_corr];
                capped_atomic_increment(value);
            }
        }
        
    }
    
    else
    {
        for (size_t i = 0; i < norm.size(); ++i)
        {
            
            if (norm[i] > 0)
            {
                indexes[index_ct ++] = i;
            }
        }
        
        for (size_t i = 0; i < index_ct; ++i)
        {
            size_t index_i = indexes[i];
            size_t index_corr = ((index_i + 1)*index_i)/2;
            for (size_t j = i; j < index_ct; ++j)
            {
                auto& value = normmatrix[index_i*colsize + indexes[j]  - index_corr];
                capped_atomic_increment(value);
            }
        }
    }
}



template <int dictsize>
void kmer_counter<dictsize>::calculatenorm(vector<std::atomic<uint16_t>>& normmatrix)
{
    size_t colsize = totalsamples;
    size_t rowsize = totalkmers;
            
    std::atomic<ull> offsite = 0;
    vector<std::atomic<int>> vec_offsite(colsize);
    for (size_t i = 0; i < colsize; ++i)
        vec_offsite[i].store(0, std::memory_order_relaxed);
    

    omp_set_num_threads(nthreads);
    #pragma omp parallel for schedule(dynamic)
    for (size_t k = 0; k < norms.size(); ++k)
    {
        const auto& norm = norms[k];

        if (*std::max_element(norm.begin(), norm.end()) == 0) continue;
        if (std::accumulate(norm.begin(), norm.end(), 0u) <= 1) continue;

        offsite_operation(offsite, norm, vec_offsite, normmatrix, colsize);
    }

    size_t lastindex = 0;
    for (size_t i =0; i< colsize ; ++i)
    {
        size_t rowend =  lastindex + colsize - i;
        
        for (;lastindex < rowend; ++lastindex)
        {
            normmatrix[lastindex] += offsite + vec_offsite[i] + vec_offsite[colsize - (rowend - lastindex)] ;
        }
    }
}

template <int dictsize>
void kmer_counter<dictsize>::determinegroup(vector<uint> &group_edges,  vector<std::atomic<uint16_t>>& matrix)
{
    size_t size = group_edges.size();
    // Compute scaled diagonal values
    vector<float> diags(size);
    for (size_t i = 0; i < size; ++i)
    {
        size_t diag_index = i * size + i - (i * (i + 1)) / 2;
        diags[i] = matrix[diag_index] * similar;
    }

    for (size_t i = 0; i < size; ++i)
    {
        group_edges[i] = i;
    }

    // Build group edges
    for (size_t i = 0; i < size; ++i)
    {
        size_t indexstart = i * size - (i * (i + 1)) / 2;
        float diag_cutoff = diags[i];

        for (size_t j = i + 1; j < size; ++j)
        {
            uint16_t edge = matrix[indexstart + j];
            float diag_cutoff2 = diags[j];

            
            if (edge > std::max((float)cutoff, std::min(diag_cutoff, diag_cutoff2)))
            {
                if (group_edges[i] < group_edges[j])
                {
                    group_edges[j] = group_edges[i];
                } else
                {
                    group_edges[i] = group_edges[j];
                }
            }
        }
    }

    // Finalize group assignments
    for (size_t i = 0; i < size; ++i)
    {
        uint group_index = group_edges[i];
        uint last_group_index = i;

        while (group_index < last_group_index) {
            last_group_index = group_index;
            group_index = group_edges[group_index];
        }

        group_edges[i] = group_index;
    }
    
    for (size_t i = 0; i < size; ++i) group_edges[i]++;
    
}

template <int dictsize>
void kmer_counter<dictsize>::getnorm(string &inputfile, string &kmerfile, string &outputfile)
{
    fasta targetfile(kmerfile.c_str());
    
    read_target(targetfile);
    
    fasta readsfile(inputfile.c_str());
    
    read_file(readsfile);
    
    vector<std::atomic<uint16_t>> normmatrix((totalsamples*(totalsamples+1))/2);
    
    calculatenorm(normmatrix);
    
    vector<uint> groups(totalsamples,0);
    
    determinegroup(groups, normmatrix);
    
    write(inputfile, outputfile, groups, norms);
    
}

template <int dictsize>
void kmer_counter<dictsize>::gettargetpartition(string &inputfile, string &kmerfile, string &outputfile, const std::unordered_set<std::string> &target_names)
{
    fasta targetfile(kmerfile.c_str());
    
    read_target(targetfile);
    
    fasta readsfile(inputfile.c_str());
    
    read_file(readsfile);
    
    write_target_partition(inputfile, outputfile, target_names);
}

/*
inline static uint compute_group_order(const vector<uint>& groups, vector<uint>& groupsort)
{
    size_t N = groups.size();
    vector<uint> counts(N + 1, 0);

    // Step 1: Count occurrences
    for (uint g : groups)
    {
        counts[g]++;
    }

    // Step 2: Gather group IDs with nonzero count
    vector<uint> unique_gids;
    for (uint g = 0; g <= N; ++g)
    {
        if (counts[g] > 0)
            unique_gids.push_back(g);
    }

    // Step 3: Sort group IDs by count descending
    sort(unique_gids.begin(), unique_gids.end(), [&](uint a, uint b) {
        return counts[a] > counts[b];
    });

    // Step 4: Fill in group ranks
    groupsort.assign(N + 1, 0); // Initialize with 0 (unused groups)
    for (size_t i = 0; i < unique_gids.size(); ++i)
    {
        groupsort[unique_gids[i]] = i + 1; // rank = 1-based
    }
    
    return unique_gids.size() + 1;
}
*/

#define MAX_OPEN_GROUP_FILES_PER_PASS 1000

inline static uint compute_group_order(const vector<uint>& groups, vector<uint>& groupsort)
{
    size_t N = groups.size();
    vector<uint> counts(N + 1, 0);

    for (uint g : groups)
    {
        counts[g]++;
    }

    vector<uint> unique_gids;
    unique_gids.reserve(N + 1);
    for (uint g = 0; g <= N; ++g)
    {
        if (counts[g] > 0)
            unique_gids.push_back(g);
    }

    sort(unique_gids.begin(), unique_gids.end(), [&](uint a, uint b) {
        return counts[a] > counts[b];
    });

    groupsort.assign(N + 1, 0); // 0 = filtered
    for (size_t i = 0; i < unique_gids.size(); ++i)
    {
        groupsort[unique_gids[i]] = static_cast<uint>(i + 1);
    }

    return static_cast<uint>(unique_gids.size() + 1); // +1 for filtered bucket 0
}

inline static void update_target_partition_counts(
    const vector<vector<bool>> &source_norms,
    const vector<bool> &is_target_sample,
    vector<bool> &keep_kmer,
    vector<int> &sample_counts)
{
    keep_kmer.assign(source_norms.size(), false);
    
    for (size_t kmer_index = 0; kmer_index < source_norms.size(); ++kmer_index)
    {
        const vector<bool> &currnorm = source_norms[kmer_index];
        bool in_target = false;
        
        for (size_t sample_index = 0; sample_index < currnorm.size() && sample_index < is_target_sample.size(); ++sample_index)
        {
            if (is_target_sample[sample_index] && currnorm[sample_index])
            {
                in_target = true;
                break;
            }
        }
        
        if (!in_target) continue;
        
        keep_kmer[kmer_index] = true;
        for (size_t sample_index = 0; sample_index < currnorm.size() && sample_index < sample_counts.size(); ++sample_index)
        {
            if (currnorm[sample_index]) sample_counts[sample_index]++;
        }
    }
}

inline static void prune_target_partition_kmers(
    const vector<vector<bool>> &source_norms,
    const vector<bool> &selected_samples,
    vector<bool> &keep_kmer)
{
    for (size_t kmer_index = 0; kmer_index < source_norms.size() && kmer_index < keep_kmer.size(); ++kmer_index)
    {
        if (!keep_kmer[kmer_index]) continue;
        
        const vector<bool> &currnorm = source_norms[kmer_index];
        bool in_selected = false;
        
        for (size_t sample_index = 0; sample_index < currnorm.size() && sample_index < selected_samples.size(); ++sample_index)
        {
            if (selected_samples[sample_index] && currnorm[sample_index])
            {
                in_selected = true;
                break;
            }
        }
        
        keep_kmer[kmer_index] = in_selected;
    }
}

template <int dictsize>
void kmer_counter<dictsize>::write_target_partition(string &inputfile, string &outputfile, const std::unordered_set<std::string> &target_names)
{
    vector<bool> is_target_sample(totalsamples, false);
    size_t matched_targets = 0;
    
    for (size_t i = 0; i < sample_names.size() && i < is_target_sample.size(); ++i)
    {
        if (target_names.find(sample_names[i]) != target_names.end())
        {
            is_target_sample[i] = true;
            matched_targets++;
        }
    }
    
    if (matched_targets == 0)
    {
        std::cerr << "WARNING: No -t/--targets names matched input sequence names.\n";
    }
    
    vector<int> sample_counts(totalsamples, 0);
    vector<bool> keep_kmer, keep_masked_kmer;
    
    update_target_partition_counts(norms, is_target_sample, keep_kmer, sample_counts);
    update_target_partition_counts(masknorms, is_target_sample, keep_masked_kmer, sample_counts);
    
    vector<bool> selected_samples(totalsamples, false);
    for (size_t i = 0; i < sample_counts.size(); ++i)
    {
        selected_samples[i] = sample_counts[i] >= cutoff;
    }
    
    prune_target_partition_kmers(norms, selected_samples, keep_kmer);
    prune_target_partition_kmers(masknorms, selected_samples, keep_masked_kmer);
    
    std::string kmer_output = outputfile + "_kmer.list";
    std::ofstream kmerfile(kmer_output);
    if (!kmerfile.is_open())
    {
        std::cerr << "Failed to open file: " << kmer_output << "\n";
        return;
    }
    
    for (const auto &pair : kmer_hash)
    {
        if (pair.second < keep_kmer.size() && keep_kmer[pair.second])
        {
            if (pair.second < kmer_texts.size() && !kmer_texts[pair.second].empty())
            {
                kmerfile << ">\n" << kmer_texts[pair.second] << "\n";
            }
            else
            {
                kmerfile << ">\n" << kmer_int_toatcg(pair.first) << "\n";
            }
        }
    }
    
    for (const auto &pair : masked_hash)
    {
        if (pair.second < keep_masked_kmer.size() && keep_masked_kmer[pair.second])
        {
            if (pair.second < masked_kmer_texts.size() && !masked_kmer_texts[pair.second].empty())
            {
                kmerfile << ">\n" << masked_kmer_texts[pair.second] << "\n";
            }
            else
            {
                kmerfile << ">\n" << kmer_int_toatcg_l(pair.first) << "\n";
            }
        }
    }
    
    if (ifbed)
    {
        write_bed_selected_assignments(inputfile, outputfile, selected_samples, "1");
        return;
    }
    
    std::string fasta_output = outputfile;
    std::ofstream seqfile(fasta_output);
    if (!seqfile.is_open())
    {
        std::cerr << "Failed to open file: " << fasta_output << "\n";
        return;
    }
    
    fasta fastafile(inputfile.c_str());
    std::string StrLine;
    int sample_index = -1;
    bool write_this_sample = false;
    
    while (fastafile.nextLine(StrLine))
    {
        if (!StrLine.empty() && StrLine[0] == '>')
        {
            sample_index++;
            write_this_sample = sample_index >= 0 &&
                                static_cast<size_t>(sample_index) < selected_samples.size() &&
                                selected_samples[sample_index];
        }
        
        if (write_this_sample)
        {
            seqfile << StrLine << "\n";
        }
    }
    
    fastafile.Close();
}

template <int dictsize>
void kmer_counter<dictsize>::write(string &inputfile, string &outputfile, vector<uint> & groups, vector<vector<bool>> &norms)
{
    vector<uint> groupsort;
    uint counter = compute_group_order(groups, groupsort);   // 0..counter-1

    // Precompute sample -> output bucket
    vector<uint> sample_writeindex(groups.size(), 0);
    for (size_t i = 0; i < groups.size(); ++i)
    {
        sample_writeindex[i] = groupsort[groups[i]];
    }

    if (ifbed)
    {
        vector<std::string> bed_scores(sample_writeindex.size(), "0");
        for (size_t i = 0; i < sample_writeindex.size(); ++i)
        {
            bed_scores[i] = std::to_string(sample_writeindex[i]);
        }
        write_bed_assignments(inputfile, outputfile, bed_scores);
        return;
    }

    auto resolve_writeindex = [&](const vector<bool>& currnorm) -> uint
    {
        uint groupindex = 0;
        for (size_t i = 0; i < groups.size(); ++i)
        {
            if (currnorm[i])
            {
                if (groupindex > 0 && groupindex != groups[i])
                {
                    return 0; // mixed groups => filtered
                }
                groupindex = groups[i];
            }
        }
        return groupsort[groupindex];
    };

    uint total_group_outputs = (counter > 1) ? (counter - 1) : 0;
    uint num_batches = std::max<uint>(
        1,
        (total_group_outputs + MAX_OPEN_GROUP_FILES_PER_PASS - 1) / MAX_OPEN_GROUP_FILES_PER_PASS
    );

    // -----------------------------
    // Pass 1: write *.fa_kmer.list
    // -----------------------------
    for (uint batch_id = 0; batch_id < num_batches; ++batch_id)
    {
        uint batch_begin = 1 + batch_id * MAX_OPEN_GROUP_FILES_PER_PASS;
        uint batch_end   = std::min(counter, batch_begin + MAX_OPEN_GROUP_FILES_PER_PASS);
        bool include_filtered = (batch_id == 0);

        std::ofstream filtered_kmerfile;
        if (include_filtered)
        {
            filtered_kmerfile.open(outputfile + "_filtered.fa_kmer.list");
            if (!filtered_kmerfile.is_open())
            {
                std::cerr << "Failed to open file: " << outputfile + "_filtered.fa_kmer.list" << "\n";
                return;
            }
        }

        std::vector<std::ofstream> kmerfiles(batch_end > batch_begin ? batch_end - batch_begin : 0);
        for (uint writeindex = batch_begin; writeindex < batch_end; ++writeindex)
        {
            std::string fname = outputfile + "p" + std::to_string(writeindex) + ".fa_kmer.list";
            kmerfiles[writeindex - batch_begin].open(fname);
            if (!kmerfiles[writeindex - batch_begin].is_open())
            {
                std::cerr << "Failed to open file: " << fname << "\n";
            }
        }

        for (const auto& pair : kmer_hash)
        {
            const vector<bool>& currnorm = norms[pair.second];
            uint writeindex = resolve_writeindex(currnorm);

            if (writeindex == 0)
            {
                if (include_filtered)
                    filtered_kmerfile << ">\n" << kmer_int_toatcg(pair.first) << "\n";
            }
            else if (writeindex >= batch_begin && writeindex < batch_end)
            {
                auto& out = kmerfiles[writeindex - batch_begin];
                if (out.is_open())
                    out << ">\n" << kmer_int_toatcg(pair.first) << "\n";
            }
        }

        for (const auto& pair : masked_hash)
        {
            const vector<bool>& currnorm = masknorms[pair.second];
            uint writeindex = resolve_writeindex(currnorm);

            if (writeindex == 0)
            {
                if (include_filtered)
                    filtered_kmerfile << ">\n" << kmer_int_toatcg_l(pair.first) << "\n";
            }
            else if (writeindex >= batch_begin && writeindex < batch_end)
            {
                auto& out = kmerfiles[writeindex - batch_begin];
                if (out.is_open())
                    out << ">\n" << kmer_int_toatcg_l(pair.first) << "\n";
            }
        }
    }

    // -----------------------------
    // Pass 2: write *.fa
    // -----------------------------
    for (uint batch_id = 0; batch_id < num_batches; ++batch_id)
    {
        uint batch_begin = 1 + batch_id * MAX_OPEN_GROUP_FILES_PER_PASS;
        uint batch_end   = std::min(counter, batch_begin + MAX_OPEN_GROUP_FILES_PER_PASS);
        bool include_filtered = (batch_id == 0);

        std::ofstream filtered_seqfile;
        if (include_filtered)
        {
            filtered_seqfile.open(outputfile + "_filtered.fa");
            if (!filtered_seqfile.is_open())
            {
                std::cerr << "Failed to open file: " << outputfile + "_filtered.fa" << "\n";
                return;
            }
        }

        std::vector<std::ofstream> seqfiles(batch_end > batch_begin ? batch_end - batch_begin : 0);
        for (uint writeindex = batch_begin; writeindex < batch_end; ++writeindex)
        {
            std::string fname = outputfile + "p" + std::to_string(writeindex) + ".fa";
            seqfiles[writeindex - batch_begin].open(fname);
            if (!seqfiles[writeindex - batch_begin].is_open())
            {
                std::cerr << "Failed to open file: " << fname << "\n";
            }
        }

        fasta fastafile(inputfile.c_str());
        std::string StrLine;
        int sample_index = -1;
        int active_writeindex = -1;
        bool write_this_sample = false;

        while (fastafile.nextLine(StrLine))
        {
            if (!StrLine.empty() && StrLine[0] == '>')
            {
                sample_index++;
                active_writeindex = static_cast<int>(sample_writeindex[sample_index]);

                if (active_writeindex == 0)
                    write_this_sample = include_filtered;
                else
                    write_this_sample = (active_writeindex >= (int)batch_begin &&
                                         active_writeindex <  (int)batch_end);
            }

            if (!write_this_sample) continue;

            if (active_writeindex == 0)
            {
                filtered_seqfile << StrLine << "\n";
            }
            else
            {
                auto& out = seqfiles[active_writeindex - batch_begin];
                if (out.is_open())
                    out << StrLine << "\n";
            }
        }

        fastafile.Close();
    }
}


#endif /* KmerCounter_hpp */
