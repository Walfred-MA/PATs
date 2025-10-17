//
//  KmerCounter.hpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 2/2/23.
//  Copyright © 2023 USC_Mark. All rights reserved.
//

#ifndef KmerCounter_hpp
#define KmerCounter_hpp

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
    
    void getnorm(string &inputfile, string &kmerfile, string &outputfile);
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


template <int dictsize>
template <class typefile>
void kmer_counter<dictsize>::read_target(typefile &fastafile)
{
        
    int current_size = 0;
    
    ull current_kmer = 0;
    ull reverse_kmer = 0;

    std::string StrLine;
    
    vector<vector<uint8_t>> norms;
    
    kmer32_dict* usehash = &kmer_hash;
    ull* usecounter = &totalkmers;
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '>': case '+':
                current_size = 0;
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
                               
            if (! (ifmask && base >= 'a'))
            {
                usehash = &kmer_hash;
                usecounter = &totalkmers;
            }
            else
            {
                usehash = &masked_hash;
                usecounter = &totalmasked;
            }

            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
            
            auto find = usehash->find(larger_kmer);
            
            if (find == usehash->end())
            {
                (*usehash)[larger_kmer] = (*usecounter)++;
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
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '>':
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
    if (old < MAX_UINT16)
        atom.compare_exchange_weak(old, old + 1, std::memory_order_relaxed);
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

template <int dictsize>
void kmer_counter<dictsize>::write(string &inputfile, string &outputfile, vector<uint> & groups, vector<vector<bool>> &norms)
{
    vector<uint> groupsort;
    uint counter = compute_group_order(groups, groupsort);
    
    std::vector<bool> ifopen(counter, 0) ;
    std::vector<string> filenames(counter, "") ;
    std::vector<std::ofstream> kmerfiles(counter);
    kmerfiles[0].open(outputfile + "_filtered.fa_kmer.list");
    ifopen[0] = 1;
    
    for (int i = 1; i < counter; ++i)
    {
        filenames[i] = outputfile + "p" + std::to_string(i)+".fa_kmer.list";
    }
    
    for (const auto& pair : kmer_hash)
    {
        vector<bool> & currnorm = norms[pair.second];
        uint groupindex = 0;
        for (int i = 0 ; i < groups.size(); ++i)
        {
            if (currnorm[i] == 1 )
            {
                if (groupindex > 0 && groupindex != groups[i])
                {
                    groupindex = 0;
                    break;
                }
                else
                {
                    groupindex = groups[i];
                }
                
            }
        }
        
        uint writeindex = groupsort[groupindex];
        
        if (not ifopen[writeindex])
        {
            kmerfiles[writeindex].open(filenames[writeindex]);
            if (!kmerfiles[writeindex].is_open()) {
                std::cerr << "Failed to open file: " << filenames[writeindex] << "\n";
                continue;
            }
            ifopen[writeindex] = 1;
        }

        kmerfiles[writeindex] << ">\n"+kmer_int_toatcg(pair.first)+"\n";
    }
    
    
    for (const auto& pair : masked_hash)
    {
        vector<bool> & currnorm = masknorms[pair.second];
        uint groupindex = 0;
        for (int i = 0 ; i < groups.size(); ++i)
        {
            if (currnorm[i] == 1 )
            {
                if (groupindex > 0 && groupindex != groups[i])
                {
                    groupindex = 0;
                    break;
                }
                else
                {
                    groupindex = groups[i];
                }
                
            }
        }
        
        uint writeindex = groupsort[groupindex];
        
        if (not ifopen[writeindex])
        {
            kmerfiles[writeindex].open(filenames[writeindex]);
            if (!kmerfiles[writeindex].is_open()) {
                std::cerr << "Failed to open file: " << filenames[writeindex] << "\n";
                continue;
            }
            ifopen[writeindex] = 1;
        }

        kmerfiles[writeindex] << ">\n"+kmer_int_toatcg_l(pair.first)+"\n";
    }
    kmerfiles.clear();
    
    std::vector<std::ofstream> seqfiles(counter);

    seqfiles[0].open(outputfile + "_filtered.fa");
    for (int i = 1; i < counter; ++i)
    {
        if (ifopen[i] == 0) continue;
        
        auto filename = outputfile + "p" + std::to_string(i)+".fa";
        seqfiles[i].open(filename);

        if (!seqfiles[i].is_open()) {
            std::cerr << "Failed to open file: " << filename << "\n";
            continue;
        }
    }
    
    fasta fastafile(inputfile.c_str());
    std::string StrLine;
    int sample_index = -1;
    int sample_group = 0;
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '>':
                sample_index ++;
                sample_group = groupsort[groups[sample_index]];
                if (not ifopen[sample_group]) sample_group = 0;
                break;
            default:
                break;
        }
        seqfiles[sample_group] << StrLine + "\n";
    }
    seqfiles.clear();
    fastafile.Close();
    
    
    
    

    return ;
}


#endif /* KmerCounter_hpp */

