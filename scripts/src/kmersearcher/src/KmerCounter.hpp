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

#include "fasta.hpp"
#include "fastq.hpp"
#include "KmerStruct.hpp"
#include "CramReader.hpp"
#include "config.hpp"
#include "KmerHash.hpp"

using namespace std;

extern bool ifmask;

using kmer64_dict = std::unordered_map<u128, uint, hash_128> ;
using kmer32_dict = std::unordered_map<ull, uint > ;

using kmer64_dict_nt = std::unordered_map<u128, countint, hash_128> ;
using kmer32_dict_nt = std::unordered_map<ull, countint>;

using kmer64_dict_mul = std::unordered_map<u128, uint*, hash_128> ;
using kmer32_dict_mul = std::unordered_map<ull, uint* > ;

template <int dictsize>
class kmer_counter
{
    using kmer_int = typename std::conditional<(dictsize>32), u128, ull>::type;
    using kmer_dict_type = typename std::conditional<(dictsize>32), kmer64_dict, kmer32_dict>::type;
    using kmer_dict_type_nt = typename std::conditional<(dictsize>32), kmer64_dict_nt, kmer32_dict_nt>::type;
    using kmer_dict_type_mul = typename std::conditional<(dictsize>32), kmer64_dict_mul, kmer32_dict_mul>::type;
    
    kmer_dict_type target_map;
    kmer_dict_type_nt target_map_nt;
    kmer_dict_type_mul target_map_mul;
    
    
    ull totalkmers = 1;
    kmer_int *kmer_records;
    
    int klen = 31 , knum = 0;
    bool iftarget = 0;
    vector<pair<ull,ull>> targetranges;
    uint targetindex = 0;
    std::atomic_uint restfileindex ;
    std::mutex Threads_lock;
    
    std::vector<std::thread*> threads;
    std::vector<std::string> inputfiles;
    std::vector<std::string> outputfiles;
    std::vector<std::string> prefixes;

public:
    
    kmer_counter (int kmersize): klen(31)
    {
        kmer_records = (kmer_int *) malloc(1000);
    };
    ~kmer_counter()
    {
        free(kmer_records);
    };
    template <class typefile>
    void read_counttarget(typefile &fastafile);
    
    void read_hashtarget(fasta &fastafile);
    
    void read_target(const char* infile);
        
    template <class typefile, class typeint>
    void count_kmer(typefile &fastafile, typeint* samplevecs);

    void read_files(std::vector<std::string>& inputfiles, std::vector<std::string>& outputfiles, std::vector<std::string>& prefixes ,int numthread);
    
    void read_file();
            
    template <class typeint>
    void write(const char* outfile, typeint*, pair<ull,ull> range);

};


inline static std::string get_basename(const std::string& filepath)
{
    size_t pos = filepath.find_last_of("/\\");
    if (pos == std::string::npos)
        return filepath;  // No slash found, it's already a filename
    return filepath.substr(pos + 1);
}

//storage additional positions
template <typename T1, typename T2>
inline static void initiate_counter_mul(T2 &target_map, T1 &larger_kmer, uint index)
{
    
    typename T2::iterator map_find = target_map.find(larger_kmer);

    if (map_find == target_map.end() )
    {
        target_map[larger_kmer] = (uint*) malloc(sizeof(uint)*2);
        target_map[larger_kmer][0] = 1;
        target_map[larger_kmer][1] = index ;
    }
    
    else
    {
        uint* &data = map_find->second;
        if (data[0]%5 == 1)
        {
            data = (uint*) realloc(data, sizeof(uint)*(data[0]+1 + 5));
        }
            
        data[0]++;
        data[data[0]] = index ;
    }
}

template <typename T1, typename T2, class typeint>
inline static void update_counter_mul(T2 &target_map_mul, T1 &larger_kmer, typeint* vec, const uint val = 1)
{
    
    typename T2::iterator map_find_mul = target_map_mul.find(larger_kmer);
    
    if (map_find_mul != target_map_mul.end())
    {
        uint* index_mul = map_find_mul->second;
        uint num_num = index_mul[0];
        
        for (uint i = 1 ; i < num_num + 1; ++i)
        {
            if (vec[index_mul[i]] + val > MAXCOUNT)
            {
                vec[index_mul[i]] = MAXCOUNT;
            } else
            {
                vec[index_mul[i]] += val;
            }
        }
    }
}


template <typename T1, typename T2, typename T3>
inline static void initiate_counter(T2 &target_map, T3 &target_map_mul, T1 &larger_kmer, ull &kindex, ull &start, T1 * &kmer_records)
{
    
    typename T2::iterator map_find = target_map.find(larger_kmer);

    if (map_find == target_map.end())
    {
        if ((kindex % 1000000) == 1)
        {
            kmer_records = (T1 *) realloc(kmer_records, sizeof(T1) * (kindex+1000000));
        }
        kmer_records[kindex] = larger_kmer;
        target_map[larger_kmer] = (int) kindex++;
        
        
    }
    
    else if (map_find->second < start)
    {
        initiate_counter_mul(target_map_mul, larger_kmer, map_find->second);
        
        if ((kindex % 1000000) == 1)
        {
            kmer_records = (T1 *) realloc(kmer_records, sizeof(T1) * (kindex+1000000));
        }
        kmer_records[kindex] = larger_kmer;
        target_map[larger_kmer] = (int) kindex++;
        
    }
    
}


template <typename T1, typename T2, typename T3, class typeint>
inline static void update_counter(T2 &target_map,T3 &target_map_mul, T1 &larger_kmer, typeint* vec, const uint val = 1)
{
    
    typename T2::iterator map_find = target_map.find(larger_kmer);

    if (map_find != target_map.end())
    {
        uint index =map_find->second;
        
        if ( index == 0) return;
        
        if (vec[index] + val > MAXCOUNT)
        {
            vec[index] = MAXCOUNT;
        } else
        {
            vec[index] += val;
        }
        
        update_counter_mul(target_map_mul, larger_kmer, vec);
    }
}


template <typename T>
static void kmer_read_c(char base, int klen, int &current_size, T &current_kmer, T &reverse_kmer)
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
        
    }
    
    else
    {
        current_size = 0;
        current_kmer = 0;
        reverse_kmer = 0;
    }
}

template <typename T>
static void kmer_read_f_64(char base, int &current_size, T &current_kmer, T &reverse_kmer)
{
    #define nbase 3
    int converted = 0;
    T reverse_converted;
    
    if (base == '\n' || base == ' ') return;
    
    if (base_to_int_64(base, converted))
    {
        current_size += nbase;
        current_kmer <<= 2*nbase;
        current_kmer += converted;
        
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
void kmer_counter<dictsize>::read_counttarget(typefile &fastafile)
{
        
    int current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    iftarget = 1;
    targetranges.emplace_back();
    auto & targetrange = targetranges[targetranges.size()-1];
    targetrange.first = totalkmers;
    
    std::string StrLine;
    
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '@':  case '+': case '>':
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

            kmer_read_c(base, klen, ++current_size, current_kmer, reverse_kmer);
           
	        if (current_size < klen || (ifmask && base >= 'a') ) continue;
                                
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
            
            initiate_counter(target_map, target_map_mul, larger_kmer, totalkmers, targetrange.first, kmer_records);
                
        }
    }
    
    fastafile.Close();
    
    targetrange.second = totalkmers ;
    
};
 
template <int dictsize>
void kmer_counter<dictsize>::read_hashtarget(fasta &fastafile)
{
    int current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    iftarget = 1;
    targetranges.emplace_back();
    auto & targetrange = targetranges[targetranges.size()-1];
    targetrange.first = totalkmers;
    
    std::string StrLine (10000, '\0');
    std::string Header;
    
    while (fastafile.nextLine(StrLine))
    {
        
        for (auto base: StrLine)
        {
            if (base == '\0') break;
                        
            if (base == '\n' || base == ' ' || base == '\t') break;
 
            kmer_read_f_64(base, current_size, current_kmer, reverse_kmer);
            
            if (current_size < klen) continue;
            //auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
            auto larger_kmer  = current_kmer;
            
            initiate_counter(target_map, target_map_mul, larger_kmer, totalkmers, targetrange.first, kmer_records);
        }
    }
    
    targetrange.second = totalkmers ;
            
}


template <int dictsize>
void kmer_counter<dictsize>::read_target(const char* targetfile)
{
    
    if (
        (std::strlen(targetfile) >= 3 && std::strcmp(targetfile + std::strlen(targetfile) - 3, ".fa") == 0 )
        ||
        (std::strlen(targetfile) >= 6 && std::strcmp(targetfile + std::strlen(targetfile) - 6, ".fasta") == 0)
        ||
        (std::strlen(targetfile) >= 5 && std::strcmp(targetfile + std::strlen(targetfile) - 5, ".list") == 0)
        )
    {
        
        fasta readsfile(targetfile);
        read_counttarget(readsfile);
        
    }
    else
    {
        fasta fastafile(targetfile);
        read_hashtarget(fastafile);
        
    }
    
}



template <int dictsize>
template <class typefile, class typeint>
void kmer_counter<dictsize>::count_kmer(typefile &fastafile, typeint* samplevecs)
{
    
    int current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    //uint64_t ifmasked = 0;
    //int num_masked = 0 ;
    std::string StrLine;
    std::string Header;
    int line = 0;
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '+': case '>':
                current_size = 0;
                continue;
            case '@':
                current_size = 0;
                break;
            case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }
        
        for (auto base: StrLine)
        {
            if (base == '\0') break;
                        
            if (base == '\n' || base == ' ' || base == '@') continue;
 
            
            
            kmer_read_c(base, klen, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen) continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
           
            if (iftarget)
            {
                update_counter(target_map, target_map_mul,larger_kmer, samplevecs);
                
                
            }
            else if (ifmask == 0 || base < 'a')
            {
                totalkmers++;
                if (target_map_nt[larger_kmer]<MAXCOUNT) target_map_nt[larger_kmer] ++;
                
            }
        }
    }
        
    fastafile.Close();
    
};

template <int dictsize>
void kmer_counter<dictsize>::read_files(std::vector<std::string>& inputs, std::vector<std::string>& outputs, std::vector<std::string>& prefs, int nthreads)
{
    
    inputfiles = inputs;
    outputfiles = outputs;
    prefixes = prefs;
    restfileindex = 0;
    
    std::vector<std::thread*> threads;
    
    for(int i=0; i< nthreads; ++i)
    {
        std::thread *newthread_ = new std::thread(&kmer_counter<dictsize>::read_file, this);
        threads.push_back(newthread_);
    }
    
    
    for(int i=0; i< nthreads; ++i)
    {
        threads[i]->join();
    }
    
}


template <int dictsize>
void kmer_counter<dictsize>::read_file()
{
    
    while (restfileindex < inputfiles.size())
    {
        
        Threads_lock.lock();
        
        int inputindex = restfileindex++ ;
        
        Threads_lock.unlock();
        
        if (inputindex >= inputfiles.size()) break;
        
        const char* inputfile = inputfiles[inputindex ].c_str();

        string prefix = "";
        if (prefixes.size() > inputindex)
        {
            prefix = prefixes[inputindex];
        }
        else if (inputfiles.size() > 1)
        {
            prefix = get_basename(inputfile);
        }
        
        int pathlen = (int)strlen(inputfile);
        
        if (
            (pathlen >= 5 && std::strcmp(inputfile + pathlen - 5, ".cram") == 0 )
            ||
            (pathlen >= 4 && std::strcmp(inputfile + pathlen - 4, ".bam") == 0)
            ||
            (pathlen >= 4 && std::strcmp(inputfile + pathlen - 4, ".sam") == 0)
            )
        {
            countint* samplevecs = (countint* )malloc(sizeof(countint) * totalkmers);
            
            memset(samplevecs, 0, sizeof(countint) * totalkmers);
            
            CramReader readsfile(inputfile);
            
            count_kmer(readsfile, samplevecs);
            
            for (int j = 0; j <outputfiles.size(); ++j)
            {
                auto outputfile = outputfiles[j] + prefix;
                pair<ull,ull> range;
                if (iftarget)
                {
                    range = targetranges[j];
                }
                else
                {
                    range = make_pair(1,totalkmers);
                }
            
                write(outputfile.c_str(), samplevecs, range);
            }
            
            free(samplevecs);
            
        }
        else if ( pathlen > 2 && strcmp(inputfile+(pathlen-3),".gz") == 0 )
        {
            countint* samplevecs = (countint* )malloc(sizeof(countint) * totalkmers);
            
            memset(samplevecs, 0, sizeof(countint) * totalkmers);
            
            fastq readsfile(inputfile);
            
            count_kmer(readsfile, samplevecs);
            
            for (int j = 0; j <outputfiles.size(); ++j)
            {
                auto outputfile = outputfiles[j] + prefix;
                pair<ull,ull> range;
                if (iftarget)
                {
                    range = targetranges[j];
                }
                else
                {
                    range = make_pair(1,totalkmers);
                }
            
                write(outputfile.c_str(), samplevecs, range);
            }
            
            free(samplevecs);
            
        }
            
        else
        {
            countint* samplevecs = (countint* )malloc(sizeof(countint) * totalkmers);
            
            memset(samplevecs, 0, sizeof(countint) * totalkmers);
            
            fasta readsfile(inputfile);
            
            count_kmer(readsfile, samplevecs);
            
            for (int j = 0; j <outputfiles.size(); ++j)
            {
                auto outputfile = outputfiles[j] + prefix;
                pair<ull,ull> range;
                if (iftarget)
                {
                    range = targetranges[j];
                }
                else
                {
                    range = make_pair(1,totalkmers);
                }
            
                write(outputfile.c_str(), samplevecs, range);
            }
            
            free(samplevecs);
        }
    }
}




template <int dictsize>
template <class typeint>
void kmer_counter<dictsize>::write(const char * outputfile, typeint* counts, pair<ull,ull> range)
{
    
    FILE *fwrite=fopen(outputfile, "w");
    
    
    if (fwrite==NULL)
    {
        std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
        
        std::_Exit(EXIT_FAILURE);
    }
    
    int code_bit = 1;
    
    int digit = (int)floor(1.0*klen/code_bit);
    
    char kmer_seq[digit+1];
    kmer_seq[digit] = '\0';
    
    if (iftarget)
    {
        
        fprintf(fwrite,"@targetsize\t%llu\n" , range.second -  range.first);
        
        for (auto i = range.first; i < range.second; ++i)
        {
            countint count = counts[i];
            kmer_int seq = kmer_records[i];
            
            if (count == 0) continue;
            
            for (int index = digit-1; index >= 0 ; --index)
            {
                kmer_seq[index] = "ACGT"[seq%4];
                seq /= 4;
            }
            fprintf(fwrite,"%s\t%d\n", kmer_seq, count);
        }
    }
    else
    {
        
        fprintf(fwrite,"@targetsize\t%llu\n", range.second -  range.first);
        
        for (const auto &kmer: target_map_nt)
        {
            
            auto seq = kmer.first;
            auto count = kmer.second;
            
            for (int index = digit-1; index >= 0 ; --index)
            {
                kmer_seq[index] = "ACGT"[seq%4];
                seq /= 4;
            }
            
            fprintf(fwrite,"%s\t%d\n", kmer_seq, count);
        }
    }
    
    fclose(fwrite);
    
    return ;
}



#endif /* KmerCounter_hpp */
