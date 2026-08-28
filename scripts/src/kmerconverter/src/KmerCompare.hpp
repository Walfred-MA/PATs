//
//  KmerCompare.hpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 9/20/23.
//  Copyright © 2023 USC_Mark. All rights reserved.
//

#ifndef KmerCompare_hpp
#define KmerCompare_hpp

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

using namespace std;


using kmer64_comp_nt = std::unordered_map<u128, pair<uint8, uint16>, hash_128> ;
using kmer32_comp_nt = std::unordered_map<ull, pair<uint8, uint16> >;

template <int dictsize>
class kmer_compare
{
    using kmer_int = typename std::conditional<(dictsize>32), u128, ull>::type;
    using kmer_compare_nt = typename std::conditional<(dictsize>32), kmer64_comp_nt, kmer32_comp_nt>::type;
    
    kmer_compare_nt target_map;
    
    
    ull totalkmers = 1;
    
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
    
    kmer_compare (int kmersize): klen(31)
    {};
    ~kmer_compare()
    {};
    template <class typefile>
    void read_counttarget(typefile &fastafile);
    
    void read_target(const char* infile);
        
    template <class typefile>
    void count_kmer(typefile &fastafile, uint8* samplevecs);
    
    void read_files(std::vector<std::string>& inputfiles, std::vector<std::string>& outputfiles, std::vector<std::string>& prefixes ,int numthread);
    
    void read_file();
            
    void write(const char* outfile);

};




template <typename T>
static void kmer_read_compare(char base, int klen, int &current_size, T &current_kmer, T &reverse_kmer)
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
        current_size = -1;
        current_kmer = 0;
        reverse_kmer = 0;
    }
}

template <int dictsize>
template <class typefile>
void kmer_compare<dictsize>::read_counttarget(typefile &fastafile)
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
            case '\n': case ' ': case '\t':
                continue;
            default:
                break;
        }
        
        for (auto base: StrLine)
        {
            if (base == '\0') break;
            
            if (base == '\n' || base == ' ' || base == '\t') continue;

            kmer_read_compare(base, klen, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen) continue;
                                
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
            
            if (target_map[larger_kmer].first<255) target_map[larger_kmer].first ++;
                
        }
    }
    
    fastafile.Close();
    
    targetrange.second = totalkmers ;
    
};
 
template <int dictsize>
void kmer_compare<dictsize>::read_target(const char* inputfile)
{
    
    fasta readsfile(inputfile);
    
    read_counttarget(readsfile);
    
}



template <int dictsize>
template <class typefile>
void kmer_compare<dictsize>::count_kmer(typefile &fastafile, uint8* samplevecs)
{
    
    int current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    //uint64_t ifmasked = 0;
    //int num_masked = 0 ;
    std::string StrLine;
    std::string Header;
    
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '@':  case '+': case '>':
                current_size = 0;
                continue;
            case '\n': case ' ': case '\t':
                continue;
            default:
                break;
        }
        
        for (auto base: StrLine)
        {
            if (base == '\0') break;
                        
            if (base == '\n' || base == ' ' || base == '\t' ) continue;
 
            kmer_read_compare(base, klen, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen) continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
            
            if (target_map.find(larger_kmer) != target_map.end() && target_map[larger_kmer].second < 0xFFFF ) target_map[larger_kmer].second ++;
        }
    }
        
    fastafile.Close();
};

template <int dictsize>
void kmer_compare<dictsize>::read_files(std::vector<std::string>& inputs, std::vector<std::string>& outputs, std::vector<std::string>& prefs, int nthreads)
{
    
    inputfiles = inputs;
    outputfiles = outputs;
    prefixes = prefs;
    restfileindex = 0;
    
    std::vector<std::thread*> threads;
    
    for(int i=0; i< nthreads; ++i)
    {
        std::thread *newthread_ = new std::thread(&kmer_compare<dictsize>::read_file, this);
        threads.push_back(newthread_);
    }
    
    
    for(int i=0; i< nthreads; ++i)
    {
        threads[i]->join();
    }
    
}


template <int dictsize>
void kmer_compare<dictsize>::read_file()
{
    
    while (restfileindex < inputfiles.size())
    {
        
        Threads_lock.lock();
        
        int inputindex = restfileindex++ ;
        
        Threads_lock.unlock();
        
        if (inputindex >= inputfiles.size()) break;
        
        uint8* samplevecs = (uint8* )malloc(sizeof(uint8) * 0);
        
        memset(samplevecs, 0, sizeof(uint8) * totalkmers);
        
        const char* inputfile = inputfiles[inputindex ].c_str();
        
        int pathlen = (int)strlen(inputfile);
        
        if ( pathlen > 2 && strcmp(inputfile+(pathlen-3),".gz") == 0 )
        {
            
            fastq readsfile(inputfile);
            
            count_kmer(readsfile, samplevecs);
            
        }
            
        else
        {
            
            fasta readsfile(inputfile);
            
            count_kmer(readsfile, samplevecs);
        }
        
        string prefix = "";
        
        if (prefixes.size() > inputindex) prefix = prefixes[inputindex];
        
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
            
            write(outputfile.c_str());
        }
        
        free(samplevecs);
        
    }
        
}




template <int dictsize>
void kmer_compare<dictsize>::write(const char * outputfile)
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
    
    fprintf(fwrite,"@targetsize\t%llu\n" , target_map.size());
    
    for (auto &kmer: target_map)
    {
        uint8 count1 = kmer.second.first;
        uint16 count2 = kmer.second.second;
        kmer_int seq = kmer.first;
        
        if ( count1 == 0 || count1 == 255 ) continue;
        
        for (int index = digit-1; index >= 0 ; --index)
        {
            kmer_seq[index] = "ACGT"[seq%4];
            seq /= 4;
        }
        
        fprintf(fwrite,"%s\t%d\t%d\n", kmer_seq, count1, count2);
    }
    
    
    fclose(fwrite);
    
    return ;
}


#endif /* KmerCompare_hpp */

