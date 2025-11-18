//
//  KmerFilter.cpp
//  merge_kmerlist
//
//  Created by Wangfei MA on 9/6/23.
//

#ifndef KmerFilter_hpp
#define KmerFilter_hpp

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


using kmer64_set_nt = std::unordered_set<u128, hash_128> ;
using kmer32_set_nt = std::unordered_set<ull>;

template <int dictsize>
class kmer_filter
{
    using kmer_int = typename std::conditional<(dictsize>32), u128, ull>::type;
    using kmer_set_type_nt = typename std::conditional<(dictsize>32), kmer64_set_nt, kmer32_set_nt>::type;
    
    kmer_set_type_nt target_map_nt;
    kmer_set_type_nt exclude_map_nt;
    
    ull totalkmers = 1;
    kmer_int *kmer_records;
    
    int klen = 31 , knum = 0;
    bool iftarget = 0;
    vector<pair<ull,ull>> targetranges;
    uint targetindex = 0;
    std::atomic_uint restfileindex ;
    std::mutex Threads_lock;
    std::mutex Map_lock;   
 
    std::vector<std::thread*> threads;
    std::vector<std::string> inputfiles;
    std::vector<std::string> outputfiles;
    std::vector<std::string> prefixes;

public:
    
    kmer_filter (int kmersize): klen(31)
    {
        kmer_records = (kmer_int *) malloc(1000);
    };
    ~kmer_filter()
    {
        free(kmer_records);
    };
    template <class typefile>
    void read_counttarget(typefile &fastafile);
    
    void read_target();
        
    template <class typefile>
    void filter_kmer(typefile &fastafile);
    
    void read_files(std::vector<std::string>& inputfiles, std::vector<std::string>& outputfiles, std::vector<std::string>& prefixes ,int numthread);
    
    void read_targets(std::vector<std::string>& inputfiles,int numthread);
    
    void read_file();
            
    void write(const char* outfile);

};



template <typename T>
static void kmer_read_c_(char base, int klen, int &current_size, T &current_kmer, T &reverse_kmer)
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
void kmer_filter<dictsize>::read_counttarget(typefile &fastafile)
{
        
    int current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    iftarget = 1;
    
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
            
            if (base == '\n' || base == ' ' || base == '\t') continue;

            kmer_read_c_(base, klen, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen) continue;
                                
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
            Map_lock.lock();
            target_map_nt.insert(larger_kmer);
            Map_lock.unlock();
        }
    }
    
    fastafile.Close();
    
    
};
 

template <int dictsize>
void kmer_filter<dictsize>::read_targets(std::vector<std::string>& targets,int nthreads)
{
    
    restfileindex = 0;
    inputfiles = targets;
    
    std::vector<std::thread*> threads;
    
    for(int i=0; i< nthreads; ++i)
    {
        std::thread *newthread_ = new std::thread(&kmer_filter<dictsize>::read_target, this);
        threads.push_back(newthread_);
    }
    
    
    for(int i=0; i< nthreads; ++i)
    {
        threads[i]->join();
    }
        
}

template <int dictsize>
void kmer_filter<dictsize>::read_target()
{
    
    while (restfileindex < inputfiles.size())
    {
        
        Threads_lock.lock();
        
        int inputindex = restfileindex++ ;
        
        Threads_lock.unlock();
        
        if (inputindex >= inputfiles.size()) break;
                        
        const char* inputfile = inputfiles[inputindex ].c_str();
        
        fasta readsfile(inputfile);
        
        read_counttarget(readsfile);
        
    }
    
}





template <int dictsize>
template <class typefile>
void kmer_filter<dictsize>::filter_kmer(typefile &fastafile)
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
            case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }
        
        for (auto base: StrLine)
        {
            if (base == '\0') break;
                        
            if (base == '\n' || base == ' ' || base == '\t') continue;
 
            kmer_read_c_(base, klen, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen) continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;

            Map_lock.lock();
            if ( target_map_nt.find(larger_kmer) != target_map_nt.end() )
            {
                    target_map_nt.erase(larger_kmer);
                    exclude_map_nt.insert(larger_kmer);
            }
            Map_lock.unlock();
        }
    }
        
    fastafile.Close();
};

template <int dictsize>
void kmer_filter<dictsize>::read_files(std::vector<std::string>& inputs, std::vector<std::string>& outputs, std::vector<std::string>& prefs, int nthreads)
{
    
    inputfiles = inputs;
    outputfiles = outputs;
    prefixes = prefs;
    restfileindex = 0;
    
    std::vector<std::thread*> threads;
    
    for(int i=0; i< nthreads; ++i)
    {
        std::thread *newthread_ = new std::thread(&kmer_filter<dictsize>::read_file, this);
        threads.push_back(newthread_);
    }
    
    
    for(int i=0; i< nthreads; ++i)
    {
        threads[i]->join();
    }
    
    write(outputs[0].c_str());
    
}


template <int dictsize>
void kmer_filter<dictsize>::read_file()
{
    while (restfileindex < inputfiles.size())
    {
        
        Threads_lock.lock();
        
        int inputindex = restfileindex++ ;
        
        Threads_lock.unlock();
        
        if (inputindex >= inputfiles.size()) break;
                        
        const char* inputfile = inputfiles[inputindex ].c_str();
        
        fasta readsfile(inputfile);
        
        filter_kmer(readsfile);
        
    }
        
}




template <int dictsize>
void kmer_filter<dictsize>::write(const char * outputfile)
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
   
    for (auto seq: target_map_nt)
    {
        for (int index = digit-1; index >= 0 ; --index)
        {
            kmer_seq[index] = "ACGT"[seq%4];
            seq /= 4;
        }
        
        fprintf(fwrite,">\n%s\n", kmer_seq);
    }
    
    fclose(fwrite);
    
    string excludefile = string(outputfile)+"_exclude.txt";
    FILE *fwrite2=fopen(excludefile.c_str(), "w");
    
    if (fwrite2==NULL)
    {
        std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
        
        std::_Exit(EXIT_FAILURE);
    }
    
    for (auto seq: exclude_map_nt)
    {
        for (int index = digit-1; index >= 0 ; --index)
        {
            kmer_seq[index] = "ACGT"[seq%4];
            seq /= 4;
        }
        
        fprintf(fwrite2,">\n%s\n", kmer_seq);
    }

    fclose(fwrite2);
    
    
    return ;
}

#endif /* KmerFilter_hpp */
