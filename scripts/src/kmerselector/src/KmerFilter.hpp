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
#include "KmerStruct.hpp"
#include "gzfile.hpp"
#include "Blacklist.hpp"

using namespace std;
extern bool singletarget;
#define nbase 3
using kmer64_set_nt = std::unordered_map<u128, bool, hash_128> ;
//using kmer32_set_nt = std::unordered_set<ull>;


template <int dictsize>
class kmer_filter
{
    using kmer_int = typename std::conditional<(dictsize>32), u128, ull>::type;
    using kmer_set_type_nt = typename std::conditional<(dictsize>32), kmer64_set_nt, kmer32_set_nt>::type;
    
    //kmer_set_type_nt target_map_nt;
    //kmer_set_type_nt exclude_map_nt;
    
    ull totalkmers = 1;
    kmer_int *kmer_records;
    
    int knum = 0;
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
    std::vector<std::string> targetfiles;

public:
    
    kmer_filter (int kmersize)
    {
        kmer_records = (kmer_int *) malloc(1000);
    };
    ~kmer_filter()
    {
        free(kmer_records);
    };
    template <class typefile>
    void read_counttarget(typefile &fastafile,  kmer_set_type_nt &target_map_nt, kmer_set_type_nt &exclude_map_nt);
    
    void read_target();
        
    template <class typefile>
    void filter_kmer(typefile &fastafile, kmer_set_type_nt &target_map_nt, kmer_set_type_nt &exclude_map_nt);
    
    void read_files(std::vector<std::string>& inputfiles, std::vector<std::string>& outputfiles, std::vector<std::string>& prefixes, std::vector<std::string>& targetfiles,int numthread);
    
    void read_targets(std::vector<std::string>& inputfiles,int numthread);
    
    void read_file();
            
    void write(const char* outfile, kmer_set_type_nt &target_map_nt, kmer_set_type_nt &exclude_map_nt);

};



template <typename T>
static void kmer_read_f_(char base, int &current_size, T &current_kmer, T &reverse_kmer)
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


template <typename T>
static void kmer_read_f_64(char base, int &current_size, T &current_kmer, T &reverse_kmer)
{
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
void kmer_filter<dictsize>::read_counttarget(typefile &fastafile, kmer_set_type_nt &target_map_nt, kmer_set_type_nt &exclude_map_nt)
{
        
    int current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    iftarget = 1;
    
    std::string StrLine(10000,'\0');
    BlacklistIntervals::HeaderRegion blacklist_region;
    unsigned long long sequence_offset = 0;
    
    while (fastafile.nextLine(StrLine))
    {
        if (StrLine.empty()) continue;

        switch (StrLine[0])
        {
            case '@':  case '+': case '>':
                current_size = 0;
                sequence_offset = 0;
                blacklist_region = blacklist_intervals.parse_header_region(StrLine);
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

            sequence_offset++;
            kmer_read_f_(base, current_size, current_kmer, reverse_kmer);
            
            if ( current_size < klen || (ifmask && base >= 'a') ) continue;
            if (blacklist_intervals.overlaps_sequence_span(blacklist_region, sequence_offset - klen, sequence_offset)) continue;
                                
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
           
            if ( exclude_map_nt.find(larger_kmer) != exclude_map_nt.end() ) continue;
 
            //target_map_nt.insert(larger_kmer);
            target_map_nt.try_emplace(larger_kmer, false);
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
        
        //read_counttarget(readsfile);
        
    }
    
}





template <int dictsize>
template <class typefile>
void kmer_filter<dictsize>::filter_kmer(typefile &fastafile, kmer_set_type_nt &target_map_nt, kmer_set_type_nt &exclude_map_nt)
{
    
    int current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    std::string StrLine (1000000, '\0');
    std::string Header;
    while (fastafile.nextLine(StrLine))
    {
        /*
        switch (StrLine[0])
        {
            case '@':  case '+': case '>':
                current_size = 0;
                break;
            case ' ': case '\n': case '\t':
                break;
            default:
                break;
        }
        */
        current_size = 0;
        for (auto base: StrLine)
        {
            if (base == '\0') break;
                        
            if (base == '\n' || base == ' ' || base == '\t') continue;
 
            kmer_read_f_64(base, current_size, current_kmer, reverse_kmer);
            
            if (current_size < klen) continue;
            //auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
            //auto larger_kmer  = current_kmer;
            target_map_nt[current_kmer] = 1;

        }
    }
        
    fastafile.Close();
};

template <int dictsize>
void kmer_filter<dictsize>::read_files(std::vector<std::string>& inputs, std::vector<std::string>& outputs, std::vector<std::string>& prefs, std::vector<std::string>& targets, int nthreads)
{
    
    inputfiles = inputs;
    outputfiles = outputs;
    targetfiles = targets;
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
        
}


template <int dictsize>
void kmer_filter<dictsize>::read_file()
{
    
    while (restfileindex < targetfiles.size())
    {
        
        Threads_lock.lock();
        
        int inputindex = restfileindex++ ;
        
        Threads_lock.unlock();
        
        if (inputindex >= targetfiles.size()) break;
    
        std::string targetfile = targetfiles[inputindex];
        std::string outputfolder = outputfiles[inputindex];
       
        cout << "merging results on: "<< targetfile << endl;
 
        kmer_set_type_nt target_map_nt;
        kmer_set_type_nt exclude_map_nt;
                
        string cachefile = targetfile + "_allkmer.cache";
        if (FILE *f = fopen(cachefile.c_str(), "r"); f && !blacklist_intervals.enabled())
        {
            fclose(f);
            read_cache_to_map(cachefile.c_str(), target_map_nt);
            std::remove(cachefile.c_str());
        }
        else
        {
            if (f) fclose(f);
            if (blacklist_intervals.enabled()) std::remove(cachefile.c_str());
            if (targetfile.size())
            {
                fasta targetfile_(targetfile.c_str());
                read_counttarget(targetfile_, target_map_nt, exclude_map_nt);
            }
        }
        
        for (std::string &prefix: prefixes)
        {
            auto outputfile_ = outputfolder + prefix;
            gzfile outputfile(outputfile_.c_str());
            filter_kmer(outputfile, target_map_nt, exclude_map_nt);
        }
        
        if (outputfiles.size() == 1 && singletarget)
        {
            targetfile = outputfiles[0];
        }
        else
        {
            targetfile = targetfile + "_kmers.txt";
        }
        
    
        write(targetfile.c_str(), target_map_nt, exclude_map_nt);
        
    }
        
}




template <int dictsize>
void kmer_filter<dictsize>::write(const char * outputfile, kmer_set_type_nt &target_map_nt, kmer_set_type_nt &exclude_map_nt)
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
    
 
    for (auto [seq0, iffilter]: target_map_nt)
    {
        //if ( exclude_map_nt.find(seq) != exclude_map_nt.end() ) continue;
        ull seq = seq0;
        
        if (iffilter) continue;
        for (int index = digit-1; index >= 0 ; --index)
        {
            kmer_seq[index] = "ACGT"[seq%4];
            seq /= 4;
        }
        
        fprintf(fwrite,">\n%s\n", kmer_seq);
    }
    
    fclose(fwrite);

    /*
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
    */
    
    return ;
}

#endif /* KmerFilter_hpp */
