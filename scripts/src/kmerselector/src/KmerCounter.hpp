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

#include "fasta.hpp"
#include "KmerStruct.hpp"
#include "KmerHash.hpp"

using namespace std;

extern bool ifmask;

using kmer64_dict = std::unordered_map<u128, uint, hash_128> ;
using kmer32_dict = std::unordered_map<ull, uint > ;

using kmer64_dict_nt = std::unordered_map<u128, uint8, hash_128> ;
using kmer32_dict_nt = std::unordered_map<ull, uint8>;

using kmer64_dict_mul = std::unordered_map<u128, uint*, hash_128> ;
using kmer32_dict_mul = std::unordered_map<ull, uint* > ;

template <int dictsize>
class kmer_counter
{
    using kmer_int = typename std::conditional<(dictsize>32), u128, ull>::type;
    using kmer_dict_type_nt = typename std::conditional<(dictsize>32), kmer64_dict_nt, kmer32_dict_nt>::type;
    
    //kmer32_dict target_map;
    Kmer32_hash Kmer_hash;
    ull totalkmers = 1;
    
    int knum = 0;
    bool iftarget = 0;
    vector<pair<ull,ull>> targetranges;
    uint targetindex = 0;
    std::atomic_uint restfileindex ;
    std::mutex Threads_lock;
    
    std::vector<std::thread*> threads;
    std::vector<std::string> inputfiles;
    std::vector<std::string> targetfiles;
    std::vector<std::string> outputfiles;
    std::vector<std::string> prefixes;

public:
    
    kmer_counter (int kmersize): Kmer_hash(large_prime)
    {
        
    };
    ~kmer_counter()
    {
    };
    template <class typefile>
    void read_counttarget(typefile &fastafile);
    
    template <class typefile>
    void count_target(typefile &fastafile, string prefix, kmer32_dict_nt &target_map_nt);
    
    void read_target(const char* infile);
        
    template <class typefile>
    void count_kmer(typefile &fastafile, uint8* samplevecs);
    
    void read_files(std::vector<std::string>& inputfiles, std::vector<std::string>& outputfiles, std::vector<std::string>& prefixes, std::vector<std::string>& targets,int numthread);
    
    void read_file();
            
    //void write(const char *outputfile, kmer32_dict &target_map, uint8* samplevecs, kmer_dict_type_nt &target_map_nt);
    void write(const char *outputfile, Kmer32_hash &target_map, uint8* samplevecs, kmer_dict_type_nt &target_map_nt);
};

inline static void update_counter(kmer32_dict &target_map, ull &larger_kmer, uint8* samplevecs)
{
    auto map_find = target_map.find(larger_kmer);

    if (map_find != target_map.end())
    {
        if (samplevecs[map_find->second] < 255) samplevecs[map_find->second] ++;
    }
}

inline static void update_counter(Kmer32_hash &target_map, ull &larger_kmer, uint8* samplevecs)
{
    auto map_find = target_map.find(larger_kmer);

    if (map_find != NULL)
    {
        if (samplevecs[*map_find] < 255) (samplevecs[*map_find]) ++;
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

            kmer_read_c(base, current_size, current_kmer, reverse_kmer);
           
            if (current_size < klen || (ifmask && base >= 'a') ) continue;
                                
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
            
            if (Kmer_hash.ifadd(larger_kmer, totalkmers)) totalkmers++;
                
        }
    }
    
    fastafile.Close();

    cout << "finishing reading file: "<< string(fastafile.filepath) << endl;
    cout << "total kmers: "<< totalkmers << endl;   

 
    targetrange.second = totalkmers ;
    
};

bool static inline isTarget(const std::string& strLine, const std::string& prefix)
{
    std::istringstream iss(strLine);
    std::vector<std::string> tokens;
    std::string token;
    
    while (std::getline(iss, token, '_')) {
        tokens.push_back(token);
    }

    if (tokens.size() < 3) {
        return false;
    }

    std::string extractedPrefix = tokens[1] + "_" + tokens[2];

    return extractedPrefix == prefix;
}

template <int dictsize>
template <class typefile>
void kmer_counter<dictsize>::count_target(typefile &fastafile, string prefix, kmer32_dict_nt &target_map_nt)
{
    
    int current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    //uint64_t ifmasked = 0;
    //int num_masked = 0 ;
    std::string StrLine;
    std::string Header;
    bool istarget = 0;
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '@':  case '+': case '>':
                current_size = 0;
                istarget = isTarget(StrLine, prefix);
                continue;
            case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }
        
        if (istarget == 0) continue;
        
        for (auto base: StrLine)
        {
            if (base == '\0') break;
                        
            if (base == '\n' || base == ' ') continue;
 
            kmer_read_c(base, current_size, current_kmer, reverse_kmer);
            
            if (current_size < klen) continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
           
            if (ifmask == 0 || base < 'a')
            {
                if (target_map_nt[larger_kmer]<255) target_map_nt[larger_kmer] ++;
                
            }
        }
    }
    
    
        
    fastafile.Close();
    
};
 
template <int dictsize>
void kmer_counter<dictsize>::read_target(const char* inputfile)
{
    
    fasta readsfile(inputfile);
    
    read_counttarget(readsfile);
    
}



template <int dictsize>
template <class typefile>
void kmer_counter<dictsize>::count_kmer(typefile &fastafile, uint8* samplevecs)
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
                        
            if (base == '\n' || base == ' ') continue;
 
            kmer_read_c(base, current_size, current_kmer, reverse_kmer);
            
            if (current_size < klen) continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
           
            update_counter(Kmer_hash, larger_kmer, samplevecs);
        }
    }
        
    fastafile.Close();
};

template <int dictsize>
void kmer_counter<dictsize>::read_files(std::vector<std::string>& inputs, std::vector<std::string>& outputs, std::vector<std::string>& prefs, std::vector<std::string>& targets, int nthreads)
{
    
    inputfiles = inputs;
    outputfiles = outputs;
    prefixes = prefs;
    targetfiles = targets;
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
        
        uint8* samplevecs = (uint8* )malloc(sizeof(uint8) * totalkmers);
        
        memset(samplevecs, 0, sizeof(uint8) * totalkmers);
        
        const char* inputfile = inputfiles[inputindex ].c_str();

        cout << "Working on sample: "<< string(inputfile) << endl;

        string prefix = "";
        if (prefixes.size() > inputindex) prefix = prefixes[inputindex];
        
        int pathlen = (int)strlen(inputfile);
        
        fasta readsfile(inputfile);
        
        count_kmer(readsfile, samplevecs);

        cout << "Outputing on sample: "<< string(inputfile) << endl;        

        for (int j = 0; j <outputfiles.size(); ++j)
        {
            kmer_dict_type_nt target_map_nt;
            
            fasta targetfile(targetfiles[j].c_str());
            
            count_target(targetfile,prefix, target_map_nt);
            
            auto outputfile = outputfiles[j] + prefix;
            pair<ull,ull> range;
            range = targetranges[j];
            
            write(outputfile.c_str(), Kmer_hash, samplevecs, target_map_nt);
        }
        
        free(samplevecs);
        
    }
        
}




template <int dictsize>
void kmer_counter<dictsize>::write(const char * outputfile, Kmer32_hash &target_map, uint8* samplevecs ,kmer_dict_type_nt &target_map_nt)
{
    
    gzFile gz_out = gzopen(outputfile, "wb6");

    if (fwrite==NULL)
    {
        std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
        
        std::_Exit(EXIT_FAILURE);
    }
   
    if (target_map_nt.size() == 0)
    {
        gzclose(gz_out);
        return;
    }


    for (const auto& [kmer, count] : target_map_nt)
    {
        
        auto loc = target_map.findvalue(kmer);
        
        if (count == samplevecs[loc] && count < 255 ) continue;
        
        std::string output = kmer_int_totext(kmer) + "\t" + std::to_string(count) + "\t" + std::to_string(samplevecs[loc] ) + "\n";
        gzwrite(gz_out, output.c_str(), output.size());
    }
    
    gzclose(gz_out);

    return ;
}


#endif /* KmerCounter_hpp */

