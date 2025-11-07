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
    Kmer32_hash kmer_hash;
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
    
    kmer_counter()
    {};
    ~kmer_counter()
    {
    };

    //ull read_counttarget(vector<std::string> &targetfiles);
    
    void read_counttarget(std::string &infile);
    
    template <class typefile>
    bool count_target(typefile &fastafile, string prefix, kmer32_dict_nt &target_map_nt);
    
    void read_targets(std::vector<std::string>& infiles);
        
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
    auto map_find = target_map.findhash(larger_kmer);

    if (map_find > 0)
    {
        if (samplevecs[map_find] < 255) (samplevecs[map_find]) ++;
    }
}

inline static void update_counter(Kmer_hash &Kmer_hash,ull &larger_kmer, uint8* samplevecs)
{
    uint index = Kmer_hash.findhash(larger_kmer);
    
    if (index > 0 && samplevecs[index] < MAX_COUNT) samplevecs[index]++;
    
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
/*
template <int dictsize>
ull kmer_counter<dictsize>::read_counttarget(vector<std::string> &targetfiles)
{
    
    int current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
        
    std::string StrLine;
    StrLine.resize(MAX_LINE);
   
    for (std::string &targetfile: targetfiles)
    {
        fasta fastafile(targetfile.c_str());

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
            
            for (char base: StrLine)
            {
                if (base == '\0') break;
                
                if (base == '\n' || base == ' ') continue;

                kmer_read_c(base, current_size, current_kmer, reverse_kmer);
                
                if (current_size < klen || (ifmask && base >= 'a')) continue;
                                    
                auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
                
                kmer_hash.preview(larger_kmer);
            }
        }
        
    }
    
    totalkmers = kmer_hash.initiate();
    
    for (std::string &targetfile: targetfiles)
    {
        fasta fastafile(targetfile.c_str());
        
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
            
            for (char base: StrLine)
            {
                if (base == '\0') break;
                
                if (base == '\n' || base == ' ') continue;

                kmer_read_c(base, current_size, current_kmer, reverse_kmer);
                
                if (current_size < klen || (ifmask && base >= 'a')) continue;
                                    
                auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
                
                kmer_hash.add(larger_kmer);
            }
        }
    }
    
    kmer_hash.finalize();
    
    return totalkmers;

};
*/


inline static void write_cache(const char* outputfile, const std::unordered_set<ull>& allkmers)
{
    std::ofstream ofs(outputfile, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error(std::string("Cannot open file for writing: ") + outputfile);
    }

    // 1) write number of elements
    uint64_t n = allkmers.size();
    ofs.write(reinterpret_cast<const char*>(&n), sizeof(n));

    // 2) write each ull
    for (ull k : allkmers) {
        ofs.write(reinterpret_cast<const char*>(&k), sizeof(k));
    }

    if (!ofs) {
        throw std::runtime_error(std::string("Error while writing cache file: ") + outputfile);
    }
}


inline static void read_cache_to_map(const char* inputfile, kmer32_dict_nt& dict)
{
    std::ifstream ifs(inputfile, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error(std::string("Cannot open file for reading: ") + inputfile);
    }

    uint64_t n = 0;
    ifs.read(reinterpret_cast<char*>(&n), sizeof(n));
    if (!ifs) {
        throw std::runtime_error(std::string("Error reading size from cache file: ") + inputfile);
    }

    dict.clear();
    dict.reserve(static_cast<size_t>(n * 1.3));  // small slack

    for (uint64_t i = 0; i < n; ++i) {
        ull k;
        ifs.read(reinterpret_cast<char*>(&k), sizeof(k));
        if (!ifs) {
            throw std::runtime_error(std::string("Error reading element from cache file: ") + inputfile);
        }
        dict.emplace(k, 0);   // value initialized as 0
    }
}

template <int dictsize>
void kmer_counter<dictsize>::read_counttarget(std::string &infile)
{
    
    fasta fastafile(infile.c_str());
    
    int current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    iftarget = 1;
    targetranges.emplace_back();
    auto & targetrange = targetranges[targetranges.size()-1];
    targetrange.first = totalkmers;
    
    unordered_set<ull> allkmers;
    
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
            
            allkmers.insert(larger_kmer);
        }
    }
    
    fastafile.Close();
    
    for (auto& larger_kmer: allkmers)
    {
        if (kmer_hash.ifadd(larger_kmer, totalkmers)) totalkmers++;
    }
    
    write_cache((infile+"_allkmer.cache").c_str(), allkmers);
    
    cout << "finishing loading kmers: "<< string(fastafile.filepath) << endl;
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
bool kmer_counter<dictsize>::count_target(typefile &fastafile, string prefix, kmer32_dict_nt &target_map_nt)
{
    
    int current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    bool iffind = 0;
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
        
        if (!istarget) continue;
        
        iffind = 1;
        
        for (auto base: StrLine)
        {
            if (base == '\0') break;
                        
            if (base == '\n' || base == ' ') continue;
 
            kmer_read_c(base, current_size, current_kmer, reverse_kmer);
            
            if (current_size < klen) continue;
            
            ull larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
             
            if (ifmask == 0 || base < 'a')
            {
                if (target_map_nt[larger_kmer]<255) target_map_nt[larger_kmer] ++;
            }
        }
    }
        
    fastafile.Close();
    
    return iffind;
};
 
template <int dictsize>
void kmer_counter<dictsize>::read_targets(std::vector<std::string>& infiles)
{
    
    for (std::string &infile: infiles)
    {
        cout << "loading kmers in: " << infile << endl;
        
        read_counttarget(infile);
    }
    
    kmer_hash.initiate();
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
            
            update_counter(kmer_hash, larger_kmer, samplevecs);
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
            
            read_cache_to_map((targetfiles[j]+"_allkmer.cache").c_str(), target_map_nt);
            
            fasta targetfile(targetfiles[j].c_str());
            
            if (!count_target(targetfile,prefix, target_map_nt))
            {
                target_map_nt.clear();
            }
            
            auto outputfile = outputfiles[j] + prefix;
            
            write(outputfile.c_str(), kmer_hash, samplevecs, target_map_nt);
            
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
        
        auto loc = target_map.findhash(kmer);
        
        if (count == samplevecs[loc] && count < 255 ) continue;
        
        std::string output = kmer_int_totext(kmer) + "\t" + std::to_string(count) + "\t" + std::to_string(samplevecs[loc] ) + "\n";
        gzwrite(gz_out, output.c_str(), output.size());
    }
    
    gzclose(gz_out);

    return ;
}


#endif /* KmerCounter_hpp */

