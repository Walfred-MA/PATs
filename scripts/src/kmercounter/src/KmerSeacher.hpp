//
//  kmers.hpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 10/14/21.
//  Copyright © 2021 USC_Mark. All rights reserved.
//

#ifndef KmerSeacher
#define KmerSeacher


#include <stdio.h>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <thread>
#include <atomic>
#include <mutex>
#include <memory>

#include "fasta.hpp"
#include "fastq.hpp"
#include "KmerStruct.hpp"

using namespace std;

using kmer64_dict_s = std::unordered_map<u128, uint16*, hash_128> ;
using kmer32_dict_s = std::unordered_map<ull, uint16* > ;

template <int dictsize>
class kmer_map
{
    using kmer_int = typename std::conditional<(dictsize>32), u128, ull>::type;
    //using target_maps_type = typename std::conditional<(dictsize>32), std::unordered_set<u128, hash_128>, std::unordered_set<ull>>::type;
    using kmer_dict_type = typename std::conditional<(dictsize>32), kmer64_dict_s, kmer32_dict_s>::type;
    using hotspot = std::tuple<int,int,int>;
    
    kmer_dict_type target_maps;
    
    int klen = 31 , knum = 0;
    bool iftarget = 0;
        vector<unsigned long long> targetlengths;
    ull totallength;
    uint targetindex = 0;
    const int hotspot_cutoff;
    
    std::atomic_uint restfileindex ;
    std::mutex Threads_lock;
    
    std::vector<std::thread*> threads;
    std::vector<std::string> inputfiles;
    std::vector<std::string> outputfiles;
    std::vector<std::string> prefixes;
        std::vector<std::string> allheaders;
public:
    
    kmer_map (int kmersize, int cutoff = 50): klen(31), hotspot_cutoff(cutoff)
    {
        allheaders.reserve(1000000);
        allheaders.push_back(string(""));
    };
    
    void read_target(const char* infile);
    
    template <class typefile>
    void read_searchtarget(typefile &fastafile);
    
    void read_files(std::vector<std::string>& inputs, std::vector<std::string>& outputs, std::vector<std::string>& prefs, int nthreads);
    
    void read_file();
            
    template <class typefile>
    void locate_kmer(typefile &fastafile, vector<hotspot*> &allhotspots, vector<ull>& sizes);
    
    void write(const char* outfile,  hotspot* allhotspots, ull sizes);
};


template <typename T>
static void kmer_read_s(char base, int klen, int &current_size, T &current_kmer, T &reverse_kmer)
{
    int converted = 0;
    T reverse_converted;
    
    if (base == '\n' || base == ' ' || base == '\t') return;
    
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

template <typename T1, typename T2>
inline static void initiate_counter_s(T2 &target_map, T1 &larger_kmer, uint16 index)
{
    
    typename T2::iterator map_find = target_map.find(larger_kmer);

    if (map_find == target_map.end())
    {
        target_map[larger_kmer] = (uint16*) malloc(sizeof(uint16)*2);
        target_map[larger_kmer][0] = 1;
        target_map[larger_kmer][1] = index ;
    }
    
    else
    {
        uint16* &data = map_find->second;

        if (std::find(data + 1, data + 1 + data[0], index) == map_find->second + 1 + data[0] )

        {
            if (data[0]%5 == 1)
            {
                data = (uint16*) realloc(data, sizeof(uint16)*(data[0]+1 + 5));
            }
            
            data[0]++;
            data[data[0]] = index ;
        }

    }
}

template <int dictsize>
template <class typefile>
void kmer_map<dictsize>::read_searchtarget(typefile &fastafile)
{
    
    int current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    
    targetindex ++;
    
    
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

            kmer_read_s(base, klen, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen || (ifmask && base >= 'a') ) continue;

            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
            
            initiate_counter_s(target_maps,larger_kmer,targetindex - 1);
                
        }

    }
    
    fastafile.Close();

    
};






template <int dictsize>
void kmer_map<dictsize>::read_target(const char* inputfile)
{

    fasta readsfile(inputfile);
    
    read_searchtarget(readsfile);
    
}

template <int dictsize>
template <class typefile>
void kmer_map<dictsize>::locate_kmer(typefile &fastafile, vector<hotspot*> &allhotspots, vector<ull>& allhotspot_sizes)
{
    int current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer=0;
    
    //const int hotspot_cutoff = 50;
    const int num_targets = (int) allhotspots.size();
    
    std::vector<std::vector<int>> kmer_posis;
    kmer_posis.resize(num_targets);
    for (auto &kmer_posi: kmer_posis)
    {
        kmer_posi.resize(hotspot_cutoff+1);
        std::fill(kmer_posi.begin(), kmer_posi.end(), -10000000);
    }
    
    std::vector<int> num_kmers(num_targets, 0);
    int posi = 0;
    int posi_start = 0;
    //uint64_t ifmasked = 0;
    //int num_masked = 0 ;
    std::string StrLine;


    int header = 0;
    
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '@':  case '+':
                current_size = 0;
                continue;
            case ' ': case '\n': case '\t':
                continue;
            case '>':
                current_size = 0;
                Threads_lock.lock();
                allheaders.push_back(StrLine.substr(1));
                header = allheaders.size() - 1;
                Threads_lock.unlock();
                std::fill(num_kmers.begin(), num_kmers.end(), 0);
                for (auto &kmer_posi: kmer_posis)
                {
                    std::fill(kmer_posi.begin(), kmer_posi.end(), -10000000);
                }
                posi = 0;
                continue;
            default:
                break;
        }
        
        for (auto base: StrLine)
        {
            if (base == '\0')
            {
                current_size = 0;
                header = 0;
                std::fill(num_kmers.begin(), num_kmers.end(), 0);
                for (auto &kmer_posi: kmer_posis)
                {
                    std::fill(kmer_posi.begin(), kmer_posi.end(), -10000000);
                }
                posi = 0;
                
                break;
            }
                        
            if (base == '\n' || base == ' ' || base == '\t') continue;
            
            posi++;
            
            kmer_read_s(base, klen, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen) continue;
            //num_masked=!!ifmasked;while(ifmasked&=ifmasked-1)num_masked++;
           
            //if (2*num_masked > klen) continue;
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
            
            auto find = target_maps.find(larger_kmer);
            if (find == target_maps.end()) continue;
           
 
            uint16* allindex = find->second;

            for (size_t j = 0; j < allindex[0]; ++j)
            {
                uint16 i = allindex[j+1] ;
                hotspot*& hotspots = allhotspots[i];
                ull& hotspots_size = allhotspot_sizes[i];
                auto& kmer_posi = kmer_posis[i];
                auto& num_kmer = num_kmers[i];
                
                posi_start = kmer_posi[num_kmer];
                
                kmer_posi[num_kmer] = posi;
                
                num_kmer = (num_kmer + 1) % hotspot_cutoff;
               
 
                if (abs(posi - posi_start) < 1000)
                {
                    if (hotspots_size &&
                        header == std::get<0>(hotspots[hotspots_size-1]) && abs(posi_start - std::get<2>(hotspots[hotspots_size-1]) ) < 1000 )
                    {
                        std::get<2>(hotspots[hotspots_size-1]) = posi ;
                    }
                    
                    else
                    {
    
                        if (hotspots_size%1000 == 1)
                        {
                            hotspots = (tuple<int,int,int>*) realloc(hotspots, (hotspots_size+1000)*sizeof(tuple<int,int,int>));
                        }
                        
                        hotspots[hotspots_size++] = std::make_tuple(header, posi_start, posi);
 
                    }
                }
            }
        }
    }
        
    fastafile.Close();
};


template <int dictsize>
void kmer_map<dictsize>::read_files(std::vector<std::string>& inputs, std::vector<std::string>& outputs, std::vector<std::string>& prefs, int nthreads)
{
    
    inputfiles = inputs;
    outputfiles = outputs;
    prefixes = prefs;
    restfileindex = 0;
        
    std::vector<std::thread*> threads;
    
    for(int i=0; i< nthreads; ++i)
    {
        std::thread *newthread_ = new std::thread(&kmer_map<dictsize>::read_file, this);
        threads.push_back(newthread_);
    }
    
    
    for(int i=0; i< nthreads; ++i)
    {
        threads[i]->join();
    }
    
}

template <int dictsize>
void kmer_map<dictsize>::read_file()
{
    
    while (restfileindex < inputfiles.size())
    {
        Threads_lock.lock();
        
        int inputindex = restfileindex ++;
        
        Threads_lock.unlock();
        
        if (inputindex >= inputfiles.size()) break;
        
        vector<hotspot*> allhotspots;
        allhotspots.resize(targetindex);
        
        vector<ull> allhotspots_size;
        allhotspots_size.resize(targetindex);
        
        for(hotspot* &hotspots:allhotspots)
        {
            hotspots = (tuple<int,int,int>*) malloc(sizeof(tuple<int,int,int>));
        }
        
        const char* inputfile = inputfiles[inputindex].c_str();
               
	string prefix = "";
	if (prefixes.size() > inputindex ) prefix = prefixes[inputindex];
 
        int pathlen = (int)strlen(inputfile);
        
        fasta readsfile(inputfile);
        
        if ( pathlen > 2 && strcmp(inputfile+(pathlen-3),".gz") == 0 )
        {
            
            fastq readsfile(inputfile);
            
            locate_kmer(readsfile, allhotspots, allhotspots_size);
            
        }
            
        else
        {
            
            fasta readsfile(inputfile);
            
            locate_kmer(readsfile, allhotspots, allhotspots_size);
        }
        
        
        for (int j = 0; j <outputfiles.size(); ++j)
        {
            auto outputfile = outputfiles[j] + prefix + "_hotspot.txt";
            
            write(outputfile.c_str(), allhotspots[j], allhotspots_size[j]);
        }
        
        for(hotspot* &hotspots:allhotspots)
        {
            free(hotspots);
        }
        
        allhotspots.clear();
        
    }
        
}




template <int dictsize>
void kmer_map<dictsize>::write(const char * outputfile, hotspot* hotspots, ull size)
{
    
    FILE *fwrite=fopen(outputfile, "w");
    if (fwrite==NULL)
    {
        std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
        
        std::_Exit(EXIT_FAILURE);
    }
            
    
    for (size_t i = 0; i < size; ++i)
    {
        hotspot hotspot = hotspots[i];
        fprintf(fwrite,"%s\t%d\t%d\n", allheaders[std::get<0>(hotspot)].c_str(), std::get<1>(hotspot), std::get<2>(hotspot));
        
    }

    fclose(fwrite);
    
    return ;
}


#endif /* kmers_hpp */
