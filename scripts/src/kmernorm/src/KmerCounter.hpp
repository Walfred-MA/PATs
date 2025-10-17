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

#include "fasta.hpp"
#include "gzfile.hpp"
#include "KmerStruct.hpp"

using namespace std;

extern bool ifweight;
extern bool ifmask;
extern bool ifhighres;

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
    kmer32_dict kmer_hash;
    
    vector<vector<uint8_t>> norms;
    
    ull totalkmers = 1;
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
   
        void calculatenorm_highres(vector<uint>& norm, vector<float> &singletons);
    void calculatenorm( vector<uint16_t>& norm, vector<float> &singletons);

    void getnorm(string &inputfile, string &kmerfile, string &outputfile);
    void write(string &outputfile, vector<uint> & norm, vector<float> &singletons);
    void write(string &outputfile, vector<uint16_t> & norm, vector<float> &singletons);
    
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
           
            if (current_size < klen || (ifmask && base >= 'a') ) continue;
                                
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
            
            auto find = kmer_hash.find(larger_kmer);
            
            if (find == kmer_hash.end())
            {
                kmer_hash[larger_kmer] = totalkmers++;
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
    
    size_t sample_index = -1;
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
            
            if (find != kmer_hash.end() && norms[find->second][sample_index] < 255)
            {
                norms[find->second][sample_index] ++;
            }
        }
    }
    
    fastafile.Close();

};

inline void offsite_operation(ull &offsite, const vector<uint8_t> norm, vector<int>& vec_offsite, vector<uint16_t>&normmatrix, vector<int>& indexes, size_t colsize)
{
    size_t total = std::accumulate(norm.begin(), norm.end(), 0u);
    size_t index_ct = 0;
    if (2*total > norm.size())
    {
        offsite ++ ;
        
        for (size_t i = 0; i < norm.size(); ++i)
        {
            if (norm[i] == 0)
            {
                vec_offsite[i] --;
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
                if (value<MAX_UINT16) value ++;
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
                if (value<MAX_UINT16) value ++;
            }
        }
    }
}


template <int dictsize>
void kmer_counter<dictsize>::calculatenorm(vector<uint16_t>& normmatrix, vector<float> &singletons)
{
    size_t colsize = totalsamples;
    size_t rowsize = totalkmers;
    normmatrix.resize((colsize*(colsize+1))/2,0);
    singletons.resize(colsize,0.0);
    
    ull offsite = 0;
    vector<int> vec_offsite(colsize,0);
    vector<int> indexes(colsize, 0);
    
    for (size_t k =0; k<norms.size() ; ++k)
    {
        auto &norm = norms[k];
        auto max_index =std::max_element(norm.begin(), norm.end());
        auto max_value = (int)*max_index;
        
        size_t total = std::accumulate(norm.begin(), norm.end(), 0u);
         
        if (max_value == 0)
        {
            continue;
        }
        else if (ifweight && total == 1)
        {
            singletons[max_index - norm.begin()] += sington_weight;
        }
        
        else if (max_value == 1 )
        {
            offsite_operation(offsite, norm, vec_offsite, normmatrix, indexes, colsize);
        }
        else
        {
            for (size_t i =0; i< colsize ; ++i)
            {
                size_t index_corr = ((i + 1)*i)/2;
                normmatrix[colsize*i + i - index_corr] += norm[i] * norm[i] ;
                
                for (size_t j =i+1; j< colsize ; ++j)
                {
                    normmatrix[colsize*i + j - index_corr] += norm[i] * norm[j];
                }
                
            }
            
        }
    }
 
    size_t lastindex = 0;
    for (size_t i =0; i< colsize ; ++i)
    {
        size_t rowend =  lastindex + colsize - i;
        normmatrix[lastindex] += (uint) (singletons[i]+0.5);
        
        for (;lastindex < rowend; ++lastindex)
        {
            normmatrix[lastindex] += offsite + vec_offsite[i] + vec_offsite[colsize - (rowend - lastindex)] ;
        }
    }
    
}


inline void offsite_operation(ull &offsite, const vector<uint8_t> norm, vector<int>& vec_offsite, vector<uint>&normmatrix, vector<int>& indexes, size_t colsize)
{
    size_t total = std::accumulate(norm.begin(), norm.end(), 0u);
    size_t index_ct = 0;
    if (2*total > norm.size())
    {
        offsite ++ ;
        
        for (size_t i = 0; i < norm.size(); ++i)
        {
            if (norm[i] == 0)
            {
                vec_offsite[i] --;
                indexes[index_ct ++] = i;
            }
        }
        
        for (size_t i = 0; i < index_ct; ++i)
        {
            size_t index_i = indexes[i];
            size_t index_corr = ((index_i + 1)*index_i)/2;
            for (size_t j = i; j < index_ct; ++j)
            {
                normmatrix[index_i*colsize + indexes[j]  - index_corr] ++;
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
                normmatrix[index_i*colsize + indexes[j]  - index_corr] ++;
            }
        }
    }
}


template <int dictsize>
void kmer_counter<dictsize>::calculatenorm_highres(vector<uint> &normmatrix, vector<float> &singletons)
{
    size_t colsize = totalsamples;
    size_t rowsize = totalkmers;
    normmatrix.resize((colsize*(colsize+1))/2,0);
    singletons.resize(colsize,0.0);
    
    ull offsite = 0;
    vector<int> vec_offsite(colsize,0);
    vector<int> indexes(colsize, 0);
    
    for (size_t k =0; k<norms.size() ; ++k)
    {
        auto &norm = norms[k];
        auto max_index =std::max_element(norm.begin(), norm.end());
        auto max_value = (int)*max_index;
        
        size_t total = std::accumulate(norm.begin(), norm.end(), 0u);
         
        if (max_value == 0)
        {
            continue;
        }
        else if (ifweight && total == 1)
        {
            singletons[max_index - norm.begin()] += sington_weight;
        }
        
        else if (max_value == 1 )
        {
            offsite_operation(offsite, norm, vec_offsite, normmatrix, indexes, colsize);
        }
        else
        {
            for (size_t i =0; i< colsize ; ++i)
            {
                size_t index_corr = ((i + 1)*i)/2;
                normmatrix[colsize*i + i - index_corr] += norm[i] * norm[i] ;
                
                for (size_t j =i+1; j< colsize ; ++j)
                {
                    normmatrix[colsize*i + j - index_corr] += norm[i] * norm[j];
                }
                
            }
            
        }
    }
 
    size_t lastindex = 0;
    for (size_t i =0; i< colsize ; ++i)
    {
        size_t rowend =  lastindex + colsize - i;
        normmatrix[lastindex] += (int) (singletons[i]+0.5);
        
        for (;lastindex < rowend; ++lastindex)
        {
            normmatrix[lastindex] += offsite + vec_offsite[i] + vec_offsite[colsize - (rowend - lastindex)] ;
        }
    }
    
}

template <int dictsize>
void kmer_counter<dictsize>::getnorm(string &inputfile, string &kmerfile, string &outputfile)
{

    if (kmerfile.size() == 0)
    {
        kmerfile = inputfile;
    }


    fasta targetfile(kmerfile.c_str());
    
    read_target(targetfile);
    
    fasta readsfile(inputfile.c_str());
    
    read_file(readsfile);
    
    if (totalsamples < 200000 || ifhighres)
    {
        vector<float> singletons;
        vector<uint> thenorm;
        
        calculatenorm_highres(thenorm, singletons);
        write(outputfile, thenorm, singletons);
    }
    else
    {
        vector<float> singletons;
        vector<uint16_t> thenorm;
        
        calculatenorm(thenorm, singletons);
        write(outputfile, thenorm, singletons);
    }
}

template <int dictsize>
void kmer_counter<dictsize>::write(string &outputfile, vector<uint16_t> & norm, vector<float> &singletons)
{
    
    uint size = (int) floor(sqrt (norm.size() * 2 ));
    
    gzFile gz_out = gzopen(outputfile.c_str(), "wb6");
    
    if (fwrite==NULL)
    {
        std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
        
        std::_Exit(EXIT_FAILURE);
    }
    
    uint nextrowsize = size;
    uint index = 0;
    for (auto value: norm)
    {
        auto out = to_string(value)+",";
        gzwrite(gz_out, out.c_str(), out.size());
        if (++index == nextrowsize)
        {
            gzwrite(gz_out,"\n", 1);
            nextrowsize --;
            index = 0;
        }
    }
    
    gzclose(gz_out);
    
    if (ifweight)
    {
        gzFile gz_out = gzopen((outputfile+"_singleton.gz").c_str(), "wb6");
        
        if (fwrite==NULL)
        {
            std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
            
            std::_Exit(EXIT_FAILURE);
        }
        
        for (auto value: singletons)
        {
            auto out = to_string(int(sington_weight_r*value+0.5))+",";
            gzwrite(gz_out, out.c_str(), out.size());
        }
        
        gzclose(gz_out);
    }

    return ;
}

template <int dictsize>
void kmer_counter<dictsize>::write(string &outputfile, vector<uint> & norm, vector<float> &singletons)
{
    
    uint size = (int) floor(sqrt (norm.size() * 2 ));
    
    gzFile gz_out = gzopen(outputfile.c_str(), "wb6");
    
    if (fwrite==NULL)
    {
        std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
        
        std::_Exit(EXIT_FAILURE);
    }
    
    uint nextrowsize = size;
    uint index = 0;
    for (auto value: norm)
    {
        auto out = to_string(value)+",";
        gzwrite(gz_out, out.c_str(), out.size());
        if (++index == nextrowsize)
        {
            gzwrite(gz_out,"\n", 1);
            nextrowsize --;
            index = 0;
        }
    }
    
    gzclose(gz_out);
    
    if (ifweight)
    {
        gzFile gz_out = gzopen((outputfile+"_singleton.gz").c_str(), "wb6");
        
        if (fwrite==NULL)
        {
            std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
            
            std::_Exit(EXIT_FAILURE);
        }
        
        for (auto value: singletons)
        {
            auto out = to_string(int(sington_weight_r*value+0.5))+",";
            gzwrite(gz_out, out.c_str(), out.size());
        }
        
        gzclose(gz_out);
    }
    

    return ;
}


#endif /* KmerCounter_hpp */
