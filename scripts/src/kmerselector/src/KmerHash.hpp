//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License.
//  If you use this code, please cite our work.
//

#ifndef KmerHash_hpp
#define KmerHash_hpp

#include <stdio.h>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <utility>
#include <functional>
#include <cstdint>
#include <memory>
#include <cstring>
#include <unordered_map>
#include <cmath>
#include <iostream>
#include <fstream>
#include <atomic>
#include <mutex>
#include <thread>
#include "config.hpp"


#define item_t uint
#define uint unsigned int
#define uint16 unsigned short

typedef unsigned long long ull;
typedef unsigned char  uint8;

using namespace std;



class Kmer32_hash
{
public:
    Kmer32_hash(): modsize(large_prime)
    {
        memset(key_sizes, 0 , modsize);
        memset(key_starts, 0 , modsize);
    };
    
    ~Kmer32_hash()
    {
        for (int i = 0 ; i < modsize; ++i)
        {
            if (key_sizes[i])
            {
                free(items[i]);
            }
        }
        free(items);
        free(key_sizes);
        free(key_starts);
    }
    
    uint* add(const ull kmer_int, const uint index);
    bool ifadd(const ull kmer_int, const uint index);
    uint* find(const ull kmer_int);
    uint findhash(const ull kmer_int);
    uint initiate();
    const size_t modsize;

private:
    item_t** items = new item_t*[modsize];
    uint16* key_sizes = new uint16[modsize];
    uint* key_starts = new uint[modsize];
};


class Kmer_hash
{
public:
    Kmer_hash()
        : Hash1(make_unique<uint[]>(large_prime)), Hash2_size(make_unique<uint8_t[]>(large_prime))
    {
        std::fill_n(Hash1.get(), large_prime, 0);
        std::fill_n(Hash2_size.get(), large_prime, 0);
    };
    ~Kmer_hash()
    {
        
    };

    void savehash(const std::string& outputfile, ull totalkmers);
    ull loadhash(const std::string& outputfile);
    void preview(const ull kmer_int);
    bool iterate(uint &hash1, uint8_t &hash2index, ull &kmer);
    uint initiate();
    uint finalize();
    bool add(const ull kmer_int);
    uint findhash(ull kmer_int);
    ull totalkmer = 0;
    uint32_t totaleles = 0;
    uint32_t bucketcollision = 0;
    unique_ptr<uint[]> Hash1;
    unique_ptr<uint8_t[]>Hash2_size;
    vector<uint> Hash2;
};




#endif /* KmerHash_hpp */
