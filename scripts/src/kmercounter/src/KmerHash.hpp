//
//  KmerHash.hpp
//  KmerHasher
//
//  Created by walfred on 12/15/24.
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
#include <unordered_map>
#include <cmath>
#include <iostream>
#include <cstring>

#include "config.hpp"

#define MAX_UINT22 4194303
#define uint unsigned int

struct uint40
{
    uint40() : first(0), second(0) {}
    uint40 (const ull x): first( x % MAX_UINT32 ) , second( x / MAX_UINT32 )
    {};
    uint40 (const uint32_t x, const uint8_t y): first( x ) , second(y)
    {};
    
    auto operator==(const uint40& other) const
    {
        return (this->first == other.first && this->second == other.second);
    }
 
    uint first;
    uint8 second;
};

typedef std::pair<uint40,uint16> item40_t ;
#define MAX_UINT22 4194303
#define uint unsigned int
#define item_default make_pair<item40_t>(MAX_INT16, 0)


using namespace std;


class Small_hash
{
public:
    Small_hash(): modsize(primeint22)
    {
         memset(key_sizes, 0 , modsize);
    };
    
    ~Small_hash()
    {
        free(key_sizes);
        free(items);
    }
    
    uint16* add(const ull kmer_int, const uint index);
    uint16* find(const ull kmer_int);
    const size_t modsize;

private:
    item40_t* items = new item40_t[modsize * 10];
    uint16_t* key_sizes = new uint16_t[modsize];
    uint32_t* key_index = new uint32_t[modsize];
    uint32_t totalsize = 0;
};

class Kmer_hash
{
public:
    Kmer_hash()
    {
        std::fill_n(Hash1.get(), large_prime, MAX_UINT32);
    };
    ~Kmer_hash()
    {
        
    }
    
    void preview(const ull kmer_int);
    void initiate();
    bool add(const ull kmer_int);
    uint* find(const ull kmer_int);
    ull totalkmer = 0;

private:
    unique_ptr<uint[]> Hash1 = make_unique<uint[]>(large_prime );
    vector<uint> Hash2;
    vector<uint8_t> Hash2_size;
    uint32_t totalcollision = 0;
    uint32_t bucketcollision = 0;
};




#endif /* KmerHash_hpp */

