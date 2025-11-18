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
#include <fstream>

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

class Kmer_hash
{
public:
    Kmer_hash()
        : Hash1(make_unique<uint[]>(large_prime)), Hash2_size(make_unique<uint8_t[]>(large_prime))
    {
        std::fill_n(Hash1.get(), large_prime, 0);
        std::fill_n(Hash2_size.get(), large_prime, 0);
        values_mul.resize(MAX_UINT16+1);
        values_mul_bucksize.resize(MAX_UINT16+1,0);
    };
    ~Kmer_hash()
    {};

    void savehash(const std::string& outputfile);
    ull loadhash(const std::string& outputfile);
    void preview(const ull kmer_int);
    uint initiate();
    uint finalize();
    bool add(const ull kmer_int);
    uint findhash(ull kmer_int);
    ull totalkmers = 1,totalvalues=0,totaltargets = 0;
    uint32_t totaleles = 0;
    uint32_t bucketcollision = 0;
    unique_ptr<uint[]> Hash1;
    unique_ptr<uint8_t[]>Hash2_size;
    vector<uint> Hash2;
    
    void previewvalue(const ull kmer_int);
    void initiatevalue();
    bool addvalue(const ull kmer_int, const uint16 value);
    bool finalizevalue();
    uint findvalue(const ull kmer_int, vector<uint16>& values);
    
    
    vector<uint16> values;
    
    vector<uint> values_mul_bucksize;
    vector<vector<pair<uint16,uint16>>> values_mul;
};


#endif /* KmerHash_hpp */

