//
//  struct.hpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 2/2/23.
//  Copyright © 2023 USC_Mark. All rights reserved.
//

#ifndef KmerStruct_hpp
#define KmerStruct_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstring>
#include <unordered_set>
#include <string>
#include <thread>
#include <atomic>
#include <zlib.h>
#include "config.hpp"

using namespace std;

typedef unsigned int uint;
typedef unsigned __int128 u128;
typedef unsigned long long ull;
typedef unsigned short  uint16;
typedef unsigned char  uint8;

extern bool ifmask;

using kmer32_dict_nt = std::unordered_map<ull, uint8>;
using kmer32_set_nt = std::unordered_map<ull,bool>;

struct hash_128 {
    
    size_t operator()(const u128& num128) const
    {
        auto hash1 = hash<ull>{}((ull)(num128%large_prime));
        return hash1;
    }
};

static std::string toString128(u128 num)
{
    std::string str;
    do {
        int digit = num % 2;
        str = std::to_string(digit) + str;
        num = num / 2;
    } while (num != 0);
    return str;
}

static std::string toString128(ull num)
{
    std::string str;
    do {
        int digit = num % 2;
        str = std::to_string(digit) + str;
        num = num / 2;
    } while (num != 0);
    return str;
}

static bool base_to_int(char base, int &converted)
{
    
    if (base >= 'a' )
    {
        base -= 32;
    }

    switch (base)
    {
        case 'A':
            converted=0b00;
            break;
        case 'T':
            converted=0b11;
            break;
        case 'C':
            converted=0b01;
            break;
        case 'G':
            converted=0b10;
            break;
        default:
            return 0;
    }
    
    return 1;
}

static bool base_to_int_64(char base, int &converted)
{
    
    if (base < '0' ) return 0;
    
    converted = base - '0';
    
    return 1;
}

static std::string kmer_int_totext(ull kmer, const int nbase = 64)
{
    int code_bit = 3;
    
    const int digit = (int)ceil(1.0*klen/code_bit);
    
    char kmer_seq[digit+1];
    kmer_seq[digit] = '\0';
    
    for (int index = digit-1; index >= 0 ; --index)
    {
        kmer_seq[index] = '0' + kmer%nbase;
        kmer /= nbase;
    }
    
    return string(kmer_seq);
    
}

static std::string kmer_int_toatcg(ull kmer)
{
    int code_bit = 1;
    
    const int digit = (int)ceil(1.0*klen/code_bit);
    
    char kmer_seq[digit+1];
    kmer_seq[digit] = '\0';
    
    for (int index = digit-1; index >= 0 ; --index)
    {
        kmer_seq[index] = "ACGT"[kmer%4];
        kmer /= 4;
    }
    
    return string(kmer_seq);
    
}


static void write_cache(const char* outputfile, const std::unordered_set<ull>& allkmers)
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


static void read_cache_to_map(const char* inputfile, kmer32_dict_nt& dict)
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

static void read_cache_to_map(const char* inputfile, kmer32_set_nt& dict)
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


#endif /* struct_hpp */
