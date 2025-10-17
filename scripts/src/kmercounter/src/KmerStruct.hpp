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

using namespace std;

typedef unsigned int uint;
typedef unsigned __int128 u128;
typedef unsigned long long ull;
typedef unsigned short  uint16;
typedef unsigned char  uint8;

extern bool ifmask;

#define large_prime 2147483647


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
    
    if (base < ':' ) return 0;
    
    converted = base - ':';
    
    return 1;
}




#endif /* struct_hpp */
