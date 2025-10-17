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
#include <functional>
#include <cstdint>
#include <cstring>

#define MIN( A , B ) ( A <= B ) ? A : B
#define MAX( A , B ) ( A >= B ) ? A : B
#define large_prime 2147483647
#define MAX_UINT32 16777215
#define item_t std::pair<uint,uint>
#define uint unsigned int
#define uint16 unsigned short

typedef unsigned long long ull;
typedef unsigned char  uint8;

using namespace std;



class Kmer32_hash
{
public:
    Kmer32_hash(int size): modsize(MAX(large_prime,size))
    {
         memset(key_sizes, 0 , modsize);
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
    }
    
    uint* add(const ull kmer_int, const uint index);
    bool ifadd(const ull kmer_int, const uint index);
    uint* find(const ull kmer_int);
    uint findvalue(const ull kmer_int);
    const size_t modsize;

private:
    item_t** items = new item_t*[modsize];
    uint16* key_sizes = new uint16[modsize];
};




#endif /* KmerHash_hpp */
