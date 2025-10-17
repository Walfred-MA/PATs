//
//  KmerHash.cpp
//  KmerHasher
//
//  Created by walfred on 12/15/24.
//

#include "KmerHash.hpp"

using namespace std;

void Kmer_hash::preview(const ull kmer_int)
{
    uint32_t hash1 = (uint32_t) kmer_int % large_prime ;
    uint32_t hash2 = (uint32_t) kmer_int / large_prime;
    uint32_t &thekey = Hash1[hash1];
        
    if (thekey == MAX_UINT32)  //initiate
    {
        thekey = hash2;
    }
    
    else if ((thekey & 0x80000000) == 0) //find collision
    {
        thekey = 0x80000002;
        bucketcollision ++;  //record for number of unique collision events
        totalcollision += 2;  //record for number of total collision
    }
    
    else //N > 1 + 1
    {
        thekey ++ ;  //collision number
        totalcollision ++;
    }
    
    
}


void Kmer_hash::initiate()
{
    //flat all collisions
    Hash2.resize(bucketcollision + 10, MAX_UINT32);
    Hash2_size.resize(bucketcollision/2 + 10);
    
    uint32_t hash2index = 0;
    for (size_t hash1 = 0; hash1 < large_prime ; ++hash1)
    {
        uint32_t &thekey = Hash1[hash1];
        uint32_t numelement = 0;
        if (thekey == MAX_UINT32 )
        {
            continue;
        }
        else if (thekey > 0x80000000)
        {
            numelement = thekey - 0x80000000;
            thekey = hash2index + 0x80000000;
            Hash2_size[hash2index/2] = (uint8_t) MIN( MAX_UINT8 , numelement);     //and record the pointer size, condense index by 2
            hash2index += (uint8_t) MIN( MAX_UINT8 , numelement);           //get next local pointer
        }
    }
}

bool Kmer_hash::add(const ull kmer_int)
{
    uint32_t hash1 = (uint32_t) (kmer_int % large_prime) ;
    uint32_t hash2 = (uint32_t) (kmer_int / large_prime);
    uint32_t& thekey = Hash1[hash1];
    
    if (thekey != MAX_UINT32)      //redirected
    {
        uint32_t redirect = thekey - 0x80000000;
        for (uint32_t index = redirect;  index < redirect + Hash2_size[redirect/2]; ++index)
        {
            if (Hash2[index] == MAX_UINT32)
            {
                Hash2[index] = hash2 ;    //unintiated
                return 1;
            }
            else if (Hash2[index] == thekey)
            {
                return 0;
            }
        }
    }
    else
    {
        thekey = hash2;
        return 1;
    }
    
    return 0;
}

uint *Kmer_hash::find(const ull kmer_int)
{
    uint32_t hash1 = (uint32_t) (kmer_int % large_prime) ;
    uint32_t hash2 = (uint32_t) (kmer_int / large_prime);
    uint32_t thekey = Hash1[hash1];
    
    uint *posi = NULL;
    if (thekey == MAX_UINT32)
    {
        return NULL;
    }
        
    else if (thekey >= 0x80000000)             //check redirect
    {
        uint32_t redirect = thekey - 0x80000000;
        
        for (uint32_t index = redirect;  index < redirect + Hash2_size[redirect]; ++index)
        {
            if (Hash2[index] == hash2) return &Hash2[index];
        }
    }
    
    else
    {
        return &Hash1[hash1] ;
    }
    
    return posi;
}
