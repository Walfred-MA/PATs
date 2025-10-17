//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License.
//  If you use this code, please cite our work.
//

#include "KmerHash.hpp"
#define uint unsigned int

/*
static int CompareItem40(const void *a, const void *b)
{
    const uint40 *key = (const uint40 *)a;
    const item40_t *item = (const item40_t *)b;
    if (key->first != item->first.first)
        return (key->first > item->first.first) ? 1 : -1;
    if (key->second != item->first.second)
        return (key->second > item->first.second) ? 1 : -1;
    return 0;
}

static int BSearch(const item40_t *arr, const uint size, const uint40 x)
{
    item40_t *item = bsearch(&x, arr, size, sizeof(item40_t), CompareItem40);
    if (item != NULL) {
        return (int)(item - arr); // The index is the difference between pointers
    }
    return -1;
}
*/
static int Search(const item_t *arr, const uint size, const uint x)
{
    for (uint i = 0 ; i < size; ++i)
    {
        if (arr[i].first == x) return i;
    }
    return -1;
}



uint* Kmer32_hash::add(const ull kmer_int, const uint val)
{
    ull key = kmer_int;
    size_t hash_ = key % modsize;
    uint reminder_ =  (uint) (key / modsize);

    auto &size = key_sizes[hash_];
    item_t *&bucket = items[hash_];

    int loc = Search(bucket, size, reminder_);

    if ( loc == -1 )
    {
        size ++;
        
        bucket = (item_t *) realloc(bucket, sizeof(item_t)*(size));

        bucket[size - 1] = make_pair(reminder_, val);
                    
        return &bucket[size - 1].second;
    }
    
    return &bucket[loc].second;

}


uint* Kmer32_hash::find(const ull kmer_int)
{
    ull key = kmer_int;
    size_t hash_ = key % modsize;
    uint reminder_ =  (uint) (key / modsize);
    
    auto size = key_sizes[hash_];
    item_t *bucket = items[hash_];
        
    int loc = Search(bucket, size, reminder_);
    
    if (loc < 0)
    {
        return NULL;
    }
    else
    {
        return &bucket[loc].second;
    }
}

uint Kmer32_hash::findvalue(const ull kmer_int)
{
    ull key = kmer_int;
    size_t hash_ = key % modsize;
    uint reminder_ =  (uint) (key / modsize);
    
    auto size = key_sizes[hash_];
    item_t *bucket = items[hash_];
        
    int loc = Search(bucket, size, reminder_);
    
    if (loc < 0)
    {
        return 0;
    }
    else
    {
        return bucket[loc].second;
    }
}


bool Kmer32_hash::ifadd(const ull kmer_int, const uint val)
{
    ull key = kmer_int;
    size_t hash_ = key % modsize;
    uint reminder_ =  (uint) (key / modsize);

    auto &size = key_sizes[hash_];
    item_t *&bucket = items[hash_];

    int loc = Search(bucket, size, reminder_);

    if ( loc == -1 )
    {
        size ++;
        
        bucket = (item_t *) realloc(bucket, sizeof(item_t)*(size));

        bucket[size - 1] = make_pair(reminder_, val);
                    
        return 1;
    }
    
    return 0;

}
