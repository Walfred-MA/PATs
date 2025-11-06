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

static constexpr uint32_t MASK30 = (1u << 30) - 1;

inline uint32_t rotr32(uint32_t x, unsigned r)
{
    return (x >> r) | (x << ((32 - r) & 31));
}

inline uint32_t g32to30(uint32_t v)
{
    uint32_t t = rotr32(v, 7) ^ rotr32(v, 13) ^ rotr32(v, 19);
    return t & MASK30;
}

inline void hashfunc(ull kmer_int, uint32_t &hash1,uint32_t &hash2)
{
    hash1 = kmer_int & MASK30;
    hash2 = (uint32_t)(kmer_int >> 30);
    hash1 = hash1 ^ g32to30(hash2);
}

inline ull unhashfunc(uint32_t hash1, uint32_t hash2)
{
    uint32_t hi = hash2;

    uint32_t lo = hash1 ^ g32to30(hi);

    ull kmer_int = (static_cast<ull>(hi) << 30) | static_cast<ull>(lo);
    return kmer_int;
}


static int Search(const item_t *arr, const uint size, const uint x)
{
    for (uint i = 0 ; i < size; ++i)
    {
        //if (arr[i].first == x) return i;
        if (arr[i] == x) return i;
    }
    return -1;
}



uint* Kmer32_hash::add(const ull kmer_int, const uint val)
{
    ull key = kmer_int;
    uint hash_, reminder_;
    hashfunc(kmer_int, hash_, reminder_);

    auto &size = key_sizes[hash_];
    item_t *&bucket = items[hash_];

    int loc = Search(bucket, size, reminder_);

    if ( loc == -1 )
    {
        size ++;
        
        bucket = (item_t *) realloc(bucket, sizeof(item_t)*(size));

        //bucket[size - 1] = make_pair(reminder_, val);
        bucket[size - 1] = reminder_;
        
        return &bucket[size - 1];
    }
    
    return &bucket[loc];

}

uint Kmer32_hash::initiate()
{
    uint32_t hash2index = 1;
    for (size_t hash1 = 0; hash1 < large_prime ; ++hash1)
    {
        key_starts[hash1] = hash2index;
        auto &size = key_sizes[hash1];
        hash2index += size;
    }
    
    return hash2index;
}


uint* Kmer32_hash::find(const ull kmer_int)
{
    ull key = kmer_int;
    uint hash_, reminder_;
    hashfunc(kmer_int, hash_, reminder_);
    
    auto size = key_sizes[hash_];
    item_t *bucket = items[hash_];
        
    int loc = Search(bucket, size, reminder_);
    
    if (loc < 0)
    {
        return NULL;
    }
    else
    {
        return &bucket[loc];
    }
}

uint Kmer32_hash::findhash(const ull kmer_int)
{
    ull key = kmer_int;
    uint hash_, reminder_;
    hashfunc(kmer_int, hash_, reminder_);
    
    auto size = key_sizes[hash_];
    item_t *bucket = items[hash_];
        
    int loc = Search(bucket, size, reminder_);
    
    if (loc < 0)
    {
        return 0;
    }
    else
    {
        return key_starts[hash_] + loc;
    }
}


bool Kmer32_hash::ifadd(const ull kmer_int, const uint val)
{
    ull key = kmer_int;
    uint hash_, reminder_;
    hashfunc(kmer_int, hash_, reminder_);

    auto &size = key_sizes[hash_];
    item_t *&bucket = items[hash_];

    int loc = Search(bucket, size, reminder_);

    if ( loc == -1 )
    {
        size ++;
        
        bucket = (item_t *) realloc(bucket, sizeof(item_t)*(size));

        bucket[size - 1] = reminder_;
                    
        return 1;
    }
    
    return 0;

}





void Kmer_hash::preview(const ull kmer_int)
{
    uint32_t hash1,hash2;
    hashfunc(kmer_int, hash1, hash2);
    
    uint32_t &thekey = Hash1[hash1];
    
    totaleles ++;
    thekey ++;
    
}


uint Kmer_hash::initiate()
{
    //flat all collisions
    Hash2.resize(totaleles + 10, 0);
    
    uint32_t hash2index = 1;
    for (size_t hash1 = 0; hash1 < large_prime ; ++hash1)
    {
        uint32_t &thekey = Hash1[hash1];
        thekey = (uint8_t) MIN( MAX_UINT8 , thekey);
        auto &counter = Hash2_size[hash1];
        counter = thekey;
        
        if (thekey == 0 )
        {
            continue;
        }
        else
        {
            thekey = hash2index;
            hash2index += counter;           //get next local pointer
        }
    }
    
    return hash2index;
}

bool Kmer_hash::iterate(uint &hash1, uint8_t &hash2index, ull &kmer)
{
    while (hash1 < large_prime) {

        uint posi    = Hash1[hash1];
        uint8_t counter = Hash2_size[hash1];

        if (hash2index >= counter) {
            hash2index = 0;
            ++hash1;
            continue;
        }

        uint hash2_ = Hash2[posi + hash2index];
        ++hash2index;

        kmer = unhashfunc(hash1, hash2_);
        return true;
    }

    return false;
}

bool Kmer_hash::add(const ull kmer_int)
{
    uint32_t hash1,hash2;
    hashfunc(kmer_int, hash1, hash2);
    
    uint32_t& thekey = Hash1[hash1];
    
    uint32_t redirect = thekey;
    for (uint32_t index = redirect;  index < redirect + Hash2_size[hash1]; ++index)
    {
        if (Hash2[index] == 0)
        {
            Hash2[index] = hash2 ;    //unintiated
            return 1;
        }
        else if (Hash2[index] == hash2)
        {
            return 0;
        }
    }
    
    return 0;
}

uint Kmer_hash::finalize()
{
    uint32_t sortcounter = 0;
    for (size_t hash1 = 0; hash1 < large_prime ; ++hash1)
    {
        const uint8_t numelement = Hash2_size[hash1];
        
        if (numelement >= 8 )
        {
            sortcounter ++;
            uint32_t redirect  = Hash1[hash1];
            uint32_t* base = Hash2.data() + redirect;
            std::sort(base, base + numelement);//get next local pointer
        }
    }
    
    return sortcounter;
}


uint Kmer_hash::findhash(const ull kmer_int)
{
    uint32_t hash1,hash2;
    hashfunc(kmer_int, hash1, hash2);
    
    uint32_t thenum = Hash2_size[hash1];
    if (thenum == 0)
    {
        return 0;
    }
    else if (thenum < 8)
    {
        uint32_t redirect = Hash1[hash1];
        for (uint32_t index = redirect;  index < redirect + thenum; ++index)
        {
            if (Hash2[index] == hash2) return index;
        }
    }
    else             //check redirect
    {
        uint32_t redirect = Hash1[hash1];
        uint32_t lo = redirect, hi = redirect + thenum;
        while (lo < hi)
        {
            uint32_t mid = lo + ((hi - lo) >> 1);
            uint32_t v = Hash2[mid];
            if (v < hash2) lo = mid + 1;
            else hi = mid;
        }
        if (lo < redirect + thenum && Hash2[lo] == hash2) return lo;
        return 0;
    }
    
    return 0;
}


void Kmer_hash::savehash(const std::string& outputfile, ull totalkmers)
{
    std::string tmpfile = outputfile + ".tmp";

    std::ofstream ofs(tmpfile, std::ios::binary);
    if (!ofs) {
        std::cerr << "Error: Cannot open temp file: " << tmpfile << std::endl;
        return;
    }

    uint64_t size = Hash2.size();
    uint64_t hash1_size = large_prime;

    // Save external totalkmers
    ofs.write(reinterpret_cast<const char*>(&totalkmers), sizeof(totalkmers));

    // Save internal stats
    ofs.write(reinterpret_cast<const char*>(&totaleles), sizeof(totaleles));
    ofs.write(reinterpret_cast<const char*>(&bucketcollision), sizeof(bucketcollision));

    // Sizes
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
    ofs.write(reinterpret_cast<const char*>(&hash1_size), sizeof(hash1_size));

    // Data
    ofs.write(reinterpret_cast<const char*>(Hash1.get()), sizeof(uint) * hash1_size);
    ofs.write(reinterpret_cast<const char*>(Hash2.data()), sizeof(uint) * size);
    ofs.write(reinterpret_cast<const char*>(Hash2_size.get()), sizeof(uint8_t) * hash1_size);

    ofs.close();
    std::rename(tmpfile.c_str(), outputfile.c_str());
}

ull Kmer_hash::loadhash(const std::string& inputfile)
{
    std::ifstream ifs(inputfile, std::ios::binary);
    if (!ifs) {
        std::cerr << "Error: Cannot open file for reading: " << inputfile << std::endl;
        return 0;
    }

    ull totalkmers = 0;

    // Load external totalkmers
    ifs.read(reinterpret_cast<char*>(&totalkmers), sizeof(totalkmers));

    // Load internal stats
    ifs.read(reinterpret_cast<char*>(&totaleles), sizeof(totaleles));
    ifs.read(reinterpret_cast<char*>(&bucketcollision), sizeof(bucketcollision));

    uint64_t size = 0;
    uint64_t hash1_size = 0;
    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
    ifs.read(reinterpret_cast<char*>(&hash1_size), sizeof(hash1_size));

    if (hash1_size != large_prime) {
        std::cerr << "[loadhash] ERROR: matrix .bin file mismatch." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    Hash1 = std::make_unique<uint[]>(hash1_size);
    Hash2.resize(size);
    Hash2_size = std::make_unique<uint8_t[]>(hash1_size);

    ifs.read(reinterpret_cast<char*>(Hash1.get()), sizeof(uint) * hash1_size);
    ifs.read(reinterpret_cast<char*>(Hash2.data()), sizeof(uint) * size);
    ifs.read(reinterpret_cast<char*>(Hash2_size.get()), sizeof(uint8_t) * hash1_size);

    ifs.close();

    std::cout << "[loadhash] Successfully loaded "
              << "Hash1[" << hash1_size << "], Hash2[" << size << "], "
              << "totalkmers = " << totalkmers << std::endl;

    return totalkmers;
}
