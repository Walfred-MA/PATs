//
//  KmerHash.cpp
//  KmerHasher
//
//  Created by walfred on 12/15/24.
//

#include "KmerHash.hpp"

using namespace std;

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

// ---------- 32-bit reversible mix ----------
static inline uint32_t mix32(uint32_t x)
{
    x ^= x >> 16;
    x *= 0x85EBCA6Bu;           // odd => invertible mod 2^32
    x ^= x >> 13;
    x *= 0xC2B2AE35u;           // odd => invertible mod 2^32
    x ^= x >> 16;
    return x;
}

static inline std::pair<uint16_t,uint16_t> pack_u32(uint32_t x)
{
    uint32_t y = mix32(x);
    return { (uint16_t)(y & 0xFFFFu), (uint16_t)(y >> 16) };
}

static inline uint16_t u32(uint32_t x)
{
    uint32_t y = mix32(x);
    return (uint16_t)(y & 0xFFFFu);
}


int find_first_occurrence(const std::vector<std::pair<uint16_t,uint16_t>>& v, uint size,uint16_t key)
{
    uint r = size;
    uint l = 0;        // search in [l, r)
    while (l < r) {
        uint m = l + ((r - l) >> 1); // avoid overflow
        if (v[m].first < key) l = m + 1;
        else                  r = m;   // move left even when equal to find first
    }
    if (l < size && v[l].first == key) return static_cast<int>(l);
    return -1;
}


void Kmer_hash::preview(const ull kmer_int)
{
    uint32_t hash1,hash2;
    hashfunc(kmer_int, hash1, hash2);
    
    uint32_t &thekey = Hash1[hash1];
    
    totalkmers ++;
    thekey ++;
    /*
    if (thekey == 0)  //initiate
    {
        thekey = 1;
    }
    
    else if (thekey == 1) //find collision
    {
        thekey = 2;
        bucketcollision ++;  //record for number of unique collision events
    }
    
    else //N > 1 + 1
    {
        thekey ++ ;  //collision number
    }
    */
    
}

void Kmer_hash::previewvalue(const ull kmer_int)
{
    uint hash = findhash(kmer_int);
    if (hash == 0) return;
    totalvalues ++;
    values[hash]++;
    
}

void Kmer_hash::initiatevalue()
{
    uint32_t hash2index = 1;
    for (uint hash = 0; hash < totalkmers+1 ; ++hash)
    {
        auto& thekey = values[hash];
        
        if (thekey > 1 )
        {
            auto hash_ = u32(hash);
            values_mul_bucksize[hash_] += thekey;
            thekey = MAX_UINT16;
            
        }
        else
        {
            thekey = 0;
        }
        
    }
    
    for (uint hash2 = 0; hash2 < MAX_UINT16+1 ; ++hash2)
    {
        values_mul[hash2].resize(values_mul_bucksize[hash2]);
        values_mul_bucksize[hash2] = 0;
    }
    
    fill_n(values_mul_bucksize.begin(), MAX_UINT16+1, 0);
}


bool Kmer_hash::addvalue(const ull kmer_int, const uint16 value)
{
    uint hash = findhash(kmer_int);
    if (hash == 0) return 0;
    if (values[hash] == MAX_UINT16)
    {
        auto [hash1,hash2]  = pack_u32(hash);
        values_mul[hash1][values_mul_bucksize[hash1]++] = {hash2,value};
    }
    else
    {
        values[hash] = value;
    }
    
    return 0;
}

bool Kmer_hash::finalizevalue()
{
    
    for (uint hash2 = 0; hash2 < MAX_UINT16+1 ; ++hash2)
    {
        std::sort(values_mul[hash2].begin(), values_mul[hash2].end());
    }
    
    return 0;
}

uint Kmer_hash::findvalue(const ull kmer_int, vector<uint16>& results)
{
    uint hash = findhash(kmer_int);
    
    if (hash == 0) return 0;
    
    if (values[hash] != MAX_UINT16)
    {
        results[0] = values[hash];
        return 1;
    }
    else
    {
        auto [hash1,hash2]  = pack_u32(hash);
        
        auto index = find_first_occurrence(values_mul[hash1], values_mul_bucksize[hash1], hash2);
        
        if (index == -1) return 0;
        
        int values_i = 0;
        for (;index < values_mul_bucksize[hash1]; ++index)
        {
            auto& thevalue = values_mul[hash1][index];
            
            if (thevalue.first == hash2)
            {
                results[values_i++] = thevalue.second;
            }
            else
            {
                break;
            }
            
        }
        return values_i;
    }
    
    return 0;
}


uint Kmer_hash::initiate()
{
    //flat all collisions
    Hash2.resize(totalkmers + 10, 0);
    values.resize(totalkmers + 10, 0);
    
    
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


void Kmer_hash::savehash(const std::string& outputfile)
{
    std::string tmpfile = outputfile + ".tmp";

    std::ofstream ofs(tmpfile, std::ios::binary);
    if (!ofs) {
        std::cerr << "Error: Cannot open temp file: " << tmpfile << std::endl;
        return;
    }

    uint64_t size = Hash2.size();
    uint64_t hash1_size = large_prime;
    
    // Sizes
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));
    ofs.write(reinterpret_cast<const char*>(&hash1_size), sizeof(hash1_size));
    
    // Save external totalkmers
    ofs.write(reinterpret_cast<const char*>(&totalkmers), sizeof(totalkmers));
    ofs.write(reinterpret_cast<const char*>(&totaltargets), sizeof(totaltargets));

    // Save internal stats
    ofs.write(reinterpret_cast<const char*>(&bucketcollision), sizeof(bucketcollision));

    // Data
    ofs.write(reinterpret_cast<const char*>(Hash1.get()), sizeof(uint) * hash1_size);
    ofs.write(reinterpret_cast<const char*>(Hash2.data()), sizeof(uint) * size);
    ofs.write(reinterpret_cast<const char*>(Hash2_size.get()), sizeof(uint8_t) * hash1_size);
    
    ofs.write(reinterpret_cast<const char*>(values.data()), sizeof(uint16) * size);
    ofs.write(reinterpret_cast<const char*>(values_mul_bucksize.data()), sizeof(uint) * (MAX_UINT16+1));
    
    for (uint64_t b = 0; b < MAX_UINT16+1; ++b)
    {
        const auto& bucket = values_mul[b];
        uint64_t n = static_cast<uint64_t>(bucket.size());
        ofs.write(reinterpret_cast<const char*>(&n), sizeof(n));
        if (n)
        {
            ofs.write(reinterpret_cast<const char*>(bucket.data()),
                      n * sizeof(std::pair<uint16_t,uint16_t>));
        }
    }
    
    

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

    uint64_t size = 0;
    uint64_t hash1_size = 0;
    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
    ifs.read(reinterpret_cast<char*>(&hash1_size), sizeof(hash1_size));

    // Load external totalkmers
    ifs.read(reinterpret_cast<char*>(&totalkmers), sizeof(totalkmers));
    ifs.read(reinterpret_cast<char*>(&totaltargets), sizeof(totaltargets));
    // Load internal stats
    ifs.read(reinterpret_cast<char*>(&bucketcollision), sizeof(bucketcollision));

    if (hash1_size != large_prime) {
        std::cerr << "[loadhash] ERROR: matrix .bin file mismatch." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    Hash1 = std::make_unique<uint[]>(hash1_size);
    Hash2.resize(size);
    Hash2_size = std::make_unique<uint8_t[]>(hash1_size);
    values.resize(size);
    
    ifs.read(reinterpret_cast<char*>(Hash1.get()), sizeof(uint) * hash1_size);
    ifs.read(reinterpret_cast<char*>(Hash2.data()), sizeof(uint) * size);
    ifs.read(reinterpret_cast<char*>(Hash2_size.get()), sizeof(uint8_t) * hash1_size);
    
    ifs.read(reinterpret_cast<char*>(values.data()), sizeof(uint16) * size);
    values_mul_bucksize.resize(MAX_UINT16+1, 0);
    ifs.read(reinterpret_cast<char*>(values_mul_bucksize.data()), sizeof(uint) * (MAX_UINT16+1) );
    
    uint64_t nbuckets = MAX_UINT16+1;
    values_mul.clear();
    values_mul.resize(nbuckets);

    for (uint64_t b = 0; b < nbuckets; ++b)
    {
        uint64_t n = 0;
        ifs.read(reinterpret_cast<char*>(&n), sizeof(n));
        auto& bucket = values_mul[b];
        bucket.resize(n);
        if (n)
        {
            ifs.read(reinterpret_cast<char*>(bucket.data()),
                     n * sizeof(std::pair<uint16_t,uint16_t>));
        }
    }
    
    
    ifs.close();

    std::cout << "[loadhash] Successfully loaded "
              << "Hash1[" << hash1_size << "], Hash2[" << size << "], "
              << "totalkmers = " << totalkmers << std::endl;

    return totalkmers;
}

