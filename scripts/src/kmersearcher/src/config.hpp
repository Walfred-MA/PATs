#ifndef config_hpp
#define config_hpp

#include <stdio.h>
#include <atomic>
#include <limits>
typedef unsigned int uint;
typedef unsigned __int128 u128;
typedef unsigned long long ull;
typedef unsigned short  uint16;
typedef unsigned char  uint8;
typedef unsigned long long ull_atom;

#ifndef countint
#define countint uint16_t
#endif

constexpr int klen = 31;
constexpr countint MAXCOUNT = std::numeric_limits<countint>::max();
//large_prime will determine the hash space and memory usage: recommend use 536870909, 1073741827, or 2147483647
#define large_prime 0x3FFFFFFF
#define large_size 1573741823

#define primeint30 1073741827
#define primeint31 2147483647
#define primeint22 4194319
#define MAX_UINT30 1073741823
#define MAX_UINT32 0xFFFFFFFF
#define MAX_UINT16 0xFFFF
#define MAX_UINT10 0x4FF
#define MAX_UINT8 0xFF
#define prime10M 9999991
#define MAX_LINE 10000000
#define Comb2( size ) (size + 1) * size / 2
#define MIN( A , B ) ( A <= B ) ? A : B
#define MAX( A , B ) ( A >= B ) ? A : B
#define MAXABS( A , B ) ( abs(A) >= abs(B) ) ? A : B
#define spair std::pair<std::string,std::string>
typedef double FLOAT_T;


#define FLOAT_T double
#define FIXCOL 6

#define genomesize_mean 6320012150.0
#define genomesize_male 6270605410.0
#define genomesize_female 6320012150.0

extern bool optioncorr;

#define errorcutoff1 7
#define errorcutoff2 11


#define sufficient 1000
#define corrstartpoint 0.3

#define windowmerge 15
#define minwindowcutoff 3
#define windowunit 30

#endif /* config_hpp */
