# cython: boundscheck=False, wraparound=False, cdivision=True
import numpy as np
cimport numpy as np
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdio cimport FILE, fopen, fread, fclose
from libc.math cimport round,log
from libc.stdlib cimport malloc, free

cdef uint64_t large_prime = 1073741827

cdef struct Stats:
    uint16_t mean_assem
    uint16_t valid
    uint16_t mean_ratio
    uint16_t var_ratio
    uint8_t overbias
    uint8_t lowerbias
    
cdef class KmerDB:
    cdef uint32_t* Hash1
    cdef uint32_t* Hash2
    cdef uint8_t* Hash2_size
    cdef size_t size1, size2
    
    # 4-field output: mean_ratio, var_ratio, bias, valid
    cdef uint8_t* mean_ratios
    cdef uint8_t* var_ratios
    cdef uint8_t* biases
    cdef uint8_t* valids
    cdef size_t num_kmers
    
    def __cinit__(self):
        self.Hash1 = NULL
        self.Hash2 = NULL
        self.Hash2_size = NULL
        self.mean_ratios = NULL
        self.var_ratios = NULL
        self.biases = NULL
        self.valids = NULL
        
    def loadhash(self, str hashfile):
        # Load hash
        cdef FILE* fp = fopen(hashfile.encode('utf-8'), "rb")
        if not fp:
            raise IOError("Cannot open hash file")
            
        cdef uint64_t totalkmers, totaleles, bucketcollision
        fread(&totalkmers, sizeof(uint64_t), 1, fp)
        fread(&totaleles, sizeof(uint64_t), 1, fp)
        fread(&bucketcollision, sizeof(uint64_t), 1, fp)
        fread(&self.size2, sizeof(uint64_t), 1, fp)
        fread(&self.size1, sizeof(uint64_t), 1, fp)
        
        self.Hash1 = <uint32_t*> malloc(self.size1 * sizeof(uint32_t))
        self.Hash2 = <uint32_t*> malloc(self.size2 * sizeof(uint32_t))
        self.Hash2_size = <uint8_t*> malloc(self.size2 * sizeof(uint8_t))
        
        if not self.Hash1 or not self.Hash2 or not self.Hash2_size:
            raise MemoryError("Hash table allocation failed")
            
        fread(self.Hash1, sizeof(uint32_t), self.size1, fp)
        fread(self.Hash2, sizeof(uint32_t), self.size2, fp)
        fread(self.Hash2_size, sizeof(uint8_t), self.size2, fp)
        fclose(fp)
        
    def loadstat(self, str statsfile, int samplesize = 300):
        # Load Stats
        fp = fopen(statsfile.encode('utf-8'), "rb")
        if not fp:
            raise IOError("Cannot open stats file")
            
        self.num_kmers = self.size2
        self.mean_ratios = <uint8_t*> malloc(self.num_kmers * sizeof(uint16_t))
        self.var_ratios = <uint8_t*> malloc(self.num_kmers * sizeof(uint8_t))
        self.biases = <uint8_t*> malloc(self.num_kmers * sizeof(uint8_t))
        self.valids = <uint8_t*> malloc(self.num_kmers * sizeof(uint8_t))
        
        if not self.mean_ratios or not self.var_ratios or not self.biases or not self.valids:
            raise MemoryError("Stats buffer allocation failed")
        cdef Stats s
        cdef size_t i

        for i in range(self.num_kmers):
            if fread(&s, sizeof(Stats), 1, fp) != 1:
                break
            
            self.valids[i] = <uint8_t> min(255, s.valid * 100/samplesize) #0-100
            self.mean_ratios[i] = <uint8_t> ( max(50,min( s.mean_ratio, 150))-50 if s.valid >= 16 else 50 )#ratio adjust, need to be balance on both sides of 100
            self.var_ratios[i] = <uint8_t> (max(10, int(round((5000.0 / max(500, <int>s.var_ratio)) ** 2))) if s.valid >= 16 else 100)
            self.biases[i] = <uint8_t> (min(90, ((<int>s.overbias + <int>s.lowerbias) * 100 // samplesize)) if s.valid >= 16 else 0)

        fclose(fp)
    
    def load(self, str file, int samplesize = 300):
        
        self.loadhash(file+".hash")
        self.loadstat(file+".stat", samplesize)
    
    cpdef tuple find(self, uint64_t kmer):
        cdef uint32_t hash1 = <uint32_t> (kmer % large_prime)
        cdef uint32_t hash2 = <uint32_t> (kmer // large_prime)
        cdef uint32_t redirect = self.Hash1[hash1]
        cdef uint8_t len = self.Hash2_size[redirect]
        cdef uint32_t idx
        for idx in range(redirect, redirect + len):
            if self.Hash2[idx] == hash2:
                return ( self.valids[idx], self.biases[idx], self.mean_ratios[idx], self.var_ratios[idx])
        return None  # Not found
    
    def __dealloc__(self):
        if self.Hash1: free(self.Hash1)
        if self.Hash2: free(self.Hash2)
        if self.Hash2_size: free(self.Hash2_size)
        if self.mean_ratios: free(self.mean_ratios)
        if self.var_ratios: free(self.var_ratios)
        if self.biases: free(self.biases)
        if self.valids: free(self.valids)

def create_kmerdb():
    return KmerDB()
