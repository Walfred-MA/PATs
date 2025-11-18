//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License.
//  If you use this code, please cite our work.
//

#ifndef CramReader_hpp
#define CramReader_hpp

#include "htslib/sam.h"
#include <zlib.h>
#include <mutex>
#include <string.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstring>

using namespace std;

class CramReader
{
   
public:
    
    CramReader(const char* inputfile)
    {
        filepath = inputfile;
        kstring = new kstring_t();
        SRread = bam_init1();
        Load();
    };
    ~CramReader()
    {
        ks_free(kstring);
        bam_destroy1(SRread);
        hts_idx_destroy(indexdata);
        sam_hdr_destroy(header);
        hts_itr_multi_destroy(iter);
    }
    
    bool nextLine(std::string &StrLine);
    bool nextLine(uint8_t * &StrLine, size_t &rlen, bam1_t* &SRread);
    bool nextLine_prt(uint8_t * &Str, size_t &rlen);
    void LoadRegion(std::vector<char *>& regions );
    void Load();
    void TotalReads(){};
    
    void Close(){ sam_close(samfile); };

    htsFile *samfile = NULL;
    
    bam1_t * SRread = NULL;
    
    kstring_t* kstring = NULL;
    
    hts_itr_t * iter = NULL;
    
    sam_hdr_t *header = NULL;
    
    hts_idx_t *indexdata = NULL;
    
    std::vector<char *> bedregions;
    std::vector<char*>* workregions=NULL;
    bool useunmap = 1;
    string refpath = "";
    const char *filepath;
    uint16_t filtertag = 0x900;
private:
    
    std::mutex IO_lock;
    bool startumap = 0;
    bool ifindex = 0;
    unsigned long long nreads = 0;
};



#endif

