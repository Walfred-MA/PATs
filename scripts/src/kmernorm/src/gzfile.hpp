//
//  fasta.hpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 10/17/21.
//  Copyright © 2021 USC_Mark. All rights reserved.
//

#ifndef gzfile_hpp
#define gzfile_hpp

#include <zlib.h>
#include <memory>

#include "readsfile.hpp"

using namespace std;

class gzfile: public readsfile
{
    
public:
    
    gzfile(const char* inputfile):readsfile(inputfile){Load();};
    
    bool nextLine(std::string &StrLine);
    
    void Load();
    
    void Close();

private:
    
    std::unique_ptr<gzFile> fafile = nullptr;
        
};


#endif /* fasta_hpp */
