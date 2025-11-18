//
//  fasta.hpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 10/17/21.
//  Copyright © 2021 USC_Mark. All rights reserved.
//

#ifndef fasta_hpp
#define fasta_hpp

#include "readsfile.hpp"

using namespace std;

class fasta: public readsfile
{
    
public:
    
    fasta(const char* inputfile):readsfile(inputfile){Load();};
    
    bool nextLine(std::string &StrLine);
    
    void Load();
    
    void Close();

    void Reset();
private:
    
    std::fstream fafile;
        
};


#endif /* fasta_hpp */
