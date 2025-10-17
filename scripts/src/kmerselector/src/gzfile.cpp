//
//  fasta.cpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 10/17/21.
//  Copyright © 2021 USC_Mark. All rights reserved.
//

#include "gzfile.hpp"

using namespace std;
//reading reference
void gzfile::Load()
{
    if (filepath.length()<2) return;
    
    fafile = std::make_unique<gzFile>(gzopen(filepath.c_str(), "rb"));
    if (!fafile) {
        
        std::cerr << "ERROR: Could not open " << filepath << " for reading.\n" << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
    
    

}

bool gzfile::nextLine(std::string &StrLine)
{
    return(bool)gzgets(*fafile, &StrLine[0], StrLine.size());
}

void gzfile::Close()
{
    gzclose(*fafile);
}
