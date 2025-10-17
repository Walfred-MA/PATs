//
//  fasta.cpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 10/17/21.
//  Copyright © 2021 USC_Mark. All rights reserved.
//

#include "fasta.hpp"
//reading reference
void fasta::Load()
{
    if (filepath.length()<2) return;
    
    fafile.open(filepath, ios_base::in);
    
    if (!fafile) {
        
        std::cerr << "ERROR: Could not open " << filepath << " for reading.\n" << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
    
    

}


bool fasta::nextLine(std::string &StrLine)
{
    if (std::getline(fafile, StrLine)) {
        return true;
    }

    // If we reached the end of fafile, check for extension file
    if (fafile.eof()) {
        
        // Only attempt once — close main file and open the extension if it exists
        if (std::filesystem::exists(filepath + "._app_")) {
            filepath = filepath + "._app_";
            fafile.close();
            fafile.open(filepath, std::ios_base::in);

            if (!fafile) {
                std::cerr << "ERROR: Could not open " << filepath << " for reading.\n";
                std::_Exit(EXIT_FAILURE);
            }

            // Try reading again from the new file
            return static_cast<bool>(std::getline(fafile, StrLine));
        }
    }

    return false;
}


void fasta::Close()
{
    
    fafile.close();
    
}
