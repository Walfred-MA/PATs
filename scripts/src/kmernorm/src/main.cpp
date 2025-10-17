//
//  main.cpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 10/13/21.
//  Copyright © 2021 USC_MarkLab. All rights reserved.
//

#include <iostream>
#include <string>
#include <filesystem>
#include <algorithm>
#include <random>

#include "KmerCounter.hpp"
#include "fasta.hpp"

extern bool ifmask;
bool ifmask = 0;

extern bool ifweight;
bool ifweight = 1;

extern bool ifhighres;
bool ifhighres = 1;

template <int dictsize>
void run(std::string &inputfile, std::string &kmerfile,std::string &outputfile)
{
    kmer_counter <32>counter;
    
    counter.getnorm(inputfile, kmerfile, outputfile);
    
    return;
}


int main(int argc, const char * argv[]) {
    
    
    std::string inputfile = "";
    std::string outputfile = "";
    std::string kmerfile = "";
    
    const char* Argument = "";
        
    int mode = 0, kmer_size = 31,  nthreads = 1, cutoff =50;
    
    for (int i = 1; i < argc ; i++)
    {
        if (argv[i][0] == '-')
        {
            Argument = argv[i];
        }
        else if (strcmp(Argument, "-i")==0 or strcmp(Argument, "--input")==0)
        {
            inputfile = argv[i];
        }
        
        else if (strcmp(Argument, "-o")==0 or strcmp(Argument, "--output")==0)
        {
            outputfile = argv[i];
        }
        else if (strcmp(Argument, "-k")==0 or strcmp(Argument, "--kmer")==0)
        {
            kmerfile = argv[i];
        }
        else if (strcmp(Argument, "-w")==0 or strcmp(Argument, "--weight")==0)
        {
            ifweight = atoi(argv[i]);
        }
        else if (strcmp(Argument, "-h")==0 or strcmp(Argument, "--highres")==0)
        {
            ifhighres = atoi(argv[i]);
        }
        else if (strcmp(Argument, "-m")==0 or strcmp(Argument, "--mask")==0)
        {
            ifmask = atoi(argv[i]);
        }
        
    }
        
    run<32>(inputfile, kmerfile, outputfile);
    
    return 0;
}
