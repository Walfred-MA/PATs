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
bool ifmask = 1;

extern bool ifweight;
bool ifweight = 0;

extern bool ifhighres;
bool ifhighres = 0;

extern int cutoff;
int cutoff = 1000;

extern int nthreads;
int nthreads = 1;

extern float similar;
float similar = 0.2;

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
        
    int mode = 0, kmer_size = 31;
    
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
        else if (strcmp(Argument, "-c")==0 or strcmp(Argument, "--cut")==0)
        {
            cutoff = atoi(argv[i]);
        }
        else if (strcmp(Argument, "-s")==0 or strcmp(Argument, "--simi")==0)
        {
            similar = std::stof(argv[i]);
        }
        else if (strcmp(Argument, "-n")==0 or strcmp(Argument, "--nthread")==0)
        {
            nthreads = std::stof(argv[i]);
        }
        else if (strcmp(Argument, "-m")==0 or strcmp(Argument, "--mask")==0)
        {
            ifmask = atoi(argv[i]);
        }
        
    }
        
    run<32>(inputfile, kmerfile, outputfile);
    
    return 0;
}
