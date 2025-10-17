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

#include "KmerTree.hpp"
#include "fasta.hpp"

extern bool ifmask;
bool ifmask = 1;

extern bool ifweight;
bool ifweight = 0;

extern bool ifhighres;
bool ifhighres = 0;

extern int nthreads;
int nthreads = 1;

void run(std::string &inputfile, std::string &normfile, std::string &outputfile)
{
    kmer_tree tree;
    
    tree.run(inputfile, normfile, outputfile);
    
    return;
}


int main(int argc, const char * argv[]) {
    
    
    std::string inputfile = "";
    std::string outputfile = "";
    std::string kmerfile = "";
    std::string normfile = "";
    const char* Argument = "";
        
    int mode = 0, kmer_size = 31, cutoff =50;
    
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
        else if (strcmp(Argument, "-n")==0 or strcmp(Argument, "--norm")==0)
        {
            normfile = argv[i];
        }
        else if (strcmp(Argument, "-t")==0 or strcmp(Argument, "--nthreads")==0)
        {
            nthreads = atoi(argv[i]);
        }
        
    }
        
    run(inputfile, normfile, outputfile);
    
    return 0;
}
