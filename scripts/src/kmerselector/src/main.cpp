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
#include "KmerFilter.hpp"
#include "Blacklist.hpp"
#include "fasta.hpp"
#include "gzfile.hpp"

extern bool ifmask;
bool ifmask = 0;

extern bool singletarget;
bool singletarget = 0;

BlacklistIntervals blacklist_intervals;

// Default: assemblies absent from a matrix do not contribute exclusion kmers.
// --strict-mod restores the historical all-assemblies/all-matrices behavior.
bool strict_mod = false;

using namespace std;

void shuffle_pair(std::vector<std::string> &vec1, std::vector<std::string> &vec2)
{
    std::random_device rd;
    std::mt19937 g(rd());
    
    std::vector<int> indices(vec1.size());
    for (size_t i = 0; i < indices.size(); i++) indices[i] = i;
    std::shuffle(std::begin(indices), std::end(indices), g);
    
    std::vector<std::string> shuffled_vec1(vec1.size());
    std::vector<std::string> shuffled_vec2(vec2.size());

    

    for (size_t i = 0; i < indices.size(); i++)
    {
        shuffled_vec1[i] = vec1[indices[i]];
        shuffled_vec2[i] = vec2[indices[i]];
    }
    
    vec1 = shuffled_vec1;
    vec2 = shuffled_vec2;
    
}

template <int dictsize>
void run(std::vector<std::string> &inputfiles, std::vector<std::string>& targetfiles, std::vector<std::string>& kmerfiles, std::vector<std::string> &outputfiles, std::vector<std::string> &prefixes, const int kmer_size, const int nthreads, const int mode, const int cutoff)
{

    shuffle_pair(inputfiles, prefixes);

    if (kmerfiles.size() == 0)
        kmerfiles = targetfiles;

    if (mode == 1 || mode == 0)
    {

        kmer_counter<dictsize> counter;

        if (kmerfiles.size())
        {
            counter.read_targets(kmerfiles);
        }
        
        counter.read_files(inputfiles, outputfiles, prefixes, targetfiles, nthreads);
    }
    
    if (mode == 2 || mode == 0)
    {

        kmer_filter<dictsize> filter(kmer_size);
        
        filter.read_files(inputfiles, outputfiles, prefixes, kmerfiles, nthreads);

    }
    
    return;
}


int main(int argc, const char * argv[]) {
    
    
    std::vector<std::string> inputfiles;
    
    std::vector<std::string> prefixes;
    
    std::vector<std::string> targetfiles;
    
    std::vector<std::string> kmerfiles;
    
    std::vector<std::string> outputfiles;
    
    std::vector<std::string> blacklistfiles;

    const char* Argument="";
        
    int mode = 0, kmer_size = 31,  nthreads = 1, cutoff =50;
    
    for (int i = 1; i < argc; i++)
    {
        if (argv[i][0] == '-')
            {
                if (strcmp(argv[i], "--strict-mod") == 0)
                    {
                        strict_mod = true;
                        Argument = "";
                        continue;
                    }
                Argument = argv[i];
            }
        else
            {
                if (strcmp(Argument, "-i") == 0 || strcmp(Argument, "--input") == 0)
                    {
                        if (std::filesystem::is_directory(argv[i]))
                            {
                                for (const auto& entry : std::filesystem::directory_iterator(argv[i]))
                                    {
                                        if (std::filesystem::is_regular_file(entry.path()))
                                            inputfiles.push_back(entry.path().string());
                                    }
                            }
                        else
                            {
                                inputfiles.push_back(argv[i]);
                            }
                    }
                else if (strcmp(Argument, "-I")==0 or strcmp(Argument, "--Inputs")==0)
                {
                    std::ifstream pathfile(argv[i]);
                    std::string line;
                    if(!pathfile)
                    {
                        std::cout<<"Error opening target file"<<std::endl;
                        return -1;
                    }

                    while (std::getline(pathfile, line))
                    {
                        if (line.empty()) continue;

                        std::size_t tabpos = line.find('\t');
                        if (tabpos != std::string::npos)
                        {
                            std::string pref = line.substr(0, tabpos);
                            std::string infile = line.substr(tabpos + 1);
                            prefixes.push_back(pref);
                            inputfiles.push_back(infile);
                        }
                        else
                        {
                            inputfiles.push_back(line);
                        }
                    }

                }
                else if (strcmp(Argument, "-t") == 0 || strcmp(Argument, "--target") == 0)
                    {
                        if (std::filesystem::is_directory(argv[i]))
                            {
                                for (const auto& entry : std::filesystem::directory_iterator(argv[i]))
                                    {
                                        if (std::filesystem::is_regular_file(entry.path()))
                                            targetfiles.push_back(entry.path().string());
                                    }
                            }
                        else
                            {
                                targetfiles.push_back(argv[i]);
                            }
                    }
                else if (strcmp(Argument, "-T") == 0 || strcmp(Argument, "--Targets") == 0)
                    {
                        std::ifstream pathfile(argv[i]);
                        if (!pathfile)
                            {
                                std::cout << "Error opening target list file: " << argv[i] << std::endl;
                                return -1;
                            }
                        std::string line;
                        while (std::getline(pathfile, line))
                            {
                                targetfiles.push_back(line);
                            }
                    }
                else if (strcmp(Argument, "-p") == 0 || strcmp(Argument, "--pref") == 0)
                    {
                        prefixes.push_back(argv[i]);
                    }
                else if (strcmp(Argument, "-P") == 0 || strcmp(Argument, "--Prefs") == 0)
                    {
                        std::ifstream pathfile(argv[i]);
                        if (!pathfile)
                            {
                                std::cout << "Error opening prefixes list file: " << argv[i] << std::endl;
                                return -1;
                            }
                        std::string line;
                        while (std::getline(pathfile, line))
                            {
                                prefixes.push_back(line);
                            }
                    }
                else if (strcmp(Argument, "-o") == 0 || strcmp(Argument, "--output") == 0)
                    {
                        singletarget = 1;
                        outputfiles.push_back(argv[i]);
                    }
                else if (strcmp(Argument, "-O") == 0 || strcmp(Argument, "--Outputs") == 0)
                    {
                        std::ifstream pathfile(argv[i]);
                        if (!pathfile)
                            {
                                std::cout << "Error opening output list file: " << argv[i] << std::endl;
                                return -1;
                            }
                        std::string line;
                        while (std::getline(pathfile, line))
                            {
                                outputfiles.push_back(line);
                            }
                    }
                else if (strcmp(Argument, "-k") == 0 || strcmp(Argument, "--kmer") == 0)
                    {
                        kmerfiles.push_back(argv[i]);
                    }
                else if (strcmp(Argument, "-b") == 0 || strcmp(Argument, "--blacklist") == 0)
                    {
                        if (std::filesystem::is_directory(argv[i]))
                            {
                                for (const auto& entry : std::filesystem::directory_iterator(argv[i]))
                                    {
                                        if (std::filesystem::is_regular_file(entry.path()))
                                            blacklistfiles.push_back(entry.path().string());
                                    }
                            }
                        else
                            {
                                blacklistfiles.push_back(argv[i]);
                            }
                    }
                else if (strcmp(Argument, "-n") == 0 || strcmp(Argument, "--nthreads") == 0)
                    {
                        nthreads = atoi(argv[i]);
                    }
                else if (strcmp(Argument, "-c") == 0 || strcmp(Argument, "--cutoff") == 0)
                    {
                        cutoff = atoi(argv[i]);
                    }
                else if (strcmp(Argument, "-m") == 0 || strcmp(Argument, "--mode") == 0)
                    {
                        mode = atoi(argv[i]);
                    }
                else if (strcmp(Argument, "-M") == 0 || strcmp(Argument, "--mask") == 0)
                    {
                        ifmask = 1;
                        i--; // <- Important: no following value, so don't skip next arg
                    }
                else
                    {
                        std::cout << "Unknown argument: " << Argument << std::endl;
                        return -1;
                    }
                Argument = ""; // Reset argument after consuming
            }
    }
    
    
    
    if (!inputfiles.size()) return 1;
    
    if (prefixes.size() < 2) prefixes.push_back("") ;

    if (blacklistfiles.size())
    {
        blacklist_intervals.load_files(blacklistfiles);
    }

    cout << "cross-matrix filtering: "
         << (strict_mod ? "strict (--strict-mod enabled)" : "off") << endl;
    
    cout <<"start running\n"<<endl;
    run<32>(inputfiles, targetfiles, kmerfiles, outputfiles, prefixes,kmer_size, nthreads,mode,cutoff);
    
    return 0;
}
