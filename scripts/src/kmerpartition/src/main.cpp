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
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <sstream>

#include "KmerAssigner.hpp"
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

extern bool ifbed;
bool ifbed = 0;

extern std::string queryfilepath;
std::string queryfilepath = "";

extern bool ifmultipleKmer;
bool ifmultipleKmer = 0;

static void print_usage(const char *program)
{
    std::cerr
        << "Usage:\n"
        << "  " << program << " -i INPUT -k KMER -o OUTPUT [options]\n"
        << "  " << program << " -i INPUT.bed -q QUERY.fa -k KMER -o OUTPUT.bed --predefined 1 [options]\n"
        << "  " << program << " -i INPUT.bed -q QUERY.fa -K KMERFILELIST -o OUTPUT.bed --predefined 1 [options]\n\n"
        << "Required/common arguments:\n"
        << "  -i, --input       Input FASTA or BED file. BED uses 0-based half-open coordinates.\n"
        << "  -o, --output      Output prefix for FASTA mode, or output BED path for BED mode.\n"
        << "  -k, --kmer        Single kmer file/list.\n"
        << "  -K, --kmer-list   Multi-matrix kmer-file list. Each line is a kmer filepath; the matrix\n"
        << "                    prefix is the true filename text before the first '_'. Only valid with\n"
        << "                    BED input and --predefined 1. Different matrices run in parallel.\n"
        << "  -q, --query       Required for BED input. If it ends with .fa/.fasta, it is one indexed\n"
        << "                    FASTA and is loaded into RAM. Otherwise it is a path file; column 2\n"
        << "                    is used as each FASTA path and intervals are fetched through .fai.\n\n"
        << "Modes:\n"
        << "  default           Build partitions from a target kmer file, then write FASTA/kmer outputs.\n"
        << "  --predefined 1    Assign sequences/intervals to predefined groups from >partition/kmer records.\n"
        << "  BED + -k          Reads BED intervals as virtual FASTA via -q and writes partition score to BED col 5.\n"
        << "  BED + -K          Multi-matrix predefined BED mode. BED rows are matched to matrix prefixes from -K,\n"
        << "                    assigned with the corresponding kmer file, and BED col 5 receives group number(s).\n\n"
        << "Other options:\n"
        << "  -t, --targets     Comma-separated target names for one-partition target mode. Not compatible\n"
        << "                    with --predefined. FASTA names use header text before tab/space; BED uses col 4.\n"
        << "  -c, --cut         Kmer count cutoff (default 1000).\n"
        << "  -s, --simi        Similarity threshold for non-predefined mode.\n"
        << "  -n, --nthread     Number of threads/matrix workers.\n"
        << "  -m, --mask        Ignore lowercase target kmers when enabled.\n"
        << "  -h, --help        Show this help.\n";
}

static std::vector<std::string> read_kmerfile_list(const std::string &listfile)
{
    std::ifstream input(listfile);
    if (!input)
    {
        std::cerr << "ERROR: Could not open kmer file list " << listfile << " for reading.\n";
        std::_Exit(EXIT_FAILURE);
    }
    
    std::unordered_map<std::string, std::string> kmer_by_prefix;
    std::vector<std::string> kmerfiles;
    
    std::string line;
    while (std::getline(input, line))
    {
        std::vector<std::string> fields = split_text_fields(line);
        if (fields.empty() || fields[0][0] == '#') continue;
        
        std::string filepath = fields[0];
        std::string prefix = assigner_prefix_from_kmer_path(filepath);
        auto inserted = kmer_by_prefix.emplace(prefix, filepath);
        if (!inserted.second && inserted.first->second != filepath)
        {
            std::cerr << "ERROR: Duplicate matrix prefix in -K list: " << prefix
                      << " maps to both " << inserted.first->second << " and " << filepath << "\n";
            std::_Exit(EXIT_FAILURE);
        }
        
        if (inserted.second) kmerfiles.push_back(filepath);
    }
    
    return kmerfiles;
}

static std::unordered_set<std::string> parse_target_names(const std::string &text)
{
    std::unordered_set<std::string> targets;
    std::stringstream ss(text);
    std::string item;
    
    while (std::getline(ss, item, ','))
    {
        std::string trimmed = assigner_trim(item);
        if (!trimmed.empty()) targets.insert(trimmed);
    }
    
    return targets;
}

template <int dictsize>
void run(std::string &inputfile, std::string &kmerfile,std::string &outputfile)
{
    kmer_counter <32>counter;
    
    counter.getnorm(inputfile, kmerfile, outputfile);
    
    return;
}

template <int dictsize>
void run_targets(std::string &inputfile, std::string &kmerfile, std::string &outputfile, const std::unordered_set<std::string> &target_names)
{
    kmer_counter <32>counter;
    
    counter.gettargetpartition(inputfile, kmerfile, outputfile, target_names);
    
    return;
}

template <int dictsize>
void run_predefined(std::string &inputfile, std::string &kmerfile,std::string &outputfile)
{
    kmer_assigner <32> assigner;
    
    assigner.assigngroup(inputfile, kmerfile, outputfile);
    
    return;
}

template <int dictsize>
void run_predefined_mul(std::string &inputfile, vector<std::string> &kmerfiles, std::string &outputfile)
{
    kmer_assigner <32> assigner;
    
    assigner.assigngroup_mul(inputfile, kmerfiles, outputfile);
    
    return;
}


int main(int argc, const char * argv[]) {
    
    
    std::string inputfile = "";
    std::string outputfile = "";
    std::string kmerfile = "";
    std::string kmerfilelist = "";
    std::string target_names_text = "";
    std::unordered_set<std::string> target_names;
    std::vector<std::string> kmerfiles;
    
    const char* Argument = "";
        
    int mode = 0, kmer_size = 31;
    bool ifpredefined = 0;
    
    for (int i = 1; i < argc ; i++)
    {
        if (argv[i][0] == '-')
        {
            if (strcmp(argv[i], "-h")==0 or strcmp(argv[i], "--help")==0)
            {
                print_usage(argv[0]);
                return 0;
            }
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
        else if (strcmp(Argument, "-K")==0 or strcmp(Argument, "--kmer-list")==0)
        {
            kmerfilelist = argv[i];
            ifmultipleKmer = 1;
            kmerfiles = read_kmerfile_list(kmerfilelist);
        }
        else if (strcmp(Argument, "-t")==0 or strcmp(Argument, "--targets")==0)
        {
            target_names_text = argv[i];
            target_names = parse_target_names(target_names_text);
        }
        else if (strcmp(Argument, "-q")==0 or strcmp(Argument, "--query")==0)
        {
            queryfilepath = argv[i];
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
        
        else if (strcmp(Argument, "--predefined")==0)
        {
            ifpredefined = atoi(argv[i]);
        }
        
    }
    
    ifbed = has_file_suffix(inputfile, ".bed");
    if (ifbed)
    {
        if (queryfilepath.empty())
        {
            std::cerr << "ERROR: BED input requires -q/--query.\n";
            return 1;
        }
        load_query_contigs(queryfilepath);
    }
    
    if (ifmultipleKmer)
    {
        if (!ifpredefined || !ifbed)
        {
            std::cerr << "ERROR: -K/--kmer-list is only supported with BED input and --predefined 1.\n";
            return 1;
        }
        if (kmerfiles.empty())
        {
            std::cerr << "ERROR: -K/--kmer-list did not load any kmer files.\n";
            return 1;
        }
    }
    
    if (!target_names.empty() && ifpredefined)
    {
        std::cerr << "ERROR: -t/--targets is not compatible with --predefined.\n";
        return 1;
    }
    
    if (!target_names.empty())
    {
        run_targets<32>(inputfile, kmerfile, outputfile, target_names);
    }
    
    else if (ifpredefined && ifmultipleKmer)
    {
        run_predefined_mul<32>(inputfile, kmerfiles, outputfile);
    }
    
    else if (ifpredefined)
    {
        run_predefined<32>(inputfile, kmerfile, outputfile);
    }
    
    else
    {
        run<32>(inputfile, kmerfile, outputfile);
    }
    
    
    return 0;
}
