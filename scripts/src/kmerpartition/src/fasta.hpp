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
#include <unordered_map>

using namespace std;

extern bool ifbed;
extern std::string queryfilepath;
extern std::unordered_map<std::string, std::string> allquerycontigs;
extern bool ifqueryloaded;

void load_query_contigs(const std::string &querypath);
void write_bed_assignments(const std::string &inputfile, const std::string &outputfile, const std::vector<std::string> &scores);
void write_bed_selected_assignments(const std::string &inputfile, const std::string &outputfile, const std::vector<bool> &selected, const std::string &score);
bool has_file_suffix(const std::string &path, const std::string &suffix);
std::vector<std::string> split_text_fields(const std::string &line);
bool parse_bed_record_fields(const std::string &line, std::vector<std::string> &fields);
std::string fetch_query_sequence(const std::string &contig, unsigned long long start, unsigned long long end);

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
    bool bed_mode = false;
    std::vector<std::string> pending_lines;
    size_t pending_index = 0;
    
    bool nextBedLine(std::string &StrLine);
        
};


#endif /* fasta_hpp */
