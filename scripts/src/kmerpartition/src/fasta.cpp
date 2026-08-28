//
//  fasta.cpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 10/17/21.
//  Copyright © 2021 USC_Mark. All rights reserved.
//

#include "fasta.hpp"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <sstream>

struct faidx_record
{
    unsigned long long length = 0;
    unsigned long long offset = 0;
    unsigned long long line_bases = 0;
    unsigned long long line_bytes = 0;
};

std::unordered_map<std::string, std::string> allquerycontigs;
bool ifqueryloaded = false;
static std::unordered_map<std::string, faidx_record> allqueryfaidx;
static std::unordered_map<std::string, std::string> allqueryseqs;

bool has_file_suffix(const std::string &path, const std::string &suffix)
{
    if (path.size() < suffix.size()) return false;
    return std::equal(suffix.rbegin(), suffix.rend(), path.rbegin(),
                      [](char a, char b)
                      {
                          return std::tolower(static_cast<unsigned char>(a)) ==
                                 std::tolower(static_cast<unsigned char>(b));
                      });
}

static bool is_fasta_path(const std::string &path)
{
    return has_file_suffix(path, ".fasta") || has_file_suffix(path, ".fa");
}

std::vector<std::string> split_text_fields(const std::string &line)
{
    std::stringstream ss(line);
    std::vector<std::string> fields;
    std::string field;
    
    while (ss >> field) fields.push_back(field);
    
    return fields;
}

bool parse_bed_record_fields(const std::string &line, std::vector<std::string> &fields)
{
    if (line.empty()) return false;
    if (line[0] == '#') return false;
    if (line.rfind("track", 0) == 0 || line.rfind("browser", 0) == 0) return false;
    
    fields = split_text_fields(line);
    if (fields.size() < 3) return false;
    
    try
    {
        std::stoull(fields[1]);
        std::stoull(fields[2]);
    }
    catch (...)
    {
        return false;
    }
    
    return true;
}

static void read_one_fai(const std::string &fastapath)
{
    std::string faipath = fastapath + ".fai";
    std::ifstream fai(faipath);
    
    if (!fai)
    {
        std::cerr << "ERROR: Could not open FASTA index " << faipath << " for reading.\n";
        std::_Exit(EXIT_FAILURE);
    }
    
    std::string line;
    while (std::getline(fai, line))
    {
        if (line.empty()) continue;
        
        std::vector<std::string> fields = split_text_fields(line);
        if (fields.size() < 5) continue;
        
        faidx_record record;
        record.length = std::stoull(fields[1]);
        record.offset = std::stoull(fields[2]);
        record.line_bases = std::stoull(fields[3]);
        record.line_bytes = std::stoull(fields[4]);
        
        allquerycontigs[fields[0]] = fastapath;
        allqueryfaidx[fields[0]] = record;
    }
}

static void load_one_fasta_to_ram(const std::string &fastapath)
{
    std::ifstream fasta_file(fastapath);
    if (!fasta_file)
    {
        std::cerr << "ERROR: Could not open query FASTA " << fastapath << " for reading.\n";
        std::_Exit(EXIT_FAILURE);
    }
    
    std::string line;
    std::string contig;
    
    while (std::getline(fasta_file, line))
    {
        if (line.empty()) continue;
        
        if (line[0] == '>')
        {
            std::vector<std::string> fields = split_text_fields(line.substr(1));
            contig = fields.empty() ? "" : fields[0];
            
            if (!contig.empty())
            {
                auto find_idx = allqueryfaidx.find(contig);
                if (find_idx != allqueryfaidx.end())
                {
                    allqueryseqs[contig].reserve(static_cast<size_t>(find_idx->second.length));
                }
                else
                {
                    allqueryseqs[contig];
                }
            }
            
            continue;
        }
        
        if (contig.empty()) continue;
        
        if (!line.empty() && line.back() == '\r') line.pop_back();
        allqueryseqs[contig] += line;
    }
    
    ifqueryloaded = true;
}

void load_query_contigs(const std::string &querypath)
{
    allquerycontigs.clear();
    allqueryfaidx.clear();
    allqueryseqs.clear();
    ifqueryloaded = false;
    
    if (querypath.empty())
    {
        std::cerr << "ERROR: BED input requires -q/--query.\n";
        std::_Exit(EXIT_FAILURE);
    }
    
    if (is_fasta_path(querypath))
    {
        read_one_fai(querypath);
        load_one_fasta_to_ram(querypath);
    }
    else
    {
        std::ifstream pathfile(querypath);
        if (!pathfile)
        {
            std::cerr << "ERROR: Could not open query path file " << querypath << " for reading.\n";
            std::_Exit(EXIT_FAILURE);
        }
        
        std::string line;
        while (std::getline(pathfile, line))
        {
            if (line.empty() || line[0] == '#') continue;
            
            std::vector<std::string> fields = split_text_fields(line);
            if (fields.empty()) continue;
            
            std::string fastapath = (fields.size() >= 2) ? fields[1] : fields[0];
            read_one_fai(fastapath);
        }
    }
    
    if (allquerycontigs.empty())
    {
        std::cerr << "ERROR: No contigs were loaded from query FASTA index input " << querypath << ".\n";
        std::_Exit(EXIT_FAILURE);
    }
}

std::string fetch_query_sequence(const std::string &contig, unsigned long long start, unsigned long long end)
{
    auto find_path = allquerycontigs.find(contig);
    auto find_idx = allqueryfaidx.find(contig);
    
    if (find_path == allquerycontigs.end() || find_idx == allqueryfaidx.end())
    {
        std::cerr << "ERROR: Contig " << contig << " was not found in any query FASTA index.\n";
        std::_Exit(EXIT_FAILURE);
    }
    
    const faidx_record &record = find_idx->second;
    start = std::min(start, record.length);
    end = std::min(end, record.length);
    if (end < start) std::swap(start, end);
    
    if (ifqueryloaded)
    {
        auto find_seq = allqueryseqs.find(contig);
        if (find_seq == allqueryseqs.end())
        {
            std::cerr << "ERROR: Contig " << contig << " was not loaded into query FASTA RAM cache.\n";
            std::_Exit(EXIT_FAILURE);
        }
        
        const std::string &seq = find_seq->second;
        start = std::min<unsigned long long>(start, seq.size());
        end = std::min<unsigned long long>(end, seq.size());
        if (end < start) std::swap(start, end);
        
        return seq.substr(static_cast<size_t>(start), static_cast<size_t>(end - start));
    }
    
    std::ifstream fasta_file(find_path->second, std::ios::binary);
    if (!fasta_file)
    {
        std::cerr << "ERROR: Could not open query FASTA " << find_path->second << " for reading.\n";
        std::_Exit(EXIT_FAILURE);
    }
    
    std::string seq;
    seq.reserve(static_cast<size_t>(end - start));
    
    unsigned long long pos = start;
    while (pos < end)
    {
        unsigned long long line_no = pos / record.line_bases;
        unsigned long long line_offset = pos % record.line_bases;
        unsigned long long take = std::min(end - pos, record.line_bases - line_offset);
        unsigned long long file_offset = record.offset + line_no * record.line_bytes + line_offset;
        
        std::string chunk(static_cast<size_t>(take), '\0');
        fasta_file.seekg(static_cast<std::streamoff>(file_offset), std::ios::beg);
        fasta_file.read(&chunk[0], static_cast<std::streamsize>(take));
        chunk.resize(static_cast<size_t>(fasta_file.gcount()));
        seq += chunk;
        pos += take;
    }
    
    return seq;
}

static std::string bed_score_from_assignment(const std::string &assignment)
{
    if (assignment.empty()) return "0";
    return assignment;
}

void write_bed_assignments(const std::string &inputfile, const std::string &outputfile, const std::vector<std::string> &scores)
{
    std::ifstream bed(inputfile);
    if (!bed)
    {
        std::cerr << "ERROR: Could not open BED input " << inputfile << " for reading.\n";
        std::_Exit(EXIT_FAILURE);
    }
    
    std::ofstream out(outputfile);
    if (!out)
    {
        std::cerr << "ERROR: Could not open BED output " << outputfile << " for writing.\n";
        std::_Exit(EXIT_FAILURE);
    }
    
    std::string line;
    size_t score_index = 0;
    
    while (std::getline(bed, line))
    {
        std::vector<std::string> fields;
        if (!parse_bed_record_fields(line, fields))
        {
            out << line << "\n";
            continue;
        }
        
        std::string score = (score_index < scores.size()) ? bed_score_from_assignment(scores[score_index]) : "0";
        score_index++;
        
        while (fields.size() < 4) fields.push_back(".");
        if (fields.size() < 5) fields.push_back(score);
        else fields[4] = score;
        
        for (size_t i = 0; i < fields.size(); ++i)
        {
            if (i) out << "\t";
            out << fields[i];
        }
        out << "\n";
    }
}

void write_bed_selected_assignments(const std::string &inputfile, const std::string &outputfile, const std::vector<bool> &selected, const std::string &score)
{
    std::ifstream bed(inputfile);
    if (!bed)
    {
        std::cerr << "ERROR: Could not open BED input " << inputfile << " for reading.\n";
        std::_Exit(EXIT_FAILURE);
    }
    
    std::ofstream out(outputfile);
    if (!out)
    {
        std::cerr << "ERROR: Could not open BED output " << outputfile << " for writing.\n";
        std::_Exit(EXIT_FAILURE);
    }
    
    std::string line;
    size_t bed_index = 0;
    
    while (std::getline(bed, line))
    {
        std::vector<std::string> fields;
        if (!parse_bed_record_fields(line, fields))
        {
            continue;
        }
        
        bool keep = bed_index < selected.size() && selected[bed_index];
        bed_index++;
        if (!keep) continue;
        
        while (fields.size() < 4) fields.push_back(".");
        if (fields.size() < 5) fields.push_back(score);
        else fields[4] = score;
        
        for (size_t i = 0; i < fields.size(); ++i)
        {
            if (i) out << "\t";
            out << fields[i];
        }
        out << "\n";
    }
}

void fasta::Load()
{
    if (strlen(filepath)<2) return;
    
    bed_mode = ifbed && has_file_suffix(filepath, ".bed");
    if (bed_mode && allquerycontigs.empty()) load_query_contigs(queryfilepath);
    
    fafile.open(filepath, ios_base::in);
    
    if (!fafile) {
        
        std::cerr << "ERROR: Could not open " << filepath << " for reading.\n" << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
}

bool fasta::nextBedLine(std::string &StrLine)
{
    if (pending_index < pending_lines.size())
    {
        StrLine = pending_lines[pending_index++];
        return true;
    }
    
    pending_lines.clear();
    pending_index = 0;
    
    std::string line;
    while (std::getline(fafile, line))
    {
        std::vector<std::string> fields;
        if (!parse_bed_record_fields(line, fields)) continue;
        
        unsigned long long start = std::stoull(fields[1]);
        unsigned long long end = std::stoull(fields[2]);
        std::string seq = fetch_query_sequence(fields[0], start, end);
        
        pending_lines.push_back(">" + line);
        for (size_t i = 0; i < seq.size(); i += 80)
        {
            pending_lines.push_back(seq.substr(i, 80));
        }
        
        StrLine = pending_lines[pending_index++];
        return true;
    }
    
    return false;
}

bool fasta::nextLine(std::string &StrLine)
{
    if (bed_mode) return nextBedLine(StrLine);
    return (bool)getline(fafile,StrLine);
}

void fasta::Close()
{
    fafile.close();
    pending_lines.clear();
    pending_index = 0;
}

void fasta::Reset()
{
    fafile.clear();
    fafile.seekg(0, std::ios::beg);
    pending_lines.clear();
    pending_index = 0;
}
