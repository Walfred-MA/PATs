//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License.
//  If you use this code, please cite our work.
//

#include "CramReader.hpp"
#include <string>
#include <chrono>
#include <thread>
#include <unistd.h>
#include <fcntl.h>

using namespace std;

int suppress_stderr() {
    int old_stderr = dup(STDERR_FILENO);  // Save current stderr
    int dev_null = open("/dev/null", O_WRONLY);
    dup2(dev_null, STDERR_FILENO);        // Redirect stderr to /dev/null
    close(dev_null);
    return old_stderr;                    // Return saved fd to restore later
}

void restore_stderr(int old_stderr) {
    dup2(old_stderr, STDERR_FILENO);      // Restore original stderr
    close(old_stderr);
}

void LoadHLARegion(hts_idx_t* index, bam_hdr_t* header, htsFile* samfile, std::vector<char *>* workregions)
{
    std::vector<char*> hla_regions;

    for (int i = 0; i < header->n_targets; ++i)
    {
        std::string contig(header->target_name[i]);

        if (contig.find("HLA") != std::string::npos ||
            contig.find("hla") != std::string::npos ||
            contig.find("MHC") != std::string::npos ||
            contig.find("chr6") != std::string::npos ||
            contig.find("HSCHR6") != std::string::npos ||
            contig.find("GL000") != std::string::npos)
        {
            workregions->push_back(strdup(header->target_name[i]));
        }
    }
}

bool hasHLARegion(const std::vector<char*>* workregions)
{
    for (const char* region : *workregions)
    {
        std::string reg(region);
        if (reg.find("HLA") != std::string::npos)
        {
            return true;
        }
    }
    return false;
}

void CramReader::Load()
{
    if (strlen(filepath)<2) return;
    
    std::string indexpath;
    if (strlen(filepath) >= 5 && std::string(filepath).substr(strlen(filepath) - 5) == ".cram")
    {
        indexpath = std::string(filepath) + ".crai";
    }
    else if (strlen(filepath) >= 4 && std::string(filepath).substr(strlen(filepath) - 4) == ".bam")
    {
        indexpath = std::string(filepath) + ".bai";
    }
    else if (strlen(filepath) >= 4 && std::string(filepath).substr(strlen(filepath) - 4) == ".sam")
    {
        indexpath = std::string(filepath) + ".sai";  // placeholder
    }
    else
    {
        std::cerr << "Unsupported file type: " << filepath << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
            
    int attemp = 200;
    do
    {
    try
    {
        samfile = hts_open(filepath, "r");
        
        if (!refpath.empty())
        {
            if (hts_set_fai_filename(samfile, refpath.c_str()) != 0)
            {
                throw std::runtime_error("Failed to set reference FASTA: " + refpath);
            }
        }
        
        indexdata = sam_index_load2(samfile, filepath, indexpath.c_str());
    }
        catch (const std::bad_alloc&)
    {
         samfile = NULL;
         indexdata = NULL;
    }
        
        
        if (!samfile || !indexdata)
        {
            std::cerr << "ERROR: Could not open " << filepath << " for reading:  attemp: "<< attemp <<" \n" << std::endl;
            std::this_thread::sleep_for(std::chrono::seconds(10));
        }
        
    } while ( (!samfile || !indexdata) && attemp-- > 0) ;

    
    indexdata = sam_index_load2(samfile, filepath, indexpath.c_str());
    
    header = sam_hdr_read(samfile);
}

void CramReader::LoadRegion(std::vector<char *>& bedregions)
{
    workregions = &bedregions;
    
    if (workregions->size() )
    {
        ifindex = 1;
        
        if (hasHLARegion(workregions)) LoadHLARegion(indexdata, header, samfile, workregions);
        
        int saved = suppress_stderr();
        iter = sam_itr_regarray(indexdata, header , workregions->data(), (int)workregions->size());
        restore_stderr(saved);
    }
        
    //strcpy(  regions[0], "chr1:1209512-1409512");
    
    //hts_reglist_t * reglist = hts_reglist_create_cram(regions, 1, &rcount, indexdata);
    
        
}

bool CramReader::nextLine_prt(uint8_t * &StrLine, size_t &rlen)
{
    //lock_guard<mutex> IO(IO_lock);

    while (true)
    {
        if (ifindex)
        {
            if (startumap)
            {
                iter = sam_itr_queryi(indexdata, HTS_IDX_NOCOOR, 0, 0);
                useunmap = 0;
                if (!iter || sam_itr_next(samfile, iter, SRread) < 0)
                {
                    return false;
                }
            }
            else
            {
                auto get = sam_itr_multi_next(samfile, iter, SRread);
                
                if ( get < 0)
                {
                    
                    if (useunmap)
                    {
                        hts_itr_destroy(iter);
                        iter = sam_itr_queryi(indexdata, HTS_IDX_NOCOOR, 0, 0);
                        useunmap = 0;
                        if (!iter || sam_itr_next(samfile, iter, SRread) < 0)
                        {
                            return false;
                        }
                        
                    }
                    else
                    {
                        return false;
                    }
                }
            }
        }
        else
        {
            if (sam_read1(samfile,header ,SRread)<0)
            {
                return false;
            }
        }
        
        if ((SRread->core.flag & filtertag) == 0)
        {
            break;
        }
    }

    rlen = SRread->core.l_qseq;


    //if (StrLine.size() <= readLength ) StrLine.resize(readLength+1);
    //uint8_t *q = bam_get_seq(SRread);

    //for (int i=0; i < readLength; i++) {StrLine[i]=seq_nt16_str[bam_seqi(q,i)];    }

    //StrLine[readLength] = '\0';
    StrLine = bam_get_seq(SRread);

    return true;
}
 
bool CramReader::nextLine(uint8_t * &StrLine, size_t &rlen, bam1_t* &SRread)
{
    std::lock_guard<std::mutex> IO(IO_lock);
    
    while (true)
    {
        if (ifindex)
        {
            if (startumap)
            {
                iter = sam_itr_queryi(indexdata, HTS_IDX_NOCOOR, 0, 0);
                useunmap = 0;
                if (!iter || sam_itr_next(samfile, iter, SRread) < 0)
                {
                    return false;
                }
            }
            else
            {
                auto get = sam_itr_multi_next(samfile, iter, SRread);
                
                if ( get < 0)
                {
                    
                    if (useunmap)
                    {
                        hts_itr_destroy(iter);
                        iter = sam_itr_queryi(indexdata, HTS_IDX_NOCOOR, 0, 0);
                        useunmap = 0;
                        if (!iter || sam_itr_next(samfile, iter, SRread) < 0)
                        {
                            return false;
                        }
                        
                    }
                    else
                    {
                        return false;
                    }
                }
            }
        }
        else
        {
            if (sam_read1(samfile,header ,SRread)<0)
            {
                return false;
            }
        }
        
        if ((SRread->core.flag & filtertag) == 0)
        {
            break;
        }
    }
        
    rlen = SRread->core.l_qseq;
    StrLine = bam_get_seq(SRread);

    return true;
}


bool CramReader::nextLine(std::string &StrLine)
{
    //lock_guard<mutex> IO(IO_lock);

    while (true)
    {
        if (ifindex)
        {
            if (startumap)
            {
                iter = sam_itr_queryi(indexdata, HTS_IDX_NOCOOR, 0, 0);
                useunmap = 0;
                if (!iter || sam_itr_next(samfile, iter, SRread) < 0)
                {
                    return false;
                }
            }
            else
            {
                auto get = sam_itr_multi_next(samfile, iter, SRread);
                
                if ( get < 0)
                {
                    
                    if (useunmap)
                    {
                        hts_itr_destroy(iter);
                        iter = sam_itr_queryi(indexdata, HTS_IDX_NOCOOR, 0, 0);
                        useunmap = 0;
                        if (!iter || sam_itr_next(samfile, iter, SRread) < 0)
                        {
                            return false;
                        }
                        
                    }
                    else
                    {
                        return false;
                    }
                }
            }
        }
        else
        {
            if (sam_read1(samfile,header ,SRread)<0)
            {
                return false;
            }
        }
        
        if ((SRread->core.flag & filtertag) == 0)
        {
            break;
        }
    }
    
    size_t readLength = SRread->core.l_qseq;

    if (StrLine.size() <= readLength ) StrLine.resize(readLength+2);
    uint8_t *q = bam_get_seq(SRread);

    for (int i=0; i < readLength; i++) {StrLine[i+1]=seq_nt16_str[bam_seqi(q,i)];    }

    StrLine[0] = '@';
    StrLine[readLength] = '\0';

    return true;
}

