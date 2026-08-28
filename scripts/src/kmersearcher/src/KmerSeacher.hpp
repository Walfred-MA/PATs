//
//  kmers.hpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 10/14/21.
//  Copyright © 2021 USC_Mark. All rights reserved.
//

#ifndef KmerSeacher
#define KmerSeacher


#include <stdio.h>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <exception>
#include <memory>
#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <limits>

#include "htslib/faidx.h"

#include "fasta.hpp"
#include "fastq.hpp"
#include "KmerStruct.hpp"
#include "KmerHash.hpp"

using namespace std;

using kmer64_dict_s = std::unordered_map<u128, uint16*, hash_128> ;
using kmer32_dict_s = std::unordered_map<ull, uint16* > ;

extern bool ifcache;
extern string cachefile;


template <typename T>
static void kmer_read_s(char base, int klen, int &current_size, T &current_kmer, T &reverse_kmer)
{
    int converted = 0;
    T reverse_converted;
    
    if (base == '\n' || base == ' ' || base == '\t') return;
    
    if (base_to_int(base, converted))
    {
        
        current_kmer <<= ( 8*sizeof(current_kmer) - 2*klen + 2 );
        current_kmer >>=  ( 8*sizeof(current_kmer) - 2*klen );
        current_kmer += converted;
        
        reverse_kmer >>= 2;
        reverse_converted = 0b11-converted;
        reverse_converted <<= (2*klen-2);
        reverse_kmer += reverse_converted;
        
    }
    
    else
    {
        current_size = 0;
        current_kmer = 0;
        reverse_kmer = 0;
    }
}

template <int dictsize>
class kmer_map
{
    using kmer_int = typename std::conditional<(dictsize>32), u128, ull>::type;
    //using target_maps_type = typename std::conditional<(dictsize>32), std::unordered_set<u128, hash_128>, std::unordered_set<ull>>::type;
    using kmer_dict_type = typename std::conditional<(dictsize>32), kmer64_dict_s, kmer32_dict_s>::type;
    using hotspot = std::tuple<int,int,int>;
    
    size_t totalkmers = 1;
    kmer_dict_type target_maps;
    Kmer_hash kmer_hash;
    
    int klen = 31 , knum = 0;
    bool iftarget = 0;
    vector<unsigned long long> targetlengths;
    ull totallength;
    uint targetindex = 0;
    const int hotspot_cutoff;
    const int kmer_window_size;
    int sample_nthreads = 1;
    
    std::atomic_uint restfileindex ;
    std::mutex Threads_lock;
    
    std::vector<std::thread*> threads;
    std::vector<std::string> inputfiles;
    std::vector<std::string> outputfiles;
    std::vector<std::string> prefixes;
    std::vector<std::string> allheaders;
public:
    
    kmer_map (int kmersize, int cutoff = 50, int kwindowsize = 1000):
        klen(31), hotspot_cutoff(cutoff), kmer_window_size(kwindowsize)
    {
        allheaders.reserve(1000000);
        allheaders.push_back(string(""));
    };
    
    ull read_targets(vector<string> &inputfiles, int nthreads);
    ull hash_targets(vector<string> &fastafiles, int nthreads);
    vector<kmer_int> read_unique_target_kmers(const string &targetfile);
    template <class consume_type>
    void process_target_kmers(
        vector<string> &targetfiles,
        int nthreads,
        consume_type consume);
    ull read_searchtarget(fasta &fastafile);
    void read_target(vector<string> &inputfiles, int nthreads);
    
    void read_files(std::vector<std::string>& inputs, std::vector<std::string>& outputs, std::vector<std::string>& prefs, int nthreads, int threads_per_sample = 1);
    
    void read_file();
            
    template <class typefile>
    void locate_kmer(typefile &fastafile, vector<hotspot*> &allhotspots, vector<ull>& sizes);

    void locate_kmer_indexed(const std::string &inputfile,
                             vector<hotspot*> &allhotspots,
                             vector<ull>& sizes);
    
    void write(const char* outfile, hotspot* allhotspots, ull sizes, uint index);
};



template <typename T1, typename T2>
inline static void initiate_counter_s(T2 &target_map, T1 &larger_kmer, uint16 index)
{
    
    typename T2::iterator map_find = target_map.find(larger_kmer);

    if (map_find == target_map.end())
    {
        target_map[larger_kmer] = (uint16*) malloc(sizeof(uint16)*2);
        target_map[larger_kmer][0] = 1;
        target_map[larger_kmer][1] = index ;
    }
    
    else
    {
        uint16* &data = map_find->second;

        if (std::find(data + 1, data + 1 + data[0], index) == map_find->second + 1 + data[0] )

        {
            if (data[0]%5 == 1)
            {
                data = (uint16*) realloc(data, sizeof(uint16)*(data[0]+1 + 5));
            }
            
            data[0]++;
            data[data[0]] = index ;
        }

    }
}


template <int dictsize>
vector<typename kmer_map<dictsize>::kmer_int>
kmer_map<dictsize>::read_unique_target_kmers(const string &targetfile)
{
    int current_size = 0;
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    std::string StrLine;
    StrLine.resize(MAX_LINE);
    vector<kmer_int> kmers;

    fasta fastafile(targetfile.c_str());
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '@': case '+': case '>':
                current_size = 0;
                continue;
            case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }

        for (char base: StrLine)
        {
            if (base == '\0') break;
            if (base == '\n' || base == ' ') continue;

            kmer_read_s(base, klen, ++current_size,
                        current_kmer, reverse_kmer);
            if (current_size < klen || base >= 'a') continue;

            kmers.push_back((current_kmer >= reverse_kmer)
                                ? current_kmer : reverse_kmer);
        }
    }

    // Preserve every non-FASTA target exactly as before.  Only raw .fa and
    // .fasta sequence targets require within-target canonical-k-mer
    // deduplication; generated .list targets are already unique.
    const bool is_fasta_target =
        (targetfile.size() >= 3 &&
         targetfile.compare(targetfile.size() - 3, 3, ".fa") == 0) ||
        (targetfile.size() >= 6 &&
         targetfile.compare(targetfile.size() - 6, 6, ".fasta") == 0);
    if (is_fasta_target)
    {
        sort(kmers.begin(), kmers.end());
        kmers.erase(unique(kmers.begin(), kmers.end()), kmers.end());
    }
    return kmers;
}


template <int dictsize>
template <class consume_type>
void kmer_map<dictsize>::process_target_kmers(
    vector<string> &targetfiles,
    int nthreads,
    consume_type consume)
{
    const size_t target_count = targetfiles.size();
    if (target_count == 0) return;

    const size_t worker_count = std::min(
        target_count,
        static_cast<size_t>(std::max(1, nthreads)));

    if (worker_count == 1)
    {
        for (size_t index = 0; index < target_count; ++index)
        {
            auto kmers = read_unique_target_kmers(targetfiles[index]);
            consume(index + 1, kmers);
        }
        return;
    }

    // A fixed-size ring bounds memory to at most one completed target vector
    // per worker.  Workers may parse/sort targets concurrently, while the
    // coordinator consumes them strictly in target-list order.
    vector<vector<kmer_int>> result_slots(worker_count);
    vector<size_t> result_indices(worker_count, target_count);
    vector<bool> result_ready(worker_count, false);
    std::mutex result_mutex;
    std::condition_variable result_changed;
    size_t next_to_assign = 0;
    size_t next_to_consume = 0;
    bool stop = false;
    std::exception_ptr worker_error;

    auto worker = [this, &targetfiles, target_count, worker_count,
                   &result_slots, &result_indices, &result_ready,
                   &result_mutex, &result_changed, &next_to_assign,
                   &next_to_consume, &stop, &worker_error]()
    {
        while (true)
        {
            size_t index = 0;
            {
                std::unique_lock<std::mutex> lock(result_mutex);
                result_changed.wait(lock, [&]() {
                    return stop || worker_error ||
                           next_to_assign >= target_count ||
                           next_to_assign < next_to_consume + worker_count;
                });
                if (stop || worker_error || next_to_assign >= target_count)
                    return;
                index = next_to_assign++;
            }

            try
            {
                auto kmers = read_unique_target_kmers(targetfiles[index]);
                {
                    std::lock_guard<std::mutex> lock(result_mutex);
                    const size_t slot = index % worker_count;
                    result_slots[slot] = std::move(kmers);
                    result_indices[slot] = index;
                    result_ready[slot] = true;
                }
                result_changed.notify_all();
            }
            catch (...)
            {
                {
                    std::lock_guard<std::mutex> lock(result_mutex);
                    if (!worker_error) worker_error = std::current_exception();
                }
                result_changed.notify_all();
                return;
            }
        }
    };

    vector<std::thread> workers;
    workers.reserve(worker_count);
    for (size_t worker_index = 0; worker_index < worker_count; ++worker_index)
        workers.emplace_back(worker);

    try
    {
        for (size_t index = 0; index < target_count; ++index)
        {
            vector<kmer_int> kmers;
            {
                std::unique_lock<std::mutex> lock(result_mutex);
                const size_t slot = index % worker_count;
                result_changed.wait(lock, [&]() {
                    return worker_error ||
                           (result_ready[slot] && result_indices[slot] == index);
                });
                if (worker_error) std::rethrow_exception(worker_error);

                kmers = std::move(result_slots[slot]);
                result_slots[slot].clear();
                result_ready[slot] = false;
                next_to_consume = index + 1;
            }
            result_changed.notify_all();
            consume(index + 1, kmers);
        }
    }
    catch (...)
    {
        {
            std::lock_guard<std::mutex> lock(result_mutex);
            stop = true;
        }
        result_changed.notify_all();
        for (auto &thread: workers) thread.join();
        throw;
    }

    {
        std::lock_guard<std::mutex> lock(result_mutex);
        stop = true;
    }
    result_changed.notify_all();
    for (auto &thread: workers) thread.join();
}


template <int dictsize>
ull kmer_map<dictsize>::read_targets(
    vector<std::string> &targetfiles,
    int nthreads)
{
    process_target_kmers(
        targetfiles, nthreads,
        [this](size_t, const vector<kmer_int> &kmers) {
        for (const auto &kmer: kmers)
            kmer_hash.previewvalue(kmer);
    });
    
    kmer_hash.initiatevalue();
    
    process_target_kmers(
        targetfiles, nthreads,
        [this](size_t index, const vector<kmer_int> &kmers) {
        const uint16 target_index = static_cast<uint16>(index);
        for (const auto &kmer: kmers)
            kmer_hash.addvalue(kmer, target_index);
    });
    
    kmer_hash.finalizevalue();
    
    return totalkmers;

};




template <int dictsize>
ull kmer_map<dictsize>::hash_targets(
    vector<std::string> &targetfiles,
    int nthreads)
{
    process_target_kmers(
        targetfiles, nthreads,
        [this](size_t, const vector<kmer_int> &kmers) {
        for (const auto &kmer: kmers)
            kmer_hash.preview(kmer);
    });
    
    totalkmers = kmer_hash.initiate();
    
    process_target_kmers(
        targetfiles, nthreads,
        [this](size_t, const vector<kmer_int> &kmers) {
        for (const auto &kmer: kmers)
            kmer_hash.add(kmer);
    });
    
    kmer_hash.finalize();
    
    return totalkmers;

};

template <int dictsize>
void kmer_map<dictsize>::read_target(
    vector<string> &inputfiles,
    int nthreads)
{
    const string cache_marker = cachefile + ".fasta_target_unique_v2";
    const bool has_fasta_target = std::any_of(
        inputfiles.begin(), inputfiles.end(), [](const string &targetfile) {
            return
                (targetfile.size() >= 3 &&
                 targetfile.compare(targetfile.size() - 3, 3, ".fa") == 0) ||
                (targetfile.size() >= 6 &&
                 targetfile.compare(
                     targetfile.size() - 6, 6, ".fasta") == 0);
        });
    const bool cache_exists =
        !cachefile.empty() && std::filesystem::exists(cachefile);
    const bool cache_is_current =
        cache_exists &&
        (!has_fasta_target || std::filesystem::exists(cache_marker));

    if (cache_is_current)
    {
        kmer_hash.loadhash(cachefile);
        
        targetindex = kmer_hash.totaltargets;
    }
    else
    {
        if (cache_exists && has_fasta_target)
        {
            cout << "Rebuilding legacy FASTA-target cache without "
                    "per-target k-mer deduplication: " << cachefile << endl;
        }

        hash_targets(inputfiles, nthreads);
        
        read_targets(inputfiles, nthreads);
        
        targetindex = inputfiles.size() + 1;
        kmer_hash.totaltargets = targetindex;
        
        if (ifcache == 1)
        {
            kmer_hash.savehash(cachefile);
            if (has_fasta_target && std::filesystem::exists(cachefile))
            {
                std::ofstream marker(cache_marker);
                marker << "fasta_target_unique_v2\n";
            }
        }
    }
}

template <int dictsize>
template <class typefile>
void kmer_map<dictsize>::locate_kmer(typefile &fastafile, vector<hotspot*> &allhotspots, vector<ull>& allhotspot_sizes)
{
    
    cout << "start reading: "<< fastafile.filepath << endl;
    
    int current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer=0;
    
    unordered_map<ull,uint16> counter;
    //const int hotspot_cutoff = 50;
    const int num_targets = (int) allhotspots.size();
    
    std::vector<std::vector<int>> kmer_posis;
    kmer_posis.resize(num_targets);
    for (auto &kmer_posi: kmer_posis)
    {
        kmer_posi.resize(hotspot_cutoff+1);
        std::fill(kmer_posi.begin(), kmer_posi.end(), -10000000);
    }
    
    std::vector<int> num_kmers(num_targets, 0);
    int posi = 0;
    int posi_start = 0;
    //uint64_t ifmasked = 0;
    //int num_masked = 0 ;
    std::string StrLine;

    
    vector<uint16> allindex(MAX_UINT16-1,0);
    int header = 0;
    
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '@':  case '+':
                current_size = 0;
                continue;
            case ' ': case '\n': case '\t':
                continue;
            case '>':
                current_size = 0;
                Threads_lock.lock();
                allheaders.push_back(StrLine.substr(1));
                header = allheaders.size() - 1;
                Threads_lock.unlock();
                std::fill(num_kmers.begin(), num_kmers.end(), 0);
                for (auto &kmer_posi: kmer_posis)
                {
                    std::fill(kmer_posi.begin(), kmer_posi.end(), -10000000);
                }
                posi = 0;
                continue;
            default:
                break;
        }
        
        for (auto base: StrLine)
        {
            if (base == '\0')
            {
                current_size = 0;
                header = 0;
                std::fill(num_kmers.begin(), num_kmers.end(), 0);
                for (auto &kmer_posi: kmer_posis)
                {
                    std::fill(kmer_posi.begin(), kmer_posi.end(), -10000000);
                }
                posi = 0;
                
                break;
            }
                        
            if (base == '\n' || base == ' ' || base == '\t') continue;
            
            posi++;
            
            kmer_read_s(base, klen, ++current_size, current_kmer, reverse_kmer);
            
            if (current_size < klen) continue;
            //num_masked=!!ifmasked;while(ifmasked&=ifmasked-1)num_masked++;
           
            //if (2*num_masked > klen) continue;
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
            
            auto find = kmer_hash.findvalue(larger_kmer, allindex);
            if (find == 0) continue;
            
            for (size_t j = 0; j < find; ++j)
            {
                uint16 i = allindex[j] ;
                
                hotspot*& hotspots = allhotspots[i];
                ull& hotspots_size = allhotspot_sizes[i];
                auto& kmer_posi = kmer_posis[i];
                auto& num_kmer = num_kmers[i];
                
                kmer_posi[num_kmer] = posi;
                
                num_kmer = (num_kmer + 1) % hotspot_cutoff;

                // After inserting the current hit, num_kmer points to the
                // oldest of the last hotspot_cutoff distinct query-position
                // hits.  The previous implementation read the slot before
                // insertion and therefore required cutoff + 1 hits.  With a
                // cutoff of 1000 and a <1000-bp window, that was impossible
                // once repeated target k-mers were correctly deduplicated.
                posi_start = kmer_posi[num_kmer];
               
 
                const int hotspot_window = kmer_window_size;
                if (abs(posi - posi_start) < hotspot_window)
                {
                    if (hotspots_size &&
                        header == std::get<0>(hotspots[hotspots_size-1]) &&
                        abs(posi_start - std::get<2>(hotspots[hotspots_size-1]))
                            < hotspot_window)
                    {
                        std::get<2>(hotspots[hotspots_size-1]) = posi ;
                    }
                    
                    else
                    {
    
                        if (hotspots_size%1000 == 1)
                        {
                            hotspots = (tuple<int,int,int>*) realloc(hotspots, (hotspots_size+1000)*sizeof(tuple<int,int,int>));
                        }
                        
                        hotspots[hotspots_size++] = std::make_tuple(header, posi_start, posi);
 
                    }
                }
            }
        }
    }
        
    fastafile.Close();
};


template <int dictsize>
void kmer_map<dictsize>::locate_kmer_indexed(
    const std::string &inputfile,
    vector<hotspot*> &allhotspots,
    vector<ull>& allhotspot_sizes)
{
    cout << "start reading indexed FASTA: " << inputfile << endl;

    const std::string fai_path = inputfile + ".fai";
    std::unique_ptr<faidx_t, decltype(&fai_destroy)> metadata(
        fai_load3(inputfile.c_str(), fai_path.c_str(), nullptr, 0),
        &fai_destroy);
    if (!metadata)
    {
        throw std::runtime_error("Cannot load FASTA index: " + fai_path);
    }

    const int num_contigs = faidx_nseq(metadata.get());
    if (num_contigs < 0)
    {
        throw std::runtime_error("Cannot read contigs from FASTA index: " +
                                 fai_path);
    }

    vector<string> contig_names(static_cast<size_t>(num_contigs));
    vector<hts_pos_t> contig_lengths(static_cast<size_t>(num_contigs));
    vector<int> header_ids(static_cast<size_t>(num_contigs));

    for (int i = 0; i < num_contigs; ++i)
    {
        const char *name = faidx_iseq(metadata.get(), i);
        if (name == nullptr)
        {
            throw std::runtime_error("Invalid contig entry in FASTA index: " +
                                     fai_path);
        }

        contig_names[static_cast<size_t>(i)] = name;
        const hts_pos_t length = faidx_seq_len64(metadata.get(), name);
        if (length < 0 ||
            length > static_cast<hts_pos_t>(std::numeric_limits<int>::max()))
        {
            throw std::runtime_error(
                "Unsupported contig length for " + std::string(name) +
                " in " + inputfile);
        }
        contig_lengths[static_cast<size_t>(i)] = length;
    }
    metadata.reset();

    // Allocate header IDs in FASTA/.fai order.  Results are merged in this
    // same order after workers finish, so -N does not change output ordering.
    {
        std::lock_guard<std::mutex> lock(Threads_lock);
        for (int i = 0; i < num_contigs; ++i)
        {
            allheaders.push_back(contig_names[static_cast<size_t>(i)]);
            header_ids[static_cast<size_t>(i)] =
                static_cast<int>(allheaders.size() - 1);
        }
    }

    struct target_state
    {
        explicit target_state(size_t cutoff):
            kmer_positions(cutoff, -10000000)
        {}

        vector<int> kmer_positions;
        size_t next_position = 0;
        vector<pair<int, int>> hotspots;
    };

    struct contig_result
    {
        unordered_map<uint16, vector<pair<int, int>>> hotspots_by_target;
    };

    vector<contig_result> results(static_cast<size_t>(num_contigs));
    if (num_contigs == 0) return;

    const size_t worker_count = std::min(
        static_cast<size_t>(num_contigs),
        static_cast<size_t>(std::max(1, sample_nthreads)));
    std::atomic_size_t next_contig(0);
    std::atomic_bool stop(false);
    std::mutex indexed_error_mutex;
    std::exception_ptr indexed_error;

    auto worker = [this, &inputfile, &fai_path, &contig_names,
                   &contig_lengths, &results, &next_contig, &stop,
                   &indexed_error_mutex, &indexed_error, num_contigs]()
    {
        try
        {
            std::unique_ptr<faidx_t, decltype(&fai_destroy)> fai(
                fai_load3(inputfile.c_str(), fai_path.c_str(), nullptr, 0),
                &fai_destroy);
            if (!fai)
            {
                throw std::runtime_error("Cannot load FASTA index: " +
                                         fai_path);
            }

            vector<uint16> allindex(static_cast<size_t>(MAX_UINT16) + 1, 0);
            constexpr hts_pos_t sequence_chunk_size = 4 * 1024 * 1024;

            while (!stop.load())
            {
                const size_t contig_index = next_contig.fetch_add(1);
                if (contig_index >= static_cast<size_t>(num_contigs)) break;

                const string &contig_name = contig_names[contig_index];
                const hts_pos_t contig_length = contig_lengths[contig_index];
                unordered_map<uint16, target_state> states;

                int current_size = 0;
                kmer_int current_kmer = 0;
                kmer_int reverse_kmer = 0;
                int position = 0;

                for (hts_pos_t start = 0; start < contig_length;
                     start += sequence_chunk_size)
                {
                    const hts_pos_t end = std::min(
                        contig_length - 1, start + sequence_chunk_size - 1);
                    hts_pos_t fetched_length = 0;
                    std::unique_ptr<char, decltype(&std::free)> sequence(
                        faidx_fetch_seq64(fai.get(), contig_name.c_str(),
                                          start, end, &fetched_length),
                        &std::free);
                    const hts_pos_t expected_length = end - start + 1;
                    if (!sequence || fetched_length != expected_length)
                    {
                        throw std::runtime_error(
                            "Cannot fetch contig " + contig_name + " from " +
                            inputfile);
                    }

                    for (hts_pos_t offset = 0; offset < fetched_length;
                         ++offset)
                    {
                        const char base = sequence.get()[offset];
                        ++position;
                        kmer_read_s(base, klen, ++current_size,
                                    current_kmer, reverse_kmer);
                        if (current_size < klen) continue;

                        const auto larger_kmer =
                            (current_kmer >= reverse_kmer)
                                ? current_kmer : reverse_kmer;
                        const uint found =
                            kmer_hash.findvalue(larger_kmer, allindex);
                        for (size_t j = 0; j < found; ++j)
                        {
                            const uint16 target = allindex[j];
                            auto inserted = states.try_emplace(
                                target,
                                static_cast<size_t>(hotspot_cutoff));
                            target_state &state = inserted.first->second;

                            state.kmer_positions[state.next_position] = position;
                            state.next_position =
                                (state.next_position + 1) %
                                static_cast<size_t>(hotspot_cutoff);
                            const int position_start =
                                state.kmer_positions[state.next_position];

                            if (std::abs(position - position_start) <
                                kmer_window_size)
                            {
                                if (!state.hotspots.empty() &&
                                    std::abs(
                                        position_start -
                                        state.hotspots.back().second) <
                                        kmer_window_size)
                                {
                                    state.hotspots.back().second = position;
                                }
                                else
                                {
                                    state.hotspots.emplace_back(position_start,
                                                                position);
                                }
                            }
                        }
                    }
                }

                auto &contig_hotspots =
                    results[contig_index].hotspots_by_target;
                contig_hotspots.reserve(states.size());
                for (auto &entry: states)
                {
                    if (!entry.second.hotspots.empty())
                    {
                        contig_hotspots.emplace(
                            entry.first, std::move(entry.second.hotspots));
                    }
                }
            }
        }
        catch (...)
        {
            {
                std::lock_guard<std::mutex> lock(indexed_error_mutex);
                if (!indexed_error) indexed_error = std::current_exception();
            }
            stop.store(true);
        }
    };

    vector<std::thread> workers;
    workers.reserve(worker_count);
    for (size_t i = 0; i < worker_count; ++i)
        workers.emplace_back(worker);
    for (auto &thread: workers) thread.join();

    if (indexed_error) std::rethrow_exception(indexed_error);

    // Merge per-contig results in index order.  Each target therefore retains
    // the same contig and coordinate order as the serial FASTA reader.
    for (size_t contig_index = 0;
         contig_index < static_cast<size_t>(num_contigs); ++contig_index)
    {
        for (auto &target_result:
             results[contig_index].hotspots_by_target)
        {
            const size_t target = static_cast<size_t>(target_result.first);
            if (target >= allhotspots.size())
            {
                throw std::runtime_error(
                    "Invalid target index while processing " + inputfile);
            }

            hotspot*& hotspots = allhotspots[target];
            ull& hotspot_count = allhotspot_sizes[target];
            for (const auto &interval: target_result.second)
            {
                if (hotspot_count % 1000 == 0)
                {
                    void *resized = std::realloc(
                        hotspots,
                        static_cast<size_t>(hotspot_count + 1000) *
                            sizeof(hotspot));
                    if (resized == nullptr) throw std::bad_alloc();
                    hotspots = static_cast<hotspot*>(resized);
                }
                hotspots[hotspot_count++] = std::make_tuple(
                    header_ids[contig_index], interval.first,
                    interval.second);
            }
        }
    }
}


template <int dictsize>
void kmer_map<dictsize>::read_files(std::vector<std::string>& inputs, std::vector<std::string>& outputs, std::vector<std::string>& prefs, int nthreads, int threads_per_sample)
{
    
    inputfiles = inputs;
    outputfiles = outputs;
    prefixes = prefs;
    sample_nthreads = threads_per_sample;
    restfileindex = 0;

    if (sample_nthreads > 1)
    {
        for (const string &inputfile: inputfiles)
        {
            const string fai_path = inputfile + ".fai";
            std::unique_ptr<faidx_t, decltype(&fai_destroy)> fai(
                fai_load3(inputfile.c_str(), fai_path.c_str(), nullptr, 0),
                &fai_destroy);
            if (!fai)
            {
                throw std::runtime_error("Cannot load FASTA index: " +
                                         fai_path);
            }
        }
    }

    std::mutex worker_error_mutex;
    std::exception_ptr worker_error;
    auto worker = [this, &worker_error_mutex, &worker_error]()
    {
        try
        {
            read_file();
        }
        catch (...)
        {
            std::lock_guard<std::mutex> lock(worker_error_mutex);
            if (!worker_error) worker_error = std::current_exception();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(static_cast<size_t>(nthreads));
    for(int i = 0; i < nthreads; ++i) threads.emplace_back(worker);
    for(auto &thread: threads) thread.join();

    if (worker_error) std::rethrow_exception(worker_error);
    
}

template <int dictsize>
void kmer_map<dictsize>::read_file()
{
    
    while (restfileindex < inputfiles.size())
    {
        Threads_lock.lock();
        
        int inputindex = restfileindex ++;
        
        Threads_lock.unlock();
        
        if (inputindex >= inputfiles.size()) break;
        
        vector<hotspot*> allhotspots;
        allhotspots.resize(targetindex);
        
        vector<ull> allhotspots_size;
        allhotspots_size.resize(targetindex);
        
        for(hotspot* &hotspots:allhotspots)
        {
            hotspots = (tuple<int,int,int>*) malloc(sizeof(tuple<int,int,int>));
        }
        
        const char* inputfile = inputfiles[inputindex].c_str();
        string outputfile_ = outputfiles[MIN(inputindex-1,outputfiles.size()-1)];
        if (outputfile_[outputfile_.length() -1] == '/') outputfile_ = "";
        
        
        string prefix = "";
        if (prefixes.size() > inputindex ) prefix = prefixes[inputindex];
 
        int pathlen = (int)strlen(inputfile);

        try
        {
            if (sample_nthreads > 1)
            {
                locate_kmer_indexed(inputfiles[inputindex], allhotspots,
                                    allhotspots_size);
            }
            else if (pathlen > 2 &&
                     strcmp(inputfile + (pathlen - 3), ".gz") == 0)
            {
                fastq readsfile(inputfile);
                locate_kmer(readsfile, allhotspots, allhotspots_size);
            }
            else
            {
                fasta readsfile(inputfile);
                locate_kmer(readsfile, allhotspots, allhotspots_size);
            }

            for (int j = 1; j < allhotspots.size(); ++j)
            {
                string outputfile = outputfile_;
                if (outputfile == "")
                {
                    outputfile = outputfiles[MIN(j-1,outputfiles.size()-1)] ;
                }
                outputfile = outputfile + prefix + "_hotspot.txt" ;

                write(outputfile.c_str(), allhotspots[j],
                      allhotspots_size[j], j);
            }
        }
        catch (...)
        {
            for (hotspot* &hotspots: allhotspots) free(hotspots);
            throw;
        }

        for (hotspot* &hotspots: allhotspots) free(hotspots);
        
    }
        
}




template <int dictsize>
void kmer_map<dictsize>::write(const char * outputfile, hotspot* hotspots, ull size, uint index)
{
    
    FILE *fwrite;
    
    if (index == 1)
    {
        fwrite=fopen(outputfile, "w");
    }
    else
    {
        fwrite=fopen(outputfile, "a");
    }
    
    if (fwrite==NULL)
    {
        std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
        
        std::_Exit(EXIT_FAILURE);
    }
            
    
    for (size_t i = 0; i < size; ++i)
    {
        hotspot hotspot = hotspots[i];
        
        std::string header;
        {
            std::lock_guard<std::mutex> lock(Threads_lock);
            header = allheaders[std::get<0>(hotspot)];
        }
        // get first whitespace-separated token
        std::string first;
        {
            std::istringstream iss(header);
            iss >> first;   // reads up to first whitespace
        }
        
        fprintf(fwrite,"%s\t%d\t%d\t%d\n", first.c_str(), index, std::get<1>(hotspot), std::get<2>(hotspot));
        
    }

    fclose(fwrite);
        
    return ;
}


#endif /* kmers_hpp */
