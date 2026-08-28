#ifndef Blacklist_hpp
#define Blacklist_hpp

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

class BlacklistIntervals
{
    using interval_type = std::pair<unsigned long long, unsigned long long>;
    std::unordered_map<std::string, std::vector<interval_type>> intervals;
    bool finalized = false;

    static std::string trim_header_prefix(const std::string &line)
    {
        if (!line.empty() && line[0] == '>') return line.substr(1);
        return line;
    }

    static bool parse_coord_text(std::string coord, std::string &contig, unsigned long long &start, unsigned long long &end)
    {
        while (!coord.empty() && std::isspace(static_cast<unsigned char>(coord.back()))) coord.pop_back();
        if (!coord.empty() && (coord.back() == '+' || coord.back() == '-')) coord.pop_back();

        size_t colon = coord.rfind(':');
        size_t dash = coord.rfind('-');
        if (colon == std::string::npos || dash == std::string::npos || dash <= colon + 1) return false;

        contig = coord.substr(0, colon);
        try
        {
            start = std::stoull(coord.substr(colon + 1, dash - colon - 1));
            end = std::stoull(coord.substr(dash + 1));
        }
        catch (...)
        {
            return false;
        }

        if (end < start) std::swap(start, end);
        return !contig.empty() && end > start;
    }

    static bool parse_fasta_header_coord(const std::string &line, std::string &contig, unsigned long long &start, unsigned long long &end)
    {
        std::istringstream iss(trim_header_prefix(line));
        std::string first;
        std::string coord;
        if (!(iss >> first)) return false;
        if (!(iss >> coord)) return false;
        return parse_coord_text(coord, contig, start, end);
    }

    static bool parse_bed_line(const std::string &line, std::string &contig, unsigned long long &start, unsigned long long &end)
    {
        if (line.empty() || line[0] == '#') return false;

        std::istringstream iss(line);
        std::string start_text;
        std::string end_text;
        if (!(iss >> contig >> start_text >> end_text)) return false;
        if (contig == "track" || contig == "browser") return false;

        try
        {
            start = std::stoull(start_text);
            end = std::stoull(end_text);
        }
        catch (...)
        {
            return false;
        }

        if (end < start) std::swap(start, end);
        return !contig.empty() && end > start;
    }

    void add_interval(const std::string &contig, unsigned long long start, unsigned long long end)
    {
        if (contig.empty() || end <= start) return;
        intervals[contig].emplace_back(start, end);
        finalized = false;
    }

public:
    struct HeaderRegion
    {
        std::string contig;
        unsigned long long start = 0;
        unsigned long long end = 0;
        bool valid = false;
    };

    bool enabled() const
    {
        return !intervals.empty();
    }

    void load_file(const std::string &path)
    {
        std::ifstream fh(path);
        if (!fh)
        {
            std::cerr << "ERROR: Could not open blacklist file: " << path << std::endl;
            std::_Exit(EXIT_FAILURE);
        }

        std::string line;
        bool saw_fasta_header = false;
        size_t loaded = 0;

        while (std::getline(fh, line))
        {
            if (line.empty()) continue;

            std::string contig;
            unsigned long long start = 0;
            unsigned long long end = 0;

            if (line[0] == '>')
            {
                saw_fasta_header = true;
                if (parse_fasta_header_coord(line, contig, start, end))
                {
                    add_interval(contig, start, end);
                    loaded++;
                }
                continue;
            }

            if (saw_fasta_header) continue;

            if (parse_bed_line(line, contig, start, end))
            {
                add_interval(contig, start, end);
                loaded++;
            }
        }

        std::cout << "loaded blacklist intervals from " << path << ": " << loaded << std::endl;
    }

    void load_files(const std::vector<std::string> &paths)
    {
        for (const std::string &path: paths) load_file(path);
        finalize();
    }

    void finalize()
    {
        if (finalized) return;

        for (auto &[contig, contig_intervals]: intervals)
        {
            if (contig_intervals.empty()) continue;

            std::sort(contig_intervals.begin(), contig_intervals.end());
            std::vector<interval_type> merged;
            merged.reserve(contig_intervals.size());

            for (const auto &span: contig_intervals)
            {
                if (merged.empty() || span.first > merged.back().second)
                {
                    merged.push_back(span);
                }
                else
                {
                    merged.back().second = std::max(merged.back().second, span.second);
                }
            }

            contig_intervals.swap(merged);
        }

        finalized = true;
    }

    bool overlaps(const std::string &contig, unsigned long long start, unsigned long long end)
    {
        if (end < start) std::swap(start, end);
        if (end <= start) return false;

        finalize();
        auto found = intervals.find(contig);
        if (found == intervals.end()) return false;

        const auto &contig_intervals = found->second;
        auto it = std::lower_bound(
            contig_intervals.begin(),
            contig_intervals.end(),
            interval_type(start, start),
            [](const interval_type &lhs, const interval_type &rhs)
            {
                return lhs.second <= rhs.first;
            }
        );

        return it != contig_intervals.end() && it->first < end && start < it->second;
    }

    bool overlaps_header(const std::string &line)
    {
        std::string contig;
        unsigned long long start = 0;
        unsigned long long end = 0;
        if (!parse_fasta_header_coord(line, contig, start, end)) return false;
        return overlaps(contig, start, end);
    }

    HeaderRegion parse_header_region(const std::string &line) const
    {
        HeaderRegion region;
        region.valid = parse_fasta_header_coord(line, region.contig, region.start, region.end);
        return region;
    }

    bool overlaps_sequence_span(const HeaderRegion &region, unsigned long long seq_start, unsigned long long seq_end)
    {
        if (!region.valid || seq_end <= seq_start) return false;
        return overlaps(region.contig, region.start + seq_start, region.start + seq_end);
    }
};

extern BlacklistIntervals blacklist_intervals;

#endif /* Blacklist_hpp */
