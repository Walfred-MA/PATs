//
//  KmerTree.cpp
//  KmerTree
//
//  Created by walfred on 6/18/25.
//

#include "KmerTree.hpp"

inline void strsplit(const std::string& str, std::vector<normvalue>& eles, int dim)
{
    const char deli = ',';
    std::string::size_type start = 0;
    size_t len = std::strlen(str.c_str());
    size_t end = str.find(deli, start);

    for (int i = 0 ; i < dim ; ++i)
    {
        eles[i] = std::atof(str.substr(start, end - start).c_str());
        start = end + 1;
        end = str.find(deli, start);
        if (end == std::string::npos) end = len;
    }
}

void find_identical_rows(const std::vector<std::vector<normvalue>>& matrix, vector<int> &redundant_indices)
{
    const double tol = 1e-6;
    int n = matrix.size();
    redundant_indices.resize(n, 0);
    int samecount = 0;
    for (int i = 0; i < n; ++i)
    {
        const auto& row_i = matrix[i];
        for (int j = 0; j < i; ++j)
        {  // only check earlier rows
            const auto& row_j = matrix[j];
            bool ifsame = 1;
            for (size_t k = 0; k < row_i.size(); ++k)
            {
                if (std::abs(row_i[k] - row_j[k]) > tol)
                {
                    ifsame = 0;
                    break;
                }
            }
            if (ifsame)
            {
                samecount ++;
                redundant_indices[i] = j+1;
                break;
            }
        }
    }
    
}

 void compute_distance_matrix(
    const std::vector<std::vector<normvalue>>& matrix, const std::vector<int>& redundant_flags, std::vector<std::vector<normvalue>> &dm)
{
     int N = matrix.size();
     std::vector<int> useindex;
     useindex.reserve(N);
     for (int i = 0; i < N; ++i)
     {
         if (redundant_flags[i] == 0) useindex.push_back(i);
     }
     
     int M = useindex.size();
     dm.resize(M, std::vector<normvalue>(M, 0));
     
     
    std::vector<double> vars(N);
    for (int i = 0; i < N; ++i)
        vars[i] = std::sqrt((double)matrix[i][i]);

    for (int i = 0; i < M; ++i)
    {
        int x = useindex[i];
        double var1 = std::max(1.0, vars[x]);

        for (int j = 0; j < i; ++j)
        {
            int y = useindex[j];
            double var2 = std::max((double)1.0, vars[y]);

            double norm = var1 * var2;
            double sim = matrix[x][y];
            double dist = (double)1.0 - sim / norm;
            dm[i][j] = dm[j][i] = dist;
            
        }
    }
     
     
}


void kmer_tree::loadnorm(string &normfilepath, vector<vector<normvalue>> &data, const int dim)
{
    
    gzfile normfile(normfilepath.c_str());

    const int buffer_size = 10000000;
    std::string line(buffer_size, '\0');
    vector<normvalue> values(dim,0);
    data.resize(dim, vector<normvalue> (dim,0));
    // Read remaining lines
    int row_index = 0;
    while (normfile.nextLine(line))
    {

        strsplit(line, values, dim - row_index);
        for (int i = row_index; i < dim; ++i)
        {
            data[row_index][i] = values[i-row_index];
            data[i][row_index] = data[row_index][i];
        }
            
        ++row_index;
    }

    normfile.Close();

}

void kmer_tree::loadcontigs(const std::string& inputfile, std::vector<std::string>& names)
{
    std::ifstream file(inputfile);
    std::string line;

    while (std::getline(file, line))
    {
        if (!line.empty() && line[0] == '>')
        {
            size_t end = line.find_first_of(" \t\n\r", 1);
            names.push_back(line.substr(1, end - 1));
        }
    }
}

std::vector<int> ExtractTreeOrder(const std::string& treetext, const std::vector<std::string>& samplelist)
{
    std::vector<std::string> names;
    std::stringstream ss(treetext.substr(0, treetext.size() - 1)); // remove final ';'
    std::string token;
    
    while (std::getline(ss, token, ','))
    {
        size_t colon_pos = token.find(':');
        std::string name = token.substr(0, colon_pos);
        name.erase(std::remove(name.begin(), name.end(), '('), name.end());
        names.push_back(name);
    }
    
    std::unordered_map<std::string, int> nametoindex;
    for (int i = 0; i < samplelist.size(); ++i)
    {
        nametoindex[samplelist[i]] = i;
    }
    
    std::vector<int> treeorder;
    for (const auto& name : names)
    {
        treeorder.push_back(nametoindex[name]);
    }
    
    return treeorder;
}


void ReorderFile(const std::string& infile, const std::string& outfile,
                 const std::vector<int>& treeorder, const std::string& name = "")
{
    std::ifstream fin(infile);
    if (!fin)
    {
        throw std::runtime_error("Failed to open input file: " + infile);
    }

    std::vector<std::string> raw_reads;
    std::string line, block, full_content;
    while (std::getline(fin, line))
    {
        full_content += line + "\n";
    }
    fin.close();

    std::stringstream ss(full_content);
    std::string segment = "";
    while (std::getline(ss, segment, '>'))
    {
        if (!segment.empty())
            raw_reads.push_back(segment);
    }
    raw_reads.push_back(segment);
    
    std::ofstream fout(outfile);
    if (!fout)
    {
        throw std::runtime_error("Failed to open output file: " + outfile);
    }

    for (int index : treeorder)
    {
        if (index < 0 || index >= static_cast<int>(raw_reads.size()))
        {
            throw std::out_of_range("Index out of bounds in treeorder.");
        }

        std::istringstream entry_stream(raw_reads[index]);
        std::string header, line, entry;
        std::getline(entry_stream, header);

        if (!name.empty())
        {
            header = name + header;
        }

        entry = ">" + header;

        while (std::getline(entry_stream, line))
        {
            entry += "\n" + line;
        }

        fout << entry << "\n";
    }

    fout.close();
}

void ReorderMatrix(const std::vector<std::vector<normvalue>>& matrix, const std::vector<int>& order, const std::string& filename)
{
    gzFile gzfile = gzopen(filename.c_str(), "wb");
    if (!gzfile)
    {
        throw std::runtime_error("Failed to open gzip file for writing: " + filename);
    }

    for (int i = 0; i < order.size(); ++i)
    {
        const auto& row = matrix[order[i]];
        for (int j = i; j < order.size(); ++j)
        {
            gzprintf(gzfile, "%.4f", row[order[j]]);
            gzputc(gzfile, ',');
        }
        gzputc(gzfile, '\n');
    }

    gzclose(gzfile);
}

void WriteTreeText(const std::string& treetext, const std::string& outputfilepath)
{
    std::ofstream out(outputfilepath );
    if (!out)
    {
        throw std::runtime_error("Failed to open file for writing tree: " + outputfilepath + "_tree.ph");
    }
    out << treetext << std::endl;
    out.close();
}


void kmer_tree::run( string &inputfilepath, string &normfilepath, string &outputfilepath)
{
    vector<string> headers, header_used;
    loadcontigs(inputfilepath, headers);
    
    vector<vector<normvalue>> norm;
    loadnorm(normfilepath, norm, headers.size());
    
    vector<int> redundant_indices;
    find_identical_rows(norm,redundant_indices);
    
    vector<vector<normvalue>> dm;
    compute_distance_matrix(norm, redundant_indices, dm);
    
    int newindex = 0;
    vector<int> condense_index(headers.size(), 0);
    vector<vector<string>> header_groups(headers.size());
    for (int i = 0; i < headers.size(); ++i)
    {
        condense_index[i] = newindex;
        if (redundant_indices[i] == 0)
        {
            condense_index[i] = newindex;
            header_used.push_back(headers[i]);
            header_groups[newindex++].push_back(headers[i]);
        }
        else
        {
            header_groups[condense_index[redundant_indices[i] - 1]].push_back(headers[i]);
        }
    }
    
    UPGMA umpma (dm, header_used, header_groups);
    auto &tree = umpma.tree;
    tree->reorder(dm);
    std::string treetext = tree->to_newick(header_groups) + ";";
        
    std::vector<int> treeorder = ExtractTreeOrder( treetext, headers);
    
    ReorderFile(inputfilepath,outputfilepath,treeorder);
    
    ReorderMatrix(norm, treeorder, outputfilepath+"_norm.txt.gz");
    
    WriteTreeText(treetext, outputfilepath+"_tree.ph");
    
}

