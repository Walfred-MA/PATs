//
//  KmerTree.hpp
//  KmerTree
//
//  Created by walfred on 6/18/25.
//

#ifndef KmerTree_hpp
#define KmerTree_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstring>
#include <unordered_set>
#include <algorithm>
#include <thread>
#include <atomic>
#include <mutex>
#include <sstream>
#include <math.h>
#include <numeric>
#include <cstdint>
#include <sstream>
#include <iomanip>

#include "KmerStruct.hpp"
#include "fasta.hpp"
#include "gzfile.hpp"

using normvalue = double;

class UPGMANode
{
public:
    std::string name;
    UPGMANode *left = NULL;
    UPGMANode *right = NULL;
    normvalue up_dist = 0.0;
    normvalue down_dist = 0.0;
    int index = -1;
    int nleaves = 1;
    UPGMANode(){};
    UPGMANode(std::string name_, int idx) : name(std::move(name_)), index(idx) {}

    UPGMANode(UPGMANode *l, UPGMANode *r)
        : left(std::move(l)), right(std::move(r))
    {
        nleaves = (left ? left->nleaves : 0) + (right ? right->nleaves : 0);
    }

    std::vector<std::string> leaves() const
    {
        if (!right) return {name};
        auto left_leaves = left->leaves();
        auto right_leaves = right->leaves();
        left_leaves.insert(left_leaves.end(), right_leaves.begin(), right_leaves.end());
        return left_leaves;
    }

    std::string to_newick(const std::vector<std::vector<std::string>>& group_dict) const
    {
        
        
        if (!right)
        {
            if (index >= 0 && index < group_dict.size() && group_dict[index].size() > 1)
            {
                return build_balanced_newick(group_dict[index], up_dist);
            }
            else
            {
                std::ostringstream oss;
                oss << name << ":" << std::fixed << std::setprecision(7) << up_dist;
                return oss.str();
            }
        }
        else
        {
            std::ostringstream oss;
            oss << "(" << left->to_newick(group_dict) << "," << right->to_newick(group_dict)
                << "):" << std::fixed << std::setprecision(6) << up_dist;
            return oss.str();
        }
    }

    std::vector<const UPGMANode*> allleaves() const
    {
        if (!left || !right) {
            return {this};  // 'this' is safe here because it's not null and points to a leaf
        }

        std::vector<const UPGMANode*> l, r;
        if (left)  l = left->allleaves();
        if (right) r = right->allleaves();

        l.insert(l.end(), r.begin(), r.end());
        return l;
    }

    void reorder(const std::vector<std::vector<normvalue>>& dist_matrix, const std::vector<int>* sibling_indexes = nullptr, int sign = 0)
    {
        if (!left || !right) return;

        auto leftleaves = left->allleaves();
        auto rightleaves = right->allleaves();

        std::vector<int> left_indexes, right_indexes;
        for (auto& l : leftleaves) left_indexes.push_back(l->index);
        for (auto& r : rightleaves) right_indexes.push_back(r->index);

        left->reorder(dist_matrix, &right_indexes, 1);
        right->reorder(dist_matrix, &left_indexes, -1);

        if (!sibling_indexes || sign == 0) return;

        normvalue left_mean = 0.0, right_mean = 0.0;
        int count = 0;
        for (int i : sampleleaves(left_indexes))
        {
            for (int j : sampleleaves(*sibling_indexes))
            {
                left_mean += dist_matrix[i][j];
                ++count;
            }
        }
        left_mean /= count;

        count = 0;
        for (int i : sampleleaves(right_indexes))
        {
            for (int j : sampleleaves(*sibling_indexes))
            {
                right_mean += dist_matrix[i][j];
                ++count;
            }
        }
        right_mean /= count;

        if (sign * right_mean > sign * left_mean)
        {
            std::swap(left, right);
        }
    }

private:
    static std::vector<int> sampleleaves(const std::vector<int>& vec)
    {
        if (vec.size() <= 100) return vec;
        std::vector<int> sampled;
        float ratio = 100.0 / vec.size();
        int lastdigit = -1;
        for (size_t i = 0; i < vec.size(); ++i)
        {
            int currdigt = (int)(i * ratio + 0.5);
            if ( currdigt != lastdigit )
            {
                sampled.push_back(vec[i]);
            }
            lastdigit = currdigt;
        }
        return sampled;
    }

    static std::string build_balanced_newick(const std::vector<std::string>& leaves, normvalue up_dist)
    {
        if (leaves.size() == 1)
        {
            std::ostringstream oss;
            oss << leaves[0] << ":" << std::fixed << std::setprecision(7) << up_dist;
            return oss.str();
        }
        else if (leaves.size() == 2)
        {
            std::ostringstream oss;
            oss << "(" << leaves[0] << ":0.0000000," << leaves[1] << ":0.0000000):" << std::fixed << std::setprecision(7) << up_dist;
            return oss.str();
        }
        else
        {
            size_t mid = leaves.size() / 2;
            std::string left = build_balanced_newick(std::vector<std::string>(leaves.begin(), leaves.begin() + mid), 0.0);
            std::string right = build_balanced_newick(std::vector<std::string>(leaves.begin() + mid, leaves.end()), 0.0);
            std::ostringstream oss;
            oss << "(" << left << "," << right << "):" << std::fixed << std::setprecision(7) << up_dist;
            return oss.str();
        }
    }
};



class UPGMA {
public:
    std::vector<std::vector<normvalue>> distances;
    std::vector<std::string> header;
    std::vector<int> exclude;
    std::vector<std::vector<normvalue>> work_matrix;
    UPGMANode *tree;
    std::vector<UPGMANode> nodes;
    int size = 0;
    
    UPGMA(const std::vector<std::vector<normvalue>>& dist_matrix, const std::vector<std::string>& header_, const  std::vector<std::vector<std::string>>& group_dict)
        : distances(dist_matrix), header(header_)
    {
        size = header.size();
        exclude.resize(2 * size, 1);
        for (int i = 0; i < size; ++i) exclude[i] = 0;
        build_tree(dist_matrix, header, group_dict);
    }

private:
    std::pair<int, int> getmindist(int currindex)
    {
        int min_row = -1, min_col = -1;
        normvalue min_val = std::numeric_limits<normvalue>::infinity();

        for (int i = 0; i < currindex; ++i)
        {
            if (exclude[i]) continue;
            for (int j = i+1; j < currindex; ++j)
            {
                if (exclude[j]) continue;
                if (work_matrix[i][j] < min_val)
                {
                    min_val = work_matrix[i][j];
                    min_row = i;
                    min_col = j;
                }
            }
        }
        //cout << "value: " << min_val << " ";
        return {min_row, min_col};
    }

    void build_tree(const std::vector<std::vector<normvalue>>& dist_matrix, const std::vector<std::string>& header, const  std::vector<std::vector<std::string>>& group_dict)
    {
        nodes.resize(2*size);
        for (int i = 0; i < size; ++i)
        {
            nodes[i] = UPGMANode(header[i], i);
            nodes[i].nleaves = group_dict[i].size();
        }

        work_matrix = std::vector<std::vector<normvalue>>(2 * size, std::vector<normvalue>(2 * size, std::numeric_limits<normvalue>::infinity()));

        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                work_matrix[i][j] = dist_matrix[i][j];
            }
        }

        for (int step = 0; step < size - 1 ; ++step)
        {
            auto[min_i, min_j] = getmindist(step+size);
            if (min_i == -1 || min_j == -1) break;

            normvalue dist = work_matrix[min_i][min_j];
            auto &node1 = nodes[min_i];
            auto &node2 = nodes[min_j];

            exclude[min_i] = exclude[min_j] = 1;
            UPGMANode new_node (&node1, &node2);
            new_node.index = -(size+step);
            nodes[step+size] = new_node;
            exclude[step+size] = 0;
            
            //cout <<"adding: " <<min_i << " and " << min_j << " to: "<< step+size << endl;
            node1.up_dist = dist / 2.0 - node1.down_dist;
            node2.up_dist = dist / 2.0 - node2.down_dist;
            new_node.down_dist = dist / 2.0;
            
            update_distance(step+size,{min_i, min_j});
        }

        tree = &nodes[2*size - 2];
    }

    void update_distance(int currindex, std::pair<int, int> least_id)
    {
        int i1 = least_id.first, i2 = least_id.second;
        int n1 = nodes[i1].nleaves, n2 = nodes[i2].nleaves;
        int total = n1 + n2;

        for (int i = 0; i < currindex; ++i)
        {
            if (exclude[i]) continue;
            
            normvalue new_dist1 = work_matrix[i][i1];
            normvalue new_dist2 = work_matrix[i][i2];
            normvalue new_dist = (new_dist1 * n1 + new_dist2 * n2) / total;
            work_matrix[i][currindex] = work_matrix[currindex][i] = new_dist;
                
        }
    }
};



class kmer_tree
{
public:
    
    kmer_tree()
    {};
    ~kmer_tree()
    {};
    
    void readfile(fasta &fastafile);
    
    void loadnorm(string &normfilepath, vector<vector<normvalue>> &data, const int dim);
    
    void loadcontigs(const std::string& inputfile, std::vector<std::string>& names);
    
    void determinetree();
    
    void write();
    
    void run(string &normfile, string &inputfile, string &outputfile);
    
    
};



#endif /* KmerTree_hpp */
