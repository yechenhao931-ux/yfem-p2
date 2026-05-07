#ifndef DUAL_GRAPH_HPP
#define DUAL_GRAPH_HPP

#include <vector>
#include <unordered_map>
#include <algorithm>
#include "../common/mesh.hpp"


struct FaceHash
{
    typedef std::vector<size_t> Face;
    size_t operator()(const Face &f) const
    {
        size_t h = 0;
        for (size_t x : f)
        {
            h ^= std::hash<size_t>()(x) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};

std::vector<std::vector<size_t>> build_dual_graph(const std::vector<std::vector<size_t>>& elements, int dim);

#endif