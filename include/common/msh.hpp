#ifndef MSH_HPP
#define MSH_HPP
#include <vector>
#include <unordered_map>
#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include "common/type.hpp"

#ifdef USE_GMSH
struct GNode
{
    size_t tag;
    real_t x, y, z;
};

struct GNodeBlock
{
    int32_t entityDim;
    int32_t entityTag;
    int32_t parametri;
    size_t numNodesInBlock;
    std::vector<GNode> nodes;
};

struct GElement
{
    size_t elementTag;
    size_t numNodesinElement;
    std::vector<size_t> nodes;
};

struct GEntityBlock
{
    int32_t entityDim;
    int32_t entityTag;
    int32_t elementType;
    size_t numElementsInBlock;
    std::vector<GElement> elements;
};

struct GEntityBlocks
{
    std::vector<GEntityBlock> blocks;
};


struct GPhysicalName
{
    int dim;
    size_t tag;
    std::string name;

    void print()
    {
        std::cout << dim << " " << tag << " " << name << std::endl;
    }
};

#endif

#endif