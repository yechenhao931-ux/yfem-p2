#ifndef READ_MSH_HPP
#define READ_MSH_HPP


#ifdef USE_GMSH
#include <vector>
#include <unordered_map>
#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include "../common/type.hpp"
#include "../common/msh.hpp"


enum ElementType{
    LINE = 1,
    TRIANGLE = 2,
    QUADRANGLE = 3,
    TRTRAHEDRON = 4,
    HEXAHEDRON = 5,
    PRISM = 6,
    PYRAMID = 7,
    POINT = 15,
};
int num_nodes_per_element(ElementType elementType);




/* ===============================
   Gmsh reader class
   =============================== */
class Gmsh41
{
public:
    std::vector<size_t> vmap;
    std::vector<GPhysicalName> physical_names;
    std::vector<GNodeBlock> node_blocks;
    std::vector<GNode> nodes;
    std::vector<GEntityBlock> entity_blocks;
    std::vector<std::vector<size_t>> basic_elements;

    void read_gmsh(const char* str);
};
#endif

#endif