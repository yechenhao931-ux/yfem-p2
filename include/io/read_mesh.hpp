#ifndef READ_MSH_HPP
#define READ_MSH_HPP

#include <vector>
#include <unordered_map>
#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include "../common/type.hpp"
#include "../common/mesh.hpp"

class MFEMMesh10{
public:
    // 单元结构体：存储属性、几何类型和顶点索引
    struct Element
    {
        int attribute;
        int geom_type;
        std::vector<size_t> vertex_indices;
    };

    struct Vertex
    {
        real_t coords[3];
    };
    
    int dimension = 0;

    std::vector<Element> elements;
    std::vector<Element> boundary;

    // vertices[i] 存储第 i 个顶点的坐标数组 [x, y, z]
    std::vector<Vertex> vertices;

    void read_mesh(const char* str);
};


#endif