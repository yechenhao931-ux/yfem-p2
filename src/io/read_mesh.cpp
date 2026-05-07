
#include <io/read_mesh.hpp>

#include <vector>
#include <map>
#include <string>
#include <fstream>


// 根据 MFEM 类型定义获取顶点数量
int GetVerticesCount(int geom_int)
{
    // 将输入的 int 强转为枚举类型进行匹配
    GeometryType geom = static_cast<GeometryType>(geom_int);

    switch (geom) {
        case GeometryType::POINT:       return 1;
        case GeometryType::SEGMENT:     return 2;
        case GeometryType::TRIANGLE:    return 3;
        case GeometryType::SQUARE:      return 4;
        case GeometryType::TETRAHEDRON: return 4;
        case GeometryType::CUBE:        return 8;
        case GeometryType::PRISM:       return 6;
        default:                        return 0; // 未知类型
    }
}

void MFEMMesh10::read_mesh(const char* str){
    std::ifstream file(str);
    if (!file.is_open()){
        std::cerr << std::string("can't open ") + str << std::endl;
        throw std::runtime_error(std::string("can't open ") + str);
    }
    int dim;
    std::string word;

    while (file >> word)
    {
        if (word == "dimension")
        {
            file >> dim;
        }
        else if (word == "elements")
        {
            int num_el;
            file >> num_el;
            elements.reserve(num_el);
            for (int i = 0; i < num_el; ++i)
            {
                Element el;
                file >> el.attribute >> el.geom_type;
                int nv = GetVerticesCount(el.geom_type);
                for (int j = 0; j < nv; ++j)
                {
                    int idx;
                    file >> idx;
                    el.vertex_indices.push_back(idx);
                }
                elements.push_back(el);
            }
        }
        else if (word == "boundary")
        {
            int num_be;
            file >> num_be;
            boundary.reserve(num_be);
            for (int i = 0; i < num_be; ++i)
            {
                Element be;
                file >> be.attribute >> be.geom_type;
                int nv = GetVerticesCount(be.geom_type);
                for (int j = 0; j < nv; ++j)
                {
                    int idx;
                    file >> idx;
                    be.vertex_indices.push_back(idx);
                }
                boundary.push_back(be);
            }
        }
        else if (word == "vertices")
        {
            int num_v, v_dim;
            file >> num_v >> v_dim;
            // 预分配空间：num_v 行，v_dim 列
            vertices.resize(num_v);
            for (int i = 0; i < num_v; ++i)
            {
                for (int d = 0; d < v_dim; ++d)
                {
                    file >> vertices[i].coords[d];
                }
            }
        }
    }
}