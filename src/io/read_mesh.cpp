
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
    std::string word;

    while (file >> word)
    {
        if (word == "dimension")
        {
            file >> dimension;
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
            int num_v;
            file >> num_v;
            // resize 会对坐标做值初始化（清零），缺失的维度保持为 0
            vertices.resize(num_v);

            // 探测下一个非空白字符：
            //   数字/正负号 -> 普通网格，紧跟着 v_dim 与坐标
            //   否则（如 "nodes"）-> 高阶网格，坐标由后面的 nodes GridFunction 给出
            file >> std::ws;
            int next_ch = file.peek();
            if (next_ch == '+' || next_ch == '-' ||
                (next_ch >= '0' && next_ch <= '9'))
            {
                int v_dim;
                file >> v_dim;
                for (int i = 0; i < num_v; ++i)
                {
                    for (int d = 0; d < v_dim; ++d)
                    {
                        file >> vertices[i].coords[d];
                    }
                }
            }
            // 否则坐标在下面的 "nodes" 分支中填充
        }
        else if (word == "nodes")
        {
            // 高阶（曲边）网格：节点坐标以 GridFunction 形式存储
            //   FiniteElementSpace
            //   FiniteElementCollection: <name>
            //   VDim: <vd>
            //   Ordering: <ord>      (0 = byNODES, 1 = byVDIM)
            //   <values...>
            // MFEM 约定前 nv 个自由度即为网格顶点（顶点自由度排在最前），
            // 据此恢复每个顶点的坐标。
            int vd = 3;
            int ordering = 0;
            std::string tok;
            while (file >> tok)
            {
                if (tok == "VDim:")
                {
                    file >> vd;
                }
                else if (tok == "Ordering:")
                {
                    file >> ordering;
                }
                else if (tok == "FiniteElementSpace")
                {
                    // 仅是分节标题
                }
                else if (tok == "FiniteElementCollection:")
                {
                    std::string name;
                    std::getline(file, name); // 跳过该行剩余的集合名称
                }
                else
                {
                    // 到达第一个数值，停止解析头部（tok 即第一个值）
                    break;
                }
            }

            std::vector<real_t> vals;
            if (!tok.empty())
            {
                try
                {
                    vals.push_back(static_cast<real_t>(std::stod(tok)));
                }
                catch (...)
                {
                }
            }
            real_t v;
            while (file >> v)
            {
                vals.push_back(v);
            }

            if (vd > 0 && !vals.empty())
            {
                size_t ndof = vals.size() / static_cast<size_t>(vd);
                for (size_t i = 0; i < vertices.size(); ++i)
                {
                    for (int d = 0; d < vd && d < 3; ++d)
                    {
                        size_t idx = (ordering == 0)
                                         ? (static_cast<size_t>(d) * ndof + i) // byNODES
                                         : (i * static_cast<size_t>(vd) + d);  // byVDIM
                        if (idx < vals.size())
                        {
                            vertices[i].coords[d] = vals[idx];
                        }
                    }
                }
            }
        }
    }
}