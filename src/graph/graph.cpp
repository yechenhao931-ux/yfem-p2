#include <graph/graph.hpp>

int compute_element_dof(int element_type, int p)
{
    switch (element_type)
    {
    case 2: // triangle
        return (p + 1) * (p + 2) / 2;

    case 3: // quad
        return (p + 1) * (p + 1);

    case 4: // tetra
        return (p + 1) * (p + 2) * (p + 3) / 6;

    case 5: // hex
        return (p + 1) * (p + 1) * (p + 1);

    default:
        throw std::runtime_error("Unsupported element type");
    }
}

Graph::Graph(const MFEMMesh10 &mfem_mesh, int dim)
{
    std::vector<std::vector<size_t>> basic_elements;
    basic_elements.resize(mfem_mesh.elements.size());

    for (int i = 0; i < mfem_mesh.elements.size(); i++)
    {
        if (mfem_mesh.elements[i].geom_type == 4)
        {
            basic_elements[i] = mfem_mesh.elements[i].vertex_indices;
        }
        else
        {
            assert(0);
        }
    }
    auto dual_graph = build_dual_graph(basic_elements, dim);

    nvtxs = basic_elements.size();

    nedges = 0;

    ncon = 1;
    // 计算总边数（每条边会被计算两次）
    for (const auto &neighbors : dual_graph)
    {
        nedges += static_cast<int>(neighbors.size());
    }

    // 分配内存
    xadj.resize(nvtxs + 1);
    adjncy.resize(nedges);
    adjwgt.resize(nedges, 1);

    // 构建 CSR (xadj,adjncy)
    xadj[0] = 0;
    for (int i = 0; i < nvtxs; ++i)
    {
        auto &neighbors = dual_graph[i];
        int start_idx = xadj[i];
        std::sort(neighbors.begin(), neighbors.end());
        int degree = static_cast<int>(neighbors.size());

        // 复制邻接顶点
        for (int_t j = 0; j < degree; ++j)
        {
            adjncy[start_idx + j] = static_cast<int>(neighbors[j]);
        }

        // 设置下一行的起始位置
        xadj[i + 1] = xadj[i] + degree;
    }

    // adjwgt
    adjwgt.assign(nedges, 1);

    // 可选：分配权重数组（初始化为1）
    vwgt.resize(nvtxs * ncon, 1);

    // 设置 Cood (目前只有一种类型的有限元)
    coordinates.resize(nvtxs);
    for (int i = 0; i < nvtxs; i++)
    {
        const auto &e = mfem_mesh.elements[i];
        real_t x = 0, y = 0, z = 0;
        int n = 0;
        for (auto &c : e.vertex_indices)
        {
            x += mfem_mesh.vertices[c].coords[0];
            y += mfem_mesh.vertices[c].coords[1];
            z += mfem_mesh.vertices[c].coords[2];
            n++;
        }

        coordinates[i].x = x / n;
        coordinates[i].y = y / n;
        coordinates[i].z = z / n;
        // printf("%d:%f %f %f\n",e.elementTag,coordinates[i].x,coordinates[i].y,coordinates[i].z);
    }
}

#ifdef USE_GMSH
Graph::Graph(const Gmsh41 &gmsh, int dim)
{
    auto dual_graph = build_dual_graph(gmsh.basic_elements, dim);

    /*初始化基础图*/

    // nvtxs nedges ncon
    nvtxs = gmsh.basic_elements.size();

    nedges = 0;

    ncon = 1;
    // 计算总边数（每条边会被计算两次）
    for (const auto &neighbors : dual_graph)
    {
        nedges += static_cast<int>(neighbors.size());
    }

    // 分配内存
    xadj.resize(nvtxs + 1);
    adjncy.resize(nedges);
    adjwgt.resize(nedges, 1);

    // 构建 CSR (xadj,adjncy)
    xadj[0] = 0;
    for (int i = 0; i < nvtxs; ++i)
    {
        auto &neighbors = dual_graph[i];
        int start_idx = xadj[i];
        std::sort(neighbors.begin(), neighbors.end());
        int degree = static_cast<int>(neighbors.size());

        // 复制邻接顶点
        for (int_t j = 0; j < degree; ++j)
        {
            adjncy[start_idx + j] = static_cast<int>(neighbors[j]);
        }

        // 设置下一行的起始位置
        xadj[i + 1] = xadj[i] + degree;
    }

    // adjwgt
    adjwgt.assign(nedges, 1);

    // 可选：分配权重数组（初始化为1）
    vwgt.resize(nvtxs * ncon, 1);

    // 设置 Cood (目前只有一种类型的有限元)

    auto &basic = gmsh.entity_blocks.back();
    int basic_e_size = basic.numElementsInBlock;
    coordinates.resize(basic_e_size);
    for (int i = 0; i < basic_e_size; i++)
    {
        const auto &e = basic.elements[i];
        real_t x = 0, y = 0, z = 0;
        int n = 0;
        for (auto &c : e.nodes)
        {
            x += gmsh.nodes[c].x;
            y += gmsh.nodes[c].y;
            z += gmsh.nodes[c].z;
            n++;
        }

        coordinates[i].x = x / n;
        coordinates[i].y = y / n;
        coordinates[i].z = z / n;
        // printf("%d:%f %f %f\n",e.elementTag,coordinates[i].x,coordinates[i].y,coordinates[i].z);
    }
}
#endif

void Graph::printDenseMatrix() const
{
    std::cout << "\n========== Dense Matrix Representation ==========" << std::endl;

    // 创建稠密矩阵
    std::vector<std::vector<int>> dense(nvtxs, std::vector<int>(nvtxs, 0));

    for (int_t i = 0; i < nvtxs; ++i)
    {
        dense[i][i] = 0; // 对角线标记为2

        int_t start = xadj[i];
        int_t end = xadj[i + 1];

        for (int_t idx = start; idx < end; ++idx)
        {
            int_t j = adjncy[idx];
            dense[i][j] = 1;
        }
    }

    // 打印矩阵
    std::cout << "\n    ";
    for (int_t j = 0; j < nvtxs; ++j)
    {
        std::cout << std::setw(2) << j << " ";
    }
    std::cout << "\n   ";
    for (int_t j = 0; j < nvtxs; ++j)
    {
        std::cout << "---";
    }
    std::cout << std::endl;

    for (int_t i = 0; i < nvtxs; ++i)
    {
        std::cout << std::setw(2) << i << " |";
        for (int_t j = 0; j < nvtxs; ++j)
        {
            std::cout << std::setw(2) << dense[i][j] << " ";
        }
        std::cout << std::endl;
    }
}