#ifndef GRAPH_HPP
#define GRAPH_HPP

#include "dual_graph.hpp"
#include <vector>
#include <memory>
#include <iostream>
#include <iomanip>
#include "../common/type.hpp"
#ifdef USE_GMSH
#include "../io/read_msh.hpp"
#endif
#include "../io/read_mesh.hpp"
typedef struct Coord
{
    real_t x, y, z;
} Coord;

struct Cnbr
{
    int pid = -1; // 相邻分区的ID
    int ed = 0;   // 连接到该分区的边数 (edge count)
};
struct Ckrinfo
{
    int id = 0;    //  internal degree
    int ed = 0;    //  external degree
    int nnbrs = 0; // 相邻分区的数量
    int inbr = -1; // 相邻分区数组的起始索引
};
struct Vnbr
{
    int pid = -1; // 相邻分区的ID
    int ed = 0;   // 连接到该分区的边数
    int gv = 0;   // 增益值 (gain value)
};
struct Vkrinfo
{
    int nid = 0;   // node internal degree
    int ned = 0;   // node external degree
    int gv = 0;    // 增益值
    int nnbrs = 0; // 相邻分区数量
    int inbr = -1; // 相邻分区数组起始索引
};
int compute_element_dof(int element_type, int p);
Coord compute_centroid(const std::vector<Coord> &nodes, const std::vector<size_t> &elem);

class Graph
{
public:
    /* ----基础图结构---- */
    int nvtxs = 0;  // number of vertices
    int nedges = 0; // number of edges * 2

    int ncon = 1; // number if constraints (default 1)

    std::vector<int> xadj;   // CSR row ptr
    std::vector<int> adjncy; // CSR col
    std::vector<int> adjwgt; // edge weights

    std::vector<int> vwgt;  // vertex weights
    std::vector<int> vsize; // volume weights (optional)

    std::vector<int> tvwgt;       // total weight per constraint
    std::vector<real_t> invtvwgt; // inverse total weight

    /* ----基于物理结构----*/
    std::vector<Coord> coordinates; // coordinate of vertexs

    /* ----多级粗化---- */

    std::vector<int> cmap;    // fine -> coarse map
    Graph *coarser = nullptr; // pointer to coarser graph
    Graph *finer = nullptr;   // pointer to finer graph

    /*----分区结果---- */

    std::vector<int> where; // partition id per vertex
    std::vector<int> pwgts; // partition weights

    int mincut = 0;
    int minvol = 0;

    /* ----边界信息----*/
    int nbnd = 0;
    std::vector<int> bndptr; // boundary pointer
    std::vector<int> bndind; // boundary vertices

    /*------------------ Bisection Refinement ------------------*/

    std::vector<int> id; // internal degree
    std::vector<int> ed; // external degree

    /*------------------ K-way Refinement ------------------*/

    std::vector<Ckrinfo> ckrinfo;
    std::vector<Cnbr> cnbrPool;
    void AllocateCkrinfo(int hint = 6)
    {
        ckrinfo.assign(nvtxs, {});
        cnbrPool.clear();
        cnbrPool.reserve((size_t)nvtxs * hint);
    }

    std::vector<Vkrinfo> vkrinfo;
    std::vector<Vnbr> vnbrPool;
    void AllocateVkrinfo(int hint = 6)
    {
        vkrinfo.assign(nvtxs, {});
        vnbrPool.clear();
        vnbrPool.reserve((size_t)nvtxs * hint);
    }

    int Vwgt(int v) const { return vwgt.empty() ? 1 : vwgt[v]; }
    int Vsize(int v) const { return vsize.empty() ? 1 : vsize[v]; }
    int Ewgt(int ei) const { return adjwgt.empty() ? 1 : adjwgt[ei]; }
    void InitTvwgt()
    {
        tvwgt.assign(ncon, 0);
        for (int v = 0; v < nvtxs; ++v)
        {
            tvwgt[0] += Vwgt(v);
        }

        invtvwgt.resize(ncon);
        for (int c = 0; c < ncon; ++c)
        {
            invtvwgt[c] = tvwgt[c] > 0 ? 1.0 / tvwgt[c] : 0.0;
        }
    }
    /* ----Function----*/
    Graph() = default;

    Graph(const MFEMMesh10 &mfem_mesh, int dim);
#ifdef USE_GMSH
    Graph(const Gmsh41 &gmsh, int dim);
#endif
    void printDenseMatrix() const;
    ~Graph() = default;
};

#endif