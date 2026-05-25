/**
 * metis_test.cpp  ──  原生 METIS baseline 与 GeoKwayPartition 对比
 *
 * 在完全相同的 Graph 上分别运行：
 *   - METIS Kway (cut)        : METIS_PartGraphKway, 最小化切边
 *   - METIS Kway (vol)        : METIS_PartGraphKway, 最小化通信量
 *   - METIS Recursive (cut)   : METIS_PartGraphRecursive
 *   - GeoKwayPartition        : 本项目几何引导多级 k-way
 *
 * 指标（口径与 GeoKwayResult 一致）：
 *   cut  切边数 | vol 通信量 | imb 最大不平衡 | comps 连通子图数 | time
 *
 * 用法：
 *   metis_test -m <mesh> -p <nparts> -s <seed>
 */
#include "graph/graph.hpp"
#include "partition/metispart.hpp"
#include "partition/geokwaypartition.hpp"
#include "common/optparser.hpp"

#include <chrono>
#include <cstdio>
#include <functional>
#include <queue>
#include <vector>

// 返回 (孤岛分区数, 总连通子图数)；理想总数 = nparts。
static std::pair<int, int> CountComponents(const Graph &g, int nparts)
{
    const int n = g.nvtxs;
    std::vector<bool> visited(n, false);
    int totalComps = 0, isolatedParts = 0;
    std::vector<int> compCountPerPart(nparts, 0);

    for (int v = 0; v < n; ++v)
    {
        if (visited[v])
            continue;
        int part = g.where[v];
        std::queue<int> q;
        q.push(v);
        visited[v] = true;
        while (!q.empty())
        {
            int u = q.front();
            q.pop();
            for (int ei = g.xadj[u]; ei < g.xadj[u + 1]; ++ei)
            {
                int nb = g.adjncy[ei];
                if (!visited[nb] && g.where[nb] == part)
                {
                    visited[nb] = true;
                    q.push(nb);
                }
            }
        }
        compCountPerPart[part]++;
        totalComps++;
    }
    for (int p = 0; p < nparts; ++p)
        if (compCountPerPart[p] > 1)
            isolatedParts++;
    return {isolatedParts, totalComps};
}

// 每次都在新建的 Graph 上运行（划分会改写 where/粗化指针等）。
static void RunOne(const char *tag, MFEMMesh10 &mmesh, int dim, int nparts,
                   const std::function<GeoKwayResult(Graph &)> &fn)
{
    Graph g(mmesh, dim);
    auto t0 = std::chrono::steady_clock::now();
    GeoKwayResult r = fn(g);
    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    auto [isoParts, totalComps] = CountComponents(g, nparts);
    std::printf("[%-22s] cut=%-6d vol=%-6d imb=%5.2f%%  bal=%s  comps=%d/%d  iso=%d  time=%.1fms\n",
                tag, r.mincut, r.minvol, r.maxImbalance * 100,
                r.balanced ? "Y" : "N", totalComps, nparts, isoParts, ms);
}

int main(int argc, char *argv[])
{
    const char *mesh_file = "../resource/box.mesh";
    int nparts = 8;
    int seed = 42;
    real_t alpha = 0.7;

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
    args.AddOption(&nparts, "-p", "--part", "Number of partitions.");
    args.AddOption(&seed, "-s", "--seed", "Random seed.");
    args.AddOption(&alpha, "-a", "--alpha", "GeoKway topology weight.");
    args.Parse();
    if (!args.Good())
    {
        args.PrintUsage(std::cout);
        return 1;
    }

    MFEMMesh10 mmesh;
    mmesh.read_mesh(mesh_file);
    int dim = mmesh.dimension > 0 ? mmesh.dimension : 3;

    std::printf("\n=== Mesh=%s  parts=%d  seed=%d  alpha=%.2f  metis=%s ===\n",
                mesh_file, nparts, seed, alpha,
                MetisAvailable() ? "available" : "NOT COMPILED");

    if (!MetisAvailable())
        std::printf("(提示：用 cmake -DENABLE_METIS=ON 重新配置以启用 METIS 对照)\n");

    // ── METIS Kway，最小化切边 ────────────────────────────────────
    RunOne("METIS Kway (cut)", mmesh, dim, nparts, [&](Graph &g) {
        MetisOptions o;
        o.nparts = nparts;
        o.seed = seed;
        o.algo = MetisAlgo::Kway;
        o.objective = MetisObjective::EdgeCut;
        return MetisPartition(g, o);
    });

    // ── METIS Kway，最小化通信量 ──────────────────────────────────
    RunOne("METIS Kway (vol)", mmesh, dim, nparts, [&](Graph &g) {
        MetisOptions o;
        o.nparts = nparts;
        o.seed = seed;
        o.algo = MetisAlgo::Kway;
        o.objective = MetisObjective::Volume;
        return MetisPartition(g, o);
    });

    // ── METIS 递归二分 ────────────────────────────────────────────
    RunOne("METIS Recursive (cut)", mmesh, dim, nparts, [&](Graph &g) {
        MetisOptions o;
        o.nparts = nparts;
        o.seed = seed;
        o.algo = MetisAlgo::RecursiveBisection;
        return MetisPartition(g, o);
    });

    // ── METIS Kway + 强制连通 ─────────────────────────────────────
    RunOne("METIS Kway (contig)", mmesh, dim, nparts, [&](Graph &g) {
        MetisOptions o;
        o.nparts = nparts;
        o.seed = seed;
        o.algo = MetisAlgo::Kway;
        o.contig = true;
        return MetisPartition(g, o);
    });

    // ── 本项目几何引导 k-way ──────────────────────────────────────
    RunOne("GeoKwayPartition", mmesh, dim, nparts, [&](Graph &g) {
        GeoKwayOptions o;
        o.nparts = nparts;
        o.seed = seed;
        o.alpha = alpha;
        return GeoKwayPartition(g, o);
    });

    return 0;
}
