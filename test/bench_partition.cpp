/**
 * bench_partition.cpp  ──  对比基线 vs 改进后的划分质量
 *
 * 评估指标：
 *   - mincut       : 切边数（越小越好）
 *   - imbalance    : 最大不平衡比（越小越好）
 *   - components   : 各分区的连通子图数量（理想值 1）
 *                    平均值 = sum(components_in_part_p) / nparts
 *
 * 用法：
 *   bench_partition -m <mesh> -p <nparts> -s <seed>
 */
#include "graph/graph.hpp"
#include "partition/geokwaypartition.hpp"
#include "partition/connectrepair.hpp"
#include "common/optparser.hpp"
#include "io/save.hpp"

#include <chrono>
#include <cstdio>
#include <queue>
#include <vector>

// 计算 (孤岛分区数, 总连通子图数)
static std::pair<int,int> CountComponents(const Graph& g, int nparts)
{
    const int n = g.nvtxs;
    std::vector<bool> visited(n, false);
    int totalComps = 0;
    int isolatedParts = 0;
    std::vector<int> compCountPerPart(nparts, 0);

    for (int v = 0; v < n; ++v) {
        if (visited[v]) continue;
        int part = g.where[v];
        std::queue<int> q;
        q.push(v); visited[v] = true;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int ei = g.xadj[u]; ei < g.xadj[u + 1]; ++ei) {
                int nb = g.adjncy[ei];
                if (!visited[nb] && g.where[nb] == part) {
                    visited[nb] = true; q.push(nb);
                }
            }
        }
        compCountPerPart[part]++;
        totalComps++;
    }
    for (int p = 0; p < nparts; ++p)
        if (compCountPerPart[p] > 1) isolatedParts++;
    return {isolatedParts, totalComps};
}

static void RunOne(const char* tag, MFEMMesh10& mmesh, int nparts, int seed,
                   const GeoKwayOptions& o)
{
    Graph g(mmesh, 3);
    auto t0 = std::chrono::steady_clock::now();
    GeoKwayResult r = GeoKwayPartition(g, o);
    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    auto [isoParts, totalComps] = CountComponents(g, nparts);
    int idealComps = nparts;
    std::printf("[%-20s] cut=%-5d  vol=%-5d  imb=%5.2f%%  comps=%d/%d  isoParts=%d  time=%.1fms\n",
                tag, r.mincut, r.minvol, r.maxImbalance * 100,
                totalComps, idealComps, isoParts, ms);
}

int main(int argc, char* argv[])
{
    const char* mesh_file = "../resource/box.mesh";
    int nparts = 8;
    int seed = 42;
    real_t alpha = 0.7;
    int useVol = 0; // 1 → 以通信量(vol)为优化目标，0 → 以切边(cut)为目标

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
    args.AddOption(&nparts, "-p", "--part", "Number of partitions.");
    args.AddOption(&seed, "-s", "--seed", "Random seed.");
    args.AddOption(&alpha, "-a", "--alpha", "Topology weight.");
    args.AddOption(&useVol, "-v", "--vol", "Objective: 1=communication volume, 0=edge cut.");
    args.Parse();
    if (!args.Good()) { args.PrintUsage(std::cout); return 1; }

    MFEMMesh10 mmesh;
    mmesh.read_mesh(mesh_file);

    std::printf("\n=== Mesh=%s  parts=%d  seed=%d  alpha=%.2f  obj=%s ===\n",
                mesh_file, nparts, seed, alpha, useVol ? "vol" : "cut");

    // 配置 1：原始基线（关闭所有改进）
    GeoKwayOptions baseline;
    baseline.nparts = nparts;
    baseline.seed = seed;
    baseline.alpha = alpha;
    baseline.useVolume = (useVol != 0);
    baseline.verbose = false;
    baseline.enforceConnectivity = false;
    baseline.coarsenWeightCap = 0;        // 禁用 HEM 权重上限
    baseline.kmeansBfsSeed = false;        // 禁用 BFS 种子
    baseline.kmeansTopoBalanceFix = false; // 用旧的 BalanceRepair（纯几何）
    baseline.kmeansEnforceConn = false;    // 禁用 K-means 内部连通性修复
    RunOne("baseline", mmesh, nparts, seed, baseline);

    // 配置 2：基线 + BFS 种子初始化
    GeoKwayOptions cfgBfs = baseline;
    cfgBfs.kmeansBfsSeed = true;
    RunOne("+ bfsSeed", mmesh, nparts, seed, cfgBfs);

    // 配置 3：+ 拓扑感知 BalanceFix
    GeoKwayOptions cfgTopo = cfgBfs;
    cfgTopo.kmeansTopoBalanceFix = true;
    RunOne("+ topoBalanceFix", mmesh, nparts, seed, cfgTopo);

    // 配置 4：+ HEM 权重上限
    GeoKwayOptions cfgHEM = cfgTopo;
    cfgHEM.coarsenWeightCap = 1.5f;
    RunOne("+ hemWeightCap", mmesh, nparts, seed, cfgHEM);

    // 配置 5：+ 连通性修复（K-means 内 + 顶层）= 默认
    GeoKwayOptions cfgAll = cfgHEM;
    cfgAll.kmeansEnforceConn = true;
    cfgAll.enforceConnectivity = true;
    RunOne("all-on (default)", mmesh, nparts, seed, cfgAll);

    return 0;
}
