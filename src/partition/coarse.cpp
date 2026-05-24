/**
 * Coarsen.cpp  ──  Heavy Edge Matching 粗化完整实现
 *
 * 核心函数：
 *   HEM_Match      : 随机遍历，贪心选最重邻居配对
 *   BuildCoarseGraph: 根据匹配结果构造粗图（合并平行边）
 *   CoarsenGraph   : 递归粗化到停止条件
 */
#include <partition/coarse.hpp>
#include <algorithm>
#include <random>
#include <numeric>
#include <map>
#include <vector>
#include <unordered_map>
#include <climits>
#ifdef _OPENMP
#include <omp.h>
#endif

// ================================================================
//  HEM_Match
//  输出：match[v] = v 配对的超级顶点编号（粗图顶点号）
//        nc = 粗化后顶点总数
// ================================================================
static int HEM_Match(const Graph &g, std::vector<int> &match, int seed)
{
    const int n = g.nvtxs;
    match.assign(n, -1);

    // 随机化访问顺序（减少顺序偏差）
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::mt19937 rng(seed);
    std::shuffle(perm.begin(), perm.end(), rng);

    int nc = 0; // 粗化后顶点计数

    for (int pi = 0; pi < n; ++pi)
    {
        int v = perm[pi];
        if (match[v] != -1)
            continue; // 已匹配

        // 找边权最重的未匹配邻居
        int bestU = -1, bestW = -1;
        for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
        {
            int u = g.adjncy[ei], ew = g.Ewgt(ei);
            if (match[u] == -1 && ew > bestW)
            {
                bestW = ew;
                bestU = u;
            }
        }

        if (bestU != -1)
        {
            // 配对 (v, bestU) → 同一超级顶点 nc
            match[v] = nc;
            match[bestU] = nc;
            ++nc;
        }
        else
        {
            // 自匹配（无未匹配邻居）
            match[v] = nc++;
        }
    }
    return nc;
}

// ================================================================
//  UnionFindMatch  ──  A-2: 带打分的 union-find 粗化（G-kway, DAC'24）
//
//  与 HEM 的区别：HEM 只能两两配对；union-find 可把多点并入同一子集，
//  层数更少。两步：
//    1) 选邻居（可并行）：每个顶点取分数最高的邻居
//         s(u,v) = c·W(u,v) − degree(v)，c = maxdeg+1（边权优先，平局选低度数邻居）
//    2) union-find（串行 DSU，O(n·α)）：沿 selected_nbr 合并连通分量，
//       可选权重上限 maxVwgt 防止合出过重子集
//  输出：match[v] = 连续粗顶点编号(0..nc-1)；返回 nc
//  注：结果由 sel[] 唯一确定（与线程数无关，可复现）
// ================================================================
static int find_root(std::vector<int> &parent, int x)
{
    while (parent[x] != x)
    {
        parent[x] = parent[parent[x]]; // 路径折半
        x = parent[x];
    }
    return x;
}

static int UnionFindMatch(const Graph &g, std::vector<int> &match, int maxVwgt)
{
    const int n = g.nvtxs;
    match.assign(n, -1);

    // c = maxdeg + 1（保证边权项主导度数项）
    int maxdeg = 1;
#pragma omp parallel for reduction(max : maxdeg) if (n > 16384)
    for (int v = 0; v < n; ++v)
    {
        int d = g.xadj[v + 1] - g.xadj[v];
        if (d > maxdeg)
            maxdeg = d;
    }
    const long c = (long)maxdeg + 1;

    // 步骤1：选邻居（并行；每个顶点独立，平局取最先出现的最高分邻居）
    std::vector<int> sel(n, -1);
#pragma omp parallel for schedule(static) if (n > 16384)
    for (int v = 0; v < n; ++v)
    {
        long best = LLONG_MIN;
        int bu = -1;
        for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
        {
            int u = g.adjncy[ei];
            int du = g.xadj[u + 1] - g.xadj[u];
            long s = c * (long)g.Ewgt(ei) - du;
            if (s > best)
            {
                best = s;
                bu = u;
            }
        }
        sel[v] = bu;
    }

    // 步骤2：union-find（串行，确定性；按 sel 边合并，较小根胜出）
    std::vector<int> parent(n), subw(n);
    std::iota(parent.begin(), parent.end(), 0);
    for (int v = 0; v < n; ++v)
        subw[v] = g.Vwgt(v);

    for (int v = 0; v < n; ++v)
    {
        int u = sel[v];
        if (u < 0)
            continue;
        int rv = find_root(parent, v), ru = find_root(parent, u);
        if (rv == ru)
            continue;
        if (maxVwgt > 0 && subw[rv] + subw[ru] > maxVwgt)
            continue; // 权重上限：避免过重子集
        if (ru < rv)
            std::swap(rv, ru); // 较小根作为新根（确定性）
        parent[ru] = rv;
        subw[rv] += subw[ru];
    }

    // 步骤3：把根重标号成连续粗顶点编号
    std::vector<int> label(n, -1);
    int nc = 0;
    for (int v = 0; v < n; ++v)
    {
        int r = find_root(parent, v);
        if (label[r] == -1)
            label[r] = nc++;
        match[v] = label[r];
    }
    return nc;
}

// ================================================================
//  BuildCoarseGraph
//  根据 match[] 构造粗图 coarse：
//    · 超级顶点权重 = 被收缩顶点权重之和
//    · 超级边权重  = 所有被收缩边权之和（去掉自环，合并平行边）
//  同时设置 fine.cmap[v] = match[v]（细图顶点 → 粗图顶点）
// ================================================================
static Graph *BuildCoarseGraph(Graph &fine, const std::vector<int> &match, int nc)
{
    const int n = fine.nvtxs;
    fine.cmap = match; // 保存映射：细图 → 粗图

    Graph *c = new Graph();
    c->nvtxs = nc;
    c->ncon = fine.ncon;
    c->finer = &fine;
    fine.coarser = c;

    // ── 超级顶点权重 ─────────────────────────────────────────
    c->vwgt.assign(nc, 0);
    c->vsize.assign(nc, 0);
    for (int v = 0; v < n; ++v)
    {
        int cv = match[v];
        c->vwgt[cv] += fine.Vwgt(v);
        c->vsize[cv] += fine.Vsize(v);
    }

    // ── 构造粗图 CSR ─────────────────────────────────────────
    // 对每个粗顶点，扫描所有细顶点的邻边，累加到粗边
    // 使用 hash map 合并平行边（临时，每个粗顶点处理后清空）
    std::vector<int> xadj(nc + 1, 0);
    // 先统计每个粗顶点的邻居数（粗图度数）
    // 方法：双遍扫描
    //   Pass1: 对每条细边 (v,u)，若 cv != cu，标记 cv 有邻居 cu
    //   用 marker 数组避免重复计数

    // Pass1：统计粗图中每个超级顶点的邻居分区数（用于 xadj）
    // 先建立粗顶点 → 细顶点列表
    std::vector<std::vector<int>> c2f(nc);
    for (int v = 0; v < n; ++v)
        c2f[match[v]].push_back(v);

    // 并行化（OpenMP）：各粗顶点 cv 的统计相互独立，仅写自己的 xadj[cv+1]。
    //   marker 用"!=cv"技巧去重；每线程独立的 marker 避免竞争，且因全局 cv 互异，
    //   该技巧在每个线程内部依然成立。结果与串行逐字节一致。
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    std::vector<std::vector<int>> markerTL(nthreads, std::vector<int>(nc, -1));

#pragma omp parallel for schedule(dynamic, 256) if (nc > 16384)
    for (int cv = 0; cv < nc; ++cv)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        std::vector<int> &marker = markerTL[tid];
        int deg = 0;
        for (int fv : c2f[cv])
        {
            for (int ei = fine.xadj[fv]; ei < fine.xadj[fv + 1]; ++ei)
            {
                int cu = match[fine.adjncy[ei]];
                if (cu == cv)
                    continue; // 自环
                if (marker[cu] != cv)
                {
                    marker[cu] = cv;
                    ++deg;
                }
            }
        }
        xadj[cv + 1] = deg;
    }
    // 前缀和
    for (int cv = 0; cv < nc; ++cv)
        xadj[cv + 1] += xadj[cv];
    c->xadj = xadj;
    c->adjncy.resize(xadj[nc]);
    c->adjwgt.resize(xadj[nc], 0);
    c->nedges = xadj[nc];

    // Pass2：填充边列表和边权（并行：各 cv 仅写自己的 [xadj[cv],xadj[cv+1]) 区间）
    std::vector<int> pos(nc);
    for (int cv = 0; cv < nc; ++cv)
        pos[cv] = xadj[cv];

    // 每线程独立的 edgeIdx（记录 cv→cu 在 adjncy 中的槽位；每个 cv 处理完复位）
    std::vector<std::vector<int>> edgeIdxTL(nthreads, std::vector<int>(nc, -1));

#pragma omp parallel for schedule(dynamic, 256) if (nc > 16384)
    for (int cv = 0; cv < nc; ++cv)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        std::vector<int> &edgeIdx = edgeIdxTL[tid];
        std::vector<int> curNbrs;
        for (int fv : c2f[cv])
        {
            for (int ei = fine.xadj[fv]; ei < fine.xadj[fv + 1]; ++ei)
            {
                int fu = fine.adjncy[ei];
                int cu = match[fu];
                int ew = fine.Ewgt(ei);
                if (cu == cv)
                    continue;

                if (edgeIdx[cu] == -1)
                {
                    // 首次遇到 cu：分配槽位
                    int slot = pos[cv]++;
                    c->adjncy[slot] = cu;
                    c->adjwgt[slot] = ew;
                    edgeIdx[cu] = slot;
                    curNbrs.push_back(cu);
                }
                else
                {
                    c->adjwgt[edgeIdx[cu]] += ew;
                }
            }
        }
        // 复位 edgeIdx，供本线程下一个 cv 使用
        for (int cu : curNbrs)
            edgeIdx[cu] = -1;
    }

    c->InitTvwgt();

    // set coord
    c->coordinates.assign(c->nvtxs,{0,0,0});

    for (int v = 0; v < n; ++v)
    {
        int cv = match[v];
        auto& cc = c->coordinates[cv];
        const auto& fv = fine.coordinates[v];
        auto k = fine.vwgt[v];
        cc.x += fv.x * k;
        cc.y += fv.y * k;
        cc.z += fv.z * k;
    }

    for (int cv = 0; cv < c->nvtxs; ++cv)
    {
        auto& cc = c->coordinates[cv];
        cc.x /= c->vwgt[cv];
        cc.y /= c->vwgt[cv];
        cc.z /= c->vwgt[cv];

        //printf("%f %f %f\n",cc.x,cc.y,cc.z);
    }


    return c;
}

// ================================================================
//  CoarsenGraph  ──  递归粗化到停止条件
// ================================================================
Graph *CoarsenGraph(Graph &fine, int nparts, const CoarsenOpts &opts)
{
    // 粗化停止阈值：至少保留 coarseLimit * nparts 个顶点
    int stopAt = std::max(opts.coarseLimit * nparts, 40);

    Graph *cur = &fine;
    int seed = opts.seed;

    // union-find 的超级顶点权重上限（总权重在粗化中不变，按 fine 计算一次）
    int totalVwgt = fine.tvwgt.empty() ? fine.nvtxs : fine.tvwgt[0];
    int maxVwgt = (opts.useUnionFind && opts.coarsenWeightCap > 0)
                      ? (int)((real_t)totalVwgt / nparts * opts.coarsenWeightCap + 0.5f)
                      : 0;

    for (int lv = 0; lv < opts.maxLevels; ++lv)
    {
        if (cur->nvtxs <= stopAt)
            break;

        // 匹配：union-find（A-2）或经典 HEM
        std::vector<int> match;
        int nc = opts.useUnionFind ? UnionFindMatch(*cur, match, maxVwgt)
                                   : HEM_Match(*cur, match, seed++);

        // 粗化效果不足：若顶点缩减比 < minCoarseRation 才停止
        // (nc/cur->nvtxs >= minCoarseRation 说明收缩很少)
        real_t ratio = (real_t)nc / cur->nvtxs;
        if (lv > 0 && ratio >= opts.minCoarseRation)
            break;

        // 构造粗图
        Graph *coarse = BuildCoarseGraph(*cur, match, nc);
        cur = coarse;
    }
    return cur; // 返回最粗层
}

// ================================================================
//  FreeCoarseGraphs
// ================================================================
void FreeCoarseGraphs(Graph &g)
{
    Graph *cur = g.coarser;
    g.coarser = nullptr;
    while (cur)
    {
        Graph *next = cur->coarser;
        cur->finer = nullptr;
        cur->coarser = nullptr;
        delete cur;
        cur = next;
    }
}
