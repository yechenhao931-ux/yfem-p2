/**
 * KwayPartition.cpp  ──  多级 k-way 划分完整实现
 *
 * 三阶段：
 *   1. CoarsenGraph     (Coarsen.h)
 *   2. InitPartition    (RecursiveBisect，多次取最优)
 *   3. ProjectAndRefine (逐层投影 + KwayFMCut/KwayFMVol)
 */
#include "partition/kwaypartition.hpp"
#include <algorithm>
#include <cstdio>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

// ================================================================
//  ProjectPartition
//  将粗图 c 的 where[] 投影回细图 f：
//    f.where[v] = c.where[f.cmap[v]]
//  然后重建 f.pwgts[]
// ================================================================
static void ProjectPartition(Graph& f, const Graph& c, int nparts) {
    f.where.resize(f.nvtxs);
    for (int v = 0; v < f.nvtxs; ++v)
        f.where[v] = c.where[f.cmap[v]];

    f.pwgts.assign(nparts, 0);
    for (int v = 0; v < f.nvtxs; ++v)
        f.pwgts[f.where[v]] += f.Vwgt(v);
}

// ================================================================
//  InitPartition
//  对最粗图用递归二分做初始划分，重复 initTrials 次取最优。
//  最优标准：平衡解优先，次比切边数。
// ================================================================
static void InitPartition(Graph& gc, int nparts, const KwayOptions& opts,
                          std::mt19937& rng) {
    int    bestCut = std::numeric_limits<int>::max();
    bool   bestBal = false;
    std::vector<int> bestWhere, bestPwgts;

    RBOptions rbo;
    rbo.nparts       = nparts;
    rbo.nTrials      = opts.rbTrials;
    rbo.nFMPasses    = opts.rbFMPasses;
    rbo.ubFactor     = opts.ubFactor;
    rbo.postRefine   = PostRefineMode::None; // 粗图上不做后精化，留给反粗化
    rbo.verbose      = false;

    for (int t = 0; t < opts.initTrials; ++t) {
        rbo.seed = (int)rng();
        RBResult r = RecursiveBisect(gc, rbo);

        bool better = (r.balanced && !bestBal) ||
                      (r.balanced == bestBal && r.totalCut < bestCut);
        if (better) {
            bestCut   = r.totalCut;
            bestBal   = r.balanced;
            bestWhere = gc.where;
            bestPwgts = gc.pwgts;
        }
    }

    gc.where  = bestWhere;
    gc.pwgts  = bestPwgts;
    // 重建 ckrinfo 供后续精化使用
    ComputeCkrinfo(gc, nparts);

    if (opts.verbose)
        std::printf("[Init] coarsest n=%d  cut=%d  %s\n",
                    gc.nvtxs, bestCut, bestBal ? "OK" : "UNBAL");
}

// ================================================================
//  KwayPartition  ──  主入口
// ================================================================
KwayResult KwayPartition(Graph& graph, const KwayOptions& opts) {
    const int K = opts.nparts;
    if (K < 1 || graph.nvtxs == 0) return {};
    if (K == 1) {
        graph.where.assign(graph.nvtxs, 0);
        graph.pwgts.assign(1, graph.tvwgt.empty() ? graph.nvtxs : graph.tvwgt[0]);
        KwayResult r; r.nlevels=0; r.balanced=true; r.partWeights=graph.pwgts;
        return r;
    }

    if (graph.tvwgt.empty()) graph.InitTvwgt();
    std::mt19937 rng(opts.seed);

    // ── 1. 粗化阶段 ───────────────────────────────────────────
    if (opts.verbose)
        std::printf("===== KwayPartition n=%d nparts=%d =====\n",
                    graph.nvtxs, K);

    CoarsenOpts co;
    co.coarseLimit    = opts.coarseLimit;
    co.minCoarseRation = opts.minCoarseRation;
    co.maxLevels      = opts.maxLevels;
    co.seed           = (int)rng();
    co.useUnionFind   = opts.useUnionFind;
    co.coarsenWeightCap = opts.coarsenWeightCap;

    Graph* gc = CoarsenGraph(graph, K, co);

    // 统计粗化层数
    int nlevels = 0;
    for (Graph* p = &graph; p->coarser; p = p->coarser) ++nlevels;

    if (opts.verbose) {
        std::printf("[Coarsen] %d levels  coarsest n=%d\n", nlevels, gc->nvtxs);
    }

    // ── 2. 初始划分（最粗图）─────────────────────────────────
    InitPartition(*gc, K, opts, rng);

    // ── 3. 反粗化 + 精化 ─────────────────────────────────────
    KwayFMOpts fmo;
    fmo.nparts   = K;
    fmo.nIter    = opts.nFMIter;
    fmo.ubFactor = opts.ubFactor;
    fmo.verbose  = false;

    // 从最粗层往回走到原图
    Graph* cur = gc;
    while (cur->finer != nullptr) {
        Graph* fine = cur->finer;

        // 投影分区到细层
        ProjectPartition(*fine, *cur, K);

        if (opts.verbose)
            std::printf("[Uncoarsen] n=%d  (from coarse n=%d)\n",
                        fine->nvtxs, cur->nvtxs);

        // 精化
        if (opts.objective == KwayObjective::EdgeCut) {
            ComputeCkrinfo(*fine, K);
            int gain = opts.useIndepSetRefine ? IndepSetRefineCut(*fine, fmo)
                                              : KwayFMCut(*fine, fmo);
            if (opts.verbose)
                std::printf("  %s gain=%d  cut=%d\n",
                            opts.useIndepSetRefine ? "IndepSet" : "KwayFMCut",
                            gain, fine->mincut);
        } else {
            ComputeCkrinfo(*fine, K);
            ComputeVkrinfo(*fine, K);
            int gain = KwayFMVol(*fine, fmo);
            if (opts.verbose)
                std::printf("  KwayFMVol gain=%d  vol=%d\n", gain, fine->minvol);
        }

        cur = fine;
    }
    // cur 现在是原图（graph）

    // ── 最终精化（额外几轮确保质量）─────────────────────────
    if (opts.objective == KwayObjective::EdgeCut) {
        ComputeCkrinfo(graph, K);
        if (opts.useIndepSetRefine) IndepSetRefineCut(graph, fmo);
        else                        KwayFMCut(graph, fmo);
    } else {
        ComputeCkrinfo(graph, K);
        ComputeVkrinfo(graph, K);
        KwayFMVol(graph, fmo);
    }

    // ── 4. 清理粗化层 ─────────────────────────────────────────
    FreeCoarseGraphs(graph);

    // ── 5. 最终重算 mincut / minvol ───────────────────────────
    {
        int c = 0;
        for (int v = 0; v < graph.nvtxs; ++v)
            for (int ei = graph.xadj[v]; ei < graph.xadj[v+1]; ++ei)
                if (graph.where[v] != graph.where[graph.adjncy[ei]])
                    c += graph.Ewgt(ei);
        graph.mincut = c / 2;
    }

    // ── 6. 负载均衡统计 ───────────────────────────────────────
    real_t ideal = (real_t)graph.tvwgt[0] / K;
    real_t maxImb = 0, sumImb = 0;
    for (int p = 0; p < K; ++p) {
        real_t d = std::abs((real_t)graph.pwgts[p] / ideal - 1.f);
        maxImb = std::max(maxImb, d);
        sumImb += d;
    }

    KwayResult res;
    res.mincut       = graph.mincut;
    res.minvol       = graph.minvol;
    res.maxImbalance = maxImb;
    res.avgImbalance = sumImb / K;
    res.balanced     = (maxImb + 1.f <= opts.ubFactor);
    res.nlevels      = nlevels;
    res.partWeights  = graph.pwgts;

    if (opts.verbose) {
        std::printf("\n===== KwayPartition 完成 =====\n");
        std::printf("粗化层数   : %d\n", nlevels);
        std::printf("总切边数   : %d\n", res.mincut);
        std::printf("最大不平衡 : %.2f%%\n", res.maxImbalance * 100);
        std::printf("平均不平衡 : %.2f%%\n", res.avgImbalance * 100);
        std::printf("满足约束   : %s\n",     res.balanced ? "是" : "否");
        std::printf("理想权重   : %.1f\n",   ideal);
        for (int p = 0; p < K; ++p) {
            real_t r = (real_t)graph.pwgts[p] / ideal;
            std::printf("  [%2d] wgt=%-6d ratio=%.3f%s\n",
                        p, graph.pwgts[p], r,
                        r > opts.ubFactor ? " ←超标" : "");
        }
        std::printf("\n");
    }
    return res;
}