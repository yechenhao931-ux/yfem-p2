#include <partition/geokwaypartition.hpp>
#include <partition/kwaypartition.hpp>
#include <partition/connectrepair.hpp>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

// ================================================================
//  [创新1] 几何感知 HEM 匹配（带权重平衡约束）
//  匹配分数 = λ × ewgt + (1-λ) × ewgt_max / (1 + dist²(u,v) / dref²)
//  使用平方距离避免开方；添加超级顶点权重上限：
//    若 vwgt[v] + vwgt[u] > maxVwgt 则禁止匹配，
//    避免出现极端不平衡的超级顶点。
//  结果：优先合并边权重且空间近邻的顶点对，且粗图各顶点权重均衡。
// ================================================================

static int GeoHEM_Match(const Graph &g, std::vector<int> &match, int seed,
                        real_t lambda, int maxVwgt)
{
    const int n = g.nvtxs;
    match.assign(n, -1);

    // 计算参考距离的平方 dref²（采样平均边长平方）
    // 优势：避免每次匹配做 sqrt
    double sumD2 = 0;
    int cnt = 0;
    int sampleStop = std::min(n, 2000);
    for (int v = 0; v < sampleStop; ++v)
    {
        for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
        {
            int u = g.adjncy[ei];
            real_t dx = g.coordinates[v].x - g.coordinates[u].x;
            real_t dy = g.coordinates[v].y - g.coordinates[u].y;
            real_t dz = g.coordinates[v].z - g.coordinates[u].z;
            sumD2 += dx * dx + dy * dy + dz * dz;
            ++cnt;
        }
    }

    real_t dref2 = (cnt > 0) ? (real_t)(sumD2 / cnt) : 1.;
    if (dref2 < 1e-20) dref2 = 1.;

    // 边权上界（用于归一化)
    int maxEw = 1;
    for (int ei = 0; ei < (int)g.adjwgt.size(); ++ei)
    {
        maxEw = std::max(maxEw, g.adjwgt[ei]);
    }
    // 随机访问顺序
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::mt19937 rng(seed);
    std::shuffle(perm.begin(), perm.end(), rng);

    int nc = 0;
    for (int pi = 0; pi < n; ++pi)
    {
        int v = perm[pi];
        if (match[v] != -1)
            continue;

        real_t bestScore = -1.0;
        int bestU = -1;
        int wv = g.Vwgt(v);

        for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
        {
            int u = g.adjncy[ei];
            if (match[u] != -1)
                continue;
            // 权重上限：避免合并出过重超级顶点
            if (maxVwgt > 0 && wv + g.Vwgt(u) > maxVwgt)
                continue;

            real_t ew_norm = (real_t)g.Ewgt(ei) / maxEw;
            real_t dx = g.coordinates[v].x - g.coordinates[u].x;
            real_t dy = g.coordinates[v].y - g.coordinates[u].y;
            real_t dz = g.coordinates[v].z - g.coordinates[u].z;
            real_t d2 = dx * dx + dy * dy + dz * dz;
            // 平方距离归一化（替代原版 1/(1+dist/dref)）
            real_t geo_score = 1. / (1. + d2 / dref2);

            real_t score = lambda * ew_norm + (1. - lambda) * geo_score;
            if (score > bestScore)
            {
                bestScore = score;
                bestU = u;
            }
        }

        if (bestU != -1)
        {
            match[v] = nc;
            match[bestU] = nc;
            ++nc;
        }
        else
        {
            match[v] = nc++;
        }
    }
    return nc;
}

// ── 复用 KwayPartition 的 BuildCoarseGraph（通过 Coarsen.cpp）──
// 注：Coarsen.cpp 中已有 BuildCoarseGraph 为 static，
//     这里我们重新调用 CoarsenGraph 的内部逻辑。
//     为避免重复代码，我们重写一个外部可调用的版本。

#include <partition/coarse.hpp>

// 利用 Coarsen.cpp 的接口：
// 我们通过 monkey-patching cmap 来指定自定义匹配，
// 然后调用 BuildCoarseGraph 的等价实现。
// 实际上复用 Coarsen.cpp 内的 BuildCoarseGraph 逻辑：
// 将 cmap 设为我们的 match，直接构造粗图。
static Graph *BuildCoarseFromMatch(Graph &fine,
                                   const std::vector<int> &match, int nc)
{
    // 与 Coarsen.cpp 的 BuildCoarseGraph 完全相同的逻辑
    const int n = fine.nvtxs;
    fine.cmap = match;

    Graph *c = new Graph();
    c->nvtxs = nc;
    c->ncon = fine.ncon;
    c->finer = &fine;
    fine.coarser = c;

    c->vwgt.assign(nc, 0);
    c->vsize.assign(nc, 0);
    // 复制坐标（加权平均）
    c->coordinates.resize(nc);
    std::vector<real_t> wsum(nc, 0.f);
    for (int v = 0; v < n; ++v)
    {
        int cv = match[v];
        real_t w = (real_t)fine.Vwgt(v);
        c->vwgt[cv] += fine.Vwgt(v);
        c->vsize[cv] += fine.Vsize(v);
        c->coordinates[cv].x += w * fine.coordinates[v].x;
        c->coordinates[cv].y += w * fine.coordinates[v].y;
        c->coordinates[cv].z += w * fine.coordinates[v].z;
        wsum[cv] += w;
    }
    for (int cv = 0; cv < nc; ++cv)
    {
        if (wsum[cv] > 0.f)
        {
            c->coordinates[cv].x /= wsum[cv];
            c->coordinates[cv].y /= wsum[cv];
            c->coordinates[cv].z /= wsum[cv];
        }
    }

    // 构造粗图 CSR（与 Coarsen.cpp 完全相同）
    std::vector<int> xadj(nc + 1, 0);
    std::vector<std::vector<int>> c2f(nc);
    for (int v = 0; v < n; ++v)
        c2f[match[v]].push_back(v);

    std::vector<int> marker(nc, -1), adjList;
    for (int cv = 0; cv < nc; ++cv)
    {
        adjList.clear();
        for (int fv : c2f[cv])
            for (int ei = fine.xadj[fv]; ei < fine.xadj[fv + 1]; ++ei)
            {
                int cu = match[fine.adjncy[ei]];
                if (cu == cv)
                    continue;
                if (marker[cu] != cv)
                {
                    marker[cu] = cv;
                    adjList.push_back(cu);
                }
            }
        xadj[cv + 1] = (int)adjList.size();
    }
    for (int cv = 0; cv < nc; ++cv)
        xadj[cv + 1] += xadj[cv];
    c->xadj = xadj;
    c->adjncy.resize(xadj[nc]);
    c->adjwgt.resize(xadj[nc], 0);
    c->nedges = xadj[nc];

    std::fill(marker.begin(), marker.end(), -1);
    std::vector<int> pos(nc), edgeIdx(nc, -1), curNbrs;
    for (int cv = 0; cv < nc; ++cv)
        pos[cv] = xadj[cv];

    for (int cv = 0; cv < nc; ++cv)
    {
        curNbrs.clear();
        for (int fv : c2f[cv])
            for (int ei = fine.xadj[fv]; ei < fine.xadj[fv + 1]; ++ei)
            {
                int fu = fine.adjncy[ei], cu = match[fu], ew = fine.Ewgt(ei);
                if (cu == cv)
                    continue;
                if (edgeIdx[cu] == -1)
                {
                    int slot = pos[cv]++;
                    c->adjncy[slot] = cu;
                    c->adjwgt[slot] = ew;
                    edgeIdx[cu] = slot;
                    curNbrs.push_back(cu);
                }
                else
                    c->adjwgt[edgeIdx[cu]] += ew;
            }
        for (int cu : curNbrs)
            edgeIdx[cu] = -1;
    }
    c->InitTvwgt();
    return c;
}
// ================================================================
//  ProjectPartition（与 KwayPartition.cpp 相同）
// ================================================================
static void ProjectPartition(Graph &f, const Graph &c, int nparts)
{
    f.where.resize(f.nvtxs);
    for (int v = 0; v < f.nvtxs; ++v)
        f.where[v] = c.where[f.cmap[v]];
    f.pwgts.assign(nparts, 0);
    for (int v = 0; v < f.nvtxs; ++v)
        f.pwgts[f.where[v]] += f.Vwgt(v);
}

// ================================================================
//  [创新2] K-means++ 初始划分
//  多次随机试验，取 inertia 最低（几何最内聚）的结果
// ================================================================
static KmeansResult InitByKmeans(Graph &gc, int nparts,
                                 const GeoKwayOptions &opts,
                                 std::mt19937 &rng)
{
    KmeansResult best;
    best.inertia = std::numeric_limits<real_t>::max();
    std::vector<int> bestWhere;
    std::vector<int> bestPwgts;

    KmeansOpts ko;
    ko.maxIter = opts.kmeansMaxIter;
    ko.ubFactor = opts.ubFactor * 1.5; // 几何初始化允许更大不平衡，FM 后修复
    ko.verbose = false;
    ko.bfsSeedExpansion   = opts.kmeansBfsSeed;
    ko.useTopoBalanceFix  = opts.kmeansTopoBalanceFix;
    ko.enforceConnectivity = opts.kmeansEnforceConn;

    for (int t = 0; t < opts.initTrials; ++t)
    {
        ko.seed = (int)rng();
        KmeansResult r = KmeansPartition(gc, nparts, ko);
        if (r.inertia < best.inertia)
        {
            best = r;
            bestWhere = gc.where;
            bestPwgts = gc.pwgts;
        }
    }
    gc.where = bestWhere;
    gc.pwgts = bestPwgts;
    return best;
}

#include <common/macro.hpp>
// ================================================================
//  GeoKwayPartition  ──  主入口
// ================================================================

GeoKwayResult GeoKwayPartition(Graph &graph, const GeoKwayOptions &opts)
{
    YFEM_ASSERT((int)graph.coordinates.size() == graph.nvtxs,
                "GeoKwayPartition: 需要 graph.coordinates[]");

    const int K = opts.nparts;
    if (K < 1 || graph.nvtxs == 0)
        return {};
    if (graph.tvwgt.empty())
        graph.InitTvwgt();

    std::mt19937 rng(opts.seed);

    if (opts.verbose)
        std::printf("===== GeoKwayPartition n=%d nparts=%d alpha=%.2f =====\n",
                    graph.nvtxs, K, opts.alpha);

    // ── 1. [创新1] 几何感知粗化（带权重平衡约束） ────────────
    int stopAt = std::max(opts.coarseLimit * K, 40);
    Graph *cur = &graph;
    int seed = opts.seed;
    int nlevels = 0;

    // 计算粗化权重上限：每个超级顶点最多承载 (totalVwgt/K) * cap
    // 防止 HEM 合并出极端不平衡的超级顶点（METIS 论文 BHEM 思路）
    int totalVwgt = graph.tvwgt[0];
    int maxVwgtCap = (opts.coarsenWeightCap > 0)
                     ? (int)((real_t)totalVwgt / K * opts.coarsenWeightCap + 0.5f)
                     : 0;

    for (int lv = 0; lv < opts.maxLevels; ++lv)
    {
        if (cur->nvtxs <= stopAt)
            break;
        std::vector<int> match;
        // 当前图的权重相对总权重而言已减小，等比例放宽 cap
        int curCap = maxVwgtCap; // 维持全局 cap，保证最粗图也均衡
        int nc = GeoHEM_Match(*cur, match, seed++, opts.geoCoarsenLambda, curCap);
        real_t ratio = (real_t)nc / cur->nvtxs;
        if (lv > 0 && ratio >= opts.minCoarseRatio)
            break;
        cur = BuildCoarseFromMatch(*cur, match, nc);
        nlevels++;
    }
    Graph *gc = cur;

    if (opts.verbose)
        std::printf("[GeoCoarsen] %d levels  coarsest n=%d\n", nlevels, gc->nvtxs);

    // ── 2. [创新2] 几何初始划分 ──────────────────────────────
    GeoKwayResult res;
    res.nlevels = nlevels;
 
    if (opts.initMethod == InitMethod::KmeansPP ||
        opts.initMethod == InitMethod::Hybrid) {
 
        KmeansResult kr = InitByKmeans(*gc, K, opts, rng);
        res.kmeansInertia = kr.inertia;
        ComputeCkrinfo(*gc, K);
 
        if (opts.verbose)
            std::printf("[KmeansPP] inertia=%.4f  maxImb=%.2f%%  %s\n",
                        kr.inertia, kr.maxImbalance*100,
                        (kr.maxImbalance+1.f<=opts.ubFactor*1.5f)?"OK":"UNBAL");
 
        // Hybrid 模式：K-means 后再跑一轮标准 KwayFMCut 预热
        if (opts.initMethod == InitMethod::Hybrid) {
            KwayFMOpts fmo; fmo.nparts=K; fmo.nIter=3; fmo.ubFactor=opts.ubFactor;
            KwayFMCut(*gc, fmo);
        }
    } else {
        // 退化到 RB
        RBOptions rbo; rbo.nparts=K; rbo.nTrials=opts.rbTrials;
        rbo.nFMPasses=6; rbo.ubFactor=opts.ubFactor;
        rbo.postRefine=PostRefineMode::None; rbo.seed=(int)rng();
        RecursiveBisect(*gc, rbo);
        ComputeCkrinfo(*gc, K);
    }

    // ── 3. [创新3] 反粗化 + 几何感知 FM 精化 ─────────────────
    GeoFMOpts gfo;
    gfo.nparts=K; gfo.nIter=opts.nFMIter; gfo.ubFactor=opts.ubFactor;
    gfo.alpha=opts.alpha; gfo.autoBeta=opts.autoBeta; gfo.beta=opts.beta;
    gfo.verbose=opts.verbose;
 
    cur = gc;
    while (cur->finer != nullptr) {
        Graph* fine = cur->finer;
        ProjectPartition(*fine, *cur, K);
 
        if (opts.verbose)
            std::printf("[Uncoarsen] n=%d  (from coarse n=%d)\n",
                        fine->nvtxs, cur->nvtxs);
 
        ComputeCkrinfo(*fine, K);
        auto cen = ComputeCentroids(*fine, K);
 
        GeoFMResult gfr;
        if (!opts.useVolume)
            gfr = GeoKwayFMCut(*fine, gfo, cen);
        else {
            ComputeVkrinfo(*fine, K);
            gfr = GeoKwayFMVol(*fine, gfo, cen);
        }
        res.betaUsed = gfr.betaUsed;
 
        if (opts.verbose)
            std::printf("  GeoFM: topoCut=%d  nMoves=%d  beta=%.4f\n",
                        gfr.gainTopo, gfr.nMoves, gfr.betaUsed);
        cur = fine;
    }
    // 最终精化（原图层）
    {
        if (graph.ckrinfo.empty()) ComputeCkrinfo(graph, K);
        auto cen = ComputeCentroids(graph, K);
        GeoFMResult gfr;
        if (!opts.useVolume)
            gfr = GeoKwayFMCut(graph, gfo, cen);
        else {
            ComputeVkrinfo(graph, K);
            gfr = GeoKwayFMVol(graph, gfo, cen);
        }
        res.betaUsed = gfr.betaUsed;
    }
 
    // ── 4. 清理 ────────────────────────────────────────────────
    FreeCoarseGraphs(graph);

    // ── 4.5 连通性修复（FE/CFD 网格关键步骤） ─────────────────
    //   K-means 初始化 + FM 精化结束后，仍可能存在分区孤岛。
    //   对 MPI 并行求解，孤岛会显著增加通信成本。
    //   策略：交替 EnforceConnectivity / PostBalanceFix，
    //   直到两者都不再修改（或迭代上限）
    if (opts.enforceConnectivity)
    {
        bool didRepair = false;
        for (int rep = 0; rep < 3; ++rep)
        {
            ConnectRepairResult cr = EnforceConnectivity(graph, K,
                                                         opts.ubFactor * 1.10f, // 修连通时容许稍宽
                                                         0.95f);
            int moved = PostBalanceFix(graph, K, opts.ubFactor);
            if (cr.verticesMoved > 0 || moved > 0) didRepair = true;
            if (opts.verbose && (cr.isolatedCount > 0 || moved > 0))
                std::printf("[ConnFix rep%d] iso=%d connMoved=%d balMoved=%d cut=%d\n",
                            rep, cr.isolatedCount, cr.verticesMoved, moved, graph.mincut);
            if (cr.isolatedCount == 0 && moved == 0) break;
        }
        // 连通+平衡修复后做一轮额外 FM 精化恢复切边质量
        // FM 自身有 ubFactor 约束，不会破坏前面 PostBalanceFix 达成的平衡
        if (didRepair)
        {
            KwayFMOpts fmo;
            fmo.nparts = K;
            fmo.nIter = std::max(2, opts.nFMIter / 2);
            fmo.ubFactor = opts.ubFactor;
            fmo.verbose = false;
            ComputeCkrinfo(graph, K);
            int g_extra = KwayFMCut(graph, fmo);
            if (opts.verbose)
                std::printf("[PostFM] gain=%d  cut=%d\n", g_extra, graph.mincut);
        }
    }

    // ── 5. 最终统计 ────────────────────────────────────────────
    { int c=0;
      for(int v=0;v<graph.nvtxs;++v)
          for(int ei=graph.xadj[v];ei<graph.xadj[v+1];++ei)
              if(graph.where[v]!=graph.where[graph.adjncy[ei]]) c+=graph.Ewgt(ei);
      graph.mincut=c/2; }
 
    real_t ideal=(real_t)graph.tvwgt[0]/K;
    real_t maxImb=0, sumImb=0;
    for(int p=0;p<K;++p){
        real_t d=std::abs((real_t)graph.pwgts[p]/ideal-1.f);
        maxImb=std::max(maxImb,d); sumImb+=d;
    }
    res.mincut=graph.mincut; res.minvol=graph.minvol;
    res.maxImbalance=maxImb; res.avgImbalance=sumImb/K;
    res.balanced=(maxImb+1.f<=opts.ubFactor);
    res.partWeights=graph.pwgts;
 
    if (opts.verbose) {
        std::printf("\n===== GeoKwayPartition 完成 =====\n");
        std::printf("粗化层数    : %d\n", nlevels);
        std::printf("总切边数    : %d\n", res.mincut);
        std::printf("K-means惯性 : %.4f\n", res.kmeansInertia);
        std::printf("beta(自动)  : %.4f\n", res.betaUsed);
        std::printf("最大不平衡  : %.2f%%\n", res.maxImbalance*100);
        std::printf("满足约束    : %s\n", res.balanced?"是":"否");
        std::printf("理想权重    : %.1f\n", ideal);
        for (int p=0;p<K;++p){
            real_t r=(real_t)graph.pwgts[p]/ideal;
            std::printf("  [%2d] wgt=%-6d ratio=%.3f%s\n",
                        p, graph.pwgts[p], r, r>opts.ubFactor?" ←超标":"");
        }
        std::printf("\n");
    }
    return res;
}

