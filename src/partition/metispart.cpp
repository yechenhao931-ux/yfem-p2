// ════════════════════════════════════════════════════════════════
//  metispart.cpp  ──  原生 METIS 划分接口封装
//
//  ⚠ 头文件包含顺序很关键：
//    项目 type.hpp 里有  `#define real_t double`（宏），
//    而 metis.h 里是      `typedef float real_t;`。
//    若先包含项目头，real_t 变成宏，metis.h 的 typedef 会被展开成
//    `typedef float double;` 直接编译失败。
//    因此先包含 <metis.h>，让它的 typedef 正常解析（METIS 的函数
//    原型此时即按 idx_t / float 完成解析）；随后项目宏只影响后续
//    文本，不再波及已解析的 METIS 声明。
// ════════════════════════════════════════════════════════════════
#ifdef USE_METIS
#include <metis.h>
namespace
{
// 在项目宏 `#define real_t double` 生效前，先把 METIS 的整型别名固定下来。
// （idx_t 不与项目冲突，这里取别名仅为可读性与可移植性。）
using metis_idx_t = ::idx_t;
} // namespace
#endif

#include <partition/metispart.hpp>
#include <partition/kwayrefine.hpp> // ComputeCkrinfo / ComputeVkrinfo

#include <cmath>
#include <cstdio>
#include <vector>

bool MetisAvailable()
{
#ifdef USE_METIS
    return true;
#else
    return false;
#endif
}

// 统计填充：与 GeoKwayPartition 第 5 节完全相同的口径，保证可比性。
static void FillStats(Graph &g, int K, real_t ubFactor, GeoKwayResult &res)
{
    if (g.tvwgt.empty())
        g.InitTvwgt();

    ComputeCkrinfo(g, K); // 填充 pwgts / mincut / ckrinfo
    ComputeVkrinfo(g, K); // 填充 minvol（依赖 ckrinfo）

    real_t ideal = (real_t)g.tvwgt[0] / K;
    real_t maxImb = 0, sumImb = 0;
    for (int p = 0; p < K; ++p)
    {
        real_t d = std::abs((real_t)g.pwgts[p] / ideal - 1.f);
        maxImb = std::max(maxImb, d);
        sumImb += d;
    }

    res.mincut = g.mincut;
    res.minvol = g.minvol;
    res.maxImbalance = maxImb;
    res.avgImbalance = sumImb / K;
    res.balanced = (maxImb + 1.f <= ubFactor);
    res.partWeights = g.pwgts;
    // nlevels / betaUsed / kmeansInertia 为本项目几何方法特有，METIS 不暴露，保持 0。
}

GeoKwayResult MetisPartition(Graph &graph, const MetisOptions &opts)
{
    const int K = opts.nparts;
    GeoKwayResult res;

    if (K < 1 || graph.nvtxs == 0)
        return res;
    if (graph.tvwgt.empty())
        graph.InitTvwgt();

    // 平凡情形：单分区，全部归 0。
    if (K == 1)
    {
        graph.where.assign(graph.nvtxs, 0);
        FillStats(graph, K, opts.ubFactor, res);
        return res;
    }

#ifndef USE_METIS
    (void)opts;
    std::fprintf(stderr,
                 "[MetisPartition] 未编译 METIS 支持（USE_METIS 未定义）。"
                 "请用 cmake -DENABLE_METIS=ON 重新配置后再编译。\n");
    return res; // 空结果
#else
    const int n = graph.nvtxs;

    // ── 1. Graph CSR → METIS idx_t 数组 ───────────────────────────
    metis_idx_t nvtxs = n;
    metis_idx_t ncon = graph.ncon > 0 ? graph.ncon : 1;

    std::vector<metis_idx_t> xadj(graph.xadj.begin(), graph.xadj.end());
    std::vector<metis_idx_t> adjncy(graph.adjncy.begin(), graph.adjncy.end());

    // 顶点权重：vwgt 为空 / 长度不匹配 → 传 NULL（METIS 视为全 1，
    // 与 Graph::Vwgt() 在空时返回 1 一致）。
    std::vector<metis_idx_t> vwgt;
    metis_idx_t *vwgtPtr = nullptr;
    if ((int)graph.vwgt.size() == n * (int)ncon)
    {
        vwgt.assign(graph.vwgt.begin(), graph.vwgt.end());
        vwgtPtr = vwgt.data();
    }

    // 边权：adjwgt 为空 → NULL（全 1，与 Graph::Ewgt() 一致）。
    std::vector<metis_idx_t> adjwgt;
    metis_idx_t *adjwgtPtr = nullptr;
    if (!graph.adjwgt.empty() && graph.adjwgt.size() == graph.adjncy.size())
    {
        adjwgt.assign(graph.adjwgt.begin(), graph.adjwgt.end());
        adjwgtPtr = adjwgt.data();
    }

    // 通信量目标需要 vsize；为空 → NULL（全 1，与 Graph::Vsize() 一致）。
    std::vector<metis_idx_t> vsize;
    metis_idx_t *vsizePtr = nullptr;
    if (opts.objective == MetisObjective::Volume && (int)graph.vsize.size() == n)
    {
        vsize.assign(graph.vsize.begin(), graph.vsize.end());
        vsizePtr = vsize.data();
    }

    metis_idx_t nparts = K;
    metis_idx_t objval = 0;
    std::vector<metis_idx_t> part(n, 0);

    // ── 2. METIS options ──────────────────────────────────────────
    metis_idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0; // 0-based 编号
    options[METIS_OPTION_SEED] = opts.seed;
    if (opts.nIter > 0)
        options[METIS_OPTION_NITER] = opts.nIter;
    if (opts.nCuts > 0)
        options[METIS_OPTION_NCUTS] = opts.nCuts;

    // ufactor：允许的最大不平衡 = (1 + ufactor/1000)
    metis_idx_t uf = (metis_idx_t)std::lround((opts.ubFactor - 1.0) * 1000.0);
    if (uf < 1)
        uf = 1;
    options[METIS_OPTION_UFACTOR] = uf;
    options[METIS_OPTION_DBGLVL] = 0;

    const bool useKway = (opts.algo == MetisAlgo::Kway);

    // OBJTYPE / CONTIG / MINCONN 仅 Kway 有意义；Recursive 仅做切边。
    if (useKway)
    {
        options[METIS_OPTION_OBJTYPE] =
            (opts.objective == MetisObjective::Volume) ? METIS_OBJTYPE_VOL
                                                       : METIS_OBJTYPE_CUT;
        options[METIS_OPTION_CONTIG] = opts.contig ? 1 : 0;
        options[METIS_OPTION_MINCONN] = opts.minconn ? 1 : 0;
    }

    // ── 3. 调用 METIS ─────────────────────────────────────────────
    int rc;
    if (useKway)
    {
        rc = METIS_PartGraphKway(&nvtxs, &ncon, xadj.data(), adjncy.data(),
                                 vwgtPtr, vsizePtr, adjwgtPtr, &nparts,
                                 /*tpwgts*/ nullptr, /*ubvec*/ nullptr,
                                 options, &objval, part.data());
    }
    else
    {
        rc = METIS_PartGraphRecursive(&nvtxs, &ncon, xadj.data(), adjncy.data(),
                                      vwgtPtr, vsizePtr, adjwgtPtr, &nparts,
                                      /*tpwgts*/ nullptr, /*ubvec*/ nullptr,
                                      options, &objval, part.data());
    }

    if (rc != METIS_OK)
    {
        std::fprintf(stderr, "[MetisPartition] METIS 返回错误码 %d\n", rc);
        return res; // 空结果
    }

    // ── 4. 写回 where[] 并统计 ────────────────────────────────────
    graph.where.assign(n, 0);
    for (int v = 0; v < n; ++v)
        graph.where[v] = (int)part[v];

    FillStats(graph, K, opts.ubFactor, res);

    if (opts.verbose)
    {
        const char *algoName = useKway ? "Kway" : "RecursiveBisection";
        const char *objName =
            (useKway && opts.objective == MetisObjective::Volume) ? "vol" : "cut";
        std::printf("\n===== MetisPartition (%s, obj=%s) =====\n", algoName, objName);
        std::printf("顶点数      : %d\n", n);
        std::printf("分区数      : %d\n", K);
        std::printf("METIS objval: %d\n", (int)objval);
        std::printf("总切边数    : %d\n", res.mincut);
        std::printf("通信量      : %d\n", res.minvol);
        std::printf("最大不平衡  : %.2f%%\n", res.maxImbalance * 100);
        std::printf("满足约束    : %s\n", res.balanced ? "是" : "否");
        real_t ideal = (real_t)graph.tvwgt[0] / K;
        for (int p = 0; p < K; ++p)
        {
            real_t r = (real_t)graph.pwgts[p] / ideal;
            std::printf("  [%2d] wgt=%-6d ratio=%.3f%s\n", p, graph.pwgts[p], r,
                        r > opts.ubFactor ? " ←超标" : "");
        }
        std::printf("\n");
    }
    return res;
#endif
}
