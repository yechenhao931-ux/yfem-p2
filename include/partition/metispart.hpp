#ifndef METISPART_HPP
#define METISPART_HPP

/**
 * metispart.h  ──  原生 METIS 划分接口封装
 *
 * 直接调用 METIS 5.x 的官方接口：
 *   · METIS_PartGraphKway       —— 多级 k-way 划分
 *   · METIS_PartGraphRecursive  —— 多级递归二分
 *
 * 作用：作为论文实验中的 **baseline**，与本项目的 GeoKwayPartition
 *       在完全相同的图（Graph CSR + 顶点/边权）上对比。
 *       因此本接口复用 GeoKwayResult 作为返回类型，统计口径与
 *       GeoKwayPartition 完全一致（mincut / minvol / maxImbalance /
 *       avgImbalance / balanced / partWeights）。
 *
 * 说明：betaUsed / kmeansInertia / nlevels 是本项目几何方法特有的指标，
 *       METIS 不暴露，保持默认值 0。
 *
 * 结果写入 graph.where[] / graph.pwgts[] / graph.mincut / graph.minvol。
 */

#include "../graph/graph.hpp"
#include "geokwaypartition.hpp" // 复用 GeoKwayResult

// METIS 的两种划分算法
enum class MetisAlgo
{
    Kway,              // METIS_PartGraphKway（多级 k-way，默认）
    RecursiveBisection // METIS_PartGraphRecursive（多级递归二分）
};

// 优化目标
enum class MetisObjective
{
    EdgeCut, // 最小化切边数（OBJTYPE_CUT；两种算法均支持）
    Volume   // 最小化通信量（OBJTYPE_VOL；仅 Kway 支持，Recursive 时回退为切边）
};

struct MetisOptions
{
    int nparts = 4;

    MetisAlgo algo = MetisAlgo::Kway;
    MetisObjective objective = MetisObjective::EdgeCut;

    // 允许的最大不平衡比，映射到 METIS ufactor = (ubFactor-1)*1000
    // 例：1.03 → ufactor=30（与 METIS k-way 默认值一致）
    real_t ubFactor = 1.03;

    int nIter = 10; // 每层精化迭代数      (METIS_OPTION_NITER)
    int nCuts = 1;  // 不同初始划分试验数  (METIS_OPTION_NCUTS)

    // 以下两项仅 Kway 有效
    bool contig = false;  // 强制每个分区连通 (METIS_OPTION_CONTIG)
    bool minconn = false; // 最小化子域连接度 (METIS_OPTION_MINCONN)

    int seed = 42;
    bool verbose = false;
};

// 是否编译进了真实 METIS（即配置时 USE_METIS 已定义）。
// 未编译时 MetisPartition 不做划分、返回空结果并打印提示。
bool MetisAvailable();

/**
 * MetisPartition
 *   用原生 METIS 对 graph 做 k-way 划分。
 *   结果写入 graph.where[]，并返回与 GeoKwayPartition 同口径的统计。
 */
GeoKwayResult MetisPartition(Graph &graph, const MetisOptions &opts = {});

#endif
