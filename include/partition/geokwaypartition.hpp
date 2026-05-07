#ifndef GEOKWAYPARTITION_HPP
#define GEOKWAYPARTITION_HPP
/**
 * GeoKwayPartition.h  ──  几何引导多级 k-way 划分（GGKP）
 *
 * ════════════════════════════════════════════════════════════════
* 论文创新点：
 *
 *   传统 KwayPartition 流程（Metis）：
 *     粗化(HEM) → 初始划分(RB) → 反粗化(KwayFM)
 *
 *   本方法 GeoKwayPartition 三层创新：
 *
 *   [创新 1] 几何感知粗化 (Geometry-Aware Coarsening, GAC)
 *     HEM 匹配分数 = λ × edge_weight + (1-λ) × 1/(1 + dist(u,v))
 *     使粗化时优先合并空间相邻的顶点对，提高粗图的几何质量
 *
 *   [创新 2] K-means++ 几何初始划分 (Geometric Initial Partitioning)
 *     在最粗图上用 K-means++ 而非 RB 做初始划分
 *     利用坐标信息得到空间连续的分区，减少后续精化的迭代次数
 *     对多物理场网格、非结构有限元网格效果显著
 *
 *   [创新 3] 几何感知 FM 精化 (Geo-FM, GeoKwayRefine.h)
 *     混合增益 = α × 拓扑增益 + (1-α) × β × 几何增益
 *     几何增益 = dist²(v,centroid[src]) - dist²(v,centroid[dst])
 *     带 autoBeta 自动量级标定，带历史回滚保证质量
 *
 * 对比实验建议（论文第5章）：
 *   A. Metis-style KwayPartition  (baseline)
 *   B. GGKP with alpha=1.0        (仅 K-means 初始化，FM 不变)
 *   C. GGKP with alpha=0.7        (完整混合增益，推荐)
 *   D. GGKP with alpha=0.3        (几何主导，验证极端情况)
 *
 * 评价指标：
 *   · mincut（切边数）
 *   · maxImbalance（最大不平衡率）
 *   · 运行时间（粗化层数 × 每层 FM 迭代数）
 *   · 分区连通性（每分区是否连通）
 * ════════════════════════════════════════════════════════════════
 */
#include "../graph/graph.hpp"
#include "coarse.hpp"
#include "kwayrefine.hpp"
#include "kmeanspart.hpp"
#include "geokwayrefine.hpp"
#include "recursivebisect.hpp"


enum class InitMethod{
    KmeansPP, // 【创新】K-means++ 几何初始划分（默认）
    RecursiveBisect, // 传统 RB 初始划分（对照组）
    Hybrid          // K-means++ 结果 + 一轮标准 KwayFMCut 预热
};

struct GeoKwayOptions
{
    int nparts = 4;
    // ── 粗化参数 ──────────────────────────────────────────────
    int coarseLimit = 20;
    real_t minCoarseRatio = 0.75;
    int maxLevels = 128;
    real_t geoCoarsenLambda = 0.5f; // [创新1] 粗化中边权/空间距离的混合比
                                    // 0=纯空间 1=纯边权（退化为标准HEM）

    // ── 初始划分参数 ──────────────────────────────────────────
    InitMethod initMethod  = InitMethod::KmeansPP; // [创新2]
    int   kmeansMaxIter    = 100;
    int   initTrials       = 3;    // K-means 多次随机试验取最优
    int   rbTrials         = 6;    // 退化为 RB 时的试验数
 
    // ── 精化参数 ──────────────────────────────────────────────
    real_t alpha = 0.7;   // [创新3] 拓扑权重（0~1）
    bool autoBeta = true;
    real_t beta = 1.0;
    // 【新增】软平衡感知项（GeoFM 内）
    //   gain_total = α*topo + (1-α)*β*geo + γ*balance
    //   γ=0 退化为旧行为；γ>0 提供"重→轻"软推力，常同时改善 mincut 与
    //   maxImbalance。推荐 0.3~0.6。
    real_t balanceGamma = 0.4;
    int nFMIter = 10;
    real_t ubFactor = 1.03;

    // ── 目标 ──────────────────────────────────────────────────
    bool  useVolume        = false; // true → 优化通信量（GeoKwayFMVol）

    // ── 连通性 / 平衡感知 ─────────────────────────────────────
    bool  enforceConnectivity = true;  // 修复后置阶段产生的孤岛
    real_t coarsenWeightCap   = 1.5f;  // 粗化权重上限 = (W/K) * coarsenWeightCap
                                       // 设为 0 关闭，> 0 时 HEM 不合并使
                                       // vwgt 超过 (totalVwgt / nparts) * cap 的顶点对
    // K-means 内部行为（透传到 KmeansOpts，便于消融实验）
    bool  kmeansBfsSeed       = true;  // BFS 邻接扩展 K-means++ 种子
    bool  kmeansTopoBalanceFix = true; // 用拓扑感知 BalanceFix 取代纯几何 BalanceRepair
    bool  kmeansEnforceConn   = true;  // K-means 内部的孤岛修复

    int   seed             = 42;
    bool  verbose          = false;
};


struct GeoKwayResult {
    int   mincut           = 0;
    int   minvol           = 0;
    real_t maxImbalance     = 0.;
    real_t avgImbalance     = 0.;
    bool  balanced         = false;
    int   nlevels          = 0;
    real_t betaUsed         = 0.;   // 实际使用的 beta 值（auto 时）
    real_t kmeansInertia    = 0.;   // K-means 初始化的簇内方差
    std::vector<int> partWeights;
};

/**
 * GeoKwayPartition
 *   主入口：几何引导多级 k-way 划分
 *   要求 graph.coordinates[] 已填充
 */
GeoKwayResult GeoKwayPartition(Graph& graph, const GeoKwayOptions& opts = {});

#endif