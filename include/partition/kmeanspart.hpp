#ifndef KMEANSPART_HPP
#define KMEANSPART_HPP
/**
 * KmeansPart.h  ──  基于坐标的 K-means++ 几何划分
 *
 * 算法：
 *   1. K-means++ 初始化：以概率正比于 dist²(v, 已选中心) 选种子
 *   2. Lloyd 迭代：分配 + 重心更新，直到收敛或达到最大迭代数
 *   3. 空间平衡修复：若某簇顶点权重严重超标，
 *      将边界顶点迁移到最近的轻簇
 *
 * 作为创新点的作用：
 *   · 替代纯随机的 RB 初始划分，提供"几何合理"的起始分区
 *   · 对坐标质量好的网格（有限元、CFD），可大幅减少后续 FM 迭代轮数
 *   · 可与拓扑 FM 混合，构成 GeoKway 的初始化子模块
 */

#include "../graph/graph.hpp"
#include <vector>

struct KmeansOpts
{
    int maxIter = 100;  // Lloyd 最大迭代数
    real_t tol = 1e-4;  // 重心移动收敛阈值
    int seed = 42;
    real_t ubFactor  = 1.10;  // 空间平衡修复的不平衡上限
    bool  verbose   = false;
    // ── 拓扑感知扩展 ─────────────────────────────────────────
    bool useTopoBalanceFix = true; // 用 BalanceFix(连通性优先) 代替 BalanceRepair
    bool enforceConnectivity = true; // 修复 K-means 结果中可能出现的孤岛
    bool bfsSeedExpansion = true; // 用 BFS 邻接扩展 K-means++ 种子，提升连通性
};

struct KmeansResult {
    real_t inertia   = 0.;    // 总簇内方差（越小越好）
    int   iters     = 0;      // 实际迭代次数
    bool  converged = false;
    real_t maxImbalance = 0.;
    std::vector<std::vector<real_t>> centroids; // K×3，最终重心坐标
};

/**
 * KmeansPartition
 *   对 graph 按坐标做 K-means 划分，结果写入 graph.where[]
 *   要求 graph.coordinates 已填充（大小 == nvtxs）
 *   nparts: 目标分区数
 */
KmeansResult KmeansPartition(Graph& graph, int nparts, const KmeansOpts& opts = {});
 
/**
 * ComputeCentroids
 *   根据当前 graph.where[] 和 graph.coordinates 计算每个分区的加权重心
 *   返回 nparts × 3 的向量
 */
std::vector<std::vector<real_t>> ComputeCentroids(const Graph& graph, int nparts);

                                                                                                                                                                 
// OPT1
KmeansResult KmeansPartition_Opt1(Graph& graph, int nparts, const KmeansOpts& opts = {});

#endif