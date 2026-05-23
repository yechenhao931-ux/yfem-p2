#ifndef GEO_KWAY_REFINE_HPP
#define GEO_KWAY_REFINE_HPP
/**
 * GeoKwayRefine.h  ──  几何感知 k-way FM 精化  【论文核心创新点】
 *
 * ════════════════════════════════════════════════════════════════
 * 创新点：混合增益函数
 *
 *   标准 k-way FM 的移动增益仅考虑图拓扑：
 *     gain_topo(v → p) = Cnbr[p].ed - ckrinfo[v].id
 *
 *   本方法引入空间坐标项，构造混合增益：
 *     gain_geo(v → p)  = dist²(v, centroid[src]) - dist²(v, centroid[dst])
 *                       （移动后 v 与新分区重心更近 → 正增益）
 *
 *     gain_hybrid(v→p) = α × gain_topo + (1-α) × β × gain_geo
 *
 *   其中 β 为尺度归一化因子（将 geo_gain 映射到与 topo_gain 相当的量级）
 *   α ∈ [0,1] 是拓扑权重，由用户控制：
 *     α = 1.0 → 退化为标准 KwayFMCut
 *     α = 0.0 → 纯几何驱动（类 K-means 移动）
 *     α ≈ 0.7 → 推荐值：拓扑主导，几何辅助打破平局
 *
 * 物理直觉：
 *   · 空间相邻的顶点更可能属于同一分区（局部连接更紧密）
 *   · 在切边数相近的候选中，选择空间上更内聚的移动方向
 *   · 对不规则网格（局部加密、非结构网格）效果显著
 *
 * 实现细节：
 *   · 每次移动后增量更新受影响分区的重心
 *     （避免每轮全量重算 O(n) 的 ComputeCentroids）
 *   · 带历史回滚（同 KwayFMCut）保证质量不退化
 *   · β 自动标定：初始时对所有边界点采样，
 *     计算 topo_gain 和 geo_gain 的 RMS 比值
 * ════════════════════════════════════════════════════════════════
 */

#include "../graph/graph.hpp"
#include "kwayrefine.hpp"
#include "kmeanspart.hpp"
#include <vector>

struct GeoFMOpts
{
    int nparts = 4;
    int nIter = 10;
    real_t ubFactor = 1.03;
    real_t alpha = 0.7;     // 拓扑权重；1-alpha 为几何权重
    bool autoBeta = true;   // 自动标定 beta（topo/geo 量级比）
    real_t beta = 1.0;   // 手动 beta（autoBeta=false 时有效）
    // 【新增 软平衡感知项】
    //   FM 混合增益新增 balance 项：当 src 偏重 / dst 偏轻时，
    //   给该方向的移动一个正向奖励（与 topo 同量级，自动标定）。
    //   作用：在 ubFactor 硬约束之外提供"主动均衡"软推力，
    //         避免 FM 在硬上限边缘震荡，从而同时改善 mincut 与
    //         maxImbalance。0 = 关闭（旧行为）；推荐 0.3~0.6。
    real_t balanceGamma = 0.4;
    bool verbose = false;
};

struct GeoFMResult {
    int    gainTopo    = 0;  // 本次精化的纯拓扑增益（切边减少数）
    int    gainGeo     = 0;  // 几何内聚增益（加权距离减少量，仅用于统计）
    real_t betaUsed    = 0;  // 实际使用的 beta 值
    int    nMoves      = 0;  // 有效移动次数
};

/**
 * GeoKwayFMCut
 *   几何感知 FM 精化（最小化切边 + 几何内聚）
 *   需要 graph.ckrinfo[] 已由 ComputeCkrinfo 填充
 *   需要 graph.coordinates[] 已填充
 *   centroids: 当前各分区加权重心（由 ComputeCentroids 获得），
 *              函数内部会随移动增量更新
 */

GeoFMResult GeoKwayFMCut(Graph& graph, const GeoFMOpts& opts,
                          std::vector<std::vector<real_t>>& centroids);


/**
 * GeoKwayFMVol
 *   几何感知 FM 精化（最小化通信量 + 几何内聚）
 */

GeoFMResult GeoKwayFMVol(Graph& graph, const GeoFMOpts& opts,
                          std::vector<std::vector<real_t>>& centroids);

// 优化的GeoKwayRefine
/**
 * GeoKwayFMCut
 *   几何感知 FM 精化（最小化切边 + 几何内聚）
 *   需要 graph.ckrinfo[] 已由 ComputeCkrinfo 填充
 *   需要 graph.coordinates[] 已填充
 *   centroids: 当前各分区加权重心（由 ComputeCentroids 获得），
 *              函数内部会随移动增量更新
 */
GeoFMResult GeoKwayFMCut_Opt1(Graph& graph, const GeoFMOpts& opts,
                          std::vector<std::vector<real_t>>& centroids);

/**
 * GeoKwayFMVol
 *   几何感知 FM 精化（最小化通信量 + 几何内聚）
 */
GeoFMResult GeoKwayVol_Opt1(Graph& graph, const GeoFMOpts& opts,
                          std::vector<std::vector<real_t>>& centroids);

#endif