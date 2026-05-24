#ifndef KWAYPARTITION_HPP
#define KWAYPARTITION_HPP

/**
 * KwayPartition.h  ──  多级 k-way 划分主接口
 *
 * 完整流程：
 *   ┌─ 粗化阶段  ──────────────────────────────────────────┐
 *   │  G0 (原图) → G1 → G2 → ... → Gc (最粗图)           │
 *   │  HEM 匹配 + BuildCoarseGraph                         │
 *   └──────────────────────────────────────────────────────┘
 *   ┌─ 初始划分  ──────────────────────────────────────────┐
 *   │  对 Gc 用递归二分（RB）得到初始 k-way 分区           │
 *   │  可重复多次随机种子取最优                            │
 *   └──────────────────────────────────────────────────────┘
 *   ┌─ 反粗化 + 精化  ────────────────────────────────────┐
 *   │  for lv = c..0:                                      │
 *   │    ProjectPartition: where[fine] ← cmap + where[coarse]│
 *   │    KwayFMCut / KwayFMVol: 在 fine 层精化            │
 *   └──────────────────────────────────────────────────────┘
 *
 * 结果写入 graph.where[] / graph.pwgts[] /
 *          graph.ckrinfo[] / graph.cnbrPool[]
 */

#include "../graph/graph.hpp"
#include "coarse.hpp"
#include "kwayrefine.hpp"
#include "recursivebisect.hpp"
enum class KwayObjective {
    EdgeCut,  // 最小化切边数（默认）
    Volume,   // 最小化通信量
};

struct KwayOptions
{
    int    nparts        = 4;

    // ── 粗化参数 ──────────────────────────────────────────────
    int    coarseLimit   = 20;    // 最粗图每分区顶点数下限（实际阈值=coarseLimit×nparts）
    real_t minCoarseRation= 0.75; // 粗化率低于此值时停止粗化
    int    maxLevels     = 128;
    bool   useUnionFind  = false; // A-2: union-find 粗化（替代 HEM）
    real_t coarsenWeightCap = 0;  // union-find 超级顶点权重上限系数（0=不限）

    // ── 初始划分参数 ──────────────────────────────────────────
    int   initTrials    = 5;     // 初始划分随机试验次数（每次 RB 内部也有多次试验）
    int   rbTrials      = 8;     // 每次 RB 的随机试验次数
    int   rbFMPasses    = 6;     // 每次 RB 的 FM 轮数

    // ── 精化参数 ──────────────────────────────────────────────
    KwayObjective objective = KwayObjective::EdgeCut;
    int   nFMIter       = 10;    // 每层 k-way FM 最大迭代轮数
    real_t ubFactor      = 1.03; // 允许最大不平衡比
    bool  useIndepSetRefine = false; // A-3: 独立集并行精化（替代串行 FM，仅 EdgeCut）
 
    int   seed          = 42;
    bool  verbose       = false;
};


struct KwayResult
{
    int   mincut        = 0;
    int   minvol        = 0;
    real_t maxImbalance  = 0.0;
    real_t avgImbalance  = 0.0;
    bool  balanced      = false;
    int   nlevels       = 0;     // 实际粗化层数
    std::vector<int>  partWeights;
};

// 主入口：对 graph 进行多级 k-way 划分，结果写入 graph.where[]
KwayResult KwayPartition(Graph& graph, const KwayOptions& opts = {});

#endif