#ifndef COARSE_HPP
#define COARSE_HPP

/**
 * Coarsen.h  ──  多级粗化接口
 *
 * 算法：Heavy Edge Matching (HEM) + 随机打乱访问顺序
 *
 *   1. 随机打乱顶点访问顺序
 *   2. 对每个未匹配顶点 v，选其未匹配邻居中边权最重的 u 配对
 *      若无未匹配邻居则自匹配（孤立点）
 *   3. 每对匹配 (v,u) 收缩为一个超级顶点：
 *        vwgt_new  = vwgt[v] + vwgt[u]
 *        vsize_new = vsize[v] + vsize[u]
 *   4. 构造粗化图：合并平行边（累加权重），去掉自环
 *
 * 粗化到以下任一条件停止：
 *   - nvtxs <= coarseLimit（顶点数足够少）
 *   - 粗化比 < minCoarseRatio（粗化效果不佳，不再继续）
 *   - 达到最大层数 maxLevels
 *
 * 粗化后每层图通过 Graph::coarser / finer 双向链表连接。
 * 调用方负责用 FreeCoarseGraphs() 释放内存。
 */


#include "../graph/graph.hpp"

struct CoarsenOpts
{
    int coarseLimit = 20; // 粗化停止的顶点数下限（相对于 nparts 的倍数）
    real_t minCoarseRation = 0.75; //粗化比 >= 此值则停止（收益不足）

    int maxLevels = 128;   // 最大粗化层数
    int seed = 42;

    // ── A-2: union-find 粗化（G-kway）──
    //   true  : 用带打分的 union-find 粗化（一趟选邻居可并行，子集可合并多点，层数更少）
    //   false : 经典 Heavy-Edge-Matching（默认，保持原行为）
    bool useUnionFind = false;
    //   超级顶点权重上限系数：>0 时 maxVwgt = coarsenWeightCap × totalVwgt/nparts，
    //   防止 union-find 合出过重子集；0 表示不限制（遵循论文仅靠打分）
    real_t coarsenWeightCap = 0;
};

// 从 fine 出发递归向上粗化，返回最粗一层的指针（已 new）
// fine.coarser 链会被填充；fine.cmap 会被设置
// nparts 用于计算 coarseLimit 的实际阈值
Graph* CoarsenGraph(Graph& fine, int nparts, const CoarsenOpts& opts = {});

// 释放从 g.coarser 开始的所有粗化层
void FreeCoarseGraphs(Graph& g);





#endif