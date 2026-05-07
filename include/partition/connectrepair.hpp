#ifndef CONNECTREPAIR_HPP
#define CONNECTREPAIR_HPP
/**
 * ConnectRepair.h  ──  分区连通性修复
 *
 * 背景：
 *   K-means 等几何方法只看坐标，忽略图拓扑，
 *   常常生成包含若干"孤岛"的分区（顶点在该分区内但
 *   在图上与该分区其余顶点不连通）。
 *   对有限元并行求解，孤岛会显著增加 MPI 通信成本，
 *   并破坏代数多重网格 (BoomerAMG) 的层次构造。
 *
 * 算法：
 *   1. 对每个分区做 BFS 找连通子图
 *   2. 保留每个分区中权重最大的连通子图作为"主体"
 *   3. 将其他"孤岛"按以下优先级重新分配：
 *        a. 选与孤岛切边最多的邻居分区
 *        b. 若候选分区即将超重，则选次重的
 *        c. 若所有邻居都超重，回到原分区（保持现状）
 *   4. 重复直到所有分区都连通（或迭代上限）
 *
 * 复杂度：O(n + m) 每轮，最多 nparts 轮。
 */
#include "../graph/graph.hpp"
#include <vector>

struct ConnectRepairResult
{
    int isolatedCount = 0;     // 检测到的孤岛个数
    int verticesMoved = 0;     // 实际重分配的顶点数
    int finalCut = 0;          // 修复后切边数
    bool fullyConnected = false; // 是否所有分区都连通
};

/**
 * EnforceConnectivity
 *   修复 graph.where[] 中的连通性问题。
 *   - nparts: 分区数
 *   - ubFactor: 平衡上限（不允许超过 ubFactor*ideal 才迁移）
 *   - keepLargestRatio: 主体连通子图最少占该分区权重的比例
 *                       低于此比例的连通子图视为孤岛 (默认 0.95)
 *   要求 graph.where[] 已填充。
 *   函数会更新 graph.where[]、graph.pwgts[]、graph.mincut。
 */
ConnectRepairResult EnforceConnectivity(Graph& graph,
                                        int nparts,
                                        real_t ubFactor = 1.05,
                                        real_t keepLargestRatio = 0.95);

/**
 * PostBalanceFix
 *   双向平衡修复（重→轻 / 轻←重），仅在边界点上做迁移，
 *   保留分区连通性。设计目标：在 EnforceConnectivity 引发的
 *   小幅不平衡后做兜底修复，使最终结果满足 ubFactor 约束。
 *
 *   返回成功迁移的顶点数。
 */
int PostBalanceFix(Graph& graph, int nparts, real_t ubFactor);

#endif
