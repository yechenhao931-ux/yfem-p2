/**
 * ConnectRepair.cpp  ──  分区连通性修复实现
 *
 * 关键步骤：
 *   1. 用 BFS 在每个分区内部找所有连通子图
 *   2. 每个分区保留"主体"（权重最大的连通子图）
 *   3. 对其他"孤岛"，扫描其所有跨分区边，
 *      选切边权重最大且接收后不严重超重的邻居分区
 *   4. 整体迁移孤岛（保持其内部连通性）
 *   5. 至多 nparts 轮，直到不再有孤岛或无法继续修复
 */
#include <partition/connectrepair.hpp>
#include <algorithm>
#include <queue>
#include <vector>

// ── 单个分区内 BFS 找连通子图 ─────────────────────────────────
//   返回每个连通子图的 (顶点列表, 总权重)
//   只在 part == targetPart 的顶点上做 BFS
static void FindComponentsInPart(const Graph& g, int targetPart,
                                 std::vector<bool>& visited,
                                 std::vector<std::vector<int>>& components,
                                 std::vector<int>& compWeight)
{
    const int n = g.nvtxs;
    for (int v = 0; v < n; ++v)
    {
        if (visited[v] || g.where[v] != targetPart) continue;

        std::vector<int> comp;
        int w = 0;
        std::queue<int> q;
        q.push(v);
        visited[v] = true;
        while (!q.empty())
        {
            int u = q.front(); q.pop();
            comp.push_back(u);
            w += g.Vwgt(u);
            for (int ei = g.xadj[u]; ei < g.xadj[u + 1]; ++ei)
            {
                int nb = g.adjncy[ei];
                if (!visited[nb] && g.where[nb] == targetPart)
                {
                    visited[nb] = true;
                    q.push(nb);
                }
            }
        }
        components.push_back(std::move(comp));
        compWeight.push_back(w);
    }
}

// ── 计算孤岛与各邻居分区的边权和 ──────────────────────────────
//   nbrEdgeWeight[p] = sum of edge weights from island vertices into part p
//   不包括到原分区的边
static void ComputeIslandNeighborWeight(const Graph& g,
                                        const std::vector<int>& island,
                                        int srcPart,
                                        std::vector<int>& nbrEdgeWeight)
{
    std::fill(nbrEdgeWeight.begin(), nbrEdgeWeight.end(), 0);
    for (int v : island)
    {
        for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
        {
            int u = g.adjncy[ei];
            int pu = g.where[u];
            if (pu == srcPart) continue;
            nbrEdgeWeight[pu] += g.Ewgt(ei);
        }
    }
}

ConnectRepairResult EnforceConnectivity(Graph& g, int nparts,
                                        real_t ubFactor,
                                        real_t keepLargestRatio)
{
    ConnectRepairResult res;
    if (nparts <= 1 || g.nvtxs == 0) {
        res.fullyConnected = true;
        return res;
    }

    // 重建 pwgts（防止外部状态不一致）
    g.pwgts.assign(nparts, 0);
    for (int v = 0; v < g.nvtxs; ++v)
        g.pwgts[g.where[v]] += g.Vwgt(v);

    int totalW = 0;
    for (int p = 0; p < nparts; ++p) totalW += g.pwgts[p];
    real_t ideal = (real_t)totalW / nparts;
    int maxAllow = (int)(ubFactor * ideal + 0.5);

    // 至多 nparts 轮（每轮可修复多个分区）
    for (int round = 0; round < nparts; ++round)
    {
        std::vector<bool> visited(g.nvtxs, false);
        bool anyMoved = false;
        int roundIsolated = 0;

        for (int p = 0; p < nparts; ++p)
        {
            std::vector<std::vector<int>> components;
            std::vector<int> compWeight;
            // 复位 visited 至 part p（其他分区仍标记，避免重访）
            // 简化：用一个独立 visited 数组每轮
            std::fill(visited.begin(), visited.end(), false);
            FindComponentsInPart(g, p, visited, components, compWeight);

            if (components.size() <= 1) continue;

            // 找权重最大的主体
            int mainIdx = 0;
            for (int i = 1; i < (int)components.size(); ++i)
                if (compWeight[i] > compWeight[mainIdx]) mainIdx = i;

            int mainW = compWeight[mainIdx];

            // 主体已占绝大多数 → 把孤岛迁走
            std::vector<int> nbrW(nparts, 0);
            for (int i = 0; i < (int)components.size(); ++i)
            {
                if (i == mainIdx) continue;
                if ((real_t)compWeight[i] / mainW > (1.0f - keepLargestRatio))
                {
                    // 该子图本身规模也较大，不算孤岛
                    if ((real_t)mainW >= keepLargestRatio * g.pwgts[p])
                        continue; // 主体已足够大，跳过判定
                }

                ++roundIsolated;
                ComputeIslandNeighborWeight(g, components[i], p, nbrW);

                // 选切边最多且不导致严重超重的邻居分区
                int bestP = -1, bestW = -1;
                int islW = compWeight[i];
                for (int q = 0; q < nparts; ++q)
                {
                    if (q == p || nbrW[q] == 0) continue;
                    if (g.pwgts[q] + islW > maxAllow) continue;
                    if (nbrW[q] > bestW)
                    {
                        bestW = nbrW[q];
                        bestP = q;
                    }
                }
                // 若所有切边邻居都超重，放宽限制选切边最多的（损失平衡换连通）
                if (bestP == -1)
                {
                    for (int q = 0; q < nparts; ++q)
                    {
                        if (q == p || nbrW[q] == 0) continue;
                        if (nbrW[q] > bestW)
                        {
                            bestW = nbrW[q];
                            bestP = q;
                        }
                    }
                }
                if (bestP == -1) continue; // 极少情况：孤岛与外部全无边

                // 整体迁移
                for (int v : components[i])
                {
                    g.pwgts[p] -= g.Vwgt(v);
                    g.pwgts[bestP] += g.Vwgt(v);
                    g.where[v] = bestP;
                    ++res.verticesMoved;
                }
                anyMoved = true;
            }
        }
        res.isolatedCount += roundIsolated;
        if (!anyMoved) { res.fullyConnected = true; break; }
    }

    // 重算 mincut
    int c = 0;
    for (int v = 0; v < g.nvtxs; ++v)
        for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
            if (g.where[v] != g.where[g.adjncy[ei]])
                c += g.Ewgt(ei);
    g.mincut = c / 2;
    res.finalCut = g.mincut;
    return res;
}

// ────────────────────────────────────────────────────────────────
//  PostBalanceFix  ──  双向边界点平衡修复
//  逻辑：
//    · 找超重 (>maxPW) 与超轻 (<minPW) 分区
//    · 重→轻：从超重分区找最优可移边界点（cut gain + 目标分区轻奖励）
//    · 轻→轻：若最重邻居有边界点能移过来则拉取
//    · 仅在边界顶点上移动，保留分区连通性
//  返回：迁移顶点数
// ────────────────────────────────────────────────────────────────
int PostBalanceFix(Graph& g, int K, real_t ubFactor)
{
    if (K <= 1 || g.nvtxs == 0) return 0;
    if (g.tvwgt.empty()) g.InitTvwgt();

    // 同步 pwgts（防止外部状态过期）
    g.pwgts.assign(K, 0);
    for (int v = 0; v < g.nvtxs; ++v) g.pwgts[g.where[v]] += g.Vwgt(v);

    real_t ideal = (real_t)g.tvwgt[0] / K;
    int maxPW = (int)(ubFactor * ideal + 0.5);
    int minPW = (int)((2.f - ubFactor) * ideal - 0.5);

    // 同分区边权和（id）
    std::vector<int> vId(g.nvtxs, 0);
    for (int v = 0; v < g.nvtxs; ++v)
    {
        int pv = g.where[v];
        for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
            if (g.where[g.adjncy[ei]] == pv) vId[v] += g.Ewgt(ei);
    }

    int moved = 0;
    bool changed = true;
    int guard = g.nvtxs;
    while (changed && --guard > 0)
    {
        changed = false;

        int heavy = -1, light = -1;
        for (int p = 0; p < K; ++p)
        {
            if (g.pwgts[p] > maxPW && (heavy == -1 || g.pwgts[p] > g.pwgts[heavy])) heavy = p;
            if (g.pwgts[p] < minPW && (light == -1 || g.pwgts[p] < g.pwgts[light])) light = p;
        }
        if (heavy == -1 && light == -1) break;

        // 重→邻居（重分区找最佳可移边界点）
        if (heavy != -1)
        {
            int bestV = -1, bestDst = -1;
            real_t bestScore = -1e30;
            for (int v = 0; v < g.nvtxs; ++v)
            {
                if (g.where[v] != heavy) continue;
                for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
                {
                    int u = g.adjncy[ei], pu = g.where[u];
                    if (pu == heavy) continue;
                    if (g.pwgts[pu] + g.Vwgt(v) > maxPW) continue;
                    int edToDst = 0;
                    for (int ei2 = g.xadj[v]; ei2 < g.xadj[v + 1]; ++ei2)
                        if (g.where[g.adjncy[ei2]] == pu) edToDst += g.Ewgt(ei2);
                    real_t gain = (real_t)(edToDst - vId[v]);
                    real_t score = gain + (ideal - (real_t)g.pwgts[pu]) / (ideal + 1.f) * 0.5f;
                    if (score > bestScore) { bestScore = score; bestV = v; bestDst = pu; }
                }
            }
            if (bestV != -1)
            {
                int wv = g.Vwgt(bestV);
                for (int ei = g.xadj[bestV]; ei < g.xadj[bestV + 1]; ++ei)
                {
                    int nb = g.adjncy[ei], ew = g.Ewgt(ei), pnb = g.where[nb];
                    if (pnb == heavy) vId[nb] -= ew;
                    if (pnb == bestDst) vId[nb] += ew;
                }
                vId[bestV] = 0;
                for (int ei = g.xadj[bestV]; ei < g.xadj[bestV + 1]; ++ei)
                    if (g.where[g.adjncy[ei]] == bestDst) vId[bestV] += g.Ewgt(ei);
                g.pwgts[heavy] -= wv;
                g.pwgts[bestDst] += wv;
                g.where[bestV] = bestDst;
                ++moved;
                changed = true;
                continue;
            }
        }

        // 轻←邻居（轻分区从最重邻居拉取一个边界点）
        if (light != -1)
        {
            int bestV = -1, bestSrc = -1;
            real_t bestScore = -1e30;
            for (int v = 0; v < g.nvtxs; ++v)
            {
                int pv = g.where[v];
                if (pv == light) continue;
                if (g.pwgts[pv] - g.Vwgt(v) < minPW) continue;
                bool adjLight = false;
                for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
                    if (g.where[g.adjncy[ei]] == light) { adjLight = true; break; }
                if (!adjLight) continue;
                if (g.pwgts[light] + g.Vwgt(v) > maxPW) continue;
                real_t score = (real_t)g.pwgts[pv]; // 越重越优先
                if (score > bestScore) { bestScore = score; bestV = v; bestSrc = pv; }
            }
            if (bestV != -1)
            {
                int wv = g.Vwgt(bestV);
                for (int ei = g.xadj[bestV]; ei < g.xadj[bestV + 1]; ++ei)
                {
                    int nb = g.adjncy[ei], ew = g.Ewgt(ei), pnb = g.where[nb];
                    if (pnb == bestSrc) vId[nb] -= ew;
                    if (pnb == light) vId[nb] += ew;
                }
                vId[bestV] = 0;
                for (int ei = g.xadj[bestV]; ei < g.xadj[bestV + 1]; ++ei)
                    if (g.where[g.adjncy[ei]] == light) vId[bestV] += g.Ewgt(ei);
                g.pwgts[bestSrc] -= wv;
                g.pwgts[light] += wv;
                g.where[bestV] = light;
                ++moved;
                changed = true;
            }
        }
    }

    // 重算 mincut
    int cut = 0;
    for (int v = 0; v < g.nvtxs; ++v)
        for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
            if (g.where[v] != g.where[g.adjncy[ei]]) cut += g.Ewgt(ei);
    g.mincut = cut / 2;

    return moved;
}
