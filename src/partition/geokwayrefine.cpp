/**
 * GeoKwayRefine.cpp  ──  几何感知 k-way FM 精化完整实现
 *
 * 核心算法见 GeoKwayRefine.h 注释。
 * 本文件实现：
 *   1. autoBeta 标定（RMS 比值法）
 *   2. 混合增益 Bucket 队列管理
 *   3. 增量重心更新
 *   4. 历史回滚（同 KwayRefine.cpp）
 */
#include <common/macro.hpp>
#include <partition/geokwayrefine.hpp>
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <vector>

// ── 几何辅助函数 ──────────────────────────────────────────────
static real_t dist2C(const Coord &a, const std::vector<real_t> &c)
{
    real_t dx = a.x - c[0], dy = a.y - c[1], dz = a.z - c[2];
    return dx * dx + dy * dy + dz * dz;
}

// ── Bucket 队列（支持浮点增益映射到整型桶）─────────────────────
// 将浮点混合增益离散化：乘以 SCALE 后取 int
// SCALE 决定精度；取 1000 可区分 0.001 的差异
static constexpr int GEO_SCALE = 1000;

struct GeoBQ
{
    int R, top;
    std::vector<std::vector<int>> bkt;
    std::vector<int> pos;
    std::vector<real_t> gval; // 原始浮点增益（用于回滚时记录精确值）
    std::vector<int> gint;    // 离散化整型增益（用于桶索引）

    GeoBQ() {}
    GeoBQ(int n, int R_)
        : R(R_), top(INT_MIN), bkt(2 * R_ + 1), pos(n, -1),
          gval(n, 0.0), gint(n, 0) {}

    bool inQ(int v) const { return pos[v] != -1; }

    void push(int v, real_t gain)
    {
        gval[v] = gain;
        int gi = std::max(-R, std::min(R, (int)(gain * GEO_SCALE)));
        gint[v] = gi;
        pos[v] = (int)bkt[gi + R].size();
        bkt[gi + R].push_back(v);
        if (gi > top)
        {
            top = gi;
        }
    }

    void erase(int v)
    {
        if (pos[v] == -1)
            return;
        int idx = gint[v] + R, p = pos[v];
        int last = bkt[idx].back();
        bkt[idx][p] = last;
        pos[last] = p;
        bkt[idx].pop_back();
        pos[v] = -1;
    }

    void update(int v, real_t newg)
    {
        if (pos[v] == -1)
        {
            return;
        }
        erase(v);
        push(v, newg);
    }

    int popMax()
    {
        while (top >= -R && bkt[top + R].empty())
        {
            top--;
        }
        if (top < -R)
            return -1;
        auto &b = bkt[top + R];
        int v = b.back();
        b.pop_back();
        pos[v] = -1;
        return v;
    }

    bool empty()
    {
        while (top >= -R && bkt[top + R].empty())
        {
            top--;
        }

        return top < -R;
    }
};

// 增量更新 src/dst 两个分区的重心（v 从 src 移到 dst）
static void UpdateCentroid(const Graph &g, std::vector<std::vector<real_t>> &cen, const std::vector<int> &pwgts, int v, int src, int dst)
{
    real_t wv = (real_t)g.Vwgt(v);
    real_t wsrc = (real_t)pwgts[src]; // 移动后的 src 权重
    real_t wdst = (real_t)pwgts[dst]; // 移动后的 dst 权重

    // src 重心：移除 v 的贡献
    real_t ws_old = wsrc + wv;
    if (ws_old > 0.0)
    {
        for (int d = 0; d < 3; d++)
        {
            real_t coord = (d == 0) ? g.coordinates[v].x : (d == 1) ? g.coordinates[v].y
                                                                    : g.coordinates[v].z;
            cen[src][d] = (wsrc > 0.0)
                              ? (cen[src][d] * ws_old - wv * coord) / wsrc
                              : 0.0;
        }
    }

    // dst 重心：加入 v 的贡献
    real_t wd_old = wdst - wv;
    for (int d = 0; d < 3; ++d)
    {
        real_t coord = (d == 0) ? g.coordinates[v].x : (d == 1) ? g.coordinates[v].y
                                                                : g.coordinates[v].z;
        cen[dst][d] = (wdst > 0.0)
                          ? (cen[dst][d] * wd_old + wv * coord) / wdst
                          : coord;
    }
}

// ================================================================
//  AutoBeta 标定
//  对所有边界顶点，采样计算 topo_gain 和 geo_gain 的 RMS，
//  令 beta = rms(topo) / rms(geo)，使两者量级相当。
// ================================================================
static real_t CalibrateBeta(const Graph &g, int K, const std::vector<std::vector<real_t>> &cen)
{
    real_t sumT2 = 0, sumG2 = 0;
    int cnt = 0;
    for (int i = 0; i < g.nbnd; i++)
    {
        int v = g.bndind[i];
        const Ckrinfo &ci = g.ckrinfo[v];
        if (ci.inbr < 0 || ci.nnbrs == 0)
            continue;

        int src = g.where[v];

        for (int j = 0; j < ci.nnbrs; j++)
        {
            int p = g.cnbrPool[ci.inbr + j].pid;
            real_t tg = (real_t)(g.cnbrPool[ci.inbr + j].ed - ci.id);
            real_t gg = dist2C(g.coordinates[v], cen[src]) - dist2C(g.coordinates[v], cen[p]);

            sumT2 += tg * tg;
            sumG2 += gg * gg;
            cnt++;
        }
    }

    if (cnt == 0 || sumG2 < 1e-12)
        return 1.0;
    return (real_t)std::sqrt(sumT2 / sumG2);
}

// ================================================================
//  GeoKwayFMCut  ──  几何感知 k-way FM 切边精化
// ================================================================
GeoFMResult GeoKwayFMCut(Graph &g, const GeoFMOpts &opts, std::vector<std::vector<real_t>> &centroids)
{
    YFEM_ASSERT(g.coordinates.size() == g.nvtxs, "GeoKwayFMCut: 需要 graph.coordinates[]");

    const int n = g.nvtxs, K = opts.nparts;
    real_t alpha = opts.alpha;

    // 平衡约束
    int W = 0;
    for (int p = 0; p < K; p++)
    {
        W += g.pwgts[p];
    }

    real_t ideal = (real_t)W / K;

    std::vector<int> maxPW(K);

    for (int p = 0; p < K; ++p)
    {
        maxPW[p] = (int)(opts.ubFactor * ideal + 0.5);
    }

    // Bucket 范围（混合增益上界估计）
    // topo 项最大 = 最大加权度；geo 项最大 = bbox 对角线²
    int maxTopoR = 1;
    for (int v = 0; v < n; v++)
    {
        int d = 0;
        for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ei++)
        {
            d += g.Ewgt(ei);
        }
        maxTopoR = std::max(maxTopoR, d);
    }
    // 几何增益的安全上界：用顶点坐标 bbox 对角线²。
    // 这是 O(n) 而非 O(K²)，且作为更紧的上界
    // （centroid 可能位于 bbox 内部，差异远小于 bbox 对角线）。
    real_t minX = std::numeric_limits<real_t>::max(), maxX = -minX;
    real_t minY = minX, maxY = -minX;
    real_t minZ = minX, maxZ = -minX;
    for (int v = 0; v < n; ++v)
    {
        const Coord &c = g.coordinates[v];
        if (c.x < minX) minX = c.x;
        if (c.x > maxX) maxX = c.x;
        if (c.y < minY) minY = c.y;
        if (c.y > maxY) maxY = c.y;
        if (c.z < minZ) minZ = c.z;
        if (c.z > maxZ) maxZ = c.z;
    }
    real_t bboxD2 = (maxX - minX) * (maxX - minX)
                  + (maxY - minY) * (maxY - minY)
                  + (maxZ - minZ) * (maxZ - minZ);
    if (bboxD2 < 1.0) bboxD2 = 1.0;

    // 自动标定 beta
    real_t beta = opts.autoBeta ? CalibrateBeta(g, K, centroids) : opts.beta;
    if (opts.verbose)
        std::printf("  [GeoFM] alpha=%.2f  beta=%.4f\n", alpha, beta);

    // Bucket 范围（以整型桶为单位，混合增益乘以 GEO_SCALE 后的最大值）
    int R = (int)((alpha * maxTopoR + (1 - alpha) * beta * std::sqrt(bboxD2) + 1.0) * GEO_SCALE + 1);
    R = std::max(R, 1);

    // Cnbr Pool 操作（与 KwayRefine 相同）
    auto findC = [&](int u, int p) -> int
    {
        int inbr = g.ckrinfo[u].inbr, nn = g.ckrinfo[u].nnbrs;
        if (inbr < 0)
            return -1;
        for (int i = 0; i < nn; ++i)
            if (g.cnbrPool[inbr + i].pid == p)
                return i;
        return -1;
    };
    auto removeC = [&](int u, int idx)
    {
        int inbr = g.ckrinfo[u].inbr;
        int &nn = g.ckrinfo[u].nnbrs;
        g.cnbrPool[inbr + idx] = g.cnbrPool[inbr + nn - 1];
        --nn;
    };
    auto addC = [&](int u, int pid, int ev)
    {
        Ckrinfo &cu = g.ckrinfo[u];
        int newInbr = (int)g.cnbrPool.size();
        for (int i = 0; i < cu.nnbrs; ++i)
            g.cnbrPool.push_back(g.cnbrPool[cu.inbr + i]);
        g.cnbrPool.push_back({pid, ev});
        cu.inbr = newInbr;
        cu.nnbrs++;
    };

    // 混合增益计算（核心公式）
    // gain_hybrid(v→p) = alpha*topo_gain + (1-alpha)*beta*geo_gain
    auto hybridGain = [&](int v, int dst) -> real_t
    {
        const Ckrinfo &ci = g.ckrinfo[v];
        int src = g.where[v];
        // 拓扑增益：找 dst 的 Cnbr
        real_t topo = 0.;
        for (int i = 0; i < ci.nnbrs; i++)
        {
            if (g.cnbrPool[ci.inbr + i].pid == dst)
            {
                topo = (real_t)(g.cnbrPool[ci.inbr + i].ed - ci.id);
                break;
            }
        }

        // 几何增益：移到 dst 后距离减少
        real_t geo = dist2C(g.coordinates[v], centroids[src]) - dist2C(g.coordinates[v], centroids[dst]);

        return alpha * topo + (1.0 - alpha) * beta * geo;
    };

    // 计算每顶点最优可行混合增益
    std::vector<int> bestPart(n, -1);
    std::vector<real_t> bestGain(n, -1e30);
    auto calcBest = [&](int v)
    {
        bestPart[v] = -1;
        bestGain[v] = -1e30;
        const Ckrinfo &ci = g.ckrinfo[v];
        if (ci.ed == 0 || ci.inbr < 0)
        {
            return;
        }
        int wv = g.Vwgt(v);
        // 拓扑：只考虑有拓扑邻居的分区（保证图连通性）
        for (int i = 0; i < ci.nnbrs; i++)
        {
            int p = g.cnbrPool[ci.inbr + i].pid;
            if (g.pwgts[p] + wv > maxPW[p])
                continue;
            real_t hg = hybridGain(v, p);
            if (hg > bestGain[v])
            {
                bestGain[v] = hg;
                bestPart[v] = p;
            }
        }
    };
    // 移动历史（用于回滚
    struct MoveRec
    {
        int v, src, dst;
    };
    std::vector<MoveRec> hist;
    hist.reserve(n);

    GeoFMResult res;
    res.betaUsed = beta;
    int totalTopoCut = 0;

    for (int iter = 0; iter < opts.nIter; ++iter)
    {
        GeoBQ bq(n, R);
        std::vector<bool> locked(n, false);
        for (int i = 0; i < g.nbnd; ++i)
        {
            int v = g.bndind[i];
            calcBest(v);
            if (bestPart[v] != -1)
                bq.push(v, bestGain[v]);
        }

        hist.clear();
        real_t cum = 0.f, best = 0.f;
        int bestStep = -1;
        int stepTopoCut = 0, bestTopoCut = 0;

        while (!bq.empty())
        {
            int v = bq.popMax();
            if (locked[v])
                continue;
            int dst = bestPart[v];
            if (dst == -1)
                continue;
            int wv = g.Vwgt(v);
            if (g.pwgts[dst] + wv > maxPW[dst])
            {
                calcBest(v);
                if (bestPart[v] != -1)
                    bq.push(v, bestGain[v]);
                continue;
            }

            int src = g.where[v];
            real_t hg = hybridGain(v, dst);

            // 精确拓扑增益
            int topoGain = 0;
            {
                int idx = findC(v, dst);
                if (idx >= 0)
                    topoGain = g.cnbrPool[g.ckrinfo[v].inbr + idx].ed - g.ckrinfo[v].id;
            }

            // 执行移动
            g.where[v] = dst;
            g.pwgts[src] -= wv;
            g.pwgts[dst] += wv;
            locked[v] = true;
            hist.push_back({v, src, dst});
            // 增量更新重心
            UpdateCentroid(g, centroids, g.pwgts, v, src, dst);
            cum += hg;
            stepTopoCut += topoGain;
            if (cum > best)
            {
                best = cum;
                bestStep = (int)hist.size() - 1;
                bestTopoCut = stepTopoCut;
            }
            ++res.nMoves;

            // 重建 v 的 Ckrinfo
            {
                Ckrinfo &cv = g.ckrinfo[v];
                std::vector<int> tmpE(K, 0), sp;
                sp.reserve(K);
                cv.id = cv.ed = 0;
                for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
                {
                    int u2 = g.adjncy[ei], ew = g.Ewgt(ei), pu = g.where[u2];
                    if (pu == dst)
                        cv.id += ew;
                    else
                    {
                        cv.ed += ew;
                        if (tmpE[pu] == 0)
                            sp.push_back(pu);
                        tmpE[pu] += ew;
                    }
                }
                cv.inbr = (int)g.cnbrPool.size();
                cv.nnbrs = (int)sp.size();
                for (int p : sp)
                {
                    g.cnbrPool.push_back({p, tmpE[p]});
                    tmpE[p] = 0;
                }
            }

            // 增量更新邻居 Ckrinfo
            for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
            {
                int u = g.adjncy[ei], ew = g.Ewgt(ei), pu = g.where[u];
                if (locked[u])
                    continue;
                Ckrinfo &cu = g.ckrinfo[u];
                if (pu == src)
                {
                    cu.id -= ew;
                    cu.ed += ew;
                    int di = findC(u, dst);
                    if (di >= 0)
                        g.cnbrPool[cu.inbr + di].ed += ew;
                    else
                        addC(u, dst, ew);
                }
                else if (pu == dst)
                {
                    cu.id += ew;
                    cu.ed -= ew;
                    int si = findC(u, src);
                    if (si >= 0)
                    {
                        g.cnbrPool[cu.inbr + si].ed -= ew;
                        if (g.cnbrPool[cu.inbr + si].ed == 0)
                            removeC(u, si);
                    }
                }
                else
                {
                    int si = findC(u, src);
                    if (si >= 0)
                    {
                        g.cnbrPool[cu.inbr + si].ed -= ew;
                        if (g.cnbrPool[cu.inbr + si].ed == 0)
                            removeC(u, si);
                    }
                    int di = findC(u, dst);
                    if (di >= 0)
                        g.cnbrPool[cu.inbr + di].ed += ew;
                    else
                        addC(u, dst, ew);
                }
                calcBest(u);
                bool was = bq.inQ(u), need = (cu.ed > 0 && bestPart[u] != -1);
                if (was && need)
                    bq.update(u, bestGain[u]);
                else if (was)
                    bq.erase(u);
                else if (need)
                    bq.push(u, bestGain[u]);
            }
        } // FM loop

        // 回滚到历史最优
        for (int i = (int)hist.size() - 1; i > bestStep; --i)
        {
            auto &m = hist[i];
            int wv = g.Vwgt(m.v);
            // 回滚重心（逆向更新）
            g.where[m.v] = m.src;
            g.pwgts[m.src] += wv;
            g.pwgts[m.dst] -= wv;
            UpdateCentroid(g, centroids, g.pwgts, m.v, m.dst, m.src);
        }
        // 全量重建（保证一致性）
        ComputeCkrinfo(g, K);
        // 同步重建重心
        centroids = ComputeCentroids(g, K);

        totalTopoCut += bestTopoCut;
        if (opts.verbose)
            std::printf("  [GeoFM iter%d] hybridGain=%.3f  topoCut=%d  cut=%d\n",
                        iter, best, bestTopoCut, g.mincut);
        if (best <= 0.0)
            break;
    }

    res.gainTopo = totalTopoCut;
    return res;
}

// ================================================================
//  GeoKwayFMVol  ──  几何感知 Volume FM
//  在 KwayFMVol 的通信量增益基础上叠加几何项
// ================================================================
GeoFMResult GeoKwayFMVol(Graph &g, const GeoFMOpts &opts,
                         std::vector<std::vector<real_t>> &centroids)
{
    assert((int)g.coordinates.size() == g.nvtxs);
    const int n = g.nvtxs, K = opts.nparts;
    real_t alpha = opts.alpha;

    int W = 0;
    for (int p = 0; p < K; ++p)
        W += g.pwgts[p];
    real_t ideal = (real_t)W / K;
    std::vector<int> maxPW(K);
    for (int p = 0; p < K; ++p)
        maxPW[p] = (int)(opts.ubFactor * ideal + 0.5f);

    int maxVol = 0;
    for (int v = 0; v < n; ++v)
        maxVol += g.Vsize(v);

    real_t beta = opts.autoBeta ? CalibrateBeta(g, K, centroids) : opts.beta;
    int R = std::max((int)((alpha * maxVol + (1 - alpha) * beta * 100.f) * GEO_SCALE + 1), 1);

    auto geoGain = [&](int v, int dst) -> real_t
    {
        int src = g.where[v];
        return dist2C(g.coordinates[v], centroids[src]) - dist2C(g.coordinates[v], centroids[dst]);
    };

    struct MoveRec
    {
        int v, src, dst;
    };
    std::vector<MoveRec> hist;
    GeoFMResult res;
    res.betaUsed = beta;

    for (int iter = 0; iter < opts.nIter; ++iter)
    {
        ComputeCkrinfo(g, K);
        ComputeVkrinfo(g, K);
        centroids = ComputeCentroids(g, K);

        GeoBQ bq(n, R);
        std::vector<bool> locked(n, false);
        std::vector<int> bPart(n, -1);
        std::vector<real_t> bGv(n, -1e30f);

        auto calcV = [&](int v)
        {
            bPart[v] = -1;
            bGv[v] = -1e30f;
            Vkrinfo &vi = g.vkrinfo[v];
            if (vi.ned == 0 || vi.inbr < 0)
                return;
            int wv = g.Vwgt(v);
            int src = g.where[v];
            for (int i = 0; i < vi.nnbrs; ++i)
            {
                const Vnbr &nb = g.vnbrPool[vi.inbr + i];
                if (g.pwgts[nb.pid] + wv > maxPW[nb.pid])
                    continue;
                real_t hg = alpha * nb.gv + (1 - alpha) * beta * geoGain(v, nb.pid);
                if (hg > bGv[v])
                {
                    bGv[v] = hg;
                    bPart[v] = nb.pid;
                }
            }
        };
        for (int i = 0; i < g.nbnd; ++i)
        {
            int v = g.bndind[i];
            calcV(v);
            if (bPart[v] != -1)
                bq.push(v, bGv[v]);
        }

        hist.clear();
        real_t cum = 0, best = 0;
        int bestStep = -1;
        while (!bq.empty())
        {
            int v = bq.popMax();
            if (locked[v])
                continue;
            int dst = bPart[v];
            if (dst == -1)
                continue;
            int wv = g.Vwgt(v);
            if (g.pwgts[dst] + wv > maxPW[dst])
            {
                calcV(v);
                if (bPart[v] != -1)
                    bq.push(v, bGv[v]);
                continue;
            }
            int src = g.where[v];
            real_t hg = bGv[v];
            g.where[v] = dst;
            g.pwgts[src] -= wv;
            g.pwgts[dst] += wv;
            locked[v] = true;
            hist.push_back({v, src, dst});
            UpdateCentroid(g, centroids, g.pwgts, v, src, dst);
            cum += hg;
            if (cum > best)
            {
                best = cum;
                bestStep = (int)hist.size() - 1;
            }
            ++res.nMoves;
        }
        for (int i = (int)hist.size() - 1; i > bestStep; --i)
        {
            auto &m = hist[i];
            int wv = g.Vwgt(m.v);
            g.where[m.v] = m.src;
            g.pwgts[m.src] += wv;
            g.pwgts[m.dst] -= wv;
            UpdateCentroid(g, centroids, g.pwgts, m.v, m.dst, m.src);
        }
        ComputeCkrinfo(g, K);
        ComputeVkrinfo(g, K);
        centroids = ComputeCentroids(g, K);
        if (opts.verbose)
            std::printf("  [GeoFMVol iter%d] hybridGain=%.3f  vol=%d\n", iter, best, g.minvol);
        if (best <= 0.f)
            break;
    }
    return res;
}

// OPT1
// ================================================================
//  GeoKwayFMCut_OPT1  ──  几何感知 k-way FM 切边精化
// ================================================================
GeoFMResult GeoKwayFMCut_Opt1(Graph &g, const GeoFMOpts &opts, std::vector<std::vector<real_t>> &centroids)
{

    YFEM_ASSERT(g.coordinates.size() == g.nvtxs, "GeoKwayFMCut: 需要 graph.coordinates[]");
    const int n = g.nvtxs, K = opts.nparts;
    real_t alpha = opts.alpha;

    // 平衡约束
    int W = 0;
    for (int p = 0; p < K; p++)
    {
        W += g.pwgts[p];
    }

    real_t ideal = (real_t)W / K;

    std::vector<int> maxPW(K);

    for (int p = 0; p < K; ++p)
    {
        maxPW[p] = (int)(opts.ubFactor * ideal + 0.5);
    }

    // Bucket 范围（混合增益上界估计）
    // topo 项最大 = 最大加权度；geo 项最大 = bbox 对角线²（O(n)）
    int maxTopoR = 1;
    for (int v = 0; v < n; v++)
    {
        int d = 0;
        for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ei++)
        {
            d += g.Ewgt(ei);
        }
        maxTopoR = std::max(maxTopoR, d);
    }
    // 用顶点 bbox 而不是 centroid bbox（更紧的上界，O(n) 而非 O(K²)）
    real_t minX = std::numeric_limits<real_t>::max(), maxX = -minX;
    real_t minY = minX, maxY = -minX;
    real_t minZ = minX, maxZ = -minX;
    for (int v = 0; v < n; ++v)
    {
        const Coord &c = g.coordinates[v];
        if (c.x < minX) minX = c.x;
        if (c.x > maxX) maxX = c.x;
        if (c.y < minY) minY = c.y;
        if (c.y > maxY) maxY = c.y;
        if (c.z < minZ) minZ = c.z;
        if (c.z > maxZ) maxZ = c.z;
    }
    real_t bboxD2 = (maxX - minX) * (maxX - minX)
                  + (maxY - minY) * (maxY - minY)
                  + (maxZ - minZ) * (maxZ - minZ);
    if (bboxD2 < 1.0) bboxD2 = 1.0;

    // 自动标定 beta
    real_t beta = opts.autoBeta ? CalibrateBeta(g, K, centroids) : opts.beta;
    if (opts.verbose)
        std::printf("  [GeoFM] alpha=%.2f  beta=%.4f\n", alpha, beta);

    // Bucket 范围（以整型桶为单位，混合增益乘以 GEO_SCALE 后的最大值）
    int R = (int)((alpha * maxTopoR + (1 - alpha) * beta * std::sqrt(bboxD2) + 1.0) * GEO_SCALE + 1);
    R = std::max(R, 1);

    // Cnbr Pool 操作（与 KwayRefine 相同）
    auto findC = [&](int u, int p) -> int
    {
        int inbr = g.ckrinfo[u].inbr, nn = g.ckrinfo[u].nnbrs;
        if (inbr < 0)
            return -1;
        for (int i = 0; i < nn; ++i)
            if (g.cnbrPool[inbr + i].pid == p)
                return i;
        return -1;
    };
    auto removeC = [&](int u, int idx)
    {
        int inbr = g.ckrinfo[u].inbr;
        int &nn = g.ckrinfo[u].nnbrs;
        g.cnbrPool[inbr + idx] = g.cnbrPool[inbr + nn - 1];
        --nn;
    };
    auto addC = [&](int u, int pid, int ev)
    {
        Ckrinfo &cu = g.ckrinfo[u];
        int newInbr = (int)g.cnbrPool.size();
        for (int i = 0; i < cu.nnbrs; ++i)
            g.cnbrPool.push_back(g.cnbrPool[cu.inbr + i]);
        g.cnbrPool.push_back({pid, ev});
        cu.inbr = newInbr;
        cu.nnbrs++;
    };

    // 混合增益计算（核心公式）
    // gain_hybrid(v→p) = alpha*topo_gain + (1-alpha)*beta*geo_gain
    auto hybridGain = [&](int v, int dst) -> real_t
    {
        const Ckrinfo &ci = g.ckrinfo[v];
        int src = g.where[v];
        // 拓扑增益：找 dst 的 Cnbr
        real_t topo = 0.;
        for (int i = 0; i < ci.nnbrs; i++)
        {
            if (g.cnbrPool[ci.inbr + i].pid == dst)
            {
                topo = (real_t)(g.cnbrPool[ci.inbr + i].ed - ci.id);
                break;
            }
        }

        // 几何增益：移到 dst 后距离减少
        real_t geo = dist2C(g.coordinates[v], centroids[src]) - dist2C(g.coordinates[v], centroids[dst]);

        return alpha * topo + (1.0 - alpha) * beta * geo;
    };

    // 计算每顶点最优可行混合增益
    std::vector<int> bestPart(n, -1);
    std::vector<real_t> bestGain(n, -1e30);
    auto calcBest = [&](int v)
    {
        bestPart[v] = -1;
        bestGain[v] = -1e30;
        const Ckrinfo &ci = g.ckrinfo[v];
        if (ci.ed == 0 || ci.inbr < 0)
        {
            return;
        }
        int wv = g.Vwgt(v);
        // 拓扑：只考虑有拓扑邻居的分区（保证图连通性）
        for (int i = 0; i < ci.nnbrs; i++)
        {
            int p = g.cnbrPool[ci.inbr + i].pid;
            if (g.pwgts[p] + wv > maxPW[p])
                continue;
            real_t hg = hybridGain(v, p);
            if (hg > bestGain[v])
            {
                bestGain[v] = hg;
                bestPart[v] = p;
            }
        }
    };
    // 移动历史（用于回滚
    struct MoveRec
    {
        int v, src, dst;
    };
    std::vector<MoveRec> hist;
    hist.reserve(n);

    GeoFMResult res;
    res.betaUsed = beta;
    int totalTopoCut = 0;

    for (int iter = 0; iter < opts.nIter; ++iter)
    {
        int cutBefore = g.mincut; // 记录本轮开始时的 cut，用于检测振荡
        GeoBQ bq(n, R);
        std::vector<bool> locked(n, false);
        for (int i = 0; i < g.nbnd; ++i)
        {
            int v = g.bndind[i];
            calcBest(v);
            if (bestPart[v] != -1)
                bq.push(v, bestGain[v]);
        }

        hist.clear();
        real_t cum = 0., best = 0.;
        int bestStep = -1;
        int stepTopoCut = 0, bestTopoCut = 0;

        while (!bq.empty())
        {
            int v = bq.popMax();
            if (locked[v])
                continue;
            int dst = bestPart[v];
            if (dst == -1)
                continue;
            int wv = g.Vwgt(v);
            if (g.pwgts[dst] + wv > maxPW[dst])
            {
                calcBest(v);
                if (bestPart[v] != -1)
                    bq.push(v, bestGain[v]);
                continue;
            }

            int src = g.where[v];
            float hg = hybridGain(v, dst);

            // 精确拓扑增益
            int topoGain = 0;
            {
                int idx = findC(v, dst);
                if (idx >= 0)
                    topoGain = g.cnbrPool[g.ckrinfo[v].inbr + idx].ed - g.ckrinfo[v].id;
            }

            // 执行移动
            g.where[v] = dst;
            g.pwgts[src] -= wv;
            g.pwgts[dst] += wv;
            locked[v] = true;
            hist.push_back({v, src, dst});
            // 增量更新重心
            UpdateCentroid(g, centroids, g.pwgts, v, src, dst);
            cum += hg;
            stepTopoCut += topoGain;
            if (cum > best)
            {
                best = cum;
                bestStep = (int)hist.size() - 1;
                bestTopoCut = stepTopoCut;
            }
            ++res.nMoves;

            // 重建 v 的 Ckrinfo
            {
                Ckrinfo &cv = g.ckrinfo[v];
                std::vector<int> tmpE(K, 0), sp;
                sp.reserve(K);
                cv.id = cv.ed = 0;
                for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
                {
                    int u2 = g.adjncy[ei], ew = g.Ewgt(ei), pu = g.where[u2];
                    if (pu == dst)
                        cv.id += ew;
                    else
                    {
                        cv.ed += ew;
                        if (tmpE[pu] == 0)
                            sp.push_back(pu);
                        tmpE[pu] += ew;
                    }
                }
                cv.inbr = (int)g.cnbrPool.size();
                cv.nnbrs = (int)sp.size();
                for (int p : sp)
                {
                    g.cnbrPool.push_back({p, tmpE[p]});
                    tmpE[p] = 0;
                }
            }

            // 增量更新邻居 Ckrinfo
            for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
            {
                int u = g.adjncy[ei], ew = g.Ewgt(ei), pu = g.where[u];
                if (locked[u])
                    continue;
                Ckrinfo &cu = g.ckrinfo[u];
                if (pu == src)
                {
                    cu.id -= ew;
                    cu.ed += ew;
                    int di = findC(u, dst);
                    if (di >= 0)
                        g.cnbrPool[cu.inbr + di].ed += ew;
                    else
                        addC(u, dst, ew);
                }
                else if (pu == dst)
                {
                    cu.id += ew;
                    cu.ed -= ew;
                    int si = findC(u, src);
                    if (si >= 0)
                    {
                        g.cnbrPool[cu.inbr + si].ed -= ew;
                        if (g.cnbrPool[cu.inbr + si].ed == 0)
                            removeC(u, si);
                    }
                }
                else
                {
                    int si = findC(u, src);
                    if (si >= 0)
                    {
                        g.cnbrPool[cu.inbr + si].ed -= ew;
                        if (g.cnbrPool[cu.inbr + si].ed == 0)
                            removeC(u, si);
                    }
                    int di = findC(u, dst);
                    if (di >= 0)
                        g.cnbrPool[cu.inbr + di].ed += ew;
                    else
                        addC(u, dst, ew);
                }
                calcBest(u);
                bool was = bq.inQ(u), need = (cu.ed > 0 && bestPart[u] != -1);
                if (was && need)
                    bq.update(u, bestGain[u]);
                else if (was)
                    bq.erase(u);
                else if (need)
                    bq.push(u, bestGain[u]);
            }
        } // FM loop

        // 回滚到历史最优
        for (int i = (int)hist.size() - 1; i > bestStep; --i)
        {
            auto &m = hist[i];
            int wv = g.Vwgt(m.v);
            // 回滚重心（逆向更新）
            g.where[m.v] = m.src;
            g.pwgts[m.src] += wv;
            g.pwgts[m.dst] -= wv;
            UpdateCentroid(g, centroids, g.pwgts, m.v, m.dst, m.src);
        }

        // 全量重建（保证一致性）
        ComputeCkrinfo(g, K);
        // 同步重建重心
        centroids = ComputeCentroids(g, K);

        totalTopoCut += bestTopoCut;
        int cutAfter = g.mincut;
        if (opts.verbose)
            std::printf("  [GeoFM iter%d] hybridGain=%.3f  topoCut=%d  cut=%d  Δcut=%d\n",
                        iter, best, bestTopoCut, cutAfter, cutBefore - cutAfter);
        // 停止条件：混合增益 ≤ 0 或实际 cut 无改善
        if (best <= 0.f || cutAfter >= cutBefore)
            break;
    }

    res.gainTopo = totalTopoCut;
    return res;
}