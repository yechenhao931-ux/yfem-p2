#include "partition/kwayrefine.hpp"
#include <algorithm>
#include <climits>
#include <cstdio>
#include <numeric>
#include <vector>

// ── Bucket 队列 ───────────────────────────────────────────────
struct KBQ
{
    int R, top;
    std::vector<std::vector<int>> bkt;
    std::vector<int> pos, gval;
    KBQ() {}
    KBQ(int n, int R_) : R(R_), top(INT_MIN), bkt(2 * R_ + 1), pos(n, -1), gval(n, 0) {}
    bool inQ(int v) const { return pos[v] != -1; }
    void push(int v, int g)
    {
        gval[v] = g;
        pos[v] = (int)bkt[g + R].size();
        bkt[g + R].push_back(v);
        if (g > top)
            top = g;
    }
    void erase(int v)
    {
        if (pos[v] == -1)
            return;
        int idx = gval[v] + R, p = pos[v], last = bkt[idx].back();
        bkt[idx][p] = last;
        pos[last] = p;
        bkt[idx].pop_back();
        pos[v] = -1;
    }
    void update(int v, int dg)
    {
        if (pos[v] == -1)
            return;
        erase(v);
        push(v, gval[v] + dg);
    }
    int popMax()
    {
        while (top >= -R && bkt[top + R].empty())
            --top;
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
            --top;
        return top < -R;
    }
};

// ================================================================
//  ComputeCkrinfo
// ================================================================
void ComputeCkrinfo(Graph &g, int nparts)
{
    const int n = g.nvtxs;
    g.AllocateCkrinfo(std::min(nparts, 16));
    g.pwgts.assign(nparts, 0);
    g.bndptr.assign(n, -1);
    g.bndind.clear();
    g.nbnd = 0;
    g.mincut = 0;
    std::vector<int> tmpEd(nparts, 0), seen;
    seen.reserve(nparts);

    for (int v = 0; v < n; ++v)
        g.pwgts[g.where[v]] += g.Vwgt(v);

    for (int v = 0; v < n; ++v)
    {
        int pv = g.where[v];
        Ckrinfo &ci = g.ckrinfo[v];
        ci.id = ci.ed = ci.nnbrs = 0;
        seen.clear();
        for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
        {
            int u = g.adjncy[ei], ew = g.Ewgt(ei), pu = g.where[u];
            if (pu == pv)
                ci.id += ew;
            else
            {
                ci.ed += ew;
                if (tmpEd[pu] == 0)
                    seen.push_back(pu);
                tmpEd[pu] += ew;
            }
        }
        ci.inbr = (int)g.cnbrPool.size();
        ci.nnbrs = (int)seen.size();
        for (int pu : seen)
        {
            g.cnbrPool.push_back({pu, tmpEd[pu]});
            tmpEd[pu] = 0;
        }
        if (ci.ed > 0)
        {
            g.mincut += ci.ed;
            g.bndptr[v] = g.nbnd++;
            g.bndind.push_back(v);
        }
    }
    g.mincut /= 2;
}

// ================================================================
//  ComputeVkrinfo
// ================================================================
void ComputeVkrinfo(Graph &g, int nparts)
{
    const int n = g.nvtxs;
    g.AllocateVkrinfo(std::min(nparts, 16));
    g.minvol = 0;
    std::vector<int> tmpEd(nparts, 0), tmpGv(nparts, 0), seen;
    seen.reserve(nparts);

    for (int v = 0; v < n; ++v)
    {
        int pv = g.where[v];
        Vkrinfo &vi = g.vkrinfo[v];
        const Ckrinfo &ci = g.ckrinfo[v];
        vi.nid = ci.id;
        vi.ned = (ci.ed > 0) ? g.Vsize(v) : 0;
        vi.gv = 0;
        seen.clear();
        if (ci.ed == 0)
        {
            vi.inbr = (int)g.vnbrPool.size();
            vi.nnbrs = 0;
            continue;
        }
        for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
        {
            int u = g.adjncy[ei], ew = g.Ewgt(ei), pu = g.where[u];
            if (pu == pv)
                continue;
            if (tmpEd[pu] == 0)
                seen.push_back(pu);
            tmpEd[pu] += ew;
            if (g.ckrinfo[u].ed == ew)
                tmpGv[pu] += g.Vsize(u);
        }
        bool becomes = (ci.id == 0);
        vi.inbr = (int)g.vnbrPool.size();
        vi.nnbrs = (int)seen.size();
        // 标准 METIS TOTALV 口径：顶点 v 需把自身数据发给它触及的每个外部
        // 分区，故通信量贡献 = vsize(v) × (相邻分区数)。旧实现只 +vsize(v)
        // 一次（等价于"界面顶点数"），会低估真实通信量。
        g.minvol += g.Vsize(v) * vi.nnbrs;
        for (int pu : seen)
        {
            int gv = tmpGv[pu] + (becomes ? g.Vsize(v) : 0);
            g.vnbrPool.push_back({pu, tmpEd[pu], gv});
            vi.gv = std::max(vi.gv, gv);
            tmpEd[pu] = 0;
            tmpGv[pu] = 0;
        }
    }
}

// ================================================================
//  KwayFMCut  ──  带回滚的 k-way FM 切边精化
//
//  每轮（iter）：
//    · 将所有边界点加入 Bucket 队列
//    · 贪心选最高增益点移动（锁定，不重复移动）
//    · 记录移动历史和累计增益，追踪历史最优位置
//    · 本轮结束后：回滚到历史最优点，全量重建 Ckrinfo
//    · 净增益 <= 0 时停止
// ================================================================
int KwayFMCut(Graph &g, const KwayFMOpts &opts)
{
    const int n = g.nvtxs, K = opts.nparts;
    int W = 0;
    for (int p = 0; p < K; ++p)
        W += g.pwgts[p];
    real_t ideal = (real_t)W / K;
    std::vector<int> maxPW(K);
    for (int p = 0; p < K; ++p)
        maxPW[p] = (int)(opts.ubFactor * ideal + 0.5f);

    // Bucket 范围
    int R = 1;
    for (int v = 0; v < n; ++v)
    {
        int d = 0;
        for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
            d += g.Ewgt(ei);
        if (d > R)
            R = d;
    }

    // ── Cnbr Pool 操作（不缓存跨 push_back 的指针）────────────
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

    // ── 计算每顶点最优可行移动 ────────────────────────────────
    std::vector<int> bestPart(n, -1), bestGain(n, INT_MIN);
    auto calcBest = [&](int v)
    {
        bestPart[v] = -1;
        bestGain[v] = INT_MIN;
        const Ckrinfo &ci = g.ckrinfo[v];
        if (ci.ed == 0 || ci.inbr < 0)
            return;
        int wv = g.Vwgt(v);
        for (int i = 0; i < ci.nnbrs; ++i)
        {
            int p = g.cnbrPool[ci.inbr + i].pid;
            int gv = g.cnbrPool[ci.inbr + i].ed - ci.id;
            if (g.pwgts[p] + wv <= maxPW[p] && gv > bestGain[v])
            {
                bestGain[v] = gv;
                bestPart[v] = p;
            }
        }
    };

    // ── 移动历史（用于回滚）──────────────────────────────────
    struct MoveRec
    {
        int v, src, dst;
    };
    std::vector<MoveRec> hist;
    hist.reserve(n);

    int totalGain = 0;

    for (int iter = 0; iter < opts.nIter; ++iter)
    {
        // 初始化队列
        KBQ bq(n, R);
        std::vector<bool> locked(n, false);
        for (int i = 0; i < g.nbnd; ++i)
        {
            int v = g.bndind[i];
            calcBest(v);
            if (bestPart[v] != -1)
                bq.push(v, bestGain[v]);
        }

        hist.clear();
        int cum = 0, best = 0, bestStep = -1;

        // ── FM 贪心移动 ───────────────────────────────────────
        while (!bq.empty())
        {
            int v = bq.popMax();
            if (locked[v])
                continue;
            int dst = bestPart[v];
            if (dst == -1)
                continue;
            int wv = g.Vwgt(v);
            // 平衡约束可能因他人移动而失效
            if (g.pwgts[dst] + wv > maxPW[dst])
            {
                calcBest(v);
                if (bestPart[v] != -1)
                    bq.push(v, bestGain[v]);
                continue;
            }

            int src = g.where[v];
            // 精确增益（从 Cnbr 查）
            int moveGain = 0;
            {
                int idx = findC(v, dst);
                if (idx >= 0)
                    moveGain = g.cnbrPool[g.ckrinfo[v].inbr + idx].ed - g.ckrinfo[v].id;
            }

            // 执行移动
            g.where[v] = dst;
            g.pwgts[src] -= wv;
            g.pwgts[dst] += wv;
            locked[v] = true;
            hist.push_back({v, src, dst});
            cum += moveGain;
            if (cum > best)
            {
                best = cum;
                bestStep = (int)hist.size() - 1;
            }

            // 重建 v 的 Ckrinfo（重扫描，写入 pool 末尾）
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
                    bq.update(u, bestGain[u] - bq.gval[u]);
                else if (was)
                    bq.erase(u);
                else if (need)
                    bq.push(u, bestGain[u]);
            }
        } // end FM loop

        // ── 回滚到历史最优点 ──────────────────────────────────
        for (int i = (int)hist.size() - 1; i > bestStep; --i)
        {
            auto &m = hist[i];
            int wv = g.Vwgt(m.v);
            g.where[m.v] = m.src;
            g.pwgts[m.src] += wv;
            g.pwgts[m.dst] -= wv;
        }

        // 全量重建（保证一致性）
        ComputeCkrinfo(g, K);

        totalGain += best;
        if (opts.verbose)
            std::printf("  [KwayFMCut iter%d] gain=%d  cut=%d\n", iter, best, g.mincut);
        if (best <= 0)
            break;
    }
    return totalGain;
}

// ================================================================
//  KwayFMVol  ──  每轮全量重建 + 贪心移动（带回滚）
// ================================================================
int KwayFMVol(Graph &g, const KwayFMOpts &opts)
{
    const int n = g.nvtxs, K = opts.nparts;
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
    maxVol = std::max(maxVol, 1);

    struct MoveRec
    {
        int v, src, dst;
    };
    std::vector<MoveRec> hist;
    int totalGain = 0;

    for (int iter = 0; iter < opts.nIter; ++iter)
    {
        ComputeCkrinfo(g, K);
        ComputeVkrinfo(g, K);

        KBQ bq(n, maxVol);
        std::vector<bool> locked(n, false);
        std::vector<int> bPart(n, -1), bGv(n, INT_MIN);

        auto calcV = [&](int v)
        {
            bPart[v] = -1;
            bGv[v] = INT_MIN;
            Vkrinfo &vi = g.vkrinfo[v];
            if (vi.ned == 0 || vi.inbr < 0)
                return;
            int wv = g.Vwgt(v);
            for (int i = 0; i < vi.nnbrs; ++i)
            {
                const Vnbr &nb = g.vnbrPool[vi.inbr + i];
                if (g.pwgts[nb.pid] + wv <= maxPW[nb.pid] && nb.gv > bGv[v])
                {
                    bGv[v] = nb.gv;
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
        int cum = 0, best = 0, bestStep = -1;

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
            int src = g.where[v], gain = bGv[v];
            g.where[v] = dst;
            g.pwgts[src] -= wv;
            g.pwgts[dst] += wv;
            locked[v] = true;
            hist.push_back({v, src, dst});
            cum += gain;
            if (cum > best)
            {
                best = cum;
                bestStep = (int)hist.size() - 1;
            }
        }

        // 回滚
        for (int i = (int)hist.size() - 1; i > bestStep; --i)
        {
            auto &m = hist[i];
            int wv = g.Vwgt(m.v);
            g.where[m.v] = m.src;
            g.pwgts[m.src] += wv;
            g.pwgts[m.dst] -= wv;
        }

        ComputeCkrinfo(g, K);
        ComputeVkrinfo(g, K);
        totalGain += best;
        if (opts.verbose)
            std::printf("  [KwayFMVol iter%d] gain=%d  vol=%d\n", iter, best, g.minvol);
        if (best <= 0)
            break;
    }
    return totalGain;
}
