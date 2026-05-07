/**
 * BisectFM.cpp
 * 锁定后的顶点不再加回队列，避免 FMPass 无限循环。
 */
#include <partition/bisectfm.hpp>
#include <algorithm>
#include <climits>
#include <cstdio>
#include <numeric>
#include <queue>
#include <random>
#include <vector>

// ─── O(1) Bucket 队列 ─────────────────────────────────────────
struct BQ
{
    int R, top;
    std::vector<std::vector<int>> bkt;
    std::vector<int> pos, gval;
    BQ() {}
    BQ(int n, int R_) : R(R_), top(INT_MIN), bkt(2 * R_ + 1), pos(n, -1), gval(n, 0) {}
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

namespace bisect
{

    // ── 全量重建 id/ed/bndind/pwgts/mincut ────────────────────────
    void ComputeParams(Graph &g)
    {
        const int n = g.nvtxs;
        g.id.assign(n, 0);
        g.ed.assign(n, 0);
        g.bndptr.assign(n, -1);
        g.bndind.clear();
        g.nbnd = 0;
        g.pwgts.assign(2, 0);
        g.mincut = 0;
        for (int v = 0; v < n; ++v)
            g.pwgts[g.where[v]] += g.Vwgt(v);
        for (int v = 0; v < n; ++v)
        {
            int pv = g.where[v];
            for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
            {
                int ew = g.Ewgt(ei), pu = g.where[g.adjncy[ei]];
                if (pu == pv)
                    g.id[v] += ew;
                else
                    g.ed[v] += ew;
            }
            if (g.ed[v] > 0)
            {
                g.mincut += g.ed[v];
                g.bndptr[v] = g.nbnd++;
                g.bndind.push_back(v);
            }
        }
        g.mincut /= 2;
    }

    // ── BFS 双种子初始化 ──────────────────────────────────────────
    void GrowInit(Graph &g, int seed, int t0, int t1)
    {
        const int n = g.nvtxs;
        g.where.assign(n, -1);
        g.pwgts.assign(2, 0);
        std::mt19937 rng(seed);
        auto rv = [&]()
        { return (int)(rng() % (unsigned)n); };
        int s0 = rv(), s1 = rv();
        while (s1 == s0)
            s1 = rv();
        g.where[s0] = 0;
        g.where[s1] = 1;
        int w0 = g.Vwgt(s0), w1 = g.Vwgt(s1);
        std::queue<int> q0, q1;
        q0.push(s0);
        q1.push(s1);
        while (!q0.empty() || !q1.empty())
        {
            auto step = [&](std::queue<int> &q, int part, int &pw, int lim)
            {
                if (q.empty())
                    return;
                int u = q.front();
                q.pop();
                for (int ei = g.xadj[u]; ei < g.xadj[u + 1]; ++ei)
                {
                    int v = g.adjncy[ei];
                    if (g.where[v] == -1 && pw < lim)
                    {
                        g.where[v] = part;
                        pw += g.Vwgt(v);
                        q.push(v);
                    }
                }
            };
            step(q0, 0, w0, t0);
            step(q1, 1, w1, t1);
        }
        // 未触达顶点放到偏差更大的一侧
        for (int v = 0; v < n; ++v)
            if (g.where[v] == -1)
            {
                int wv = g.Vwgt(v);
                real_t d0 = (real_t)(t0 - w0) / std::max(t0, 1), d1 = (real_t)(t1 - w1) / std::max(t1, 1);
                int part = (d0 >= d1) ? 0 : 1;
                g.where[v] = part;
                if (part == 0)
                    w0 += wv;
                else
                    w1 += wv;
            }
        g.pwgts[0] = w0;
        g.pwgts[1] = w1;
    }

    // ── FM 单轮精化 ───────────────────────────────────────────────
    // 关键修复：
    //   1. 锁定顶点 cand 后不加回队列（避免无限循环）
    //   2. 邻居 u 非锁定时才更新队列
    int FMPass(Graph &g, int maxPW0, int maxPW1)
    {
        const int n = g.nvtxs;
        int R = 1;
        for (int v = 0; v < n; ++v)
        {
            int d = 0;
            for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
                d += g.Ewgt(ei);
            if (d > R)
                R = d;
        }

        BQ bq0(n, R), bq1(n, R);
        std::vector<bool> locked(n, false);
        auto gain = [&](int v)
        { return g.ed[v] - g.id[v]; };

        // 初始化：所有边界点入队（未锁定）
        for (int i = 0; i < g.nbnd; ++i)
        {
            int v = g.bndind[i];
            if (g.where[v] == 0)
                bq0.push(v, gain(v));
            else
                bq1.push(v, gain(v));
        }

        struct Mv
        {
            int v, from;
        };
        std::vector<Mv> hist;
        hist.reserve(n / 2);
        int cum = 0, best = 0, bestStep = -1;

        for (;;)
        {
            // 从每侧最高桶找第一个满足平衡约束的未锁定候选
            int cand = -1, cFrom = -1, cGain = INT_MIN;
            auto tryPick = [&](BQ &bq, int from, int maxTo)
            {
                // 推进 top 到最高非空桶
                while (bq.top >= -bq.R && bq.bkt[bq.top + bq.R].empty())
                    --bq.top;
                if (bq.top < -bq.R)
                    return;
                auto &b = bq.bkt[bq.top + bq.R];
                for (int i = (int)b.size() - 1; i >= 0; --i)
                {
                    int v = b[i];
                    if (!locked[v] && g.pwgts[1 - from] + g.Vwgt(v) <= maxTo && gain(v) > cGain)
                    {
                        cand = v;
                        cFrom = from;
                        cGain = gain(v);
                        return;
                    }
                }
            };
            tryPick(bq0, 0, maxPW1);
            tryPick(bq1, 1, maxPW0);
            if (cand == -1)
                break;

            // 从队列删除 cand 并锁定（锁定后不再重入队列）
            if (cFrom == 0)
                bq0.erase(cand);
            else
                bq1.erase(cand);
            locked[cand] = true;

            int src = cFrom, dst = 1 - src, wv = g.Vwgt(cand);
            g.pwgts[src] -= wv;
            g.pwgts[dst] += wv;
            g.where[cand] = dst;
            std::swap(g.id[cand], g.ed[cand]);

            // 增量更新邻居（仅未锁定的）
            for (int ei = g.xadj[cand]; ei < g.xadj[cand + 1]; ++ei)
            {
                int u = g.adjncy[ei], ew = g.Ewgt(ei);
                if (locked[u])
                    continue;
                int old = gain(u);
                if (g.where[u] == dst)
                {
                    g.id[u] += ew;
                    g.ed[u] -= ew;
                }
                else
                {
                    g.id[u] -= ew;
                    g.ed[u] += ew;
                }
                int dg = gain(u) - old;
                BQ &qu = (g.where[u] == 0) ? bq0 : bq1;
                bool was = qu.inQ(u), need = (g.ed[u] > 0);
                if (was && need)
                    qu.update(u, dg);
                else if (was)
                    qu.erase(u);
                else if (need)
                    qu.push(u, gain(u));
            }
            // ← 不把 cand 加回队列（锁定）

            hist.push_back({cand, src});
            cum += cGain;
            if (cum > best)
            {
                best = cum;
                bestStep = (int)hist.size() - 1;
            }
        }

        // 回滚到历史最优
        for (int i = (int)hist.size() - 1; i > bestStep; --i)
        {
            int v = hist[i].v, fp = hist[i].from, wv = g.Vwgt(v);
            g.where[v] = fp;
            g.pwgts[fp] += wv;
            g.pwgts[1 - fp] -= wv;
        }
        ComputeParams(g);
        return best;
    }

    // ── 平衡修复 Pass ─────────────────────────────────────────────
    int BalanceFix(Graph &g, int t0, int t1, real_t ub)
    {
        if (Imbalance(g, t0, t1) <= ub)
            return 0;
        int mPW0 = (int)(ub * t0 + 0.5f), mPW1 = (int)(ub * t1 + 0.5f);
        int moved = 0, maxMoves = g.nvtxs; // 防止死循环上限
        while (Imbalance(g, t0, t1) > ub && moved < maxMoves)
        {
            bool h0 = (g.pwgts[0] * t1 > g.pwgts[1] * t0);
            int src = h0 ? 0 : 1, dst = 1 - src, maxDst = h0 ? mPW1 : mPW0;
            int bestV = -1, bestG = INT_MIN;
            for (int i = 0; i < g.nbnd; ++i)
            {
                int v = g.bndind[i];
                if (g.where[v] != src || g.pwgts[dst] + g.Vwgt(v) > maxDst * 2)
                    continue;
                int gv = g.ed[v] - g.id[v];
                if (gv > bestG)
                {
                    bestG = gv;
                    bestV = v;
                }
            }
            if (bestV == -1)
                break;
            int wv = g.Vwgt(bestV);
            g.pwgts[src] -= wv;
            g.pwgts[dst] += wv;
            g.where[bestV] = dst;
            for (int ei = g.xadj[bestV]; ei < g.xadj[bestV + 1]; ++ei)
            {
                int u = g.adjncy[ei], ew = g.Ewgt(ei);
                if (g.where[u] == dst)
                {
                    g.id[u] += ew;
                    g.ed[u] -= ew;
                }
                else
                {
                    g.id[u] -= ew;
                    g.ed[u] += ew;
                }
            }
            std::swap(g.id[bestV], g.ed[bestV]);
            ++moved;
        }
        if (moved > 0)
        {
            ComputeParams(g);
            FMPass(g, mPW0, mPW1);
        }
        return moved;
    }

    real_t Imbalance(const Graph &g, int t0, int t1)
    {
        if (t0 <= 0 || t1 <= 0)
            return 1.0;
        return std::max((real_t)g.pwgts[0] / t0, (real_t)g.pwgts[1] / t1);
    }

} // namespace bisect

// ─────────────────────────────────────────────────────────────
BisectResult Bisect(Graph &g, const BisectOptions &opts)
{
    using namespace bisect;
    if (g.nvtxs == 0)
        return {};
    if (g.tvwgt.empty())
        g.InitTvwgt();
    int W = g.tvwgt[0];
    int t0 = opts.target0 > 0 ? opts.target0 : W / 2;
    int t1 = std::max(W - t0, 1);
    t0 = W - t1;
    int mPW0 = (int)(opts.ubFactor * t0 + 0.5f), mPW1 = (int)(opts.ubFactor * t1 + 0.5f);

    BisectResult best;
    std::vector<int> bWhere, bPwgts;
    std::mt19937 rng(opts.seed);

    for (int tr = 0; tr < opts.nTrials; ++tr)
    {
        GrowInit(g, (int)rng(), t0, t1);
        ComputeParams(g);
        if (opts.verbose)
            std::printf("  [trial%2d] init cut=%d pw=[%d,%d]\n", tr, g.mincut, g.pwgts[0], g.pwgts[1]);
        for (int p = 0; p < opts.nFMPasses; ++p)
        {
            int gain = FMPass(g, mPW0, mPW1);
            if (opts.verbose)
                std::printf("    FM%d gain=%d cut=%d\n", p, gain, g.mincut);
            if (gain <= 0)
                break;
        }
        real_t imb = Imbalance(g, t0, t1);
        if (imb > opts.ubFactor)
        {
            BalanceFix(g, t0, t1, opts.ubFactor);
            imb = Imbalance(g, t0, t1);
        }
        bool bal = (imb <= opts.ubFactor);
        if (opts.verbose)
            std::printf("  [trial%2d] cut=%d imb=%.4f %s\n", tr, g.mincut, imb, bal ? "OK" : "UNBAL");
        bool better = (bal && !best.balanced) || (bal == best.balanced && g.mincut < best.mincut);
        if (better)
        {
            best = {g.mincut, imb, bal, tr};
            bWhere = g.where;
            bPwgts = g.pwgts;
        }
    }
    g.where = bWhere;
    g.pwgts = bPwgts;
    ComputeParams(g);
    if (opts.verbose)
        std::printf("  ==> best=%d cut=%d imb=%.4f %s\n\n", best.bestTrial, best.mincut, best.imbalance, best.balanced ? "BAL" : "UNBAL");
    return best;
}