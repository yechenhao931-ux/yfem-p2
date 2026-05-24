#include "partition/kwayrefine.hpp"
#include <algorithm>
#include <climits>
#include <cstdio>
#include <numeric>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

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
// 并行化说明（OpenMP）：
//   邻边扫描按顶点天然独立，是 O(V+E) 的热点（每轮 FM、每层反粗化都会调用）。
//   为保证结果与串行版本逐字节一致：
//     · 每线程使用独立的 tmpEd/seen 临时数组，避免共享写竞争；
//     · cnbrPool 的偏移用前缀和按 v 递增顺序预先算好（= 串行追加顺序）；
//     · 每个顶点内邻接分区对仍按"首次出现"顺序写入（= 串行 seen 顺序）；
//     · pwgts / mincut / 边界表在 v 递增的串行小循环中累加（O(V)，顺序与串行同）。
void ComputeCkrinfo(Graph &g, int nparts)
{
    const int n = g.nvtxs;
    g.AllocateCkrinfo(std::min(nparts, 16));
    g.pwgts.assign(nparts, 0);
    g.bndptr.assign(n, -1);
    g.bndind.clear();
    g.nbnd = 0;
    g.mincut = 0;

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    std::vector<std::vector<int>> tmpEdTL(nthreads, std::vector<int>(nparts, 0));
    std::vector<std::vector<int>> seenTL(nthreads);

    // ── Pass A（并行）：每个顶点的 id / ed / nnbrs（工作量均匀 → static）──
    // if(n>阈值)：小图直接串行，避免线程组开销（box.mesh 等小网格不退化）
#pragma omp parallel for schedule(static) if (n > 16384)
    for (int v = 0; v < n; ++v)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        std::vector<int> &tmpEd = tmpEdTL[tid];
        std::vector<int> &seen = seenTL[tid];
        seen.clear();
        int pv = g.where[v];
        int id = 0, ed = 0;
        for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
        {
            int u = g.adjncy[ei], ew = g.Ewgt(ei), pu = g.where[u];
            if (pu == pv)
                id += ew;
            else
            {
                ed += ew;
                if (tmpEd[pu] == 0)
                    seen.push_back(pu);
                tmpEd[pu] += ew;
            }
        }
        Ckrinfo &ci = g.ckrinfo[v];
        ci.id = id;
        ci.ed = ed;
        ci.nnbrs = (int)seen.size();
        for (int pu : seen)
            tmpEd[pu] = 0; // 复位本线程临时数组
    }

    // ── 前缀和（串行）：inbr 偏移 == 串行追加顺序 ──
    size_t total = 0;
    for (int v = 0; v < n; ++v)
    {
        g.ckrinfo[v].inbr = (int)total;
        total += (size_t)g.ckrinfo[v].nnbrs;
    }
    g.cnbrPool.assign(total, Cnbr{});

    // ── Pass B（并行）：仅边界顶点(nnbrs>0)需写 cnbrPool；
    //    内部顶点直接跳过，省去重复的邻边扫描（FE 网格内部点占绝大多数）。
    //    工作量不均匀 → dynamic 负载均衡。──
#pragma omp parallel for schedule(dynamic, 256) if (n > 16384)
    for (int v = 0; v < n; ++v)
    {
        if (g.ckrinfo[v].nnbrs == 0)
            continue;
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        std::vector<int> &tmpEd = tmpEdTL[tid];
        std::vector<int> &seen = seenTL[tid];
        seen.clear();
        int pv = g.where[v];
        for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
        {
            int u = g.adjncy[ei], ew = g.Ewgt(ei), pu = g.where[u];
            if (pu != pv)
            {
                if (tmpEd[pu] == 0)
                    seen.push_back(pu);
                tmpEd[pu] += ew;
            }
        }
        int base = g.ckrinfo[v].inbr;
        for (size_t i = 0; i < seen.size(); ++i)
        {
            int pu = seen[i];
            g.cnbrPool[base + i] = {pu, tmpEd[pu]};
            tmpEd[pu] = 0;
        }
    }

    // ── pwgts / mincut / 边界（串行 v 递增，结果与串行一致）──
    for (int v = 0; v < n; ++v)
        g.pwgts[g.where[v]] += g.Vwgt(v);
    for (int v = 0; v < n; ++v)
    {
        int ed = g.ckrinfo[v].ed;
        if (ed > 0)
        {
            g.mincut += ed;
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
        g.minvol += g.Vsize(v);
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

// ================================================================
//  IndepSetRefineCut  ──  A-3: 独立集并行精化（G-kway, DAC'24）
//
//  与串行 FM 的区别：FM 用优先队列每次只移动一个顶点；本算法每轮并行地
//  找到一个"独立集"的正增益移动（互不相邻），一次性应用，迭代度更高。
//  每轮三步（对应论文 Algorithm 2/3）：
//    1) 找移动：并行计算每个边界点的最佳合法目标(gain>0 且不撑爆目标分区)；
//       仅当相邻"想动"的顶点里本点 ID 最小时才入选 → 保证独立集（移动相加性）
//    2) 排序：按 gain 降序（平局按顶点 ID 升序，确定性）
//    3) 选择：取"最长不破坏平衡前缀"并行应用（论文的 scan + 最长平衡子序列）
//  迭代到无正增益移动或无法在平衡下推进为止。
//  注：移动集为独立集 → 总增益 = 各移动增益之和（边不被两端同时移动）。
// ================================================================
int IndepSetRefineCut(Graph &g, const KwayFMOpts &opts)
{
    const int n = g.nvtxs, K = opts.nparts;

    int W = 0;
    for (int v = 0; v < n; ++v)
        W += g.Vwgt(v);
    real_t ideal = (real_t)W / K;
    std::vector<int> maxPW(K);
    for (int p = 0; p < K; ++p)
        maxPW[p] = (int)(opts.ubFactor * ideal + 0.5f);

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif

    struct Move
    {
        int v, src, dst, gain;
    };
    int totalGain = 0;

    for (int iter = 0; iter < opts.nIter; ++iter)
    {
        ComputeCkrinfo(g, K); // 刷新 ckrinfo / cnbrPool / pwgts / bndind

        std::vector<int> bestDst(n, -1), bestGain(n, 0);

        // ── 步骤1a：每个边界点的最佳合法目标（并行）──
#pragma omp parallel for schedule(dynamic, 256) if (g.nbnd > 4096)
        for (int bi = 0; bi < g.nbnd; ++bi)
        {
            int v = g.bndind[bi];
            const Ckrinfo &ci = g.ckrinfo[v];
            if (ci.ed == 0 || ci.inbr < 0)
                continue;
            int wv = g.Vwgt(v), bd = -1, bg = 0;
            for (int i = 0; i < ci.nnbrs; ++i)
            {
                int p = g.cnbrPool[ci.inbr + i].pid;
                int gain = g.cnbrPool[ci.inbr + i].ed - ci.id;
                if (gain > 0 && g.pwgts[p] + wv <= maxPW[p] &&
                    (gain > bg || (gain == bg && (bd == -1 || p < bd))))
                {
                    bg = gain;
                    bd = p;
                }
            }
            bestDst[v] = bd;
            bestGain[v] = bg;
        }

        // ── 步骤1b：独立集选择（并行）──
        //   v 入选 ⇔ 有合法目标 且 相邻"想动"顶点中 v 的 ID 最小
        std::vector<std::vector<Move>> moveTL(nthreads);
#pragma omp parallel for schedule(dynamic, 256) if (g.nbnd > 4096)
        for (int bi = 0; bi < g.nbnd; ++bi)
        {
            int v = g.bndind[bi];
            if (bestDst[v] < 0)
                continue;
            bool selected = true;
            for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
            {
                int u = g.adjncy[ei];
                if (bestDst[u] >= 0 && u < v)
                {
                    selected = false;
                    break;
                }
            }
            if (!selected)
                continue;
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            moveTL[tid].push_back({v, g.where[v], bestDst[v], bestGain[v]});
        }

        std::vector<Move> moves;
        for (auto &mv : moveTL)
            moves.insert(moves.end(), mv.begin(), mv.end());
        if (moves.empty())
            break;

        // ── 步骤2：按增益降序（平局按 ID 升序）──
        std::sort(moves.begin(), moves.end(), [](const Move &a, const Move &b)
                  { return a.gain != b.gain ? a.gain > b.gain : a.v < b.v; });

        // ── 步骤3：最长"不破坏平衡"前缀 ──
        //   limit[p] = max(maxPW[p], startPW[p])：欠载分区可填到 maxPW；
        //   已超载分区(初始不平衡)不允许再增长。nViol 跟踪越界分区数。
        std::vector<int> cumPW = g.pwgts;
        std::vector<int> lim(K);
        for (int p = 0; p < K; ++p)
            lim[p] = std::max(maxPW[p], g.pwgts[p]);
        auto over = [&](int p)
        { return cumPW[p] > lim[p]; };

        int nViol = 0, bestLen = 0, prefixGain = 0, runGain = 0;
        for (size_t j = 0; j < moves.size(); ++j)
        {
            const Move &m = moves[j];
            int wv = g.Vwgt(m.v);
            bool sB = over(m.src);
            cumPW[m.src] -= wv;
            nViol += (int)over(m.src) - (int)sB; // src 减少：可能 +越界→不越界
            bool dB = over(m.dst);
            cumPW[m.dst] += wv;
            nViol += (int)over(m.dst) - (int)dB; // dst 增加：可能 不越界→越界
            runGain += m.gain;
            if (nViol == 0)
            {
                bestLen = (int)j + 1;
                prefixGain = runGain;
            }
        }

        if (bestLen == 0)
            break; // 平衡约束下无法推进

        // 应用前缀（移动互不相邻 → 可并行；这里直接顺序应用，O(bestLen)）
        for (int j = 0; j < bestLen; ++j)
        {
            const Move &m = moves[j];
            int wv = g.Vwgt(m.v);
            g.where[m.v] = m.dst;
            g.pwgts[m.src] -= wv;
            g.pwgts[m.dst] += wv;
        }
        totalGain += prefixGain;

        if (opts.verbose)
            std::printf("  [IndepSet iter%d] moves=%zu applied=%d gain=%d\n",
                        iter, moves.size(), bestLen, prefixGain);
    }

    ComputeCkrinfo(g, K); // 刷新 mincut / ckrinfo 供下游使用
    return totalGain;
}
