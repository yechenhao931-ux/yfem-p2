#include <partition/kmeanspart.hpp>
#include <partition/connectrepair.hpp>
#include <common/macro.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <vector>

// ── 欧氏距离平方 ─────────────────────────────────────────────
static real_t dist2(const Coord &a, const Coord &b)
{
    real_t dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

static real_t dist2C(const Coord &a, const std::vector<real_t> &c)
{
    real_t dx = a.x - c[0], dy = a.y - c[1], dz = a.z - c[2];
    return dx * dx + dy * dy + dz * dz;
}

// ================================================================
//  K-means++ 种子初始化
//  第一个种子随机选，后续每个种子以 dist²(v, 最近已选中心)
//  为概率权重抽样，使初始中心尽量分散
// ================================================================

static std::vector<std::vector<real_t>> KmeansPPInit(const Graph &g, int K, std::mt19937 &rng)
{
    const int n = g.nvtxs;
    std::vector<std::vector<real_t>> centers;
    centers.reserve(K);

    // 第一个中心：随机选一个顶点
    std::uniform_int_distribution<int> uni(0, n - 1);
    int first = uni(rng);
    centers.push_back({g.coordinates[first].x, g.coordinates[first].y, g.coordinates[first].z});

    std::vector<real_t> minDist(n, std::numeric_limits<real_t>::max());

    for (int k = 1; k < K; k++)
    {
        // 更新每个顶点到最近已选中心的距离平方
        const auto &last = centers.back();
        real_t total = 0.;
        for (int v = 0; v < n; v++)
        {
            real_t d = dist2C(g.coordinates[v], last);
            if (d < minDist[v])
            {
                minDist[v] = d;
            }
            total += minDist[v];
        }

        // 按 dist² 比例抽样（轮盘赌）
        std::uniform_real_distribution<real_t> draw(0., total);
        real_t r = draw(rng);
        real_t cum = 0.;
        int chosen = n - 1;
        for (int v = 0; v < n; ++v)
        {
            cum += minDist[v];
            if (cum >= r)
            {
                chosen = v;
                break;
            }
        }
        centers.push_back({g.coordinates[chosen].x,
                           g.coordinates[chosen].y,
                           g.coordinates[chosen].z});
    }

    return centers;
}

// ================================================================
//  ComputeCentroids
//  加权重心：centroid[p] = Σ(vwgt[v] * coord[v]) / Σ vwgt[v]
// ================================================================

std::vector<std::vector<real_t>> ComputeCentroids(const Graph &g, int nparts)
{
    std::vector<std::vector<real_t>> c(nparts, std::vector<real_t>(3, 0.));

    std::vector<real_t> wsum(nparts, 0.);

    for (int v = 0; v < g.nvtxs; v++)
    {
        int p = g.where[v];
        real_t w = (real_t)g.Vwgt(v);
        c[p][0] += w * g.coordinates[v].x;
        c[p][1] += w * g.coordinates[v].y;
        c[p][2] += w * g.coordinates[v].z;
        wsum[p] += w;
    }
    for (int p = 0; p < nparts; ++p)
    {
        if (wsum[p] > 0.f)
        {
            c[p][0] /= wsum[p];
            c[p][1] /= wsum[p];
            c[p][2] /= wsum[p];
        }
    }

    return c;
}
// ================================================================
//  BalanceFix  ──  双向平衡修复（可推可拉）
//
//  双向策略：
//    · 重→轻：重分区的边界顶点移到轻邻居分区
//    · 轻→轻：若轻分区的邻居有更重的，优先从最重邻居处拉取
//  连通性保证：只移动在源分区有其他图邻居的顶点
//             （叶节点允许移动，因为移走后源分区仍有其他顶点）
// ================================================================
static void BalanceFix(Graph& g,int K,real_t ubFactor){
    real_t ideal=(real_t)g.tvwgt[0]/K;
    int maxPW=(int)(ubFactor*ideal+0.5);
    int minPW=(int)((2.f-ubFactor)*ideal-0.5); // 允许的最小权重
 
    // 预计算 id（同分区边权和）
    std::vector<int> vId(g.nvtxs,0);
    for(int v=0;v<g.nvtxs;++v){
        int pv=g.where[v];
        for(int ei=g.xadj[v];ei<g.xadj[v+1];++ei)
            if(g.where[g.adjncy[ei]]==pv)vId[v]+=g.Ewgt(ei);
    }
 
    bool changed=true;
    int guard=g.nvtxs;
    while(changed&&--guard>0){
        changed=false;
 
        // 找最重和最轻分区
        int heavy=-1,light=-1;
        for(int p=0;p<K;++p){
            if(g.pwgts[p]>maxPW&&(heavy==-1||g.pwgts[p]>g.pwgts[heavy]))heavy=p;
            if(g.pwgts[p]<minPW&&(light==-1||g.pwgts[p]<g.pwgts[light]))light=p;
        }
        if(heavy==-1&&light==-1)break;
 
        // 优先处理最重分区：找其边界顶点中最优的可移动点
        if(heavy!=-1){
            int bestV=-1,bestDst=-1;real_t bestScore=-1e30;
            for(int v=0;v<g.nvtxs;++v){
                if(g.where[v]!=heavy)continue;
                for(int ei=g.xadj[v];ei<g.xadj[v+1];++ei){
                    int u=g.adjncy[ei],pu=g.where[u];
                    if(pu==heavy)continue;
                    if(g.pwgts[pu]+g.Vwgt(v)>maxPW)continue;
                    // 增益：到目标分区的边减少切边；离开源分区增加切边
                    int edToDst=0;
                    for(int ei2=g.xadj[v];ei2<g.xadj[v+1];++ei2)
                        if(g.where[g.adjncy[ei2]]==pu)edToDst+=g.Ewgt(ei2);
                    real_t gain=(real_t)(edToDst-vId[v]);
                    // 奖励目标分区越轻越好
                    real_t score=gain+(ideal-(real_t)g.pwgts[pu])/(ideal+1.f)*0.5f;
                    if(score>bestScore){bestScore=score;bestV=v;bestDst=pu;}
                }
            }
            if(bestV!=-1){
                int wv=g.Vwgt(bestV);
                // 更新 vId：v 移走后，其邻居的 id/ed 发生变化
                for(int ei=g.xadj[bestV];ei<g.xadj[bestV+1];++ei){
                    int nb=g.adjncy[ei],ew=g.Ewgt(ei),pnb=g.where[nb];
                    if(pnb==heavy)vId[nb]-=ew; // nb 失去同分区邻居 v
                    if(pnb==bestDst)vId[nb]+=ew;// nb 获得新同分区邻居 v
                }
                vId[bestV]=0;// 重算 v 的 id
                for(int ei=g.xadj[bestV];ei<g.xadj[bestV+1];++ei)
                    if(g.where[g.adjncy[ei]]==bestDst)vId[bestV]+=g.Ewgt(ei);
                g.pwgts[heavy]-=wv;g.pwgts[bestDst]+=wv;g.where[bestV]=bestDst;
                changed=true;continue;
            }
        }
 
        // 处理最轻分区：从最重邻居处拉取一个边界顶点
        if(light!=-1){
            int bestV=-1,bestSrc=-1;real_t bestScore=-1e30f;
            for(int v=0;v<g.nvtxs;++v){
                int pv=g.where[v];
                if(pv==light)continue;
                if(g.pwgts[pv]-g.Vwgt(v)<minPW)continue; // 源分区不能变太轻
                // v 是否与 light 分区相邻
                bool adjLight=false;
                for(int ei=g.xadj[v];ei<g.xadj[v+1];++ei)
                    if(g.where[g.adjncy[ei]]==light){adjLight=true;break;}
                if(!adjLight)continue;
                if(g.pwgts[light]+g.Vwgt(v)>maxPW)continue;
                real_t score=(real_t)g.pwgts[pv]; // 越重越优先被拉取
                if(score>bestScore){bestScore=score;bestV=v;bestSrc=pv;}
            }
            if(bestV!=-1){
                int wv=g.Vwgt(bestV);
                for(int ei=g.xadj[bestV];ei<g.xadj[bestV+1];++ei){
                    int nb=g.adjncy[ei],ew=g.Ewgt(ei),pnb=g.where[nb];
                    if(pnb==bestSrc)vId[nb]-=ew;
                    if(pnb==light)vId[nb]+=ew;
                }
                vId[bestV]=0;
                for(int ei=g.xadj[bestV];ei<g.xadj[bestV+1];++ei)
                    if(g.where[g.adjncy[ei]]==light)vId[bestV]+=g.Ewgt(ei);
                g.pwgts[bestSrc]-=wv;g.pwgts[light]+=wv;g.where[bestV]=light;
                changed=true;
            }
        }
    }
}
 
static int EvalCut(const Graph& g){
    int cut=0;
    for(int v=0;v<g.nvtxs;++v)
        for(int ei=g.xadj[v];ei<g.xadj[v+1];++ei)
            if(g.where[v]!=g.where[g.adjncy[ei]])cut+=g.Ewgt(ei);
    return cut/2;
}
// ================================================================
//  BalanceRepair
//  若某分区权重超过 ubFactor * ideal，将其边界顶点
//  （即在 graph 拓扑中有跨分区邻居的顶点）迁移到最近的轻分区，
//  直到不平衡消除。
//  注意：此处"边界"用简单定义——有任意邻居在不同分区。
// ================================================================
static void BalanceRepair(Graph &g, int nparts, std::vector<std::vector<real_t>> &centers, real_t ubFactor)
{
    real_t ideal = (real_t)g.tvwgt[0] / nparts;
    std::vector<int> pwgts(nparts, 0);
    for (int v = 0; v < g.nvtxs; v++)
    {
        pwgts[g.where[v]] += g.Vwgt(v);
    }

    bool changed = true;
    int guard = g.nvtxs; // 防止死循环

    while (changed && --guard > 0)
    {
        changed = false;
        for (int v = 0; v < g.nvtxs; v++)
        {
            int src = g.where[v];
            if ((real_t)pwgts[src] <= ubFactor * ideal)
            {
                continue;
            }
            // 检查是否是边界点（有邻居在不同分区）
            bool isBnd = false;
            for (int ei = g.xadj[v]; ei < g.xadj[v + 1]; ++ei)
            {
                if (g.where[g.adjncy[ei]] != src)
                {
                    isBnd = true;
                    break;
                }
            }

            if (!isBnd)
            {
                continue;
            }

            // 找最近的、接收后不超重的分区
            real_t bestD = std::numeric_limits<real_t>::max();
            int bestP = -1;
            int wv = g.Vwgt(v);

            for (int p = 0; p < nparts; p++)
            {
                if (p == src)
                    continue;
                if ((real_t)(pwgts[p] + wv) > ubFactor * ideal * 1.05)
                {
                    continue;
                }
                real_t d = dist2C(g.coordinates[v], centers[p]);

                if (d < bestD)
                {
                    bestD = d;
                    bestP = p;
                }
            }
            if (bestP == -1)
                continue;

            // 迁移
            pwgts[src] -= wv;
            pwgts[bestP] += wv;
            g.where[v] = bestP;
            changed = true;

            // 更新重心（增量更新）

            real_t ws = (real_t)pwgts[src] + wv;   // 迁移前 src 的权重
            real_t wd = (real_t)pwgts[bestP] - wv; // 迁移前 dst 的权重

            if (ws > 0)
            {
                centers[src][0] = (centers[src][0] * ws - (real_t)wv * g.coordinates[v].x) / (ws - wv);
                centers[src][1] = (centers[src][1] * ws - (real_t)wv * g.coordinates[v].y) / (ws - wv);
                centers[src][2] = (centers[src][2] * ws - (real_t)wv * g.coordinates[v].z) / (ws - wv);
            }
            centers[bestP][0] = (centers[bestP][0] * wd + (real_t)wv * g.coordinates[v].x) / (wd + wv);
            centers[bestP][1] = (centers[bestP][1] * wd + (real_t)wv * g.coordinates[v].y) / (wd + wv);
            centers[bestP][2] = (centers[bestP][2] * wd + (real_t)wv * g.coordinates[v].z) / (wd + wv);
        }
    }

    g.pwgts = pwgts;
}

// ================================================================
//  BFSSeedExpansion
//  从 K-means++ 几何种子出发，沿图邻接做 K-轮交替 BFS，
//  得到一个"图上连通"的初始划分。
//
//  关键差异：相比纯几何分配，本方法保证每个分区在图上
//  连通（最多有 0..1 个未触达顶点回退到几何最近）。
//  对有限元/CFD 网格，这极大降低了初始切边数和后续 FM 迭代次数。
// ================================================================
static void BFSSeedExpansion(Graph &g, int K,
                             const std::vector<std::vector<real_t>> &centers,
                             std::vector<int> &where)
{
    const int n = g.nvtxs;
    where.assign(n, -1);

    // 找最接近每个 centroid 的图顶点作为 BFS 源
    std::vector<int> sources(K, -1);
    std::vector<real_t> bestD(K, std::numeric_limits<real_t>::max());
    for (int v = 0; v < n; ++v)
    {
        for (int k = 0; k < K; ++k)
        {
            real_t d = dist2C(g.coordinates[v], centers[k]);
            if (d < bestD[k])
            {
                bestD[k] = d;
                sources[k] = v;
            }
        }
    }

    // 多源 BFS，按权重轮转扩展（保持平衡）
    std::vector<std::queue<int>> q(K);
    std::vector<int> partW(K, 0);
    int totW = 0;
    for (int v = 0; v < n; ++v) totW += g.Vwgt(v);
    real_t target = (real_t)totW / K;

    for (int k = 0; k < K; ++k)
    {
        int s = sources[k];
        if (s < 0 || where[s] != -1) {
            // 种子冲突或无效：用未分配的最近顶点
            real_t bd = std::numeric_limits<real_t>::max();
            int bv = -1;
            for (int v = 0; v < n; ++v)
            {
                if (where[v] != -1) continue;
                real_t d = dist2C(g.coordinates[v], centers[k]);
                if (d < bd) { bd = d; bv = v; }
            }
            if (bv < 0) continue;
            s = bv;
        }
        where[s] = k;
        partW[k] += g.Vwgt(s);
        q[k].push(s);
    }

    // 轮转 BFS：当前权重最轻的分区优先扩展，趋向平衡
    bool any = true;
    while (any)
    {
        any = false;
        // 选最轻且队列非空的分区
        int pick = -1;
        for (int k = 0; k < K; ++k)
        {
            if (q[k].empty()) continue;
            if (partW[k] >= 1.05f * target) continue; // 不让单分区严重超重
            if (pick == -1 || partW[k] < partW[pick]) pick = k;
        }
        if (pick == -1)
        {
            // 所有"未超标"队列空了：兜底用任何非空队列
            for (int k = 0; k < K; ++k)
                if (!q[k].empty()) { pick = k; break; }
            if (pick == -1) break;
        }
        any = true;
        int u = q[pick].front();
        q[pick].pop();
        for (int ei = g.xadj[u]; ei < g.xadj[u + 1]; ++ei)
        {
            int v = g.adjncy[ei];
            if (where[v] == -1)
            {
                where[v] = pick;
                partW[pick] += g.Vwgt(v);
                q[pick].push(v);
            }
        }
    }

    // 仍未触达的顶点（如孤立点 / 不连通分量）→ 选最近 centroid
    for (int v = 0; v < n; ++v)
    {
        if (where[v] == -1)
        {
            real_t bd = std::numeric_limits<real_t>::max();
            int bp = 0;
            for (int k = 0; k < K; ++k)
            {
                real_t d = dist2C(g.coordinates[v], centers[k]);
                if (d < bd) { bd = d; bp = k; }
            }
            where[v] = bp;
        }
    }
}

// ================================================================
//  KmeansPartition  ──  主入口
// ================================================================s
KmeansResult KmeansPartition(Graph &graph, int nparts, const KmeansOpts &opts)
{
    YFEM_ASSERT((int)graph.coordinates.size() == graph.nvtxs, "KmeansPartition: coordinates 未填充");
    if (graph.tvwgt.empty())
    {
        graph.InitTvwgt();
    }

    const int n = graph.nvtxs;
    std::mt19937 rng(opts.seed);

    // ── K-means++ 初始化 ──────────────────────────────────────
    auto centers = KmeansPPInit(graph, nparts, rng);

    graph.where.assign(n, 0);

    // ── BFS 种子扩展（拓扑感知初始化） ────────────────────────
    if (opts.bfsSeedExpansion)
    {
        BFSSeedExpansion(graph, nparts, centers, graph.where);
        // 用 BFS 结果重新计算 centroid，使 Lloyd 从拓扑合理的位置出发
        centers = ComputeCentroids(graph, nparts);
    }

    std::vector<real_t> minD(n);

    KmeansResult res;

    // 中心间距离² 矩阵（用于 Elkan/Phillips 三角不等式剪枝）
    //   引理：若 d²(v, c_curBest) ≤ 0.25 * d²(c_curBest, c_q)，则 q 不可能更近，
    //         可跳过 d²(v, c_q) 计算。dist² 单调，无需开方。
    //   收益：对 K 较大（≥8）时显著，相同时间内可跑更多 K-means trial，
    //         从而提升初始划分质量。
    std::vector<std::vector<real_t>> dC2(nparts, std::vector<real_t>(nparts, 0.));
    auto recomputeDC2 = [&]()
    {
        for (int i = 0; i < nparts; ++i)
        {
            dC2[i][i] = 0.0;
            for (int j = i + 1; j < nparts; ++j)
            {
                real_t dx = centers[i][0] - centers[j][0];
                real_t dy = centers[i][1] - centers[j][1];
                real_t dz = centers[i][2] - centers[j][2];
                real_t d2 = dx * dx + dy * dy + dz * dz;
                dC2[i][j] = dC2[j][i] = d2;
            }
        }
    };
    recomputeDC2();

    // ── Lloyd 迭代 ────────────────────────────────────────────
    for (int iter = 0; iter < opts.maxIter; iter++)
    {
        bool anyChange = false;
        real_t inertia = 0.0;

        // 分配阶段：每个顶点分配到最近中心（带三角不等式剪枝）
        for (int v = 0; v < n; v++)
        {
            // 从当前分配开始，最大化剪枝命中率
            int bestP = graph.where[v];
            if (bestP < 0 || bestP >= nparts) bestP = 0;
            real_t bestD = dist2C(graph.coordinates[v], centers[bestP]);
            for (int p = 0; p < nparts; p++)
            {
                if (p == bestP) continue;
                // 剪枝：4*d²(v,c_bestP) ≤ d²(c_bestP, c_p) ⇒ c_p 不可能更近
                if (4.0 * bestD <= dC2[bestP][p]) continue;
                real_t d = dist2C(graph.coordinates[v], centers[p]);
                if (d < bestD)
                {
                    bestD = d;
                    bestP = p;
                }
            }
            minD[v] = bestD;
            inertia += bestD;

            if (graph.where[v] != bestP)
            {
                graph.where[v] = bestP;
                anyChange = true;
            }
        }

        res.inertia = inertia;
        res.iters = iter + 1;

        // 更新阶段：重新计算加权重心
        std::vector<std::vector<real_t>> newC(nparts, std::vector<real_t>(3, 0.0));
        std::vector<real_t> wsum(nparts, 0.0);

        for (int v = 0; v < n; v++)
        {
            int p = graph.where[v];
            real_t w = (real_t)graph.Vwgt(v);
            newC[p][0] += w * graph.coordinates[v].x;
            newC[p][1] += w * graph.coordinates[v].y;
            newC[p][2] += w * graph.coordinates[v].z;
            wsum[p] += w;
        }

        real_t maxShift = 0.0;
        for (int p = 0; p < nparts; p++)
        {
            if (wsum[p] > 0.0)
            {
                newC[p][0] /= wsum[p];
                newC[p][1] /= wsum[p];
                newC[p][2] /= wsum[p];
            }
            real_t shift = std::sqrt(
                (newC[p][0] - centers[p][0]) * (newC[p][0] - centers[p][0]) +
                (newC[p][1] - centers[p][1]) * (newC[p][1] - centers[p][1]) +
                (newC[p][2] - centers[p][2]) * (newC[p][2] - centers[p][2]));
            maxShift = std::max(maxShift, shift);
        }
        centers = newC;
        // 中心移动后必须重算 dC2 才能继续安全剪枝
        recomputeDC2();

        if (opts.verbose)
            std::printf("  [Kmeans iter%3d] inertia=%.4f  maxShift=%.6f\n",
                        iter, inertia, maxShift);

        if (!anyChange || maxShift < opts.tol)
        {
            res.converged = true;
            break;
        }
    }

    // ── 负载均衡修复 ──────────────────────────────────────────
    graph.pwgts.assign(nparts, 0);
    for (int v = 0; v < n; v++)
    {
        graph.pwgts[graph.where[v]] += graph.Vwgt(v);
    }
    if (opts.useTopoBalanceFix)
    {
        // 拓扑感知：仅在边界点上换 src/dst，
        // 配合 BFS 初始化保留分区连通性
        BalanceFix(graph, nparts, opts.ubFactor);
    }
    else
    {
        // 旧路径：纯几何 BalanceRepair（可能破坏连通性）
        BalanceRepair(graph, nparts, centers, opts.ubFactor);
    }

    // ── 连通性修复（K-means 几何分配可能产生孤岛） ──────────
    if (opts.enforceConnectivity)
    {
        EnforceConnectivity(graph, nparts, opts.ubFactor, 0.95f);
        // 重新计算 centroid 以反映迁移
        centers = ComputeCentroids(graph, nparts);
    }

    res.centroids = centers;

    // 统计最终不平衡
    real_t ideal = (real_t)graph.tvwgt[0] / nparts;
    real_t maxImb = 0.;

    for (int p = 0; p < nparts; p++)
    {
        maxImb = std::max(maxImb, std::abs((real_t)graph.pwgts[p] / ideal - 1.));
    }
    res.maxImbalance = maxImb;

    if (opts.verbose)
        std::printf("  [Kmeans] converged=%d iters=%d inertia=%.4f maxImb=%.2f%%\n",
                    res.converged, res.iters, res.inertia, res.maxImbalance * 100);

    return res;
}

// OPT1
static real_t dist2v(const Coord &a, const Coord &b)
{
    real_t dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}
static real_t dist2c(const Coord &a, const std::vector<real_t> &c)
{
    real_t dx = a.x - c[0], dy = a.y - c[1], dz = a.z - c[2];
    return dx * dx + dy * dy + dz * dz;
}

// K-means++ 种子顶点
static std::vector<int> SelectSeeds(const Graph &g, int K, std::mt19937 &rng)
{
    const int n = g.nvtxs;
    std::vector<int> seeds;
    seeds.reserve(K);
    std::vector<real_t> minD(n, std::numeric_limits<real_t>::max());
    std::uniform_int_distribution<int> uni(0, n - 1);
    seeds.push_back(uni(rng));
    for (int k = 1; k < K; ++k)
    {
        int last = seeds.back();
        real_t total = 0.0;
        for (int v = 0; v < n; ++v)
        {
            real_t d = dist2v(g.coordinates[v], g.coordinates[last]);
            if (d < minD[v])
                minD[v] = d;
            total += minD[v];
        }
        if (total < 1e-12f)
        {
            seeds.push_back(uni(rng));
            continue;
        }
        std::uniform_real_distribution<real_t> draw(0.0, total);
        real_t r = draw(rng), cum = 0.0;
        int chosen = n - 1;
        for (int v = 0; v < n; ++v)
        {
            cum += minD[v];
            if (cum >= r)
            {
                chosen = v;
                break;
            }
        }
        seeds.push_back(chosen);
    }
    return seeds;
}
