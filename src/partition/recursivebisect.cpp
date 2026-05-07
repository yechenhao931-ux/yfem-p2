#include "partition/recursivebisect.hpp"
#include <algorithm>
#include <cstdio>
#include <numeric>
#include <random>
#include <vector>

static void ExtractSubgraph(const Graph &par, const std::vector<int> &vtx,
                            Graph &sub, std::vector<int> &l2o)
{
    const int ns = (int)vtx.size();
    l2o = vtx;
    std::vector<int> g2l(par.nvtxs, -1);
    for (int i = 0; i < ns; ++i)
        g2l[vtx[i]] = i;
    sub.nvtxs = ns;
    sub.ncon = par.ncon;
    sub.xadj.assign(ns + 1, 0);
    sub.adjncy.clear();
    sub.adjwgt.clear();
    sub.vwgt.resize(ns);
    sub.vsize.resize(ns);
    for (int i = 0; i < ns; ++i)
    {
        int oi = vtx[i];
        sub.vwgt[i] = par.Vwgt(oi);
        sub.vsize[i] = par.Vsize(oi);
        for (int ei = par.xadj[oi]; ei < par.xadj[oi + 1]; ++ei)
        {
            int lnb = g2l[par.adjncy[ei]];
            if (lnb < 0)
                continue;
            sub.adjncy.push_back(lnb);
            sub.adjwgt.push_back(par.Ewgt(ei));
            ++sub.xadj[i + 1];
        }
    }
    for (int i = 0; i < ns; ++i)
        sub.xadj[i + 1] += sub.xadj[i];
    sub.nedges = (int)sub.adjncy.size();
    sub.InitTvwgt();
}

static void RBRecurse(const Graph &orig, const std::vector<int> &vtx,
                      int partOff, int nparts, const RBOptions &opts,
                      std::vector<int> &whereOut, std::mt19937 &rng, int &totCut)
{
    if (nparts == 1)
    {
        for (int v : vtx)
            whereOut[v] = partOff;
        return;
    }
    if ((int)vtx.size() <= 1)
    {
        if (!vtx.empty())
            whereOut[vtx[0]] = partOff;
        return;
    }

    Graph sub;
    std::vector<int> l2o;
    ExtractSubgraph(orig, vtx, sub, l2o);

    int n0 = nparts / 2, n1 = nparts - n0;
    int W = sub.tvwgt[0];
    int t0 = std::max(1, (int)((real_t)n0 / nparts * W + 0.5f));
    int t1 = std::max(1, W - t0);

    if (opts.verbose)
        std::printf("[RB] nparts=%d vtxs=%d t=[%d,%d]\n", nparts, (int)vtx.size(), t0, t1);

    BisectOptions bo;
    bo.target0 = t0;
    bo.target1 = t1;
    bo.nTrials = opts.nTrials;
    bo.nFMPasses = opts.nFMPasses;
    bo.ubFactor = opts.ubFactor;
    bo.seed = (int)rng();
    bo.verbose = false;

    BisectResult br = Bisect(sub, bo);
    totCut += br.mincut;
    if (opts.verbose)
        std::printf("  cut=%d imb=%.4f %s\n", br.mincut, br.imbalance, br.balanced ? "OK" : "UNBAL");

    std::vector<int> p0, p1;
    for (int li = 0; li < sub.nvtxs; ++li)
    {
        if (sub.where[li] == 0)
            p0.push_back(l2o[li]);
        else
            p1.push_back(l2o[li]);
    }
    RBRecurse(orig, p0, partOff, n0, opts, whereOut, rng, totCut);
    RBRecurse(orig, p1, partOff + n0, n1, opts, whereOut, rng, totCut);
}

RBResult RecursiveBisect(Graph &graph, const RBOptions &opts)
{
    if (opts.nparts < 1 || graph.nvtxs == 0)
        return {};
    if (graph.tvwgt.empty())
        graph.InitTvwgt();

    graph.where.assign(graph.nvtxs, 0);
    int totCut = 0;
    std::mt19937 rng(opts.seed);
    std::vector<int> all(graph.nvtxs);
    std::iota(all.begin(), all.end(), 0);

    if (opts.verbose)
        std::printf("===== RecursiveBisect n=%d nparts=%d =====\n", graph.nvtxs, opts.nparts);

    RBRecurse(graph, all, 0, opts.nparts, opts, graph.where, rng, totCut);

    graph.pwgts.assign(opts.nparts, 0);
    for (int v = 0; v < graph.nvtxs; ++v)
        graph.pwgts[graph.where[v]] += graph.Vwgt(v);

    const int K = opts.nparts;
    KwayFMOpts ko;
    ko.nparts = K;
    ko.nIter = opts.postKwayIter;
    ko.ubFactor = opts.ubFactor;
    ko.verbose = opts.verbose;

    if (opts.postRefine == PostRefineMode::EdgeCut || opts.postRefine == PostRefineMode::Both)
    {
        ComputeCkrinfo(graph, K);
        int g = KwayFMCut(graph, ko);
        if (opts.verbose)
            std::printf("[PostRefine Cut] gain=%d cut=%d\n", g, graph.mincut);
    }
    if (opts.postRefine == PostRefineMode::Volume || opts.postRefine == PostRefineMode::Both)
    {
        if (graph.ckrinfo.empty())
            ComputeCkrinfo(graph, K);
        ComputeVkrinfo(graph, K);
        int g = KwayFMVol(graph, ko);
        if (opts.verbose)
            std::printf("[PostRefine Vol] gain=%d vol=%d\n", g, graph.minvol);
    }
    if (opts.postRefine == PostRefineMode::None)
        ComputeCkrinfo(graph, K);

    // 重算 mincut
    {
        int c = 0;
        for (int v = 0; v < graph.nvtxs; ++v)
            for (int ei = graph.xadj[v]; ei < graph.xadj[v + 1]; ++ei)
                if (graph.where[v] != graph.where[graph.adjncy[ei]])
                    c += graph.Ewgt(ei);
        graph.mincut = c / 2;
    }

    real_t ideal = (real_t)graph.tvwgt[0] / K;
    real_t maxImb = 0, sumImb = 0;
    for (int p = 0; p < K; ++p)
    {
        real_t r = (real_t)graph.pwgts[p] / ideal;
        real_t d = r - 1.f;
        maxImb = std::max(maxImb, std::abs(d));
        sumImb += std::abs(d);
    }

    RBResult res;
    res.totalCut = graph.mincut;
    res.totalVol = graph.minvol;
    res.maxImbalance = maxImb;
    res.avgImbalance = sumImb / K;
    res.balanced = (maxImb + 1.f <= opts.ubFactor);
    res.partWeights = graph.pwgts;

    if (opts.verbose)
    {
        std::printf("\n===== 完成 =====\n总切边数   : %d\n最大不平衡 : %.2f%%\n平均不平衡 : %.2f%%\n满足约束   : %s\n理想权重   : %.1f\n",
                    res.totalCut, res.maxImbalance * 100, res.avgImbalance * 100, res.balanced ? "是" : "否", ideal);
        for (int p = 0; p < K; ++p)
        {
            real_t r = (real_t)graph.pwgts[p] / ideal;
            std::printf("  [%2d] wgt=%-5d ratio=%.3f%s\n", p, graph.pwgts[p], r, (r > opts.ubFactor) ? " ←超标" : "");
        }
    }
    return res;
}