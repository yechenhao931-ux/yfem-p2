#ifndef BISECTFM_HPP
#define BISECTFM_HPP
#include "../graph/graph.hpp"
#include <limits>

struct BisectOptions
{
    int target0 = 0, target1 = 0;
    int nTrials = 10, nFMPasses = 8;
    real_t ubFactor = 1.03;
    int seed = 42;
    bool verbose = false;
};

struct BisectResult
{
    int mincut = std::numeric_limits<int>::max();
    real_t imbalance = 0.0;
    bool balanced = false;
    int bestTrial = -1;
};

BisectResult Bisect(Graph& g, const BisectOptions& opts={});

namespace bisect {
    void  ComputeParams(Graph& g);
    void  GrowInit(Graph& g, int seed, int t0, int t1);
    int   FMPass(Graph& g, int maxPW0, int maxPW1);
    int   BalanceFix(Graph& g, int t0, int t1, real_t ub);
    real_t Imbalance(const Graph& g, int t0, int t1);
}

#endif