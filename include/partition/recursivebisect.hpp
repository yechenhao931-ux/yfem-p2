#ifndef RECURSIVEBISECT_HPP
#define RECURSIVEBISECT_HPP

#include "../graph/graph.hpp"
#include "bisectfm.hpp"
#include "kwayrefine.hpp"

enum class PostRefineMode
{
    None,
    EdgeCut,
    Volume,
    Both
};
struct RBOptions
{
    int nparts = 4, nTrials = 10, nFMPasses = 8;
    real_t ubFactor = 1.03;
    int seed = 42;
    PostRefineMode postRefine = PostRefineMode::EdgeCut;
    int postKwayIter = 10;
    bool verbose = false;
};
struct RBResult
{
    int totalCut = 0, totalVol = 0;
    real_t maxImbalance = 0., avgImbalance = 0.;
    bool balanced = false;
    std::vector<int> partWeights;
};
RBResult RecursiveBisect(Graph &graph, const RBOptions &opts = {});
#endif