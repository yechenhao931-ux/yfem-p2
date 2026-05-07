#ifndef KWAYREFINE_HPP
#define KWAYREFINE_HPP
#include "../graph/graph.hpp"
void ComputeCkrinfo(Graph &g, int nparts);
void ComputeVkrinfo(Graph &g, int nparts);
struct KwayFMOpts
{
    int nparts = 4;
    int nIter = 10;
    real_t ubFactor = 1.03;
    bool verbose = false;
};
int KwayFMCut(Graph &g, const KwayFMOpts &opts);
int KwayFMVol(Graph &g, const KwayFMOpts &opts);

#endif