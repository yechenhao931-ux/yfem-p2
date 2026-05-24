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
// A-3: 独立集并行精化（G-kway）——替代串行 FM，仅优化 EdgeCut
int IndepSetRefineCut(Graph &g, const KwayFMOpts &opts);

#endif