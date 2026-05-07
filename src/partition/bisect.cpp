#include "partition/bisect.hpp"
#include "graph/graph.hpp"
#include <limits>
int ComputeCut(const Graph& g){
    if(g.where.size() == 0){
        return 0;
    }

    int cut = 0;
    for(int v = 0; v < g.nvtxs; v++){
        for(int e = g.xadj[v]; e < g.xadj[v + 1]; e++){
            int u = g.adjncy[e];

            if(g.where[v] != g.where[u]){
                cut ++;
            }
        }
    }
    return cut / 2;
}

void BestBisection(Graph& g, int ncuts = 8){
    int best_cut = std::numeric_limits<int>::max();
    std::vector<int> best_where;

    for(int i = 0; i < ncuts; i++){
        Graph tmp = g;
        assert(0);
    }


    g.where = best_where;
    
}

Graph BuildSubGraph(const Graph&g, int part){
    std::vector<int> map(g.nvtxs, -1);

    int nvtxs_sub = 0;
    for(int v = 0; v < g.nvtxs; v++){
        if(g.where[v] == part){
            map[v] = nvtxs_sub++;
        }
    }

    int edges = 0;
    for(int v = 0; v < g.nvtxs; v++){
        if(map[v] == -1){
            continue;
        }

        for(int e = g.xadj[v]; e < g.xadj[v + 1]; e++){
            int u = g.adjncy[e];
            if(map[u] != -1){
                edges++;
            }
        }
    }

    // 初始化子图
    Graph sub;

    sub.nvtxs = nvtxs_sub;
    sub.nedges = edges;
    sub.ncon = g.ncon;
    
    sub.xadj.resize(nvtxs_sub + 1);
    sub.adjncy.resize(edges);
    sub.adjwgt.assign(edges, 1);

    int edgeptr = 0;

    for(int v = 0; v < g.nvtxs; v++){
        if(map[v] == -1){
            continue;
        }

        int newv = map[v];

        sub.xadj[newv] = edgeptr;

        for(int e = g.xadj[v]; e < g.xadj[v + 1]; e++){
            int u = g.adjncy[e];

            if(map[u] != -1){
                sub.adjncy[edgeptr++] = map[u];
            }
        }
    }

    sub.xadj[nvtxs_sub] = edgeptr;

    return sub;
}

void RecursivePartition(Graph &g, int nparts){
    if(nparts == 1){
        return;
    }

    BestBisection(g);

    Graph g0 = BuildSubGraph(g,0);
    Graph g1 = BuildSubGraph(g,1);

    RecursivePartition(g0, nparts/2);
    RecursivePartition(g1, nparts/2);

    int offset = nparts / 2;

    for(int v=0; v<g.nvtxs; v++)
    {
        if(g.where[v]==0)
            g.where[v] = g0.where[v];
        else
            g.where[v] = g1.where[v] + offset;
    }
}


