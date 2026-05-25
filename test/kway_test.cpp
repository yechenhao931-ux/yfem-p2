#include <iostream>
#include <graph/graph.hpp>
#include <partition/kwaypartition.hpp>
#include <partition/geokwaypartition.hpp>
#include <common/macro.hpp>
#include <io/save.hpp>
#include <common/optparser.hpp>
#include <cstdio>


int main(int argc,char* argv[]){
    const char *mesh_file = "../resource/box.mesh";
    int n_part = 4;

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
    args.AddOption(&n_part, "-p", "--part",
                  "Number of partitions.");
    args.Parse();
    if (!args.Good())
    {
        args.PrintUsage(std::cout);
        return 1;
    }
    args.PrintOptions(std::cout);

    std::string path = mesh_file;
    // 1. 最後の '/' を見つける
    size_t last_slash = path.find_last_of("/\\");
    size_t start = (last_slash == std::string::npos) ? 0 : last_slash + 1;

    // 2. 最後の '.' を見つける
    size_t last_dot = path.find_last_of('.');
    
    // 3. 切り出し
    std::string file_name = path.substr(start, last_dot - start);
    MFEMMesh10 mfemmesh;
    mfemmesh.read_mesh(path.c_str());
    Graph g = Graph(mfemmesh, 3);
    //RBOptions o; o.nparts=8; o.nTrials=10; o.nFMPasses=6; o.ubFactor=1.03; o.seed=42; o.postRefine=PostRefineMode::EdgeCut; o.postKwayIter=5; o.verbose=true;
    KwayOptions o;
    o.objective = KwayObjective::EdgeCut;
    o.nparts = n_part;
    o.verbose = true;



    std::stringstream ss;
    ss  << file_name << ".part";
    std::string s = ss.str(); // "Score: 95.57"

    auto r = KwayPartition(g, o);
    SavePartitionResult(g,n_part,s);


    // Graph g1 = Graph(reader);
    // GeoKwayOptions o1; 
    // o1.nparts=parts;
    // o1.verbose =true;

    // GeoKwayResult r1=GeoKwayPartition(g1,o1);
    // std::printf("\n[GeoKway Volume α=0.7] cut=%d  vol=%d  imb=%.2f%%  beta=%.4f\n",
    //                 r1.mincut, r1.minvol, r1.maxImbalance*100, r1.betaUsed);
    // SavePartitionResult(g1,"outputg.data");
    printf("finish\n");
    
    return 0;
}