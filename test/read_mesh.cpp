#include "graph/graph.hpp"
#include "partition/kwaypartition.hpp"
#include "partition/geokwaypartition.hpp"
#include "io/save.hpp"
#include "common/optparser.hpp"

#include <random>

int main(int argc, char *argv[])
{

    const char *mesh_file_path = "../resource/box.mesh";
    int n_part = 8;
    char k_or_g = 'g';
    double alpha = 0.7;

    // 乱数を生成する
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1024);
    int seed = dis(gen);

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file_path, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&n_part, "-p", "--part",
                   "Number of partitions.");
    args.AddOption(&k_or_g, "-w", "--way",
                   "Way to divide");
    args.AddOption(&seed, "-s", "--seed",
                   "set the seed");
    args.AddOption(&alpha, "-a", "--alpha",
                   " 拓扑权重（0~1）");
    args.Parse();
    if (!args.Good())
    {
        args.PrintUsage(std::cout);
        return 1;
    }
    args.PrintOptions(std::cout);

    std::string path = mesh_file_path;
    // 1. 最後の '/' を見つける
    size_t last_slash = path.find_last_of("/\\");
    size_t start = (last_slash == std::string::npos) ? 0 : last_slash + 1;

    // 2. 最後の '.' を見つける
    size_t last_dot = path.find_last_of('.');

    // 3. 切り出し
    std::string file_name = path.substr(start, last_dot - start);

    // 4. 保存文件的名字
    std::stringstream ss;
    ss << file_name << ".part";
    std::string s = ss.str();

    int npart = 8;
    MFEMMesh10 mmesh;

    mmesh.read_mesh(mesh_file_path);
    if (k_or_g == 'g')
    {
        // GKWAY划分

        Graph g1(mmesh, 3);
        GeoKwayOptions o1;
        o1.nparts = npart;
        o1.verbose = false;
        o1.seed = seed;
        o1.alpha = alpha;

        GeoKwayResult r1 = GeoKwayPartition(g1, o1);
        std::printf("\n[GeoKway Volume α=%f] cut=%d  vol=%d  imb=%.2f%%  beta=%.4f\n",
                    alpha,r1.mincut, r1.minvol, r1.maxImbalance * 100, r1.betaUsed);


    int count[16] = {0};

        for (int x : g1.where)
        {
            count[x]++;
        }
        for(int i =0; i < npart ;i++){
            printf("p[%d]%d\t",i,count[i]);
        }
        SavePartitionResult(g1, npart, s);
    }
    else if (k_or_g == 'k')
    {

        // Kway划分
        Graph g2 = Graph(mmesh, 3);
        KwayOptions o2;
        o2.objective = KwayObjective::EdgeCut;
        o2.nparts = npart;
        o2.verbose = false;
        o2.seed = seed;

        KwayResult r2 = KwayPartition(g2, o2);
        std::printf("\n[Kway] cut=%d  vol=%d  imb=%.2f%% \n",
                    r2.mincut, r2.minvol, r2.maxImbalance * 100);
    }

    std::cout << "finished\n";
    return 0;
}
