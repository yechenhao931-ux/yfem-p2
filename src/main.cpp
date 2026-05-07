#include <mpi.h>
#include <iostream>
#include <parallel/par.hpp>
#include <common/optparser.hpp>
char buf[64];


int main(int argc, char** argv){
    const char *mesh_file = "../resource/mashbin.msh";
    int parts = 4;


    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
    args.AddOption(&parts, "-p", "--part",
                  "Number of partitions.");
    args.Parse();
    int rank, size;
    int tag = 0;
    MPI_Status status;
    MPI_Comm comm = MPI_COMM_WORLD;
    // 初始化MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 构造一个网格
    Mesh mesh;
    mesh.num_vertices = 4;
    Element e1 = Element(0,1,2);
    Element e2 = Element(1,3,2);
    mesh.elements.push_back(std::move(e1));
    mesh.elements.push_back(std::move(e2));
    mesh.Print();

    // 并行网格
    ParMesh pmesh(comm,mesh);

    // 并行FESpace
    ParFiniteElementSpace fes(pmesh);

    // 并行装配
    ParBilinearForm a(fes);
    a.Assemble();

    if (rank == 1) {
        int send_data = 12;
        MPI_Send(&send_data , 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
        // 阻塞直到接收方开始接收（或缓冲区可用）
    } else if (rank == 0) {
        int recv_data;
        double wait_start = MPI_Wtime();

        MPI_Recv(&recv_data, 1, MPI_INT, 1, tag, MPI_COMM_WORLD, &status);
        // 阻塞直到收到数据
        printf("sync finished\n");
    }
    if(rank == 0){
        std::cout << "Local assembled entries:\n";
        for(auto &kv : a.mat){
            std::cout << "("
                      << kv.first.first << ","
                      << kv.first.second << ") = "
                      << kv.second << "\n";
         }


    }


    MPI_Finalize();

    std::cout << "コリアタウン"  << std::endl;
    return 0;
}