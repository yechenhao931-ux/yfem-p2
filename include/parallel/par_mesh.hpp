#ifndef PARMESH_HPP
#define PARMESH_HPP
#include <mpi.h>
#include <set>
#include "mesh.hpp"
class ParMesh{
public:
    MPI_Comm comm;
    int rank, size;
    std::vector<Element> local_elements;
    std::set<int> shared_vertices;

    ParMesh() = default;

    ParMesh(MPI_Comm c, const Mesh &mesh):comm(c){
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);
        std::cout << "rank:" <<  rank << std::endl;
        Partition(mesh);
        BuildSharedVertices();
    }
    ~ParMesh() = default;

private:
    void Partition(const Mesh &mesh)
    {
        // 最简单的按 element id 划分
        for (int i = 0; i < (int)mesh.elements.size(); i++)
        {
            
            if (i % size == rank){
                local_elements.push_back(mesh.elements[i]);
                std::cout << "rank[" << rank << "] local_elements:" <<  &local_elements << "push" << i  << std::endl ;
               }   
        }
    }

    void BuildSharedVertices()
    {
        for (auto &e : local_elements){
            
            for (int i = 0; i < 3; i++){
                shared_vertices.insert(e.vertice[i]);
                std::cout << "rank [" << rank << "] shared_vertice insert:"   << e.vertice[i] << std::endl;
            }     
        }
        for(auto &v : shared_vertices){
            std::cout << "rank [" << rank << "] shared_vertice :"   << v << std::endl;
        }
    }
    
};

#endif