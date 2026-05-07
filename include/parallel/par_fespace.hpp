#ifndef PAR_FESPACE_HPP
#define PAR_FESPACE_HPP
#include "par_mesh.hpp"
#include <map>
class ParFiniteElementSpace{
public:
    ParMesh &pmesh;

    std::map<int, int> local_dof;
    std::map<int, int> global_dof;

    int local_ndofs;
    int global_offset;

    ParFiniteElementSpace(ParMesh &pmesh_):pmesh(pmesh_){
        BuildLocalDofs();
        BuildGlobalDofs();
    }
    ~ParFiniteElementSpace() = default;

private:
    void BuildLocalDofs()
    {
        int cnt = 0;
        for (int v : pmesh.shared_vertices)
            local_dof[v] = cnt++;

        local_ndofs = cnt;
    }

    void BuildGlobalDofs()
    {
        MPI_Scan(&local_ndofs, &global_offset,
                 1, MPI_INT, MPI_SUM, pmesh.comm);

        global_offset -= local_ndofs;

        for (auto &kv : local_dof){
            global_dof[kv.first] = global_offset + kv.second;
        }
            
    }
};

#endif