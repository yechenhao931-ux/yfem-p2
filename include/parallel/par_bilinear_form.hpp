#ifndef PAR_BILINEAR_FORM_HPP
#define PAR_BILINEAR_FORM_HPP
#include "par_fespace.hpp"
#include <map>

class ParBilinearForm{
public:
    ParFiniteElementSpace &fes;
    std::map<std::pair<int,int>,double> mat;

public:
    ParBilinearForm(ParFiniteElementSpace &f):fes(f){

    }
    ~ParBilinearForm() = default;

    void Assemble(){
        for(auto &e : fes.pmesh.local_elements){
            AssembleElement(e);
        }
    }
private:
    void AssembleElement(const Element &e){
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                int gi = fes.global_dof[e.vertice[i]];
                int gj = fes.global_dof[e.vertice[j]];

                mat[{gi,gj}] += (i==j ? 2.0 : -1.0);
            }
        }
    }
};
#endif