#ifndef PAR_LINALG_HPP
#define PAR_LINALG_HPP

#include <mpi.h>
#include <vector>

inline double ParDot(MPI_Comm comm,
                     const std::vector<double> &x,
                     const std::vector<double> &y)
{
    double local = 0.0;
    for (size_t i = 0; i < x.size(); i++)
        local += x[i] * y[i];

    double global;
    MPI_Allreduce(&local, &global, 1,
                  MPI_DOUBLE, MPI_SUM, comm);
    return global;
}
#endif