// header file which declares the util functions
#ifndef VERTEX_H //inlcude guard
#define VERTEX_H
#include "Index.h"
#include "mpi.h"
#include "globals.h"
#include <complex>
#include <vector>

class Vertex {

    const int Nl_internal_;
    std::vector<std::complex<double>> vertex_;

public:
    Vertex(const int Nl);
    Vertex(std::complex<double> initial_value, const int Nl);

    const std::complex<double> operator()(const IndexF indexk1, const IndexF indexk2, const IndexIBZ indexq1) const noexcept
    {

        return vertex_[indexq1.ind() * Nw * Nl_internal_ * Nw * NxIBZ + indexk2.ind() * Nw * NxIBZ + indexk1.ind()];

    }
};


#endif // vertex.h
