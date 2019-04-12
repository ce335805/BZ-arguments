#include "FF.h"
#include "GreenFunction.h"
#include "Index.h"
#include "globals.h"
#include "mpi.h"
#include "output.inc"
#include "util.h"
#include <assert.h>
#include <chrono>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>

#define MKL_Complex16 std::complex<double>
#include "mkl.h"

int main()
{

    static_assert(N > 0, "specified 0 sites or negative");
    static_assert(Nw > 0, "specified 0 prositive frequencies or negative");
    static_assert(N == Nx * Ny, "Nx * Ny neq N");

    MPI_Init(nullptr, nullptr);

    int worldsize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldsize);
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    std::cout << "I am rank " << myrank << " of " << worldsize << '\n';

    std::string file_name = "data/inidata/pick_list_10.init";
    //new FF object
    FF ff(file_name);

    //extract Nl from FF class
    const int Nl { ff.get_NL() };

    std::cout << "model parameters: " << '\n'
              << "Nx x Ny = " << Nx << " x " << Ny << '\n'
              << "Nw = " << Nw << '\n'
              << "Nl = " << Nl << '\n'
              << "Beta = " << beta << '\n';

    IndexF ::Nl = Nl;
    IndexB ::Nl = Nl;

    std::vector<std::complex<double>> sigma(Nx * Ny * 2 * Nw, std::complex<double>(0.0, 0.0));

    //simple GreenFunction object
    GreenFunction Gk(2, 0.0, sigma.data(), 1);

    //output Green's function to look at it
    if (myrank == 0)
        output::output_BZ(Gk.data(), Nx * Nx * 2 * Nw, "Green_file.h5", "Gk");

    MPI_Finalize();
}
