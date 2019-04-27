#include "Index.h"
#include "mpi.h"
#include "gtest/gtest.h"

namespace {

TEST(IndexIBZTest, initialization)
{

    int Nx = 10;
    int Ny = 10;
    int NxIBZ = (Nx / 2 + 1) * (Nx / 2 + 2);
    int Nw = 16;

    IndexIBZ::Nx_ = Nx;
    IndexIBZ::Ny_ = Ny;
    IndexIBZ::NxIBZ_ = (IndexIBZ::Nx_ / 2 + 1) * (IndexIBZ::Nx_ / 2 + 2) ;
    IndexIBZ::Nw = Nw;


    int worldsize;
    int myrank;

    MPI_Comm_size(MPI_COMM_WORLD, &worldsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    start dummy_start;
    end dummy_end;

    //IndexIBZ
    IndexIBZ ibz1;
    EXPECT_EQ(ibz1.ind(), (Nw * NxIBZ) / worldsize * myrank) << "IndexIBZ: initialization with no arguments";

    IndexIBZ ibz2(dummy_start);
    EXPECT_EQ(ibz2.ind(), (Nw * NxIBZ) / worldsize * myrank) << "IndexIBZ: initialization with start argument";

    IndexIBZ ibz3(dummy_end);
    EXPECT_EQ(ibz3.ind(), (Nw * NxIBZ / worldsize) * (myrank + 1)) << "IndexIBZ: initialization with end argument";
}

TEST(IndexIBZTest, operators)
{

    int Nx = 10;
    int Ny = 10;
    int NxIBZ = (Nx / 2 + 1) * (Nx / 2 + 2);
    int Nw = 16;

    IndexIBZ::Nx_ = Nx;
    IndexIBZ::Ny_ = Ny;
    IndexIBZ::NxIBZ_ = (IndexIBZ::Nx_ / 2 + 1) * (IndexIBZ::Nx_ / 2 + 2) ;
    IndexIBZ::Nw = Nw;

    int worldsize;
    int myrank;

    MPI_Comm_size(MPI_COMM_WORLD, &worldsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    start dummy_start;
    end dummy_end;

    //IndexIBZ
    IndexIBZ ibz1;
    ++ibz1;
    EXPECT_EQ(ibz1.ind(), (Nw * NxIBZ) / worldsize * myrank + 1) << "IndexIBZ: incrementation operator '++()' - 1";
    ++(++ibz1);
    EXPECT_EQ(ibz1.ind(), (Nw * NxIBZ) / worldsize * myrank + 3) << "IndexIBZ: incrementation operator '++()' - 3";
    ++(++ibz1);
    EXPECT_EQ(ibz1.ind(), (Nw * NxIBZ) / worldsize * myrank + 5) << "IndexIBZ: incrementation operator '++()' - 5";

    IndexIBZ ibz2;
    EXPECT_GT(ibz1, ibz2) << "IndexBK: comparison operator '>'";
    EXPECT_LT(ibz2, ibz1) << "IndexBK: comparison operator '<'";
}

} //namespace
