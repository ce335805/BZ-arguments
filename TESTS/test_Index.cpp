#include "Index.h"
#include "gtest/gtest.h"

namespace {

TEST(IndexTest, initialization)
{

    int Nx = 10;
    int Ny = 10;
    int Nl = Nx * Ny;
    int Nw = 16;

    IndexF::Nl = Nl;
    IndexB::Nl = Nl;

    IndexFK::Nx_ = Nx;
    IndexFK::Ny_ = Ny;
    IndexBK::Nx_ = Nx;
    IndexBK::Ny_ = Ny;

    IndexF::Nw = Nw;
    IndexB::Nw = Nw;
    IndexFK::Nw = Nw;
    IndexBK::Nw = Nw;


    //IndexF
    IndexF f1;
    EXPECT_EQ(f1.ind(), 0) << "IndexF: initialization with no arguments";
    start dummy_start;

    IndexF f2(dummy_start);
    EXPECT_EQ(f2.ind(), 0) << "IndexF: initialization with start argument";
    end dummy_end;

    IndexF f3(dummy_end);
    EXPECT_EQ(f3.ind(), 2 * Nw * Nl) << "IndexF: initialization with end argument";

    //initialization with extended frequency range
    IndexF f4(dummy_end, 2);
    EXPECT_EQ(f4.ind(), 4 * Nw * Nl) << "IndexF: initialization with end argument";

    //initialization via copy constructor
    IndexF f5(f4);
    EXPECT_EQ(f5.ind(), 4 * Nw * Nl) << "IndexF: initialization via copy constructor";

    //initialization via copy constructor with changed frequency range - trivial case
    IndexF f6(f1, 1);
    EXPECT_EQ(f6.ind(), 0) << "IndexF: initialization via copy constructor with changed frequency range, Nw_range = 1";

    //initialization via copy constructor with changed frequency range
    IndexF f7(f1, 2);
    EXPECT_EQ(f7.ind(), (Nw - 1) * Nl + Nl) << "IndexF: initialization via copy constructor with changed frequency range, Nw_range = 2";

    //check narrowed down index
    EXPECT_EQ(f7.ind(1), 0) << "IndexF: check narrowed down index";

    //IndexB
    IndexB b1;
    EXPECT_EQ(b1.ind(), 0) << "IndexB: initialization with no arguments";

    IndexB b2(dummy_start);
    EXPECT_EQ(b2.ind(), 0) << "IndexB: initialization with start argument";

    IndexB b3(dummy_end);
    EXPECT_EQ(b3.ind(), Nw * Nl) << "IndexB: initialization with end argument";

    //initialization with extended frequency range
    IndexB b4(dummy_end, 2);
    EXPECT_EQ(b4.ind(), 2 * Nw * Nl) << "IndexB: initialization with end argument";

    //initialization via copy constructor
    IndexB b5(b4);
    EXPECT_EQ(b5.ind(), 2 * Nw * Nl) << "IndexB: initialization via copy constructor";

    //initialization via copy constructor with changed frequency range - trivial case
    IndexB b6(b1, 1);
    EXPECT_EQ(b6.ind(), 0) << "IndexB: initialization via copy constructor with changed frequency range, Nw_range = 1";

    //initialization via copy constructor with changed frequency range - also tricial here
    IndexB b7(b1, 2);
    EXPECT_EQ(b7.ind(), 0) << "IndexB: initialization via copy constructor with changed frequency range, Nw_range = 2";

    //check narrowed down index
    EXPECT_EQ(b7.ind(1), 0) << "IndexF: check narrowed down index";

    //IndexFK
    IndexFK fk1;
    EXPECT_EQ(fk1.ind(), 0) << "IndexFK: initialization with no arguments";

    IndexFK fk2(dummy_start);
    EXPECT_EQ(fk2.ind(), 0) << "IndexFK: initialization with start argument";

    IndexFK fk3(dummy_end);
    EXPECT_EQ(fk3.ind(), 2 * Nw * Nx * Ny) << "IndexFK: initialization with end argument";

    //initialization with extended frequency range
    IndexFK fk4(dummy_end, 2);
    EXPECT_EQ(fk4.ind(), 4 * Nw * Nx * Ny) << "IndexFK: initialization with end argument";

    //initialization via copy constructor
    IndexFK fk5(fk4);
    EXPECT_EQ(fk5.ind(), 4 * Nw * Nx * Ny) << "IndexFK: initialization via copy constructor";

    //initialization via copy constructor with changed frequency range - trivial case
    IndexFK fk6(fk1, 1);
    EXPECT_EQ(fk6.ind(), 0) << "IndexFK: initialization via copy constructor with changed frequency range, Nw_range = 1";

    //initialization via copy constructor with changed frequency range
    IndexFK fk7(fk1, 2);
    EXPECT_EQ(fk7.ind(), (Nw - 1) * Nx * Ny + (Ny - 1) * Nx + Nx) << "IndexFK: initialization via copy constructor with changed frequency range, Nw_range = 2";

    //check narrowed down index
    EXPECT_EQ(f7.ind(1), 0) << "IndexFK: check narrowed down index";

    //IndexBK
    IndexBK bk1;
    EXPECT_EQ(bk1.ind(), 0) << "IndexBK: initialization with no arguments";

    IndexBK bk2(dummy_start);
    EXPECT_EQ(bk2.ind(), 0) << "IndexBK: initialization with start argument";

    IndexBK bk3(dummy_end);
    EXPECT_EQ(bk3.ind(), Nw * Nx * Ny) << "IndexBK: initialization with end argument";

    //initialization with extended frequency range
    IndexBK bk4(dummy_end, 2);
    EXPECT_EQ(bk4.ind(), 2 * Nw * Nx * Ny) << "IndexBK: initialization with end argument";

    //initialization via copy constructor
    IndexBK bk5(bk4);
    EXPECT_EQ(bk5.ind(), 2 * Nw * Nx * Ny) << "IndexB: initialization via copy constructor";

    //initialization via copy constructor with changed frequency range - trivial case
    IndexBK bk6(bk1, 1);
    EXPECT_EQ(bk6.ind(), 0) << "IndexB: initialization via copy constructor with changed frequency range, Nw_range = 1";

    //initialization via copy constructor with changed frequency range - also tricial here
    IndexBK bk7(bk1, 2);
    EXPECT_EQ(bk7.ind(), 0) << "IndexB: initialization via copy constructor with changed frequency range, Nw_range = 2";

    //check narrowed down index
    EXPECT_EQ(bk7.ind(1), 0) << "IndexF: check narrowed down index";
}

TEST(IndexTest, operators)
{

    int Nx = 10;
    int Ny = 10;
    int Nl = Nx * Ny;
    int Nw = 16;

    IndexF::Nl = Nl;
    IndexB::Nl = Nl;

    IndexFK::Nx_ = Nx;
    IndexFK::Ny_ = Ny;
    IndexBK::Nx_ = Nx;
    IndexBK::Ny_ = Ny;

    IndexF::Nw = Nw;
    IndexB::Nw = Nw;
    IndexFK::Nw = Nw;
    IndexBK::Nw = Nw;

    //IndexF
    IndexF f1;
    ++f1;
    EXPECT_EQ(f1.ind(), 1) << "IndexF: incrementation operator '++()'";
    ++(++f1);
    EXPECT_EQ(f1.ind(), 3) << "IndexF: incrementation operator '++()'";

    IndexF f2;
    EXPECT_GT(f1, f2) << "IndexF: comparison operator '>'";
    EXPECT_LT(f2, f1) << "IndexF: comparison operator '<'";

    //IndexB
    IndexB b1;
    ++b1;
    EXPECT_EQ(b1.ind(), 1) << "IndexB: incrementation operator '++()'";
    ++(++b1);
    EXPECT_EQ(b1.ind(), 3) << "IndexB: incrementation operator '++()'";

    IndexB b2;
    EXPECT_GT(b1, b2) << "IndexB: comparison operator '>'";
    EXPECT_LT(b2, b1) << "IndexB: comparison operator '<'";

    //IndexFK
    IndexFK fk1;
    ++fk1;
    EXPECT_EQ(fk1.ind(), 1) << "IndexFK: incrementation operator '++()'";
    ++(++fk1);
    EXPECT_EQ(fk1.ind(), 3) << "IndexFK: incrementation operator '++()'";

    IndexFK fk2;
    EXPECT_GT(fk1, fk2) << "IndexFK: comparison operator '>'";
    EXPECT_LT(fk2, fk1) << "IndexFK: comparison operator '<'";

    //IndexBK
    IndexBK bk1;
    ++bk1;
    EXPECT_EQ(bk1.ind(), 1) << "IndexBK: incrementation operator '++()'";
    ++(++bk1);
    EXPECT_EQ(bk1.ind(), 3) << "IndexBK: incrementation operator '++()'";

    IndexBK bk2;
    EXPECT_GT(bk1, bk2) << "IndexBK: comparison operator '>'";
    EXPECT_LT(bk2, bk1) << "IndexBK: comparison operator '<'";
}
}
