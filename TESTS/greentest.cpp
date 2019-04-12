#include "GreenFunction.h"
#include "Index.h"
#include "globals.h"
#include "util.h"
#include "gtest/gtest.h"

#include <complex>
#include <vector>

namespace {

TEST(GreenTest, initialization)
{

    std::vector<std::complex<double>> sigma(Nx * Ny * 4 * Nw, std::complex<double> (0.0, 0.0));

    //to check decent initialization check symmetries of bare Green's function
    GreenFunction Gk1(1, 0.0, sigma.data(), 1);

    std::vector<int> k { 2, 1 };

    for (int sym = 0; sym < 8; ++sym) {
        int f = 1;
        std::vector<int> k_sym = sym_sqr(sym, k);
        IndexFK index(k[0], k[1], f);
        IndexFK index_sym(k_sym[0], k_sym[1], f);

        EXPECT_NEAR(Gk1.val(index).real(), Gk1.val(index_sym).real(), 1e-15) << "RE(GK) != RE(Sym(Gk)); sym = " << sym << " in small box";
        EXPECT_NEAR(Gk1.val(index).imag(), Gk1.val(index_sym).imag(), 1e-15) << "IM(GK) != IM(Sym(Gk)); sym = " << sym << " in small box";
    }

    GreenFunction Gk2(2, 0.0, sigma.data(), 1);

    for (int sym = 0; sym < 8; ++sym) {
        int f = Nw + 3;
        std::vector<int> k_sym = sym_sqr(sym, k);
        IndexFK index(k[0], k[1], f);
        IndexFK index_sym(k_sym[0], k_sym[1], f);

        EXPECT_NEAR(Gk2.val(index).real(), Gk2.val(index_sym).real(), 1e-15) << "RE(GK) != RE(Sym(Gk)); sym = " << sym << "in extended box";
        EXPECT_NEAR(Gk2.val(index).imag(), Gk2.val(index_sym).imag(), 1e-15) << "IM(GK) != IM(Sym(Gk)); sym = " << sym << "in extended box";
    }

    IndexFK index_neg(0, 0, -3, 1);
    IndexFK index_pos(0, 0, 3, 1);

		//std::cout << '\n';
		//std::cout << "ind() = " << index_neg.ind() << '\n' 
		//					<< "f_ = " << index_neg.get_f() << '\n' 
		//					<< "ky_ = " << index_neg.get_ky() << '\n' 
		//					<< "kx_ = " << index_neg.get_kx() << '\n';
		//std::cout << '\n';
		//std::cout << "ind() = " << index_pos.ind() << '\n' 
		//					<< "f_ = " << index_pos.get_f() << '\n' 
		//					<< "ky_ = " << index_pos.get_ky() << '\n' 
		//					<< "kx_ = " << index_pos.get_kx() << '\n';

		//std::cout << Gk1.data()[index_neg.ind()] << '\n';
		//std::cout << Gk1.data()[index_pos.ind()] << '\n';

    //check symmetry and antisymmetry in frequency arguments
    EXPECT_NEAR(Gk1.val(index_neg).real(), Gk1.val(index_pos).real(), 1e-15) << "RE(GK(-f)) != RE(GK(f)) - in small box";
    EXPECT_NEAR(Gk1.val(index_neg).imag(), -Gk1.val(index_pos).imag(), 1e-15)<< "IM(GK(-f)) != -IM(GK(f)) - in small box";

    IndexFK index_neg2(0, 0, -3, 2);
    IndexFK index_pos2(0, 0, 3, 2);

    EXPECT_NEAR(Gk2.val(index_neg2).real(), Gk2.val(index_pos2).real(), 1e-15)<< "RE(GK(-f)) != RE(GK(f)) - in extended box";
    EXPECT_NEAR(Gk2.val(index_neg2).imag(), -Gk2.val(index_pos2).imag(), 1e-15)<< "IM(GK(-f)) != -IM(GK(f)) - in extended box";

}
}
