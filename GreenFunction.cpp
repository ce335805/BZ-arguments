#include "GreenFunction.h"
#include "globals.h"

#include <assert.h>
#include <cmath>

GreenFunction::GreenFunction(const int frequencyextension_init, const double mu_init, const std::complex<double>* selfenergy_ptr, const int sig_frequencyextension_init)
    : frequencyextension_ { frequencyextension_init }
    , mu_ { mu_init }
    , selfenergy_ { selfenergy_ptr }
    , sig_frequencyextension_ { sig_frequencyextension_init }
{

    //cannot use bigger window for selfenergy than for Green's function
    assert(frequencyextension_ >= sig_frequencyextension_);
    //allocate data for the actual Greens function
    Gk = std::vector<std::complex<double>>(Nx * Ny * 2 * Nw * frequencyextension_, std::complex<double>(0.0, 0.0));

    //case that Green's function is to be computed in the same frequency box as the selfenergy is available
    if (frequencyextension_ == sig_frequencyextension_) {

        for (IndexFK index(dummy_start_, frequencyextension_); index < IndexFK(dummy_end_, frequencyextension_); ++index) {

            //            std::cout << "ind() = " << index.ind() << '\n'
            //                      << "f_ = " << index.get_f() << '\n'
            //                      << "ky_ = " << index.get_ky() << '\n'
            //                      << "kx_ = " << index.get_kx() << '\n';

            double px = index.get_kx() * 2.0 * PI / double(Nx);
            double py = index.get_ky() * 2.0 * PI / double(Ny);
            std::complex<double> freq(0.0, PI / beta * index.get_f());

            Gk[index.ind()] = 1.0 / (freq + mu_ - 2.0 * (cos(px) + cos(py)) - selfenergy_ptr[index.ind()]);
        }
    } else {

        //        std::cout << "constructing Green's function in bigger box than selfenergy" << '\n';

        //fill up values where selfenergy is not known - first below known selfenergy ...
        for (IndexFK index(dummy_start_, frequencyextension_);
             index < IndexFK(dummy_start_, sig_frequencyextension_);
             ++index) {

            double px = index.get_kx() * 2.0 * PI / double(Nx);
            double py = index.get_ky() * 2.0 * PI / double(Ny);
            std::complex<double> freq(0.0, PI / beta * index.get_f());

            Gk[index.ind()] = 1.0 / (freq + mu_ - 2.0 * (cos(px) + cos(py)));
        }

        // ... then above
        for (IndexFK index(IndexFK(dummy_end_, sig_frequencyextension_), frequencyextension_);
             index < IndexFK(dummy_end_, frequencyextension_);
             ++index) {

            double px = index.get_kx() * 2.0 * PI / double(Nx);
            double py = index.get_ky() * 2.0 * PI / double(Ny);
            std::complex<double> freq(0.0, PI / beta * index.get_f());

            Gk[index.ind()] = 1.0 / (freq + mu_ - 2.0 * (cos(px) + cos(py)));
        }

        //then compute part where selfergy is known
        for (IndexFK index(IndexFK(dummy_start_, sig_frequencyextension_), frequencyextension_);
             index < IndexFK(dummy_end_, sig_frequencyextension_);
             ++index) {

            double px = index.get_kx() * 2.0 * PI / double(Nx);
            double py = index.get_ky() * 2.0 * PI / double(Ny);
            std::complex<double> freq(0.0, PI / beta * index.get_f());

            Gk[index.ind()]
                = 1.0 / (freq + mu_ - 2.0 * (cos(px) + cos(py)) - selfenergy_ptr[index.ind(sig_frequencyextension_)]);
        }
    }
}

//function to return G(k + q) in a vertain k-range
//to have G(k + q) data contiguously in memory
const std::vector<std::complex<double>> GreenFunction::Gkpq(const IndexBK indexq, const int Nw_range) const
{

    //check that boxsizes are ok
    assert(Nw_range <= frequencyextension_);
    assert(2 * indexq.get_f() + 2 * Nw_range - 1 <= frequencyextension_ * Nw - 1);

    std::vector<std::complex<double>> Gkpq(2 * Nw_range * Nw * Nx * Ny, std::complex<double>(0.0, 0.0));

    //case that one wants G(k+q) in the same box as sigma is known
    for (IndexFK indexk(dummy_start_, Nw_range); indexk < IndexFK(dummy_end_, Nw_range); ++indexk) {

        Gkpq[indexk.ind()] = Gk[(indexk + indexq).ind()];
    }

    return Gkpq;
}

//function to return G(k + q) in a vertain k-range
//to have G(k + q) data contiguously in memory
const std::vector<std::complex<double>> GreenFunction::Gkpq(const IndexIBZ indexq, const int Nw_range) const
{

    //check that boxsizes are ok
    assert(Nw_range <= frequencyextension_);
    assert(2 * indexq.get_f() + 2 * Nw_range - 1 <= frequencyextension_ * Nw - 1);

    std::vector<std::complex<double>> Gkpq(2 * Nw_range * Nw * Nx * Ny, std::complex<double>(0.0, 0.0));

    //case that one wants G(k+q) in the same box as sigma is known
    for (IndexFK indexk(dummy_start_, Nw_range); indexk < IndexFK(dummy_end_, Nw_range); ++indexk) {

        Gkpq[indexk.ind()] = Gk[(indexq + indexk).ind()];
    }

    return Gkpq;
}
