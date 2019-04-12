// header file which declares the util functions
#ifndef GREENFUNCTION_H //inlcude guard
#define GREENFUNCTION_H

#include "Index.h"
#include "globals.h"

#include <complex>
#include <vector>

class GreenFunction {

private:
    //multiples of Nw for which Gk should be saved
    const int frequencyextension_;
    //chemical potential
    const double mu_;
    //for initialization of IndexX
    const start dummy_start_ = start();
    const end dummy_end_ = end();

    //pointer to selfenergy
    const std::complex<double>* selfenergy_;
    //multiples of frequencies for which selfenergy is known
    const int sig_frequencyextension_;

    //actual Greens function
    std::vector<std::complex<double>> Gk;

public:
    GreenFunction(const int frequencyextesion_init,
        const double mu_init,
        const std::complex<double>* selfenergy_ptr,
        const int sig_frequencyextension_init);

    const std::vector<std::complex<double>> Gkpq(const IndexBK indexq, const int Nw_range) const;
    const std::vector<std::complex<double>> Gkpq(const IndexIBZ indexq, const int Nw_range) const;

    inline const std::complex<double>* data() const noexcept { return Gk.data(); }
    inline const std::complex<double>* data(IndexFK index) const noexcept
    {
        //check that frequency is inside box
        assert(std::abs(index.get_f()) < 2 * frequencyextension_ * Nw - 1);
        return &(Gk.data()[index.ind()]);
    }
    inline const std::complex<double> val(IndexFK index) const noexcept { return Gk[index.ind()]; }
};

#endif // GreenFuncion.h
