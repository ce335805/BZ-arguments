#include "Index.h"
#include "mpi.h"
#include <assert.h>

int IndexF::Nl;
int IndexB::Nl;

int IndexFK::Nx_;
int IndexFK::Ny_;
int IndexBK::Nx_;
int IndexBK::Ny_;
int IndexIBZ::Nx_;
int IndexIBZ::Ny_;
int IndexIBZ::NxIBZ_;

int IndexF::Nw;
int IndexB::Nw;
int IndexFK::Nw;
int IndexBK::Nw;
int IndexIBZ::Nw;

///////////////////////////////////
//definition of functions in IndexF
///////////////////////////////////

//no arguments - call fermionic start constructor
IndexF::IndexF()
    : IndexF(start())
{
}

//constructor with explicitly determined arguments
IndexF::IndexF(const int l_init, const int f_init)
    : Nw_internal_ { Nw }
{
    //assert that passsed frequency in uneven
    assert(std::abs(f_init % 2) == 1);
    assert(Nl > 0);
    assert(Nw > 0);
    this->l_ = l_init;
    this->f_ = f_init;
}
//start with first fermionic frequency
IndexF::IndexF(const start dummy_start)
    : Nw_internal_ { Nw }
{
    assert(Nl > 0);
    assert(Nw > 0);
    this->l_ = 0;
    this->f_ = -2 * Nw_internal_ + 1;
}
//returns the IndexF PAST the end
IndexF::IndexF(const end dummy_end)
    : Nw_internal_ { Nw }
{
    assert(Nl > 0);
    assert(Nw > 0);
    this->l_ = 0;
    this->f_ = 2 * Nw_internal_ + 1;
}

//no arguments - call fermionic start constructor
IndexF::IndexF(const int Nw_range)
    : IndexF(start(), Nw_range)
{
    assert(Nl > 0);
    assert(Nw > 0);
}

//constructor with explicitly determined arguments
IndexF::IndexF(const int l_init, const int f_init, const int Nw_range)
    : Nw_internal_ { Nw_range * Nw }
{
    //assert that passsed frequency in uneven
    assert(std::abs(f_init % 2) == 1);
    assert(Nl > 0);
    assert(Nw > 0);
    this->l_ = l_init;
    this->f_ = f_init;
}
//start with first fermionic frequency
IndexF::IndexF(const start dummy_start, const int Nw_range)
    : Nw_internal_ { Nw_range * Nw }
{
    assert(Nl > 0);
    assert(Nw > 0);
    this->l_ = 0;
    this->f_ = -2 * Nw_internal_ + 1;
}
//returns the IndexF PAST the end
IndexF::IndexF(const end dummy_end, const int Nw_range)
    : Nw_internal_ { Nw_range * Nw }
{
    assert(Nl > 0);
    assert(Nw > 0);
    this->l_ = 0;
    this->f_ = 2 * Nw_internal_ + 1;
}

//obvious(=default) copy constructor
IndexF::IndexF(const IndexF& to_copy)
    : Nw_internal_ { to_copy.Nw_internal_ }
{
    this->l_ = to_copy.get_l();
    this->f_ = to_copy.get_f();
}

//copy contructor with different number of frequencies
IndexF::IndexF(const IndexF& to_copy, const int Nw_range)
    : Nw_internal_ { Nw * Nw_range }
{
    this->l_ = to_copy.get_l();
    this->f_ = to_copy.get_f();
}

//addition/subtraction Fermion +- Fermion
//add only frequencies but keep first l-IndexF
IndexB IndexF::operator+(const IndexF& to_add)
{
    //assert that one adds indeces with the same frequency range
    assert(this->Nw_internal_ == to_add.get_Nw_int());
    return IndexB(this->l_, this->f_ + to_add.get_f(), int(this->Nw_internal_ / Nw));
}
//subtract only frequencies but keep first l-IndexF
IndexB IndexF::operator-(const IndexF& to_sub)
{
    //assert that one adds indeces with the same frequency range
    assert(this->Nw_internal_ == to_sub.get_Nw_int());
    return IndexB(this->l_, this->f_ - to_sub.get_f(), int(this->Nw_internal_ / Nw));
}

//addition/subtraction Fermion +- Boson
//add only frequencies but keep first l-IndexF
IndexF IndexF::operator+(const IndexB& to_add)
{
    //assert that one adds indeces with the same frequency range
    assert(this->Nw_internal_ == to_add.get_Nw_int());
    return IndexF(this->l_, this->f_ + to_add.get_f(), int(this->Nw_internal_ / Nw));
}
//add only frequencies but keep first l-IndexF
IndexF IndexF::operator-(const IndexB& to_sub)
{
    //assert that one adds indeces with the same frequency range
    assert(this->Nw_internal_ == to_sub.get_Nw_int());
    return IndexF(this->l_, this->f_ - to_sub.get_f(), int(this->Nw_internal_ / Nw));
}

//define output of indices
std::ostream& operator<<(std::ostream& os, const IndexF& idx)
{
    os << "idx(l, f, i) = (" << idx.get_l() << ", " << idx.get_f() << ", " << idx.ind() << ")";
    return os;
}

///////////////////////////////////
//definition of functions in IndexB
///////////////////////////////////

//no arguments - call fermionic start constructor
IndexB::IndexB()
    : IndexB(start())
{
}

//constructor with explicitl determined arguments
IndexB::IndexB(const int l_init, const int f_init)
    : Nw_internal_ { Nw }
{
    //check that passed frequency is even
    assert(f_init % 2 == 0);
    assert(Nl > 0);
    assert(Nw > 0);
    this->l_ = l_init;
    this->f_ = f_init;
}
//start with first fermionic frequency
IndexB::IndexB(const start dummy_start)
    : Nw_internal_ { Nw }
{
    assert(Nl > 0);
    assert(Nw > 0);
    this->l_ = 0;
    this->f_ = 0;
}
//returns the IndexB PAST the end
IndexB::IndexB(const end dummy_end)
    : Nw_internal_ { Nw }
{
    assert(Nl > 0);
    assert(Nw > 0);
    this->l_ = 0;
    this->f_ = 2 * Nw_internal_;
}

//no arguments - call fermionic start constructor
IndexB::IndexB(const int Nw_range)
    : IndexB(start(), Nw_range)
{
    assert(Nl > 0);
    assert(Nw > 0);
}

//constructor with explicitl determined arguments
IndexB::IndexB(const int l_init, const int f_init, const int Nw_range)
    : Nw_internal_ { Nw_range * Nw }
{
    //check that passed frequency is even
    assert(f_init % 2 == 0);
    assert(Nl > 0);
    assert(Nw > 0);
    this->l_ = l_init;
    this->f_ = f_init;
}
//start with first fermionic frequency
IndexB::IndexB(const start dummy_start, const int Nw_range)
    : Nw_internal_ { Nw_range * Nw }
{
    assert(Nl > 0);
    assert(Nw > 0);
    this->l_ = 0;
    this->f_ = 0;
}
//returns the IndexB PAST the end
IndexB::IndexB(const end dummy_end, const int Nw_range)
    : Nw_internal_ { Nw_range * Nw }
{
    assert(Nl > 0);
    assert(Nw > 0);
    this->l_ = 0;
    this->f_ = 2 * Nw_internal_;
}

//obvious(=default) copy constructor
IndexB::IndexB(const IndexB& to_copy)
    : Nw_internal_ { to_copy.Nw_internal_ }
{
    this->l_ = to_copy.get_l();
    this->f_ = to_copy.get_f();
}

//copy contructor with different number of frequencies
IndexB::IndexB(const IndexB& to_copy, const int Nw_range)
    : Nw_internal_ { Nw * Nw_range }
{
    this->l_ = to_copy.get_l();
    this->f_ = to_copy.get_f();
}

//addition/subtraction Boson +- Boson
//add only frequencies but keep first l-IndexF
IndexB IndexB::operator+(const IndexB& to_add)
{
    //assert that one adds indeces with the same frequency range
    assert(this->Nw_internal_ == to_add.get_Nw_int());
    return IndexB(this->l_, this->f_ + to_add.f_, int(this->Nw_internal_ / Nw));
}
//subtract only frequencies but keep first l-IndexF
IndexB IndexB::operator-(const IndexB& to_sub)
{
    //assert that one adds indeces with the same frequency range
    assert(this->Nw_internal_ == to_sub.get_Nw_int());
    return IndexB(this->l_, this->f_ - to_sub.f_, int(this->Nw_internal_ / Nw));
}

//addition/subtraction Boson +- Fermion
//add only frequencies but keep first l-IndexF
IndexF IndexB::operator+(const IndexF& to_add)
{
    //assert that one adds indeces with the same frequency range
    assert(this->Nw_internal_ == to_add.get_Nw_int());
    return IndexF(this->l_, this->f_ + to_add.get_f(), int(this->Nw_internal_ / Nw));
}
//add only frequencies but keep first l-IndexF
IndexF IndexB::operator-(const IndexF& to_sub)
{
    //assert that one adds indeces with the same frequency range
    assert(this->Nw_internal_ == to_sub.get_Nw_int());
    return IndexF(this->l_, this->f_ - to_sub.get_f(), int(this->Nw_internal_ / Nw));
}

//define output of indices
std::ostream& operator<<(std::ostream& os, const IndexB& idx)
{
    os << "idx(l, f, i) = (" << idx.get_l() << ", " << idx.get_f() << ", " << idx.ind() << ")";
    return os;
}

///////////////////////////////////
//definition of functions in IndexFK
///////////////////////////////////

//no arguments - call fermionic start constructor
IndexFK::IndexFK()
    : IndexFK(start())
{
}

//constructor with explicitl determined arguments
IndexFK::IndexFK(const int kx_init, const int ky_init, const int f_init)
    : Nw_internal_ { Nw }
{
    //check that passed frequency is uneven
    assert(std::abs(f_init % 2) == 1);
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(Nw > 0);

    this->kx_ = kx_init % Nx_;
    this->ky_ = ky_init % Ny_;
    this->f_ = f_init;
}
//start with first fermionic frequency
IndexFK::IndexFK(const start dummy_start)
    : Nw_internal_ { Nw }
{
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(Nw > 0);

    this->kx_ = 0;
    this->ky_ = 0;
    this->f_ = -2 * Nw_internal_ + 1;
}
//returns the IndexFK PAST the end
IndexFK::IndexFK(const end dummy_end)
    : Nw_internal_ { Nw }
{
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(Nw > 0);

    this->kx_ = 0;
    this->ky_ = 0;
    this->f_ = 2 * Nw_internal_ + 1;
}

//no arguments - call fermionic start constructor
IndexFK::IndexFK(const int Nw_range)
    : IndexFK(start(), Nw_range)
{
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(Nw > 0);
}

//constructor with explicitl determined arguments
IndexFK::IndexFK(const int kx_init, const int ky_init, const int f_init, const int Nw_range)
    : Nw_internal_ { Nw_range * Nw }
{
    //check that passed frequency is uneven
    assert(std::abs(f_init % 2) == 1);
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(Nw > 0);

    this->kx_ = kx_init % Nx_;
    this->ky_ = ky_init % Ny_;
    this->f_ = f_init;
}
//start with first fermionic frequency
IndexFK::IndexFK(const start dummy_start, const int Nw_range)
    : Nw_internal_ { Nw_range * Nw }
{
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(Nw > 0);

    this->kx_ = 0;
    this->ky_ = 0;
    this->f_ = -2 * Nw_internal_ + 1;
}
//returns the IndexFK PAST the end
IndexFK::IndexFK(const end dummy_end, const int Nw_range)
    : Nw_internal_ { Nw_range * Nw }
{
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(Nw > 0);

    this->kx_ = 0;
    this->ky_ = 0;
    this->f_ = 2 * Nw_internal_ + 1;
}

//obvious(=default) copy constructor
IndexFK::IndexFK(const IndexFK& to_copy)
    : Nw_internal_ { to_copy.Nw_internal_ }
{
    this->kx_ = to_copy.get_kx();
    this->ky_ = to_copy.get_ky();
    this->f_ = to_copy.get_f();
}

//copy contructor with different number of frequencies
IndexFK::IndexFK(const IndexFK& to_copy, const int Nw_range)
    : Nw_internal_ { Nw * Nw_range }
{
    this->kx_ = to_copy.get_kx();
    this->ky_ = to_copy.get_ky();
    this->f_ = to_copy.get_f();
}

//addition/subtraction Fermion +- Fermion
IndexBK IndexFK::operator+(const IndexFK& to_add)
{
    return IndexBK(this->kx_ + to_add.get_kx(), this->ky_ + to_add.get_ky(), this->f_ + to_add.get_f());
}

IndexBK IndexFK::operator-(const IndexFK& to_sub)
{
    return IndexBK(this->kx_ - to_sub.get_kx() + Nx_, this->ky_ - to_sub.get_ky() + Ny_, this->f_ - to_sub.get_f());
}

//addition/subtraction Fermion +- Boson
IndexFK IndexFK::operator+(const IndexBK& to_add)
{
    return IndexFK(this->kx_ + to_add.get_kx(), this->ky_ + to_add.get_ky(), this->f_ + to_add.get_f());
}

IndexFK IndexFK::operator-(const IndexBK& to_sub)
{
    return IndexFK(this->kx_ - to_sub.get_kx() + Nx_, this->ky_ - to_sub.get_ky() + Ny_, this->f_ - to_sub.get_f());
}

//define output of indices
std::ostream& operator<<(std::ostream& os, const IndexFK& idx)
{
    os << "idx(kx, ky, f, i) = (" << idx.get_kx() << ", " << idx.get_ky() << ", " << idx.get_f() << ", " << idx.ind() << ")";
    return os;
}

///////////////////////////////////
//definition of functions in IndexBK
///////////////////////////////////

//no arguments - call fermionic start constructor
IndexBK::IndexBK()
    : IndexBK(start())
{
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(Nw > 0);
}

//constructor with explicitl determined arguments
IndexBK::IndexBK(const int kx_init, const int ky_init, const int f_init)
    : Nw_internal_ { Nw }
{
    //check that passed frequency is uneven
    assert(f_init % 2 == 0);
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(Nw > 0);

    this->kx_ = kx_init % Nx_;
    this->ky_ = ky_init % Ny_;
    this->f_ = f_init;
}
//start with first fermionic frequency
IndexBK::IndexBK(const start dummy_start)
    : Nw_internal_ { Nw }
{
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(Nw > 0);

    this->kx_ = 0;
    this->ky_ = 0;
    this->f_ = 0;
}
//returns the IndexBK PAST the end
IndexBK::IndexBK(const end dummy_end)
    : Nw_internal_ { Nw }
{
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(Nw > 0);

    this->kx_ = 0;
    this->ky_ = 0;
    this->f_ = 2 * Nw_internal_;
}

//no arguments - call fermionic start constructor
IndexBK::IndexBK(const int Nw_range)
    : IndexBK(start(), Nw_range)
{
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(Nw > 0);
}

//constructor with explicitl determined arguments
IndexBK::IndexBK(const int kx_init, const int ky_init, const int f_init, const int Nw_range)
    : Nw_internal_ { Nw_range * Nw }
{
    //check that passed frequency is uneven
    assert(f_init % 2 == 0);
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(Nw > 0);

    this->kx_ = kx_init % Nx_;
    this->ky_ = ky_init % Ny_;
    this->f_ = f_init;
}
//start with first fermionic frequency
IndexBK::IndexBK(const start dummy_start, const int Nw_range)
    : Nw_internal_ { Nw_range * Nw }
{
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(Nw > 0);

    this->kx_ = 0;
    this->ky_ = 0;
    this->f_ = 0;
}
//returns the IndexBK PAST the end
IndexBK::IndexBK(const end dummy_end, const int Nw_range)
    : Nw_internal_ { Nw_range * Nw }
{
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(Nw > 0);

    this->kx_ = 0;
    this->ky_ = 0;
    this->f_ = 2 * Nw_internal_;
}

//obvious(=default) copy constructor
IndexBK::IndexBK(const IndexBK& to_copy)
    : Nw_internal_ { to_copy.Nw_internal_ }
{
    this->kx_ = to_copy.get_kx();
    this->ky_ = to_copy.get_ky();
    this->f_ = to_copy.get_f();
}

//copy contructor with different number of frequencies
IndexBK::IndexBK(const IndexBK& to_copy, const int Nw_range)
    : Nw_internal_ { Nw * Nw_range }
{
    this->kx_ = to_copy.get_kx();
    this->ky_ = to_copy.get_ky();
    this->f_ = to_copy.get_f();
}

//addition/subtraction Boson +- Boson
IndexBK IndexBK::operator+(const IndexBK& to_add)
{
    return IndexBK(this->kx_ + to_add.get_kx(), this->ky_ + to_add.get_ky(), this->f_ + to_add.get_f());
}
//subtract only frequencies but keep first l-IndexBK
IndexBK IndexBK::operator-(const IndexBK& to_sub)
{
    return IndexBK(this->kx_ - to_sub.get_kx() + Nx_, this->ky_ - to_sub.get_ky() + Ny_, this->f_ - to_sub.get_f());
}

//addition/subtraction Boson +- Fermion
IndexFK IndexBK::operator+(const IndexFK& to_add)
{
    return IndexFK(this->kx_ + to_add.get_kx(), this->ky_ + to_add.get_ky(), this->f_ + to_add.get_f());
}

IndexFK IndexBK::operator-(const IndexFK& to_sub)
{
    return IndexFK(this->kx_ - to_sub.get_kx() + Nx_, this->ky_ - to_sub.get_ky() + Ny_, this->f_ - to_sub.get_f());
}

//define output of indices
std::ostream& operator<<(std::ostream& os, const IndexBK& idx)
{
    os << "idx(kx, ky, f, i) = (" << idx.get_kx() << ", " << idx.get_ky() << ", " << idx.get_f() << ", " << idx.ind() << ")";
    return os;
}

///////////////////////////////////
//definition of functions in IndexIBZ
///////////////////////////////////

//no arguments - call fermionic start constructor
IndexIBZ::IndexIBZ()
    : IndexIBZ(start())
{
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(NxIBZ_ > 0);
    assert(Nw > 0);
}

//constructor with explicitl determined arguments
IndexIBZ::IndexIBZ(const int kx_init, const int ky_init, const int f_init)
    : Nw_internal_ { Nw }
{
    //check that passed frequency is uneven
    assert(f_init % 2 == 0);
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(NxIBZ_ > 0);
    assert(Nw > 0);

    this->kx_ = kx_init % Nx_;
    this->ky_ = ky_init % Ny_;
    this->f_ = f_init;
}
//start with first fermionic frequency
IndexIBZ::IndexIBZ(const start dummy_start)
    : Nw_internal_ { Nw }
{
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(NxIBZ_ > 0);
    assert(Nw > 0);

    // get the correct index at own task
    int worldsize;
    int myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &worldsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    const int index_start = (Nw * NxIBZ_) / worldsize * myrank;
    this->f_ = 2 * index_start / NxIBZ_;
    const int index_momentum = index_start % NxIBZ_;

    kx_ = 0;
    ky_ = 0;

    //determine ky_ by looping through 'wedge' of k_points
    for (int pos = 0; pos < Nx_ / 2 + 1; ++pos) {

        this->ky_ += int(index_momentum > (pos * (pos + 1)) / 2);
    }

    //kx_ then is the 'rest'
    kx_ = index_momentum - (ky_ * (ky_ + 1)) / 2;
}
//returns the IndexIBZ PAST the end
IndexIBZ::IndexIBZ(const end dummy_end)
    : Nw_internal_ { Nw }
{
    assert(Nx_ > 0);
    assert(Ny_ > 0);
    assert(NxIBZ_ > 0);
    assert(Nw > 0);

    // get the correct index at own task
    int worldsize;
    int myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &worldsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    const int index_start = (Nw * NxIBZ_) / worldsize * (myrank + 1);
    this->f_ =  2 * index_start / NxIBZ_;
    const int index_momentum = index_start % NxIBZ_;

    kx_ = 0;
    ky_ = 0;

    //determine ky_ by looping through 'wedge' of k_points
    for (int pos = 0; pos < Nx_ / 2 + 1; ++pos) {

        this->ky_ += int(index_momentum > (pos * (pos + 1)) / 2);
    }

    //kx_ then is the 'rest'
    kx_ = index_momentum - (ky_ * (ky_ + 1)) / 2;
}

//obvious(=default) copy constructor
IndexIBZ::IndexIBZ(const IndexIBZ& to_copy)
    : Nw_internal_ { to_copy.Nw_internal_ }
{
    this->kx_ = to_copy.get_kx();
    this->ky_ = to_copy.get_ky();
    this->f_ = to_copy.get_f();
}

//addition/subtraction Boson +- Fermion
IndexFK IndexIBZ::operator+(const IndexFK& to_add) const noexcept
{
    return IndexFK(this->kx_ + to_add.get_kx(), this->ky_ + to_add.get_ky(), this->f_ + to_add.get_f());
}

IndexFK IndexIBZ::operator-(const IndexFK& to_sub) const noexcept
{
    return IndexFK(this->kx_ - to_sub.get_kx() + Nx_, this->ky_ - to_sub.get_ky() + Ny_, this->f_ - to_sub.get_f());
}

//define output of indices
std::ostream& operator<<(std::ostream& os, const IndexIBZ& idx)
{
    os << "idx(kx, ky, f, i) = (" << idx.get_kx() << ", " << idx.get_ky() << ", " << idx.get_f() << ", " << idx.ind() << ")";
    return os;
}
