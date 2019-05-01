// header file which declares the util functions
#ifndef INDEX_H //inlcude guard
#define INDEX_H

#include <assert.h>
#include <iostream>

//to differentiate which constructor to use
struct start {
};
struct end {
};
//forward decleration of classes
class IndexF;
class IndexB;
class IndexFK;
class IndexBK;

class IndexF {

    int l_;
    int f_;
    const int Nw_internal_;

public:
    static int Nl;
		static int Nw;


    //no arguments - call fermionic start constructor
    IndexF();
    //constructor with explicitl determined arguments
    IndexF(const int l_init, const int f_init);
    //start with first fermionic frequency
    IndexF(const start dummy_start);
    //returns the IndexF PAST the end
    IndexF(const end dummy_end);

    //no arguments - call fermionic start constructor
    IndexF(const int Nw_range);
    //constructor with explicitl determined arguments
    IndexF(const int l_init, const int f_init, const int Nw_range);
    //start with first fermionic frequency
    IndexF(const start dummy_start, const int Nw_range);
    //returns the IndexF PAST the end
    IndexF(const end dummy_end, const int Nw_range);

    //obvious(=default) copy constructor
    IndexF(const IndexF& to_copy);
    //copy constructor with changed frequency range
    IndexF(const IndexF& to_copy, const int Nw_range);

    //default destructor
    ~IndexF(void) {}

    //prefix operator - no postfix to avoid unnecessary copying
    inline IndexF& operator++() noexcept
    {
        //increments f only when l is at last value
        f_ += 2 * ((l_ + 1) / Nl);
        l_ = (l_ + 1) % Nl;
        return *this;
    }

    //equality operator
    inline bool operator==(const IndexF& to_comp) const noexcept { return (this->get_f() == to_comp.get_f() && this->get_l() == to_comp.get_l()); }

    //compares to indices based on actual arguments
    inline bool operator>(const IndexF& to_comp) const noexcept
    {
        return (this->f_ > to_comp.get_f()
            || (this->f_ == to_comp.get_f()) && this->l_ > to_comp.get_l());
    }
    inline bool operator<(const IndexF& to_comp) const noexcept
    {
        return (this->f_ < to_comp.get_f()
            || (this->f_ == to_comp.get_f()) && this->l_ < to_comp.get_l());
    }

    //add only frequencies but keep first l-IndexF
    IndexB operator+(const IndexF& to_add);
    //subtract only frequencies but keep first l-IndexF
    IndexB operator-(const IndexF& to_sub);

    //add only frequencies but keep first l-IndexF
    IndexF operator+(const IndexB& to_add);
    //subtract only frequencies but keep first l-IndexF
    IndexF operator-(const IndexB& to_sub);

    //Index for array Indexing
    inline int ind() const noexcept { return ((f_ - 1) / 2 + Nw_internal_) * Nl + l_; }
    //Index for narrower frequency box
    inline int ind(const int Nw_range) const noexcept
    {
        //does not make sense to 'enlarge' box
        assert(Nw_range * Nw <= Nw_internal_);
        //check that narrowing box does not yield an invalid frequency
        assert(std::abs(this->f_) <= 2.0 * Nw_range * Nw - 1);

        return ((f_ - 1) / 2 + Nw_range * Nw) * Nl + l_;
    }
    //getters
    inline int get_f() const noexcept { return f_; }
    inline int get_l() const noexcept { return l_; }
    inline int get_Nw_int() const noexcept { return Nw_internal_; }

    friend std::ostream& operator<<(std::ostream& os, const IndexF& idx);
};

///////////////////////
//Bosonic index class
//////////////////////

class IndexB {

    int l_;
    int f_;
    const int Nw_internal_;

public:
    static int Nl;
		static int Nw;

    //no arguments - call fermionic start constructor
    IndexB();
    //constructor with explicitl determined arguments
    IndexB(const int l_init, const int f_init);
    //start with first fermionic frequency
    IndexB(const start dummy_start);
    //returns the IndexB PAST the end
    IndexB(const end dummy_end);

    //no arguments - call fermionic start constructor
    IndexB(const int Nw_range);
    //constructor with explicitl determined arguments
    IndexB(const int l_init, const int f_init, const int Nw_range);
    //start with first fermionic frequency
    IndexB(const start dummy_start, const int Nw_range);
    //returns the IndexB PAST the end
    IndexB(const end dummy_end, const int Nw_range);

    //obvious(=default) copy constructor
    IndexB(const IndexB& to_copy);
    //copy constructor with changed frequency range
    IndexB(const IndexB& to_copy, const int Nw_range);

    //default destructor
    ~IndexB(void) {}

    //prefix operator - no postfix to avoid unnecessary copying
    inline IndexB& operator++() noexcept
    {
        //increments f only when l is at last value
        f_ += 2 * ((l_ + 1) / Nl);
        l_ = (l_ + 1) % Nl;
        return *this;
    }

    //equality operator
    inline bool operator==(const IndexB& to_comp) const noexcept { return (this->get_f() == to_comp.get_f() && this->get_l() == to_comp.get_l()); }

    //compares to indices based on actual arguments
    inline bool operator>(const IndexB& to_comp) const noexcept
    {
        return (this->f_ > to_comp.get_f()
            || (this->f_ == to_comp.get_f()) && this->l_ > to_comp.get_l());
    }
    inline bool operator<(const IndexB& to_comp) const noexcept
    {
        return (this->f_ < to_comp.get_f()
            || (this->f_ == to_comp.get_f()) && this->l_ < to_comp.get_l());
    }

    //add only frequencies but keep first l-IndexF
    IndexB operator+(const IndexB& to_add);
    //subtract only frequencies but keep first l-IndexF
    IndexB operator-(const IndexB& to_sub);

    //add only frequencies but keep first l-IndexF
    IndexF operator+(const IndexF& to_add);
    //subtract only frequencies but keep first l-IndexF
    IndexF operator-(const IndexF& to_sub);

    inline int ind() const noexcept { return (f_ / 2) * Nl + l_; }
    //narrowed down index - nothing todo in this case
    inline int ind(const int Nw_range) const noexcept
    {
        assert(this->f_ <= Nw_range * Nw);
        return (f_ / 2) * Nl + l_;
    }

    inline int get_f() const noexcept { return f_; }
    inline int get_l() const noexcept { return l_; }
    inline int get_Nw_int() const noexcept { return Nw_internal_; }

    friend std::ostream& operator<<(std::ostream& os, const IndexB& idx);
};

///////////////////////
//Fermionic index class in k-space
//////////////////////


class IndexFK {

    int kx_;
    int ky_;
    int f_;
    const int Nw_internal_;

public:
    static int Nl;
        static int Nx_;
        static int Ny_;
		static int Nw;

    //no arguments - call fermionic start constructor
    IndexFK();
    //constructor with explicitl determined arguments
    IndexFK(const int kx_init, const int ky_init, const int f_init);
    //start with first fermionic frequency
    IndexFK(const start dummy_start);
    //returns the IndexFK PAST the end
    IndexFK(const end dummy_end);

    //no arguments - call fermionic start constructor
    IndexFK(const int Nw_range);
    //constructor with explicitl determined arguments
    IndexFK(const int kx_init, const int ky_init, const int f_init, const int Nw_range);
    //start with first fermionic frequency
    IndexFK(const start dummy_start, const int Nw_range);
    //returns the IndexFK PAST the end
    IndexFK(const end dummy_end, const int Nw_range);

    //obvious(=default) copy constructor
    IndexFK(const IndexFK& to_copy);
    //copy constructor with changed frequency range
    IndexFK(const IndexFK& to_copy, const int Nw_range);

    //default destructor
    ~IndexFK(void) {}

    //prefix operator - no postfix to avoid unnecessary copying
    inline IndexFK& operator++() noexcept
    {
        //        //increments f only when ky and kx have reached 'largest' value
        //        f_ += 2 * ((ky_ * Nx_ + kx_ + 1) / (Nx_ * Ny_));
        //        ky_ += ((kx_ + 1) / Nx_);
        //        ky_ = ky_ % Ny_;
        //        kx_ = (kx_ + 1) % Nx_;

        f_ += 2 * static_cast<int>(kx_ == Nx_ - 1 && ky_ == Ny_ - 1);
        ky_ = (ky_ + static_cast<int>(kx_ == Nx_ - 1)) % Ny_;
        kx_ = ++kx_ % Nx_;

        return *this;
    }

    //compares to indices based on actual arguments
    inline bool operator>(const IndexFK& to_comp) const noexcept
    {
        return (this->f_ > to_comp.get_f()
            || (this->f_ == to_comp.get_f()) && this->ky_ > to_comp.get_ky()
            || (this->f_ == to_comp.get_f()) && this->ky_ == to_comp.get_ky() && this->kx_ > to_comp.get_kx());
    }
    inline bool operator<(const IndexFK& to_comp) const noexcept
    {
        return (this->f_ < to_comp.get_f()
            || (this->f_ == to_comp.get_f()) && this->ky_ < to_comp.get_ky()
            || (this->f_ == to_comp.get_f()) && this->ky_ == to_comp.get_ky() && this->kx_ < to_comp.get_kx());
    }

    //add only frequencies but keep first l-IndexFK
    IndexBK operator+(const IndexFK& to_add);
    //subtract only frequencies but keep first l-IndexFK
    IndexBK operator-(const IndexFK& to_sub);

    //add only frequencies but keep first l-IndexFK
    IndexFK operator+(const IndexBK& to_add);
    //subtract only frequencies but keep first l-IndexFK
    IndexFK operator-(const IndexBK& to_sub);

    //Index for array Indexing
    inline int ind() const noexcept { return ((f_ - 1) / 2 + Nw_internal_) * Nx_ * Ny_ + ky_ * Nx_ + kx_; }
    //narrowed down index
    inline int ind(const int Nw_range) const noexcept
    {
        assert(Nw_range * Nw <= Nw_internal_);
        assert(std::abs(f_) <= 2.0 * (Nw_range * Nw) - 1);
        return ((f_ - 1) / 2 + Nw_range * Nw) * Nx_ * Ny_ + ky_ * Nx_ + kx_;
    }

    inline int get_kx() const noexcept { return kx_; }
    inline int get_ky() const noexcept { return ky_; }
    inline int get_f() const noexcept { return f_; }
    inline int get_Nw_int() const noexcept { return Nw_internal_; }

    friend std::ostream& operator<<(std::ostream& os, const IndexFK& idx);
};

//////////////////////////
//bosonic index class in k-space
//////////////////////////

class IndexBK {

    int kx_;
    int ky_;
    int f_;
    const int Nw_internal_;

public:
    static int Nl;
        static int Nx_;
        static int Ny_;
		static int Nw;

    //no arguments - call fermionic start constructor
    IndexBK();
    //constructor with explicitl determined arguments
    IndexBK(const int kx_init, const int ky_init, const int f_init);
    //start with first fermionic frequency
    IndexBK(const start dummy_start);
    //returns the IndexBK PAST the end
    IndexBK(const end dummy_end);

    //no arguments - call fermionic start constructor
    IndexBK(const int Nw_range);
    //constructor with explicitl determined arguments
    IndexBK(const int kx_init, const int ky_init, const int f_init, const int Nw_range);
    //start with first fermionic frequency
    IndexBK(const start dummy_start, const int Nw_range);
    //returns the IndexBK PAST the end
    IndexBK(const end dummy_end, const int Nw_range);

    //obvious(=default) copy constructor
    IndexBK(const IndexBK& to_copy);
    //copy constructor with changed frequency range
    IndexBK(const IndexBK& to_copy, const int Nw_range);

    //default destructor
    ~IndexBK(void) {}

    //prefix operator - no postfix to avoid unnecessary copying
    inline IndexBK& operator++() noexcept
    {
        //increments f only when kx and ky are at 'largest value'
        f_ += 2 * ((ky_ * Nx_ + kx_ + 1) / (Nx_ * Ny_));
        ky_ += ((kx_ + 1) / Nx_);
        ky_ = ky_ % Ny_;
        kx_ = (kx_ + 1) % Nx_;

        return *this;
    }

    //compares to indices based on actual arguments
    inline bool operator>(const IndexBK& to_comp) const noexcept
    {
        return (this->f_ > to_comp.get_f()
            || (this->f_ == to_comp.get_f()) && this->ky_ > to_comp.get_ky()
            || (this->f_ == to_comp.get_f()) && this->ky_ == to_comp.get_ky() && this->kx_ > to_comp.get_kx());
    }
    inline bool operator<(const IndexBK& to_comp) const noexcept
    {
        return (this->f_ < to_comp.get_f()
            || (this->f_ == to_comp.get_f()) && this->ky_ < to_comp.get_ky()
            || (this->f_ == to_comp.get_f()) && this->ky_ == to_comp.get_ky() && this->kx_ < to_comp.get_kx());
    }

    //add only frequencies but keep first l-IndexBK
    IndexBK operator+(const IndexBK& to_add);
    //subtract only frequencies but keep first l-IndexBK
    IndexBK operator-(const IndexBK& to_sub);

    //add only frequencies but keep first l-IndexBK
    IndexFK operator+(const IndexFK& to_add);
    //subtract only frequencies but keep first l-IndexFK
    IndexFK operator-(const IndexFK& to_sub);

    //Index for array Indexing
    inline int ind() const noexcept { return (f_ / 2) * Nx_ * Ny_ + ky_ * Nx_ + kx_; }
    //narrowed index - nothing to do here
    inline int ind(const int Nw_range) const noexcept
    {
        //check that function is used as intended
        assert(Nw_range * Nw <= Nw_internal_);
        assert(f_ <= Nw_range * Nw);
        return (f_ / 2) * Nx_ * Ny_ + ky_ * Nx_ + kx_;
    }

    inline int get_kx() const noexcept { return kx_; }
    inline int get_ky() const noexcept { return ky_; }
    inline int get_f() const noexcept { return f_; }
    inline int get_Nw_int() const noexcept { return Nw_internal_; }

    friend std::ostream& operator<<(std::ostream& os, const IndexBK& idx);
};

////////////////////////////
//bosonic index class with k-argument in irreducible brillouin zone
////////////////////////////

class IndexIBZ {

    int kx_;
    int ky_;
    int f_;
    const int Nw_internal_;

public:
    static int Nl;
        static int Nx_;
        static int Ny_;
        static int NxIBZ_;
		static int Nw;

    //no arguments - call fermionic start constructor
    IndexIBZ();
    //constructor with explicitl determined arguments
    IndexIBZ(const int kx_init, const int ky_init, const int f_init);
    //start with first fermionic frequency
    IndexIBZ(const start dummy_start);
    //returns the IndexIBZ PAST the end
    IndexIBZ(const end dummy_end);

    //obvious(=default) copy constructor
    IndexIBZ(const IndexIBZ& to_copy);

    //default destructor
    ~IndexIBZ(void) {}

    //prefix operator - no postfix to avoid unnecessary copying
    inline IndexIBZ& operator++() noexcept
    {

        f_ += 2 * static_cast<int>(kx_ == Nx_ / 2 && ky_ == Ny_ / 2);
        ++kx_;
        kx_ = kx_ % (ky_ + 1);
        ky_ += static_cast<int>(kx_ == 0);
        ky_ = (ky_ % ((Ny_ / 2) + 1));

        return *this;
    }

    //compares to indices based on actual arguments
    inline bool operator>(const IndexIBZ& to_comp) const noexcept
    {
        return (this->f_ > to_comp.get_f()
            || ((this->f_ == to_comp.get_f()) && this->ky_ > to_comp.get_ky())
            || ((this->f_ == to_comp.get_f()) && this->ky_ == to_comp.get_ky() && this->kx_ > to_comp.get_kx()));
    }
    inline bool operator<(const IndexIBZ& to_comp) const noexcept
    {
        return (this->f_ < to_comp.get_f()
            || ((this->f_ == to_comp.get_f()) && this->ky_ < to_comp.get_ky())
            || ((this->f_ == to_comp.get_f()) && this->ky_ == to_comp.get_ky() && this->kx_ < to_comp.get_kx()));
    }

    //add only frequencies but keep first l-IndexIBZ
    IndexFK operator+(const IndexFK& to_add) const noexcept;
    //subtract only frequencies but keep first l-IndexFK
    IndexFK operator-(const IndexFK& to_sub) const noexcept;

    //Index for array Indexing
    inline int ind() const noexcept { return (f_ / 2) * NxIBZ_ + (ky_ * (ky_ + 1))/2 + kx_; }

    inline int get_kx() const noexcept { return kx_; }
    inline int get_ky() const noexcept { return ky_; }
    inline int get_f() const noexcept { return f_; }
    inline int get_Nw_int() const noexcept { return Nw_internal_; }

    friend std::ostream& operator<<(std::ostream& os, const IndexIBZ& idx);
};

#endif //INDEX_H
