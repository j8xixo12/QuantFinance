#ifndef SDE_HPP_
#define SDE_HPP_
#include <functional>
#include "Option.hpp"

// Functions of arity 2 (two input arguments)
template <typename T> using FunctionType = std::function<T (const T& arg1, const T& arg2)>;

// Interface to simulate any SDE with drift/diffusion
template <typename T> using ISde = std::tuple<FunctionType<T>, FunctionType<T>>;

template <typename T = double> class Sde {
    protected:
        FunctionType<T> dr_;
        FunctionType<T> diff_;
        T mu; // Drift
        T vol; // Constant volatility
        T d; // Constant dividend yield
        T ic; // Initial condition
        T exp; // Expiry
    public:
        Sde() {};
        Sde(const FunctionType<T>& drift, const FunctionType<T>& diffusion,
            const Option&opt, const T& initialcCondition)
            : dr_(drift), diff_(diffusion), mu(opt.r_), vol(opt.sig_), 
            d(opt.D_), ic(initialcCondition), exp(opt.T_) {}
        Sde(const Sde<T>& sde2)
            : dr_(sde2.dr_), diff_(sde2.diff_), mu(sde2.mu), vol(sde2.vol), 
            d(sde2.d), ic(sde2.ic), exp(sde2.exp) {}

        virtual ~Sde() {};

        T drift(const T& S, const T& t) const { return dr_(S, t); }
        T diffusion(const T& S, const T& t) const { return diff_(S, t); }

        T Expiry() { return exp; }
        T InitialCondition() { return ic; }
};

template<typename T = double> class GBMSDE : public Sde<T> {    
    public:
        GBMSDE() {};
        virtual ~GBMSDE() {};
        GBMSDE(T driftCoefficient, T diffusionCoefficient,
            T dividendYield, T initialCondition,
            T expiry) {
                this->mu = driftCoefficient;
                this->vol = diffusionCoefficient;
                this->d = dividendYield;
                this->ic = initialCondition;
                this->exp = expiry;

                auto dr = [&] (T x, T t) -> T { return (this->mu - this->d) * x;};
                auto diff = [&] (T x, T t) -> T { return this->vol * x;};

                FunctionType<T> dr_fp = dr;
                FunctionType<T> diff_fp = diff;

                this->dr_ = dr_fp;
                this->diff_ = diff_fp;
        }

        GBMSDE(const Option &opt, T initialCondition) {
            this->mu = opt.r_;
            this->vol = opt.sig_;
            this->d = opt.D_;
            this->ic = initialCondition;
            this->exp = opt.T_;
            
            auto dr = [&] (T x, T t) -> T { return (this->mu - this->d) * x;};
            auto diff = [&] (T x, T t) -> T { return this->vol * x;};

            FunctionType<T> dr_fp = dr;
            FunctionType<T> diff_fp = diff;

            this->dr_ = dr_fp;
            this->diff_ = diff_fp;
        }
};

template<typename T = double> class CEVSDE : public Sde<T> {    
    private:
        T Beta;
    public:
        CEVSDE() {};

        CEVSDE(const Option &opt, T initialCondition, T beta_in) {
            this->mu = opt.r_;
            this->vol = opt.sig_;
            this->d = opt.D_;
            this->ic = initialCondition;
            this->exp = opt.T_;
            this->Beta = beta_in;
            
            auto dr = [&](T x, T t) { return this->mu * std::pow(x, Beta); }; 
            auto diff = [&](T x, T t) { return this->vol * x; };

            FunctionType<T> dr_fp = dr;
            FunctionType<T> diff_fp = diff;

            this->dr_ = dr_fp;
            this->diff_ = diff_fp;
        }
};

#endif // SDE_HPP_