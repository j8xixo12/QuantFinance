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
    public:
        T ic;
        T B; // SDE on interval [0,B]
        Sde() {};
        Sde(const FunctionType<T>& drift, const FunctionType<T>& diffusion,
            const T& initialcCondition, const T& expiration)
            : dr_(drift), diff_(diffusion), ic(initialcCondition),
            B(expiration) {}
        Sde(const Sde<T>& sde2, const T& initialcCondition, const
        T& expiration)
            : dr_(sde2.dr_), diff_(sde2.diff_), ic(initialcCondition),
            B(expiration) {}
        Sde(const ISde<T>& functions, const T& initialcCondition, const T& expiration)
            : dr_(std::get<0>(functions)), diff_(std::get<1>(functions)),
                ic(initialcCondition), B(expiration) {}
        T drift(const T& S, const T& t) const { return dr_(S, t); }
        T diffusion(const T& S, const T& t) const { return diff_(S, t); }
};

template<typename T = double> class GBM : public Sde<T> {
    private:
        double mu; // Drift
        double vol; // Constant volatility
        double d; // Constant dividend yield
        double ic; // Initial condition
        double exp; // Expiry
    
    public:
        GBM() {};
        GBM(double driftCoefficient, double diffusionCoefficient,
            double dividendYield, double initialCondition,
            double expiry) {
                mu = driftCoefficient;
                vol = diffusionCoefficient;
                d = dividendYield;
                ic = initialCondition;
                exp = expiry;

                auto dr = [&] (double x, double t) -> double { return (mu - d) * x;};
                auto diff = [&] (double x, double t) -> double { return vol * x;};

                FunctionType<T> dr_fp = dr;
                FunctionType<T> diff_fp = diff;

                this->dr_ = dr_fp;
                this->diff_ = diff_fp;
        }

        double Expiry() { return exp; }
        double InitialCondition() { return ic; }

        GBM(Option &opt, double initialCondition) {
            mu = opt.r_;
            vol = opt.sig_;
            d = opt.D_;
            ic = initialCondition;
            exp = opt.T_;
            
            auto dr = [&] (double x, double t) -> double { return (mu - d) * x;};
            auto diff = [&] (double x, double t) -> double { return vol * x;};

            FunctionType<T> dr_fp = dr;
            FunctionType<T> diff_fp = diff;

            this->dr_ = dr_fp;
            this->diff_ = diff_fp;
        }
};

#endif // SDE_HPP_