#ifndef FDM_HPP_
#define FDM_HPP_

#include <random>
#include <memory>
#include "Sde.hpp"

template<typename T> class Fdm {
    protected:
        std::shared_ptr<Sde<T>> sde;
        int NT;

    public:
        std::vector<T> x;
        T dt;
        T dtSqrt;

        Fdm() {}
        virtual ~Fdm() {};
        Fdm(const std::shared_ptr<Sde<T>>& oneFactorProcess)
        : sde(oneFactorProcess), NT(0), dt(0.0), dtSqrt(0.0) {}

        // Compute x(t_n + dt) in terms of x(t_n)
        virtual T advance(T xn, T tn, T dt, T normalVar) const = 0;
};

#endif // FDM_HPP_