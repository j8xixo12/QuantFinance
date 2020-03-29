#ifndef FDM_EULER_HPP_
#define FDM_EULER_HPP_
#include <cmath>
#include <random>
#include "Fdm.hpp"
#include "Sde.hpp"

template <typename T>
class FdmEuler : public Fdm<FdmEuler<T>, T> {
    public:
        FdmEuler(const std::shared_ptr<Sde<T>>& OneFactorProcess,
                const std::normal_distribution<T>& normalDist,
                const std::default_random_engine& engine)
        : Fdm<FdmEuler<T>, T>(OneFactorProcess, normalDist, engine) {}

    // Compute x(t_n + dt) in terms of x(t_n)
    T advance(T xn, T tn, T dt) {
        T normalVar = this->generateRN();
        return xn + this->sde->drift(xn, tn) * dt + this->sde->diffusion(xn, tn) * std::sqrt(dt) * normalVar;
    }
};

#endif // FDM_EULER_HPP_