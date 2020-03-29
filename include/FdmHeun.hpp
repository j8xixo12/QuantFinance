#ifndef FDM_HEUN_HPP_
#define FDM_HEUN_HPP_
#include <random>
#include "Fdm.hpp"

template <class T>
class FdmHeun : public Fdm<FdmHeun<T>, T> {
    public:
        FdmHeun(const std::shared_ptr<Sde<T>>& OneFactorProcess,
                const std::normal_distribution<T>& normalDist,
                const std::default_random_engine& engine)
        : Fdm<FdmHeun<T>, T>(OneFactorProcess, normalDist, engine) {}
    
    T advance(T xn, T tn, T dt) {
        T a = this->sde->drift(xn, tn);
        T b = this->sde->diffusion(xn, tn);
        T normalVar = this->generateRN();
        T suppValue = xn + a * dt + b * std::sqrt(dt) * normalVar;
        return xn + 0.5 * (this->sde->drift(suppValue, tn) + a) * dt + 0.5 * (this->sde->diffusion(suppValue, tn) + b) * std::sqrt(dt) * normalVar;
    }
};
#endif // FDM_HEUN_HPP_