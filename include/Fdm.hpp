#ifndef FDM_HPP_
#define FDM_HPP_

#include <random>
#include "Sde.hpp"

template<class Derived, typename T> class Fdm {
    protected:
        std::normal_distribution<T> dist;
        std::default_random_engine eng;
        std::shared_ptr<Sde<T>> sde;

    public:
        Fdm(const std::shared_ptr<Sde<T>>& oneFactorProcess,
            const std::normal_distribution<T>& normalDist,
            const std::default_random_engine& engine)
        : dist(normalDist), eng(engine), sde(oneFactorProcess) {}

        // Compute x(t_n + dt) in terms of x(t_n)
        T advance(T xn, double tn, double dt) {
            return static_cast<Derived*> (this)->advance(xn, tn, dt); 
        }

        T generateRN() {
            return dist(eng);
        }
};

#endif // FDM_HPP_