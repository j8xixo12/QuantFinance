#ifndef CRRLATTICEALGORITHM_HPP_
#define CRRLATTICEALGORITHM_HPP_

#include "Option.hpp"
#include <cmath>

class CRRLatticeAlgorithms {
    public:
        explicit CRRLatticeAlgorithms(Option& option, const double dt) {
            double s = option.sig_;
            double r = option.r_;
            double R1 = (r - 0.5 * s * s) * dt;
            double R2 = s * std::sqrt(dt);
            u = std::exp(R1 + R2);
            d = std::exp(R1 - R2);
            discounting = std::exp(- r*dt);
            p = 0.5;
        }

        std::tuple<double, double> operator() (const double val) {
            return std::make_tuple(u * val, d * val);
        }

        double operator() (const double upper, const double lower) {
            return discounting * (p * upper + (1.0 - p) * lower); 
        }

    private:
        double u;
        double d;
        double p;
        double discounting;
};

#endif // CRRLATTICEALGORITHM_HPP_