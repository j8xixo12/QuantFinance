#ifndef FDM_EULER_HPP_
#define FDM_EULER_HPP_
#include <cmath>
#include <random>
#include "Fdm.hpp"
#include "Sde.hpp"

template <typename T>
class FdmEuler : public Fdm<T> {
    public:
        FdmEuler() {}
        ~FdmEuler() {}
        FdmEuler(const std::shared_ptr<Sde<T>>& OneFactorProcess, int numSubdivisions)
        : Fdm<T>(OneFactorProcess) {
            this->NT = numSubdivisions;
            this->dt = this->sde->Expiry() / static_cast<double>(this->NT); 
            this->dtSqrt = std::sqrt(this->dt);
            this->x = std::vector<double>(this->NT + 1);
            // Create the mesh array
            this->x[0] = 0.0;
            for (std::size_t n = 1; n < this->x.size(); ++n) {
                this->x[n] = this->x[n - 1] + this->dt; 
            }
        }

        // Compute x(t_n + dt) in terms of x(t_n)
        virtual T advance(T xn, T tn, T dt, T normalVar) const {
            return xn + this->sde->drift(xn, tn) * dt + this->sde->diffusion(xn, tn) * this->dtSqrt * normalVar;
        }
};


// template <typename T> class EulerFdm { 
//     private:
//         std::shared_ptr<Sde<T>> sde;
//         int NT;
//     public:
//         std::vector<double> x;
//         double dt;
//         double dtSqrt;
//     // The mesh array
//     // Mesh size
//     EulerFdm() = default;
//     EulerFdm(const std::shared_ptr<Sde<T>>& stochasticEquation, int numSubdivisions)
//     : sde(stochasticEquation), NT(numSubdivisions) {
//         k = sde->Expiry() / static_cast<double>(NT); 
//         dtSqrt = std::sqrt(k);
//         x = std::vector<double>(NT + 1);
//         // Create the mesh array
//         x[0] = 0.0;
//         for (std::size_t n = 1; n < x.size(); ++n) {
//             x[n] = x[n - 1] + k; 
//         }
//     }
//     double advance(double xn, double tn, double dt,
//                     double  normalVar) const {
//         return xn + sde->drift(xn, tn) * dt
//         + sde->diffusion(xn, tn) * dtSqrt * normalVar;
//     } 
// };

#endif // FDM_EULER_HPP_