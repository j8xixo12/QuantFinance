#ifndef FDM_HEUN_HPP_
#define FDM_HEUN_HPP_
#include <random>
#include "Fdm.hpp"

template <class T>
class FdmHeun : public Fdm<T> {
    public:
        FdmHeun() {}
        ~FdmHeun() {}
        FdmHeun(const std::shared_ptr<Sde<T>>& OneFactorProcess, int numSubdivisions)
        : Fdm<T>(OneFactorProcess) {
            this->NT = numSubdivisions;
            this->dt = this->sde->Expiry() / static_cast<T>(this->NT); 
            this->dtSqrt = std::sqrt(this->dt);
            this->x = std::vector<T>(this->NT + 1);
            // Create the mesh array
            this->x[0] = 0.0;
            for (std::size_t n = 1; n < this->x.size(); ++n) {
                this->x[n] = this->x[n - 1] + this->dt; 
            }
        }
    
        virtual T advance(T xn, T tn, T dt, T normalVar) const {
            T a = this->sde->drift(xn, tn);
            T b = this->sde->diffusion(xn, tn);
            T suppValue = xn + a * this->dt + b * this->dtSqrt * normalVar;
            return xn + 0.5 * (this->sde->drift(suppValue, tn) + a) * this->dt + 0.5 * (this->sde->diffusion(suppValue, tn) + b) * this->dtSqrt * normalVar;
        }
};
#endif // FDM_HEUN_HPP_