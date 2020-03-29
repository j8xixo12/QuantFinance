#ifndef SDE_FACTORY_HPP_
#define SDE_FACTORY_HPP_

#include "Sde.hpp"
#include "Option.hpp"
#include <memory>

enum SdeType { GBM, CEV, CIR };

template <class T = double>
class SdeFactory {
    public:
        static std::shared_ptr<Sde<T>> GetSde(const Option& data, T S0, SdeType type) {
            switch (type){
                SdeFactory<T> ret;
                case GBM:
                    return ret.GbmSde(data, S0);
                    break;
                
                case CEV:
                    return ret.CevSde(data, S0);
                    break;
                
                case CIR:
                {
                    // Feller condition => 2 a b >= sigmaˆ2 
                    double alpha = 0.01;
                    double b = 7.0;
                    double sigma = 0.9;
                    double B = 1.0; // Upper limit of interval 
                    double X0 = 0.4;
                    return ret.CirSde(alpha, b, sigma, X0, B);
                    break;
                }
                
                default:
                    return ret.GbmSde(data, S0);
                    break;
            }
        }
    
    private:
        std::shared_ptr<Sde<T>> GbmSde(const Option& data, T S0) {
            T r = data.r_; 
            T sig = data.sig_; 
            T B = data.T_;
            auto drift = [=](T t, T S) { return r * S; };
            auto diffusion = [=](T t, T S) { return sig * S; };
            ISde<T> Functions = std::make_tuple(drift, diffusion);
            // Create Sde
            return std::shared_ptr<Sde<T>> (std::make_shared<Sde<T>> (Functions, S0, B));
        }

        std::shared_ptr<Sde<T>> CevSde(const Option& data, T S0) {
            T r = data.r_; 
            T sig = data.sig_; 
            T B = data.T_;
            T beta = 0.5; // Hard-coded
            auto drift = [=](T t, T S) { return r * std::pow(S, beta); }; 
            auto diffusion = [=](T t, T S) { return sig * S; };
            ISde<T> Functions = std::make_tuple(drift, diffusion);
            // Create Sde
            return std::shared_ptr<Sde<T>>(std::make_shared<Sde<T>> (Functions, S0, B));
        }


        std::shared_ptr<Sde<T>> CirSde(T alpha, T b, T sigma, T S0, T B) {
            auto drift = [=](T t, T X) { return alpha*(b-X); };
            auto diffusion = [=](T t, T X) { return sigma * std::sqrt(X); }; 
            ISde<T> Functions = std::make_tuple(drift, diffusion);

            return std::shared_ptr<Sde<T>>(std::make_shared<Sde<T>> (Functions, S0, B));
        }
};

#endif // SDE_FACTORY_HPP_