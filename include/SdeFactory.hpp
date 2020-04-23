#ifndef SDE_FACTORY_HPP_
#define SDE_FACTORY_HPP_

#include "Sde.hpp"
#include "Option.hpp"
#include <memory>

enum SdeType { GBM, CEV};

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
                
                default:
                    return ret.GbmSde(data, S0);
                    break;
            }
        }
    
    private:
        std::shared_ptr<Sde<T>> GbmSde(const Option& data, T S0) {
            return std::shared_ptr<Sde<T>> (std::make_shared<GBMSDE<T>> (data, S0));
        }

        std::shared_ptr<Sde<T>> CevSde(const Option& data, T S0) {
            T beta = 1.0; // Hard-coded
            return std::shared_ptr<Sde<T>>(std::make_shared<CEVSDE<T>> (data, S0, beta));
        }
};

#endif // SDE_FACTORY_HPP_