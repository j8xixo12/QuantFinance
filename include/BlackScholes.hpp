#ifndef BLACKSCHOLES_HPP_
#define BLACKSCHOLES_HPP_

#include <algorithm>
#include "IBvp.hpp"
#include "Option.hpp"

class BlackScholesPde : public IBvpImp {
    public:
        Option opt;
        BlackScholesPde(const Option& option) : opt(option) {};
        double Diffusion(double x, double t) const { return 0.5 * opt.sig_ * opt.sig_ * x * x; }
        double Convection(double x, double t) const { return opt.b_ * x; }
        double Reaction(double x, double t) const { return -opt.r_; }
        double Rhs(double x, double t) const { return 0.0; }
        // Boundary and initial conditions
        double Bcl(double t) const { 
            if (opt.type == OptionType::Call) { 
                return 0.0;
            } else {
                return (opt.K_) * std::exp(- (opt.r_) * ((opt.T_) - t));
            }
        }
        double Bcr(double t) const {
            if (opt.type == OptionType::Call) {
                // Right boundary condition
                return opt.SMax_ - (opt.K_) * std::exp(- (opt.r_) * ((opt.T_) - t));
            } else {
                return 0.0; //P 
            }
        }
        double Ic(double x) const {
            if (opt.type == OptionType::Call) {
                return std::max<double>(x - opt.K_, 0.0);
            } else {
                return std::max<double>(opt.K_ - x, 0.0);
            }
        }
};


#endif // BLACKSCHOLES_HPP_