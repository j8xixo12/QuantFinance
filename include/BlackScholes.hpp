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

class OdeBlackScholes {
    private:
        using value_type = double;
        using Function = std::function<value_type (value_type, value_type)>;
        using bc_Function = std::function<value_type (value_type)>;

        // The type of container used to hold the state vector
        typedef std::vector<value_type> state_type;
        std::vector<value_type> mesh;
        double h;
        Function diffusion_, convection_, reaction_, penalty_;
        bc_Function bcl_, bcr_;
    public:
        OdeBlackScholes(std::size_t NX, double meshSize, Function diffusion, bc_Function bcl,
                         bc_Function bcr, Function convection, Function reaction, Function penalty)
                          :mesh(std::vector<value_type>(NX + 1, 0.0)), h(meshSize),
                          diffusion_(diffusion), convection_(convection),
                          reaction_(reaction), penalty_(penalty), bcl_(bcl), bcr_(bcr) {
            mesh[0] = 0.0;
            for (std::size_t j = 1; j < mesh.size(); ++j) {
                mesh[j] = mesh[j - 1] + h; 
            }
        }

        void operator() (const state_type& U, state_type& dUdt,
             const value_type t) {
            double xval, diff,con;
            value_type h2 = 1.0 / (h * h); 
            value_type hm1 = 1.0 / (2.0 * h);
            // Boundaries
            // Left boundary (j = 1)
            std::size_t index = 1;
            xval = mesh[index];
            dUdt[index] = diffusion_(xval,t) * h2 * (U[index + 1] - 2.0 * U[index] + bcl_(t))
                        + 0.5 * convection_(xval, t) * (U[index + 1] - bcl_(t)) / h
                        + reaction_(xval, t) * U[index]
                        + penalty_(mesh[index], U[index]);
            // Right boundary (j = J-1) 
            index = U.size() - 2;
            xval = mesh[index];
            dUdt[index]
            = diffusion_(xval, t) * h2 * (bcr_(t) - 2.0 * U[index] + U[index - 1]) + 0.5 * convection_(xval, t) * (bcr_(t) - U[index - 1]) / h
            + reaction_(xval, t) * U[index]
            + penalty_(mesh[index], U[index]);
            // Interior of domain (1 < j < J-1)
            for (std::size_t index = 2; index < U.size() - 2; ++index) {
                xval = mesh[index];
                diff = diffusion_(xval, t) * h2;
                con = convection_(xval, t) * hm1;
                dUdt[index] = (diff + con) * U[index + 1]
                + (-2.0 * diff + reaction_(xval, t)) * U[index] + (diff - con) * U[index - 1]
                + penalty_(mesh[index], U[index]);
            } 
        }

        void operator () (const state_type& U, const value_type t) {
            if (t >= 1.0) { // To avoid explosion of output at the moment
          // Your code here
            }
        }

        std::vector<value_type> get_mesh() { return mesh; }

};


#endif // BLACKSCHOLES_HPP_