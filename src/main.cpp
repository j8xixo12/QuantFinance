#include <vector>
#include <boost/numeric/odeint.hpp>
#include <fstream>
#include "BlackScholes.hpp"

using value_type = double;
typedef std::vector<value_type> state_type;

// Option Data
value_type K = 50.0;
value_type T = 0.5;
value_type r = 0.08;
value_type d = 0.0;
value_type sig = 0.4;

int main(int argc, char *argv[]) {

    namespace Bode = boost::numeric::odeint;
    // Input data
    const value_type T_0 = 0.0; 
    const value_type dt = 0.1;
    const value_type A = 0.0;
    value_type B = 1. * K;
    const int NX = 50.0;
    value_type h = (B - A) / static_cast<value_type>(NX);

    // Functions
    auto diffusion =  [&] (value_type x, value_type t) -> value_type { return 0.5 * sig * sig * x * x; }; 
    auto convection = [&] (value_type x, value_type t) -> value_type { return (r - d) * x; }; 
    auto reaction = [&] (value_type x, value_type t) -> value_type { return -r; };

    // BC
    auto bcl = [&](value_type t) -> value_type { return K * std::exp(-(r - d) * t); };
    auto bcr = [&](value_type t) -> value_type { return 0.; };
    auto payoff = [&] (value_type x) -> value_type { return std::max<value_type>(K - x, 0.0); };

    // For American options
    auto penalty = [&] (value_type x, value_type u) -> value_type {
        // The rationale is that if correction < 0 the
        // constraint is satisfied and the penalty = 0.
        // Otherise, we have to include a non-zero penalty

        value_type correction = payoff(x) - u;

        // Constraint already satisfied 
        if (correction <= 0.0) return 0.0;

        // A. (Courant) quadratic barrier function 
        value_type lambda = 1e+12;
        value_type penalty = lambda * correction * correction; 
        return penalty;
    };

    // Initial condition; discretise IC
    state_type U(NX + 1); // [0,J]
    U[0] = bcl(0.0); 
    U[U.size() - 1] = bcr(0.0); // Compatibility
    value_type x = A + h;
    for (std::size_t j = 1; j < U.size() - 1; ++j) {
        U[j] = payoff(x); 
        x += h;
    }

    // Integration_class
    OdeBlackScholes ode(NX, h, diffusion, bcl, bcr, convection, reaction, penalty);
    Bode::bulirsch_stoer<state_type, value_type> myStepper;
    // Bode::runge_kutta_cash_karp54<state_type, value_type> myStepper;

    std::ofstream output("out.csv", std::ios::out);
    std::vector<double> mesh = ode.get_mesh();

    auto write_out = [&] (const state_type& U, const value_type t ) {
        if (std::fmod(t, 0.01) == 0) {
            for (std::size_t i = 0; i < mesh.size(); ++i) {
                output << t << " " << mesh[i] << " " << U[i] << std::endl;
            }
        }
    };

    std::size_t steps = Bode::integrate_adaptive (myStepper, ode, U, T_0, T, dt, write_out);
    // std::cout << "Steps: " << steps << '\n';
    // std::size_t steps = Bode::integrate(ode, U, T_0, T, dt, write_out); 
    std::cout << "Steps: " << steps << '\n';

    output.close();

    return 0;
}