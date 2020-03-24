#include "Lattice.hpp"
#include "CRRLattice.hpp"

int main(int argc, char* argv[]) {
    Option opt(0.08, 0.3, 65.0, 0.25);

    int N = 101;
    double dt = opt.T_ / N;

    CRRLatticeAlgorithms algorithm(opt, dt);

    LatticeMechanism::Lattice<double, 2> lattice(N, 0.0);

    double rootval = 60.0;

    double K = opt.K_;
    auto PayOff = [&K](const double S) -> double { return std::max<double>(K - S, 0.0); };

    LatticeMechanism::ForwardInduction<double, 2>(lattice, algorithm, rootval);

    double res = LatticeMechanism::BackwardInduction<double, 2>(lattice, algorithm, PayOff);

    std::cout << "Plain Option price, classic #1: ";
    std::cout << '\n' << res << std::endl;
    return 0;
}
