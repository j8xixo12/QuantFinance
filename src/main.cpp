#include "Lattice.hpp"
#include "CRRLattice.hpp"

using namespace LatticeMechanism;

int main(int argc, char* argv[]) {
    Option opt(0.08, 0.3, 65.0, 0.25);

    constexpr double rootval = 60.0;
    constexpr int N = 3;
    double dt = opt.T_ / N;

    CRRLatticeAlgorithms algorithm(opt, dt);

    Lattice<double, 2> asset(N, 0.0);

    ForwardInduction<double, 2>(asset, algorithm, rootval);

    // Kinds of payoff
    double K = opt.K_;
    auto PutPayoff = [&K] (double S)-> double {return std::max<double> (K - S, 0.0);};
    
    // American early exercise constraint
    auto AmericanPutAdjuster = [&PutPayoff] (double& V, double S)->void { // e.g. early exercise
        V = std::max<double>(V, PutPayoff(S));
    };

    Lattice<double, 2> earlyPut(N, 0.0);
    ForwardInduction<double, 2>(earlyPut, algorithm, rootval);

    double earlyPutPrice = BackwardInduction<double, 2>(asset, earlyPut, algorithm, PutPayoff, AmericanPutAdjuster);
    std::cout << "Early Put: " << earlyPutPrice << std::endl;

    Lattice<std::tuple<double, double>, 2> combinedLattice = merge(asset, earlyPut);
    // print(combinedLattice);
    for (std::size_t i = 0; i < combinedLattice.Depth(); ++i) {
        for (std::size_t j = 0; j < combinedLattice[i].capacity(); ++j) {
            std::cout << std::get<0>(combinedLattice[i][j]) << ", " << std::get<1>(combinedLattice[i][j]) << " | ";
        }
        std::cout << std::endl;
    }
    return 0;
}
