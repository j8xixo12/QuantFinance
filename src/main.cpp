#include "Lattice.hpp"
#include "CRRLattice.hpp"

using namespace LatticeMechanism;

int main(int argc, char* argv[]) {
    Option opt(0.08, 0.3, 65.0, 0.25);

    constexpr double rootval = 65.0;
    constexpr int N = 100;
    double dt = opt.T_ / N;

    CRRLatticeAlgorithms algorithm(opt, dt);

    Lattice<double, 2> asset(N, 0.0);

    ForwardInduction<double, 2>(asset, algorithm, rootval);

    // Kinds of payoff
    double K = opt.K_;
    auto PutPayoff = [&K] (double S)-> double {return std::max<double> (K - S, 0.0);};
    auto CallPayoff = [&K] (double S)-> double {return std::max<double> (S - K, 0.0);};
    
    // American early exercise constraint
    auto AmericanPutAdjuster = [&PutPayoff] (double& V, double S)->void { // e.g. early exercise
        V = std::max<double>(V, PutPayoff(S));
    };
    auto AmericanCallAdjuster = [&CallPayoff] (double& V, double S)->void { // e.g. early exercise
        V = std::max<double>(V, CallPayoff(S));
    };

    Lattice<double, 2> euroPut(N, 0.0);
    ForwardInduction<double, 2>(euroPut, algorithm, rootval);

    Lattice<double, 2> euroCall(N, 0.0);
    ForwardInduction<double, 2>(euroCall, algorithm, rootval);
    Lattice<double, 2> earlyPut(N, 0.0);
    ForwardInduction<double, 2>(earlyPut, algorithm, rootval);
    Lattice<double, 2> earlyCall(N, 0.0);
    ForwardInduction<double, 2>(earlyCall, algorithm, rootval);

    double euroPutPrice = BackwardInduction<double, 2>(asset, euroPut, algorithm, PutPayoff);
    std::cout << "Euro Put: " << euroPutPrice << std::endl;

    double euroCallPrice = BackwardInduction<double, 2>(asset, euroCall, algorithm, CallPayoff);
    std::cout << "Euro Call: " << euroCallPrice << std::endl;

    double earlyPutPrice = BackwardInduction<double, 2>(asset, earlyPut, algorithm, PutPayoff, AmericanPutAdjuster);
    std::cout << "Early Put: " << earlyPutPrice << std::endl;

    double earlyCallPrice = BackwardInduction<double, 2>(asset, earlyCall, algorithm, CallPayoff, AmericanCallAdjuster);
    std::cout << "Early Call: " << earlyCallPrice << std::endl;
    
    return 0;
}
