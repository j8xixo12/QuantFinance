#include <thread>
#include <iostream>
#include <chrono>
#include "CRRLattice.hpp"
#include "Lattice.hpp"
#include "Option.hpp"

using namespace LatticeMechanism;

auto elapsed(const std::chrono::steady_clock::time_point time1, const std::chrono::steady_clock::time_point time2) {
    return std::chrono::duration_cast<std::chrono::seconds>(time2 - time1).count();
}

int main(int argc, char* argv[]) {
    std::cout << std::thread::hardware_concurrency() << " processors/cores detected." <<'\n';
    Option opt;
    opt.K_ = 65.0;
    opt.T_ = 0.25;
    opt.r_ = 0.08;
    opt.sig_ = 0.3;

    int N = 7000;

    double dt = opt.T_ / static_cast<double>(N);
    std::cout << "Stress test, number of time steps: " << N << '\n';

    CRRLatticeAlgorithms algo(opt, dt);

    // Create basic asset lattice
    Lattice<double, 2> asset(N, 0.0); // init
    double rootVal = 60.0; // S(0) 
    ForwardInduction<double, 2>(asset, algo, rootVal);

    double K = opt.K_;
    auto PutPayoff = [&K] (double S)-> double { return std::max<double>(K - S, 0.0); }; 
    auto CallPayoff = [&K] (double S)-> double { return std::max<double>(S - K, 0.0); };

    // American early exercise constraint
    auto AmericanPutAdjuster = [&PutPayoff](double& V, double S)->void { // e.g. early exercise
        V = std::max<double>(V, PutPayoff(S));
    };
    auto AmericanCallAdjuster = [&CallPayoff](double& V, double S)->void { // e.g. early exercise
        V = std::max<double>(V, CallPayoff(S));
    };
    auto EmptyAdjuster = [](double& V, double S)->void { // No action executed at a node
        // Do nothing, no code generated in client code
    };

    std::chrono::steady_clock::time_point time1 = std::chrono::steady_clock::now();

    Lattice<double, 2> euroPut(N, 0.0); 
    ForwardInduction<double, 2>(euroPut, algo, rootVal); 
    Lattice<double, 2> euroCall(N, 0.0); 
    ForwardInduction<double, 2>(euroCall, algo, rootVal); 
    Lattice<double, 2> earlyPut(N, 0.0); 
    ForwardInduction<double, 2>(earlyPut, algo, rootVal); 
    Lattice<double, 2> earlyCall(N, 0.0); 
    ForwardInduction<double, 2>(earlyCall, algo, rootVal);

    double euroPutPrice = BackwardInduction<double>
            (asset, euroPut, algo, PutPayoff, EmptyAdjuster);
    std::cout << "Euro put: " << euroPutPrice << std::endl;
    double euroCallPrice = BackwardInduction<double>
              (asset, euroCall, algo, CallPayoff, EmptyAdjuster);
    std::cout << "Euro call: " << euroCallPrice << std::endl;

    // Price an early exercise option
    double earlyPutPrice = BackwardInduction<double>
           (asset, earlyPut, algo, PutPayoff,AmericanPutAdjuster);
    std::cout << "Early put: " << earlyPutPrice << std::endl;
    double earlyCallPrice = BackwardInduction<double>
           (asset, earlyCall, algo, CallPayoff,AmericanCallAdjuster);
    std::cout << "Early call: " << earlyCallPrice << std::endl;

    std::chrono::steady_clock::time_point time2 = std::chrono::steady_clock::now();
    std::cout << "Elapsed time sequential code: " << elapsed(time1, time2) << std::endl;

    double Euro_put, Euro_call, Early_put, Early_call;
    auto fn1 = [&]() {
        ForwardInduction<double> (euroPut, algo, rootVal);
        Euro_put = BackwardInduction<double> (asset, euroPut, algo, PutPayoff, EmptyAdjuster);
    };
    auto fn2 = [&]() {
        ForwardInduction<double> (euroCall, algo, rootVal);
        Euro_call = BackwardInduction<double> (asset, euroCall, algo, CallPayoff, EmptyAdjuster);
    };
    auto fn3 = [&]() {
        ForwardInduction<double> (earlyPut, algo, rootVal);
        Early_put = BackwardInduction<double> (asset, earlyPut, algo, PutPayoff, AmericanPutAdjuster);
    };
    auto fn4 = [&]() {
        ForwardInduction<double> (earlyCall, algo, rootVal);
        Early_call = BackwardInduction<double> (asset, earlyCall, algo, CallPayoff, AmericanCallAdjuster);
    };

    time1 = std::chrono::steady_clock::now();
    std::cout << '\n';
    std::thread t1(fn1);
    std::thread t2(fn2);
    std::thread t3(fn3);
    std::thread t4(fn4);
    // No shared data so we define 1 barrier/rendezvous 
    t1.join();
    t2.join();
    t3.join();
    t4.join();
    time2 = std::chrono::steady_clock::now();
    std::cout << "Euro put: " << Euro_put << std::endl;
    std::cout << "Euro call: " << Euro_call << std::endl;
    std::cout << "Early put: " << Early_put << std::endl;
    std::cout << "Early call: " << Early_call << std::endl;
    std::cout << "Elapsed time C++11 threads in seconds: " << elapsed(time1, time2) << std::endl;

    return 0;
}
