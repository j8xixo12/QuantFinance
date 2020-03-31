#include <iostream>
#include "Distributions.hpp"

template< typename T, template<typename TT> class Dist>
void TestAnyDistribution(const Distribution<T, Dist>& dist) {
    T x = 1.1;
    T p = 0.05;
    std::cout << "pdf: " << dist.pdf(x) << '\n';
    std::cout << "cdf: " << dist.cdf(x) << '\n';
    std::cout << "mean: " << dist.mean() << '\n'; std::cout << "mode: " << dist.mode() << '\n'; std::cout << "median: " << dist.median() << '\n';
    std::cout << "standard dev: "<<dist.standard_deviation() << '\n'; std::cout << "variance: " << dist.variance() << '\n';
    std::cout << "hazard: " << dist.standard_deviation() << '\n'; std::cout << "cumulative hazard: " << dist.variance() << '\n';
    std::cout << "kurtosis: " << dist.kurtosis() << '\n'; std::cout << "skewness: " << dist.skewness() << '\n';
    std::cout << "quantile " << p << ", "<< dist.quantile(p) << '\n';
}

template <typename T>
    using DistNormal = boost::math::normal_distribution<T>;
template <typename T>
    using DistBinomial = boost::math::binomial_distribution<T>;

int main(int argc, char* argv[]) {
    try {
        Distribution<double, DistNormal> normalDist;
        TestAnyDistribution(normalDist);
        long N = 100; 
        double p = 0.5; 
        Distribution<double, DistBinomial> binDist(N,p); 
        TestAnyDistribution(binDist);
    } catch(std::exception& e) {
        std::cout << e.what() << '\n'; 
    }

    return 0;
}