#ifndef UTILITY_HPP_
#define UTILITY_HPP_

std::vector<double> CreateDiscreteFunction(const std::vector<double>& xarr,
                                const std::function<double (double)>& fun) {
    std::size_t nx = xarr.size();

    std::vector<double> Vec (nx, 0.0);

    for (std::size_t i = 0; i < nx; ++i) {
        Vec[i] = fun(xarr[i]);
    }
    return Vec;
}

std::vector<double> CreateMesh(std::size_t n, double a, double b) { // Create a mesh of size n+1 on closed interval [a,b]
    std::vector<double> x(n + 1); 
    x[0] = a; 
    x[x.size()-1] = b;
    double h = (b - a) / static_cast<double>(n); 
    for (std::size_t j = 1; j < x.size() - 1; ++j) {
        x[j] = x[j - 1] + h; 
    }
    return x;
}

#endif // UTILITY_HPP_