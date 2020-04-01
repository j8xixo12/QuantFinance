#include <iostream>
#include "GoursatFdm.hpp"
#include "utility.hpp"

int main(int argc, char* argv[]) {

    using value_type = double;
    //using value_type = float;
    // Truncated lower boundaries
    value_type AL = -8.0; 
    value_type BL = -8.0;
    // Meshes
    std::size_t NX = 200; 
    std::size_t NY = 200; 
    value_type a = 6.0; 
    value_type b = 6.0; 
    auto xarr = CreateMesh(NX, AL, a);
    auto yarr = CreateMesh(NY, BL, b);
    // Construct FDM object
    value_type rho = -0.1995;
    BVNFunction<value_type> bvn(rho); 
    GoursatFdm<value_type> fdm(AL, BL, bvn, xarr, yarr);
    // Compute a value
    std::cout << "Value: " << std::setprecision(16) << fdm(a, b) << '\n';

    return 0;
}