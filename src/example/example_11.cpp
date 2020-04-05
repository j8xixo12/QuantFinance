#include "utility.hpp"
#include "GoursatFdm.hpp"

double M(double a, double b, double rho) {
    using value_type = double;
    std::size_t NX = 50; 
    std::size_t NY = 50; 
    value_type AL = -8.0; 
    value_type BL = -8.0; 
    auto xarr = CreateMesh(NX, AL, a);
    auto yarr = CreateMesh(NY, BL, b);
    BVNFunction<value_type> bvn(rho);
    GoursatFdmExtrapolation<value_type> fdm(AL, BL, bvn, NX, NY);
    return fdm(a, b);
}


int main(int argc, char *argv[]) {
    // 2-asset correlation option (Haug 2007) 
    {
        double S1 = 52.0; 
        double S2 = 65.0; 
        double T = 0.5;
        double K1 = 50.0; 
        double K2 = 70.0; 
        double s1 = 0.2; 
        double s2 = 0.3;
        double r = 0.1;
        double b1 = 0.1; 
        double b2 = 0.1; 
        double rho = 0.75;
        // Other variables
        double y1 = (std::log(S1 / K1) + (b1 - s1 * s1 / 2.) * T) / (s1 * std::sqrt(T)); 
        double y2 = (std::log(S2 / K2) + (b2 - s2 * s2 / 2.) * T) / (s2 * std::sqrt(T));
        
        // 2-asset correlation options
        double call = S2 * std::exp((b2 - r) * T) * M(y2 + s2 * std::sqrt(T), 1.
        + rho * s2 * std::sqrt(T), rho)
        - K2 * std::exp(-r * T) * M(y2, y1, rho);

        std::cout << "Call correlation: " << call << std::endl;

        double put = K2 * std::exp(-r * T) * M(-y2, -y1, rho)
        - S2 * std::exp((b2 - r) * T) * M(-y2 - s2 * std::sqrt(T),
        -y1 - rho * s2 * std::sqrt(T), rho);

        std::cout << "Put correlation: " << put << std::endl;
    }
    {
        // double S1 = 100.0; 
        // double S2 = 105.0;
        // double T = 0.5;
        // double K = 98.0;
        // double s1 = 0.11; 
        // double s2 = 0.16;
        // double r = 0.05;
        // double b1 = 0; 
        // double b2 = 0; // Investigate b < 0 */
        // double rho = 0.63;
        double S1 = 85.0; 
        double S2 = 60.0; 
        double T = 2.0;
        double K = 100.0;
        double s1 = 0.40; 
        double s2 = 0.25; 
        double r = 0.08;
        double b1 = r; 
        double b2 = r; // Investigate b < 0 
        double rho = -0.70;
        // Other variables
        double s = std::sqrt(s1 * s1 + s2 * s2 - 2.0 * rho * s1 * s2); 
        double d = (std::log(S1 / S2) + (b1 - b2 + s * s / 2.) * T)
            / (s * std::sqrt(T));
        double y1 = (std::log(S1 / K) + (b1 + s1 * s1 / 2.) * T) / (s1 * sqrt(T));
        double y2 = (std::log(S2 / K) + (b2 + s2 * s2 / 2.) * T) / (s2 * sqrt(T));
        double rho1 = (s1 - rho * s2) / s; 
        double rho2 = (s2 - rho * s1) / s;
        // As a function of S1,S2,K,T
        double cMinK =
        S1 * std::exp((b1 - r) * T) * M(y1, -d, -rho1)
        + S2 * std::exp((b2 - r) * T) * M(y2, d - s * std::sqrt(T), -rho2)
        - K * std::exp(-r * T)
        * M(y1 - s1 * std::sqrt(T), y2 - s2 * std::sqrt(T), rho);
        std::cout << "Call min K: " << cMinK << std::endl;
        double cMaxK = S1 * std::exp((b1 - r) * T) * M(y1, d, rho1)
        + S2 * std::exp((b2 - r) * T) * M(y2, -d + s * std::sqrt(T), rho2) - K * std::exp(-r * T) * (1.0 - M(-y1 + s1 * std::sqrt(T), -y2 + s2 * std::sqrt(T), rho));
        std::cout << "Call max K: " << cMaxK << std::endl;
    }
    {
        double S1 = 105.0; 
        double S2 = 100.0; 
        double T = 0.5;
        double Q = 5.0;
        double s1 = 0.3; 
        double s2 = 0.2;
        double r = 0.08;
        double b1 = 0.06; 
        double b2 = 0.02;
        double rho = 0.75; 
        double K = 100.0;
        // Other variables
        double s = std::sqrt(s1 * s1 + s2 * s2 - 2.0 * rho * s1 * s2); 
        double y = (std::log(S1 / S2) + (b1 - b2 + s * s / 2.) * T)
            / (s * std::sqrt(T));
        double z1 = (std::log(S1 / K) + (b1 + s1 * s1 / 2.) * T)
            / (s1 * std::sqrt(T));
        double z2 = (std::log(S2 / K) + (b2 + s2 * s2 / 2.) * T)
            / (s2 * std::sqrt(T));
        double rho1 = (s1 - rho * s2) / s; 
        double rho2 = (s2 - rho * s1) / s;
        // As a function of S1,S2,K,T
        double cBest = Q * std::exp(-r * T) * (M(y, z1, -rho1) + M(-y, z2, -rho2)); 
        std::cout << "Call worst cash-or-nothing best: "
                << cBest << std::endl;
    }
    {
        double S1 = 100.0; 
        double S2 = 100.0; 
        double K1 = 110.0; 
        double K2 = 90.0; 
        double T = 0.5;
        double Q = 10.0;
        double s1 = 0.2; 
        double s2 = 0.25; 
        double r = 0.1;
        double b1 = 0.05; 
        double b2 = 0.06; 
        double rho = 0.5;
        
        double d11 = (std::log(S1 / K1) + (b1 - s1 * s1 / 2.) * T) / (s1 * std::sqrt(T));
        double d22 = (std::log(S2 / K2) + (b2 - s2 * s2 / 2.) * T) / (s2 * std::sqrt(T));

        
        double cType1 = Q * std::exp(-r * T) * M(d11, d22, rho); 
        double cType2 = Q * std::exp(-r * T) * M(-d11, -d22, rho);
        double cType3 = Q  *std::exp(-r * T) * M(d11, -d22, -rho); 
        double cType4 = Q * std::exp(-r * T) * M(-d11, d22, -rho);
        std::cout << "Two-asset cash-or-nothing types: " << std::endl; 
        std::cout << "Type 1: " << cType1 << std::endl;
        std::cout << "Type 2: " << cType2 << std::endl;
        std::cout << "Type 3: " << cType3 << std::endl;
        std::cout << "Type 4: " << cType4 << std::endl;
    }
    return 0;
}