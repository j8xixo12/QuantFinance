#include <iostream>
#include <random>
#include <list>
#include "RNGenerator.hpp"
#include "utility.hpp"

int main(int argc, char* argv[]) {
    // Generate list and vector of uniform <int> and <double>
    double A = 0; 
    double B = 1.0;
    RNGenerator<double, std::uniform_real_distribution> rng(A, B); 
    std::cout << rng() << '\n';
    std::cout << generateRN(A, B) << '\n';
    int a = 1; 
    int b = 6;
    RNGenerator<int, std::uniform_int_distribution> rng2(a,b); 
    std::cout << rng2() << '\n';
    // std::size_t N = 7;
    // auto vec = rng.RandomArray(N); 
    // print(vec);
    // std::size_t M = 30;
    // auto vec2 = rng2.RandomArray(M); 
    // print(vec2);
    // auto vec3 = rng.RandomArrayII(N); 
    // print(vec3);
    // // Randomise 1st N elements of a container N =4;
    // std::vector<double> v(8, 3.1); 
    // rng.Randomise(N, v);
    // print(v);
    return 0;
}