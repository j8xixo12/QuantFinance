#ifndef GOURSATFDM_HPP_
#define GOURSATFDM_HPP_
#include <vector>
#include <functional>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/math/constants/constants.hpp>

template <typename T>
using Function = std::function<T (const T& arg1, const T& arg2)>;

template <typename T>
class GoursatFdm {
    private:
        T AL;
        T BL;
        std::vector<T> xarr;
        std::vector<T> yarr;
        Function<T> func;
    public:
        GoursatFdm(T xLower, T yLower, const Function<T>& function,
                    const std::vector<T>& xMesh, const std::vector<T>& yMesh)
        : AL(xLower), BL(yLower), xarr(xMesh), yarr(yMesh), func(function) {}

        T operator () (T x, T y) const {
            std::size_t NX = xarr.size() - 1;
            std::size_t NY = yarr.size() - 1;
            T hx = std::abs(x - AL) / static_cast<T>(NX); 
            T hy = std::abs(y - BL) / static_cast<T>(NY);
            // Matrix to hold results, boundary conditions == 0 
            boost::numeric::ublas::matrix<T> Values(NX + 1, NY + 1, 0.0);
            T hxy = hx * hy;
            for (std::size_t j = 1; j < Values.size2(); ++j) {
                for (std::size_t i = 1; i < Values.size1(); ++i) {
                    Values(i, j)= Values(i, j - 1) + Values(i - 1, j) - Values(i - 1, j - 1) + hxy * (func(xarr[i - 1] + hx / 2, yarr[j - 1] + hy / 2));
                } 
            }
            return Values(Values.size1() - 1, Values.size2() - 1); 
        }
};

template<class T>
class BVNFunction {
    private:
        T rho_;
    public:
        BVNFunction(T rho) : rho_(rho) {};
        T operator() (const T& x, const T& y) const {
            T tmp = 1. - rho_ * rho_;
            T tmp2 = 1. / (2. * boost::math::constants::pi<T>() * std::sqrt(tmp));
            return tmp2 * std::exp(-(x * x - 2. * rho_ * x * y + y * y) / (2. * tmp));
        };
};

#endif // GOURSATFDM_HPP_