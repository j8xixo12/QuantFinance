#ifndef DOUBLE_SWEEP_HPP_
#define DOUBLE_SWEEP_HPP_
#include <memory>
#include <vector>

template<typename T, template<typename T, typename Alloc> class Container = std::vector, typename Alloc = std::allocator<T>>

class DoubleSweep {
    private:
        // The vectors of length J and start index 0
        Container<T, Alloc> a, b, c, f;
        // Dirichlet boundary conditions
        T left; // Left boundary condition
        T right; // Right boundary condition

        // Work arrays
        Container<T, Alloc> L;
        Container<T, Alloc> K;

    public:
        DoubleSweep() = delete;
        DoubleSweep(const DoubleSweep<T,Container, Alloc>& s2) = delete;

        DoubleSweep(const Container<T, Alloc>& LowerDiagonal,
                    const Container<T, Alloc>& Diagonal,
                    const Container<T, Alloc>& UpperDiagonal,
                    const Container<T, Alloc>& F,
                    const T& bc_left, const T& bc_right) {
            // Vectors are copied
            a = LowerDiagonal;
            b = Diagonal;
            c = UpperDiagonal;
            f = F;

            left = bc_left;
            right = bc_right;

            std::size_t N = a.size();

            // Work arrays
            L = Container<T, Alloc>(N, 0);
            K = Container<T, Alloc>(N, 0);
        }
        virtual ~DoubleSweep() = default;

        // Operator overloading
        DoubleSweep<T, Container, Alloc>& operator= (const DoubleSweep<T, Container>& i2) = delete;

        Container<T, Alloc> solve() {
            std::size_t N = a.size();

            // eq 13.7
            L[0] = 0.0;
            K[0] = left;

            // eq 13.6
            std::size_t SZ = L.size();

            for (std::size_t j = 1; j < SZ; ++j) {
                double tmp = b[j] + (a[j] * L[j - 1]);

                L[j] = -c[j] / tmp;
                K[j] = (f[j] - (a[j] * K[j - 1])) / tmp;
            }

            // Equation 13.5. Recycle array f for u
            f[0] = left;
            f[N - 1] = right;
            for (std::size_t j = f.size() - 2; j >= 1; --j) {
                f[j] = (L[j] * f[j + 1]) + K[j];
            }

            return f;
        }

        Container<T, Alloc> operator() () { return solve(); }
};

#endif // DOUBLE_SWEEP_HPP_