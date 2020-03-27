#ifndef LUTRIDIAGONALSOLVER_HPP_
#define LUTRIDIAGONALSOLVER_HPP_

#include <vector>
#include <cmath>

template<class T> using Vector = std::vector<T>;

template <class T> class LUTridiagonalSolver {
    private:
        // Defining arrays (input)
        // V2 optimise so to work with pointers
        Vector<T> a; // The lower-diagonal array [0..J]
        Vector<T> b; // The diagonal array [0..J] "baseline array"
        Vector<T> c; // The upper-diagonal array [0..J]
        Vector<T> r; // The right-hand side of the equation
        
        // Work arrays
        Vector<T> beta;
        Vector<T> gamma;
        Vector<T> z;

        std::size_t Size;

        void calculateBetaGamma_ZU(Vector<T>& r) {
            // eq 13.10
            beta[0] = b[0];
            gamma[0] = c[0] / beta[0];

            for (std::size_t i = 1; i < Size - 1; ++i) {
                beta[i] = b[i] - (a[i] * gamma[i - 1]);
                gamma[i] = c[i] / beta[i];
            }

            beta[Size - 1] = b[Size - 1] - (a[Size - 1] * gamma[Size - 2]);

            // Calculate z and u
            // Forward direction, equation 13.11
            z[0] = r[0] / beta[0];

            for (std::size_t j = 1; j < Size; ++j) {
                z[j] = (r[j] - (a[j] * z[j - 1])) / beta[j];
            }

            // Backward direction, equation 13.12
            r[Size - 1] = z[Size - 1];
            for (long i = Size - 2; i >= 0; --i) {
                r[i] = z[i] - (gamma[i] * r[i + 1]);
            }
        }

    public:
        LUTridiagonalSolver() = delete;
        LUTridiagonalSolver(const LUTridiagonalSolver<T>& source) = delete;
        virtual ~LUTridiagonalSolver() = default;
        LUTridiagonalSolver<T>& operator = (const LUTridiagonalSolver<T>& source) = delete;

        LUTridiagonalSolver(Vector<T>& lower, Vector<T>& diagonal,
                            Vector<T>& upper, Vector<T>& RHS) {
            a = lower;
            b = diagonal;
            c = upper;
            r = RHS;

            Size = diagonal.size();

            beta = Vector<T>(Size);
            gamma = Vector<T>(Size);
            z = Vector<T>(Size);
        }

        Vector<T> solve() {
            calculateBetaGamma_ZU(r);
            return r;
        }

        Vector<T> operator() () { return solve(); }

        bool diagonalDominant() const {
            if (std::abs(b[0]) < std::abs(c[0])) return false;
            if(std::abs(b[Size - 1]) < std::abs(a[Size - 1])) return false;

            for (std::size_t i = 1; i < Size - 1; ++i) {
                if (std::abs(b[i]) < std::abs(a[i]) + std::abs(c[i])) return false;
            }

            return true;
        }
};

#endif // LUTRIDIAGONALSOLVER_HPP_