#ifndef IBVP_HPP_
#define IBVP_HPP_

#include <vector>
#include "Range.hpp"
#include "DoubleSweep.hpp"

class IBvpImp {
    // Abstract base class modelling the implementation objects in the
    // Bridge pattern. Derived classes (e.g. BS, CEV, CIR) must implement 
    // all the pure virtual functions in IBvpImp.
    public:
        // Selector functions
        // Coefficients of parabolic second order operator
        virtual double Diffusion(double x, double t) const = 0;
        virtual double Convection(double x, double t) const = 0;
        virtual double Reaction(double x, double t) const = 0;
        virtual double Rhs(double x, double t) const = 0;
        // Boundary and initial conditions
        virtual double Bcl(double t) const = 0;
        virtual double Bcr(double t) const = 0;
        virtual double Ic(double x) const = 0;
};

class IBvp {
    // Models linear convection diffusion reaction PDE with associated
    // boundary and initial conditions
    public:
        Range<double> xaxis;
        Range<double> taxis;
        // Bridge pattern: IBvp is the 'abstraction' object that is decoupled 
        // from its various implementations.
        IBvpImp* imp; // Bridge implementation

        IBvp() = delete;
        IBvp (IBvpImp& executor, const Range<double>& xrange,
          const Range<double>& trange) : xaxis(xrange), taxis(trange), imp(&executor) {};
        // Coefficients of parabolic second order operator
        // Coefficient of second derivative
        double Diffusion(double x, double t) const { return imp->Diffusion(x, t); }
        // Coefficient of first derivative
        double Convection(double x, double t) const { return imp->Convection(x, t); }
        // Coefficient of zero derivative
        double Reaction(double x, double t) const { return imp->Reaction(x, t); }
        // Inhomogeneous forcing term
        double Rhs(double x, double t) const { return imp->Rhs(x, t); }
        // Boundary and initial conditions
        // Left hand boundary condition
        double Bcl(double t) const { return imp->Bcl(t); }
        // Right hand boundary condition
        double Bcr(double t) const { return imp->Bcr(t); }
        // Initial condition
        double Ic(double x) const { return imp->Ic(x); }
        // The domain in which the PDE is 'played'
        const Range<double>& xrange() const { return xaxis; }
        const Range<double>& trange() const { return taxis; }
};

class IBvpSolver {
    private:
        // Utility functions
        void initMesh(long NSteps, long JSteps) {
            N = NSteps;
            J = JSteps;
            T = ibvp->trange().spread();
            k = T / static_cast<double>(N);
            h = ibvp ->xrange().spread() / static_cast<double>(J); 
            h2 = h * h;
            // Other numbers
            DN = static_cast<double>(N);
            DJ = static_cast<double>(J);
            DJJ = static_cast<double>(J * J);
            xarr = ibvp->xrange().mesh(J);
            // Array in t direction
            tarr = ibvp->trange().mesh(N); 
            tIndex = 0;
            vecOld = Vector(xarr.size(), 0.0); 
            vecNew = Vector(xarr.size(), 0.0);
        }
        void initIC() {
            // Utility function to initialise the payoff function
            // Initialise at the boundaries
            vecOld[0] = ibvp->Bcl(ibvp->trange().Low());
            vecOld[vecOld.size()-1] = ibvp->Bcr(ibvp->trange().High());
            // Initialise values in interior of interval using // the initial function 'IC' from the PDE
            for( std::size_t j = 1; j < xarr.size() - 1; ++j) {
                vecOld[ j ] = ibvp->Ic(xarr[j]); 
            }
        }
    protected:
        IBvp* ibvp;
        long N;
        double k;
        long J;
        double h, h2;
        double DN;
        double DJ;
        double DJJ;
        // Pointer to 'server'
        // The number of subdivisions of interval in IBVP
        // Step length; redundant data but is efficient
        // The number of subdivisions of interval in IBVP
        // Step length; redundant data but is efficient
        double tprev, tnow, T;
        long currentIndex, maxIndex, tIndex;
    public:
        using Vector = std::vector<double>;
        Vector xarr; // Useful work array
        Vector tarr; // Useful work array
        // Other data
        long n;  // Current counter
        Vector vecOld;
        Vector vecNew;

        IBvpSolver(const IBvpSolver& source) = delete;
        IBvpSolver& operator = (const IBvpSolver& source) = delete;
        IBvpSolver() {}
        IBvpSolver(IBvp& source, long NSteps, long JSteps) {
            ibvp = &source;
            // Create a mesh
            initMesh(NSteps, JSteps);
            // Initialise from the payoff function
            initIC();
        } 
        virtual ~IBvpSolver() {}
        virtual Vector& result() {
            // Initialise time.
            // The state machine, really; we march from t = 0 to t = T. 
            for(std::size_t n = 1; n < tarr.size(); ++n) {
                tnow = tarr[n];
                // The two methods that represent the variant parts // of the Template Method Pattern.
                calculate();
                tprev = tnow;
                for (std::size_t j = 0; j < vecNew.size(); ++j) { // Combine in previous loop
                    vecOld[j] = vecNew[j];
                }
            }
            return vecNew;
        }
        const Vector& XValues() const { return xarr; }
        const Vector& TValues() const { return tarr; }
        // The result of the calculation
                // Array of x values
                // Array of time values
        // Hook functions for Template Method pattern
        virtual void calculate() = 0; // Tells how to calculate sol.
        // at n+1
};

class CNIBVP : public IBvpSolver {
    private:
        // Notice that we store the data that 'belongs' to
        // this class. It is private and will not pollute the // other classes.
        Vector A, B, C; // Lower, diagonal, upper 
        Vector F; // Right-hand side of matrix
    public:
        CNIBVP(): IBvpSolver(), A(Vector(J + 1)), B(Vector(J + 1)),
                  C(Vector(J + 1)), F(Vector(J + 1)) {}
        CNIBVP(IBvp& source, long NSteps, long JSteps)
        : IBvpSolver(source, NSteps, JSteps),
        A(Vector(J + 1)), B(Vector(J + 1)), C(Vector(J + 1)),
        F(Vector(J + 1)) {}
        void calculate() {
            // Tells how to calculate sol. at n+1
            // In general we need to solve a tridiagonal system
            double t1, t2, t3, Low, Mid, Upp;
            for (std::size_t i = 0; i < F.size(); i++) {
                t1 = (0.5 * k * ibvp->Diffusion(xarr[i], tnow));
                t2 = 0.25 * k * h * ibvp->Convection(xarr[i], tnow); 
                t3 = 0.5 * k * h2 * ibvp->Reaction(xarr[i], tnow);
                // Coefficients of the U terms 
                A[i] = t1 - t2;
                B[i] = -h2 - 2.0 * t1 + t3; 
                C[i] = t1 + t2;
                // Coefficients of the U terms
                double t1A = 0.5 * k * ibvp->Diffusion(xarr[i], tnow);
                double t2A = 0.25 * k * h * ibvp->Convection(xarr[i], tprev);
                double t3A = 0.5 * k * h2 * ibvp->Reaction(xarr[i], tprev);
                Low = -t1A + t2A;
                Mid = -h2 + 2.0 * t1A - t3A; 
                Upp = -t1A - t2A;
                F[i] = Low * vecOld[i - 1] + Mid * vecOld[i] + Upp * vecOld[i + 1] + 0.5 * k * h2 * ibvp->Rhs(xarr[i], tnow)
                 + 0.5 * k * h2 * ibvp->Rhs(xarr[i], tprev);
            }
            double BCL = ibvp->Bcl(tnow);
            double BCR = ibvp->Bcr(tnow);
            DoubleSweep<double> mySolver(A, B, C, F, BCL, BCR);
            // The matrix must be diagonally dominant; we call the 
            // assert macro and the program stops
            // assert (mySolver.diagonallyDominant() == true);
            vecNew = mySolver.solve();
        }
};

template <typename D>
class ADE_CRTP: public IBvpSolver { 
    protected:
        // Derived classes contain the necessary data structures
    public:
        ADE_CRTP() : IBvpSolver() {};
        ADE_CRTP(IBvp& source, long NSteps, long JSteps) : IBvpSolver(source, NSteps, JSteps) {};
        // Hook function for Template Method pattern
        void calculate() {
            static_cast<D*>(this) -> calculate();
        }
};

class ADE_BC_CRTP: public ADE_CRTP<ADE_BC_CRTP> {
    private:
        // Intermediate values
        Vector U;
        Vector V;
        Vector UOld;
        Vector VOld;
    public:
        ADE_BC_CRTP(IBvp& source, long NSteps, long JSteps) : ADE_CRTP<ADE_BC_CRTP>(source, NSteps, JSteps),
        U(Vector(vecNew)), V(Vector(vecNew)), UOld(Vector(vecOld)),
        VOld(Vector(vecOld)) {}
        // Hook function for Template Method pattern
        void calculate() {
            // Tells how to calculate sol. at n+1, Explicit ADE_BC_CRTP schemes
            U[0] = ibvp->Bcl(tnow); 
            U[U.size()-1] = ibvp->Bcr(tnow);
            // Necessary?
            V[0] = ibvp->Bcl(tnow); 
            V[V.size()-1] = ibvp->Bcr(tnow);

            //std::cout << tnow << std::endl;
            for (std::size_t j = 1; j < U.size()-1; ++j) {
                  UOld[j] = VOld[j] = vecOld[j];
            }
            double t1, t2, t3, t4; 
            // double H = h / 2.0;
            // Upward sweep
            for (std::size_t j = 1; j < U.size()-1; ++j) {
                // Towler-Yang
                //t1 = k* fitting_factor(xarr[j], tnow)/ h2;
                t1 = k * ibvp->Diffusion(xarr[j], tnow) / h2;
                t2 = (0.5 * k * (ibvp->Convection(xarr[j], tnow))) / h; 
                t3 = (1.0 + t1 - ibvp->Reaction(xarr[j], tnow) * k); 
                t4 = -k * ibvp->Rhs(xarr[j], tnow);
                U[j] = ((t1 - t2) * U[j - 1] + (1.0 - t1) * UOld[j]
                                    + (t1 + t2) * UOld[j + 1] + t4) / t3;
            }
            // Downward sweep
            for (std::size_t j = V.size()-2; j >= 1; --j) {
                // Towler-Yang
                t1 = k * ibvp->Diffusion(xarr[j], tnow) / h2;
                t2 = (0.5 * k * (ibvp->Convection(xarr[j], tnow))) / h; 
                t3 = (1.0 + t1 - ibvp->Reaction(xarr[j], tnow) * k); 
                t4 = -k * ibvp->Rhs(xarr[j], tnow);
                V[j] = ((t1 - t2) * VOld[j - 1] + (1.0 - t1) * VOld[j] + (t1 + t2) * V[j + 1] + t4) / t3;
            }
            for (std::size_t j = 0; j < vecNew.size(); ++j) { 
                // Combine in previous loop
                vecNew[j] = 0.5 * (U[j] + V[j]); 
            }
        }
};
#endif // IBVP_HPP_
