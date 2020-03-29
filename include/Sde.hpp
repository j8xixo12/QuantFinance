#ifndef SDE_HPP_
#define SDE_HPP_

// Functions of arity 2 (two input arguments)
template <typename T> using FunctionType = std::function<T (const T& arg1, const T& arg2)>;

// Interface to simulate any SDE with drift/diffusion
template <typename T> using ISde = std::tuple<FunctionType<T>, FunctionType<T>>;

template <typename T = double> class Sde {
    private:
        FunctionType<T> dr_;
        FunctionType<T> diff_;
    public:
        T ic;
        T B; // SDE on interval [0,B]
        Sde() = delete;
        Sde(const FunctionType<T>& drift, const FunctionType<T>& diffusion,
            const T& initialcCondition, const T& expiration)
            : dr_(drift), diff_(diffusion), ic(initialcCondition),
            B(expiration) {}
        Sde(const Sde<T>& sde2, const T& initialcCondition, const
        T& expiration)
            : dr_(sde2.dr_), diff_(sde2.diff_), ic(initialcCondition),
            B(expiration) {}
        Sde(const ISde<T>& functions, const T& initialcCondition, const T& expiration)
            : dr_(std::get<0>(functions)), diff_(std::get<1>(functions)),
                ic(initialcCondition), B(expiration) {}
        T drift(const T& S, const T& t) const { return dr_(S, t); }
        T diffusion(const T& S, const T& t) const { return diff_(S, t); }
};

#endif // SDE_HPP_