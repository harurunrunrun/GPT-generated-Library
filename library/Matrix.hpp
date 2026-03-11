#include <bits/stdc++.h>
using namespace std;

template <class AddMonoid, class MulMonoid>
struct SemiringFromMonoids {
    using value_type = typename AddMonoid::value_type;
    static_assert(std::is_same_v<value_type, typename MulMonoid::value_type>);

    static value_type add(const value_type& a, const value_type& b) {
        return AddMonoid::op(a, b);
    }
    static value_type mul(const value_type& a, const value_type& b) {
        return MulMonoid::op(a, b);
    }
    static value_type zero() { return AddMonoid::id(); }
    static value_type one() { return MulMonoid::id(); }

#ifndef NDEBUG
    static vector<value_type> debug_samples() {
        return { zero(), one() };
    }
#endif
};

#ifndef NDEBUG
template <class SR>
struct SemiringValidator {
    using T = typename SR::value_type;

    static void run() {
        static const bool once = []() {
            const T z = SR::zero();
            const T o = SR::one();
            auto samples = SR::debug_samples();
            if (samples.empty()) samples = {z, o};

            for (const T& a : samples) {
                assert(SR::add(z, a) == a);
                assert(SR::add(a, z) == a);
                assert(SR::mul(o, a) == a);
                assert(SR::mul(a, o) == a);
                assert(SR::mul(z, a) == z);
                assert(SR::mul(a, z) == z);
            }
            return true;
        }();
        (void)once;
    }
};
#endif

template <class SR>
struct Matrix {
    using T = typename SR::value_type;

    int H, W;
    vector<T> dat;

    Matrix() : H(0), W(0) {
#ifndef NDEBUG
        SemiringValidator<SR>::run();
#endif
    }

    Matrix(int h, int w) : H(h), W(w), dat(h * w, SR::zero()) {
#ifndef NDEBUG
        SemiringValidator<SR>::run();
#endif
    }

    Matrix(int h, int w, const T& v) : H(h), W(w), dat(h * w, v) {
#ifndef NDEBUG
        SemiringValidator<SR>::run();
#endif
    }

    Matrix(const vector<vector<T>>& a) {
#ifndef NDEBUG
        SemiringValidator<SR>::run();
#endif
        H = (int)a.size();
        W = H ? (int)a[0].size() : 0;
        dat.assign(H * W, SR::zero());
        for (int i = 0; i < H; ++i) {
            assert((int)a[i].size() == W);
            for (int j = 0; j < W; ++j) {
                (*this)(i, j) = a[i][j];
            }
        }
    }

    static Matrix identity(int n) {
        Matrix I(n, n, SR::zero());
        for (int i = 0; i < n; ++i) I(i, i) = SR::one();
        return I;
    }

    T& operator()(int i, int j) { return dat[i * W + j]; }
    const T& operator()(int i, int j) const { return dat[i * W + j]; }

    Matrix& operator+=(const Matrix& rhs) {
        assert(H == rhs.H && W == rhs.W);
        for (int i = 0; i < H * W; ++i) dat[i] = SR::add(dat[i], rhs.dat[i]);
        return *this;
    }

    Matrix operator+(const Matrix& rhs) const {
        Matrix res = *this;
        res += rhs;
        return res;
    }

    Matrix operator*(const Matrix& rhs) const {
        assert(W == rhs.H);
        Matrix res(H, rhs.W, SR::zero());
        for (int i = 0; i < H; ++i) {
            for (int k = 0; k < W; ++k) {
                const T& aik = (*this)(i, k);
                for (int j = 0; j < rhs.W; ++j) {
                    res(i, j) = SR::add(res(i, j), SR::mul(aik, rhs(k, j)));
                }
            }
        }
        return res;
    }

    Matrix& operator*=(const Matrix& rhs) {
        return *this = (*this) * rhs;
    }

    Matrix pow(long long e) const {
        assert(H == W);
        assert(e >= 0);
        Matrix base = *this;
        Matrix res = identity(H);
        while (e > 0) {
            if (e & 1) res *= base;
            base *= base;
            e >>= 1;
        }
        return res;
    }

    vector<T> operator*(const vector<T>& x) const {
        assert(W == (int)x.size());
        vector<T> y(H, SR::zero());
        for (int i = 0; i < H; ++i) {
            T acc = SR::zero();
            for (int j = 0; j < W; ++j) {
                acc = SR::add(acc, SR::mul((*this)(i, j), x[j]));
            }
            y[i] = acc;
        }
        return y;
    }
};

// ---------- basic monoids ----------

template <class T>
struct AddMonoid {
    using value_type = T;
    static T op(const T& a, const T& b) { return a + b; }
    static T id() { return T(0); }
};

template <class T>
struct MulMonoid {
    using value_type = T;
    static T op(const T& a, const T& b) { return a * b; }
    static T id() { return T(1); }
};

template <class T, T INF>
struct MinMonoid {
    using value_type = T;
    static T op(const T& a, const T& b) { return min(a, b); }
    static T id() { return INF; }
};

template <class T, T NEG_INF>
struct MaxMonoid {
    using value_type = T;
    static T op(const T& a, const T& b) { return max(a, b); }
    static T id() { return NEG_INF; }
};

// ---------- corrected min-plus / max-plus mul monoids ----------

template <class T, T INF>
struct MinPlusMulMonoid {
    using value_type = T;
    static T op(const T& a, const T& b) {
        if (a == INF || b == INF) return INF;
        return a + b;
    }
    static T id() { return T(0); }
};

template <class T, T NEG_INF>
struct MaxPlusMulMonoid {
    using value_type = T;
    static T op(const T& a, const T& b) {
        if (a == NEG_INF || b == NEG_INF) return NEG_INF;
        return a + b;
    }
    static T id() { return T(0); }
};

struct OrMonoid {
    using value_type = bool;
    static bool op(bool a, bool b) { return a || b; }
    static bool id() { return false; }
};

struct AndMonoid {
    using value_type = bool;
    static bool op(bool a, bool b) { return a && b; }
    static bool id() { return true; }
};

// ---------- semirings ----------

using LongLongSemiring =
    SemiringFromMonoids<AddMonoid<long long>, MulMonoid<long long>>;

constexpr long long INF64 = (1LL << 60);
using MinPlusSemiring =
    SemiringFromMonoids<MinMonoid<long long, INF64>, MinPlusMulMonoid<long long, INF64>>;

constexpr long long NEG_INF64 = -(1LL << 60);
using MaxPlusSemiring =
    SemiringFromMonoids<MaxMonoid<long long, NEG_INF64>, MaxPlusMulMonoid<long long, NEG_INF64>>;

using BoolSemiring =
    SemiringFromMonoids<OrMonoid, AndMonoid>;

// ---------- optional: richer debug samples ----------

template <>
vector<long long> MinPlusSemiring::debug_samples() {
    return { zero(), one(), 1, 5, 10 };
}

template <>
vector<long long> MaxPlusSemiring::debug_samples() {
    return { zero(), one(), -3, 0, 7 };
}
