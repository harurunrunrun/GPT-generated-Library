#pragma once
#include <bits/stdc++.h>
using namespace std;

namespace cp_counting {

// ============================================================================
// ModInt
// ============================================================================

// MOD を法とする整数を扱う構造体。
// 四則演算、累乗、逆元を提供する。
// 想定: MOD は素数。
template <int MOD>
struct ModInt {
    int v;

    // 値 x を [0, MOD) に正規化して保持する。
    // 計算量: O(1)
    ModInt(long long x = 0) {
        x %= MOD;
        if (x < 0) x += MOD;
        v = (int)x;
    }

    // 法 MOD を返す。
    // 計算量: O(1)
    static constexpr int mod() { return MOD; }

    // 加算代入。
    // 計算量: O(1)
    ModInt& operator+=(const ModInt& other) {
        v += other.v;
        if (v >= MOD) v -= MOD;
        return *this;
    }

    // 減算代入。
    // 計算量: O(1)
    ModInt& operator-=(const ModInt& other) {
        v -= other.v;
        if (v < 0) v += MOD;
        return *this;
    }

    // 乗算代入。
    // 計算量: O(1)
    ModInt& operator*=(const ModInt& other) {
        v = (long long)v * other.v % MOD;
        return *this;
    }

    // 除算代入。逆元を用いる。
    // 計算量: O(log MOD)
    ModInt& operator/=(const ModInt& other) {
        return *this *= other.inv();
    }

    // 加算。
    // 計算量: O(1)
    friend ModInt operator+(ModInt a, const ModInt& b) { return a += b; }

    // 減算。
    // 計算量: O(1)
    friend ModInt operator-(ModInt a, const ModInt& b) { return a -= b; }

    // 乗算。
    // 計算量: O(1)
    friend ModInt operator*(ModInt a, const ModInt& b) { return a *= b; }

    // 除算。
    // 計算量: O(log MOD)
    friend ModInt operator/(ModInt a, const ModInt& b) { return a /= b; }

    // 等値判定。
    // 計算量: O(1)
    friend bool operator==(const ModInt& a, const ModInt& b) { return a.v == b.v; }

    // 非等値判定。
    // 計算量: O(1)
    friend bool operator!=(const ModInt& a, const ModInt& b) { return a.v != b.v; }

    // x^e を二分累乗法で求める。
    // 計算量: O(log e)
    ModInt pow(long long e) const {
        ModInt x = *this, r = 1;
        while (e > 0) {
            if (e & 1) r *= x;
            x *= x;
            e >>= 1;
        }
        return r;
    }

    // 逆元を返す。
    // Fermat の小定理 a^(MOD-2) を用いる。
    // 前提: MOD は素数、v != 0。
    // 計算量: O(log MOD)
    ModInt inv() const {
        assert(v != 0);
        return pow(MOD - 2);
    }

    // 出力演算子。
    // 計算量: O(1)
    friend ostream& operator<<(ostream& os, const ModInt& x) {
        return os << x.v;
    }

    // 入力演算子。
    // 計算量: O(1)
    friend istream& operator>>(istream& is, ModInt& x) {
        long long t;
        is >> t;
        x = ModInt(t);
        return is;
    }
};

// ============================================================================
// Comb
// ============================================================================

// 階乗 fact、逆階乗 ifact、逆元 inv を前計算して
// 組合せや順列を高速に求めるための構造体。
// 想定: MOD は素数、かつ扱う n は n < MOD。
template <int MOD>
struct Comb {
    using mint = ModInt<MOD>;
    vector<mint> fact, ifact, inv;

    // 長さ n まで前計算を行うコンストラクタ。
    // 計算量: O(n)
    Comb(int n = 0) {
        fact = {1};
        ifact = {1};
        inv = {0};
        if (n > 0) ensure(n);
    }

    // 必要に応じて fact, ifact, inv を n まで拡張する。
    // すでに十分なら何もしない。
    // 計算量: 追加で伸ばす長さを m として O(m)
    void ensure(int n) {
        if ((int)fact.size() > n) return;
        assert(n < MOD);
        int old = (int)fact.size() - 1;
        fact.resize(n + 1);
        ifact.resize(n + 1);
        inv.resize(n + 1);
        for (int i = max(1, old + 1); i <= n; ++i) {
            if (i == 1) inv[i] = 1;
            else inv[i] = mint(MOD - (long long)(MOD / i) * inv[MOD % i].v % MOD);
            fact[i] = fact[i - 1] * i;
            ifact[i] = ifact[i - 1] * inv[i];
        }
    }

    // 二項係数 C(n, k) = n! / (k! (n-k)!)
    // 条件を満たさない場合は 0 を返す。
    // 前提: n < MOD
    // 計算量: 償却 O(1)（ensure が必要な場合を除く）
    mint C(long long n, long long k) {
        if (n < 0 || k < 0 || k > n) return 0;
        assert(n < MOD);
        ensure((int)n);
        return fact[n] * ifact[k] * ifact[n - k];
    }

    // 順列数 P(n, k) = n! / (n-k)!
    // 条件を満たさない場合は 0 を返す。
    // 前提: n < MOD
    // 計算量: 償却 O(1)（ensure が必要な場合を除く）
    mint P(long long n, long long k) {
        if (n < 0 || k < 0 || k > n) return 0;
        assert(n < MOD);
        ensure((int)n);
        return fact[n] * ifact[n - k];
    }

    // 重複組合せ H(n, k) = C(n+k-1, k)
    // n 種類のものから重複を許して k 個選ぶ個数。
    // 例: 仕切り法で現れる。
    // 計算量: 償却 O(1)（ensure が必要な場合を除く）
    mint H(long long n, long long k) {
        if (n == 0 && k == 0) return 1;
        if (n <= 0 || k < 0) return 0;
        return C(n + k - 1, k);
    }

    // Catalan 数 Catalan(n) = C(2n, n) / (n+1)
    // 正しい括弧列、二分木、三角形分割などの個数に対応する。
    // 前提: 2n < MOD
    // 計算量: 償却 O(1)（ensure が必要な場合を除く）
    mint catalan(long long n) {
        if (n < 0) return 0;
        assert(2 * n < MOD);
        ensure((int)(2 * n));
        return C(2 * n, n) / mint(n + 1);
    }

    // compositions(n, k)
    // n を k 個の正整数の和で表す順序付き分割の個数。
    // いわゆる ordered partition / composition。
    // 値は C(n-1, k-1)。
    // 計算量: 償却 O(1)（ensure が必要な場合を除く）
    mint compositions(long long n, long long k) {
        if (n < 0 || k < 0) return 0;
        if (n == 0 && k == 0) return 1;
        if (n == 0 || k == 0) return 0;
        return C(n - 1, k - 1);
    }
};

// ============================================================================
// Stirling numbers of the second kind
// ============================================================================

// stirling_second_table(n, k)
// 第2種 Stirling 数 S(i, j) の表を 0<=i<=n, 0<=j<=k について作る。
// S(n, k) は n 個の区別できる要素を k 個の非空グループに分ける個数。
// 漸化式: S(i, j) = S(i-1, j-1) + j * S(i-1, j)
// 計算量: O(nk)
template <int MOD>
vector<vector<ModInt<MOD>>> stirling_second_table(int n, int k) {
    using mint = ModInt<MOD>;
    vector<vector<mint>> dp(n + 1, vector<mint>(k + 1, 0));
    dp[0][0] = 1;
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= min(i, k); ++j) {
            dp[i][j] = dp[i - 1][j - 1] + dp[i - 1][j] * j;
        }
    }
    return dp;
}

// stirling_second(n, k)
// 第2種 Stirling 数 S(n, k) を返す。
// 計算量: O(nk)
template <int MOD>
ModInt<MOD> stirling_second(int n, int k) {
    if (n < 0 || k < 0 || k > n) return 0;
    return stirling_second_table<MOD>(n, k)[n][k];
}

// ============================================================================
// Stirling numbers of the first kind (unsigned)
// ============================================================================

// stirling_first_unsigned_table(n, k)
// 符号なし第1種 Stirling 数 c(i, j) の表を作る。
// c(n, k) は n 要素の順列を k 個の巡回置換に分解する個数。
// 漸化式: c(i, j) = c(i-1, j-1) + (i-1) * c(i-1, j)
// 計算量: O(nk)
template <int MOD>
vector<vector<ModInt<MOD>>> stirling_first_unsigned_table(int n, int k) {
    using mint = ModInt<MOD>;
    vector<vector<mint>> dp(n + 1, vector<mint>(k + 1, 0));
    dp[0][0] = 1;
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= min(i, k); ++j) {
            dp[i][j] = dp[i - 1][j - 1] + dp[i - 1][j] * (i - 1);
        }
    }
    return dp;
}

// stirling_first_unsigned(n, k)
// 符号なし第1種 Stirling 数 c(n, k) を返す。
// 計算量: O(nk)
template <int MOD>
ModInt<MOD> stirling_first_unsigned(int n, int k) {
    if (n < 0 || k < 0 || k > n) return 0;
    return stirling_first_unsigned_table<MOD>(n, k)[n][k];
}

// ============================================================================
// Bell number
// ============================================================================

// bell_number(n)
// Bell 数 B(n) を返す。
// n 個の区別できる要素の集合分割総数。
// B(n) = sum_{k=0}^{n} S(n, k)
// 計算量: O(n^2)
template <int MOD>
ModInt<MOD> bell_number(int n) {
    if (n < 0) return 0;
    auto s = stirling_second_table<MOD>(n, n);
    ModInt<MOD> ans = 0;
    for (int k = 0; k <= n; ++k) ans += s[n][k];
    return ans;
}

// ============================================================================
// Derangement
// ============================================================================

// derangement(n)
// 撹乱順列 !n を返す。
// n 個の要素の順列で、どの要素も元の位置に来ないものの個数。
// 漸化式: !0=1, !1=0, !n=(n-1)(!(n-1)+!(n-2))
// 計算量: O(n)
template <int MOD>
ModInt<MOD> derangement(int n) {
    using mint = ModInt<MOD>;
    if (n < 0) return 0;
    if (n == 0) return 1;
    if (n == 1) return 0;
    mint a = 1, b = 0;
    for (int i = 2; i <= n; ++i) {
        mint c = mint(i - 1) * (a + b);
        a = b;
        b = c;
    }
    return b;
}

// ============================================================================
// Integer partition
// ============================================================================

// partition_exact_k_table(n, k)
// dp[sum][parts] = sum をちょうど parts 個の正整数の和で表す方法数
// （順序は無視）を表す表を作る。
// 漸化式: p(sum, parts) = p(sum-1, parts-1) + p(sum-parts, parts)
// 前者は 1 を 1 個使う場合、後者は各部分から 1 を引ける場合に対応する。
// 計算量: O(nk)
template <int MOD>
vector<vector<ModInt<MOD>>> partition_exact_k_table(int n, int k) {
    using mint = ModInt<MOD>;
    vector<vector<mint>> dp(n + 1, vector<mint>(k + 1, 0));
    dp[0][0] = 1;
    for (int sum = 1; sum <= n; ++sum) {
        for (int parts = 1; parts <= min(sum, k); ++parts) {
            dp[sum][parts] = dp[sum - 1][parts - 1];
            if (sum - parts >= 0) dp[sum][parts] += dp[sum - parts][parts];
        }
    }
    return dp;
}

// partition_exact_k(n, k)
// n をちょうど k 個の正整数の和で表す方法数を返す。
// 順序は無視する。
// 例: 7 = 3+2+2 と 2+3+2 は同じものとして数える。
// 計算量: O(nk)
template <int MOD>
ModInt<MOD> partition_exact_k(int n, int k) {
    if (n < 0 || k < 0 || k > n) return 0;
    return partition_exact_k_table<MOD>(n, k)[n][k];
}

// partition_number(n)
// n を正整数の和で表す方法数の総数を返す。
// 順序は無視する。いわゆる整数分割数 p(n)。
// 計算量: O(n^2)
template <int MOD>
ModInt<MOD> partition_number(int n) {
    if (n < 0) return 0;
    auto dp = partition_exact_k_table<MOD>(n, n);
    ModInt<MOD> ans = 0;
    for (int k = 0; k <= n; ++k) ans += dp[n][k];
    return ans;
}

// ============================================================================
// Surjection count
// ============================================================================

// surjection_count(n, k, comb)
// n 個の区別できる要素から k 個の区別できる箱への全射の個数を返す。
// すべての箱が少なくとも 1 回使われる写像の個数。
// 値は k! * S(n, k)。
// 計算量: O(nk)
template <int MOD>
ModInt<MOD> surjection_count(int n, int k, Comb<MOD>& comb) {
    if (n < 0 || k < 0 || k > n) return 0;
    comb.ensure(k);
    return comb.fact[k] * stirling_second<MOD>(n, k);
}

} // namespace cp_counting
