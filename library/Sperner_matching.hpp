#include <bits/stdc++.h>
using namespace std;

using Mask = unsigned long long;

/*
    Sperner の定理・対称鎖分解を用いたレベル間最大マッチング構築

    このファイルでできること
    ------------------------
    n 個の要素 {0, 1, ..., n-1} の部分集合全体 B_n を考える。
    そのうち、

        B_{n,k}     : 要素数が k の部分集合全体
        B_{n,k+1}   : 要素数が k+1 の部分集合全体

    の間で、次の条件を満たすペア (S, T) をできるだけ多く作る。

        |S| = k
        |T| = k+1
        S ⊂ T

    さらに、同じ S や同じ T は 2 回以上使わない。
    つまり、B_{n,k} と B_{n,k+1} の間の二部グラフにおける
    最大マッチングを構築する。

    集合の表し方
    ------------
    集合は unsigned long long の bitmask で表す。

        bit i が 1 である  <=>  要素 i を含む

    例:

        mask = 0b0101

    は集合 {0, 2} を表す。

    使い方
    ------
    例: n = 4, k = 1 のとき、
        1 要素集合を、それを含む 2 要素集合に重複なく対応させたい。

        int n = 4;
        int k = 1;
        auto matching = sperner_level_matching(n, k);

        for (auto [S, T] : matching) {
            cout << mask_to_set_string(S, n) << " -> "
                 << mask_to_set_string(T, n) << '\n';
        }

    出力例の一つ:

        {0} -> {0,1}
        {3} -> {0,3}
        {2} -> {2,3}
        {1} -> {1,2}

    出力されるマッチングのサイズ
    ----------------------------
    この関数が返すマッチングのサイズは

        min(C(n,k), C(n,k+1))

    になる。これは最大サイズである。

    計算量
    ------
    symmetric_chain_decomposition(n):
        時間計算量: O(n * 2^n)
        空間計算量: O(2^n)

    sperner_level_matching(n, k):
        時間計算量: O(n * 2^n)
        空間計算量: O(2^n)

    注意
    ----
    すべての部分集合を内部的に扱うため、n が大きいと使えない。
    unsigned long long を使っているので、bitmask としては n <= 64 が上限だが、
    計算量の都合上、実用的には n <= 20 程度を目安にする。
*/

/*
    B_n の対称鎖分解を構築する。

    返り値:
        vector<vector<Mask>> chains

    各 chains[i] は 1 本の鎖を表す。
    たとえば、

        chains[i] = {S0, S1, S2, ...}

    なら、

        S0 ⊂ S1 ⊂ S2 ⊂ ...

    となっている。

    また、各隣接ペアは要素を 1 つだけ追加した関係になっている。
    つまり、|S_{j+1}| = |S_j| + 1 である。

    計算量:
        時間計算量: O(n * 2^n)
        空間計算量: O(2^n)
*/
vector<vector<Mask>> symmetric_chain_decomposition(int n) {
    vector<vector<Mask>> chains;

    // B_0 は空集合だけからなる。
    // したがって、鎖は { empty set } の 1 本だけ。
    chains.push_back({0});

    // 要素 x を 0, 1, ..., n-1 の順に追加しながら、
    // B_x の対称鎖分解から B_{x+1} の対称鎖分解を作る。
    for (int x = 0; x < n; x++) {
        Mask bit = 1ULL << x;
        vector<vector<Mask>> next_chains;

        for (const auto& chain : chains) {
            int m = (int)chain.size();

            // 1 本目の鎖:
            // 元の鎖の最後に「最後の集合へ x を追加した集合」をつなげる。
            //
            // S0 ⊂ S1 ⊂ ... ⊂ S_last
            // から
            // S0 ⊂ S1 ⊂ ... ⊂ S_last ⊂ S_last ∪ {x}
            // を作る。
            vector<Mask> c1 = chain;
            c1.push_back(chain.back() | bit);
            next_chains.push_back(c1);

            // 2 本目の鎖:
            // 元の鎖の最後以外に x を追加したものを並べる。
            //
            // S0 ∪ {x} ⊂ S1 ∪ {x} ⊂ ... ⊂ S_{last-1} ∪ {x}
            //
            // m == 1 の場合は空の鎖になってしまうため作らない。
            if (m >= 2) {
                vector<Mask> c2;
                for (int i = 0; i + 1 < m; i++) {
                    c2.push_back(chain[i] | bit);
                }
                next_chains.push_back(c2);
            }
        }

        chains = move(next_chains);
    }

    return chains;
}

/*
    B_{n,k} と B_{n,k+1} の間の最大マッチングを構築する。

    引数:
        n:
            全体集合 {0, 1, ..., n-1} の要素数。

        k:
            左側に置く集合のサイズ。
            左側は B_{n,k}、右側は B_{n,k+1} になる。

    返り値:
        vector<pair<Mask, Mask>> matching

        各 pair は (S, T) を表す。

            pair.first  = S  : k 要素集合
            pair.second = T  : k+1 要素集合

        必ず S ⊂ T を満たす。
        また、同じ S や同じ T は複数回出現しない。

    構築方法:
        1. B_n の対称鎖分解を作る。
        2. 各鎖の中で、サイズ k の集合の直後にあるサイズ k+1 の集合をペアにする。
        3. 鎖どうしは互いに交わらないため、同じ集合を重複して使うことがない。

    返り値のサイズ:
        min(C(n,k), C(n,k+1))

    計算量:
        時間計算量: O(n * 2^n)
        空間計算量: O(2^n)

    注意:
        k < 0 または k >= n の場合、B_{n,k+1} が意味を持たないので空配列を返す。
*/
vector<pair<Mask, Mask>> sperner_level_matching(int n, int k) {
    vector<pair<Mask, Mask>> matching;

    if (k < 0 || k >= n) return matching;
    if (n < 0 || n > 63) return matching;

    auto chains = symmetric_chain_decomposition(n);

    for (const auto& chain : chains) {
        for (int i = 0; i + 1 < (int)chain.size(); i++) {
            Mask small = chain[i];
            Mask large = chain[i + 1];

            // 鎖の中ではサイズが 1 ずつ増える。
            // そのため、small が k 要素集合なら、large は k+1 要素集合である。
            if (__builtin_popcountll(small) == k) {
                matching.push_back({small, large});
            }
        }
    }

    return matching;
}

/*
    bitmask を {0,2,5} のような文字列に変換する補助関数。

    デバッグや出力確認用。

    計算量:
        時間計算量: O(n)
        空間計算量: O(n)
*/
string mask_to_set_string(Mask mask, int n) {
    string s = "{";
    bool first = true;

    for (int i = 0; i < n; i++) {
        if ((mask >> i) & 1ULL) {
            if (!first) s += ",";
            first = false;
            s += to_string(i);
        }
    }

    s += "}";
    return s;
}

/*
    マッチングが正しいか確認する補助関数。

    確認する条件:
        1. 各ペア (S, T) について |S| = k
        2. 各ペア (S, T) について |T| = k+1
        3. 各ペア (S, T) について S ⊂ T
        4. 同じ S が 2 回以上使われていない
        5. 同じ T が 2 回以上使われていない

    返り値:
        true  : 正しい
        false : どこかに誤りがある

    計算量:
        時間計算量: O(M)
            M は matching.size()
        空間計算量: O(M)
*/
bool verify_level_matching(int n, int k, const vector<pair<Mask, Mask>>& matching) {
    (void)n;

    unordered_set<Mask> used_left;
    unordered_set<Mask> used_right;

    for (auto [S, T] : matching) {
        if (__builtin_popcountll(S) != k) return false;
        if (__builtin_popcountll(T) != k + 1) return false;

        // S ⊂ T であることを確認する。
        // bitmask では、S の 1-bit がすべて T に含まれていれば (S & T) == S になる。
        if ((S & T) != S) return false;

        // 今回は |T| = |S| + 1 なので、(S & T) == S なら S != T も自動的に成り立つ。

        if (used_left.count(S)) return false;
        if (used_right.count(T)) return false;

        used_left.insert(S);
        used_right.insert(T);
    }

    return true;
}

/*
    使用例
    ------
    このファイルを単体で実行して動作確認したい場合は、
    コンパイル時に SPERNER_MATCHING_EXAMPLE を定義する。

    例:

        g++ -std=c++17 -O2 -DSPERNER_MATCHING_EXAMPLE sperner_level_matching_jp.cpp
        ./a.out

    ライブラリとして使う場合は、SPERNER_MATCHING_EXAMPLE を定義しなければよい。
    その場合、この main 関数はコンパイルされない。
*/
#ifdef SPERNER_MATCHING_EXAMPLE
int main() {
    int n = 4;
    int k = 1;

    auto matching = sperner_level_matching(n, k);

    cout << "n = " << n << ", k = " << k << '\n';
    cout << "matching size = " << matching.size() << '\n';
    cout << "valid = " << (verify_level_matching(n, k, matching) ? "true" : "false") << '\n';
    cout << '\n';

    for (auto [S, T] : matching) {
        cout << mask_to_set_string(S, n)
             << " -> "
             << mask_to_set_string(T, n)
             << '\n';
    }

    return 0;
}
#endif
