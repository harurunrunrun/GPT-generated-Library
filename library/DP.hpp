#pragma once

/*
================================================================================
  dp_library.hpp
  Header-only C++17 Dynamic Programming library.

  Goal:
    - Collect well-known DP patterns and DP optimizations in one place.
    - Let you write the recurrence declaratively and plug it into a fitting solver.
    - Cover classical DP, tree/bit/digit/interval DP, and major optimizations.

  Important note:
    - No library can fully "auto-solve" arbitrary DP problems from a statement alone.
    - What this library does is abstract the common DP shapes so that once you identify
      the recurrence family, the implementation becomes almost mechanical.

  Included:
    1. Generic memoized DP
    2. LIS / LCS / Edit Distance / DAG DP
    3. Knapsack / subset sum / bounded knapsack
    4. Bit DP (TSP / SOS DP)
    5. Digit DP
    6. Interval DP (generic cubic)
    7. Rerooting Tree DP
    8. Divide-and-Conquer optimization
    9. Knuth optimization
   10. Monge / Totally Monotone optimization via SMAWK
   11. Partition DP builder with automatic backend selection by hint
   12. Li Chao Tree (classic DP optimization helper)
   13. Alien DP / WQS binary search

  Standard:
    - C++17

  Basic usage style:
    - Use the ready-made functions for classic DP.
    - Use the generic frameworks when the recurrence is problem-specific.

================================================================================
*/

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <deque>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace dplib {

// -----------------------------------------------------------------------------
// Small helpers
// -----------------------------------------------------------------------------

template <class T>
constexpr bool chmin(T& a, const T& b) {
    if (b < a) {
        a = b;
        return true;
    }
    return false;
}

template <class T>
constexpr bool chmax(T& a, const T& b) {
    if (a < b) {
        a = b;
        return true;
    }
    return false;
}

template <class T>
constexpr T inf_v() {
    return std::numeric_limits<T>::max() / 4;
}

// -----------------------------------------------------------------------------
// 1. Generic Memoized DP
// -----------------------------------------------------------------------------
/*
Usage:
  - Use this when the state space is sparse or irregular.
  - You provide a recursive lambda that receives (self, state).

Example:
  struct State {
      int i, j;
      bool operator==(const State& other) const { return i == other.i && j == other.j; }
  };
  struct Hash {
      size_t operator()(const State& s) const {
          return (static_cast<size_t>(s.i) << 32) ^ static_cast<size_t>(s.j);
      }
  };

  dplib::MemoDP<State, long long, Hash> solver(
      [&](auto& self, const State& st) -> long long {
          if (st.i == 0) return st.j;
          long long ans = 0;
          ans += self.solve({st.i - 1, st.j});
          if (st.j > 0) ans += self.solve({st.i, st.j - 1});
          return ans;
      }
  );
  auto ans = solver.solve({10, 20});
*/

template <class State, class Value, class Hash = std::hash<State>>
class MemoDP {
public:
    using Recurrence = std::function<Value(MemoDP&, const State&)>;

    explicit MemoDP(Recurrence rec) : rec_(std::move(rec)) {}

    Value solve(const State& s) {
        auto it = memo_.find(s);
        if (it != memo_.end()) return it->second;
        Value val = rec_(*this, s);
        memo_.emplace(s, val);
        return val;
    }

    bool contains(const State& s) const {
        return memo_.find(s) != memo_.end();
    }

    const std::unordered_map<State, Value, Hash>& memo() const {
        return memo_;
    }

    void clear() {
        memo_.clear();
    }

private:
    std::unordered_map<State, Value, Hash> memo_;
    Recurrence rec_;
};

// -----------------------------------------------------------------------------
// 2. Classical DP utilities
// -----------------------------------------------------------------------------

/*
LIS (Longest Increasing Subsequence)
-----------------------------------
Usage:
  int len = dplib::lis_length(a);            // strict
  int len2 = dplib::lis_length(a, false);    // non-decreasing

  auto each = dplib::lis_length_ending_at(a);
  each[i] = LIS length ending exactly at i.
*/

template <class T>
int lis_length(const std::vector<T>& a, bool strict = true) {
    std::vector<T> dp;
    for (const T& x : a) {
        auto it = strict ? std::lower_bound(dp.begin(), dp.end(), x)
                         : std::upper_bound(dp.begin(), dp.end(), x);
        if (it == dp.end()) dp.push_back(x);
        else *it = x;
    }
    return static_cast<int>(dp.size());
}

template <class T>
std::vector<int> lis_length_ending_at(const std::vector<T>& a, bool strict = true) {
    std::vector<T> dp;
    std::vector<int> res(a.size());
    for (int i = 0; i < static_cast<int>(a.size()); ++i) {
        auto it = strict ? std::lower_bound(dp.begin(), dp.end(), a[i])
                         : std::upper_bound(dp.begin(), dp.end(), a[i]);
        int pos = static_cast<int>(it - dp.begin());
        if (it == dp.end()) dp.push_back(a[i]);
        else *it = a[i];
        res[i] = pos + 1;
    }
    return res;
}

/*
LCS (Longest Common Subsequence)
--------------------------------
Usage:
  int len = dplib::lcs_length(s, t);

Time:
  O(|s| * |t|)
*/
inline int lcs_length(const std::string& s, const std::string& t) {
    int n = static_cast<int>(s.size());
    int m = static_cast<int>(t.size());
    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            chmax(dp[i + 1][j + 1], dp[i][j + 1]);
            chmax(dp[i + 1][j + 1], dp[i + 1][j]);
            chmax(dp[i + 1][j + 1], dp[i][j] + (s[i] == t[j]));
        }
    }
    return dp[n][m];
}

/*
Edit Distance (Levenshtein distance)
------------------------------------
Usage:
  int dist = dplib::edit_distance(s, t);

Allowed operations:
  - insert
  - delete
  - replace

Time:
  O(|s| * |t|)
*/
inline int edit_distance(const std::string& s, const std::string& t) {
    int n = static_cast<int>(s.size());
    int m = static_cast<int>(t.size());
    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));
    for (int i = 0; i <= n; ++i) dp[i][0] = i;
    for (int j = 0; j <= m; ++j) dp[0][j] = j;
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            dp[i][j] = std::min({
                dp[i - 1][j] + 1,
                dp[i][j - 1] + 1,
                dp[i - 1][j - 1] + (s[i - 1] != t[j - 1])
            });
        }
    }
    return dp[n][m];
}

/*
DAG DP / Longest Path in DAG
----------------------------
Usage:
  std::vector<std::vector<std::pair<int,long long>>> g(n);
  // g[u].push_back({v, weight});
  auto [dist, order] = dplib::dag_longest_path<long long>(n, g, 0);

Meaning:
  - dist[v] = longest distance from source to v.
  - unreachable vertices get -INF.
  - order is one valid topological order.

Requirement:
  - graph must be a DAG.
*/

template <class T>
std::pair<std::vector<T>, std::vector<int>>
dag_longest_path(int n,
                 const std::vector<std::vector<std::pair<int, T>>>& g,
                 int source,
                 T neg_inf = -inf_v<T>()) {
    std::vector<int> indeg(n, 0);
    for (int u = 0; u < n; ++u) {
        for (auto [v, w] : g[u]) {
            (void)w;
            ++indeg[v];
        }
    }
    std::queue<int> q;
    for (int i = 0; i < n; ++i) if (indeg[i] == 0) q.push(i);
    std::vector<int> order;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        order.push_back(u);
        for (auto [v, w] : g[u]) {
            (void)w;
            if (--indeg[v] == 0) q.push(v);
        }
    }
    assert(static_cast<int>(order.size()) == n && "Graph is not a DAG.");

    std::vector<T> dist(n, neg_inf);
    dist[source] = 0;
    for (int u : order) {
        if (dist[u] == neg_inf) continue;
        for (auto [v, w] : g[u]) {
            chmax(dist[v], dist[u] + w);
        }
    }
    return {dist, order};
}

// -----------------------------------------------------------------------------
// 3. Knapsack / subset sum
// -----------------------------------------------------------------------------

/*
0/1 Knapsack (maximize value under weight limit)
------------------------------------------------
Usage:
  std::vector<std::pair<int,long long>> items = {{weight, value}, ...};
  auto dp = dplib::knapsack01_by_weight<long long>(W, items);
  long long ans = *std::max_element(dp.begin(), dp.end());

Meaning:
  - dp[w] = maximum value with total weight exactly w (or -INF if impossible).

Common pattern:
  - If you want "at most W", take max over dp[0..W].
*/

template <class T>
std::vector<T> knapsack01_by_weight(int W,
                                    const std::vector<std::pair<int, T>>& items,
                                    T neg_inf = -inf_v<T>()) {
    std::vector<T> dp(W + 1, neg_inf);
    dp[0] = 0;
    for (auto [wt, val] : items) {
        for (int w = W; w >= wt; --w) {
            if (dp[w - wt] == neg_inf) continue;
            chmax(dp[w], dp[w - wt] + val);
        }
    }
    return dp;
}

/*
Unbounded Knapsack
------------------
Usage:
  auto dp = dplib::knapsack_unbounded_by_weight<long long>(W, items);
*/

template <class T>
std::vector<T> knapsack_unbounded_by_weight(int W,
                                            const std::vector<std::pair<int, T>>& items,
                                            T neg_inf = -inf_v<T>()) {
    std::vector<T> dp(W + 1, neg_inf);
    dp[0] = 0;
    for (auto [wt, val] : items) {
        for (int w = wt; w <= W; ++w) {
            if (dp[w - wt] == neg_inf) continue;
            chmax(dp[w], dp[w - wt] + val);
        }
    }
    return dp;
}

/*
Bounded Knapsack
----------------
Usage:
  Each item is (weight, value, count).
  auto dp = dplib::knapsack_bounded_by_weight<long long>(W, items);

Implementation:
  - Binary splitting.
  - Good general-purpose implementation.
*/

template <class T>
std::vector<T> knapsack_bounded_by_weight(
    int W,
    const std::vector<std::tuple<int, T, int>>& items,
    T neg_inf = -inf_v<T>()) {
    std::vector<std::pair<int, T>> expanded;
    for (auto [wt, val, cnt] : items) {
        int k = 1;
        while (cnt > 0) {
            int take = std::min(k, cnt);
            expanded.push_back({wt * take, val * static_cast<T>(take)});
            cnt -= take;
            k <<= 1;
        }
    }
    return knapsack01_by_weight<T>(W, expanded, neg_inf);
}

/*
Subset Sum (reachable sums)
---------------------------
Usage:
  auto ok = dplib::subset_sum_reachable(a, S);
  if (ok[x]) { ... }

Meaning:
  - ok[s] = whether sum s is reachable.

Time:
  O(nS)
*/
inline std::vector<char> subset_sum_reachable(const std::vector<int>& a, int S) {
    std::vector<char> dp(S + 1, 0);
    dp[0] = 1;
    for (int x : a) {
        for (int s = S; s >= x; --s) {
            dp[s] = dp[s] || dp[s - x];
        }
    }
    return dp;
}

/*
Subset Sum count modulo MOD
---------------------------
Usage:
  auto ways = dplib::subset_sum_count_mod(a, S, MOD);
  ways[s] = number of ways to make sum s modulo MOD.
*/
inline std::vector<long long> subset_sum_count_mod(const std::vector<int>& a, int S, long long MOD) {
    std::vector<long long> dp(S + 1, 0);
    dp[0] = 1;
    for (int x : a) {
        for (int s = S; s >= x; --s) {
            dp[s] += dp[s - x];
            if (dp[s] >= MOD) dp[s] -= MOD;
        }
    }
    return dp;
}

// -----------------------------------------------------------------------------
// 4. Bit DP
// -----------------------------------------------------------------------------

/*
Traveling Salesman Problem (minimum Hamiltonian cycle)
------------------------------------------------------
Usage:
  std::vector<std::vector<long long>> dist(n, std::vector<long long>(n));
  long long ans = dplib::tsp_min_cycle<long long>(dist, 0);

Meaning:
  - Start at `start`, visit every vertex exactly once, return to `start`.

Time:
  O(n^2 2^n)

Practical range:
  - Roughly n <= 20
*/

template <class T>
T tsp_min_cycle(const std::vector<std::vector<T>>& dist,
                int start = 0,
                T INF = inf_v<T>()) {
    int n = static_cast<int>(dist.size());
    int FULL = 1 << n;
    std::vector<std::vector<T>> dp(FULL, std::vector<T>(n, INF));
    dp[1 << start][start] = 0;
    for (int mask = 0; mask < FULL; ++mask) {
        for (int u = 0; u < n; ++u) {
            if (dp[mask][u] == INF) continue;
            for (int v = 0; v < n; ++v) {
                if (mask & (1 << v)) continue;
                chmin(dp[mask | (1 << v)][v], dp[mask][u] + dist[u][v]);
            }
        }
    }
    T ans = INF;
    int all = FULL - 1;
    for (int u = 0; u < n; ++u) {
        if (dp[all][u] == INF) continue;
        chmin(ans, dp[all][u] + dist[u][start]);
    }
    return ans;
}

/*
SOS DP (subset zeta transform)
------------------------------
Usage:
  std::vector<long long> f(1 << n);
  dplib::subset_zeta_transform(f);

Meaning after call:
  - f[mask] = sum of original f[sub] over all sub ⊆ mask.

This is useful for:
  - subset DP acceleration
  - counting over all submasks
*/

template <class T>
void subset_zeta_transform(std::vector<T>& f) {
    int N = static_cast<int>(f.size());
    int n = 0;
    while ((1 << n) < N) ++n;
    assert((1 << n) == N && "Size must be a power of two.");
    for (int bit = 0; bit < n; ++bit) {
        for (int mask = 0; mask < N; ++mask) {
            if (mask & (1 << bit)) f[mask] += f[mask ^ (1 << bit)];
        }
    }
}

/*
SOS DP (subset Möbius transform)
--------------------------------
Usage:
  This is the inverse of subset_zeta_transform.
*/

template <class T>
void subset_mobius_transform(std::vector<T>& f) {
    int N = static_cast<int>(f.size());
    int n = 0;
    while ((1 << n) < N) ++n;
    assert((1 << n) == N && "Size must be a power of two.");
    for (int bit = 0; bit < n; ++bit) {
        for (int mask = 0; mask < N; ++mask) {
            if (mask & (1 << bit)) f[mask] -= f[mask ^ (1 << bit)];
        }
    }
}

// -----------------------------------------------------------------------------
// 5. Digit DP
// -----------------------------------------------------------------------------
/*
Generic Digit DP framework
--------------------------
This solves DP over decimal digits for numbers in [0, upper].

You provide:
  - State type `State`
  - Value type `Value`
  - init_state()
  - transition(state, digit, started_before, started_after) -> optional<State>
  - terminal_value(state, started) -> Value
  - merge(acc, child_value): combine answers from all digit choices
  - zero: additive identity for Value

Typical use cases:
  - count numbers with some digit property
  - sum numbers satisfying constraints
  - DP on digit automata

Important interpretation of `started`:
  - `started == false` means we have only seen leading zeros so far.
  - When digit != 0, `started_after` becomes true.
  - If the entire number is 0, then at the terminal state `started == false`.
  - So you can decide whether to count 0 by handling terminal_value appropriately.

Example: count numbers with no digit 4 in [0, N]
------------------------------------------------
  struct S {};
  dplib::DigitDP<S, long long> dp(
      []() { return S{}; },
      [](const S&, int digit, bool started_before, bool started_after) -> std::optional<S> {
          (void)started_before;
          (void)started_after;
          if (digit == 4) return std::nullopt;
          return S{};
      },
      [](const S&, bool started) -> long long {
          // count all numbers including 0
          return 1;
      },
      [](long long& acc, const long long& x) {
          acc += x;
      },
      0LL
  );
  long long ans = dp.solve_upto(123456789LL);
*/

template <class State, class Value, class Hash = std::hash<State>>
class DigitDP {
public:
    using InitState = std::function<State()>;
    using Transition = std::function<std::optional<State>(const State&, int, bool, bool)>;
    using Terminal = std::function<Value(const State&, bool)>;
    using Merge = std::function<void(Value&, const Value&)>;

    DigitDP(InitState init_state,
            Transition transition,
            Terminal terminal,
            Merge merge,
            Value zero)
        : init_state_(std::move(init_state)),
          transition_(std::move(transition)),
          terminal_(std::move(terminal)),
          merge_(std::move(merge)),
          zero_(std::move(zero)) {}

    Value solve_upto(long long upper) {
        if (upper < 0) return zero_;
        digits_ = std::to_string(upper);
        memo_.clear();
        return dfs(0, false, false, init_state_());
    }

private:
    struct Key {
        int pos;
        bool started;
        State st;
        bool operator==(const Key& other) const {
            return pos == other.pos && started == other.started && st == other.st;
        }
    };

    struct KeyHash {
        size_t operator()(const Key& k) const {
            size_t h1 = std::hash<int>{}(k.pos);
            size_t h2 = std::hash<bool>{}(k.started);
            size_t h3 = Hash{}(k.st);
            return h1 * 1315423911u ^ (h2 << 1) ^ (h3 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
        }
    };

    Value dfs(int pos, bool tight, bool started, const State& st) {
        if (pos == static_cast<int>(digits_.size())) {
            return terminal_(st, started);
        }
        if (!tight) {
            Key key{pos, started, st};
            auto it = memo_.find(key);
            if (it != memo_.end()) return it->second;
            Value ans = zero_;
            int lim = 9;
            for (int d = 0; d <= lim; ++d) {
                bool started2 = started || (d != 0);
                auto nxt = transition_(st, d, started, started2);
                if (!nxt.has_value()) continue;
                Value child = dfs(pos + 1, false, started2, *nxt);
                merge_(ans, child);
            }
            memo_.emplace(key, ans);
            return ans;
        } else {
            Value ans = zero_;
            int lim = digits_[pos] - '0';
            for (int d = 0; d <= lim; ++d) {
                bool tight2 = tight && (d == lim);
                bool started2 = started || (d != 0);
                auto nxt = transition_(st, d, started, started2);
                if (!nxt.has_value()) continue;
                Value child = dfs(pos + 1, tight2, started2, *nxt);
                merge_(ans, child);
            }
            return ans;
        }
    }

    std::string digits_;
    std::unordered_map<Key, Value, KeyHash> memo_;
    InitState init_state_;
    Transition transition_;
    Terminal terminal_;
    Merge merge_;
    Value zero_;
};

// -----------------------------------------------------------------------------
// 6. Generic Interval DP
// -----------------------------------------------------------------------------
/*
Generic interval DP of the form:
  dp[l][r] = best over split k of trans(l, k, r, dp[l][k], dp[k][r])

Half-open interval convention:
  - interval is [l, r)
  - valid range: 0 <= l <= r <= n

Usage:
  auto dp = dplib::interval_dp<long long>(
      n,
      INF,
      [](int l, int r) -> long long {
          if (r - l <= 1) return 0;
          return INF;
      },
      [&](int l, int k, int r, long long left, long long right) -> long long {
          return left + right + cost(l, r);
      },
      true // minimize
  );
  answer = dp[0][n];

This covers classic interval DP such as:
  - matrix chain multiplication
  - optimal merge pattern on intervals
  - removing intervals
  - palindrome-like interval transitions

When Knuth optimization applies, use the dedicated solver below instead.
*/

template <class T, class BaseFn, class TransFn>
std::vector<std::vector<T>> interval_dp(int n,
                                        T bad,
                                        BaseFn base_fn,
                                        TransFn trans_fn,
                                        bool minimize = true) {
    std::vector<std::vector<T>> dp(n + 1, std::vector<T>(n + 1, bad));
    for (int l = 0; l <= n; ++l) {
        dp[l][l] = base_fn(l, l);
    }
    for (int len = 1; len <= n; ++len) {
        for (int l = 0; l + len <= n; ++l) {
            int r = l + len;
            dp[l][r] = base_fn(l, r);
            for (int k = l + 1; k < r; ++k) {
                T cand = trans_fn(l, k, r, dp[l][k], dp[k][r]);
                if (minimize) chmin(dp[l][r], cand);
                else chmax(dp[l][r], cand);
            }
        }
    }
    return dp;
}

// -----------------------------------------------------------------------------
// 7. Rerooting Tree DP
// -----------------------------------------------------------------------------
/*
Generic rerooting DP
--------------------
This is the standard framework for "tree DP for all roots".

You provide:
  - DP type `DP`
  - EdgeData type `E`
  - merge(a, b): combine contributions from children
  - identity(): neutral element for merge
  - add_edge(dp_from_child, edge): converts a child's DP value when passing through an edge
  - add_root(merged_children, vertex): finalizes value at the current root vertex

Usage example skeleton:
  using DP = long long;
  struct E { int w; };

  dplib::Rerooting<DP, E> rr(
      n,
      [](const DP& a, const DP& b) { return std::max(a, b); },
      []() { return 0LL; },
      [](const DP& from_child, const E& e) { return from_child + e.w; },
      [](const DP& merged, int v) { return merged; }
  );
  rr.add_edge(u, v, {w});
  auto ans = rr.solve();

Meaning:
  - ans[v] = answer when the tree is rooted at v.
*/

template <class DP, class EdgeData>
class Rerooting {
public:
    using Merge = std::function<DP(const DP&, const DP&)>;
    using Identity = std::function<DP()>;
    using AddEdge = std::function<DP(const DP&, const EdgeData&)>;
    using AddRoot = std::function<DP(const DP&, int)>;

    Rerooting(int n, Merge merge, Identity identity, AddEdge add_edge, AddRoot add_root)
        : n_(n),
          g_(n),
          merge_(std::move(merge)),
          identity_(std::move(identity)),
          add_edge_(std::move(add_edge)),
          add_root_(std::move(add_root)),
          sub_(n),
          ans_(n) {}

    void add_edge(int u, int v, const EdgeData& data) {
        int iu = static_cast<int>(g_[u].size());
        int iv = static_cast<int>(g_[v].size());
        g_[u].push_back({v, iv, data});
        g_[v].push_back({u, iu, data});
    }

    std::vector<DP> solve(int root = 0) {
        dfs_sub(root, -1);
        dfs_all(root, -1, identity_());
        return ans_;
    }

private:
    struct Edge {
        int to;
        int rev;
        EdgeData data;
    };

    int n_;
    std::vector<std::vector<Edge>> g_;
    Merge merge_;
    Identity identity_;
    AddEdge add_edge_;
    AddRoot add_root_;
    std::vector<DP> sub_, ans_;

    void dfs_sub(int v, int p) {
        DP acc = identity_();
        for (const auto& e : g_[v]) {
            if (e.to == p) continue;
            dfs_sub(e.to, v);
            acc = merge_(acc, add_edge_(sub_[e.to], e.data));
        }
        sub_[v] = add_root_(acc, v);
    }

    void dfs_all(int v, int p, const DP& from_parent) {
        int m = static_cast<int>(g_[v].size());
        std::vector<DP> vals(m);
        for (int i = 0; i < m; ++i) {
            const auto& e = g_[v][i];
            if (e.to == p) vals[i] = add_edge_(from_parent, e.data);
            else vals[i] = add_edge_(sub_[e.to], e.data);
        }
        std::vector<DP> pref(m + 1, identity_()), suff(m + 1, identity_());
        for (int i = 0; i < m; ++i) pref[i + 1] = merge_(pref[i], vals[i]);
        for (int i = m - 1; i >= 0; --i) suff[i] = merge_(vals[i], suff[i + 1]);
        ans_[v] = add_root_(pref[m], v);
        for (int i = 0; i < m; ++i) {
            const auto& e = g_[v][i];
            if (e.to == p) continue;
            DP without_child = merge_(pref[i], suff[i + 1]);
            DP next_parent = add_root_(without_child, v);
            dfs_all(e.to, v, next_parent);
        }
    }
};

// -----------------------------------------------------------------------------
// 8. Divide-and-Conquer optimization
// -----------------------------------------------------------------------------
/*
Recurrence shape:
  cur[j] = min_{k in [optL, optR]} prev[k] + cost(k, j)
with monotone argmin:
  opt[j] <= opt[j+1]

This reduces one DP row from O(n^2) to O(n log n) or O(n) levels of recursion.

Usage:
  auto [cur, arg] = dplib::divide_and_conquer_row<long long>(
      prev,
      n,
      INF,
      [&](int k, int j) { return prev[k] + cost(k, j); },
      0,
      n - 1,
      0,
      n - 1
  );

More convenient interface:
  use PartitionDP below when your recurrence is the standard partition form.
*/

template <class T, class Eval>
void divide_and_conquer_row_impl(int l, int r,
                                 int optl, int optr,
                                 std::vector<T>& cur,
                                 std::vector<int>& arg,
                                 T INF,
                                 Eval eval) {
    if (l > r) return;
    int mid = (l + r) >> 1;
    std::pair<T, int> best = {INF, -1};
    for (int k = optl; k <= std::min(mid, optr); ++k) {
        T cand = eval(k, mid);
        if (cand < best.first) best = {cand, k};
    }
    cur[mid] = best.first;
    arg[mid] = best.second;
    divide_and_conquer_row_impl<T>(l, mid - 1, optl, best.second, cur, arg, INF, eval);
    divide_and_conquer_row_impl<T>(mid + 1, r, best.second, optr, cur, arg, INF, eval);
}

template <class T, class Eval>
std::pair<std::vector<T>, std::vector<int>>
divide_and_conquer_row(int n,
                       T INF,
                       Eval eval,
                       int l = 0,
                       int r = -1,
                       int optl = 0,
                       int optr = -1) {
    if (r == -1) r = n - 1;
    if (optr == -1) optr = n - 1;
    std::vector<T> cur(n, INF);
    std::vector<int> arg(n, -1);
    divide_and_conquer_row_impl<T>(l, r, optl, optr, cur, arg, INF, eval);
    return {cur, arg};
}

// -----------------------------------------------------------------------------
// 9. SMAWK / Monge / Totally Monotone optimization
// -----------------------------------------------------------------------------
/*
SMAWK computes row minima of a totally monotone matrix in O(rows + cols) calls
to the comparison oracle, up to constant factors.

You provide:
  - rows: list of row indices
  - cols: list of column indices
  - better(row, c1, c2): returns true iff column c1 is at least as good as c2 on row

Typical use:
  auto arg = dplib::smawk_argmin(
      rows,
      cols,
      [&](int row, int c1, int c2) {
          return value(row, c1) <= value(row, c2);
      }
  );

Precondition:
  - Matrix must be totally monotone.
  - A Monge matrix satisfies this, so Monge DP is a common use case.
*/

template <class Better>
std::map<int, int> smawk_argmin(std::vector<int> rows,
                                std::vector<int> cols,
                                Better better) {
    if (rows.empty()) return {};

    std::function<std::map<int, int>(std::vector<int>, std::vector<int>)> rec =
        [&](std::vector<int> rs, std::vector<int> cs) -> std::map<int, int> {
            if (rs.empty()) return {};

            std::vector<int> reduced;
            for (int c : cs) {
                while (!reduced.empty()) {
                    int sz = static_cast<int>(reduced.size());
                    int row = rs[sz - 1];
                    if (better(row, c, reduced.back())) reduced.pop_back();
                    else break;
                }
                if (static_cast<int>(reduced.size()) < static_cast<int>(rs.size())) {
                    reduced.push_back(c);
                }
            }

            std::vector<int> odd_rows;
            for (int i = 1; i < static_cast<int>(rs.size()); i += 2) odd_rows.push_back(rs[i]);
            auto ans = rec(odd_rows, reduced);

            int left = 0;
            for (int i = 0; i < static_cast<int>(rs.size()); i += 2) {
                int row = rs[i];
                int right = static_cast<int>(reduced.size()) - 1;
                if (i + 1 < static_cast<int>(rs.size())) {
                    right = static_cast<int>(std::find(reduced.begin(), reduced.end(), ans[rs[i + 1]]) - reduced.begin());
                }
                int best_col = reduced[left];
                for (int j = left; j <= right; ++j) {
                    if (better(row, reduced[j], best_col)) best_col = reduced[j];
                }
                ans[row] = best_col;
                left = right;
            }
            return ans;
        };

    return rec(std::move(rows), std::move(cols));
}

/*
Monge row minima wrapper
------------------------
Usage:
  auto [argmin, minval] = dplib::monge_row_minima<long long>(rows, cols, value);

Where:
  value(r, c) returns matrix[r][c].

This is convenient when you want actual values as well.
*/

template <class T, class ValueFn>
std::pair<std::map<int, int>, std::map<int, T>>
monge_row_minima(const std::vector<int>& rows,
                 const std::vector<int>& cols,
                 ValueFn value) {
    auto arg = smawk_argmin(rows, cols, [&](int r, int c1, int c2) {
        return value(r, c1) <= value(r, c2);
    });
    std::map<int, T> val;
    for (int r : rows) val[r] = value(r, arg[r]);
    return {arg, val};
}

// -----------------------------------------------------------------------------
// 10. Partition DP builder with backend hint
// -----------------------------------------------------------------------------
/*
Standard partition DP form:
  dp[t][j] = min_{0 <= k < j} dp[t-1][k] + cost(k, j)

Interpretation:
  - first j items are partitioned into t groups
  - last group is (k, j] or [k, j), depending on your cost definition

This builder supports three backends:
  1. Naive O(stages * n^2)
  2. Divide-and-Conquer optimization
  3. Monge optimization via SMAWK

Choose by hint:
  - OptimizationHint::None
  - OptimizationHint::DivideAndConquer
  - OptimizationHint::Monge

IMPORTANT:
  - The library cannot verify the mathematical preconditions for you.
  - If you choose DivideAndConquer, you are claiming argmin monotonicity.
  - If you choose Monge, you are claiming total monotonicity / Monge property.

Usage:
  dplib::PartitionDP<long long, decltype(cost)> solver(
      stages, n, INF, cost, dplib::OptimizationHint::DivideAndConquer
  );
  solver.set_base0();  // dp[0][0] = 0, dp[0][j>0] = INF
  auto ans = solver.solve();
  // ans[t][j]
*/

enum class OptimizationHint {
    None,
    DivideAndConquer,
    Monge
};

template <class T, class Cost>
class PartitionDP {
public:
    PartitionDP(int stages, int n, T INF, Cost cost, OptimizationHint hint = OptimizationHint::None)
        : stages_(stages), n_(n), INF_(INF), cost_(std::move(cost)), hint_(hint),
          dp_(stages + 1, std::vector<T>(n + 1, INF_)), arg_(stages + 1, std::vector<int>(n + 1, -1)) {}

    void set_base0() {
        dp_[0][0] = 0;
        for (int j = 1; j <= n_; ++j) dp_[0][j] = INF_;
    }

    void set_base_row(const std::vector<T>& base) {
        assert(static_cast<int>(base.size()) == n_ + 1);
        dp_[0] = base;
    }

    const std::vector<std::vector<T>>& solve() {
        for (int t = 1; t <= stages_; ++t) {
            if (hint_ == OptimizationHint::None) {
                solve_naive_row(t);
            } else if (hint_ == OptimizationHint::DivideAndConquer) {
                solve_dc_row(t);
            } else {
                solve_monge_row(t);
            }
        }
        return dp_;
    }

    const std::vector<std::vector<int>>& argmin() const {
        return arg_;
    }

private:
    int stages_, n_;
    T INF_;
    Cost cost_;
    OptimizationHint hint_;
    std::vector<std::vector<T>> dp_;
    std::vector<std::vector<int>> arg_;

    void solve_naive_row(int t) {
        dp_[t][0] = 0;
        arg_[t][0] = 0;
        for (int j = 1; j <= n_; ++j) {
            std::pair<T, int> best = {INF_, -1};
            for (int k = 0; k < j; ++k) {
                if (dp_[t - 1][k] == INF_) continue;
                T cand = dp_[t - 1][k] + cost_(k, j);
                if (cand < best.first) best = {cand, k};
            }
            dp_[t][j] = best.first;
            arg_[t][j] = best.second;
        }
    }

    void solve_dc_row(int t) {
        auto eval = [&](int k, int j) -> T {
            if (k >= j || dp_[t - 1][k] == INF_) return INF_;
            return dp_[t - 1][k] + cost_(k, j);
        };
        auto [cur, arg] = divide_and_conquer_row<T>(n_ + 1, INF_, eval, 1, n_, 0, n_ - 1);
        dp_[t][0] = 0;
        arg_[t][0] = 0;
        for (int j = 1; j <= n_; ++j) {
            dp_[t][j] = cur[j];
            arg_[t][j] = arg[j];
        }
    }

    void solve_monge_row(int t) {
        std::vector<int> rows, cols;
        for (int j = 1; j <= n_; ++j) rows.push_back(j);
        for (int k = 0; k < n_; ++k) cols.push_back(k);
        auto better = [&](int row, int c1, int c2) {
            T v1 = (c1 < row && dp_[t - 1][c1] != INF_) ? dp_[t - 1][c1] + cost_(c1, row) : INF_;
            T v2 = (c2 < row && dp_[t - 1][c2] != INF_) ? dp_[t - 1][c2] + cost_(c2, row) : INF_;
            return v1 <= v2;
        };
        auto best_col = smawk_argmin(rows, cols, better);
        dp_[t][0] = 0;
        arg_[t][0] = 0;
        for (int j = 1; j <= n_; ++j) {
            int k = best_col[j];
            arg_[t][j] = k;
            dp_[t][j] = (k < j && dp_[t - 1][k] != INF_) ? dp_[t - 1][k] + cost_(k, j) : INF_;
        }
    }
};

// -----------------------------------------------------------------------------
// 11. Knuth optimization
// -----------------------------------------------------------------------------
/*
Knuth optimization solves interval DP of the form:
  dp[l][r] = min_{l < k < r} dp[l][k] + dp[k][r] + w(l, r)
with monotone opt:
  opt[l][r-1] <= opt[l][r] <= opt[l+1][r]

Half-open intervals: [l, r)

Usage:
  auto [dp, opt] = dplib::knuth_optimization<long long>(n, cost, INF);
  answer = dp[0][n];

You provide:
  - n
  - w(l, r): interval cost on [l, r)

Base convention in this implementation:
  - dp[l][l] = 0
  - dp[l][l+1] = 0
  - transitions start from length >= 2

Typical problems:
  - optimal binary search tree
  - interval merge cost under suitable quadrangle inequality
*/

template <class T, class Cost>
std::pair<std::vector<std::vector<T>>, std::vector<std::vector<int>>>
knuth_optimization(int n, Cost w, T INF = inf_v<T>()) {
    std::vector<std::vector<T>> dp(n + 1, std::vector<T>(n + 1, INF));
    std::vector<std::vector<int>> opt(n + 1, std::vector<int>(n + 1, -1));
    for (int l = 0; l <= n; ++l) {
        dp[l][l] = 0;
        opt[l][l] = l;
        if (l + 1 <= n) {
            dp[l][l + 1] = 0;
            opt[l][l + 1] = l + 1;
        }
    }
    for (int len = 2; len <= n; ++len) {
        for (int l = 0; l + len <= n; ++l) {
            int r = l + len;
            int left = opt[l][r - 1];
            int right = opt[l + 1][r];
            if (left == -1) left = l + 1;
            if (right == -1) right = r - 1;
            left = std::max(left, l + 1);
            right = std::min(right, r - 1);
            for (int k = left; k <= right; ++k) {
                T cand = dp[l][k] + dp[k][r] + w(l, r);
                if (cand < dp[l][r]) {
                    dp[l][r] = cand;
                    opt[l][r] = k;
                }
            }
        }
    }
    return {dp, opt};
}

// -----------------------------------------------------------------------------
// 12. Li Chao Tree (DP optimization helper)
// -----------------------------------------------------------------------------
/*
Li Chao Tree for minimum of lines
---------------------------------
This is not a DP by itself, but it is one of the most common DP optimization tools.

Supports:
  - add line: y = ax + b
  - query min at x

Usage:
  dplib::LiChaoTree<long long> lichao(xs, INF);
  lichao.add_line(a, b);
  long long best = lichao.query(x);

Important:
  - Coordinate-compressed static x-set version.
  - Pass every x value you will ever query in `xs` at construction time.

Typical DP form:
  dp[i] = min_j (m_j * x_i + b_j)
*/

template <class T>
class LiChaoTree {
public:
    struct Line {
        T a, b;
        Line(T a_ = 0, T b_ = inf_v<T>()) : a(a_), b(b_) {}
        T get(T x) const { return a * x + b; }
    };

    LiChaoTree() = default;

    LiChaoTree(std::vector<T> xs, T INF = inf_v<T>())
        : xs_(std::move(xs)), n_(1), INF_(INF) {
        std::sort(xs_.begin(), xs_.end());
        xs_.erase(std::unique(xs_.begin(), xs_.end()), xs_.end());
        assert(!xs_.empty() && "LiChaoTree requires at least one query x-coordinate.");
        while (n_ < static_cast<int>(xs_.size())) n_ <<= 1;
        xs_.resize(n_, xs_.back());
        seg_.assign(2 * n_, Line(0, INF_));
    }

    void add_line(T a, T b) {
        add_line(Line(a, b), 1, 0, n_ - 1);
    }

    T query(T x) const {
        int idx = static_cast<int>(std::lower_bound(xs_.begin(), xs_.end(), x) - xs_.begin());
        assert(idx < static_cast<int>(xs_.size()) && xs_[idx] == x && "x must belong to initial coordinate set.");
        idx += n_;
        T res = INF_;
        while (idx > 0) {
            chmin(res, seg_[idx].get(x));
            idx >>= 1;
        }
        return res;
    }

private:
    std::vector<T> xs_;
    int n_ = 1;
    T INF_ = inf_v<T>();
    std::vector<Line> seg_;

    void add_line(Line nw, int node, int l, int r) {
        int mid = (l + r) >> 1;
        T xl = xs_[std::min(l, static_cast<int>(xs_.size()) - 1)];
        T xm = xs_[std::min(mid, static_cast<int>(xs_.size()) - 1)];
        T xr = xs_[std::min(r, static_cast<int>(xs_.size()) - 1)];

        Line lo = seg_[node], hi = nw;
        if (lo.get(xl) > hi.get(xl)) std::swap(lo, hi);
        if (lo.get(xr) <= hi.get(xr)) {
            seg_[node] = lo;
            return;
        }
        if (lo.get(xm) <= hi.get(xm)) {
            seg_[node] = lo;
            if (mid + 1 <= r) add_line(hi, node << 1 | 1, mid + 1, r);
        } else {
            seg_[node] = hi;
            if (l <= mid) add_line(lo, node << 1, l, mid);
        }
    }
};

// -----------------------------------------------------------------------------
// 13. Alien DP / WQS Binary Search
// -----------------------------------------------------------------------------
/*
Alien DP / WQS Binary Search
----------------------------
Use this when the DP objective depends on an exact count k of chosen operations,
segments, facilities, etc., and it is easier to solve the penalized version.

Typical setup:
  Original goal:
    maximize score with exactly K selections

  Penalized solver for penalty p:
    best(score - p * count)
  and also returns the count used in the optimum, with tie-breaking chosen so that
  count is monotone in p.

Then binary search p, and recover:
  exact_answer = penalized_value + p * K

Interface below targets the common "maximize with exactly K actions" case.

You provide `solve_penalized(p)` returning:
  std::pair<Value, long long> = {best_value_minus_p_times_count, used_count}

Tie-breaking rule you should implement inside solve_penalized:
  - If two transitions give the same penalized value, prefer the one with larger count.

Why this tie-break?
  - It makes used_count monotone non-increasing as penalty increases,
    which is the standard assumption for the binary search below.

Usage pattern:
  auto res = dplib::alien_dp_max<long long>(
      low_penalty,
      high_penalty,
      K,
      [&](long long p) -> std::pair<long long,long long> {
          // run DP for penalized problem
          // return {best_penalized_value, used_count}
      }
  );
  // res.best_exact = answer for exactly K

Notes:
  - Penalty range should cover the optimum.
  - Usually integer penalties are used.
*/

template <class T>
struct AlienDPResult {
    T best_exact;
    T penalty;
    T best_penalized;
    long long used_count;
};

template <class T, class SolvePenalized>
AlienDPResult<T> alien_dp_max(T penalty_lo,
                              T penalty_hi,
                              long long target_k,
                              SolvePenalized solve_penalized) {
    T lo = penalty_lo, hi = penalty_hi;
    while (lo < hi) {
        T mid = lo + (hi - lo + 1) / 2;
        auto [best, cnt] = solve_penalized(mid);
        (void)best;
        if (cnt >= target_k) lo = mid;
        else hi = mid - 1;
    }
    auto [penalized, cnt] = solve_penalized(lo);
    return AlienDPResult<T>{penalized + lo * static_cast<T>(target_k), lo, penalized, cnt};
}

// -----------------------------------------------------------------------------
// 14. Extra helpers for common DP modeling
// -----------------------------------------------------------------------------
/*
Prefix sums helper
------------------
This is not a DP itself, but many cost(k, j) formulas in partition / interval DP
need O(1) range sums.

Usage:
  dplib::PrefixSum<long long> ps(a);
  long long sum_lr = ps.sum(l, r); // [l, r)
*/

template <class T>
class PrefixSum {
public:
    PrefixSum() = default;
    explicit PrefixSum(const std::vector<T>& a) : pref_(a.size() + 1, 0) {
        for (int i = 0; i < static_cast<int>(a.size()); ++i) pref_[i + 1] = pref_[i] + a[i];
    }
    T sum(int l, int r) const {
        return pref_[r] - pref_[l];
    }
private:
    std::vector<T> pref_;
};

/*
Argmax pair helper for Alien DP or custom DP states
---------------------------------------------------
Stores (value, count) and compares by:
  1. larger value
  2. if tied, larger count

Usage:
  using Node = dplib::MaxWithCount<long long>;
  Node a{value, count};
  Node b{value2, count2};
  Node best = std::max(a, b);
*/

template <class T>
struct MaxWithCount {
    T value;
    long long count;
    bool operator<(const MaxWithCount& other) const {
        if (value != other.value) return value < other.value;
        return count < other.count;
    }
};

} // namespace dplib

