#include <bits/stdc++.h>
using namespace std;

// Exact 0-1 knapsack solver inspired by Pisinger's core-based approach.
// - Sort by value/weight density
// - Find the break item of the greedy fractional solution
// - Branch near the break item first (expanding core order)
// - Use fractional-knapsack upper bounds for pruning
//
// Assumptions:
//   * weights[i] >= 0
//   * values[i]  >= 0
//   * 0-weight positive-value items are handled automatically
//
// Returned "taken" is in the ORIGINAL item order.

template <class V = long long, class W = long long>
struct CoreExactKnapsack {
  struct Result {
    V value = 0;
    W weight = 0;
    vector<int> taken;  // 0/1 in original index order
  };

 private:
  struct Item {
    V v;
    W w;
    int id;
    long double dens;
  };

  int n = 0;
  W C = 0;
  vector<Item> a;             // sorted by density desc
  vector<int> order;          // branching order around break item
  vector<int> state;          // -1 undecided, 0 excluded, 1 included
  vector<int> best_state;     // best complete state in sorted order
  V best_value = 0;
  W best_weight = 0;

  static bool by_density(const Item& x, const Item& y) {
    if (x.dens != y.dens) return x.dens > y.dens;
    if (x.w != y.w) return x.w < y.w;
    if (x.v != y.v) return x.v > y.v;
    return x.id < y.id;
  }

  void update_best(const vector<int>& full_state, W w, V v) {
    if (v > best_value || (v == best_value && w < best_weight)) {
      best_value = v;
      best_weight = w;
      best_state = full_state;
    }
  }

  // Feasible completion heuristic:
  // fill remaining undecided items greedily in density order.
  void greedy_completion(W cur_w, V cur_v) {
    vector<int> cand = state;
    W w = cur_w;
    V v = cur_v;

    for (int i = 0; i < n; ++i) {
      if (cand[i] != -1) continue;
      if (w + a[i].w <= C) {
        cand[i] = 1;
        w += a[i].w;
        v += a[i].v;
      } else {
        cand[i] = 0;
      }
    }
    update_best(cand, w, v);
  }

  // Fractional-knapsack upper bound on undecided items.
  long double upper_bound(W cur_w, V cur_v) const {
    if (cur_w > C) return -1e100L;
    long double ub = static_cast<long double>(cur_v);
    W rem = C - cur_w;

    for (int i = 0; i < n; ++i) {
      if (state[i] != -1) continue;
      if (a[i].w <= rem) {
        rem -= a[i].w;
        ub += static_cast<long double>(a[i].v);
      } else {
        if (a[i].w > 0) {
          ub += static_cast<long double>(a[i].v) *
                (static_cast<long double>(rem) / static_cast<long double>(a[i].w));
        }
        break;
      }
    }
    return ub;
  }

  void dfs(int pos, W cur_w, V cur_v) {
    if (cur_w > C) return;

    // Get a fast incumbent.
    greedy_completion(cur_w, cur_v);

    // Prune by LP relaxation.
    long double ub = upper_bound(cur_w, cur_v);
    if (ub < static_cast<long double>(best_value) + 1e-18L) return;

    if (pos == (int)order.size()) {
      // All items decided.
      vector<int> full = state;
      for (int i = 0; i < n; ++i) {
        if (full[i] == -1) full[i] = 0;
      }
      update_best(full, cur_w, cur_v);
      return;
    }

    int idx = order[pos];
    if (state[idx] != -1) {
      dfs(pos + 1, cur_w, cur_v);
      return;
    }

    // Try include first, then exclude.
    if (cur_w + a[idx].w <= C) {
      state[idx] = 1;
      dfs(pos + 1, cur_w + a[idx].w, cur_v + a[idx].v);
    }

    state[idx] = 0;
    dfs(pos + 1, cur_w, cur_v);

    state[idx] = -1;
  }

 public:
  Result run(const vector<V>& values, const vector<W>& weights, W capacity) {
    const int m = (int)values.size();
    Result ans;
    ans.taken.assign(m, 0);

    if ((int)weights.size() != m) {
      throw invalid_argument("values.size() must equal weights.size()");
    }
    if (capacity < 0) {
      throw invalid_argument("capacity must be nonnegative");
    }

    // Preprocess trivial items.
    V forced_value = 0;
    W forced_weight = 0;  // stays 0 because only zero-weight items are forced
    vector<int> forced_take(m, 0);

    a.clear();
    for (int i = 0; i < m; ++i) {
      if (weights[i] < 0 || values[i] < 0) {
        throw invalid_argument("this implementation assumes nonnegative weights/values");
      }
      if (weights[i] == 0) {
        if (values[i] > 0) {
          forced_take[i] = 1;
          forced_value += values[i];
        }
        continue;
      }
      if (weights[i] > capacity) continue; // impossible to take
      if (values[i] == 0) continue;        // never useful
      a.push_back(Item{
          values[i], weights[i], i,
          static_cast<long double>(values[i]) / static_cast<long double>(weights[i])
      });
    }

    n = (int)a.size();
    C = capacity;
    sort(a.begin(), a.end(), by_density);

    ans.value = forced_value;
    ans.weight = forced_weight;
    for (int i = 0; i < m; ++i) ans.taken[i] = forced_take[i];

    if (n == 0) return ans;

    // If everything fits, take everything.
    W total_w = 0;
    V total_v = 0;
    for (const auto& it : a) {
      total_w += it.w;
      total_v += it.v;
    }
    if (total_w <= C) {
      ans.value += total_v;
      ans.weight += total_w;
      for (const auto& it : a) ans.taken[it.id] = 1;
      return ans;
    }

    // Find break item of greedy-by-density packing.
    W wsum = 0;
    int br = 0;
    while (br < n && wsum + a[br].w <= C) {
      wsum += a[br].w;
      ++br;
    }

    // Expanding-core branching order around the break item.
    order.clear();
    vector<int> used(n, 0);

    auto push_if = [&](int idx) {
      if (0 <= idx && idx < n && !used[idx]) {
        used[idx] = 1;
        order.push_back(idx);
      }
    };

    if (br < n) push_if(br);     // first excluded by greedy
    push_if(br - 1);             // last included by greedy

    for (int d = 1; (br - 1 - d) >= 0 || (br + d) < n; ++d) {
      push_if(br + d);
      push_if(br - 1 - d);
    }

    // In edge cases, ensure all items are present.
    for (int i = 0; i < n; ++i) push_if(i);

    state.assign(n, -1);
    best_state.assign(n, 0);
    best_value = 0;
    best_weight = 0;

    dfs(0, 0, 0);

    // Reconstruct.
    ans.value += best_value;
    ans.weight += best_weight;
    for (int i = 0; i < n; ++i) {
      if (best_state[i] == 1) ans.taken[a[i].id] = 1;
    }
    return ans;
  }
};

// -------------------- example --------------------
// int main() {
//   int N;
//   long long W;
//   cin >> N >> W;
//   vector<long long> v(N), w(N);
//   for (int i = 0; i < N; ++i) cin >> v[i] >> w[i];
//
//   CoreExactKnapsack<long long, long long> solver;
//   auto res = solver.run(v, w, W);
//
//   cout << res.value << '\n';
//   // chosen items in original order:
//   // for (int i = 0; i < N; ++i) if (res.taken[i]) cout << i << ' ';
//   // cout << '\n';
// }
