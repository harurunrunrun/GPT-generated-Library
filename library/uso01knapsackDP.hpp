#include <bits/stdc++.h>
using namespace std;

using ll = long long;
using i128 = __int128_t;

struct BranchAndBoundKnapsack {
  struct Result {
    ll value = 0;
    ll weight = 0;
    vector<int> taken;  // original order
  };

  struct Group {
    ll v, w;
    vector<int> ids;  // original indices of identical items
  };

  int G = 0;
  ll CAP = 0;

  vector<Group> groups;
  vector<ll> suffix_weight;
  vector<i128> suffix_value;

  vector<int> choose_cnt;
  vector<int> best_choose_cnt;

  i128 best_value = 0;
  ll best_weight = 0;

  static bool density_cmp(const Group& a, const Group& b) {
    i128 lhs = (i128)a.v * (i128)b.w;
    i128 rhs = (i128)b.v * (i128)a.w;
    if (lhs != rhs) return lhs > rhs;
    if (a.w != b.w) return a.w < b.w;
    if (a.v != b.v) return a.v > b.v;
    return a.ids.size() > b.ids.size();
  }

  static ll to_ll(i128 x) {
    return (ll)x;
  }

  void update_best(i128 value, ll used_weight, const vector<int>& cand) {
    if (value > best_value || (value == best_value && used_weight < best_weight)) {
      best_value = value;
      best_weight = used_weight;
      best_choose_cnt = cand;
    }
  }

  long double relax_upper_bound(int idx, ll rem) const {
    long double ub = 0.0L;
    for (int i = idx; i < G; ++i) {
      if (rem <= 0) break;
      ll can_take = min<ll>((ll)groups[i].ids.size(), rem / groups[i].w);
      ub += (long double)can_take * (long double)groups[i].v;
      rem -= can_take * groups[i].w;

      if (can_take < (ll)groups[i].ids.size() && rem > 0) {
        ub += (long double)rem * (long double)groups[i].v / (long double)groups[i].w;
        break;
      }
    }
    return ub;
  }

  void greedy_completion(int idx, ll rem, i128 cur_value) {
    vector<int> cand = choose_cnt;
    i128 val = cur_value;
    ll used_weight = CAP - rem;

    for (int i = idx; i < G; ++i) {
      ll can_take = min<ll>((ll)groups[i].ids.size(), rem / groups[i].w);
      cand[i] = (int)can_take;
      rem -= can_take * groups[i].w;
      used_weight += can_take * groups[i].w;
      val += (i128)can_take * (i128)groups[i].v;
    }

    update_best(val, used_weight, cand);
  }

  void dfs(int idx, ll rem, i128 cur_value) {
    if (idx == G) {
      update_best(cur_value, CAP - rem, choose_cnt);
      return;
    }

    if (rem < 0) return;

    // 全残りが入るなら即確定
    if (suffix_weight[idx] <= rem) {
      vector<int> cand = choose_cnt;
      for (int i = idx; i < G; ++i) cand[i] = (int)groups[i].ids.size();
      update_best(cur_value + suffix_value[idx], CAP - rem + suffix_weight[idx], cand);
      return;
    }

    // 分数緩和上界
    long double ub = (long double)cur_value + relax_upper_bound(idx, rem);

    // 目的値は整数なので ub < best + 1 なら改善不能
    if (ub < (long double)best_value + 1.0L - 1e-18L) return;

    // まず貪欲 completion で incumbent を強化
    greedy_completion(idx, rem, cur_value);

    ll max_take = min<ll>((ll)groups[idx].ids.size(), rem / groups[idx].w);

    // 最後のグループだけなら最大個数を取るのが最善
    if (idx == G - 1) {
      choose_cnt[idx] = (int)max_take;
      update_best(cur_value + (i128)max_take * (i128)groups[idx].v,
                  CAP - rem + max_take * groups[idx].w,
                  choose_cnt);
      choose_cnt[idx] = 0;
      return;
    }

    // 高密度グループなので、大きい個数から試す
    // x を減らすほど上界は悪くなるので、途中で break できる
    for (ll x = max_take; x >= 0; --x) {
      ll nrem = rem - x * groups[idx].w;
      i128 nval = cur_value + (i128)x * (i128)groups[idx].v;

      long double child_ub = (long double)nval + relax_upper_bound(idx + 1, nrem);
      if (child_ub < (long double)best_value + 1.0L - 1e-18L) {
        break;
      }

      choose_cnt[idx] = (int)x;
      dfs(idx + 1, nrem, nval);
    }

    choose_cnt[idx] = 0;
  }

  Result run(const vector<ll>& values, const vector<ll>& weights, ll capacity) {
    if (values.size() != weights.size()) {
      throw invalid_argument("values.size() must equal weights.size()");
    }
    if (capacity < 0) {
      throw invalid_argument("capacity must be nonnegative");
    }

    CAP = capacity;
    int N = (int)values.size();

    Result res;
    res.taken.assign(N, 0);

    i128 forced_value = 0;
    ll forced_weight = 0;

    vector<tuple<ll, ll, int>> items;
    items.reserve(N);

    for (int i = 0; i < N; ++i) {
      ll v = values[i];
      ll w = weights[i];

      if (v < 0 || w < 0) {
        throw invalid_argument("this implementation assumes nonnegative values and weights");
      }

      if (w == 0) {
        if (v > 0) {
          res.taken[i] = 1;
          forced_value += v;
        }
        continue;
      }

      if (w > CAP) continue;
      if (v == 0) continue;

      items.emplace_back(v, w, i);
    }

    if (items.empty()) {
      res.value = to_ll(forced_value);
      res.weight = forced_weight;
      return res;
    }

    // 同一 (v, w) を圧縮
    sort(items.begin(), items.end(),
         [](const auto& a, const auto& b) {
           if (get<0>(a) != get<0>(b)) return get<0>(a) < get<0>(b);
           if (get<1>(a) != get<1>(b)) return get<1>(a) < get<1>(b);
           return get<2>(a) < get<2>(b);
         });

    groups.clear();
    for (auto [v, w, id] : items) {
      if (groups.empty() || groups.back().v != v || groups.back().w != w) {
        groups.push_back(Group{v, w, {}});
      }
      groups.back().ids.push_back(id);
    }

    sort(groups.begin(), groups.end(), density_cmp);

    G = (int)groups.size();
    choose_cnt.assign(G, 0);
    best_choose_cnt.assign(G, 0);

    suffix_weight.assign(G + 1, 0);
    suffix_value.assign(G + 1, 0);

    for (int i = G - 1; i >= 0; --i) {
      suffix_weight[i] = suffix_weight[i + 1] + (ll)groups[i].ids.size() * groups[i].w;
      suffix_value[i] = suffix_value[i + 1] + (i128)groups[i].ids.size() * (i128)groups[i].v;
    }

    best_value = forced_value;
    best_weight = forced_weight;

    dfs(0, CAP, forced_value);

    res.value = to_ll(best_value);
    res.weight = best_weight;

    for (int i = 0; i < G; ++i) {
      int k = best_choose_cnt[i];
      for (int j = 0; j < k; ++j) {
        res.taken[groups[i].ids[j]] = 1;
      }
    }

    return res;
  }
};

/*
int main(){
    ll N,W;
    cin>>N>>W;
    vector<ll> V(N),WW(N);
    for(ll i=0;i<N;i++){
        cin>>V[i]>>WW[i];
    }
    BranchAndBoundKnapsack BB;
    ll ans=BB.run(V,WW,W).value;
    cout<<ans<<endl;
}
*/
