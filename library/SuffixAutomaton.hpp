#include <bits/stdc++.h>
using namespace std;

struct SuffixAutomatonLower {
    static constexpr int SIGMA = 26;

    struct Node {
        int next[SIGMA];
        int link;
        int len;
        int end_cnt;    // 非 clone 状態なら 1、clone なら 0
        int first_pos;  // この状態に属する文字列のある出現の終端位置
        long long sub;  // この状態から始まる異なる部分文字列数
        int miss_len;   // この状態からの最短未出現文字列長

        Node() : link(-1), len(0), end_cnt(0), first_pos(-1), sub(0), miss_len(1) {
            memset(next, -1, sizeof(next));
        }
    };

    vector<Node> st;
    vector<int> order;
    vector<int> cnt;
    vector<int> occ_cache;
    vector<unsigned char> terminal;
    vector<int> prefix_state;

    vector<long long> exact_distinct_by_len;
    vector<long long> pref_distinct_cnt;
    vector<long long> pref_distinct_len;

    vector<long long> multi_sub; // 多重度つき部分文字列総数 DP
    vector<int> min_end_pos;
    vector<int> max_end_pos;

    vector<vector<int>> link_children;
    vector<int> tin, tout, euler;

    int last = 0;
    string s;

    bool order_ready = false;
    bool occ_ready = false;
    bool sub_ready = false;
    bool term_ready = false;
    bool exact_len_ready = false;
    bool miss_ready = false;
    bool multi_ready = false;
    bool minmax_ready = false;
    bool link_tree_ready = false;

    // 説明: 英小文字を 0..25 に変換する。
    // 計算量: O(1)
    static inline int id(char ch) {
#ifndef NDEBUG
        assert('a' <= ch && ch <= 'z');
#endif
        return ch - 'a';
    }

    // 説明: 各種前計算キャッシュを無効化する。
    // 計算量: O(1)
    void invalidate() {
        order_ready = false;
        occ_ready = false;
        sub_ready = false;
        term_ready = false;
        exact_len_ready = false;
        miss_ready = false;
        multi_ready = false;
        minmax_ready = false;
        link_tree_ready = false;
    }

    // 説明: 空文字列用の automaton を作る。max_len を与えると 2*max_len 程度を reserve する。
    // 計算量: O(1) amortized
    explicit SuffixAutomatonLower(int max_len = 0) {
        init(max_len);
    }

    // 説明: automaton を空に初期化する。
    // 計算量: O(1) amortized
    void init(int max_len = 0) {
        st.clear();
        st.reserve(max(2, 2 * max_len));
        st.push_back(Node());

        order.clear();
        cnt.clear();
        occ_cache.clear();
        terminal.clear();
        prefix_state.clear();

        exact_distinct_by_len.clear();
        pref_distinct_cnt.clear();
        pref_distinct_len.clear();

        multi_sub.clear();
        min_end_pos.clear();
        max_end_pos.clear();

        link_children.clear();
        tin.clear();
        tout.clear();
        euler.clear();

        last = 0;
        s.clear();
        invalidate();
    }

    // 説明: 文字列 t 全体から suffix automaton を構築する。
    // 計算量: O(|t|)
    void build(const string& t) {
        init((int)t.size());
        for (char ch : t) extend(ch);
    }

    // 説明: 末尾に 1 文字追加してオンライン更新する。英小文字専用。
    // 計算量: O(1) amortized
    void extend(char ch) {
        const int c = id(ch);
        const int cur = (int)st.size();
        st.push_back(Node());
        st[cur].len = st[last].len + 1;
        st[cur].end_cnt = 1;
        st[cur].first_pos = st[cur].len - 1;

        int p = last;
        while (p != -1 && st[p].next[c] == -1) {
            st[p].next[c] = cur;
            p = st[p].link;
        }

        if (p == -1) {
            st[cur].link = 0;
        } else {
            const int q = st[p].next[c];
            if (st[p].len + 1 == st[q].len) {
                st[cur].link = q;
            } else {
                const int clone = (int)st.size();
                st.push_back(st[q]);
                st[clone].len = st[p].len + 1;
                st[clone].end_cnt = 0;
                while (p != -1 && st[p].next[c] == q) {
                    st[p].next[c] = clone;
                    p = st[p].link;
                }
                st[q].link = clone;
                st[cur].link = clone;
            }
        }

        last = cur;
        s.push_back(ch);
        prefix_state.push_back(last);
        invalidate();
    }

    // 説明: 元文字列長を返す。
    // 計算量: O(1)
    int length() const {
        return (int)s.size();
    }

    // 説明: 状態数を返す。
    // 計算量: O(1)
    int states() const {
        return (int)st.size();
    }

    // 説明: root の状態番号を返す。
    // 計算量: O(1)
    int root() const {
        return 0;
    }

    // 説明: 元文字列を返す。
    // 計算量: O(1)
    const string& original_string() const {
        return s;
    }

    // 説明: 状態 v が clone かどうかを返す。
    // 計算量: O(1)
    bool is_clone_state(int v) const {
        return st[v].end_cnt == 0;
    }

    // 説明: prefix s[0..i] を読んだ状態番号を返す。範囲外なら -1。
    // 計算量: O(1)
    int prefix_state_at(int i) const {
        if (i < 0 || i >= (int)prefix_state.size()) return -1;
        return prefix_state[i];
    }

    // 説明: 各状態を len 昇順に並べる。
    // 計算量: O(|S| + states)
    void build_order() {
        if (order_ready) return;

        const int max_len = s.empty() ? 0 : st[last].len;
        cnt.assign(max_len + 1, 0);
        for (const auto& node : st) ++cnt[node.len];
        for (int i = 1; i <= max_len; ++i) cnt[i] += cnt[i - 1];

        order.assign(st.size(), 0);
        for (int i = (int)st.size() - 1; i >= 0; --i) {
            order[--cnt[st[i].len]] = i;
        }

        order_ready = true;
    }

    // 説明: suffix-link 木を構築し、Euler tour の tin/tout も計算する。
    // 計算量: O(states)
    void prepare_link_tree() {
        if (link_tree_ready) return;

        link_children.assign(st.size(), {});
        for (int v = 1; v < (int)st.size(); ++v) {
            link_children[st[v].link].push_back(v);
        }

        tin.assign(st.size(), -1);
        tout.assign(st.size(), -1);
        euler.clear();
        euler.reserve(st.size());

        vector<pair<int, int>> stack;
        stack.reserve(st.size());

        tin[0] = 0;
        euler.push_back(0);
        stack.push_back({0, 0});

        while (!stack.empty()) {
            int v = stack.back().first;
            int& it = stack.back().second;
            if (it < (int)link_children[v].size()) {
                int u = link_children[v][it++];
                tin[u] = (int)euler.size();
                euler.push_back(u);
                stack.push_back({u, 0});
            } else {
                tout[v] = (int)euler.size();
                stack.pop_back();
            }
        }

        link_tree_ready = true;
    }

    // 説明: 各状態の endpos size を計算する。
    // 計算量: O(|S| + states)
    void prepare_occ() {
        if (occ_ready) return;
        build_order();

        occ_cache.assign(st.size(), 0);
        for (int i = 0; i < (int)st.size(); ++i) occ_cache[i] = st[i].end_cnt;

        for (int i = (int)order.size() - 1; i > 0; --i) {
            int v = order[i];
            occ_cache[st[v].link] += occ_cache[v];
        }

        occ_ready = true;
    }

    // 説明: 各状態から始まる異なる部分文字列数 sub を計算する。
    // 計算量: O(SIGMA * states)
    void prepare_sub() {
        if (sub_ready) return;
        build_order();

        for (int i = (int)order.size() - 1; i >= 0; --i) {
            int v = order[i];
            long long ways = 0;
            for (int c = 0; c < SIGMA; ++c) {
                int u = st[v].next[c];
                if (u != -1) ways += 1 + st[u].sub;
            }
            st[v].sub = ways;
        }

        sub_ready = true;
    }

    // 説明: suffix を表す terminal 状態を計算する。
    // 計算量: O(|S|)
    void prepare_terminal() {
        if (term_ready) return;
        terminal.assign(st.size(), 0);
        for (int v = last; v != -1; v = st[v].link) terminal[v] = 1;
        term_ready = true;
    }

    // 説明: 長さ L ごとの異なる部分文字列数と、その prefix sum を計算する。
    //        exact_distinct_by_len[L] = 長さ L の異なる部分文字列数
    // 計算量: O(|S| + states)
    void prepare_exact_distinct_by_len() {
        if (exact_len_ready) return;

        int n = length();
        exact_distinct_by_len.assign(n + 1, 0);
        pref_distinct_cnt.assign(n + 1, 0);
        pref_distinct_len.assign(n + 1, 0);

        vector<long long> diff(n + 2, 0);

        for (int v = 1; v < (int)st.size(); ++v) {
            int lo = st[st[v].link].len + 1;
            int hi = st[v].len;
            ++diff[lo];
            --diff[hi + 1];
        }

        for (int L = 1; L <= n; ++L) {
            diff[L] += diff[L - 1];
            exact_distinct_by_len[L] = diff[L];
            pref_distinct_cnt[L] = pref_distinct_cnt[L - 1] + exact_distinct_by_len[L];
            pref_distinct_len[L] = pref_distinct_len[L - 1] + 1LL * L * exact_distinct_by_len[L];
        }

        exact_len_ready = true;
    }

    // 説明: 各状態からの最短未出現文字列長 miss_len を計算する。
    // 計算量: O(SIGMA * states)
    void prepare_missing() {
        if (miss_ready) return;
        build_order();

        for (int i = (int)order.size() - 1; i >= 0; --i) {
            int v = order[i];

            bool has_missing = false;
            for (int c = 0; c < SIGMA; ++c) {
                if (st[v].next[c] == -1) {
                    has_missing = true;
                    break;
                }
            }

            if (has_missing) {
                st[v].miss_len = 1;
            } else {
                int best = INT_MAX;
                for (int c = 0; c < SIGMA; ++c) {
                    int u = st[v].next[c];
                    best = min(best, 1 + st[u].miss_len);
                }
                st[v].miss_len = best;
            }
        }

        miss_ready = true;
    }

    // 説明: 多重度つき部分文字列総数 DP を計算する。
    //        multi_sub[v] = 状態 v から始まる全部分文字列出現回数の総和
    // 計算量: O(SIGMA * states) + 初回のみ O(|S| + states)
    void prepare_multi_sub() {
        if (multi_ready) return;
        prepare_occ();
        build_order();

        multi_sub.assign(st.size(), 0);
        for (int i = (int)order.size() - 1; i >= 0; --i) {
            int v = order[i];
            long long ways = 0;
            for (int c = 0; c < SIGMA; ++c) {
                int u = st[v].next[c];
                if (u != -1) ways += occ_cache[u] + multi_sub[u];
            }
            multi_sub[v] = ways;
        }

        multi_ready = true;
    }

    // 説明: 各状態の endpos 集合の最小終端位置・最大終端位置を計算する。
    //        非重複反復部分文字列の判定に使う。
    // 計算量: O(|S| + states)
    void prepare_minmax_end() {
        if (minmax_ready) return;
        build_order();

        const int INF = 1e9;
        min_end_pos.assign(st.size(), INF);
        max_end_pos.assign(st.size(), -INF);

        for (int v = 0; v < (int)st.size(); ++v) {
            if (st[v].end_cnt > 0) {
                min_end_pos[v] = st[v].first_pos;
                max_end_pos[v] = st[v].first_pos;
            }
        }

        for (int i = (int)order.size() - 1; i > 0; --i) {
            int v = order[i];
            int p = st[v].link;
            min_end_pos[p] = min(min_end_pos[p], min_end_pos[v]);
            max_end_pos[p] = max(max_end_pos[p], max_end_pos[v]);
        }

        minmax_ready = true;
    }

    // 説明: 長さ制約を [1, |S|] に丸める。
    // 計算量: O(1)
    pair<int, int> normalize_length_range(int L, int R) const {
        L = max(L, 1);
        R = min(R, length());
        if (L > R) return {1, 0};
        return {L, R};
    }

    // 説明: 文字列 t を読んだ先の状態番号を返す。存在しなければ -1。
    // 計算量: O(|t|)
    int state_of(string_view t) const {
        int v = 0;
        for (char ch : t) {
            int c = id(ch);
            v = st[v].next[c];
            if (v == -1) return -1;
        }
        return v;
    }

    // 説明: 状態 v から文字 ch で遷移した先を返す。無ければ -1。
    // 計算量: O(1)
    int transition(int v, char ch) const {
        return st[v].next[id(ch)];
    }

    // 説明: 状態 v の suffix link を返す。root では -1。
    // 計算量: O(1)
    int suffix_link(int v) const {
        return st[v].link;
    }

    // 説明: 状態 v が表す部分文字列長の最大値を返す。
    // 計算量: O(1)
    int state_max_length(int v) const {
        return st[v].len;
    }

    // 説明: 状態 v が表す部分文字列長の最小値を返す。
    // 計算量: O(1)
    int state_min_length(int v) const {
        return (v == 0 ? 0 : st[st[v].link].len + 1);
    }

    // 説明: 状態 v が表す長さ範囲 [min_len, max_len] を返す。
    // 計算量: O(1)
    pair<int, int> state_length_range(int v) const {
        return {state_min_length(v), state_max_length(v)};
    }

    // 説明: 状態 v の代表となる最長文字列を返す。root なら空文字列。
    // 計算量: O(state_max_length(v))
    string state_longest_string(int v) const {
        if (v == 0) return {};
        int len = st[v].len;
        int l = st[v].first_pos - len + 1;
        return s.substr(l, len);
    }

    // 説明: 状態 v の代表となる最短文字列を返す。root なら空文字列。
    // 計算量: O(state_min_length(v))
    string state_shortest_string(int v) const {
        if (v == 0) return {};
        int len = state_min_length(v);
        int l = st[v].first_pos - len + 1;
        return s.substr(l, len);
    }

    // 説明: 状態 v に属する長さ want_len の代表文字列を返す。範囲外なら空文字列。
    // 計算量: O(want_len)
    string state_string_of_length(int v, int want_len) const {
        auto [mn, mx] = state_length_range(v);
        if (want_len < mn || want_len > mx) return {};
        int l = st[v].first_pos - want_len + 1;
        return s.substr(l, want_len);
    }

    // 説明: 文字列 t が部分文字列として存在するかを返す。
    // 計算量: O(|t|)
    bool contains(string_view t) const {
        return state_of(t) != -1;
    }

    // 説明: 文字列 t が suffix かどうかを返す。
    // 計算量: O(|t|) + 初回のみ O(|S|)
    bool is_suffix(string_view t) {
        int v = state_of(t);
        if (v == -1) return false;
        prepare_terminal();
        return terminal[v] != 0;
    }

    // 説明: suffix を表す状態番号を長い順に返す。
    // 計算量: O(|S|)
    vector<int> all_suffix_states_descending_length() const {
        vector<int> res;
        for (int v = last; v != 0; v = st[v].link) res.push_back(v);
        return res;
    }

    // 説明: 状態 v の出現回数を返す。
    //        状態 v に属する任意の長さの文字列は同じ出現回数を持つ。
    // 計算量: O(1) + 初回のみ O(|S| + states)
    int state_occurrences(int v) {
        prepare_occ();
        return occ_cache[v];
    }

    // 説明: 文字列 t の出現回数を返す。存在しなければ 0。
    // 計算量: O(|t|) + 初回のみ O(|S| + states)
    int count_occurrences(string_view t) {
        int v = state_of(t);
        if (v == -1) return 0;
        prepare_occ();
        return occ_cache[v];
    }

    // 説明: 状態 v に属する文字列のある出現の終端位置を返す。
    // 計算量: O(1)
    int state_first_end_position(int v) const {
        return st[v].first_pos;
    }

    // 説明: 文字列 t の最初の出現開始位置を返す。無ければ -1。
    //        0-indexed。
    // 計算量: O(|t|)
    int first_occurrence(string_view t) const {
        int v = state_of(t);
        if (v == -1) return -1;
        return st[v].first_pos - (int)t.size() + 1;
    }

    // 説明: 異なる部分文字列数を返す。
    // 計算量: O(1) + 初回のみ O(SIGMA * states)
    long long distinct_substrings() {
        prepare_sub();
        return st[0].sub;
    }

    // 説明: 長さちょうど L の異なる部分文字列数を返す。範囲外なら 0。
    // 計算量: O(1) + 初回のみ O(|S| + states)
    long long distinct_substrings_exact_length(int L) {
        if (L < 0 || L > length()) return 0;
        prepare_exact_distinct_by_len();
        return exact_distinct_by_len[L];
    }

    // 説明: 長さ区間 [L, R] の異なる部分文字列数を返す。
    // 計算量: O(1) + 初回のみ O(|S| + states)
    long long distinct_substrings_length_range(int L, int R) {
        prepare_exact_distinct_by_len();
        tie(L, R) = normalize_length_range(L, R);
        if (L > R) return 0;
        return pref_distinct_cnt[R] - pref_distinct_cnt[L - 1];
    }

    // 説明: 長さごとの異なる部分文字列数を返す。res[L] = 長さ L の異なる部分文字列数。
    // 計算量: O(1) + 初回のみ O(|S| + states)
    const vector<long long>& all_distinct_substrings_exact_length() {
        prepare_exact_distinct_by_len();
        return exact_distinct_by_len;
    }

    // 説明: 異なる部分文字列の総長
    //        sum_{u distinct} |u|
    // を返す。
    // 計算量: O(states)
    long long total_length_of_distinct_substrings() const {
        long long ans = 0;
        for (int v = 1; v < (int)st.size(); ++v) {
            long long lo = st[st[v].link].len + 1;
            long long hi = st[v].len;
            ans += (lo + hi) * (hi - lo + 1) / 2;
        }
        return ans;
    }

    // 説明: 長さ区間 [L, R] の異なる部分文字列の総長
    //        sum_{u distinct, L<=|u|<=R} |u|
    // を返す。
    // 計算量: O(1) + 初回のみ O(|S| + states)
    long long total_length_of_distinct_substrings_in_range(int L, int R) {
        prepare_exact_distinct_by_len();
        tie(L, R) = normalize_length_range(L, R);
        if (L > R) return 0;
        return pref_distinct_len[R] - pref_distinct_len[L - 1];
    }

    // 説明: 異なる部分文字列それぞれの出現回数の総和
    //        sum_{u distinct} occ(u)
    // を返す。
    // 計算量: O(1) + 初回のみ O(SIGMA * states) + O(|S| + states)
    long long total_occurrences_of_distinct_substrings() {
        prepare_multi_sub();
        return multi_sub[0];
    }

    // 説明: 異なる部分文字列それぞれについて (長さ × 出現回数) の総和
    //        sum_{u distinct} |u| * occ(u)
    // を返す。
    // 計算量: O(states) + 初回のみ O(|S| + states)
    long long total_length_of_distinct_substrings_with_multiplicity() {
        prepare_occ();
        long long ans = 0;
        for (int v = 1; v < (int)st.size(); ++v) {
            long long lo = st[st[v].link].len + 1;
            long long hi = st[v].len;
            long long sum_len = (lo + hi) * (hi - lo + 1) / 2;
            ans += sum_len * occ_cache[v];
        }
        return ans;
    }

    // 説明: min_occ 回以上出る異なる部分文字列数を返す。
    // 計算量: O(states) + 初回のみ O(|S| + states)
    long long count_distinct_substrings_with_min_occ(int min_occ) {
        prepare_occ();
        long long ans = 0;
        for (int v = 1; v < (int)st.size(); ++v) {
            if (occ_cache[v] >= min_occ) {
                ans += st[v].len - st[st[v].link].len;
            }
        }
        return ans;
    }

    // 説明: 出現回数がちょうど exact_occ の異なる部分文字列数を返す。
    // 計算量: O(states) + 初回のみ O(|S| + states)
    long long count_distinct_substrings_with_exact_occ(int exact_occ) {
        prepare_occ();
        long long ans = 0;
        for (int v = 1; v < (int)st.size(); ++v) {
            if (occ_cache[v] == exact_occ) {
                ans += st[v].len - st[st[v].link].len;
            }
        }
        return ans;
    }

    // 説明: 出現回数が [low_occ, high_occ] に入る異なる部分文字列数を返す。
    // 計算量: O(states) + 初回のみ O(|S| + states)
    long long count_distinct_substrings_with_occ_range(int low_occ, int high_occ) {
        prepare_occ();
        long long ans = 0;
        for (int v = 1; v < (int)st.size(); ++v) {
            if (low_occ <= occ_cache[v] && occ_cache[v] <= high_occ) {
                ans += st[v].len - st[st[v].link].len;
            }
        }
        return ans;
    }

    // 説明: 出現回数が [low_occ, high_occ] に入る異なる部分文字列の総長を返す。
    // 計算量: O(states) + 初回のみ O(|S| + states)
    long long total_length_of_distinct_substrings_with_occ_range(int low_occ, int high_occ) {
        prepare_occ();
        long long ans = 0;
        for (int v = 1; v < (int)st.size(); ++v) {
            if (low_occ <= occ_cache[v] && occ_cache[v] <= high_occ) {
                long long lo = st[st[v].link].len + 1;
                long long hi = st[v].len;
                ans += (lo + hi) * (hi - lo + 1) / 2;
            }
        }
        return ans;
    }

    // 説明: 辞書順で k 番目の異なる部分文字列を返す。1-indexed。
    //        範囲外なら空文字列。
    // 計算量: O(SIGMA * |answer|) + 初回のみ O(SIGMA * states)
    string kth_distinct_substring(long long k) {
        prepare_sub();
        if (k <= 0 || k > st[0].sub) return {};

        string res;
        int v = 0;

        while (true) {
            for (int c = 0; c < SIGMA; ++c) {
                int u = st[v].next[c];
                if (u == -1) continue;

                long long block = 1 + st[u].sub;
                if (k > block) {
                    k -= block;
                    continue;
                }

                res.push_back(char('a' + c));
                if (k == 1) return res;
                --k;
                v = u;
                break;
            }
        }
    }

    // 説明: 辞書順で k 番目の「多重度つき」部分文字列を返す。1-indexed。
    //        同じ文字列も出現回数ぶん別要素として数える。
    //        範囲外なら空文字列。
    // 計算量: O(SIGMA * |answer|) + 初回のみ O(SIGMA * states) + O(|S| + states)
    string kth_substring_with_multiplicity(long long k) {
        prepare_multi_sub();
        if (k <= 0 || k > multi_sub[0]) return {};

        string res;
        int v = 0;

        while (true) {
            for (int c = 0; c < SIGMA; ++c) {
                int u = st[v].next[c];
                if (u == -1) continue;

                long long block = occ_cache[u] + multi_sub[u];
                if (k > block) {
                    k -= block;
                    continue;
                }

                res.push_back(char('a' + c));
                if (k <= occ_cache[u]) return res;
                k -= occ_cache[u];
                v = u;
                break;
            }
        }
    }

    // 説明: 辞書順で先頭 limit 個の異なる部分文字列を列挙する。
    // 計算量: O(出力総文字数 + limit * SIGMA)
    vector<string> first_k_distinct_substrings_lex(long long limit) const {
        vector<string> res;
        if (limit <= 0) return res;

        struct Frame {
            int v;
            int next_c;
        };

        vector<Frame> stack;
        vector<char> cur;
        stack.push_back({0, 0});

        while (!stack.empty()) {
            auto& fr = stack.back();
            bool moved = false;

            for (int c = fr.next_c; c < SIGMA; ++c) {
                fr.next_c = c + 1;
                int u = st[fr.v].next[c];
                if (u == -1) continue;

                cur.push_back(char('a' + c));
                res.emplace_back(cur.begin(), cur.end());
                if ((long long)res.size() >= limit) return res;

                stack.push_back({u, 0});
                moved = true;
                break;
            }

            if (!moved) {
                stack.pop_back();
                if (!cur.empty()) cur.pop_back();
            }
        }

        return res;
    }

    // 説明: 辞書順で先頭 limit 個の異なる部分文字列と、その出現回数を列挙する。
    // 計算量: O(出力総文字数 + limit * SIGMA) + 初回のみ O(|S| + states)
    vector<pair<string, int>> first_k_distinct_substrings_lex_with_occ(long long limit) {
        prepare_occ();

        vector<pair<string, int>> res;
        if (limit <= 0) return res;

        struct Frame {
            int v;
            int next_c;
        };

        vector<Frame> stack;
        vector<char> cur;
        stack.push_back({0, 0});

        while (!stack.empty()) {
            auto& fr = stack.back();
            bool moved = false;

            for (int c = fr.next_c; c < SIGMA; ++c) {
                fr.next_c = c + 1;
                int u = st[fr.v].next[c];
                if (u == -1) continue;

                cur.push_back(char('a' + c));
                res.push_back({string(cur.begin(), cur.end()), occ_cache[u]});
                if ((long long)res.size() >= limit) return res;

                stack.push_back({u, 0});
                moved = true;
                break;
            }

            if (!moved) {
                stack.pop_back();
                if (!cur.empty()) cur.pop_back();
            }
        }

        return res;
    }

    // 説明: 最短未出現部分文字列の長さを返す。
    // 計算量: O(1) + 初回のみ O(SIGMA * states)
    int shortest_absent_substring_length() {
        prepare_missing();
        return st[0].miss_len;
    }

    // 説明: 最短未出現部分文字列のうち、辞書順最小のものを返す。
    // 計算量: O(SIGMA * |answer|) + 初回のみ O(SIGMA * states)
    string shortest_absent_substring() {
        prepare_missing();

        string res;
        int v = 0;
        while (true) {
            if (st[v].miss_len == 1) {
                for (int c = 0; c < SIGMA; ++c) {
                    if (st[v].next[c] == -1) {
                        res.push_back(char('a' + c));
                        return res;
                    }
                }
            } else {
                for (int c = 0; c < SIGMA; ++c) {
                    int u = st[v].next[c];
                    if (u != -1 && 1 + st[u].miss_len == st[v].miss_len) {
                        res.push_back(char('a' + c));
                        v = u;
                        break;
                    }
                }
            }
        }
    }

    // 説明: 文字列 t との最長共通部分文字列の位置を返す。
    //        返り値は {開始位置, 長さ} で、位置は t 側の 0-indexed。
    //        存在しなければ {-1, 0}。
    // 計算量: O(|t|)
    pair<int, int> longest_common_substring_pos(string_view t) const {
        int v = 0;
        int l = 0;
        int best_len = 0;
        int best_r = -1;

        for (int i = 0; i < (int)t.size(); ++i) {
            int c = id(t[i]);

            while (v != 0 && st[v].next[c] == -1) {
                v = st[v].link;
                l = st[v].len;
            }

            if (st[v].next[c] != -1) {
                v = st[v].next[c];
                ++l;
            } else {
                v = 0;
                l = 0;
            }

            if (l > best_len) {
                best_len = l;
                best_r = i;
            }
        }

        if (best_len == 0) return {-1, 0};
        return {best_r - best_len + 1, best_len};
    }

    // 説明: 文字列 t との最長共通部分文字列そのものを返す。
    // 計算量: O(|t|)
    string longest_common_substring(string_view t) const {
        auto [l, len] = longest_common_substring_pos(t);
        if (len == 0) return {};
        return string(t.substr(l, len));
    }

    // 説明: 文字列 t について、各状態で「t と共通な最大長」を返す。
    //        返り値 mx[v] は状態 v に属する文字列のうち、t にも現れる最大長。
    // 計算量: O(|t| + states)
    vector<int> match_limit_single(string_view t) {
        build_order();

        vector<int> mx(st.size(), 0);
        int v = 0;
        int l = 0;

        for (char ch : t) {
            int c = id(ch);

            while (v != 0 && st[v].next[c] == -1) {
                v = st[v].link;
                l = st[v].len;
            }

            if (st[v].next[c] != -1) {
                v = st[v].next[c];
                ++l;
            } else {
                v = 0;
                l = 0;
            }

            mx[v] = max(mx[v], l);
        }

        for (int i = (int)order.size() - 1; i > 0; --i) {
            int x = order[i];
            int p = st[x].link;
            mx[p] = max(mx[p], min(mx[x], st[p].len));
        }

        return mx;
    }

    // 説明: 元文字列と t の両方に現れる異なる部分文字列数を返す。
    // 計算量: O(|t| + states)
    long long count_distinct_common_substrings_with(string_view t) {
        vector<int> mx = match_limit_single(t);

        long long ans = 0;
        for (int v = 1; v < (int)st.size(); ++v) {
            int lo = st[st[v].link].len + 1;
            int hi = min(st[v].len, mx[v]);
            if (hi >= lo) ans += hi - lo + 1;
        }
        return ans;
    }

    // 説明: 元文字列と t の両方に現れる異なる部分文字列の総長
    //        sum_{u distinct, u in S and u in t} |u|
    // を返す。
    // 計算量: O(|t| + states)
    long long total_length_of_distinct_common_substrings_with(string_view t) {
        vector<int> mx = match_limit_single(t);

        long long ans = 0;
        for (int v = 1; v < (int)st.size(); ++v) {
            long long lo = st[st[v].link].len + 1;
            long long hi = min(st[v].len, mx[v]);
            if (hi >= lo) ans += (lo + hi) * (hi - lo + 1) / 2;
        }
        return ans;
    }

    // 説明: 元文字列と others の全てに共通する最長部分文字列の位置を返す。
    //        返り値は base string 側の {開始位置, 長さ}。
    //        others が空なら全文字列を返す。
    // 計算量: O(states * |others| + 総文字数)
    pair<int, int> longest_common_substring_among_pos(const vector<string>& others) {
        if (others.empty()) return {0, length()};

        build_order();
        vector<int> mn(st.size(), 0);
        for (int v = 0; v < (int)st.size(); ++v) mn[v] = st[v].len;

        for (const string& t : others) {
            vector<int> mx = match_limit_single(t);
            for (int v = 0; v < (int)st.size(); ++v) {
                mn[v] = min(mn[v], mx[v]);
            }
        }

        int best_len = 0;
        int best_r = -1;
        for (int v = 1; v < (int)st.size(); ++v) {
            if (mn[v] > best_len) {
                best_len = mn[v];
                best_r = st[v].first_pos;
            }
        }

        if (best_len == 0) return {-1, 0};
        return {best_r - best_len + 1, best_len};
    }

    // 説明: 元文字列と others の全てに共通する最長部分文字列を返す。
    //        others が空なら全文字列を返す。
    // 計算量: O(states * |others| + 総文字数)
    string longest_common_substring_among(const vector<string>& others) {
        auto [l, len] = longest_common_substring_among_pos(others);
        if (len == 0) return {};
        return s.substr(l, len);
    }

    // 説明: 元文字列と others の全てに共通する異なる部分文字列数を返す。
    //        others が空なら全異なる部分文字列数を返す。
    // 計算量: O(states * |others| + 総文字数)
    long long count_distinct_common_substrings_among(const vector<string>& others) {
        if (others.empty()) return distinct_substrings();

        build_order();
        vector<int> mn(st.size(), 0);
        for (int v = 0; v < (int)st.size(); ++v) mn[v] = st[v].len;

        for (const string& t : others) {
            vector<int> mx = match_limit_single(t);
            for (int v = 0; v < (int)st.size(); ++v) {
                mn[v] = min(mn[v], mx[v]);
            }
        }

        long long ans = 0;
        for (int v = 1; v < (int)st.size(); ++v) {
            int lo = st[st[v].link].len + 1;
            int hi = mn[v];
            if (hi >= lo) ans += hi - lo + 1;
        }
        return ans;
    }

    // 説明: 元文字列と others の全てに共通する異なる部分文字列の総長を返す。
    //        others が空なら全異なる部分文字列総長を返す。
    // 計算量: O(states * |others| + 総文字数)
    long long total_length_of_distinct_common_substrings_among(const vector<string>& others) {
        if (others.empty()) return total_length_of_distinct_substrings();

        build_order();
        vector<int> mn(st.size(), 0);
        for (int v = 0; v < (int)st.size(); ++v) mn[v] = st[v].len;

        for (const string& t : others) {
            vector<int> mx = match_limit_single(t);
            for (int v = 0; v < (int)st.size(); ++v) {
                mn[v] = min(mn[v], mx[v]);
            }
        }

        long long ans = 0;
        for (int v = 1; v < (int)st.size(); ++v) {
            long long lo = st[st[v].link].len + 1;
            long long hi = mn[v];
            if (hi >= lo) ans += (lo + hi) * (hi - lo + 1) / 2;
        }
        return ans;
    }

    // 説明: min_occ 回以上出る最長部分文字列の位置を返す。
    //        返り値は {開始位置, 長さ}。存在しなければ {-1, 0}。
    // 計算量: O(states) + 初回のみ O(|S| + states)
    pair<int, int> longest_substring_with_min_occ_pos(int min_occ = 2) {
        prepare_occ();

        int best_len = 0;
        int best_r = -1;
        for (int v = 1; v < (int)st.size(); ++v) {
            if (occ_cache[v] >= min_occ && st[v].len > best_len) {
                best_len = st[v].len;
                best_r = st[v].first_pos;
            }
        }

        if (best_len == 0) return {-1, 0};
        return {best_r - best_len + 1, best_len};
    }

    // 説明: min_occ 回以上出る最長部分文字列そのものを返す。
    // 計算量: O(states) + 初回のみ O(|S| + states)
    string longest_substring_with_min_occ(int min_occ = 2) {
        auto [l, len] = longest_substring_with_min_occ_pos(min_occ);
        if (len == 0) return {};
        return s.substr(l, len);
    }

    // 説明: 出現回数がちょうど exact_occ の最長部分文字列の位置を返す。
    //        返り値は {開始位置, 長さ}。存在しなければ {-1, 0}。
    // 計算量: O(states) + 初回のみ O(|S| + states)
    pair<int, int> longest_substring_with_exact_occ_pos(int exact_occ) {
        prepare_occ();

        int best_len = 0;
        int best_r = -1;
        for (int v = 1; v < (int)st.size(); ++v) {
            if (occ_cache[v] == exact_occ && st[v].len > best_len) {
                best_len = st[v].len;
                best_r = st[v].first_pos;
            }
        }

        if (best_len == 0) return {-1, 0};
        return {best_r - best_len + 1, best_len};
    }

    // 説明: 出現回数がちょうど exact_occ の最長部分文字列を返す。
    // 計算量: O(states) + 初回のみ O(|S| + states)
    string longest_substring_with_exact_occ(int exact_occ) {
        auto [l, len] = longest_substring_with_exact_occ_pos(exact_occ);
        if (len == 0) return {};
        return s.substr(l, len);
    }

    // 説明: 出現回数が [low_occ, high_occ] に入る最長部分文字列の位置を返す。
    //        返り値は {開始位置, 長さ}。存在しなければ {-1, 0}。
    // 計算量: O(states) + 初回のみ O(|S| + states)
    pair<int, int> longest_substring_with_occ_range_pos(int low_occ, int high_occ) {
        prepare_occ();

        int best_len = 0;
        int best_r = -1;
        for (int v = 1; v < (int)st.size(); ++v) {
            if (low_occ <= occ_cache[v] && occ_cache[v] <= high_occ && st[v].len > best_len) {
                best_len = st[v].len;
                best_r = st[v].first_pos;
            }
        }

        if (best_len == 0) return {-1, 0};
        return {best_r - best_len + 1, best_len};
    }

    // 説明: 出現回数が [low_occ, high_occ] に入る最長部分文字列を返す。
    // 計算量: O(states) + 初回のみ O(|S| + states)
    string longest_substring_with_occ_range(int low_occ, int high_occ) {
        auto [l, len] = longest_substring_with_occ_range_pos(low_occ, high_occ);
        if (len == 0) return {};
        return s.substr(l, len);
    }

    // 説明: 通常の最長反復部分文字列の位置を返す。
    // 計算量: O(states) + 初回のみ O(|S| + states)
    pair<int, int> longest_repeated_substring_pos() {
        return longest_substring_with_min_occ_pos(2);
    }

    // 説明: 通常の最長反復部分文字列を返す。
    // 計算量: O(states) + 初回のみ O(|S| + states)
    string longest_repeated_substring() {
        return longest_substring_with_min_occ(2);
    }

    // 説明: 重なりなし最長反復部分文字列の位置を返す。
    //        同一文字列でも区間が重ならなければよい。
    //        返り値は {開始位置, 長さ}。存在しなければ {-1, 0}。
    // 計算量: O(states) + 初回のみ O(|S| + states)
    pair<int, int> longest_nonoverlapping_repeated_substring_pos() {
        prepare_minmax_end();

        int best_len = 0;
        int best_r = -1;

        for (int v = 1; v < (int)st.size(); ++v) {
            if (max_end_pos[v] < 0) continue;

            int cand = min(st[v].len, max_end_pos[v] - min_end_pos[v]);
            if (cand >= state_min_length(v) && cand > best_len) {
                best_len = cand;
                best_r = max_end_pos[v];
            }
        }

        if (best_len == 0) return {-1, 0};
        return {best_r - best_len + 1, best_len};
    }

    // 説明: 重なりなし最長反復部分文字列を返す。
    // 計算量: O(states) + 初回のみ O(|S| + states)
    string longest_nonoverlapping_repeated_substring() {
        auto [l, len] = longest_nonoverlapping_repeated_substring_pos();
        if (len == 0) return {};
        return s.substr(l, len);
    }

    // 説明: 状態 v の suffix-link 木部分木に含まれる非 clone 状態の end position を昇順で返す。
    //        occurrence の列挙用であり、集計系より重い。
    // 計算量: O(subtree size + occ log occ)
    vector<int> state_occurrence_end_positions(int v) {
        prepare_link_tree();

        vector<int> ends;
        for (int i = tin[v]; i < tout[v]; ++i) {
            int x = euler[i];
            if (st[x].end_cnt > 0) ends.push_back(st[x].first_pos);
        }

        sort(ends.begin(), ends.end());
        return ends;
    }

    // 説明: 状態 v に属する長さ pat_len の文字列の出現開始位置を昇順で返す。
    //        pat_len が範囲外なら空。
    // 計算量: O(subtree size + occ log occ)
    vector<int> state_occurrence_start_positions(int v, int pat_len) {
        auto [mn, mx] = state_length_range(v);
        if (pat_len < mn || pat_len > mx) return {};

        vector<int> ends = state_occurrence_end_positions(v);
        for (int& p : ends) p = p - pat_len + 1;
        return ends;
    }

    // 説明: 文字列 t の出現終端位置を昇順で返す。存在しなければ空。
    // 計算量: O(|t| + subtree size + occ log occ)
    vector<int> occurrence_end_positions(string_view t) {
        int v = state_of(t);
        if (v == -1) return {};
        return state_occurrence_end_positions(v);
    }

    // 説明: 文字列 t の出現開始位置を昇順で返す。存在しなければ空。
    //        返り値は 0-indexed。
    // 計算量: O(|t| + subtree size + occ log occ)
    vector<int> occurrence_start_positions(string_view t) {
        int v = state_of(t);
        if (v == -1) return {};
        return state_occurrence_start_positions(v, (int)t.size());
    }
};
