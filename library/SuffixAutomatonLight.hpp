#include <bits/stdc++.h>
using namespace std;

struct SuffixAutomatonLower {
    static constexpr int SIGMA = 26;

    struct Node {
        int next[SIGMA];
        int link;
        int len;
        int end_cnt;    // 非 clone 状態なら 1、clone なら 0
        int first_pos;  // この状態の最長文字列の、ある出現の終端位置

        Node() : link(-1), len(0), end_cnt(0), first_pos(-1) {
            memset(next, -1, sizeof(next));
        }
    };

    vector<Node> st;
    vector<int> order, cnt, occ_cache;
    int last = 0;
    string s;
    bool order_ready = false;
    bool occ_ready = false;

    // 説明: 英小文字を 0..25 に変換する。
    // 計算量: O(1)
    static inline int id(char ch) {
#ifndef NDEBUG
        assert('a' <= ch && ch <= 'z');
#endif
        return ch - 'a';
    }

    // 説明: 前計算キャッシュを無効化する。
    // 計算量: O(1)
    void invalidate() {
        order_ready = false;
        occ_ready = false;
    }

    // 説明: 空文字列用に初期化する。max_len を与えると 2*max_len 程度を reserve する。
    // 計算量: O(1) amortized
    explicit SuffixAutomatonLower(int max_len = 0) {
        init(max_len);
    }

    // 説明: automaton を空に初期化する。
    // 計算量: O(1) amortized
    void init(int max_len = 0) {
        st.clear();
        st.reserve(max(2, 2 * max_len));
        st.push_back(Node()); // root
        order.clear();
        cnt.clear();
        occ_cache.clear();
        last = 0;
        s.clear();
        invalidate();
    }

    // 説明: 末尾に 1 文字追加して online に更新する。英小文字専用。
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
        invalidate();
    }

    // 説明: 文字列 t 全体から suffix automaton を構築する。
    // 計算量: O(|t|)
    void build(const string& t) {
        init((int)t.size());
        for (char ch : t) extend(ch);
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

    // 説明: 文字列 t を読んだ先の状態番号を返す。存在しなければ -1。
    // 計算量: O(|t|)
    int state_of(string_view t) const {
        int v = 0;
        for (char ch : t) {
            v = st[v].next[id(ch)];
            if (v == -1) return -1;
        }
        return v;
    }

    // 説明: 文字列 t が元文字列の部分文字列かどうかを返す。
    // 計算量: O(|t|)
    bool contains(string_view t) const {
        return state_of(t) != -1;
    }

    // 説明: 文字列 t の、ある出現の開始位置を返す。存在しなければ -1。
    //        返り値は 0-indexed。
    // 計算量: O(|t|)
    int first_occurrence(string_view t) const {
        int v = state_of(t);
        if (v == -1) return -1;
        return st[v].first_pos - (int)t.size() + 1;
    }

    // 説明: 状態を len 昇順に並べる。
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

    // 説明: 各状態の endpos size を計算する。
    //        これにより、各状態に属する文字列の出現回数が得られる。
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

    // 説明: 文字列 t の出現回数を返す。存在しなければ 0。
    // 計算量: O(|t|) + 初回のみ O(|S| + states)
    int count_occurrences(string_view t) {
        int v = state_of(t);
        if (v == -1) return 0;
        prepare_occ();
        return occ_cache[v];
    }

    // 説明: 異なる部分文字列数を返す。
    //        公式 sum_v (len[v] - len[link[v]]) を使う。
    // 計算量: O(states)
    long long distinct_substrings() const {
        long long ans = 0;
        for (int v = 1; v < (int)st.size(); ++v) {
            ans += st[v].len - st[st[v].link].len;
        }
        return ans;
    }

    // 説明: 文字列 t との最長共通部分文字列の位置を返す。
    //        返り値は {開始位置, 長さ} で、位置は t 側の 0-indexed。
    //        共通部分文字列が無ければ {-1, 0}。
    // 計算量: O(|t|)
    pair<int, int> longest_common_substring_pos(string_view t) const {
        int v = 0;
        int l = 0;
        int best_len = 0;
        int best_r = -1;

        for (int i = 0; i < (int)t.size(); ++i) {
            const int c = id(t[i]);

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

    // 説明: min_occ 回以上出る最長部分文字列を返す。
    //        既定値 min_occ = 2 なので、通常の最長反復部分文字列になる。
    // 計算量: O(states) + 初回のみ O(|S| + states)
    string longest_substring_with_min_occ(int min_occ = 2) {
        auto [l, len] = longest_substring_with_min_occ_pos(min_occ);
        if (len == 0) return {};
        return s.substr(l, len);
    }

    // 説明: 最長反復部分文字列の位置を返す。
    //        返り値は {開始位置, 長さ}。存在しなければ {-1, 0}。
    // 計算量: O(states) + 初回のみ O(|S| + states)
    pair<int, int> longest_repeated_substring_pos() {
        return longest_substring_with_min_occ_pos(2);
    }

    // 説明: 最長反復部分文字列を返す。
    // 計算量: O(states) + 初回のみ O(|S| + states)
    string longest_repeated_substring() {
        return longest_substring_with_min_occ(2);
    }
};
