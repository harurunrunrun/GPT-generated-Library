#include <bits/stdc++.h>
using namespace std;

// ============================================================
// AtCoder でよく使うモノイド
// ============================================================
namespace atcoder_monoid {

template<class T>
struct Add {
    using S = T;
    static S id() { return S(0); }
    static S op(const S& a, const S& b) { return a + b; }
};

template<class T>
struct Mul {
    using S = T;
    static S id() { return S(1); }
    static S op(const S& a, const S& b) { return a * b; }
};

template<class T, T INF>
struct Min {
    using S = T;
    static S id() { return INF; }
    static S op(const S& a, const S& b) { return min(a, b); }
};

template<class T, T NEG_INF>
struct Max {
    using S = T;
    static S id() { return NEG_INF; }
    static S op(const S& a, const S& b) { return max(a, b); }
};

template<class T>
struct BitXor {
    using S = T;
    static S id() { return S(0); }
    static S op(const S& a, const S& b) { return a ^ b; }
};

template<class T>
struct BitOr {
    using S = T;
    static S id() { return S(0); }
    static S op(const S& a, const S& b) { return a | b; }
};

template<class T>
struct BitAnd {
    using S = T;
    static S id() { return ~S(0); }
    static S op(const S& a, const S& b) { return a & b; }
};

template<class T>
struct Gcd {
    using S = T;
    static S id() { return S(0); }
    static S op(const S& a, const S& b) { return std::gcd(a, b); }
};

template<class T>
struct Lcm {
    using S = T;
    static S id() { return S(1); }
    static S op(const S& a, const S& b) { return std::lcm(a, b); }
};

// すぐ使える using
using SumLL = Add<long long>;
using SumInt = Add<int>;

using ProdLL = Mul<long long>;
using ProdInt = Mul<int>;

using MinLL = Min<long long, (1LL << 62)>;
using MaxLL = Max<long long, -(1LL << 62)>;

using MinInt = Min<int, (1 << 30)>;
using MaxInt = Max<int, -(1 << 30)>;

using XorLL = BitXor<long long>;
using XorInt = BitXor<int>;

using OrLL = BitOr<long long>;
using OrInt = BitOr<int>;

using AndLL = BitAnd<long long>;
using AndInt = BitAnd<int>;

using GcdLL = Gcd<long long>;
using GcdInt = Gcd<int>;

using LcmLL = Lcm<long long>;
using LcmInt = Lcm<int>;

} // namespace atcoder_monoid

// ============================================================
// 圧縮区間を 1 ノードとして持つ implicit treap
// - 各ノードは [連続した同値区間] を表す
// - 区間全体のモノイド積を保持
// - 隣接して同値なら自動マージ
//
// 要件:
// - Monoid は
//     using S = ...;
//     static S id();
//     static S op(const S&, const S&);
//   を持つこと
// - S は隣接マージ判定のため == 可能であること
// ============================================================
template<class Monoid, class Eq = std::equal_to<typename Monoid::S>>
struct IntervalTreap {
    using S = typename Monoid::S;
    using i64 = long long;

    struct Node {
        Node *l, *r;
        uint32_t pri;
        i64 len;   // このノードが表す区間長
        i64 sz;    // 部分木全体の長さ
        int segs;  // 部分木内の圧縮区間数
        S val;     // この区間の値
        S run;     // val を len 回並べたモノイド積
        S prod;    // 部分木全体のモノイド積

        Node(i64 len_, const S& v, uint32_t pri_)
            : l(nullptr), r(nullptr), pri(pri_),
              len(len_), sz(len_), segs(1), val(v),
              run(Monoid::id()), prod(Monoid::id()) {}
    };

    Node* root = nullptr;

    // ------------------------------------------------------------
    // O(1)
    // 空の木を作る
    // ------------------------------------------------------------
    IntervalTreap() = default;

    // ------------------------------------------------------------
    // 期待 O(1)
    // 全体 [0, n) を値 init で初期化する
    // ------------------------------------------------------------
    IntervalTreap(i64 n, const S& init) {
        if (n > 0) root = new_node(n, init);
    }

    // ------------------------------------------------------------
    // 期待 O(K log K), K = 圧縮後区間数, worst O(N log N)
    // 配列から構築する
    // 隣接して同値な部分は自動で 1 区間に圧縮する
    // ------------------------------------------------------------
    explicit IntervalTreap(const vector<S>& a) {
        build_from_vector(a);
    }

    IntervalTreap(const IntervalTreap&) = delete;
    IntervalTreap& operator=(const IntervalTreap&) = delete;

    // ------------------------------------------------------------
    // O(1)
    // ムーブ構築
    // ------------------------------------------------------------
    IntervalTreap(IntervalTreap&& other) noexcept : root(other.root) {
        other.root = nullptr;
    }

    // ------------------------------------------------------------
    // O(M)
    // ムーブ代入
    // ------------------------------------------------------------
    IntervalTreap& operator=(IntervalTreap&& other) noexcept {
        if (this != &other) {
            clear(root);
            root = other.root;
            other.root = nullptr;
        }
        return *this;
    }

    // ------------------------------------------------------------
    // O(M)
    // デストラクタ
    // ------------------------------------------------------------
    ~IntervalTreap() {
        clear(root);
    }

    // ------------------------------------------------------------
    // O(1)
    // 全体長を返す
    // ------------------------------------------------------------
    i64 size() const {
        return total_len(root);
    }

    // ------------------------------------------------------------
    // O(1)
    // 圧縮後の区間数を返す
    // ------------------------------------------------------------
    int segment_count() const {
        return total_segs(root);
    }

    // ------------------------------------------------------------
    // O(1)
    // 空かどうか
    // ------------------------------------------------------------
    bool empty() const {
        return root == nullptr;
    }

    // ------------------------------------------------------------
    // O(M)
    // 全体を破棄して空にする
    // ------------------------------------------------------------
    void clear() {
        clear(root);
        root = nullptr;
    }

    // ------------------------------------------------------------
    // 期待 O(K log K), K = 圧縮後区間数
    // 配列から再構築する
    // ------------------------------------------------------------
    void build_from_vector(const vector<S>& a) {
        clear(root);
        root = nullptr;
        if (a.empty()) return;
        i64 l = 0;
        for (i64 i = 1; i <= (i64)a.size(); i++) {
            if (i == (i64)a.size() || !same_value(a[l], a[i])) {
                root = merge_adj(root, new_node(i - l, a[l]));
                l = i;
            }
        }
    }

    // ------------------------------------------------------------
    // 期待 O(log M)
    // 末尾に長さ len, 値 v の区間を追加する
    // ------------------------------------------------------------
    void push_back(i64 len, const S& v) {
        if (len <= 0) return;
        root = merge_adj(root, new_node(len, v));
    }

    // ------------------------------------------------------------
    // 期待 O(log M)
    // pos の位置に長さ len, 値 v の区間を挿入する
    // ------------------------------------------------------------
    void insert(i64 pos, i64 len, const S& v) {
        assert(0 <= pos && pos <= size());
        if (len <= 0) return;
        Node *a, *b;
        split_by_pos(root, pos, a, b);
        root = merge_adj(merge_adj(a, new_node(len, v)), b);
    }

    // ------------------------------------------------------------
    // 期待 O(log M)
    // [l, r) を削除する
    // ------------------------------------------------------------
    void erase(i64 l, i64 r) {
        assert(0 <= l && l <= r && r <= size());
        if (l == r) return;
        Node *a, *b, *c;
        split3(root, l, r, a, b, c);
        clear(b);
        root = merge_adj(a, c);
    }

    // ------------------------------------------------------------
    // 期待 O(log M)
    // [l, r) を値 v に代入する
    // ------------------------------------------------------------
    void assign(i64 l, i64 r, const S& v) {
        assert(0 <= l && l <= r && r <= size());
        if (l == r) return;
        Node *a, *b, *c;
        split3(root, l, r, a, b, c);
        clear(b);
        root = merge_adj(merge_adj(a, new_node(r - l, v)), c);
    }

    // ------------------------------------------------------------
    // 期待 O(log M)
    // [l, r) のモノイド積を返す
    // 非可換モノイドでも左から右の順で正しく積を取る
    // ------------------------------------------------------------
    S prod(i64 l, i64 r) {
        assert(0 <= l && l <= r && r <= size());
        if (l == r) return Monoid::id();
        Node *a, *b, *c;
        split3(root, l, r, a, b, c);
        S res = prod_of(b);
        root = merge_adj(merge_adj(a, b), c);
        return res;
    }

    // ------------------------------------------------------------
    // O(1)
    // 全体のモノイド積を返す
    // ------------------------------------------------------------
    S prod_all() const {
        return prod_of(root);
    }

    // ------------------------------------------------------------
    // 期待 O(log M)
    // 点 pos の値を返す
    // ------------------------------------------------------------
    S get(i64 pos) const {
        assert(0 <= pos && pos < size());
        Node* t = root;
        while (t) {
            i64 ls = total_len(t->l);
            if (pos < ls) {
                t = t->l;
            } else if (pos < ls + t->len) {
                return t->val;
            } else {
                pos -= ls + t->len;
                t = t->r;
            }
        }
        assert(false);
        return Monoid::id();
    }

    // ------------------------------------------------------------
    // 期待 O(log M)
    // pos を含む圧縮区間 [l, r), value を返す
    // ------------------------------------------------------------
    tuple<i64, i64, S> segment_at(i64 pos) const {
        assert(0 <= pos && pos < size());
        Node* t = root;
        i64 base = 0;
        while (t) {
            i64 ls = total_len(t->l);
            if (pos < ls) {
                t = t->l;
            } else if (pos < ls + t->len) {
                i64 l = base + ls;
                i64 r = l + t->len;
                return {l, r, t->val};
            } else {
                base += ls + t->len;
                pos -= ls + t->len;
                t = t->r;
            }
        }
        assert(false);
        return {0, 0, Monoid::id()};
    }

    // ------------------------------------------------------------
    // 期待 O(log M)
    // [l, r) がちょうど 1 つの圧縮区間からなるか判定する
    // ------------------------------------------------------------
    bool same(i64 l, i64 r) {
        assert(0 <= l && l <= r && r <= size());
        if (l == r) return true;
        Node *a, *b, *c;
        split3(root, l, r, a, b, c);
        bool ok = (total_segs(b) == 1);
        root = merge_adj(merge_adj(a, b), c);
        return ok;
    }

    // ------------------------------------------------------------
    // 期待 O(log M)
    // [l, r) に含まれる圧縮区間数を返す
    // ------------------------------------------------------------
    int count_segments(i64 l, i64 r) {
        assert(0 <= l && l <= r && r <= size());
        if (l == r) return 0;
        Node *a, *b, *c;
        split3(root, l, r, a, b, c);
        int cnt = total_segs(b);
        root = merge_adj(merge_adj(a, b), c);
        return cnt;
    }

    // ------------------------------------------------------------
    // 期待 O(k + log M)
    // [l, r) 内の各圧縮区間に対して
    // f(seg_l, seg_r, value) を呼ぶ
    // k = 走査する圧縮区間数
    // ------------------------------------------------------------
    template<class F>
    void for_each_segment(i64 l, i64 r, F f) {
        assert(0 <= l && l <= r && r <= size());
        if (l == r) return;
        Node *a, *b, *c;
        split3(root, l, r, a, b, c);
        i64 pos = l;
        dfs_for_each(b, pos, f);
        root = merge_adj(merge_adj(a, b), c);
    }

    // ------------------------------------------------------------
    // 期待 O((k + 1) log M + k log k)
    // [l, r) 内の各圧縮区間の値を変換する
    // f(len, value_ref) を各区間に適用する
    // 長さは変更しない
    // k = 変換対象の圧縮区間数
    // ------------------------------------------------------------
    template<class F>
    void transform_segments(i64 l, i64 r, F f) {
        assert(0 <= l && l <= r && r <= size());
        if (l == r) return;
        Node *a, *b, *c;
        split3(root, l, r, a, b, c);

        vector<pair<i64, S>> runs;
        runs.reserve(max(1, total_segs(b)));
        collect_transformed_runs(b, runs, f);
        clear(b);
        b = build_from_runs(runs);

        root = merge_adj(merge_adj(a, b), c);
    }

    // ------------------------------------------------------------
    // 期待 O(log M * log L)
    // pred(prod[l, x)) == true を満たす最大の x を返す
    // pred は単調で、pred(id()) == true を仮定する
    // L は 1 区間の長さの最大値程度
    // ------------------------------------------------------------
    template<class F>
    i64 max_right(i64 l, F pred) {
        assert(0 <= l && l <= size());
        assert(pred(Monoid::id()));
        Node *a, *b;
        split_by_pos(root, l, a, b);
        S cur = Monoid::id();
        i64 add = max_right_sub(b, cur, pred);
        root = merge_adj(a, b);
        return l + add;
    }

    // ------------------------------------------------------------
    // 期待 O(log M * log L)
    // pred(prod[x, r)) == true を満たす最小の x を返す
    // pred は単調で、pred(id()) == true を仮定する
    // L は 1 区間の長さの最大値程度
    // ------------------------------------------------------------
    template<class F>
    i64 min_left(i64 r, F pred) {
        assert(0 <= r && r <= size());
        assert(pred(Monoid::id()));
        Node *a, *b;
        split_by_pos(root, r, a, b);
        S cur = Monoid::id();
        i64 take = min_left_sub(a, cur, pred);
        root = merge_adj(a, b);
        return r - take;
    }

    // ------------------------------------------------------------
    // 期待 O(log M)
    // 位置 k で分割し、後半を返す
    // このオブジェクトには前半 [0, k) が残る
    // ------------------------------------------------------------
    IntervalTreap split_off(i64 k) {
        assert(0 <= k && k <= size());
        Node *a, *b;
        split_by_pos(root, k, a, b);
        root = a;
        IntervalTreap res;
        res.root = b;
        return res;
    }

    // ------------------------------------------------------------
    // 期待 O(log M)
    // 他の木を末尾に連結する
    // 境界の値が同じなら自動でマージする
    // ------------------------------------------------------------
    void concat(IntervalTreap&& other) {
        root = merge_adj(root, other.root);
        other.root = nullptr;
    }

    // ------------------------------------------------------------
    // O(M)
    // デバッグ用
    // 全圧縮区間を {l, r, value} の列として返す
    // ------------------------------------------------------------
    vector<tuple<i64, i64, S>> dump() const {
        vector<tuple<i64, i64, S>> out;
        out.reserve(segment_count());
        i64 pos = 0;
        dfs_dump(root, pos, out);
        return out;
    }

private:
    // ------------------------------------------------------------
    // O(1)
    // 疑似乱数 priority を生成する
    // ------------------------------------------------------------
    static uint32_t rng32() {
        static uint64_t x = 0x9e3779b97f4a7c15ULL;
        x += 0x9e3779b97f4a7c15ULL;
        uint64_t z = x;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        z = z ^ (z >> 31);
        return (uint32_t)(z >> 32);
    }

    // ------------------------------------------------------------
    // O(1)
    // 値の等値判定
    // ------------------------------------------------------------
    static bool same_value(const S& a, const S& b) {
        return Eq{}(a, b);
    }

    // ------------------------------------------------------------
    // O(1)
    // 部分木の全長
    // ------------------------------------------------------------
    static i64 total_len(Node* t) {
        return t ? t->sz : 0;
    }

    // ------------------------------------------------------------
    // O(1)
    // 部分木の圧縮区間数
    // ------------------------------------------------------------
    static int total_segs(Node* t) {
        return t ? t->segs : 0;
    }

    // ------------------------------------------------------------
    // O(1)
    // 部分木のモノイド積
    // ------------------------------------------------------------
    static S prod_of(Node* t) {
        return t ? t->prod : Monoid::id();
    }

    // ------------------------------------------------------------
    // O(log k)
    // x を k 回並べたモノイド積を二分累乗で返す
    // ------------------------------------------------------------
    static S power(S x, i64 k) {
        S res = Monoid::id();
        while (k > 0) {
            if (k & 1) res = Monoid::op(res, x);
            x = Monoid::op(x, x);
            k >>= 1;
        }
        return res;
    }

    // ------------------------------------------------------------
    // O(size of subtree)
    // 部分木を再帰的に解放する
    // ------------------------------------------------------------
    static void clear(Node* t) {
        if (!t) return;
        clear(t->l);
        clear(t->r);
        delete t;
    }

    // ------------------------------------------------------------
    // 期待 O(1)
    // 新しい 1 区間ノードを作る
    // ------------------------------------------------------------
    static Node* new_node(i64 len, const S& v) {
        Node* t = new Node(len, v, rng32());
        update_run(t);
        pull(t);
        return t;
    }

    // ------------------------------------------------------------
    // O(log len)
    // ノード自身の run を再計算する
    // ------------------------------------------------------------
    static void update_run(Node* t) {
        t->run = power(t->val, t->len);
    }

    // ------------------------------------------------------------
    // O(1)
    // 子から集約情報を再計算する
    // ------------------------------------------------------------
    static void pull(Node* t) {
        if (!t) return;
        t->sz = total_len(t->l) + t->len + total_len(t->r);
        t->segs = total_segs(t->l) + 1 + total_segs(t->r);
        t->prod = Monoid::op(Monoid::op(prod_of(t->l), t->run), prod_of(t->r));
    }

    // ------------------------------------------------------------
    // 期待 O(log M)
    // 通常の treap merge
    // 隣接区間の値が同じでもここではマージしない
    // ------------------------------------------------------------
    static Node* raw_merge(Node* a, Node* b) {
        if (!a || !b) return a ? a : b;
        if (a->pri < b->pri) {
            a->r = raw_merge(a->r, b);
            pull(a);
            return a;
        } else {
            b->l = raw_merge(a, b->l);
            pull(b);
            return b;
        }
    }

    // ------------------------------------------------------------
    // 期待 O(log M)
    // 最左の 1 区間を切り出して first に入れる
    // 残りは rest に入れる
    // ------------------------------------------------------------
    static void split_first_seg(Node* t, Node*& first, Node*& rest) {
        if (!t) {
            first = rest = nullptr;
            return;
        }
        if (!t->l) {
            first = t;
            rest = t->r;
            first->r = nullptr;
            pull(first);
            return;
        }
        split_first_seg(t->l, first, t->l);
        pull(t);
        rest = t;
    }

    // ------------------------------------------------------------
    // 期待 O(log M)
    // 最右の 1 区間を切り出して last に入れる
    // 残りは rest に入れる
    // ------------------------------------------------------------
    static void split_last_seg(Node* t, Node*& rest, Node*& last) {
        if (!t) {
            rest = last = nullptr;
            return;
        }
        if (!t->r) {
            last = t;
            rest = t->l;
            last->l = nullptr;
            pull(last);
            return;
        }
        split_last_seg(t->r, t->r, last);
        pull(t);
        rest = t;
    }

    // ------------------------------------------------------------
    // 期待 O(log M)
    // a と b を連結する
    // 境界の値が同じなら 1 区間にマージする
    // ------------------------------------------------------------
    static Node* merge_adj(Node* a, Node* b) {
        if (!a || !b) return a ? a : b;

        Node *a0, *x, *y, *b0;
        split_last_seg(a, a0, x);
        split_first_seg(b, y, b0);

        if (same_value(x->val, y->val)) {
            x->len += y->len;
            update_run(x);
            pull(x);
            delete y;
            return raw_merge(raw_merge(a0, x), b0);
        } else {
            return raw_merge(raw_merge(a0, x), raw_merge(y, b0));
        }
    }

    // ------------------------------------------------------------
    // 期待 O(log M)
    // 長さ k で [0, k), [k, n) に split する
    // 必要なら 1 ノードを途中で 2 つに分割する
    // ------------------------------------------------------------
    static void split_by_pos(Node* t, i64 k, Node*& a, Node*& b) {
        if (!t) {
            a = b = nullptr;
            return;
        }

        i64 ls = total_len(t->l);

        if (k < ls) {
            split_by_pos(t->l, k, a, t->l);
            pull(t);
            b = t;
            return;
        }

        if (k > ls + t->len) {
            split_by_pos(t->r, k - ls - t->len, t->r, b);
            pull(t);
            a = t;
            return;
        }

        if (k == ls) {
            a = t->l;
            t->l = nullptr;
            pull(t);
            b = t;
            return;
        }

        if (k == ls + t->len) {
            b = t->r;
            t->r = nullptr;
            pull(t);
            a = t;
            return;
        }

        i64 left_len = k - ls;
        i64 right_len = t->len - left_len;

        Node* L = new_node(left_len, t->val);
        Node* R = new_node(right_len, t->val);

        a = raw_merge(t->l, L);
        b = raw_merge(R, t->r);

        t->l = t->r = nullptr;
        delete t;
    }

    // ------------------------------------------------------------
    // 期待 O(log M)
    // [0, l), [l, r), [r, n) に 3 分割する
    // ------------------------------------------------------------
    static void split3(Node* t, i64 l, i64 r, Node*& a, Node*& b, Node*& c) {
        Node* ab;
        split_by_pos(t, r, ab, c);
        split_by_pos(ab, l, a, b);
    }

    // ------------------------------------------------------------
    // 期待 O(K log K)
    // {len, value} の列から木を構築する
    // 隣接同値はマージする
    // ------------------------------------------------------------
    static Node* build_from_runs(const vector<pair<i64, S>>& runs) {
        Node* t = nullptr;
        for (auto& [len, v] : runs) {
            if (len <= 0) continue;
            t = merge_adj(t, new_node(len, v));
        }
        return t;
    }

    // ------------------------------------------------------------
    // O(k)
    // 部分木内の各圧縮区間に対して f(l, r, value) を呼ぶ
    // pos は開始位置として使い、走査後に末尾位置になる
    // ------------------------------------------------------------
    template<class F>
    static void dfs_for_each(Node* t, i64& pos, F& f) {
        if (!t) return;
        dfs_for_each(t->l, pos, f);
        i64 l = pos;
        i64 r = l + t->len;
        f(l, r, t->val);
        pos = r;
        dfs_for_each(t->r, pos, f);
    }

    // ------------------------------------------------------------
    // O(k)
    // 部分木内の各圧縮区間を走査し、変換後の run 列を作る
    // 隣接同値はここでマージする
    // ------------------------------------------------------------
    template<class F>
    static void collect_transformed_runs(Node* t, vector<pair<i64, S>>& runs, F& f) {
        if (!t) return;
        collect_transformed_runs(t->l, runs, f);

        S v = t->val;
        f(t->len, v);

        if (!runs.empty() && same_value(runs.back().second, v)) {
            runs.back().first += t->len;
        } else {
            runs.push_back({t->len, v});
        }

        collect_transformed_runs(t->r, runs, f);
    }

    // ------------------------------------------------------------
    // O(M)
    // デバッグ用に全圧縮区間をベクタへ出力する
    // ------------------------------------------------------------
    static void dfs_dump(Node* t, i64& pos, vector<tuple<i64, i64, S>>& out) {
        if (!t) return;
        dfs_dump(t->l, pos, out);
        i64 l = pos, r = l + t->len;
        out.emplace_back(l, r, t->val);
        pos = r;
        dfs_dump(t->r, pos, out);
    }

    // ------------------------------------------------------------
    // O(log len * log len)
    // 1 区間 val^len の prefix のうち、
    // pred(cur * val^x) を満たす最大 x を返す
    // ------------------------------------------------------------
    template<class F>
    static i64 max_prefix_in_run(const S& cur, const S& val, i64 len, F& pred) {
        i64 ok = 0, ng = len + 1;
        while (ng - ok > 1) {
            i64 mid = (ok + ng) >> 1;
            if (pred(Monoid::op(cur, power(val, mid)))) ok = mid;
            else ng = mid;
        }
        return ok;
    }

    // ------------------------------------------------------------
    // O(log len * log len)
    // 1 区間 val^len の suffix のうち、
    // pred(val^x * cur) を満たす最大 x を返す
    // ------------------------------------------------------------
    template<class F>
    static i64 max_suffix_in_run(const S& val, i64 len, const S& cur, F& pred) {
        i64 ok = 0, ng = len + 1;
        while (ng - ok > 1) {
            i64 mid = (ok + ng) >> 1;
            if (pred(Monoid::op(power(val, mid), cur))) ok = mid;
            else ng = mid;
        }
        return ok;
    }

    // ------------------------------------------------------------
    // 期待 O(log M * log L)
    // max_right 用の内部 DFS
    // cur はこれまでの積
    // 返り値はこの部分木から何文字分取れるか
    // ------------------------------------------------------------
    template<class F>
    static i64 max_right_sub(Node* t, S& cur, F& pred) {
        if (!t) return 0;

        if (pred(Monoid::op(cur, prod_of(t)))) {
            cur = Monoid::op(cur, prod_of(t));
            return total_len(t);
        }

        i64 res = 0;

        if (t->l) {
            if (!pred(Monoid::op(cur, prod_of(t->l)))) {
                return max_right_sub(t->l, cur, pred);
            }
            cur = Monoid::op(cur, prod_of(t->l));
            res += total_len(t->l);
        }

        i64 take = max_prefix_in_run(cur, t->val, t->len, pred);
        cur = Monoid::op(cur, power(t->val, take));
        res += take;
        if (take < t->len) return res;

        return res + max_right_sub(t->r, cur, pred);
    }

    // ------------------------------------------------------------
    // 期待 O(log M * log L)
    // min_left 用の内部 DFS
    // cur はこれまで右側に積んだ積
    // 返り値はこの部分木から左へ何文字分取れるか
    // ------------------------------------------------------------
    template<class F>
    static i64 min_left_sub(Node* t, S& cur, F& pred) {
        if (!t) return 0;

        if (pred(Monoid::op(prod_of(t), cur))) {
            cur = Monoid::op(prod_of(t), cur);
            return total_len(t);
        }

        i64 res = 0;

        if (t->r) {
            if (!pred(Monoid::op(prod_of(t->r), cur))) {
                return min_left_sub(t->r, cur, pred);
            }
            cur = Monoid::op(prod_of(t->r), cur);
            res += total_len(t->r);
        }

        i64 take = max_suffix_in_run(t->val, t->len, cur, pred);
        cur = Monoid::op(power(t->val, take), cur);
        res += take;
        if (take < t->len) return res;

        return res + min_left_sub(t->l, cur, pred);
    }
};

/*
使い方例:

using Seg = IntervalTreap<atcoder_monoid::SumLL>;

int main() {
    Seg tr(10, 0);
    tr.assign(2, 7, 5);
    tr.insert(4, 3, 1);

    cout << tr.prod(0, tr.size()) << '\n';
    cout << tr.get(5) << '\n';

    auto [l, r, v] = tr.segment_at(5);
    cout << l << ' ' << r << ' ' << v << '\n';

    tr.for_each_segment(0, tr.size(), [&](long long l, long long r, long long v) {
        cout << "[" << l << "," << r << ") = " << v << '\n';
    });

    long long x = tr.max_right(0, [&](long long s) {
        return s <= 12;
    });
    cout << x << '\n';
}
*/
