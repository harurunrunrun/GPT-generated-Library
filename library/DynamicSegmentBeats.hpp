#include <bits/stdc++.h>
using namespace std;

struct DynamicSegTreeBeats {
    static constexpr ll INF = (1LL << 60);

    struct Node {
        Node *l = nullptr, *r = nullptr;

        ll sum = 0;

        ll max_v = 0, smax_v = -INF;
        int max_c = 1;

        ll min_v = 0, smin_v = INF;
        int min_c = 1;

        ll add = 0;
    };

    ll L, R;   // 管理区間 [L, R)
    Node* root;

    /*
        概要:
            管理区間 [L, R) を持つ動的 Segment Tree Beats を構築する。
            初期値は全て 0。

        時間計算量:
            O(1)
    */
    DynamicSegTreeBeats(ll L_, ll R_) : L(L_), R(R_) {
        root = new Node();
    }

    /*
        概要:
            新しいノードを 1 つ生成する。
            初期状態は「区間内の値が全て 0」であることを表す。

        時間計算量:
            O(1)
    */
    static Node* make_node() {
        return new Node();
    }

    /*
        概要:
            区間長 r - l を返す。

        時間計算量:
            O(1)
    */
    static ll seg_len(ll l, ll r) {
        return r - l;
    }

    /*
        概要:
            ノード t が表す区間全体に +x を作用させる。
            sum / min / max / 2番目の min,max をまとめて更新する。

        時間計算量:
            O(1)
    */
    static void push_add(Node* t, ll x, ll sz) {
        if (!t) return;
        t->sum += x * sz;
        t->max_v += x;
        if (t->smax_v != -INF) t->smax_v += x;
        t->min_v += x;
        if (t->smin_v != INF) t->smin_v += x;
        t->add += x;
    }

    /*
        概要:
            ノード t が表す区間全体に chmin(x) を直接適用する。
            ただし、呼び出し側で
                t->smax_v < x < t->max_v
            またはそれに準ずる適用可能条件が満たされていることを前提とする。

        時間計算量:
            O(1)
    */
    static void push_chmin(Node* t, ll x) {
        if (!t || t->max_v <= x) return;

        t->sum += (x - t->max_v) * 1LL * t->max_c;

        if (t->min_v == t->max_v) {
            t->min_v = x;
        } else if (t->smin_v == t->max_v) {
            t->smin_v = x;
        }
        t->max_v = x;
    }

    /*
        概要:
            ノード t が表す区間全体に chmax(x) を直接適用する。
            ただし、呼び出し側で
                t->min_v < x < t->smin_v
            またはそれに準ずる適用可能条件が満たされていることを前提とする。

        時間計算量:
            O(1)
    */
    static void push_chmax(Node* t, ll x) {
        if (!t || t->min_v >= x) return;

        t->sum += (x - t->min_v) * 1LL * t->min_c;

        if (t->max_v == t->min_v) {
            t->max_v = x;
        } else if (t->smax_v == t->min_v) {
            t->smax_v = x;
        }
        t->min_v = x;
    }

    /*
        概要:
            子 2 つの情報から親ノード t の情報を再計算する。

        時間計算量:
            O(1)
    */
    static void pull(Node* t) {
        Node *a = t->l, *b = t->r;

        t->sum = a->sum + b->sum;

        // max 系
        if (a->max_v > b->max_v) {
            t->max_v = a->max_v;
            t->max_c = a->max_c;
            t->smax_v = max(a->smax_v, b->max_v);
        } else if (a->max_v < b->max_v) {
            t->max_v = b->max_v;
            t->max_c = b->max_c;
            t->smax_v = max(a->max_v, b->smax_v);
        } else {
            t->max_v = a->max_v;
            t->max_c = a->max_c + b->max_c;
            t->smax_v = max(a->smax_v, b->smax_v);
        }

        // min 系
        if (a->min_v < b->min_v) {
            t->min_v = a->min_v;
            t->min_c = a->min_c;
            t->smin_v = min(a->smin_v, b->min_v);
        } else if (a->min_v > b->min_v) {
            t->min_v = b->min_v;
            t->min_c = b->min_c;
            t->smin_v = min(a->min_v, b->smin_v);
        } else {
            t->min_v = a->min_v;
            t->min_c = a->min_c + b->min_c;
            t->smin_v = min(a->smin_v, b->smin_v);
        }

        t->add = 0;
    }

    /*
        概要:
            ノード t の子ノードを必要なら生成し、
            さらに親に載っている遅延情報(add, chmin/chmax の影響)を子へ伝播する。

        時間計算量:
            O(1)
    */
    void ensure_children(Node* t, ll l, ll r) {
        if (r - l == 1) return;

        ll m = (l + r) >> 1;
        if (!t->l) t->l = make_node();
        if (!t->r) t->r = make_node();

        ll lsz = m - l;
        ll rsz = r - m;

        if (t->add != 0) {
            push_add(t->l, t->add, lsz);
            push_add(t->r, t->add, rsz);
        }

        if (t->l->max_v > t->max_v) push_chmin(t->l, t->max_v);
        if (t->r->max_v > t->max_v) push_chmin(t->r, t->max_v);

        if (t->l->min_v < t->min_v) push_chmax(t->l, t->min_v);
        if (t->r->min_v < t->min_v) push_chmax(t->r, t->min_v);

        t->add = 0;
    }

    /*
        概要:
            区間 [ql, qr) に +x を加える。

        時間計算量:
            償却 O(log(R-L))
    */
    void range_add(ll ql, ll qr, ll x) {
        range_add(root, L, R, ql, qr, x);
    }

    /*
        概要:
            区間 [ql, qr) に対して a[i] = min(a[i], x) を適用する。

        時間計算量:
            償却 O(log(R-L))
    */
    void range_chmin(ll ql, ll qr, ll x) {
        range_chmin(root, L, R, ql, qr, x);
    }

    /*
        概要:
            区間 [ql, qr) に対して a[i] = max(a[i], x) を適用する。

        時間計算量:
            償却 O(log(R-L))
    */
    void range_chmax(ll ql, ll qr, ll x) {
        range_chmax(root, L, R, ql, qr, x);
    }

    /*
        概要:
            区間 [ql, qr) の総和を返す。

        時間計算量:
            O(log(R-L))
    */
    ll range_sum(ll ql, ll qr) {
        return range_sum(root, L, R, ql, qr);
    }

    /*
        概要:
            区間 [ql, qr) の最小値を返す。

        時間計算量:
            O(log(R-L))
    */
    ll range_min(ll ql, ll qr) {
        return range_min(root, L, R, ql, qr);
    }

    /*
        概要:
            区間 [ql, qr) の最大値を返す。

        時間計算量:
            O(log(R-L))
    */
    ll range_max(ll ql, ll qr) {
        return range_max(root, L, R, ql, qr);
    }

private:
    /*
        概要:
            内部実装。区間 [ql, qr) に +x を加える。

        時間計算量:
            償却 O(log(R-L))
    */
    void range_add(Node* t, ll l, ll r, ll ql, ll qr, ll x) {
        if (!t || qr <= l || r <= ql) return;

        if (ql <= l && r <= qr) {
            push_add(t, x, seg_len(l, r));
            return;
        }

        ensure_children(t, l, r);
        ll m = (l + r) >> 1;
        range_add(t->l, l, m, ql, qr, x);
        range_add(t->r, m, r, ql, qr, x);
        pull(t);
    }

    /*
        概要:
            内部実装。区間 [ql, qr) に chmin(x) を適用する。

        時間計算量:
            償却 O(log(R-L))
    */
    void range_chmin(Node* t, ll l, ll r, ll ql, ll qr, ll x) {
        if (!t || qr <= l || r <= ql || t->max_v <= x) return;

        if (ql <= l && r <= qr && t->smax_v < x) {
            push_chmin(t, x);
            return;
        }

        ensure_children(t, l, r);
        ll m = (l + r) >> 1;
        range_chmin(t->l, l, m, ql, qr, x);
        range_chmin(t->r, m, r, ql, qr, x);
        pull(t);
    }

    /*
        概要:
            内部実装。区間 [ql, qr) に chmax(x) を適用する。

        時間計算量:
            償却 O(log(R-L))
    */
    void range_chmax(Node* t, ll l, ll r, ll ql, ll qr, ll x) {
        if (!t || qr <= l || r <= ql || t->min_v >= x) return;

        if (ql <= l && r <= qr && t->smin_v > x) {
            push_chmax(t, x);
            return;
        }

        ensure_children(t, l, r);
        ll m = (l + r) >> 1;
        range_chmax(t->l, l, m, ql, qr, x);
        range_chmax(t->r, m, r, ql, qr, x);
        pull(t);
    }

    /*
        概要:
            内部実装。区間 [ql, qr) の総和を返す。

        時間計算量:
            O(log(R-L))
    */
    ll range_sum(Node* t, ll l, ll r, ll ql, ll qr) {
        if (!t || qr <= l || r <= ql) return 0;

        if (ql <= l && r <= qr) return t->sum;

        ensure_children(t, l, r);
        ll m = (l + r) >> 1;
        return range_sum(t->l, l, m, ql, qr) + range_sum(t->r, m, r, ql, qr);
    }

    /*
        概要:
            内部実装。区間 [ql, qr) の最小値を返す。

        時間計算量:
            O(log(R-L))
    */
    ll range_min(Node* t, ll l, ll r, ll ql, ll qr) {
        if (!t || qr <= l || r <= ql) return INF;

        if (ql <= l && r <= qr) return t->min_v;

        ensure_children(t, l, r);
        ll m = (l + r) >> 1;
        return min(range_min(t->l, l, m, ql, qr), range_min(t->r, m, r, ql, qr));
    }

    /*
        概要:
            内部実装。区間 [ql, qr) の最大値を返す。

        時間計算量:
            O(log(R-L))
    */
    ll range_max(Node* t, ll l, ll r, ll ql, ll qr) {
        if (!t || qr <= l || r <= ql) return -INF;

        if (ql <= l && r <= qr) return t->max_v;

        ensure_children(t, l, r);
        ll m = (l + r) >> 1;
        return max(range_max(t->l, l, m, ql, qr), range_max(t->r, m, r, ql, qr));
    }
};

int main() {
    DynamicSegTreeBeats seg(0, 1000000000LL);

    seg.range_add(2, 10, 5);   // [2, 10) に +5
    seg.range_chmin(0, 7, 3);  // [0, 7) を min(., 3)
    seg.range_chmax(5, 12, 4); // [5, 12) を max(., 4)

    cout << seg.range_sum(0, 12) << '\n';
    cout << seg.range_min(0, 12) << '\n';
    cout << seg.range_max(0, 12) << '\n';

    return 0;
}
