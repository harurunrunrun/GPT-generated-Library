// https://github.com/harurunrunrun/GPT-generated-Library

#include <bits/stdc++.h>
using namespace std;

// Dynamic (implicit) segtree: monoid (S, op, e)
// API: set/get/prod/all_prod/max_right/min_left  (ACL segtree 互換寄り)
template <class S, S (*op)(S, S), S (*e)()>
struct dynamic_segtree {
    struct Node {
        S val;
        Node *l, *r;
        explicit Node(const S& v) : val(v), l(nullptr), r(nullptr) {}
    };

    long long n;     // range: [0, n)
    Node* root;

    dynamic_segtree() : n(0), root(nullptr) {}
    explicit dynamic_segtree(long long n_) : n(n_), root(nullptr) {
        assert(n >= 0);
    }

    // ----- basic -----
    S all_prod() const { return root ? root->val : e(); }

    S get(long long p) const {
        assert(0 <= p && p < n);
        return prod(p, p + 1);
    }

    void set(long long p, const S& x) {
        assert(0 <= p && p < n);
        root = set_impl(root, 0, n, p, x);
    }

    S prod(long long l, long long r) const {
        assert(0 <= l && l <= r && r <= n);
        return prod_impl(root, 0, n, l, r);
    }

    // ----- binary search on segtree (ACL style) -----
    // max_right(l, f): 最大の r (l<=r<=n) で f(prod(l,r)) が true
    // ただし f(e()) == true を要求
    template <class F>
    long long max_right(long long l, F f) const {
        assert(0 <= l && l <= n);
        S sm = e();
        if (!f(sm)) return l;
        return max_right_impl(root, 0, n, l, f, sm);
    }

    // min_left(r, f): 最小の l (0<=l<=r) で f(prod(l,r)) が true
    // ただし f(e()) == true を要求
    template <class F>
    long long min_left(long long r, F f) const {
        assert(0 <= r && r <= n);
        S sm = e();
        if (!f(sm)) return r;
        return min_left_impl(root, 0, n, r, f, sm);
    }

private:
    static S node_val(Node* t) { return t ? t->val : e(); }

    Node* set_impl(Node* t, long long nl, long long nr, long long p, const S& x) {
        if (!t) t = new Node(e());
        if (nr - nl == 1) {
            t->val = x;
            return t;
        }
        long long mid = nl + (nr - nl) / 2;
        if (p < mid) t->l = set_impl(t->l, nl, mid, p, x);
        else         t->r = set_impl(t->r, mid, nr, p, x);
        t->val = op(node_val(t->l), node_val(t->r));
        return t;
    }

    S prod_impl(Node* t, long long nl, long long nr, long long ql, long long qr) const {
        if (qr <= nl || nr <= ql) return e();
        if (ql <= nl && nr <= qr) return node_val(t);
        long long mid = nl + (nr - nl) / 2;
        S lv = prod_impl(t ? t->l : nullptr, nl, mid, ql, qr);
        S rv = prod_impl(t ? t->r : nullptr, mid, nr, ql, qr);
        return op(lv, rv);
    }

    template <class F>
    long long max_right_impl(Node* t, long long nl, long long nr, long long ql, F f, S& sm) const {
        if (nr <= ql) return ql;

        if (ql <= nl) {
            S nxt = op(sm, node_val(t));
            if (f(nxt)) {
                sm = nxt;
                return nr;
            }
            if (nr - nl == 1) return nl;
        }

        long long mid = nl + (nr - nl) / 2;
        long long res = max_right_impl(t ? t->l : nullptr, nl, mid, ql, f, sm);
        if (res < mid) return res;
        return max_right_impl(t ? t->r : nullptr, mid, nr, ql, f, sm);
    }

    template <class F>
    long long min_left_impl(Node* t, long long nl, long long nr, long long qr, F f, S& sm) const {
        if (qr <= nl) return qr;

        if (nr <= qr) {
            S nxt = op(node_val(t), sm);
            if (f(nxt)) {
                sm = nxt;
                return nl;
            }
            if (nr - nl == 1) return nr;
        }

        long long mid = nl + (nr - nl) / 2;
        long long res = min_left_impl(t ? t->r : nullptr, mid, nr, qr, f, sm);
        if (mid < res) return res;
        return min_left_impl(t ? t->l : nullptr, nl, mid, qr, f, sm);
    }
};
