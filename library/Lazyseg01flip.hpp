#include <bits/stdc++.h>
using namespace std;

struct LazySegTree01 {
    using ll = long long;

    static constexpr int ID   = 0;
    static constexpr int SET0 = 1;
    static constexpr int SET1 = 2;
    static constexpr int FLIP = 3;

    int n;
    vector<ll> ones;   // 区間内の 1 の個数
    vector<int> lazy;  // 遅延作用

    LazySegTree01() : n(0) {}

    LazySegTree01(const vector<int>& a) {
        int sz = (int)a.size();
        n = 1;
        while (n < sz) n <<= 1;
        ones.assign(2 * n, 0);
        lazy.assign(2 * n, ID);

        for (int i = 0; i < sz; i++) ones[n + i] = a[i];
        for (int i = n - 1; i >= 1; i--) {
            ones[i] = ones[i << 1] + ones[i << 1 | 1];
        }
    }

    // g を先に適用し、その後 f を適用する合成 f∘g
    int compose(int f, int g) {
        if (f == ID) return g;
        if (f == SET0) return SET0;
        if (f == SET1) return SET1;
        // f == FLIP
        if (g == ID)   return FLIP;
        if (g == SET0) return SET1;
        if (g == SET1) return SET0;
        return ID; // FLIP ∘ FLIP = ID
    }

    void apply_node(int k, int len, int op) {
        if (op == ID) return;
        if (op == SET0) ones[k] = 0;
        else if (op == SET1) ones[k] = len;
        else if (op == FLIP) ones[k] = len - ones[k];
        lazy[k] = compose(op, lazy[k]);
    }

    void push(int k, int len) {
        if (lazy[k] == ID || k >= n) return;
        int half = len >> 1;
        apply_node(k << 1,     half, lazy[k]);
        apply_node(k << 1 | 1, half, lazy[k]);
        lazy[k] = ID;
    }

    void pull(int k) {
        ones[k] = ones[k << 1] + ones[k << 1 | 1];
    }

    void range_assign(int a, int b, int v) {
        range_apply(a, b, v ? SET1 : SET0, 1, 0, n);
    }

    void range_flip(int a, int b) {
        range_apply(a, b, FLIP, 1, 0, n);
    }

    void range_apply(int a, int b, int op, int k, int l, int r) {
        if (r <= a || b <= l) return;
        if (a <= l && r <= b) {
            apply_node(k, r - l, op);
            return;
        }
        push(k, r - l);
        int m = (l + r) >> 1;
        range_apply(a, b, op, k << 1, l, m);
        range_apply(a, b, op, k << 1 | 1, m, r);
        pull(k);
    }

    ll range_sum(int a, int b) {
        return range_sum(a, b, 1, 0, n);
    }

    ll range_sum(int a, int b, int k, int l, int r) {
        if (r <= a || b <= l) return 0;
        if (a <= l && r <= b) return ones[k];
        push(k, r - l);
        int m = (l + r) >> 1;
        return range_sum(a, b, k << 1, l, m)
             + range_sum(a, b, k << 1 | 1, m, r);
    }

    int get(int p) {
        return (int)range_sum(p, p + 1);
    }
};

/*
int main() {
    vector<int> a = {1, 0, 1, 1, 0, 0, 1, 0};
    LazySegTree01 seg(a);

    cout << seg.range_sum(0, 8) << '\n'; // 4

    seg.range_flip(2, 7);                // [2,7) を反転
    cout << seg.range_sum(0, 8) << '\n'; // 5

    seg.range_assign(1, 5, 0);           // [1,5) を全部 0
    cout << seg.range_sum(0, 8) << '\n'; // 2

    seg.range_assign(3, 8, 1);           // [3,8) を全部 1
    cout << seg.range_sum(0, 8) << '\n'; // 6

    cout << seg.get(4) << '\n';          // 1
}
*/
