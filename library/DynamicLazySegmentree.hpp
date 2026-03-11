


/*
  dynamic lazy segtree
  - implicit / dynamic segment tree
  - 未訪問ノードは持たない
  - 区間長 len に応じた初期値 init(len) を返す
*/
template <class S, S (*op)(S, S), S (*e)(), S (*init)(long long), class F,
          S (*mapping)(F, S, long long), F (*composition)(F, F), F (*id)()>
struct dynamic_lazy_segtree {
  private:
    struct Node {
        S val;
        F lz;
        bool has_lz;
        Node* l;
        Node* r;
        explicit Node(S v) : val(v), lz(id()), has_lz(false), l(nullptr), r(nullptr) {}
    };

    long long _n;
    long long size;
    Node* root;

    long long seg_len(long long nl, long long nr) const {
        if (nl >= _n) return 0;
        return min(nr, _n) - nl;
    }

    S init_range(long long nl, long long nr) const {
        return init(seg_len(nl, nr));
    }

    Node* make_node(long long nl, long long nr) const {
        return new Node(init_range(nl, nr));
    }

    S node_val(Node* t, long long nl, long long nr) const {
        return t ? t->val : init_range(nl, nr);
    }

    void pull(Node* t, long long nl, long long nr) {
        long long mid = (nl + nr) >> 1;
        t->val = op(node_val(t->l, nl, mid), node_val(t->r, mid, nr));
    }

    void all_apply(Node* t, F f, long long nl, long long nr) {
        long long len = seg_len(nl, nr);
        t->val = mapping(f, t->val, len);
        if (nr - nl > 1 && len > 0) {
            if (t->has_lz) t->lz = composition(f, t->lz);
            else {
                t->lz = f;
                t->has_lz = true;
            }
        }
    }

    void push(Node* t, long long nl, long long nr) {
        if (!t || !t->has_lz || nr - nl == 1) return;
        long long mid = (nl + nr) >> 1;
        if (!t->l) t->l = make_node(nl, mid);
        if (!t->r) t->r = make_node(mid, nr);
        all_apply(t->l, t->lz, nl, mid);
        all_apply(t->r, t->lz, mid, nr);
        t->lz = id();
        t->has_lz = false;
    }

    void set(Node*& t, long long nl, long long nr, long long p, S x) {
        if (!t) t = make_node(nl, nr);
        if (nr - nl == 1) {
            t->val = x;
            t->lz = id();
            t->has_lz = false;
            return;
        }
        push(t, nl, nr);
        long long mid = (nl + nr) >> 1;
        if (p < mid) set(t->l, nl, mid, p, x);
        else set(t->r, mid, nr, p, x);
        pull(t, nl, nr);
    }

    S prod(Node* t, long long nl, long long nr, long long ql, long long qr) {
        if (qr <= nl || nr <= ql) return e();
        if (ql <= nl && nr <= qr) return node_val(t, nl, nr);
        if (!t) {
            long long mid = (nl + nr) >> 1;
            return op(prod(nullptr, nl, mid, ql, qr),
                      prod(nullptr, mid, nr, ql, qr));
        }
        push(t, nl, nr);
        long long mid = (nl + nr) >> 1;
        return op(prod(t->l, nl, mid, ql, qr),
                  prod(t->r, mid, nr, ql, qr));
    }

    void apply(Node*& t, long long nl, long long nr,
               long long ql, long long qr, F f) {
        if (qr <= nl || nr <= ql || seg_len(nl, nr) == 0) return;
        if (!t) t = make_node(nl, nr);
        if (ql <= nl && nr <= qr) {
            all_apply(t, f, nl, nr);
            return;
        }
        push(t, nl, nr);
        long long mid = (nl + nr) >> 1;
        apply(t->l, nl, mid, ql, qr, f);
        apply(t->r, mid, nr, ql, qr, f);
        pull(t, nl, nr);
    }

    template <class G>
    long long max_right(Node* t, long long nl, long long nr,
                        long long ql, G& g, S& sm) {
        if (nr <= ql || seg_len(nl, nr) == 0) return ql;
        S cur = node_val(t, nl, nr);
        if (ql <= nl) {
            S nxt = op(sm, cur);
            if (g(nxt)) {
                sm = nxt;
                return min(nr, _n);
            }
            if (nr - nl == 1) return nl;
        }
        if (nr - nl == 1) return min(nr, _n);
        if (t) push(t, nl, nr);
        long long mid = (nl + nr) >> 1;
        if (ql < mid) {
            long long res = max_right(t ? t->l : nullptr, nl, mid, ql, g, sm);
            if (res < min(mid, _n)) return res;
        }
        return max_right(t ? t->r : nullptr, mid, nr, ql, g, sm);
    }

    template <class G>
    long long min_left(Node* t, long long nl, long long nr,
                       long long qr, G& g, S& sm) {
        if (qr <= nl || seg_len(nl, nr) == 0) return qr;
        S cur = node_val(t, nl, nr);
        if (nr <= qr) {
            S nxt = op(cur, sm);
            if (g(nxt)) {
                sm = nxt;
                return nl;
            }
            if (nr - nl == 1) return nr;
        }
        if (nr - nl == 1) return nl;
        if (t) push(t, nl, nr);
        long long mid = (nl + nr) >> 1;
        if (mid < qr) {
            long long res = min_left(t ? t->r : nullptr, mid, nr, qr, g, sm);
            if (mid < res) return res;
        }
        return min_left(t ? t->l : nullptr, nl, mid, qr, g, sm);
    }

  public:
    dynamic_lazy_segtree() : _n(0), size(1), root(nullptr) {}

    /*
      n 要素で初期化
      - 初期値は init(1) を各点に入れたものとみなされる
      - たとえば sum 系なら 0、min 系なら INF、max 系なら -INF
    */
    explicit dynamic_lazy_segtree(long long n) : _n(n), size(1), root(nullptr) {
        assert(0 <= n);
        while (size < _n) size <<= 1;
    }

    /*
      vector で初期化
      - v[i] をそのまま i 番目に入れる
      - O(n log n)
    */
    explicit dynamic_lazy_segtree(const vector<S>& v)
        : dynamic_lazy_segtree((long long)v.size()) {
        for (long long i = 0; i < (long long)v.size(); ++i) set(i, v[i]);
    }

    long long size_() const { return _n; }

    /* a[p] = x */
    void set(long long p, S x) {
        assert(0 <= p && p < _n);
        set(root, 0, size, p, x);
    }

    /* a[p] を取得 */
    S get(long long p) {
        assert(0 <= p && p < _n);
        return prod(p, p + 1);
    }

    /* op(a[l], ..., a[r-1]) */
    S prod(long long l, long long r) {
        assert(0 <= l && l <= r && r <= _n);
        return prod(root, 0, size, l, r);
    }

    /* op(a[0], ..., a[n-1]) */
    S all_prod() {
        return node_val(root, 0, size);
    }

    /* a[p] に 1 点作用 */
    void apply(long long p, F f) {
        assert(0 <= p && p < _n);
        apply(p, p + 1, f);
    }

    /* a[l..r) に作用 */
    void apply(long long l, long long r, F f) {
        assert(0 <= l && l <= r && r <= _n);
        apply(root, 0, size, l, r, f);
    }

    template <class G>
    long long max_right(long long l, G g) {
        assert(0 <= l && l <= _n);
        assert(g(e()));
        if (l == _n) return _n;
        S sm = e();
        return max_right(root, 0, size, l, g, sm);
    }

    template <class G>
    long long min_left(long long r, G g) {
        assert(0 <= r && r <= _n);
        assert(g(e()));
        if (r == 0) return 0;
        S sm = e();
        return min_left(root, 0, size, r, g, sm);
    }
};

namespace dynamic_lazy_segtree_presets {

using ll = long long;
constexpr ll INF64 = (1LL << 60);

/*
  ここで載せているのは
  - ACL の lazy_segtree 的に素直に書ける代表的なもの
  - 競プロでよく使うもの
  です。

  注意:
  - range chmin / range chmax / range chmin+chmax+sum などは
    通常の lazy segtree ではなく Segment Tree Beats 側です。
  - それらはここには含めていません。
*/

// ============================================================
// 1. range assign / range sum
// ============================================================

namespace range_assign_range_sum {

/*
  作用:
    区間 [l, r) をすべて x にする

  取得:
    区間和

  配列イメージ:
    a[i] : long long

  初期値:
    すべて 0

  典型操作:
    seg.apply(l, r, assign(x));   // a[l..r) = x
    seg.prod(l, r);               // sum
    seg.get(p);                   // a[p]
    seg.set(p, x);                // a[p] = x

  初期化:
    dynamic_lazy_segtree_range_assign_range_sum seg(n);
    vector<long long> v = {...};
    dynamic_lazy_segtree_range_assign_range_sum seg(v);
*/

using S = ll;
struct F {
    ll x;
    bool has;
};

inline S op(S a, S b) { return a + b; }
inline S e() { return 0; }
inline S init(long long) { return 0; }

inline S mapping(F f, S s, long long len) {
    return f.has ? f.x * (ll)len : s;
}

inline F composition(F f, F g) {
    return f.has ? f : g;
}

inline F id() { return {0, false}; }

/* 区間代入作用を作る */
inline F assign(ll x) { return {x, true}; }

using segtree = dynamic_lazy_segtree<S, op, e, init, F, mapping, composition, id>;

}  // namespace range_assign_range_sum

using dynamic_lazy_segtree_range_assign_range_sum =
    range_assign_range_sum::segtree;

// ============================================================
// 2. range add / range sum
// ============================================================

namespace range_add_range_sum {

/*
  作用:
    区間 [l, r) に x を足す

  取得:
    区間和

  配列イメージ:
    a[i] : long long

  初期値:
    すべて 0

  典型操作:
    seg.apply(l, r, add(x));      // a[i] += x
    seg.prod(l, r);               // sum
    seg.get(p);                   // a[p]
    seg.set(p, x);                // a[p] = x

  初期化:
    dynamic_lazy_segtree_range_add_range_sum seg(n);
    vector<long long> v = {...};
    dynamic_lazy_segtree_range_add_range_sum seg(v);
*/

using S = ll;
struct F {
    ll x;
};

inline S op(S a, S b) { return a + b; }
inline S e() { return 0; }
inline S init(long long) { return 0; }

inline S mapping(F f, S s, long long len) {
    return s + f.x * (ll)len;
}

inline F composition(F f, F g) {
    return {f.x + g.x};
}

inline F id() { return {0}; }

/* 区間加算作用を作る */
inline F add(ll x) { return {x}; }

using segtree = dynamic_lazy_segtree<S, op, e, init, F, mapping, composition, id>;

}  // namespace range_add_range_sum

using dynamic_lazy_segtree_range_add_range_sum =
    range_add_range_sum::segtree;

// ============================================================
// 3. range assign / range min
// ============================================================

namespace range_assign_range_min {

/*
  作用:
    区間 [l, r) をすべて x にする

  取得:
    区間最小値

  配列イメージ:
    a[i] : long long

  初期値:
    すべて INF64

  注意:
    「未設定は十分大きい値」とみなす用途向け。
    0 初期化の min をやりたいなら vector で 0 を詰めて初期化するか、
    全域に assign(0) してください。

  典型操作:
    seg.apply(l, r, assign(x));
    seg.prod(l, r);               // min
    seg.get(p);
    seg.set(p, x);

  初期化:
    dynamic_lazy_segtree_range_assign_range_min seg(n);     // 全部 INF64
    vector<long long> v(n, 0);
    dynamic_lazy_segtree_range_assign_range_min seg(v);     // 全部 0
*/

using S = ll;
struct F {
    ll x;
    bool has;
};

inline S op(S a, S b) { return min(a, b); }
inline S e() { return INF64; }
inline S init(long long) { return INF64; }

inline S mapping(F f, S s, long long) {
    return f.has ? f.x : s;
}

inline F composition(F f, F g) {
    return f.has ? f : g;
}

inline F id() { return {0, false}; }

inline F assign(ll x) { return {x, true}; }

using segtree = dynamic_lazy_segtree<S, op, e, init, F, mapping, composition, id>;

}  // namespace range_assign_range_min

using dynamic_lazy_segtree_range_assign_range_min =
    range_assign_range_min::segtree;

// ============================================================
// 4. range add / range min
// ============================================================

namespace range_add_range_min {

/*
  作用:
    区間 [l, r) に x を足す

  取得:
    区間最小値

  配列イメージ:
    a[i] : long long

  初期値:
    すべて INF64

  注意:
    これも未設定は INF64 とみなす。
    0 初期化したいなら vector<long long>(n, 0) で初期化。

  典型操作:
    seg.apply(l, r, add(x));
    seg.prod(l, r);               // min
    seg.get(p);
    seg.set(p, x);
*/

using S = ll;
struct F {
    ll x;
};

inline S op(S a, S b) { return min(a, b); }
inline S e() { return INF64; }
inline S init(long long) { return INF64; }

inline S mapping(F f, S s, long long) {
    return s + f.x;
}

inline F composition(F f, F g) {
    return {f.x + g.x};
}

inline F id() { return {0}; }

inline F add(ll x) { return {x}; }

using segtree = dynamic_lazy_segtree<S, op, e, init, F, mapping, composition, id>;

}  // namespace range_add_range_min

using dynamic_lazy_segtree_range_add_range_min =
    range_add_range_min::segtree;

// ============================================================
// 5. range assign / range max
// ============================================================

namespace range_assign_range_max {

/*
  作用:
    区間 [l, r) をすべて x にする

  取得:
    区間最大値

  配列イメージ:
    a[i] : long long

  初期値:
    すべて -INF64

  注意:
    「未設定は十分小さい値」とみなす用途向け。
    0 初期化したいなら vector で 0 を入れて初期化。
*/

using S = ll;
struct F {
    ll x;
    bool has;
};

inline S op(S a, S b) { return max(a, b); }
inline S e() { return -INF64; }
inline S init(long long) { return -INF64; }

inline S mapping(F f, S s, long long) {
    return f.has ? f.x : s;
}

inline F composition(F f, F g) {
    return f.has ? f : g;
}

inline F id() { return {0, false}; }

inline F assign(ll x) { return {x, true}; }

using segtree = dynamic_lazy_segtree<S, op, e, init, F, mapping, composition, id>;

}  // namespace range_assign_range_max

using dynamic_lazy_segtree_range_assign_range_max =
    range_assign_range_max::segtree;

// ============================================================
// 6. range add / range max
// ============================================================

namespace range_add_range_max {

/*
  作用:
    区間 [l, r) に x を足す

  取得:
    区間最大値

  配列イメージ:
    a[i] : long long

  初期値:
    すべて -INF64

  注意:
    0 初期化したいなら vector<long long>(n, 0) で初期化。
*/

using S = ll;
struct F {
    ll x;
};

inline S op(S a, S b) { return max(a, b); }
inline S e() { return -INF64; }
inline S init(long long) { return -INF64; }

inline S mapping(F f, S s, long long) {
    return s + f.x;
}

inline F composition(F f, F g) {
    return {f.x + g.x};
}

inline F id() { return {0}; }

inline F add(ll x) { return {x}; }

using segtree = dynamic_lazy_segtree<S, op, e, init, F, mapping, composition, id>;

}  // namespace range_add_range_max

using dynamic_lazy_segtree_range_add_range_max =
    range_add_range_max::segtree;

// ============================================================
// 7. range assign / range minmax
// ============================================================

namespace range_assign_range_minmax {

/*
  作用:
    区間 [l, r) をすべて x にする

  取得:
    区間の (min, max)

  配列イメージ:
    a[i] : long long

  初期値:
    未設定は (INF64, -INF64)

  用途:
    min と max を同時に見たいとき

  典型操作:
    auto res = seg.prod(l, r);
    cout << res.mn << " " << res.mx << "\n";
*/

struct S {
    ll mn, mx;
};
struct F {
    ll x;
    bool has;
};

inline S op(S a, S b) {
    return {min(a.mn, b.mn), max(a.mx, b.mx)};
}

inline S e() {
    return {INF64, -INF64};
}

inline S init(long long) {
    return e();
}

inline S mapping(F f, S s, long long) {
    return f.has ? S{f.x, f.x} : s;
}

inline F composition(F f, F g) {
    return f.has ? f : g;
}

inline F id() { return {0, false}; }

inline F assign(ll x) { return {x, true}; }

using segtree = dynamic_lazy_segtree<S, op, e, init, F, mapping, composition, id>;

}  // namespace range_assign_range_minmax

using dynamic_lazy_segtree_range_assign_range_minmax =
    range_assign_range_minmax::segtree;

// ============================================================
// 8. range add / range minmax
// ============================================================

namespace range_add_range_minmax {

/*
  作用:
    区間 [l, r) に x を足す

  取得:
    区間の (min, max)

  配列イメージ:
    a[i] : long long

  初期値:
    未設定は (INF64, -INF64)

  典型操作:
    seg.apply(l, r, add(x));
    auto res = seg.prod(l, r);
*/

struct S {
    ll mn, mx;
};
struct F {
    ll x;
};

inline S op(S a, S b) {
    return {min(a.mn, b.mn), max(a.mx, b.mx)};
}

inline S e() {
    return {INF64, -INF64};
}

inline S init(long long) {
    return e();
}

inline S mapping(F f, S s, long long) {
    if (s.mn == INF64 && s.mx == -INF64) return s;
    return {s.mn + f.x, s.mx + f.x};
}

inline F composition(F f, F g) {
    return {f.x + g.x};
}

inline F id() { return {0}; }

inline F add(ll x) { return {x}; }

using segtree = dynamic_lazy_segtree<S, op, e, init, F, mapping, composition, id>;

}  // namespace range_add_range_minmax

using dynamic_lazy_segtree_range_add_range_minmax =
    range_add_range_minmax::segtree;

// ============================================================
// 9. range affine / range sum
// ============================================================

template <class mint>
struct range_affine_range_sum {

/*
  作用:
    区間 [l, r) の各要素 x に対して x -> a*x + b

  取得:
    区間和

  配列イメージ:
    a[i] : mint

  初期値:
    すべて 0

  典型用途:
    modint 上の affine update + range sum

  初期化:
    using mint = atcoder::modint998244353;
    dynamic_lazy_segtree_range_affine_range_sum<mint> seg(n);

    vector<range_affine_range_sum<mint>::S> v(n);
    for (int i = 0; i < n; i++) v[i] = {mint(x[i]), 1};
    dynamic_lazy_segtree_range_affine_range_sum<mint> seg(v);

  典型操作:
    seg.apply(l, r, {a, b});      // x -> a*x + b
    auto res = seg.prod(l, r);
    cout << res.sum.val() << '\n';

  注意:
    S は {sum, size}
    vector で初期化するときは各点の size=1 にする
*/

    struct S {
        mint sum;
        long long size;
    };
    struct F {
        mint a, b;
    };

    static S op(S l, S r) {
        return {l.sum + r.sum, l.size + r.size};
    }

    static S e() {
        return {0, 0};
    }

    static S init(long long len) {
        return {0, len};
    }

    static S mapping(F f, S s, long long) {
        return {f.a * s.sum + f.b * s.size, s.size};
    }

    static F composition(F f, F g) {
        return {f.a * g.a, f.a * g.b + f.b};
    }

    static F id() {
        return {1, 0};
    }

    using segtree = dynamic_lazy_segtree<S, op, e, init, F, mapping, composition, id>;
};

template <class mint>
using dynamic_lazy_segtree_range_affine_range_sum =
    typename range_affine_range_sum<mint>::segtree;

// ============================================================
// 10. range flip / range inversion
// ============================================================

namespace range_flip_range_inversion {

/*
  作用:
    0/1 配列の区間反転
    0 <-> 1

  取得:
    区間転倒数

  配列イメージ:
    a[i] は 0 または 1

  S:
    zero      = 0 の個数
    one       = 1 の個数
    inversion = 転倒数

  初期値:
    すべて 0
    よって各点は {1, 0, 0}

  初期化:
    dynamic_lazy_segtree_range_flip_range_inversion seg(n);   // 全部 0

    vector<S> v(n);
    for (int i = 0; i < n; i++) {
        if (bit[i] == 0) v[i] = {1, 0, 0};
        else             v[i] = {0, 1, 0};
    }
    dynamic_lazy_segtree_range_flip_range_inversion seg(v);

  典型操作:
    seg.apply(l, r, true);        // 0/1 を反転
    auto res = seg.prod(l, r);
    cout << res.inversion << '\n';

  get(p) で返るのも S なので、
    zero==1 なら 0
    one ==1 なら 1
  と判定する
*/

struct S {
    long long zero, one, inversion;
};
using F = bool;

inline S op(S l, S r) {
    return {
        l.zero + r.zero,
        l.one + r.one,
        l.inversion + r.inversion + l.one * r.zero
    };
}

inline S e() {
    return {0, 0, 0};
}

inline S init(long long len) {
    return {len, 0, 0};
}

inline S mapping(F f, S s, long long) {
    if (!f) return s;
    return {s.one, s.zero, s.one * s.zero - s.inversion};
}

inline F composition(F f, F g) {
    return f ^ g;
}

inline F id() {
    return false;
}

using segtree = dynamic_lazy_segtree<S, op, e, init, F, mapping, composition, id>;

}  // namespace range_flip_range_inversion

using dynamic_lazy_segtree_range_flip_range_inversion =
    range_flip_range_inversion::segtree;

}  // namespace dynamic_lazy_segtree_presets

/*
------------------------------------------------------------
簡単な使用例
------------------------------------------------------------
*/
/*
int main() {
    using namespace dynamic_lazy_segtree_presets;

    // --------------------------------------------------------
    // 1) range add / range sum
    //   初期値 0 の長さ 1e9 の配列を扱う
    // --------------------------------------------------------
    dynamic_lazy_segtree_range_add_range_sum seg1((long long)1e9);

    seg1.apply(3, 10, range_add_range_sum::add(5));   // [3,10) += 5
    seg1.apply(7, 20, range_add_range_sum::add(2));   // [7,20) += 2

    cout << seg1.prod(0, 100) << '\n';                // 区間和
    cout << seg1.get(8) << '\n';                      // 1 点取得
    seg1.set(8, 100);                                 // a[8] = 100
    cout << seg1.get(8) << '\n';

    // --------------------------------------------------------
    // 2) range assign / range min
    //   vector で 0 初期化してから使う例
    // --------------------------------------------------------
    vector<long long> init_min(10, 0);
    dynamic_lazy_segtree_range_assign_range_min seg2(init_min);

    seg2.apply(2, 7, range_assign_range_min::assign(5));
    cout << seg2.prod(0, 10) << '\n';                 // min

    // --------------------------------------------------------
    // 3) range flip / range inversion
    //   0/1 列 [0,1,1,0] を作る例
    // --------------------------------------------------------
    {
        using T = range_flip_range_inversion::S;
        vector<T> v = {
            {1, 0, 0},   // 0
            {0, 1, 0},   // 1
            {0, 1, 0},   // 1
            {1, 0, 0},   // 0
        };
        dynamic_lazy_segtree_range_flip_range_inversion seg3(v);
        cout << seg3.prod(0, 4).inversion << '\n';
        seg3.apply(1, 4, true);                       // flip
        cout << seg3.prod(0, 4).inversion << '\n';
    }

    return 0;
}
*/
