#include <bits/stdc++.h>
using namespace std;

/*
    ============================================================
    Dynamic Persistent Dual Segment Tree
    ============================================================

    この実装は「区間更新・一点取得」専用の双対セグ木です。
    さらに path copying により永続化しています。

    重要:
      - 双対セグ木なので、得意なのは
            range_apply(l, r, f) : [l, r) に作用 f を適用
            point_get(p)         : 点 p の値を取得
        です。
      - 通常の意味での RMQ (range minimum query) には向きません。
        RMQ をしたい場合は通常の segtree / lazy segtree を使ってください。

    Action に要求するもの:
      using S = 点の値の型
      using F = 作用の型

      static F id()
        - 作用の単位元

      static F compose(const F& f, const F& g)
        - 「先に f を適用し、その後 g を適用する」合成
        - つまり compose(f, g) = g ∘ f

      static S apply(const F& f, const S& x)
        - 値 x に作用 f を適用した結果

      static S value_id()
        - point_get(version, p) で基底値省略時に使う初期値

    例:
      加算:
        compose(f, g) = f + g
        apply(f, x)   = x + f

      代入:
        compose(f, g) = (g が有効なら g, そうでなければ f)
        apply(f, x)   = (f が有効なら f.value, そうでなければ x)
*/

/*
    ------------------------------------------------------------
    作用例 1: Range Add / Point Get
    ------------------------------------------------------------
    区間に +delta を加える。
*/
template <class T>
struct RangeAddPointGet {
    using S = T;
    using F = T;

    // 単位作用
    static F id() {
        return T(0);
    }

    // 先に f, 後に g
    static F compose(const F& f, const F& g) {
        return f + g;
    }

    // 値 x に作用 f を適用
    static S apply(const F& f, const S& x) {
        return x + f;
    }

    // 基底値省略時の初期値
    static S value_id() {
        return T(0);
    }

    // 使用例:
    //   auto f = RangeAddPointGet<long long>::make(5);
    static F make(const T& delta) {
        return delta;
    }
};

/*
    ------------------------------------------------------------
    作用例 2: Range Update(Assign) / Point Get
    ------------------------------------------------------------
    区間を値 v で上書きする。
*/
template <class T>
struct RangeUpdatePointGet {
    using S = T;

    struct F {
        T value;
        bool has;   // true なら「代入作用あり」
    };

    static F id() {
        return F{T(), false};
    }

    // 先に f, 後に g
    // 後から来た代入が優先される
    static F compose(const F& f, const F& g) {
        return g.has ? g : f;
    }

    static S apply(const F& f, const S& x) {
        return f.has ? f.value : x;
    }

    static S value_id() {
        return T(0);
    }

    // 使用例:
    //   auto f = RangeUpdatePointGet<long long>::make(7);
    static F make(const T& v) {
        return F{v, true};
    }
};

/*
    ------------------------------------------------------------
    作用例 3: Range XOR / Point Get
    ------------------------------------------------------------
    区間に xor mask を適用する。
*/
template <class T>
struct RangeXorPointGet {
    using S = T;
    using F = T;

    static F id() {
        return T(0);
    }

    static F compose(const F& f, const F& g) {
        return f ^ g;
    }

    static S apply(const F& f, const S& x) {
        return x ^ f;
    }

    static S value_id() {
        return T(0);
    }

    // 使用例:
    //   auto f = RangeXorPointGet<int>::make(3);
    static F make(const T& mask) {
        return mask;
    }
};

/*
    ------------------------------------------------------------
    作用例 4: Range Affine / Point Get
    ------------------------------------------------------------
    各点に対して x -> a*x + b を適用する。
    先に f(x)=f.a*x+f.b, 後に g(x)=g.a*x+g.b とすると、
      compose(f, g) = g ∘ f
                    = (g.a*f.a) x + (g.a*f.b + g.b)
*/
template <class T>
struct RangeAffinePointGet {
    using S = T;

    struct F {
        T a, b;
    };

    static F id() {
        return F{T(1), T(0)};
    }

    static F compose(const F& f, const F& g) {
        return F{
            g.a * f.a,
            g.a * f.b + g.b
        };
    }

    static S apply(const F& f, const S& x) {
        return f.a * x + f.b;
    }

    static S value_id() {
        return T(0);
    }

    // 使用例:
    //   auto f = RangeAffinePointGet<long long>::make(2, 3); // x -> 2x+3
    static F make(const T& a, const T& b) {
        return F{a, b};
    }
};

/*
    ============================================================
    完全永続な核
    ============================================================

    内部的には path copying により永続化しています。
    API としては任意の root から更新できるので、理論的には
    「完全永続」に近い実装です。

    この後に示す PartiallyPersistentDynamicDualSegTree が
    「最新 version からしか更新しない」制約を付けた
    部分永続ラッパです。
*/
template <class Action>
class DynamicPersistentDualSegTreeCore {
public:
    using S = typename Action::S;
    using F = typename Action::F;

private:
    struct Node {
        F lazy;
        Node* left;
        Node* right;

        Node(const F& lazy_ = Action::id(), Node* left_ = nullptr, Node* right_ = nullptr)
            : lazy(lazy_), left(left_), right(right_) {}
    };

public:
    using root_type = Node*;

private:
    long long L0, R0;   // 管理範囲 [L0, R0)

    /*
        clone(t)
        - 使い方:
            ノード t を 1 個複製する。
            永続化のため、更新時には元ノードを壊さず複製して使う。
        - 時間計算量:
            O(1)
    */
    static Node* clone(Node* t) {
        if (!t) return new Node();
        return new Node(t->lazy, t->left, t->right);
    }

    /*
        range_apply_impl(t, nl, nr, ql, qr, f)
        - 使い方:
            現在ノード t が表す区間 [nl, nr) に対して、
            [ql, qr) に作用 f を適用した新しい root を返す。
            永続化のため、必要なノードだけ複製する。
        - 時間計算量:
            O(log(R0-L0))
          生成ノード数も O(log(R0-L0))
    */
    Node* range_apply_impl(Node* t, long long nl, long long nr,
                           long long ql, long long qr, const F& f) {
        if (qr <= nl || nr <= ql) {
            return t;
        }

        Node* res = clone(t);

        if (ql <= nl && nr <= qr) {
            res->lazy = Action::compose(res->lazy, f);
            return res;
        }

        long long mid = nl + ((nr - nl) >> 1);
        if (ql < mid) {
            res->left = range_apply_impl(res->left, nl, mid, ql, qr, f);
        }
        if (mid < qr) {
            res->right = range_apply_impl(res->right, mid, nr, ql, qr, f);
        }
        return res;
    }

    /*
        point_get_op_impl(t, nl, nr, p, acc)
        - 使い方:
            点 p に効いている作用を root から leaf までたどって集約する。
            acc は今までに集約した作用。
        - 時間計算量:
            O(log(R0-L0))
    */
    F point_get_op_impl(Node* t, long long nl, long long nr, long long p, F acc) const {
        if (t) {
            acc = Action::compose(acc, t->lazy);
        }

        if (nr - nl == 1) {
            return acc;
        }

        long long mid = nl + ((nr - nl) >> 1);
        if (p < mid) {
            return point_get_op_impl(t ? t->left : nullptr, nl, mid, p, acc);
        } else {
            return point_get_op_impl(t ? t->right : nullptr, mid, nr, p, acc);
        }
    }

public:
    /*
        DynamicPersistentDualSegTreeCore(L, R)
        - 使い方:
            管理範囲 [L, R) を持つ空の dynamic persistent dual segtree を作る。
            座標範囲が大きくても、アクセスした部分だけノード生成する。
        - 時間計算量:
            O(1)
    */
    DynamicPersistentDualSegTreeCore(long long L, long long R) : L0(L), R0(R) {
        assert(L0 < R0);
    }

    /*
        empty_root()
        - 使い方:
            何も更新していない空 version の root を返す。
            典型的には最初の version に使う。
        - 時間計算量:
            O(1)

        使用例:
            auto root0 = seg.empty_root();
    */
    root_type empty_root() const {
        return nullptr;
    }

    /*
        range_apply(root, l, r, f)
        - 使い方:
            version root に対し、区間 [l, r) に作用 f を適用した
            新しい version の root を返す。
            元の root は変更されない。
        - 時間計算量:
            O(log(R0-L0))
          追加メモリ:
            O(log(R0-L0))

        使用例:
            auto root1 = seg.range_apply(root0, 10, 20, RangeAddPointGet<ll>::make(5));
    */
    root_type range_apply(root_type root, long long l, long long r, const F& f) {
        assert(L0 <= l && l <= r && r <= R0);
        return range_apply_impl(root, L0, R0, l, r, f);
    }

    /*
        point_get_op(root, p)
        - 使い方:
            version root において、点 p にかかっている「作用そのもの」を返す。
            値ではなく作用を見たいときに使う。
        - 時間計算量:
            O(log(R0-L0))

        使用例:
            auto op = seg.point_get_op(root1, 15);
    */
    F point_get_op(root_type root, long long p) const {
        assert(L0 <= p && p < R0);
        return point_get_op_impl(root, L0, R0, p, Action::id());
    }

    /*
        point_get(root, p, base_value)
        - 使い方:
            version root における点 p の値を返す。
            初期値が一様に base_value であるとみなして、
            そこへ accumulated action を適用する。
        - 時間計算量:
            O(log(R0-L0))

        使用例:
            long long x = seg.point_get(root1, 15, 0LL);
    */
    S point_get(root_type root, long long p, const S& base_value) const {
        return Action::apply(point_get_op(root, p), base_value);
    }

    /*
        point_get(root, p)
        - 使い方:
            初期値を Action::value_id() として点 p の値を返す。
        - 時間計算量:
            O(log(R0-L0))

        使用例:
            long long x = seg.point_get(root1, 15);
    */
    S point_get(root_type root, long long p) const {
        return point_get(root, p, Action::value_id());
    }
};

/*
    ============================================================
    部分永続ラッパ
    ============================================================

    「更新は常に最新 version に対してだけ許す」ラッパです。
    これにより運用上は部分永続になります。

    部分永続の意味:
      - 過去 version は参照可能
      - 更新は最新 version からのみ行う

    内部実装は path copying なので、更新のたびに
    O(log(R-L)) 個のノードだけ新規作成し、
    触っていない部分木は古い version と共有します。
*/
template <class Action>
class PartiallyPersistentDynamicDualSegTree {
public:
    using core_type = DynamicPersistentDualSegTreeCore<Action>;
    using S = typename Action::S;
    using F = typename Action::F;
    using root_type = typename core_type::root_type;

private:
    core_type core;
    vector<root_type> roots;   // roots[i] = version i の root

public:
    /*
        PartiallyPersistentDynamicDualSegTree(L, R)
        - 使い方:
            version 0 を空状態として作成する。
        - 時間計算量:
            O(1)

        使用例:
            PartiallyPersistentDynamicDualSegTree<RangeAddPointGet<long long>> seg(0, 1LL<<30);
    */
    PartiallyPersistentDynamicDualSegTree(long long L, long long R) : core(L, R) {
        roots.push_back(core.empty_root()); // version 0
    }

    /*
        num_versions()
        - 使い方:
            現在保持している version 数を返す。
        - 時間計算量:
            O(1)

        使用例:
            int n = seg.num_versions();
    */
    int num_versions() const {
        return (int)roots.size();
    }

    /*
        latest_version()
        - 使い方:
            最新 version 番号を返す。
        - 時間計算量:
            O(1)

        使用例:
            int v = seg.latest_version();
    */
    int latest_version() const {
        return (int)roots.size() - 1;
    }

    /*
        range_apply(l, r, f)
        - 使い方:
            最新 version に対して [l, r) に作用 f を適用し、
            新しい version を 1 つ追加する。
        - 時間計算量:
            O(log(R-L))
          追加メモリ:
            O(log(R-L))

        使用例:
            seg.range_apply(10, 20, RangeAddPointGet<long long>::make(5));
            // これで新しい version が末尾に追加される
    */
    void range_apply(long long l, long long r, const F& f) {
        roots.push_back(core.range_apply(roots.back(), l, r, f));
    }

    /*
        point_get_op(version, p)
        - 使い方:
            指定 version において点 p にかかる作用を返す。
        - 時間計算量:
            O(log(R-L))

        使用例:
            auto op = seg.point_get_op(3, 15);
    */
    F point_get_op(int version, long long p) const {
        assert(0 <= version && version < (int)roots.size());
        return core.point_get_op(roots[version], p);
    }

    /*
        point_get(version, p, base_value)
        - 使い方:
            指定 version における点 p の値を返す。
            初期値が一様に base_value であるとみなす。
        - 時間計算量:
            O(log(R-L))

        使用例:
            long long x = seg.point_get(2, 15, 0LL);
    */
    S point_get(int version, long long p, const S& base_value) const {
        assert(0 <= version && version < (int)roots.size());
        return core.point_get(roots[version], p, base_value);
    }

    /*
        point_get(version, p)
        - 使い方:
            初期値を Action::value_id() として点 p の値を返す。
        - 時間計算量:
            O(log(R-L))

        使用例:
            long long x = seg.point_get(2, 15);
    */
    S point_get(int version, long long p) const {
        assert(0 <= version && version < (int)roots.size());
        return core.point_get(roots[version], p);
    }

    /*
        root(version)
        - 使い方:
            指定 version の root を取得する。
            デバッグや独自操作用。
        - 時間計算量:
            O(1)
    */
    root_type root(int version) const {
        assert(0 <= version && version < (int)roots.size());
        return roots[version];
    }
};

/*
    ============================================================
    よく使う using
    ============================================================

    双対セグ木なので、よく使うのは
      - RAQ : Range Add Query(ここでは point get)
      - RUQ : Range Update(Assign) Query(ここでは point get)
    です。

    注意:
      - 通常の RMQ は range minimum query の意味であり、
        双対セグ木の守備範囲ではありません。
*/
template <class T>
using DynamicPPDualSegTree_RAQ = PartiallyPersistentDynamicDualSegTree<RangeAddPointGet<T>>;

template <class T>
using DynamicPPDualSegTree_RUQ = PartiallyPersistentDynamicDualSegTree<RangeUpdatePointGet<T>>;

template <class T>
using DynamicPPDualSegTree_RXQ = PartiallyPersistentDynamicDualSegTree<RangeXorPointGet<T>>;

template <class T>
using DynamicPPDualSegTree_RAffineQ = PartiallyPersistentDynamicDualSegTree<RangeAffinePointGet<T>>;

/*
    ============================================================
    使用例
    ============================================================
*/
/*
int main() {
    using ll = long long;

    {
        /*
            RAQ:
              区間加算・一点取得

            version 0: 何もしていない
            version 1: [10, 20) に +5
            version 2: [15, 30) に +7

            点 16 の値:
              v0 -> 0
              v1 -> 5
              v2 -> 12
        */
        DynamicPPDualSegTree_RAQ<ll> seg(0, 1LL << 30);

        seg.range_apply(10, 20, RangeAddPointGet<ll>::make(5)); // version 1
        seg.range_apply(15, 30, RangeAddPointGet<ll>::make(7)); // version 2

        cout << seg.point_get(0, 16) << '\n'; // 0
        cout << seg.point_get(1, 16) << '\n'; // 5
        cout << seg.point_get(2, 16) << '\n'; // 12
    }

    {
        /*
            RUQ:
              区間代入・一点取得

            version 0: 全部 0
            version 1: [3, 8) を 10 にする
            version 2: [5, 7) を 4 にする
        */
        DynamicPPDualSegTree_RUQ<ll> seg(0, 16);

        seg.range_apply(3, 8, RangeUpdatePointGet<ll>::make(10)); // version 1
        seg.range_apply(5, 7, RangeUpdatePointGet<ll>::make(4));  // version 2

        cout << seg.point_get(0, 6) << '\n'; // 0
        cout << seg.point_get(1, 6) << '\n'; // 10
        cout << seg.point_get(2, 6) << '\n'; // 4
        cout << seg.point_get(2, 2) << '\n'; // 0
    }

    return 0;
}
*/
