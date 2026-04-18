#include <bits/stdc++.h>
using namespace std;

/*
    Generic Old Driver Tree / Chtholly Tree

    - 1-indexed closed intervals [l, r]
    - stores maximal constant segments
    - supports:
        * assign(l, r, x)
        * update(l, r, f)     // apply f(value) for each segment in [l, r]
        * fold(l, r)          // monoid aggregate
        * kth(l, r, k)        // k-th smallest value in [l, r] (1-indexed)
        * dump()              // all segments

    Requirements:
      S must be comparable by operator< and operator==
      Monoid must provide:
        using value_type;
        static value_type id();
        static value_type op(const value_type&, const value_type&);
        static value_type lift(const S&, long long len);
*/

template <class S, class Monoid>
class OldDriverTree {
public:
    using T = typename Monoid::value_type;

private:
    struct Node {
        int l, r;
        mutable S v;

        Node(int l_, int r_, const S& v_) : l(l_), r(r_), v(v_) {}

        bool operator<(const Node& other) const {
            return l < other.l;
        }
    };

    int n_;
    std::set<Node> st_;

    typename std::set<Node>::iterator split(int pos) {
        if (pos > n_) return st_.end();
        auto it = prev(st_.upper_bound(Node(pos, 0, S{})));
        if (it->l == pos) return it;
        if (it->r < pos) return next(it);

        int l = it->l, r = it->r;
        S v = it->v;
        st_.erase(it);
        st_.insert(Node(l, pos - 1, v));
        return st_.insert(Node(pos, r, v)).first;
    }

    void merge_around(typename std::set<Node>::iterator it) {
        if (it == st_.end()) return;

        // merge with previous
        while (it != st_.begin()) {
            auto prv = prev(it);
            if (prv->r + 1 == it->l && prv->v == it->v) {
                int nl = prv->l, nr = it->r;
                S v = it->v;
                st_.erase(prv);
                st_.erase(it);
                it = st_.insert(Node(nl, nr, v)).first;
            } else {
                break;
            }
        }

        // merge with next
        while (true) {
            auto nxt = next(it);
            if (nxt != st_.end() && it->r + 1 == nxt->l && it->v == nxt->v) {
                int nl = it->l, nr = nxt->r;
                S v = it->v;
                st_.erase(nxt);
                st_.erase(it);
                it = st_.insert(Node(nl, nr, v)).first;
            } else {
                break;
            }
        }
    }

public:
    OldDriverTree() : n_(0) {}

    explicit OldDriverTree(const vector<S>& a) {
        build(a);
    }

    void build(const vector<S>& a) {
        // a must be 0-indexed vector of size n
        n_ = (int)a.size();
        st_.clear();
        if (n_ == 0) return;

        int l = 1;
        for (int i = 2; i <= n_; ++i) {
            if (!(a[i - 1] == a[l - 1])) {
                st_.insert(Node(l, i - 1, a[l - 1]));
                l = i;
            }
        }
        st_.insert(Node(l, n_, a[l - 1]));
    }

    int size() const { return n_; }
    bool empty() const { return n_ == 0; }

    S get(int pos) const {
        auto it = prev(st_.upper_bound(Node(pos, 0, S{})));
        return it->v;
    }

    vector<tuple<int,int,S>> dump() const {
        vector<tuple<int,int,S>> res;
        res.reserve(st_.size());
        for (auto&& nd : st_) {
            res.emplace_back(nd.l, nd.r, nd.v);
        }
        return res;
    }

    vector<S> materialize() const {
        vector<S> a(n_);
        for (auto&& nd : st_) {
            for (int i = nd.l; i <= nd.r; ++i) a[i - 1] = nd.v;
        }
        return a;
    }

    void assign(int l, int r, const S& x) {
        if (l > r || l < 1 || r > n_) return;
        auto itr = split(r + 1);
        auto itl = split(l);
        st_.erase(itl, itr);
        auto it = st_.insert(Node(l, r, x)).first;
        merge_around(it);
    }

    template <class F>
    void update(int l, int r, F f) {
        if (l > r || l < 1 || r > n_) return;
        auto itr = split(r + 1);
        auto itl = split(l);

        vector<Node> buf;
        for (auto it = itl; it != itr; ++it) {
            S nv = f(it->v);
            buf.emplace_back(it->l, it->r, nv);
        }
        st_.erase(itl, itr);

        for (auto& nd : buf) {
            auto it = st_.insert(nd).first;
            merge_around(it);
        }
    }

    template <class F>
    void for_each_segment(int l, int r, F f) const {
        if (l > r || l < 1 || r > n_) return;
        auto itr = const_cast<OldDriverTree*>(this)->split(r + 1);
        auto itl = const_cast<OldDriverTree*>(this)->split(l);
        for (auto it = itl; it != itr; ++it) {
            f(it->l, it->r, it->v);
        }
    }

    T fold(int l, int r) const {
        if (l > r || l < 1 || r > n_) return Monoid::id();
        auto itr = const_cast<OldDriverTree*>(this)->split(r + 1);
        auto itl = const_cast<OldDriverTree*>(this)->split(l);

        T res = Monoid::id();
        for (auto it = itl; it != itr; ++it) {
            long long len = 1LL * it->r - it->l + 1;
            res = Monoid::op(res, Monoid::lift(it->v, len));
        }
        return res;
    }

    S kth(int l, int r, long long k) const {
        if (l > r || l < 1 || r > n_) throw std::out_of_range("invalid range");

        auto itr = const_cast<OldDriverTree*>(this)->split(r + 1);
        auto itl = const_cast<OldDriverTree*>(this)->split(l);

        vector<pair<S, long long>> segs;
        for (auto it = itl; it != itr; ++it) {
            segs.emplace_back(it->v, 1LL * it->r - it->l + 1);
        }
        sort(segs.begin(), segs.end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
        });

        for (auto&& [v, cnt] : segs) {
            if (k <= cnt) return v;
            k -= cnt;
        }
        throw std::out_of_range("k is too large");
    }
};

/* =========================
   Common monoids
   ========================= */

// Sum of values
template <class S, class T = S>
struct SumMonoid {
    using value_type = T;
    static value_type id() { return value_type(0); }
    static value_type op(const value_type& a, const value_type& b) { return a + b; }
    static value_type lift(const S& x, long long len) { return value_type(x) * value_type(len); }
};

// Minimum value on range
template <class S>
struct MinMonoid {
    using value_type = S;
    static value_type id() { return std::numeric_limits<value_type>::max(); }
    static value_type op(const value_type& a, const value_type& b) { return std::min(a, b); }
    static value_type lift(const S& x, long long) { return x; }
};

// Maximum value on range
template <class S>
struct MaxMonoid {
    using value_type = S;
    static value_type id() { return std::numeric_limits<value_type>::lowest(); }
    static value_type op(const value_type& a, const value_type& b) { return std::max(a, b); }
    static value_type lift(const S& x, long long) { return x; }
};

// Bitwise OR
template <class S>
struct OrMonoid {
    using value_type = S;
    static value_type id() { return value_type(0); }
    static value_type op(const value_type& a, const value_type& b) { return a | b; }
    static value_type lift(const S& x, long long) { return x; }
};

// Bitwise AND
template <class S>
struct AndMonoid {
    using value_type = S;
    static value_type id() { return ~value_type(0); }
    static value_type op(const value_type& a, const value_type& b) { return a & b; }
    static value_type lift(const S& x, long long) { return x; }
};

// GCD
template <class S>
struct GcdMonoid {
    using value_type = S;
    static value_type id() { return value_type(0); }
    static value_type op(const value_type& a, const value_type& b) { return std::gcd(a, b); }
    static value_type lift(const S& x, long long) { return x; }
};

// XOR
template <class S>
struct XorMonoid {
    using value_type = S;
    static value_type id() { return value_type(0); }
    static value_type op(const value_type& a, const value_type& b) { return a ^ b; }
    static value_type lift(const S& x, long long len) {
        return (len & 1LL) ? x : value_type(0);
    }
};

// Count of elements equal to One
template <class S, S One = S(1), class T = long long>
struct CountEqualMonoid {
    using value_type = T;
    static value_type id() { return value_type(0); }
    static value_type op(const value_type& a, const value_type& b) { return a + b; }
    static value_type lift(const S& x, long long len) {
        return (x == One ? value_type(len) : value_type(0));
    }
};

/* =========================
   Friendly aliases
   ========================= */

using RangeSumLL = SumMonoid<long long>;
using RangeMinLL = MinMonoid<long long>;
using RangeMaxLL = MaxMonoid<long long>;
using RangeGcdLL = GcdMonoid<long long>;
using RangeXorLL = XorMonoid<long long>;
using RangeOrULL = OrMonoid<unsigned long long>;
using RangeAndULL = AndMonoid<unsigned long long>;
using CountOnes = CountEqualMonoid<int, 1, long long>;

/* =========================
   Example
   ========================= */

/*
int main() {
    vector<long long> a = {1, 1, 2, 2, 2, 5, 5, 3, 3, 3};

    OldDriverTree<long long, RangeSumLL> odt(a);

    // [3, 7] を 10 に代入
    odt.assign(3, 7, 10);

    // [1, 5] の各値に +2
    odt.update(1, 5, [](long long x) { return x + 2; });

    // 区間和
    cout << odt.fold(1, 10) << '\n';

    // [1, 10] の 4 番目に小さい値
    cout << odt.kth(1, 10, 4) << '\n';

    // 中身確認
    for (auto [l, r, v] : odt.dump()) {
        cout << "[" << l << "," << r << "] = " << v << '\n';
    }
}
*/

/*
============================================================
ODT が向いている問題集 + 解法コメント集
============================================================

前提:
- ODT は「同じ値が連続する最大区間」を set / map で持つ
- 基本操作は split(pos), assign(l,r,x), 区間走査
- 強いのは「区間代入」「区間塗り替え」「区間削除/復元」
- 弱いのは「最悪ケース保証」「assign がほぼ無い問題」

------------------------------------------------------------
[問題例 1] 区間代入 + 区間和
------------------------------------------------------------
問題:
- a[i] を持つ
- クエリ:
  1) [l,r] を x に代入
  2) sum(l,r) を答える

解法:
- assign(l,r,x) は ODT の本領
- sum は [l,r] にかかる各区間 [L,R] を走査し、
    res += v * (R-L+1)
- これは前回の Monoid なら RangeSumLL で fold できる

コメント:
- 「代入が多い」ならかなり相性がよい
- 「加算だけ多い」なら Fenwick / Lazy SegTree の方が普通は安全

------------------------------------------------------------
[問題例 2] 区間代入 + 区間 min / max / gcd / xor
------------------------------------------------------------
問題:
- [l,r] を x に代入
- 区間 min / max / gcd / xor を答える

解法:
- 走査してモノイドで畳み込む
- min/max/gcd は各色段 1 個ずつ見ればよい
- xor は長さの偶奇だけ注意すればよい
- 前回の alias:
    RangeMinLL / RangeMaxLL / RangeGcdLL / RangeXorLL

コメント:
- 「各色段ごとの寄与」が書けるなら実装しやすい
- 可変 mod が入ると通常の monoid からは外れるが、区間走査で直接計算はできる

------------------------------------------------------------
[問題例 3] 区間加算 + 区間代入 + k 番目
------------------------------------------------------------
問題:
- クエリ:
  1) [l,r] の各要素に x を加算
  2) [l,r] を x に代入
  3) [l,r] の k 番目に小さい値
  4) [l,r] で sum(a[i]^p mod mod) を求める

解法:
- ODT の代表例
- add は [l,r] にかかる各区間の値 v に対して v += x
- assign は 1 本化
- kth は [l,r] 内の各区間について (値, 個数) を集めて sort
- べき和は
    res += len * mod_pow(v, p, mod)
  を区間ごとに足す

コメント:
- これが ODT の教科書的な問題
- 代入で区間数がつぶれるので全体が持ちやすい
- ただし「ランダム性依存」を忘れない

------------------------------------------------------------
[問題例 4] 区間塗り替え + 色数 / 出現色集合
------------------------------------------------------------
問題:
- 各位置に色 c[i]
- [l,r] を色 x に塗る
- [l,r] に何色あるか / 色集合を答える

解法:
- ODT で色段を管理
- クエリ時に走査し、
    seen.insert(color)
  するだけ
- 色数が小さいなら bitset / ビットマスクで OR
- 例えば 64 色以下なら uint64_t で高速

コメント:
- 「色を塗る」系は ODT と非常に相性がよい
- 典型的な “色段” 問題

------------------------------------------------------------
[問題例 5] 区間 0/1 管理 + 1 の個数
------------------------------------------------------------
問題:
- 配列は 0/1
- [l,r] を 0 または 1 に代入
- [l,r] の 1 の個数を答える

解法:
- 値は 0/1 だけなので色段管理がしやすい
- 1 の個数は
    sum += (v == 1) * len
- 前回の alias なら CountOnes がそのまま使える

コメント:
- メモリ使用フラグ、座席占有、道路開通/閉鎖などにそのまま使える

------------------------------------------------------------
[問題例 6] 区間 0/1 管理 + 空き区間検索
------------------------------------------------------------
問題:
- 0 = 空き, 1 = 使用中
- [l,r] を 0/1 に代入
- 長さ k 以上の最左空き区間を探す

解法:
- ODT を左から走査
- 値 0 の区間を見て、長さ >= k ならその左端が答え
- 予約したら assign(ans, ans+k-1, 1)
- 解放は assign(l,r,0)

コメント:
- 典型的なメモリアロケータ風問題
- 先頭からの線形走査で足りる制約ならかなり簡潔
- 高速な first-fit を厳密に求めるならセグ木の方が本筋

------------------------------------------------------------
[問題例 7] 区間削除 / 区間復元
------------------------------------------------------------
問題:
- 各位置が alive / dead
- [l,r] を dead にする
- [l,r] を alive に戻す
- 生存数や生存区間を答える

解法:
- 実質 0/1 区間代入と同じ
- ODT で alive 段だけを走査して数える
- dead をまとめてつぶせるので実装しやすい

コメント:
- 「区間が連続して消える/戻る」系は適性が高い

------------------------------------------------------------
[問題例 8] 文字列の区間書き換え
------------------------------------------------------------
問題:
- 文字列 s
- [l,r] を文字 ch で埋める
- [l,r] に文字 x が何個あるか
- [l,r] の distinct 文字数を答える

解法:
- 各区間の値を char にする
- 個数は走査して
    if (v == ch) ans += len;
- distinct は set<char> や bitmask OR

コメント:
- 「文字の run-length 管理」と考えると ODT は自然
- ランレングス圧縮の動的版に近い

------------------------------------------------------------
[問題例 9] 座標が巨大だが変更区間が少ない問題
------------------------------------------------------------
問題:
- 座標が 1..1e18 のように大きい
- 実際に操作される境界は少数
- 区間代入や区間加算をする

解法:
- 必要なら座標圧縮してから ODT
- あるいは map ベース ODT で境界だけ持つ
- 「触った場所だけ区切る」ので sparse に扱える

コメント:
- 巨大座標に対しても「境界の数」が小さければ十分現実的
- ただしクエリ内容によっては圧縮 + Lazy SegTree の方が安定

------------------------------------------------------------
[問題例 10] 値ごとの寄与が区間長だけで決まる問題
------------------------------------------------------------
問題:
- [l,r] を走査し、
    Σ f(value, length_of_segment_part)
  の形で答えたい

解法:
- ODT で各色段ごとに寄与を足す
- 例えば:
    res += len
    res += value * len
    res += (value == x ? len : 0)
    res += len * powmod(value, p, mod)

コメント:
- ODT が本当に強いのはここ
- 「要素単位」ではなく「色段単位」に計算できると強い

------------------------------------------------------------
[問題例 11] 区間内の値の頻度分布がほしい
------------------------------------------------------------
問題:
- [l,r] 内で値ごとの出現回数を求めたい
- 最頻値、上位何個、しきい値以上の個数など

解法:
- 走査して map<value, cnt> に len を足す
- 最頻値なら最大 cnt を取る
- 上位 t 個なら map を vector 化して sort

コメント:
- 値が同じ連続区間が大きく保たれるならかなり有効
- 値がバラける worst case には弱い

------------------------------------------------------------
[問題例 12] 動的 RLE が欲しい問題
------------------------------------------------------------
問題:
- 配列/文字列を run-length encoding 的に扱いたい
- 更新は「区間を同じ値にする」が中心
- 問い合わせも run 単位で処理できる

解法:
- まさに ODT
- split して必要部分だけ切り出し
- assign で run を置き換え
- 前後と同値なら merge

コメント:
- 「RLE をオンライン更新したい」問題は ODT の本質そのもの

------------------------------------------------------------
[問題例 13] “見た目の塗りつぶし” 系シミュレーション
------------------------------------------------------------
問題:
- 棒、テープ、道路、画面の区間を色や状態で塗る
- 最終的な色の面積や色数を求める

解法:
- 各色段を管理
- 最終走査で色ごとの合計長を集計

コメント:
- 幾何というより 1 次元ペイント問題
- ODT で実装がかなり短くなる

------------------------------------------------------------
[問題例 14] 区間の反転が “値変換” として書ける 2 値問題
------------------------------------------------------------
問題:
- 0/1 配列
- [l,r] を反転 (0<->1)
- 1 の個数を答える

解法:
- [l,r] を走査し、各区間の値 v を 1-v に変える
- その後、前後 merge
- ODT の update(l,r,f) がそのまま使える

コメント:
- 値の種類が少ないほど扱いやすい
- ただし反転だけ連発で assign が無いと色段が増えやすい

------------------------------------------------------------
[問題例 15] 区間 chmin/chmax 風に見えるが、実は値集合が小さい問題
------------------------------------------------------------
問題:
- 値集合がとても小さい
- 区間の値を規則的に写像する
- 問い合わせは色段走査で出せる

解法:
- 各区間に対して v = trans[v] を適用
- 値域が小さいほど merge が起きやすい
- 実際には assign/paint に近い問題ほど向く

コメント:
- “値変換後に色段が潰れる” なら候補
- 潰れないなら危険

------------------------------------------------------------
[問題例 16] ODT を使わない方がよい問題
------------------------------------------------------------
問題:
- 区間加算だけ
- 区間和だけ
- 最悪 O(log n) を保証したい
- hack されやすい
- 値がほとんど毎回全部バラバラになる

解法:
- Fenwick Tree
- Lazy Segment Tree
- Segment Tree Beats
- 平衡 BST + augmented info
などを使う

コメント:
- ODT を “なんとなく” で選ばない
- assign/paint が本体でないなら別解が本命

============================================================
実装上の定石
============================================================

1. split(r+1) -> split(l) の順に切る
2. 区間全消去して assign 1 本挿入
3. 更新後は前後 merge
4. 集約は「要素ごと」ではなく「色段ごと」に書く
5. 値域が小さい / 代入が多い / run が潰れやすいほど有利

============================================================
前回の抽象 ODT にそのまま乗る代表モノイド
============================================================

- RangeSumLL   : 区間和
- RangeMinLL   : 区間最小
- RangeMaxLL   : 区間最大
- RangeGcdLL   : 区間 gcd
- RangeXorLL   : 区間 xor
- RangeOrULL   : 色集合 bitmask
- CountOnes    : 1 の個数

============================================================
面接・コンテスト中の判断基準
============================================================

次の 4 条件のうち 3 つ以上を満たしたら ODT 候補:

[ ] 区間代入 / 区間塗り替えがある
[ ] 問い合わせは区間を色段ごとに見ると簡単
[ ] 値がまとまりやすい / ランダム入力寄り
[ ] 最悪ケース保証がそこまで厳しくない

2 つ以下なら、たいていは Lazy SegTree 系を先に疑う。
*/
