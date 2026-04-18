# 目次

- [Wavelet Matrix](#wavelet-matrix)
  - [概要](#概要)
  - [記号](#記号)
  - [BitVector](#bitvector)
  - [WaveletMatrix&lt;T&gt;](#waveletmatrixt)
  - [型エイリアス](#型エイリアス)
- [Weighted Wavelet Matrix](#weighted-wavelet-matrix)
  - [概要](#概要-1)
  - [記号](#記号-1)
  - [WeightedWaveletMatrix](#weightedwaveletmatrix)
  - [RectangleSum](#rectanglesum)
  - [StaticRangeMonoid](#staticrangemonoid)
  - [CommutativeMonoidWaveletMatrix](#commutativemonoidwaveletmatrix)
  - [RectangleMonoid](#rectanglemonoid)
  - [rect_monoid::Sum](#rect_monoidsum)
  - [rect_monoid::Product](#rect_monoidproduct)
  - [rect_monoid::Min](#rect_monoidmin)
  - [rect_monoid::Max](#rect_monoidmax)
  - [rect_monoid::Gcd](#rect_monoidgcd)
  - [rect_monoid::Lcm](#rect_monoidlcm)
  - [rect_monoid::Xor](#rect_monoidxor)
  - [rect_monoid::BitAnd](#rect_monoidbitand)
  - [rect_monoid::BitOr](#rect_monoidbitor)
  - [rect_monoid::AddMod998244353](#rect_monoidaddmod998244353)
  - [rect_monoid::MulMod998244353](#rect_monoidmulmod998244353)
  - [型エイリアス](#型エイリアス-1)

# Wavelet Matrix
[ライブラリのリンク](https://github.com/harurunrunrun/GPT-generated-Library/blob/main/library/WaveletMatrixMultifunction.hpp)

## 概要

- `WaveletMatrix<T>` は**静的配列向け**で、`build` 後の更新はできません。
- bitwise 系 API は **`T` が整数型のときのみ**使えます。
- sum 系 API は **`T` が `long long` に変換可能なときのみ**使えます。

## 記号

- `N`: 配列長
- `σ`: 異なる値の個数
- `W`: `unsigned(T)` のビット幅
- `M`: bitwise OR / AND 系で探索中に訪れるノード数（コメント上では最悪 `O(2^W)`）
- `m`: 出力される異なる値の個数

## BitVector

### 構築・準備

| 関数   | 時間計算量 | 説明 |
|---|---|---|
| `BitVector() = default` | O(1) | 空の BitVector を構築する。 |
| `explicit BitVector(int n)` | O(n / 64) | 長さ n の BitVector を構築し、init(n) を呼ぶ。 |
| `void init(int n_)` | O(n / 64) | 長さ n の空ビット列を作る。 |
| `void build()` | O(n / 64) | rank 用の累積情報を構築する。 set() を終えたあとに 1 回呼ぶ。 |

### 参照・集計

| 関数   | 時間計算量 | 説明 |
|---|---|---|
| `int size() const` | O(1) | ビット列の長さを返す。 |
| `void set(int i)` | O(1) | i 番目のビットを 1 にする。 |
| `bool access(int i) const` | O(1) | i 番目のビットを返す。 |
| `int rank1(int r) const` | O(1) | 区間 [0, r) に含まれる 1 の個数を返す。 |
| `int rank0(int r) const` | O(1) | 区間 [0, r) に含まれる 0 の個数を返す。 |
| `int rank(bool bit, int r) const` | O(1) | 区間 [0, r) に含まれる bit の個数を返す。 |
| `int rank1(int l, int r) const` | O(1) | 区間 [l, r) に含まれる 1 の個数を返す。 |
| `int rank0(int l, int r) const` | O(1) | 区間 [l, r) に含まれる 0 の個数を返す。 |

### 位置検索

| 関数   | 時間計算量 | 説明 |
|---|---|---|
| `int select1(int k) const` | O(log n) / ※ rank を用いた二分探索 | 0-indexed で k 番目の 1 がある位置を返す。 存在しなければ -1 を返す。 |
| `int select0(int k) const` | O(log n) | 0-indexed で k 番目の 0 がある位置を返す。 存在しなければ -1 を返す。 |
| `int select(bool bit, int k) const` | O(log n) | 0-indexed で k 番目の bit がある位置を返す。 存在しなければ -1 を返す。 |

## WaveletMatrix<T>

### 基本操作・メタ情報


| 関数   | 時間計算量 | 説明 |
|---|---|---|
| `WaveletMatrix() = default` | O(1) | 空の WaveletMatrix を構築する。 |
| `explicit WaveletMatrix(const std::vector<T>& data)` | O(N log N + N log σ) / ※ 座標圧縮の sort を含む | 配列 data から Wavelet Matrix を構築する。 |
| `void build(const std::vector<T>& data)` | O(N log N + N log σ) | 配列 data から再構築する。 |
| `int size() const` | O(1) | 元の配列長 N を返す。 |
| `bool empty() const` | O(1) | 空かどうかを返す。 |
| `int distinct_size() const` | O(1) | 異なる値の個数 σ を返す。 |
| `int bit_size() const` | O(1) | 内部で使っているビット長を返す。 |
| `const std::vector<T>& values() const` | O(1) | 座標圧縮後の「昇順に並んだ異なる値一覧」を返す。 |
| `const T& access(int i) const` | O(1) | 元配列の i 番目の値を返す。 |
| `const T& operator[](int i) const` | O(1) | access(i) の別名。 |
| `bool contains(const T& value) const` | O(log σ) | value が配列全体に 1 回以上出現するかを返す。 |
| `int index_of(const T& value) const` | O(log σ) | 座標圧縮後の ID を返す。 存在しなければ -1 を返す。 |

### 出現回数・位置

| 関数   | 時間計算量 | 説明 |
|---|---|---|
| `int rank(const T& value, int r) const` | O(log σ) | 区間 [0, r) に含まれる value の個数を返す。 |
| `int rank(const T& value, int l, int r) const` | O(log σ) | 区間 [l, r) に含まれる value の個数を返す。 |
| `int count(const T& value) const` | O(log σ) | 配列全体に含まれる value の個数を返す。 |
| `int count(const T& value, int l, int r) const` | O(log σ) | 区間 [l, r) に含まれる value の個数を返す。 |
| `int select(const T& value, int kth) const` | O(log σ log N) / ※ 各レベルで BitVector::select を使う | value の kth 回目の出現位置を返す。 kth は 0-indexed。 存在しなければ -1 を返す。 |

### 順序統計・代表値

| 関数   | 時間計算量 | 説明 |
|---|---|---|
| `const T& kth_smallest(int l, int r, int k) const` | O(log σ) | 区間 [l, r) の k 番目に小さい値を返す。 k は 0-indexed。 |
| `const T& kth_largest(int l, int r, int k) const` | O(log σ) | 区間 [l, r) の k 番目に大きい値を返す。 k は 0-indexed。 |
| `const T& quantile(int l, int r, int k) const` | O(log σ) | kth_smallest の別名。 |
| `std::optional<T> min_value(int l, int r) const` | O(log σ) | 区間 [l, r) の最小値を返す。 空区間なら std::nullopt を返す。 |
| `std::optional<T> max_value(int l, int r) const` | O(log σ) | 区間 [l, r) の最大値を返す。 空区間なら std::nullopt を返す。 |
| `std::optional<T> median_lower(int l, int r) const` | O(log σ) | 区間 [l, r) の下側中央値を返す。 要素数 m に対して floor((m-1)/2) 番目。 空区間なら std::nullopt。 |
| `std::optional<T> median_upper(int l, int r) const` | O(log σ) | 区間 [l, r) の上側中央値を返す。 要素数 m に対して floor(m/2) 番目。 空区間なら std::nullopt。 |

### 範囲カウント

| 関数   | 時間計算量 | 説明 |
|---|---|---|
| `int range_freq(int l, int r, const T& upper) const` | O(log σ) | 区間 [l, r) の要素 x について、x < upper である要素数を返す。 |
| `int range_freq(int l, int r, const T& lower, const T& upper) const` | O(log σ) | 区間 [l, r) の要素 x について、lower <= x < upper である要素数を返す。 |
| `int count_less(int l, int r, const T& upper) const` | O(log σ) | 区間 [l, r) の要素 x について、x < upper である要素数を返す。 |
| `int count_less_equal(int l, int r, const T& upper) const` | O(log σ) | 区間 [l, r) の要素 x について、x <= upper である要素数を返す。 |
| `int count_greater(int l, int r, const T& lower) const` | O(log σ) | 区間 [l, r) の要素 x について、x > lower である要素数を返す。 |
| `int count_greater_equal(int l, int r, const T& lower) const` | O(log σ) | 区間 [l, r) の要素 x について、x >= lower である要素数を返す。 |

### 前後探索

| 関数   | 時間計算量 | 説明 |
|---|---|---|
| `std::optional<T> next_value_ge(int l, int r, const T& lower) const` | O(log σ) | 区間 [l, r) で lower 以上の最小値を返す。 存在しなければ std::nullopt。 |
| `std::optional<T> next_value_gt(int l, int r, const T& lower) const` | O(log σ) | 区間 [l, r) で lower より大きい最小値を返す。 存在しなければ std::nullopt。 |
| `std::optional<T> prev_value_lt(int l, int r, const T& upper) const` | O(log σ) | 区間 [l, r) で upper 未満の最大値を返す。 存在しなければ std::nullopt。 |
| `std::optional<T> prev_value_le(int l, int r, const T& upper) const` | O(log σ) | 区間 [l, r) で upper 以下の最大値を返す。 存在しなければ std::nullopt。 |
| `std::optional<T> next_value(int l, int r, const T& lower) const` | O(log σ) | next_value_ge の別名。 |
| `std::optional<T> prev_value(int l, int r, const T& upper) const` | O(log σ) | prev_value_lt の別名。 |

### 総和クエリ

| 関数   | 時間計算量 | 説明 |
|---|---|---|
| `std::pair<int, sum_type> count_and_sum_less(int l, int r, const T& upper) const` | O(log σ) | 区間 [l, r) の要素 x について、x < upper である要素の (個数, 総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_equal(int l, int r, const T& upper) const` | O(log σ) | 区間 [l, r) の要素 x について、x <= upper である要素の (個数, 総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater(int l, int r, const T& lower) const` | O(log σ) | 区間 [l, r) の要素 x について、x > lower である要素の (個数, 総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_equal(int l, int r, const T& lower) const` | O(log σ) | 区間 [l, r) の要素 x について、x >= lower である要素の (個数, 総和) を返す。 |
| `sum_type sum_all(int l, int r) const` | O(1) | x が xl <= x < xr を満たす点の重み和を返す。 |
| `sum_type sum_less(int l, int r, const T& upper) const` | O(log σ) | 区間 [l, r) において x < upper の要素の総和を返す。 |
| `sum_type sum_less_equal(int l, int r, const T& upper) const` | O(log σ) | 区間 [l, r) において x <= upper の要素の総和を返す。 |
| `sum_type sum_range(int l, int r, const T& lower, const T& upper) const` | O(log σ) | 区間 [l, r) において lower <= x < upper の要素の総和を返す。 |
| `sum_type sum_equal(int l, int r, const T& value) const` | O(log σ) | 区間 [l, r) における value の総和を返す。 |
| `sum_type sum_greater(int l, int r, const T& lower) const` | O(log σ) | 区間 [l, r) において x > lower の要素の総和を返す。 |
| `sum_type sum_greater_equal(int l, int r, const T& lower) const` | O(log σ) | 区間 [l, r) において x >= lower の要素の総和を返す。 |
| `sum_type sum_k_smallest(int l, int r, int k) const` | O(log σ) | 区間 [l, r) の小さい方から k 個の要素の総和を返す。 k は個数であり、0 <= k <= r-l を満たす必要がある。 |
| `sum_type sum_k_largest(int l, int r, int k) const` | O(log σ) | 区間 [l, r) の大きい方から k 個の要素の総和を返す。 k は個数であり、0 <= k <= r-l を満たす必要がある。 |

### 列挙・頻度分析


| 関数   | 時間計算量 | 説明 |
|---|---|---|
| `std::vector<std::pair<T, int>> top_k_frequent(int l, int r, int k) const` | おおよそ O(k log σ log(k log σ)) / ※ 優先度付きキューを使う | 区間 [l, r) で頻度上位 k 個の (値, 出現回数) を返す。 |
| `std::optional<std::pair<T, int>> mode(int l, int r) const` | top_k_frequent(..., 1) に依存 / おおよそ O(log σ log log σ) 〜 O(log σ log σ) / 実装上は top_k_frequent を利用 | 区間 [l, r) の最頻値を (値, 出現回数) で返す。 空区間なら std::nullopt。 |
| `std::vector<std::tuple<T, int, int>> intersect(int l1, int r1, int l2, int r2) const` | O(z log σ) / z は共通する異なる値の個数 | 2 区間 [l1, r1), [l2, r2) に共通して現れる値を列挙する。 |
| `std::vector<std::pair<T, int>> list_frequencies(int l, int r) const` | おおよそ O((m + 1) log σ) / m は出力される異なる値の個数 | 区間 [l, r) に現れる異なる値を (値, 個数) で昇順列挙する。 |
| `std::vector<std::pair<T, int>> list_frequencies(int l, int r, const T& lower, const T& upper) const` | おおよそ O((m + 1) log σ) / m は出力される異なる値の個数 | 区間 [l, r) について、値域 [lower, upper) に属する 異なる値を (値, 個数) で昇順列挙する。 |
| `std::vector<T> distinct_values(int l, int r) const` | list_frequencies に依存 / おおよそ O((m + 1) log σ) | 区間 [l, r) に現れる異なる値を昇順で返す。 |

### Bitwise 共通 API


| 関数   | 時間計算量 | 説明 |
|---|---|---|
| `int count_less_bitwise(int l, int r, unsigned long long mask, unsigned long long upper, BitwiseOperation op) const` | XOR のとき O(W) / OR / AND のとき O(M) / W は unsigned(T) のビット幅、M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、y < upper である要素数を返す。 |
| `int range_freq_bitwise( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper, BitwiseOperation op ) const` | count_less_bitwise を 2 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、lower <= y < upper である要素数を返す。 |
| `int count_less_equal_bitwise(int l, int r, unsigned long long mask, unsigned long long upper, BitwiseOperation op) const` | count_less_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、y <= upper である要素数を返す。 |
| `int count_equal_bitwise(int l, int r, unsigned long long mask, unsigned long long value, BitwiseOperation op) const` | count_less_bitwise を 2 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、y == value である要素数を返す。 |
| `int count_greater_bitwise(int l, int r, unsigned long long mask, unsigned long long lower, BitwiseOperation op) const` | count_less_equal_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、y > lower である要素数を返す。 |
| `int count_greater_equal_bitwise(int l, int r, unsigned long long mask, unsigned long long lower, BitwiseOperation op) const` | count_less_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、y >= lower である要素数を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_bitwise( int l, int r, unsigned long long mask, unsigned long long upper, BitwiseOperation op ) const` | XOR のとき O(W) / OR / AND のとき O(M) / W は unsigned(T) のビット幅、M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、y < upper である要素について、(個数, 元の x の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_equal_bitwise( int l, int r, unsigned long long mask, unsigned long long upper, BitwiseOperation op ) const` | count_and_sum_less_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、y <= upper である要素について、(個数, 元の x の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_range_bitwise( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper, BitwiseOperation op ) const` | count_and_sum_less_bitwise を 2 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、lower <= y < upper である要素について、(個数, 元の x の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_bitwise( int l, int r, unsigned long long mask, unsigned long long lower, BitwiseOperation op ) const` | count_and_sum_less_equal_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、y > lower である要素について、(個数, 元の x の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_equal_bitwise( int l, int r, unsigned long long mask, unsigned long long lower, BitwiseOperation op ) const` | count_and_sum_less_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、y >= lower である要素について、(個数, 元の x の総和) を返す。 |
| `sum_type sum_less_bitwise( int l, int r, unsigned long long mask, unsigned long long upper, BitwiseOperation op ) const` | count_and_sum_less_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、y < upper である要素の元の x の総和を返す。 |
| `sum_type sum_less_equal_bitwise( int l, int r, unsigned long long mask, unsigned long long upper, BitwiseOperation op ) const` | count_and_sum_less_equal_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、y <= upper である要素の元の x の総和を返す。 |
| `sum_type sum_equal_bitwise( int l, int r, unsigned long long mask, unsigned long long value, BitwiseOperation op ) const` | count_and_sum_less_bitwise を 2 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、y == value である要素の元の x の総和を返す。 |
| `sum_type sum_range_bitwise( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper, BitwiseOperation op ) const` | count_and_sum_less_bitwise を 2 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、lower <= y < upper である要素の元の x の総和を返す。 |
| `sum_type sum_greater_bitwise( int l, int r, unsigned long long mask, unsigned long long lower, BitwiseOperation op ) const` | count_and_sum_greater_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、y > lower である要素の元の x の総和を返す。 |
| `sum_type sum_greater_equal_bitwise( int l, int r, unsigned long long mask, unsigned long long lower, BitwiseOperation op ) const` | count_and_sum_greater_equal_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、y >= lower である要素の元の x の総和を返す。 |

### Bitwise 共通 API（変換後の値に対する検索・列挙）


| 関数   | 時間計算量 | 説明 |
|---|---|---|
| `unsigned long long access_bitwise(int i, unsigned long long mask, BitwiseOperation op) const` | O(1) | i 番目の要素 value に bitwise 演算を施した結果 op(value, mask) を返す。 |
| `int rank_bitwise(unsigned long long value, int r, unsigned long long mask, BitwiseOperation op) const` | count_equal_bitwise に準ずる | prefix 区間 [0, r) の各要素 x に対して y = op(x, mask) としたとき、y == value である要素数を返す。 |
| `int rank_bitwise(unsigned long long value, int l, int r, unsigned long long mask, BitwiseOperation op) const` | count_equal_bitwise に準ずる | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、y == value である要素数を返す。 |
| `int count_bitwise(unsigned long long value, unsigned long long mask, BitwiseOperation op) const` | count_equal_bitwise に準ずる | 配列全体の各要素 x に対して y = op(x, mask) としたとき、y == value である要素数を返す。 |
| `int count_bitwise(unsigned long long value, int l, int r, unsigned long long mask, BitwiseOperation op) const` | count_equal_bitwise に準ずる | 区間 [l, r) の各要素 x に対して y = op(x, mask) としたとき、y == value である要素数を返す。 |
| `bool contains_bitwise(unsigned long long value, unsigned long long mask, BitwiseOperation op) const` | count_bitwise に準ずる | 配列全体の各要素 x に対して y = op(x, mask) としたとき、y == value である要素が存在するかを返す。 |
| `std::vector<unsigned long long> values_bitwise(unsigned long long mask, BitwiseOperation op) const` | list_frequencies_bitwise(0, n_, mask, op) に準ずる | 配列全体において現れる op(element, mask) の異なる値を昇順で返す。 |
| `int distinct_size_bitwise(unsigned long long mask, BitwiseOperation op) const` | values_bitwise(mask, op) に準ずる | 配列全体において現れる op(element, mask) の異なる値の個数を返す。 |
| `int index_of_bitwise(unsigned long long value, unsigned long long mask, BitwiseOperation op) const` | values_bitwise(mask, op) に準ずる | values_bitwise(mask, op) が返す昇順 distinct 配列の中で、 value が現れる位置を返す。存在しなければ -1 を返す。 |
| `int select_bitwise(unsigned long long value, int kth, unsigned long long mask, BitwiseOperation op) const` | O(log N * C) / C は rank_bitwise(value, r, mask, op) 1 回分 | 0-indexed で kth 番目に現れる、y = op(x, mask) が y == value となる位置を返す。存在しなければ -1 を返す。 |
| `unsigned long long kth_smallest_bitwise( int l, int r, int k, unsigned long long mask, BitwiseOperation op ) const` | O(W * C) / W は unsigned(T) のビット幅、 / C は count_less_equal_bitwise 1 回分 | 区間 [l, r) の op(element, mask) たちを昇順に見たときの k 番目の値を返す。k は 0-indexed。 |
| `unsigned long long kth_largest_bitwise( int l, int r, int k, unsigned long long mask, BitwiseOperation op ) const` | kth_smallest_bitwise に準ずる | 区間 [l, r) の op(element, mask) たちを降順に見たときの k 番目の値を返す。k は 0-indexed。 |
| `unsigned long long quantile_bitwise( int l, int r, int k, unsigned long long mask, BitwiseOperation op ) const` | kth_smallest_bitwise に準ずる | kth_smallest_bitwise の別名。 |
| `std::optional<unsigned long long> min_value_bitwise( int l, int r, unsigned long long mask, BitwiseOperation op ) const` | kth_smallest_bitwise に準ずる | 区間 [l, r) の op(element, mask) の最小値を返す。 空区間なら std::nullopt を返す。 |
| `std::optional<unsigned long long> max_value_bitwise( int l, int r, unsigned long long mask, BitwiseOperation op ) const` | kth_largest_bitwise に準ずる | 区間 [l, r) の op(element, mask) の最大値を返す。 空区間なら std::nullopt を返す。 |
| `std::optional<unsigned long long> median_lower_bitwise( int l, int r, unsigned long long mask, BitwiseOperation op ) const` | kth_smallest_bitwise に準ずる | 区間 [l, r) の op(element, mask) の下側中央値を返す。 空区間なら std::nullopt を返す。 |
| `std::optional<unsigned long long> median_upper_bitwise( int l, int r, unsigned long long mask, BitwiseOperation op ) const` | kth_smallest_bitwise に準ずる | 区間 [l, r) の op(element, mask) の上側中央値を返す。 空区間なら std::nullopt を返す。 |
| `std::optional<unsigned long long> next_value_ge_bitwise( int l, int r, unsigned long long lower, unsigned long long mask, BitwiseOperation op ) const` | count_less_bitwise と kth_smallest_bitwise に準ずる | 区間 [l, r) の op(element, mask) の中で、 lower 以上の最小値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> next_value_gt_bitwise( int l, int r, unsigned long long lower, unsigned long long mask, BitwiseOperation op ) const` | count_less_equal_bitwise と kth_smallest_bitwise に準ずる | 区間 [l, r) の op(element, mask) の中で、 lower より大きい最小値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> prev_value_lt_bitwise( int l, int r, unsigned long long upper, unsigned long long mask, BitwiseOperation op ) const` | count_less_bitwise と kth_smallest_bitwise に準ずる | 区間 [l, r) の op(element, mask) の中で、 upper 未満の最大値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> prev_value_le_bitwise( int l, int r, unsigned long long upper, unsigned long long mask, BitwiseOperation op ) const` | count_less_equal_bitwise と kth_smallest_bitwise に準ずる | 区間 [l, r) の op(element, mask) の中で、 upper 以下の最大値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> next_value_bitwise( int l, int r, unsigned long long lower, unsigned long long mask, BitwiseOperation op ) const` | next_value_ge_bitwise に準ずる | next_value_ge_bitwise の別名。 |
| `std::optional<unsigned long long> prev_value_bitwise( int l, int r, unsigned long long upper, unsigned long long mask, BitwiseOperation op ) const` | prev_value_lt_bitwise に準ずる | prev_value_lt_bitwise の別名。 |
| `std::vector<std::pair<unsigned long long, int>> list_frequencies_bitwise( int l, int r, unsigned long long mask, BitwiseOperation op ) const` | O((m + 1) log σ + m log m) / m は区間 [l, r) に現れる異なる元値の個数 | 区間 [l, r) の op(element, mask) に現れる異なる値を (値, 個数) で昇順列挙する。 |
| `std::vector<std::pair<unsigned long long, int>> list_frequencies_bitwise( int l, int r, unsigned long long lower, unsigned long long upper, unsigned long long mask, BitwiseOperation op ) const` | O((m + 1) log σ + m log m) / m は区間 [l, r) に現れる異なる元値の個数 | 区間 [l, r) の op(element, mask) について、 値域 [lower, upper) に属する異なる値を (値, 個数) で昇順列挙する。 |
| `std::vector<unsigned long long> distinct_values_bitwise( int l, int r, unsigned long long mask, BitwiseOperation op ) const` | list_frequencies_bitwise(l, r, mask, op) に準ずる | 区間 [l, r) の op(element, mask) に現れる異なる値を 昇順で返す。 |
| `std::vector<std::pair<unsigned long long, int>> top_k_frequent_bitwise( int l, int r, int k, unsigned long long mask, BitwiseOperation op ) const` | O((m + 1) log σ + m log m) / m は区間 [l, r) に現れる異なる元値の個数 | 区間 [l, r) の op(element, mask) に対し、 頻度上位 k 個の (値, 出現回数) を返す。 |
| `std::optional<std::pair<unsigned long long, int>> mode_bitwise( int l, int r, unsigned long long mask, BitwiseOperation op ) const` | top_k_frequent_bitwise に準ずる | 区間 [l, r) の op(element, mask) の最頻値を (値, 出現回数) で返す。空区間なら std::nullopt。 |
| `std::vector<std::tuple<unsigned long long, int, int>> intersect_bitwise( int l1, int r1, int l2, int r2, unsigned long long mask, BitwiseOperation op ) const` | O((m1 + m2 + 1) log σ + m1 log m1 + m2 log m2) / m1, m2 は各区間に現れる異なる元値の個数 | 2 区間 [l1, r1), [l2, r2) について、 共通して現れる op(element, mask) の値を列挙する。 |

### XOR ラッパー


| 関数   | 時間計算量 | 説明 |
|---|---|---|
| `int count_less_xor(int l, int r, unsigned long long mask, unsigned long long upper) const` | O(W) | 区間 [l, r) において (value xor mask) < upper の個数を返す。 |
| `int range_freq_xor( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper ) const` | O(W) | 区間 [l, r) において lower <= (value xor mask) < upper の個数を返す。 |
| `int count_less_equal_xor(int l, int r, unsigned long long mask, unsigned long long upper) const` | O(W) | 区間 [l, r) において (value xor mask) <= upper の個数を返す。 |
| `int count_equal_xor(int l, int r, unsigned long long mask, unsigned long long value) const` | O(W) | 区間 [l, r) の各要素 x に対して y = (x xor mask) としたとき、y == value である要素数を返す。 |
| `int count_greater_xor(int l, int r, unsigned long long mask, unsigned long long lower) const` | O(W) | 区間 [l, r) において (value xor mask) > lower の個数を返す。 |
| `int count_greater_equal_xor(int l, int r, unsigned long long mask, unsigned long long lower) const` | O(W) | 区間 [l, r) において (value xor mask) >= lower の個数を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_xor( int l, int r, unsigned long long mask, unsigned long long upper ) const` | O(W) | 区間 [l, r) の各要素 x に対して y = (x xor mask) としたとき、y < upper である要素について、(個数, 元の x の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_equal_xor( int l, int r, unsigned long long mask, unsigned long long upper ) const` | O(W) | 区間 [l, r) の各要素 x に対して y = (x xor mask) としたとき、y <= upper である要素について、(個数, 元の x の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_range_xor( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper ) const` | O(W) | 区間 [l, r) の各要素 x に対して y = (x xor mask) としたとき、lower <= y < upper である要素について、(個数, 元の x の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_xor( int l, int r, unsigned long long mask, unsigned long long lower ) const` | O(W) | 区間 [l, r) の各要素 x に対して y = (x xor mask) としたとき、y > lower である要素について、(個数, 元の x の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_equal_xor( int l, int r, unsigned long long mask, unsigned long long lower ) const` | O(W) | 区間 [l, r) の各要素 x に対して y = (x xor mask) としたとき、y >= lower である要素について、(個数, 元の x の総和) を返す。 |
| `sum_type sum_less_xor(int l, int r, unsigned long long mask, unsigned long long upper) const` | O(W) | 区間 [l, r) の各要素 x に対して y = (x xor mask) としたとき、y < upper である要素の元の x の総和を返す。 |
| `sum_type sum_less_equal_xor(int l, int r, unsigned long long mask, unsigned long long upper) const` | O(W) | 区間 [l, r) の各要素 x に対して y = (x xor mask) としたとき、y <= upper である要素の元の x の総和を返す。 |
| `sum_type sum_equal_xor(int l, int r, unsigned long long mask, unsigned long long value) const` | O(W) | 区間 [l, r) の各要素 x に対して y = (x xor mask) としたとき、y == value である要素の元の x の総和を返す。 |
| `sum_type sum_range_xor( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper ) const` | O(W) | 区間 [l, r) の各要素 x に対して y = (x xor mask) としたとき、lower <= y < upper である要素の元の x の総和を返す。 |
| `sum_type sum_greater_xor(int l, int r, unsigned long long mask, unsigned long long lower) const` | O(W) | 区間 [l, r) の各要素 x に対して y = (x xor mask) としたとき、y > lower である要素の元の x の総和を返す。 |
| `sum_type sum_greater_equal_xor(int l, int r, unsigned long long mask, unsigned long long lower) const` | O(W) | 区間 [l, r) の各要素 x に対して y = (x xor mask) としたとき、y >= lower である要素の元の x の総和を返す。 |
| `unsigned long long access_xor(int i, unsigned long long mask) const` | O(1) | i 番目の要素に (value xor mask) を施した結果を返す。 |
| `int rank_xor(unsigned long long value, int r, unsigned long long mask) const` | O(log σ) | prefix 区間 [0, r) の各要素 x に対して y = (x xor mask) としたとき、y == value である要素数を返す。 |
| `int rank_xor(unsigned long long value, int l, int r, unsigned long long mask) const` | O(log σ) | 区間 [l, r) の各要素 x に対して y = (x xor mask) としたとき、y == value である要素数を返す。 |
| `int count_xor(unsigned long long value, unsigned long long mask) const` | O(log σ) | 配列全体の各要素 x に対して y = (x xor mask) としたとき、y == value である要素数を返す。 |
| `int count_xor(unsigned long long value, int l, int r, unsigned long long mask) const` | O(log σ) | 区間 [l, r) の各要素 x に対して y = (x xor mask) としたとき、y == value である要素数を返す。 |
| `bool contains_xor(unsigned long long value, unsigned long long mask) const` | O(log σ) | 配列全体の各要素 x に対して y = (x xor mask) としたとき、y == value である要素が存在するかを返す。 |
| `std::vector<unsigned long long> values_xor(unsigned long long mask) const` | おおよそ O((σ + 1)W) / W は unsigned(T) のビット幅 | 配列全体において現れる (element xor mask) の異なる値を昇順で返す。 |
| `int distinct_size_xor(unsigned long long mask) const` | O(1) | 配列全体において現れる (element xor mask) の異なる値の個数を返す。 XOR は全単射なので元の distinct 数と一致する。 |
| `int index_of_xor(unsigned long long value, unsigned long long mask) const` | おおよそ O((σ + 1)W) / W は unsigned(T) のビット幅 | values_xor(mask) が返す昇順 distinct 配列の中で、value が現れる位置を返す。存在しなければ -1 を返す。 |
| `int select_xor(unsigned long long value, int kth, unsigned long long mask) const` | O(log σ log N) | 0-indexed で kth 番目に現れる、y = (x xor mask) が y == value となる位置を返す。存在しなければ -1 を返す。 |
| `unsigned long long kth_smallest_xor(int l, int r, int k, unsigned long long mask) const` | O(W) | 区間 [l, r) の各要素 x に対して z = (x ^ mask) とした値を昇順に見たときの k 番目の値を返す。k は 0-indexed。 |
| `unsigned long long kth_largest_xor(int l, int r, int k, unsigned long long mask) const` | O(W) | 区間 [l, r) の各要素 x に対して z = (x ^ mask) とした値を降順に見たときの k 番目の値を返す。k は 0-indexed。 |
| `unsigned long long quantile_xor(int l, int r, int k, unsigned long long mask) const` | O(W) | kth_smallest_xor の別名。 |
| `std::optional<unsigned long long> min_value_xor(int l, int r, unsigned long long mask) const` | O(W) | 区間 [l, r) の (element xor mask) の最小値を返す。空区間なら std::nullopt。 |
| `std::optional<unsigned long long> max_value_xor(int l, int r, unsigned long long mask) const` | O(W) | 区間 [l, r) の (element xor mask) の最大値を返す。空区間なら std::nullopt。 |
| `std::optional<unsigned long long> median_lower_xor(int l, int r, unsigned long long mask) const` | O(W) | 区間 [l, r) の各要素 x に対して z = (x ^ mask) とした値の下側中央値を返す。空区間なら std::nullopt を返す。 |
| `std::optional<unsigned long long> median_upper_xor(int l, int r, unsigned long long mask) const` | O(W) | 区間 [l, r) の各要素 x に対して z = (x ^ mask) とした値の上側中央値を返す。空区間なら std::nullopt を返す。 |
| `std::optional<unsigned long long> next_value_ge_xor(int l, int r, unsigned long long lower, unsigned long long mask) const` | O(W) | 区間 [l, r) の (element xor mask) の中で、lower 以上の最小値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> next_value_gt_xor(int l, int r, unsigned long long lower, unsigned long long mask) const` | O(W) | 区間 [l, r) の (element xor mask) の中で、lower より大きい最小値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> prev_value_lt_xor(int l, int r, unsigned long long upper, unsigned long long mask) const` | O(W) | 区間 [l, r) の (element xor mask) の中で、upper 未満の最大値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> prev_value_le_xor(int l, int r, unsigned long long upper, unsigned long long mask) const` | O(W) | 区間 [l, r) の (element xor mask) の中で、upper 以下の最大値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> next_value_xor(int l, int r, unsigned long long lower, unsigned long long mask) const` | O(W) | next_value_ge_xor の別名。 |
| `std::optional<unsigned long long> prev_value_xor(int l, int r, unsigned long long upper, unsigned long long mask) const` | O(W) | prev_value_lt_xor の別名。 |
| `std::vector<std::pair<unsigned long long, int>> list_frequencies_xor(int l, int r, unsigned long long mask) const` | おおよそ O((z + 1)W) / z は出力される異なる値の個数、W は unsigned(T) のビット幅 | 区間 [l, r) の各要素 x に対して z = (x ^ mask) とした値について、(値, 個数) を昇順に列挙する。変換後に同じ値になったものは集約する。 |
| `std::vector<std::pair<unsigned long long, int>> list_frequencies_xor(int l, int r, unsigned long long lower, unsigned long long upper, unsigned long long mask) const` | おおよそ O((z + 1)W) / z は出力される異なる値の個数、W は unsigned(T) のビット幅 | 区間 [l, r) の各要素 x に対して z = (x ^ mask) とした値について、(値, 個数) を昇順に列挙する。変換後に同じ値になったものは集約する。 |
| `std::vector<unsigned long long> distinct_values_xor(int l, int r, unsigned long long mask) const` | list_frequencies_xor(l, r, mask) に準ずる | 区間 [l, r) の各要素 x に対して z = (x ^ mask) としたときに現れる異なる値を昇順に列挙する。 |
| `std::vector<std::pair<unsigned long long, int>> top_k_frequent_xor(int l, int r, int k, unsigned long long mask) const` | O((m + 1) log σ + m log k + k log k) / m は区間 [l, r) に現れる異なる元値の個数 | 区間 [l, r) の (element xor mask) に対し、頻度上位 k 個の (値, 出現回数) を返す。 |
| `std::optional<std::pair<unsigned long long, int>> mode_xor(int l, int r, unsigned long long mask) const` | O((m + 1) log σ) / m は区間 [l, r) に現れる異なる元値の個数 | 区間 [l, r) の (element xor mask) の最頻値を (値, 出現回数) で返す。空区間なら std::nullopt。 |
| `std::vector<std::tuple<unsigned long long, int, int>> intersect_xor(int l1, int r1, int l2, int r2, unsigned long long mask) const` | おおよそ O((z + 1)W) / z は共通して現れる異なる値の個数、W は unsigned(T) のビット幅 | 2 区間 [l1, r1), [l2, r2) について、共通して現れる (element xor mask) の値を列挙する。 |

### OR ラッパー

| 関数   | 時間計算量 | 説明 |
|---|---|---|
| `int count_less_or(int l, int r, unsigned long long mask, unsigned long long upper) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value \| mask) < upper の個数を返す。 |
| `int range_freq_or( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper ) const` | count_less_or を 2 回呼ぶので、それに準ずる | 区間 [l, r) において lower <= (value \| mask) < upper の個数を返す。 |
| `int count_less_equal_or(int l, int r, unsigned long long mask, unsigned long long upper) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value \| mask) <= upper の個数を返す。 |
| `int count_equal_or(int l, int r, unsigned long long mask, unsigned long long value) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して y = (x \| mask) としたとき、y == value である要素数を返す。 |
| `int count_greater_or(int l, int r, unsigned long long mask, unsigned long long lower) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value \| mask) > lower の個数を返す。 |
| `int count_greater_equal_or(int l, int r, unsigned long long mask, unsigned long long lower) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value \| mask) >= lower の個数を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_or( int l, int r, unsigned long long mask, unsigned long long upper ) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して y = (x \| mask) としたとき、y < upper である要素について、(個数, 元の x の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_equal_or( int l, int r, unsigned long long mask, unsigned long long upper ) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して y = (x \| mask) としたとき、y <= upper である要素について、(個数, 元の x の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_range_or( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper ) const` | count_and_sum_less_or を 2 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = (x \| mask) としたとき、lower <= y < upper である要素について、(個数, 元の x の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_or( int l, int r, unsigned long long mask, unsigned long long lower ) const` | count_and_sum_less_equal_or を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = (x \| mask) としたとき、y > lower である要素について、(個数, 元の x の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_equal_or( int l, int r, unsigned long long mask, unsigned long long lower ) const` | count_and_sum_less_or を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = (x \| mask) としたとき、y >= lower である要素について、(個数, 元の x の総和) を返す。 |
| `sum_type sum_less_or(int l, int r, unsigned long long mask, unsigned long long upper) const` | count_and_sum_less_or を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = (x \| mask) としたとき、y < upper である要素の元の x の総和を返す。 |
| `sum_type sum_less_equal_or(int l, int r, unsigned long long mask, unsigned long long upper) const` | count_and_sum_less_equal_or を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = (x \| mask) としたとき、y <= upper である要素の元の x の総和を返す。 |
| `sum_type sum_equal_or(int l, int r, unsigned long long mask, unsigned long long value) const` | count_and_sum_less_or を 2 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = (x \| mask) としたとき、y == value である要素の元の x の総和を返す。 |
| `sum_type sum_range_or( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper ) const` | count_and_sum_less_or を 2 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = (x \| mask) としたとき、lower <= y < upper である要素の元の x の総和を返す。 |
| `sum_type sum_greater_or(int l, int r, unsigned long long mask, unsigned long long lower) const` | count_and_sum_greater_or を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = (x \| mask) としたとき、y > lower である要素の元の x の総和を返す。 |
| `sum_type sum_greater_equal_or(int l, int r, unsigned long long mask, unsigned long long lower) const` | count_and_sum_greater_equal_or を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = (x \| mask) としたとき、y >= lower である要素の元の x の総和を返す。 |
| `unsigned long long access_or(int i, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | i 番目の要素に (value \| mask) を施した結果を返す。 |
| `int rank_or(unsigned long long value, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | prefix 区間 [0, r) の各要素 x に対して y = (x \| mask) としたとき、y == value である要素数を返す。 |
| `int rank_or(unsigned long long value, int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して y = (x \| mask) としたとき、y == value である要素数を返す。 |
| `int count_or(unsigned long long value, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 配列全体の各要素 x に対して y = (x \| mask) としたとき、y == value である要素数を返す。 |
| `int count_or(unsigned long long value, int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して y = (x \| mask) としたとき、y == value である要素数を返す。 |
| `bool contains_or(unsigned long long value, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 配列全体の各要素 x に対して y = (x \| mask) としたとき、y == value である要素が存在するかを返す。 |
| `std::vector<unsigned long long> values_or(unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 配列全体において現れる (value \| mask) の異なる値を昇順で返す。 |
| `int distinct_size_or(unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 配列全体において現れる (value \| mask) の異なる値の個数を返す。 |
| `int index_of_or(unsigned long long value, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | values_or(mask) が返す昇順 distinct 配列の中で、value が現れる位置を返す。存在しなければ -1 を返す。 |
| `int select_or(unsigned long long value, int kth, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 0-indexed で kth 番目に現れる、z = (x \| mask) が z == value となる位置を返す。存在しなければ -1 を返す。 |
| `unsigned long long kth_smallest_or(int l, int r, int k, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して z = (x \| mask) とした値を昇順に見たときの k 番目の値を返す。k は 0-indexed。 |
| `unsigned long long kth_largest_or(int l, int r, int k, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して z = (x \| mask) とした値を降順に見たときの k 番目の値を返す。k は 0-indexed。 |
| `unsigned long long quantile_or(int l, int r, int k, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | kth_smallest_or の別名。 |
| `std::optional<unsigned long long> min_value_or(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) の最小値を返す。空区間なら std::nullopt。 |
| `std::optional<unsigned long long> max_value_or(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) の最大値を返す。空区間なら std::nullopt。 |
| `std::optional<unsigned long long> median_lower_or(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して z = (x \| mask) とした値の下側中央値を返す。空区間なら std::nullopt を返す。 |
| `std::optional<unsigned long long> median_upper_or(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して z = (x \| mask) とした値の上側中央値を返す。空区間なら std::nullopt を返す。 |
| `std::optional<unsigned long long> next_value_ge_or(int l, int r, unsigned long long lower, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) の中で、lower 以上の最小値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> next_value_gt_or(int l, int r, unsigned long long lower, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) の中で、lower より大きい最小値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> prev_value_lt_or(int l, int r, unsigned long long upper, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) の中で、upper 未満の最大値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> prev_value_le_or(int l, int r, unsigned long long upper, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) の中で、upper 以下の最大値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> next_value_or(int l, int r, unsigned long long lower, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | next_value_ge_or の別名。 |
| `std::optional<unsigned long long> prev_value_or(int l, int r, unsigned long long upper, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | prev_value_lt_or の別名。 |
| `std::vector<std::pair<unsigned long long, int>> list_frequencies_or(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して z = (x \| mask) とした値について、(値, 個数) を昇順に列挙する。変換後に同じ値になったものは集約する。 |
| `std::vector<std::pair<unsigned long long, int>> list_frequencies_or(int l, int r, unsigned long long lower, unsigned long long upper, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して z = (x \| mask) とした値について、(値, 個数) を昇順に列挙する。変換後に同じ値になったものは集約する。 |
| `std::vector<unsigned long long> distinct_values_or(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して z = (x \| mask) としたときに現れる異なる値を昇順に列挙する。 |
| `std::vector<std::pair<unsigned long long, int>> top_k_frequent_or(int l, int r, int k, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) に対し、頻度上位 k 個の (値, 出現回数) を返す。 |
| `std::optional<std::pair<unsigned long long, int>> mode_or(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) の最頻値を (値, 出現回数) で返す。空区間なら std::nullopt。 |
| `std::vector<std::tuple<unsigned long long, int, int>> intersect_or(int l1, int r1, int l2, int r2, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 2 区間 [l1, r1), [l2, r2) について、共通して現れる (value \| mask) の値を列挙する。 |

### AND ラッパー


| 関数   | 時間計算量 | 説明 |
|---|---|---|
| `int count_less_and(int l, int r, unsigned long long mask, unsigned long long upper) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value & mask) < upper の個数を返す。 |
| `int range_freq_and( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper ) const` | count_less_and を 2 回呼ぶので、それに準ずる | 区間 [l, r) において lower <= (value & mask) < upper の個数を返す。 |
| `int count_less_equal_and(int l, int r, unsigned long long mask, unsigned long long upper) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value & mask) <= upper の個数を返す。 |
| `int count_equal_and(int l, int r, unsigned long long mask, unsigned long long value) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して y = (x & mask) としたとき、y == value である要素数を返す。 |
| `int count_greater_and(int l, int r, unsigned long long mask, unsigned long long lower) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value & mask) > lower の個数を返す。 |
| `int count_greater_equal_and(int l, int r, unsigned long long mask, unsigned long long lower) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value & mask) >= lower の個数を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_and( int l, int r, unsigned long long mask, unsigned long long upper ) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して y = (x & mask) としたとき、y < upper である要素について、(個数, 元の x の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_equal_and( int l, int r, unsigned long long mask, unsigned long long upper ) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して y = (x & mask) としたとき、y <= upper である要素について、(個数, 元の x の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_range_and( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper ) const` | count_and_sum_less_and を 2 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = (x & mask) としたとき、lower <= y < upper である要素について、(個数, 元の x の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_and( int l, int r, unsigned long long mask, unsigned long long lower ) const` | count_and_sum_less_equal_and を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = (x & mask) としたとき、y > lower である要素について、(個数, 元の x の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_equal_and( int l, int r, unsigned long long mask, unsigned long long lower ) const` | count_and_sum_less_and を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = (x & mask) としたとき、y >= lower である要素について、(個数, 元の x の総和) を返す。 |
| `sum_type sum_less_and(int l, int r, unsigned long long mask, unsigned long long upper) const` | count_and_sum_less_and を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = (x & mask) としたとき、y < upper である要素の元の x の総和を返す。 |
| `sum_type sum_less_equal_and(int l, int r, unsigned long long mask, unsigned long long upper) const` | count_and_sum_less_equal_and を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = (x & mask) としたとき、y <= upper である要素の元の x の総和を返す。 |
| `sum_type sum_equal_and(int l, int r, unsigned long long mask, unsigned long long value) const` | count_and_sum_less_and を 2 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = (x & mask) としたとき、y == value である要素の元の x の総和を返す。 |
| `sum_type sum_range_and( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper ) const` | count_and_sum_less_and を 2 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = (x & mask) としたとき、lower <= y < upper である要素の元の x の総和を返す。 |
| `sum_type sum_greater_and(int l, int r, unsigned long long mask, unsigned long long lower) const` | count_and_sum_greater_and を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = (x & mask) としたとき、y > lower である要素の元の x の総和を返す。 |
| `sum_type sum_greater_equal_and(int l, int r, unsigned long long mask, unsigned long long lower) const` | count_and_sum_greater_equal_and を 1 回呼ぶので、それに準ずる | 区間 [l, r) の各要素 x に対して y = (x & mask) としたとき、y >= lower である要素の元の x の総和を返す。 |
| `unsigned long long access_and(int i, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | i 番目の要素に (value & mask) を施した結果を返す。 |
| `int rank_and(unsigned long long value, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | prefix 区間 [0, r) の各要素 x に対して y = (x & mask) としたとき、y == value である要素数を返す。 |
| `int rank_and(unsigned long long value, int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して y = (x & mask) としたとき、y == value である要素数を返す。 |
| `int count_and(unsigned long long value, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 配列全体の各要素 x に対して y = (x & mask) としたとき、y == value である要素数を返す。 |
| `int count_and(unsigned long long value, int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して y = (x & mask) としたとき、y == value である要素数を返す。 |
| `bool contains_and(unsigned long long value, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 配列全体の各要素 x に対して y = (x & mask) としたとき、y == value である要素が存在するかを返す。 |
| `std::vector<unsigned long long> values_and(unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 配列全体において現れる (value & mask) の異なる値を昇順で返す。 |
| `int distinct_size_and(unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 配列全体において現れる (value & mask) の異なる値の個数を返す。 |
| `int index_of_and(unsigned long long value, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | values_and(mask) が返す昇順 distinct 配列の中で、value が現れる位置を返す。存在しなければ -1 を返す。 |
| `int select_and(unsigned long long value, int kth, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 0-indexed で kth 番目に現れる、z = (x & mask) が z == value となる位置を返す。存在しなければ -1 を返す。 |
| `unsigned long long kth_smallest_and(int l, int r, int k, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して z = (x & mask) とした値を昇順に見たときの k 番目の値を返す。k は 0-indexed。 |
| `unsigned long long kth_largest_and(int l, int r, int k, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して z = (x & mask) とした値を降順に見たときの k 番目の値を返す。k は 0-indexed。 |
| `unsigned long long quantile_and(int l, int r, int k, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | kth_smallest_and の別名。 |
| `std::optional<unsigned long long> min_value_and(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) の最小値を返す。空区間なら std::nullopt。 |
| `std::optional<unsigned long long> max_value_and(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) の最大値を返す。空区間なら std::nullopt。 |
| `std::optional<unsigned long long> median_lower_and(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して z = (x & mask) とした値の下側中央値を返す。空区間なら std::nullopt を返す。 |
| `std::optional<unsigned long long> median_upper_and(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して z = (x & mask) とした値の上側中央値を返す。空区間なら std::nullopt を返す。 |
| `std::optional<unsigned long long> next_value_ge_and(int l, int r, unsigned long long lower, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) の中で、lower 以上の最小値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> next_value_gt_and(int l, int r, unsigned long long lower, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) の中で、lower より大きい最小値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> prev_value_lt_and(int l, int r, unsigned long long upper, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) の中で、upper 未満の最大値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> prev_value_le_and(int l, int r, unsigned long long upper, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) の中で、upper 以下の最大値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> next_value_and(int l, int r, unsigned long long lower, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | next_value_ge_and の別名。 |
| `std::optional<unsigned long long> prev_value_and(int l, int r, unsigned long long upper, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | prev_value_lt_and の別名。 |
| `std::vector<std::pair<unsigned long long, int>> list_frequencies_and(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して z = (x & mask) とした値について、(値, 個数) を昇順に列挙する。変換後に同じ値になったものは集約する。 |
| `std::vector<std::pair<unsigned long long, int>> list_frequencies_and(int l, int r, unsigned long long lower, unsigned long long upper, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して z = (x & mask) とした値について、(値, 個数) を昇順に列挙する。変換後に同じ値になったものは集約する。 |
| `std::vector<unsigned long long> distinct_values_and(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の各要素 x に対して z = (x & mask) としたときに現れる異なる値を昇順に列挙する。 |
| `std::vector<std::pair<unsigned long long, int>> top_k_frequent_and(int l, int r, int k, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) に対し、頻度上位 k 個の (値, 出現回数) を返す。 |
| `std::optional<std::pair<unsigned long long, int>> mode_and(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) の最頻値を (値, 出現回数) で返す。空区間なら std::nullopt。 |
| `std::vector<std::tuple<unsigned long long, int, int>> intersect_and(int l1, int r1, int l2, int r2, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 2 区間 [l1, r1), [l2, r2) について、共通して現れる (value & mask) の値を列挙する。 |

## 型エイリアス

### `WaveletMatrix<T>`

| 型エイリアス | 説明 |
|---|---|
| `value_type = T` | 値型。 |
| `sum_type = long long` | 総和系クエリで使う型。 |
| `unsigned_value_type = typename WaveletMatrixBitwiseUnsignedHelper<T>::type` | bitwise 系クエリで使う符号なし整数型。 |

# Weighted Wavelet Matrix
[ライブラリのリンク](https://github.com/harurunrunrun/GPT-generated-Library/blob/main/library/WeightedWaveletMatrixMonoid.hpp)

## 概要

- この章では、`WeightedWaveletMatrix<T, Weight>` と、その上に載る 2 次元クエリ・モノイド系の構造をまとめています。
- `WaveletMatrix<T>` と `WeightedWaveletMatrix<T, Weight>` は**別クラス**です。ここでは両者を混ぜず、重み付き・モノイド付きの側だけを扱います。
- `WeightedWaveletMatrix<T, Weight>` は、個数だけでなく**重み和**も扱います。
- `CommutativeMonoidWaveletMatrix<T, Monoid>` と `RectangleMonoid<X, Y, Monoid>` は、和だけでなく**可換モノイドの積**を扱います。

## 記号

- `N`: 要素数または点数
- `σ`: 異なる値の個数
- `m`: 出力される異なる値の個数

## WeightedWaveletMatrix

### 構築・メタ情報

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `WeightedWaveletMatrix() = default` | O(1) | 空の `WeightedWaveletMatrix` を作る。 |
| `WeightedWaveletMatrix(const std::vector<T>& data, const std::vector<Weight>& weights)` | O(N log N + N log σ) | 値列 `data` と重み列 `weights` から構築する。 |
| `void build(const std::vector<T>& data, const std::vector<Weight>& weights)` | O(N log N + N log σ) | 値列 `data` と重み列 `weights` から再構築する。 |
| `int size() const` | いずれも O(1) | 元配列の長さ `N` を返す。 |
| `bool empty() const` | O(1) | 空かどうかを返す。 |
| `int distinct_size() const` | O(1) | 異なる値の個数 `σ` を返す。 |
| `int bit_size() const` | O(1) | 内部で使うビット長を返す。 |
| `const std::vector<T>& values() const` | O(1) | 座標圧縮後の異なる値を昇順で返す。 |

### 参照・出現回数・位置

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `const T& access(int i) const` | いずれも O(1) | 元配列の `i` 番目の値を返す。 |
| `const T& operator[](int i) const` | O(1) | `access(i)` の別名。 |
| `const Weight& weight(int i) const` | O(1) | 元配列の `i` 番目の重みを返す。 |
| `bool contains(const T& value) const` | いずれも O(log σ) | 配列全体に `value` が 1 回以上現れるかを返す。 |
| `int index_of(const T& value) const` | O(log σ) | `value` の座標圧縮後の ID を返す。存在しなければ `-1`。 |
| `int rank(const T& value, int r) const` | いずれも O(log σ) | 区間 `[0, r)` に `value` が何回現れるかを返す。 |
| `int rank(const T& value, int l, int r) const` | O(log σ) | 区間 `[l, r)` に `value` が何回現れるかを返す。 |
| `int count(const T& value) const` | O(log σ) | 配列全体に `value` が何回現れるかを返す。 |
| `int count(const T& value, int l, int r) const` | O(log σ) | 区間 `[l, r)` に `value` が何回現れるかを返す。 |
| `sum_type weight_sum(const T& value, int l, int r) const` | O(log σ) | 区間 `[l, r)` にある `value` の重み和を返す。 |
| `sum_type weight_sum(const T& value) const` | O(log σ) | 配列全体にある `value` の重み和を返す。 |
| `int select(const T& value, int kth) const` | O(log σ log N) | `value` が 0-indexed で `kth` 回目に現れる位置を返す。存在しなければ `-1`。 |

### 順序統計・代表値

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `const T& kth_smallest(int l, int r, int k) const` | いずれも O(log σ) | 区間 `[l, r)` の値を昇順に見たときの `k` 番目の値を返す。 |
| `const T& kth_largest(int l, int r, int k) const` | O(log σ) | 区間 `[l, r)` の値を降順に見たときの `k` 番目の値を返す。 |
| `const T& quantile(int l, int r, int k) const` | O(log σ) | `kth_smallest(l, r, k)` の別名。 |
| `std::optional<T> min_value(int l, int r) const` | いずれも O(log σ) | 区間 `[l, r)` の最小値を返す。空区間なら `std::nullopt`。 |
| `std::optional<T> max_value(int l, int r) const` | O(log σ) | 区間 `[l, r)` の最大値を返す。空区間なら `std::nullopt`。 |
| `std::optional<T> median_lower(int l, int r) const` | O(log σ) | 区間 `[l, r)` の下側中央値を返す。空区間なら `std::nullopt`。 |
| `std::optional<T> median_upper(int l, int r) const` | O(log σ) | 区間 `[l, r)` の上側中央値を返す。空区間なら `std::nullopt`。 |

### 値域条件での個数・重み和

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `int range_freq(int l, int r, const T& upper) const` | いずれも O(log σ) | 区間 `[l, r)` で `x < upper` となる要素数を返す。 |
| `int range_freq(int l, int r, const T& lower, const T& upper) const` | O(log σ) | 区間 `[l, r)` で `lower <= x < upper` となる要素数を返す。 |
| `int count_less(int l, int r, const T& upper) const` | O(log σ) | 区間 `[l, r)` で `x < upper` となる要素数を返す。 |
| `int count_less_equal(int l, int r, const T& upper) const` | O(log σ) | 区間 `[l, r)` で `x <= upper` となる要素数を返す。 |
| `int count_greater(int l, int r, const T& lower) const` | O(log σ) | 区間 `[l, r)` で `x > lower` となる要素数を返す。 |
| `int count_greater_equal(int l, int r, const T& lower) const` | O(log σ) | 区間 `[l, r)` で `x >= lower` となる要素数を返す。 |
| `std::pair<int, sum_type> count_and_weight_less(int l, int r, const T& upper) const` | いずれも O(log σ) | 区間 `[l, r)` で `x < upper` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_less_equal(int l, int r, const T& upper) const` | O(log σ) | 区間 `[l, r)` で `x <= upper` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_greater(int l, int r, const T& lower) const` | O(log σ) | 区間 `[l, r)` で `x > lower` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_greater_equal(int l, int r, const T& lower) const` | O(log σ) | 区間 `[l, r)` で `x >= lower` となる要素の `(個数, 重み和)` を返す。 |
| `sum_type weight_sum_all(int l, int r) const` | weight_sum_all は O(1) / それ以外は O(log σ) | 区間 `[l, r)` の重み和を返す。 |
| `sum_type weight_sum_less(int l, int r, const T& upper) const` | O(log σ) | 区間 `[l, r)` で `x < upper` となる要素の重み和を返す。 |
| `sum_type weight_sum_less_equal(int l, int r, const T& upper) const` | O(log σ) | 区間 `[l, r)` で `x <= upper` となる要素の重み和を返す。 |
| `sum_type weight_sum_range(int l, int r, const T& lower, const T& upper) const` | O(log σ) | 区間 `[l, r)` で `lower <= x < upper` となる要素の重み和を返す。 |
| `sum_type weight_sum_greater(int l, int r, const T& lower) const` | O(log σ) | 区間 `[l, r)` で `x > lower` となる要素の重み和を返す。 |
| `sum_type weight_sum_greater_equal(int l, int r, const T& lower) const` | O(log σ) | 区間 `[l, r)` で `x >= lower` となる要素の重み和を返す。 |

### k 個の重み和

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `sum_type weight_sum_k_smallest(int l, int r, int k) const` | O(log σ) | 区間 `[l, r)` の値を昇順に見たとき、小さい方から `k` 個の重み和を返す。 |
| `sum_type weight_sum_k_largest(int l, int r, int k) const` | O(log σ) | 区間 `[l, r)` の値を降順に見たとき、大きい方から `k` 個の重み和を返す。 |

### 前後探索

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `std::optional<T> next_value_ge(int l, int r, const T& lower) const` | いずれも O(log σ) | 区間 `[l, r)` で `x >= lower` となる最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<T> next_value_gt(int l, int r, const T& lower) const` | O(log σ) | 区間 `[l, r)` で `x > lower` となる最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<T> prev_value_lt(int l, int r, const T& upper) const` | O(log σ) | 区間 `[l, r)` で `x < upper` となる最大値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<T> prev_value_le(int l, int r, const T& upper) const` | O(log σ) | 区間 `[l, r)` で `x <= upper` となる最大値を返す。存在しなければ `std::nullopt`。 |

### 列挙

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `std::vector<std::pair<T, int>> list_frequencies(int l, int r) const` | おおよそ O((m + 1) log σ) / m は出力される異なる値の個数 | 区間 `[l, r)` に現れる異なる値を `(値, 個数)` で昇順に返す。 |
| `std::vector<std::pair<T, int>> list_frequencies(int l, int r, const T& lower, const T& upper) const` | おおよそ O((m + 1) log σ) / m は出力される異なる値の個数 | 区間 `[l, r)` で `lower <= x < upper` となる値だけを `(値, 個数)` で昇順に返す。 |
| `std::vector<std::pair<T, sum_type>> list_weight_sums(int l, int r) const` | おおよそ O((m + 1) log σ) / m は出力される異なる値の個数 | 区間 `[l, r)` に現れる異なる値を `(値, 重み和)` で昇順に返す。 |
| `std::vector<std::pair<T, sum_type>> list_weight_sums(int l, int r, const T& lower, const T& upper) const` | おおよそ O((m + 1) log σ) / m は出力される異なる値の個数 | 区間 `[l, r)` で `lower <= x < upper` となる値だけを `(値, 重み和)` で昇順に返す。 |
| `std::vector<std::tuple<T, int, sum_type>> list_frequencies_and_weight_sums(int l, int r) const` | おおよそ O((m + 1) log σ) / m は出力される異なる値の個数 | 区間 `[l, r)` の `y' = (y & mask)` に現れる異なる値を `(値, 個数)` で昇順に返す。 |
| `std::vector<std::tuple<T, int, sum_type>> list_frequencies_and_weight_sums( int l, int r, const T& lower, const T& upper ) const` | おおよそ O((m + 1) log σ) / m は出力される異なる値の個数 | 区間 `[l, r)` の `y' = (y & mask)` に現れる異なる値を `(値, 個数)` で昇順に返す。 |
| `std::vector<T> distinct_values(int l, int r) const` | おおよそ O((m + 1) log σ) / m は出力される異なる値の個数 | 区間 `[l, r)` に現れる異なる値を昇順で返す。 |

### XOR 変換後のクエリ

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `int count_xor_less(int l, int r, const T& mask, const T& upper) const` | 以下の変換値クエリは、内部で区間内の異なる値を列挙して集約するため、 / おおよそ O((m + 1) log σ + m log m)。 / m は区間 [l, r) に現れる異なる元の値の個数。 | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' < upper` となる要素数を返す。 |
| `int count_xor_less_equal(int l, int r, const T& mask, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' <= upper` となる要素数を返す。 |
| `int count_xor_range(int l, int r, const T& mask, const T& lower, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `lower <= y' < upper` となる要素数を返す。 |
| `int count_xor_greater(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' > lower` となる要素数を返す。 |
| `int count_xor_greater_equal(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' >= lower` となる要素数を返す。 |
| `int count_xor_equal(int l, int r, const T& mask, const T& value) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' == value` となる要素数を返す。 |
| `sum_type weight_sum_xor_less(int l, int r, const T& mask, const T& upper) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' < upper` となる要素の重み和を返す。 |
| `sum_type weight_sum_xor_less_equal(int l, int r, const T& mask, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' <= upper` となる要素の重み和を返す。 |
| `sum_type weight_sum_xor_range(int l, int r, const T& mask, const T& lower, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `lower <= y' < upper` となる要素の重み和を返す。 |
| `sum_type weight_sum_xor_greater(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' > lower` となる要素の重み和を返す。 |
| `sum_type weight_sum_xor_greater_equal(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' >= lower` となる要素の重み和を返す。 |
| `sum_type weight_sum_xor_equal(int l, int r, const T& mask, const T& value) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' == value` となる要素の重み和を返す。 |
| `std::pair<int, sum_type> count_and_weight_xor_less(int l, int r, const T& mask, const T& upper) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' < upper` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_xor_less_equal(int l, int r, const T& mask, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' <= upper` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_xor_range(int l, int r, const T& mask, const T& lower, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `lower <= y' < upper` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_xor_greater(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' > lower` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_xor_greater_equal(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' >= lower` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_xor_equal(int l, int r, const T& mask, const T& value) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' == value` となる要素の `(個数, 重み和)` を返す。 |
| `T kth_smallest_xor(int l, int r, const T& mask, int k) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y ^ mask)` を昇順に見たときの `k` 番目の値を返す。 |
| `T kth_largest_xor(int l, int r, const T& mask, int k) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y ^ mask)` を降順に見たときの `k` 番目の値を返す。 |
| `T quantile_xor(int l, int r, const T& mask, int k) const` | おおよそ O((m + 1) log σ + m log m) | 上の `kth_smallest_*` の別名。 |
| `std::optional<T> min_xor(int l, int r, const T& mask) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y ^ mask)` の最小値を返す。空区間なら `std::nullopt`。 |
| `std::optional<T> max_xor(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y ^ mask)` の最大値を返す。空区間なら `std::nullopt`。 |
| `std::optional<T> median_lower_xor(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y ^ mask)` の下側中央値を返す。空区間なら `std::nullopt`。 |
| `std::optional<T> median_upper_xor(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y ^ mask)` の上側中央値を返す。空区間なら `std::nullopt`。 |
| `std::optional<T> next_xor_value_ge(int l, int r, const T& mask, const T& lower) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y ^ mask)` のうち `y' >= lower` となる最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<T> next_xor_value_gt(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y ^ mask)` のうち `y' > lower` となる最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<T> prev_xor_value_lt(int l, int r, const T& mask, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y ^ mask)` のうち `y' < upper` となる最大値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<T> prev_xor_value_le(int l, int r, const T& mask, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y ^ mask)` のうち `y' <= upper` となる最大値を返す。存在しなければ `std::nullopt`。 |
| `std::vector<std::pair<T, int>> list_frequencies_xor(int l, int r, const T& mask) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y ^ mask)` に現れる異なる値を `(値, 個数)` で昇順に返す。 |
| `std::vector<std::pair<T, sum_type>> list_weight_sums_xor(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y ^ mask)` に現れる異なる値を `(値, 重み和)` で昇順に返す。 |
| `std::vector<std::tuple<T, int, sum_type>> list_frequencies_and_weight_sums_xor(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y ^ mask)` に現れる異なる値を `(値, 個数)` で昇順に返す。 |
| `std::vector<T> distinct_values_xor(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y ^ mask)` に現れる異なる値を昇順で返す。 |

### OR 変換後のクエリ

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `int count_or_less(int l, int r, const T& mask, const T& upper) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' < upper` となる要素数を返す。 |
| `int count_or_less_equal(int l, int r, const T& mask, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' <= upper` となる要素数を返す。 |
| `int count_or_range(int l, int r, const T& mask, const T& lower, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `lower <= y' < upper` となる要素数を返す。 |
| `int count_or_greater(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' > lower` となる要素数を返す。 |
| `int count_or_greater_equal(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' >= lower` となる要素数を返す。 |
| `int count_or_equal(int l, int r, const T& mask, const T& value) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' == value` となる要素数を返す。 |
| `sum_type weight_sum_or_less(int l, int r, const T& mask, const T& upper) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' < upper` となる要素の重み和を返す。 |
| `sum_type weight_sum_or_less_equal(int l, int r, const T& mask, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' <= upper` となる要素の重み和を返す。 |
| `sum_type weight_sum_or_range(int l, int r, const T& mask, const T& lower, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `lower <= y' < upper` となる要素の重み和を返す。 |
| `sum_type weight_sum_or_greater(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' > lower` となる要素の重み和を返す。 |
| `sum_type weight_sum_or_greater_equal(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' >= lower` となる要素の重み和を返す。 |
| `sum_type weight_sum_or_equal(int l, int r, const T& mask, const T& value) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' == value` となる要素の重み和を返す。 |
| `std::pair<int, sum_type> count_and_weight_or_less(int l, int r, const T& mask, const T& upper) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' < upper` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_or_less_equal(int l, int r, const T& mask, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' <= upper` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_or_range(int l, int r, const T& mask, const T& lower, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `lower <= y' < upper` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_or_greater(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' > lower` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_or_greater_equal(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' >= lower` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_or_equal(int l, int r, const T& mask, const T& value) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' == value` となる要素の `(個数, 重み和)` を返す。 |
| `T kth_smallest_or(int l, int r, const T& mask, int k) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y \| mask)` を昇順に見たときの `k` 番目の値を返す。 |
| `T kth_largest_or(int l, int r, const T& mask, int k) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y \| mask)` を降順に見たときの `k` 番目の値を返す。 |
| `T quantile_or(int l, int r, const T& mask, int k) const` | おおよそ O((m + 1) log σ + m log m) | 上の `kth_smallest_*` の別名。 |
| `std::optional<T> min_or(int l, int r, const T& mask) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y \| mask)` の最小値を返す。空区間なら `std::nullopt`。 |
| `std::optional<T> max_or(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y \| mask)` の最大値を返す。空区間なら `std::nullopt`。 |
| `std::optional<T> median_lower_or(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y \| mask)` の下側中央値を返す。空区間なら `std::nullopt`。 |
| `std::optional<T> median_upper_or(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y \| mask)` の上側中央値を返す。空区間なら `std::nullopt`。 |
| `std::optional<T> next_or_value_ge(int l, int r, const T& mask, const T& lower) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y \| mask)` のうち `y' >= lower` となる最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<T> next_or_value_gt(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y \| mask)` のうち `y' > lower` となる最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<T> prev_or_value_lt(int l, int r, const T& mask, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y \| mask)` のうち `y' < upper` となる最大値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<T> prev_or_value_le(int l, int r, const T& mask, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y \| mask)` のうち `y' <= upper` となる最大値を返す。存在しなければ `std::nullopt`。 |
| `std::vector<std::pair<T, int>> list_frequencies_or(int l, int r, const T& mask) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y \| mask)` に現れる異なる値を `(値, 個数)` で昇順に返す。 |
| `std::vector<std::pair<T, sum_type>> list_weight_sums_or(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y \| mask)` に現れる異なる値を `(値, 重み和)` で昇順に返す。 |
| `std::vector<std::tuple<T, int, sum_type>> list_frequencies_and_weight_sums_or(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y \| mask)` に現れる異なる値を `(値, 個数)` で昇順に返す。 |
| `std::vector<T> distinct_values_or(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y \| mask)` に現れる異なる値を昇順で返す。 |

### AND 変換後のクエリ

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `std::vector<std::tuple<T, int, sum_type>> list_frequencies_and_weight_sums(int l, int r) const` | おおよそ O((m + 1) log σ) / m は出力される異なる値の個数 | 区間 `[l, r)` の `y' = (y & mask)` に現れる異なる値を `(値, 個数)` で昇順に返す。 |
| `std::vector<std::tuple<T, int, sum_type>> list_frequencies_and_weight_sums( int l, int r, const T& lower, const T& upper ) const` | おおよそ O((m + 1) log σ) / m は出力される異なる値の個数 | 区間 `[l, r)` の `y' = (y & mask)` に現れる異なる値を `(値, 個数)` で昇順に返す。 |
| `std::pair<int, sum_type> count_and_weight_xor_less(int l, int r, const T& mask, const T& upper) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' < upper` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_xor_less_equal(int l, int r, const T& mask, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' <= upper` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_xor_range(int l, int r, const T& mask, const T& lower, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `lower <= y' < upper` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_xor_greater(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' > lower` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_xor_greater_equal(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' >= lower` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_xor_equal(int l, int r, const T& mask, const T& value) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y ^ mask)` としたとき `y' == value` となる要素の `(個数, 重み和)` を返す。 |
| `std::vector<std::tuple<T, int, sum_type>> list_frequencies_and_weight_sums_xor(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y ^ mask)` に現れる異なる値を `(値, 個数)` で昇順に返す。 |
| `std::pair<int, sum_type> count_and_weight_or_less(int l, int r, const T& mask, const T& upper) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' < upper` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_or_less_equal(int l, int r, const T& mask, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' <= upper` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_or_range(int l, int r, const T& mask, const T& lower, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `lower <= y' < upper` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_or_greater(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' > lower` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_or_greater_equal(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' >= lower` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_or_equal(int l, int r, const T& mask, const T& value) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y \| mask)` としたとき `y' == value` となる要素の `(個数, 重み和)` を返す。 |
| `std::vector<std::tuple<T, int, sum_type>> list_frequencies_and_weight_sums_or(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y \| mask)` に現れる異なる値を `(値, 個数)` で昇順に返す。 |
| `int count_and_less(int l, int r, const T& mask, const T& upper) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y & mask)` としたとき `y' < upper` となる要素数を返す。 |
| `int count_and_less_equal(int l, int r, const T& mask, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y & mask)` としたとき `y' <= upper` となる要素数を返す。 |
| `int count_and_range(int l, int r, const T& mask, const T& lower, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y & mask)` としたとき `lower <= y' < upper` となる要素数を返す。 |
| `int count_and_greater(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y & mask)` としたとき `y' > lower` となる要素数を返す。 |
| `int count_and_greater_equal(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y & mask)` としたとき `y' >= lower` となる要素数を返す。 |
| `int count_and_equal(int l, int r, const T& mask, const T& value) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y & mask)` としたとき `y' == value` となる要素数を返す。 |
| `sum_type weight_sum_and_less(int l, int r, const T& mask, const T& upper) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y & mask)` としたとき `y' < upper` となる要素の重み和を返す。 |
| `sum_type weight_sum_and_less_equal(int l, int r, const T& mask, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y & mask)` としたとき `y' <= upper` となる要素の重み和を返す。 |
| `sum_type weight_sum_and_range(int l, int r, const T& mask, const T& lower, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y & mask)` としたとき `lower <= y' < upper` となる要素の重み和を返す。 |
| `sum_type weight_sum_and_greater(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y & mask)` としたとき `y' > lower` となる要素の重み和を返す。 |
| `sum_type weight_sum_and_greater_equal(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y & mask)` としたとき `y' >= lower` となる要素の重み和を返す。 |
| `sum_type weight_sum_and_equal(int l, int r, const T& mask, const T& value) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y & mask)` としたとき `y' == value` となる要素の重み和を返す。 |
| `std::pair<int, sum_type> count_and_weight_and_less(int l, int r, const T& mask, const T& upper) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y & mask)` としたとき `y' < upper` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_and_less_equal(int l, int r, const T& mask, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y & mask)` としたとき `y' <= upper` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_and_range(int l, int r, const T& mask, const T& lower, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y & mask)` としたとき `lower <= y' < upper` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_and_greater(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y & mask)` としたとき `y' > lower` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_and_greater_equal(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y & mask)` としたとき `y' >= lower` となる要素の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_weight_and_equal(int l, int r, const T& mask, const T& value) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` で `y' = (y & mask)` としたとき `y' == value` となる要素の `(個数, 重み和)` を返す。 |
| `T kth_smallest_and(int l, int r, const T& mask, int k) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y & mask)` を昇順に見たときの `k` 番目の値を返す。 |
| `T kth_largest_and(int l, int r, const T& mask, int k) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y & mask)` を降順に見たときの `k` 番目の値を返す。 |
| `T quantile_and(int l, int r, const T& mask, int k) const` | おおよそ O((m + 1) log σ + m log m) | 上の `kth_smallest_*` の別名。 |
| `std::optional<T> min_and(int l, int r, const T& mask) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y & mask)` の最小値を返す。空区間なら `std::nullopt`。 |
| `std::optional<T> max_and(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y & mask)` の最大値を返す。空区間なら `std::nullopt`。 |
| `std::optional<T> median_lower_and(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y & mask)` の下側中央値を返す。空区間なら `std::nullopt`。 |
| `std::optional<T> median_upper_and(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y & mask)` の上側中央値を返す。空区間なら `std::nullopt`。 |
| `std::optional<T> next_and_value_ge(int l, int r, const T& mask, const T& lower) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y & mask)` のうち `y' >= lower` となる最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<T> next_and_value_gt(int l, int r, const T& mask, const T& lower) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y & mask)` のうち `y' > lower` となる最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<T> prev_and_value_lt(int l, int r, const T& mask, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y & mask)` のうち `y' < upper` となる最大値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<T> prev_and_value_le(int l, int r, const T& mask, const T& upper) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y & mask)` のうち `y' <= upper` となる最大値を返す。存在しなければ `std::nullopt`。 |
| `std::vector<std::pair<T, int>> list_frequencies_and(int l, int r, const T& mask) const` | いずれもおおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y & mask)` に現れる異なる値を `(値, 個数)` で昇順に返す。 |
| `std::vector<std::pair<T, sum_type>> list_weight_sums_and(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y & mask)` に現れる異なる値を `(値, 重み和)` で昇順に返す。 |
| `std::vector<std::tuple<T, int, sum_type>> list_frequencies_and_weight_sums_and(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y & mask)` に現れる異なる値を `(値, 個数)` で昇順に返す。 |
| `std::vector<T> distinct_values_and(int l, int r, const T& mask) const` | おおよそ O((m + 1) log σ + m log m) | 区間 `[l, r)` の `y' = (y & mask)` に現れる異なる値を昇順で返す。 |

## RectangleSum

### 構築・メタ情報

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `RectangleSum() = default` | O(1) | 空の `RectangleSum` を作る。 |
| `explicit RectangleSum(const std::vector<std::tuple<X, Y, Weight>>& points)` | O(N log N + N log σ) | 点列 `(x, y, w)` から構築する。 |
| `RectangleSum( const std::vector<X>& xs, const std::vector<Y>& ys, const std::vector<Weight>& weights )` | O(N log N + N log σ) | 座標列 `xs`, `ys` と重み列 `weights` から構築する。 |
| `void build(const std::vector<std::tuple<X, Y, Weight>>& points)` | O(N log N + N log σ) | 点列 `(x, y, w)` から再構築する。 |
| `void build( const std::vector<X>& xs, const std::vector<Y>& ys, const std::vector<Weight>& weights )` | O(N log N + N log σ) | 座標列 `xs`, `ys` と重み列 `weights` から再構築する。 |
| `int size() const` | いずれも O(1) | 内部に保持している点数を返す。 |
| `bool empty() const` | O(1) | 空かどうかを返す。 |
| `const std::vector<X>& x_values() const` | いずれも O(1) | `x` で昇順に並べた `x` 座標列を返す。 |
| `const std::vector<Y>& y_values_sorted_by_x() const` | O(1) | `x` で昇順に並べた順の `y` 座標列を返す。 |
| `const std::vector<Weight>& weights_sorted_by_x() const` | O(1) | `x` で昇順に並べた順の重み列を返す。 |
| `const wm_type& wavelet_matrix() const` | O(1) | 内部で使っている `WeightedWaveletMatrix` を返す。 |

### 基本の矩形クエリ

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `int count(const X& xl, const X& xr, const Y& yl, const Y& yr) const` | O(log N + log σ) | 矩形 `[xl, xr) × [yl, yr)` に入る点数を返す。 |
| `sum_type sum(const X& xl, const X& xr, const Y& yl, const Y& yr) const` | O(log N + log σ) | 矩形 `[xl, xr) × [yl, yr)` に入る点の重み和を返す。 |
| `int count_all(const X& xl, const X& xr) const` | count_all は O(log N) / sum_all は O(log N) | `x` が `[xl, xr)` に入る点数を返す。 |
| `sum_type sum_all(const X& xl, const X& xr) const` | O(log N + log σ) | `x` が `[xl, xr)` に入る点の重み和を返す。 |

### x 範囲内での y 条件による個数・重み和

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `int count_less_y(const X& xl, const X& xr, const Y& upper) const` | いずれも O(log N + log σ) | `x` が `[xl, xr)` に入り、`y < upper` となる点数を返す。 |
| `int count_less_equal_y(const X& xl, const X& xr, const Y& upper) const` | O(log N + log σ) | `x` が `[xl, xr)` に入り、`y <= upper` となる点数を返す。 |
| `int count_greater_y(const X& xl, const X& xr, const Y& lower) const` | O(log N + log σ) | `x` が `[xl, xr)` に入り、`y > lower` となる点数を返す。 |
| `int count_greater_equal_y(const X& xl, const X& xr, const Y& lower) const` | O(log N + log σ) | `x` が `[xl, xr)` に入り、`y >= lower` となる点数を返す。 |
| `sum_type sum_less_y(const X& xl, const X& xr, const Y& upper) const` | いずれも O(log N + log σ) | `x` が `[xl, xr)` に入り、`y < upper` となる点の重み和を返す。 |
| `sum_type sum_less_equal_y(const X& xl, const X& xr, const Y& upper) const` | O(log N + log σ) | `x` が `[xl, xr)` に入り、`y <= upper` となる点の重み和を返す。 |
| `sum_type sum_greater_y(const X& xl, const X& xr, const Y& lower) const` | O(log N + log σ) | `x` が `[xl, xr)` に入り、`y > lower` となる点の重み和を返す。 |
| `sum_type sum_greater_equal_y(const X& xl, const X& xr, const Y& lower) const` | O(log N + log σ) | `x` が `[xl, xr)` に入り、`y >= lower` となる点の重み和を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_y(const X& xl, const X& xr, const Y& upper) const` | いずれも O(log N + log σ) | `x` が `[xl, xr)` に入り、`y < upper` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_equal_y(const X& xl, const X& xr, const Y& upper) const` | O(log N + log σ) | `x` が `[xl, xr)` に入り、`y <= upper` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_y(const X& xl, const X& xr, const Y& lower) const` | O(log N + log σ) | `x` が `[xl, xr)` に入り、`y > lower` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_equal_y(const X& xl, const X& xr, const Y& lower) const` | O(log N + log σ) | `x` が `[xl, xr)` に入り、`y >= lower` となる点の `(個数, 重み和)` を返す。 |

### x 範囲内での y の順序統計・探索

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `const Y& kth_smallest_y(const X& xl, const X& xr, int k) const` | O(log N + log σ) | `x` が `[xl, xr)` に入る点の `y` を昇順に見たときの `k` 番目を返す。 |
| `const Y& kth_largest_y(const X& xl, const X& xr, int k) const` | O(log N + log σ) | `x` が `[xl, xr)` に入る点の `y` を降順に見たときの `k` 番目を返す。 |
| `sum_type sum_k_smallest_y(const X& xl, const X& xr, int k) const` | O(log N + log σ) | `x` が `[xl, xr)` に入る点の `y` を昇順に見たとき、小さい方から `k` 個の重み和を返す。 |
| `sum_type sum_k_largest_y(const X& xl, const X& xr, int k) const` | O(log N + log σ) | `x` が `[xl, xr)` に入る点の `y` を降順に見たとき、大きい方から `k` 個の重み和を返す。 |
| `std::optional<Y> min_y(const X& xl, const X& xr) const` | いずれも O(log N + log σ) | `x` が `[xl, xr)` に入る点の `y` の最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> max_y(const X& xl, const X& xr) const` | O(log N + log σ) | `x` が `[xl, xr)` に入る点の `y` の最大値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> next_y_ge(const X& xl, const X& xr, const Y& lower) const` | O(log N + log σ) | `x` が `[xl, xr)` に入る点の `y` のうち `y >= lower` となる最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> next_y_gt(const X& xl, const X& xr, const Y& lower) const` | O(log N + log σ) | `x` が `[xl, xr)` に入る点の `y` のうち `y > lower` となる最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> prev_y_lt(const X& xl, const X& xr, const Y& upper) const` | O(log N + log σ) | `x` が `[xl, xr)` に入る点の `y` のうち `y < upper` となる最大値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> prev_y_le(const X& xl, const X& xr, const Y& upper) const` | O(log N + log σ) | `x` が `[xl, xr)` に入る点の `y` のうち `y <= upper` となる最大値を返す。存在しなければ `std::nullopt`。 |

### x 範囲内での y の列挙

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `std::vector<std::pair<Y, int>> list_frequencies(const X& xl, const X& xr) const` | おおよそ O(log N + (m + 1) log σ) / m は出力される異なる y の個数 | `x` が `[xl, xr)` に入る点の `y` を `(y, 個数)` で昇順に返す。 |
| `std::vector<std::pair<Y, int>> list_frequencies( const X& xl, const X& xr, const Y& yl, const Y& yr ) const` | O(log N) + おおよそ O((m + 1) log σ) | `x` が `[xl, xr)` に入り、`yl <= y < yr` となる値だけを `(y, 個数)` で昇順に返す。 |
| `std::vector<std::pair<Y, sum_type>> list_weight_sums(const X& xl, const X& xr) const` | O(log N) + おおよそ O((m + 1) log σ) | `x` が `[xl, xr)` に入る点の `y` を `(y, 重み和)` で昇順に返す。 |
| `std::vector<std::pair<Y, sum_type>> list_weight_sums( const X& xl, const X& xr, const Y& yl, const Y& yr ) const` | O(log N) + おおよそ O((m + 1) log σ) | `x` が `[xl, xr)` に入り、`yl <= y < yr` となる値だけを `(y, 重み和)` で昇順に返す。 |
| `std::vector<std::tuple<Y, int, sum_type>> list_frequencies_and_weight_sums( const X& xl, const X& xr ) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y & mask)` を `(値, 個数)` で昇順に返す。 |
| `std::vector<std::tuple<Y, int, sum_type>> list_frequencies_and_weight_sums( const X& xl, const X& xr, const Y& yl, const Y& yr ) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y & mask)` を `(値, 個数)` で昇順に返す。 |

### x 範囲内での XOR 変換 y クエリ

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `int count_xor_less_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | 以下の変換値クエリは、x 範囲の特定に O(log N)、 / さらに WeightedWaveletMatrix 側でおおよそ O((m + 1) log σ + m log m)。 | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' < upper` となる点数を返す。 |
| `int count_xor_less_equal_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' <= upper` となる点数を返す。 |
| `int count_xor_range_y(const X& xl, const X& xr, const Y& mask, const Y& lower, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `lower <= y' < upper` となる点数を返す。 |
| `int count_xor_greater_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' > lower` となる点数を返す。 |
| `int count_xor_greater_equal_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' >= lower` となる点数を返す。 |
| `int count_xor_equal_y(const X& xl, const X& xr, const Y& mask, const Y& value) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' == value` となる点数を返す。 |
| `sum_type sum_xor_less_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' < upper` となる点の重み和を返す。 |
| `sum_type sum_xor_less_equal_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' <= upper` となる点の重み和を返す。 |
| `sum_type sum_xor_range_y(const X& xl, const X& xr, const Y& mask, const Y& lower, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `lower <= y' < upper` となる点の重み和を返す。 |
| `sum_type sum_xor_greater_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' > lower` となる点の重み和を返す。 |
| `sum_type sum_xor_greater_equal_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' >= lower` となる点の重み和を返す。 |
| `sum_type sum_xor_equal_y(const X& xl, const X& xr, const Y& mask, const Y& value) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' == value` となる点の重み和を返す。 |
| `std::pair<int, sum_type> count_and_sum_xor_less_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' < upper` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_xor_less_equal_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' <= upper` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_xor_range_y(const X& xl, const X& xr, const Y& mask, const Y& lower, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `lower <= y' < upper` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_xor_greater_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' > lower` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_xor_greater_equal_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' >= lower` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_xor_equal_y(const X& xl, const X& xr, const Y& mask, const Y& value) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' == value` となる点の `(個数, 重み和)` を返す。 |
| `Y kth_smallest_xor_y(const X& xl, const X& xr, const Y& mask, int k) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y ^ mask)` を昇順に見たときの `k` 番目を返す。 |
| `Y kth_largest_xor_y(const X& xl, const X& xr, const Y& mask, int k) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y ^ mask)` を降順に見たときの `k` 番目を返す。 |
| `Y quantile_xor_y(const X& xl, const X& xr, const Y& mask, int k) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | 上の `kth_smallest_*` の別名。 |
| `std::optional<Y> min_xor_y(const X& xl, const X& xr, const Y& mask) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y ^ mask)` の最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> max_xor_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y ^ mask)` の最大値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> median_lower_xor_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y ^ mask)` の下側中央値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> median_upper_xor_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y ^ mask)` の上側中央値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> next_xor_value_ge_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y ^ mask)` のうち `y' >= lower` となる最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> next_xor_value_gt_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y ^ mask)` のうち `y' > lower` となる最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> prev_xor_value_lt_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y ^ mask)` のうち `y' < upper` となる最大値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> prev_xor_value_le_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y ^ mask)` のうち `y' <= upper` となる最大値を返す。存在しなければ `std::nullopt`。 |
| `std::vector<std::pair<Y, int>> list_frequencies_xor_y(const X& xl, const X& xr, const Y& mask) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y ^ mask)` を `(値, 個数)` で昇順に返す。 |
| `std::vector<std::pair<Y, sum_type>> list_weight_sums_xor_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y ^ mask)` を `(値, 重み和)` で昇順に返す。 |
| `std::vector<std::tuple<Y, int, sum_type>> list_frequencies_and_weight_sums_xor_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y ^ mask)` を `(値, 個数)` で昇順に返す。 |
| `std::vector<Y> distinct_values_xor_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y ^ mask)` に現れる異なる値を昇順で返す。 |

### x 範囲内での OR 変換 y クエリ

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `int count_or_less_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' < upper` となる点数を返す。 |
| `int count_or_less_equal_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' <= upper` となる点数を返す。 |
| `int count_or_range_y(const X& xl, const X& xr, const Y& mask, const Y& lower, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `lower <= y' < upper` となる点数を返す。 |
| `int count_or_greater_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' > lower` となる点数を返す。 |
| `int count_or_greater_equal_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' >= lower` となる点数を返す。 |
| `int count_or_equal_y(const X& xl, const X& xr, const Y& mask, const Y& value) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' == value` となる点数を返す。 |
| `sum_type sum_or_less_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' < upper` となる点の重み和を返す。 |
| `sum_type sum_or_less_equal_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' <= upper` となる点の重み和を返す。 |
| `sum_type sum_or_range_y(const X& xl, const X& xr, const Y& mask, const Y& lower, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `lower <= y' < upper` となる点の重み和を返す。 |
| `sum_type sum_or_greater_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' > lower` となる点の重み和を返す。 |
| `sum_type sum_or_greater_equal_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' >= lower` となる点の重み和を返す。 |
| `sum_type sum_or_equal_y(const X& xl, const X& xr, const Y& mask, const Y& value) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' == value` となる点の重み和を返す。 |
| `std::pair<int, sum_type> count_and_sum_or_less_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' < upper` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_or_less_equal_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' <= upper` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_or_range_y(const X& xl, const X& xr, const Y& mask, const Y& lower, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `lower <= y' < upper` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_or_greater_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' > lower` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_or_greater_equal_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' >= lower` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_or_equal_y(const X& xl, const X& xr, const Y& mask, const Y& value) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' == value` となる点の `(個数, 重み和)` を返す。 |
| `Y kth_smallest_or_y(const X& xl, const X& xr, const Y& mask, int k) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y \| mask)` を昇順に見たときの `k` 番目を返す。 |
| `Y kth_largest_or_y(const X& xl, const X& xr, const Y& mask, int k) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y \| mask)` を降順に見たときの `k` 番目を返す。 |
| `Y quantile_or_y(const X& xl, const X& xr, const Y& mask, int k) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | 上の `kth_smallest_*` の別名。 |
| `std::optional<Y> min_or_y(const X& xl, const X& xr, const Y& mask) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y \| mask)` の最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> max_or_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y \| mask)` の最大値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> median_lower_or_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y \| mask)` の下側中央値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> median_upper_or_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y \| mask)` の上側中央値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> next_or_value_ge_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y \| mask)` のうち `y' >= lower` となる最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> next_or_value_gt_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y \| mask)` のうち `y' > lower` となる最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> prev_or_value_lt_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y \| mask)` のうち `y' < upper` となる最大値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> prev_or_value_le_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y \| mask)` のうち `y' <= upper` となる最大値を返す。存在しなければ `std::nullopt`。 |
| `std::vector<std::pair<Y, int>> list_frequencies_or_y(const X& xl, const X& xr, const Y& mask) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y \| mask)` を `(値, 個数)` で昇順に返す。 |
| `std::vector<std::pair<Y, sum_type>> list_weight_sums_or_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y \| mask)` を `(値, 重み和)` で昇順に返す。 |
| `std::vector<std::tuple<Y, int, sum_type>> list_frequencies_and_weight_sums_or_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y \| mask)` を `(値, 個数)` で昇順に返す。 |
| `std::vector<Y> distinct_values_or_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y \| mask)` に現れる異なる値を昇順で返す。 |

### x 範囲内での AND 変換 y クエリ

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `std::vector<std::tuple<Y, int, sum_type>> list_frequencies_and_weight_sums( const X& xl, const X& xr ) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y & mask)` を `(値, 個数)` で昇順に返す。 |
| `std::vector<std::tuple<Y, int, sum_type>> list_frequencies_and_weight_sums( const X& xl, const X& xr, const Y& yl, const Y& yr ) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y & mask)` を `(値, 個数)` で昇順に返す。 |
| `std::pair<int, sum_type> count_and_sum_xor_less_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' < upper` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_xor_less_equal_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' <= upper` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_xor_range_y(const X& xl, const X& xr, const Y& mask, const Y& lower, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `lower <= y' < upper` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_xor_greater_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' > lower` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_xor_greater_equal_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' >= lower` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_xor_equal_y(const X& xl, const X& xr, const Y& mask, const Y& value) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y ^ mask)` について `y' == value` となる点の `(個数, 重み和)` を返す。 |
| `std::vector<std::tuple<Y, int, sum_type>> list_frequencies_and_weight_sums_xor_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y ^ mask)` を `(値, 個数)` で昇順に返す。 |
| `std::pair<int, sum_type> count_and_sum_or_less_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' < upper` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_or_less_equal_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' <= upper` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_or_range_y(const X& xl, const X& xr, const Y& mask, const Y& lower, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `lower <= y' < upper` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_or_greater_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' > lower` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_or_greater_equal_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' >= lower` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_or_equal_y(const X& xl, const X& xr, const Y& mask, const Y& value) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y \| mask)` について `y' == value` となる点の `(個数, 重み和)` を返す。 |
| `std::vector<std::tuple<Y, int, sum_type>> list_frequencies_and_weight_sums_or_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y \| mask)` を `(値, 個数)` で昇順に返す。 |
| `int count_and_less_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y & mask)` について `y' < upper` となる点数を返す。 |
| `int count_and_less_equal_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y & mask)` について `y' <= upper` となる点数を返す。 |
| `int count_and_range_y(const X& xl, const X& xr, const Y& mask, const Y& lower, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y & mask)` について `lower <= y' < upper` となる点数を返す。 |
| `int count_and_greater_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y & mask)` について `y' > lower` となる点数を返す。 |
| `int count_and_greater_equal_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y & mask)` について `y' >= lower` となる点数を返す。 |
| `int count_and_equal_y(const X& xl, const X& xr, const Y& mask, const Y& value) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y & mask)` について `y' == value` となる点数を返す。 |
| `sum_type sum_and_less_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y & mask)` について `y' < upper` となる点の重み和を返す。 |
| `sum_type sum_and_less_equal_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y & mask)` について `y' <= upper` となる点の重み和を返す。 |
| `sum_type sum_and_range_y(const X& xl, const X& xr, const Y& mask, const Y& lower, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y & mask)` について `lower <= y' < upper` となる点の重み和を返す。 |
| `sum_type sum_and_greater_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y & mask)` について `y' > lower` となる点の重み和を返す。 |
| `sum_type sum_and_greater_equal_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y & mask)` について `y' >= lower` となる点の重み和を返す。 |
| `sum_type sum_and_equal_y(const X& xl, const X& xr, const Y& mask, const Y& value) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y & mask)` について `y' == value` となる点の重み和を返す。 |
| `std::pair<int, sum_type> count_and_sum_and_less_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y & mask)` について `y' < upper` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_and_less_equal_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y & mask)` について `y' <= upper` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_and_range_y(const X& xl, const X& xr, const Y& mask, const Y& lower, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y & mask)` について `lower <= y' < upper` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_and_greater_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y & mask)` について `y' > lower` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_and_greater_equal_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y & mask)` について `y' >= lower` となる点の `(個数, 重み和)` を返す。 |
| `std::pair<int, sum_type> count_and_sum_and_equal_y(const X& xl, const X& xr, const Y& mask, const Y& value) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入り、変換後の `y' = (y & mask)` について `y' == value` となる点の `(個数, 重み和)` を返す。 |
| `Y kth_smallest_and_y(const X& xl, const X& xr, const Y& mask, int k) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y & mask)` を昇順に見たときの `k` 番目を返す。 |
| `Y kth_largest_and_y(const X& xl, const X& xr, const Y& mask, int k) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y & mask)` を降順に見たときの `k` 番目を返す。 |
| `Y quantile_and_y(const X& xl, const X& xr, const Y& mask, int k) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | 上の `kth_smallest_*` の別名。 |
| `std::optional<Y> min_and_y(const X& xl, const X& xr, const Y& mask) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y & mask)` の最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> max_and_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y & mask)` の最大値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> median_lower_and_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y & mask)` の下側中央値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> median_upper_and_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y & mask)` の上側中央値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> next_and_value_ge_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y & mask)` のうち `y' >= lower` となる最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> next_and_value_gt_y(const X& xl, const X& xr, const Y& mask, const Y& lower) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y & mask)` のうち `y' > lower` となる最小値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> prev_and_value_lt_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y & mask)` のうち `y' < upper` となる最大値を返す。存在しなければ `std::nullopt`。 |
| `std::optional<Y> prev_and_value_le_y(const X& xl, const X& xr, const Y& mask, const Y& upper) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y & mask)` のうち `y' <= upper` となる最大値を返す。存在しなければ `std::nullopt`。 |
| `std::vector<std::pair<Y, int>> list_frequencies_and_y(const X& xl, const X& xr, const Y& mask) const` | いずれもおおよそ O(log N + (m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y & mask)` を `(値, 個数)` で昇順に返す。 |
| `std::vector<std::pair<Y, sum_type>> list_weight_sums_and_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y & mask)` を `(値, 重み和)` で昇順に返す。 |
| `std::vector<std::tuple<Y, int, sum_type>> list_frequencies_and_weight_sums_and_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y & mask)` を `(値, 個数)` で昇順に返す。 |
| `std::vector<Y> distinct_values_and_y(const X& xl, const X& xr, const Y& mask) const` | O(log N) + おおよそ O((m + 1) log σ + m log m) | `x` が `[xl, xr)` に入る点の変換後の `y' = (y & mask)` に現れる異なる値を昇順で返す。 |

## StaticRangeMonoid

静的配列上の range fold を行うセグメント木。

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `StaticRangeMonoid() = default` | O(1) | 空のオブジェクトを作る。 |
| `explicit StaticRangeMonoid(const std::vector<value_type>& a)` | O(N) | 引数からオブジェクトを構築する。 |
| `void build(const std::vector<value_type>& a)` | O(N) | 配列から再構築する。 |
| `value_type prod(int l, int r) const` | O(log N) | 区間 [l, r) のモノイド積を返す。 |

## CommutativeMonoidWaveletMatrix

値列 data と、それぞれに対応するモノイド値 vals を持つ Wavelet Matrix。

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `CommutativeMonoidWaveletMatrix() = default` | O(1) | 空のオブジェクトを作る。 |
| `CommutativeMonoidWaveletMatrix(const std::vector<T>& data, const std::vector<prod_type>& vals)` | O(N log N + N log σ) | 引数からオブジェクトを構築する。 |
| `void build(const std::vector<T>& data, const std::vector<prod_type>& vals)` | O(N log N + N log σ) | 値列とモノイド値列から再構築する。 |
| `int size() const` | O(1) | 元の配列長を返す。 |
| `bool empty() const` | O(1) | 空かどうかを返す。 |
| `int distinct_size() const` | O(1) | 異なる値の個数を返す。 |
| `int bit_size() const` | O(1) | 内部で使うビット長を返す。 |
| `const std::vector<T>& values() const` | O(1) | 座標圧縮後の異なる値一覧を返す。 |
| `int range_freq(int l, int r, const T& lower, const T& upper) const` | O(log σ) | 区間 [l, r) の要素 x について、lower <= x < upper である要素数を返す。 |
| `prod_type prod_range(int l, int r, const T& lower, const T& upper) const` | O(log N * log σ) | 区間 [l, r) で値域 [lower, upper) に入る要素のモノイド積を返す。 |
| `prod_type prod_all(int l, int r) const` | O(log N * log σ) | 区間 [l, r) 全体のモノイド積を返す。 |

## RectangleMonoid

点集合 (x, y, value) に対して、長方形 [xl, xr) x [yl, yr) 内の 個数とモノイド積を返す静的データ構造。

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `RectangleMonoid() = default` | O(1) | 空のオブジェクトを作る。 |
| `explicit RectangleMonoid(const std::vector<std::tuple<X, Y, prod_type>>& points)` | O(N log N + N log σ) | 引数からオブジェクトを構築する。 |
| `RectangleMonoid( const std::vector<X>& xs, const std::vector<Y>& ys, const std::vector<prod_type>& vals )` | O(N log N + N log σ) | 引数からオブジェクトを構築する。 |
| `void build(const std::vector<std::tuple<X, Y, prod_type>>& points)` | O(N log N + N log σ) | 点集合から再構築する。 |
| `void build( const std::vector<X>& xs, const std::vector<Y>& ys, const std::vector<prod_type>& vals )` | O(N log N + N log σ) | 点集合から再構築する。 |
| `int size() const` | O(1) | 保持している点数を返す。 |
| `bool empty() const` | O(1) | 空かどうかを返す。 |
| `int count(const X& xl, const X& xr, const Y& yl, const Y& yr) const` | O(log N + log σ) | 矩形 [xl, xr) × [yl, yr) に入る点数を返す。 |
| `prod_type prod(const X& xl, const X& xr, const Y& yl, const Y& yr) const` | O(log N + log N * log σ) | 矩形 [xl, xr) × [yl, yr) のモノイド積を返す。 |
| `prod_type prod_all(const X& xl, const X& xr) const` | O(log N + log N * log σ) | x 範囲 [xl, xr) 全体のモノイド積を返す。 |

## rect_monoid::Sum

通常の加法モノイド。

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `static constexpr value_type id()` | O(1) | モノイドの単位元を返す。 |
| `static constexpr value_type op(const value_type& a, const value_type& b)` | O(1) | モノイド演算を適用する。 |

## rect_monoid::Product

通常の乗法モノイド。

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `static constexpr value_type id()` | O(1) | モノイドの単位元を返す。 |
| `static constexpr value_type op(const value_type& a, const value_type& b)` | O(1) | モノイド演算を適用する。 |

## rect_monoid::Min

最小値モノイド。

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `static constexpr value_type id()` | O(1) | モノイドの単位元を返す。 |
| `static constexpr value_type op(const value_type& a, const value_type& b)` | O(1) | モノイド演算を適用する。 |

## rect_monoid::Max

最大値モノイド。

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `static constexpr value_type id()` | O(1) | モノイドの単位元を返す。 |
| `static constexpr value_type op(const value_type& a, const value_type& b)` | O(1) | モノイド演算を適用する。 |

## rect_monoid::Gcd

最大公約数モノイド。

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `static constexpr value_type id()` | O(1) | モノイドの単位元を返す。 |
| `static constexpr value_type op(const value_type& a, const value_type& b)` | O(1) | モノイド演算を適用する。 |

## rect_monoid::Lcm

最小公倍数モノイド。

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `static constexpr value_type id()` | O(1) | モノイドの単位元を返す。 |
| `static constexpr value_type op(const value_type& a, const value_type& b)` | O(1) | モノイド演算を適用する。 |

## rect_monoid::Xor

xor モノイド。

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `static constexpr value_type id()` | O(1) | モノイドの単位元を返す。 |
| `static constexpr value_type op(const value_type& a, const value_type& b)` | O(1) | モノイド演算を適用する。 |

## rect_monoid::BitAnd

bitwise AND モノイド。

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `static constexpr value_type id()` | O(1) | モノイドの単位元を返す。 |
| `static constexpr value_type op(const value_type& a, const value_type& b)` | O(1) | モノイド演算を適用する。 |

## rect_monoid::BitOr

bitwise OR モノイド。

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `static constexpr value_type id()` | O(1) | モノイドの単位元を返す。 |
| `static constexpr value_type op(const value_type& a, const value_type& b)` | O(1) | モノイド演算を適用する。 |

## rect_monoid::AddMod998244353

mod 998244353 での加法モノイド。

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `static constexpr value_type id()` | O(1) | モノイドの単位元を返す。 |
| `static constexpr value_type normalize(value_type x)` | O(1) | 値を法 998244353 の範囲へ正規化する。 |
| `static constexpr value_type op(value_type a, value_type b)` | O(1) | モノイド演算を適用する。 |

## rect_monoid::MulMod998244353

mod 998244353 での乗法モノイド。

| 関数 | 時間計算量 | 説明 |
|---|---|---|
| `static constexpr value_type id()` | O(1) | モノイドの単位元を返す。 |
| `static constexpr value_type normalize(value_type x)` | O(1) | 値を法 998244353 の範囲へ正規化する。 |
| `static constexpr value_type op(value_type a, value_type b)` | O(1) | モノイド演算を適用する。 |

## 型エイリアス

### `WeightedWaveletMatrix`

| 型エイリアス | 説明 |
|---|---|
| `value_type = T` | 値型。 |
| `weight_type = Weight` | 重み型。 |
| `sum_type = std::common_type_t<Weight, long long>` | 重み和で使う型。 |

### `RectangleSum`

| 型エイリアス | 説明 |
|---|---|
| `x_type = X` | x 座標の型。 |
| `y_type = Y` | y 座標の型。 |
| `weight_type = Weight` | 重み型。 |
| `wm_type = WeightedWaveletMatrix<Y, Weight>` | 内部で使う WeightedWaveletMatrix の型。 |
| `sum_type = typename wm_type::sum_type` | 重み和で使う型。 |

### `rect_monoid::*`

| 型エイリアス | 説明 |
|---|---|
| `rect_monoid::Sum<T>::value_type = T` | 加法モノイドの値型。 |
| `rect_monoid::Product<T>::value_type = T` | 乗法モノイドの値型。 |
| `rect_monoid::Min<T>::value_type = T` | 最小値モノイドの値型。 |
| `rect_monoid::Max<T>::value_type = T` | 最大値モノイドの値型。 |
| `rect_monoid::Gcd<T>::value_type = T` | gcd モノイドの値型。 |
| `rect_monoid::Lcm<T>::value_type = T` | lcm モノイドの値型。 |
| `rect_monoid::Xor<T>::value_type = T` | xor モノイドの値型。 |
| `rect_monoid::BitAnd<T>::value_type = T` | bitwise and モノイドの値型。 |
| `rect_monoid::BitOr<T>::value_type = T` | bitwise or モノイドの値型。 |
| `rect_monoid::AddMod998244353::value_type = std::int64_t` | mod 998244353 の加法モノイドで使う型。 |
| `rect_monoid::MulMod998244353::value_type = std::int64_t` | mod 998244353 の乗法モノイドで使う型。 |

### `StaticRangeMonoid`

| 型エイリアス | 説明 |
|---|---|
| `value_type = typename Monoid::value_type` | 内部セグメント木で扱う値型。 |

### `CommutativeMonoidWaveletMatrix`

| 型エイリアス | 説明 |
|---|---|
| `value_type = T` | 値型。 |
| `prod_type = typename Monoid::value_type` | モノイド積で使う型。 |

### `RectangleMonoid`

| 型エイリアス | 説明 |
|---|---|
| `x_type = X` | x 座標の型。 |
| `y_type = Y` | y 座標の型。 |
| `prod_type = typename Monoid::value_type` | 矩形クエリで返すモノイド積の型。 |
| `wm_type = CommutativeMonoidWaveletMatrix<Y, Monoid>` | 内部で使う Wavelet Matrix の型。 |

### 便利エイリアス

| 型エイリアス | 説明 |
|---|---|
| `MonoidWaveletMatrix` | 任意の可換モノイドで、値域 `[lower, upper)` 上のモノイド積を求める用途のエイリアス。 |
| `SumWaveletMatrix` | 値域ごとの和を求める用途の Wavelet Matrix 用エイリアス。 |
| `ProductWaveletMatrix` | 値域ごとの積を求める用途の Wavelet Matrix 用エイリアス。 |
| `MinWaveletMatrix` | 値域ごとの最小値を求める用途の Wavelet Matrix 用エイリアス。 |
| `MaxWaveletMatrix` | 値域ごとの最大値を求める用途の Wavelet Matrix 用エイリアス。 |
| `GcdWaveletMatrix` | 値域ごとの gcd を求める用途の Wavelet Matrix 用エイリアス。 |
| `LcmWaveletMatrix` | 値域ごとの lcm を求める用途の Wavelet Matrix 用エイリアス。 |
| `XorWaveletMatrix` | 値域ごとの xor を求める用途の Wavelet Matrix 用エイリアス。 |
| `BitAndWaveletMatrix` | 値域ごとの bitwise and を求める用途の Wavelet Matrix 用エイリアス。 |
| `BitOrWaveletMatrix` | 値域ごとの bitwise or を求める用途の Wavelet Matrix 用エイリアス。 |
| `AddMod998244353WaveletMatrix` | 値域ごとの和を mod 998244353 で求める用途の Wavelet Matrix 用エイリアス。 |
| `MulMod998244353WaveletMatrix` | 値域ごとの積を mod 998244353 で求める用途の Wavelet Matrix 用エイリアス。 |
| `RectangleAdd` | 矩形内の和を求める用途のエイリアス。 |
| `RectangleProduct` | 矩形内の積を求める用途のエイリアス。 |
| `RectangleMin` | 矩形内の最小値を求める用途のエイリアス。 |
| `RectangleMax` | 矩形内の最大値を求める用途のエイリアス。 |
| `RectangleGcd` | 矩形内の gcd を求める用途のエイリアス。 |
| `RectangleLcm` | 矩形内の lcm を求める用途のエイリアス。 |
| `RectangleXor` | 矩形内の xor を求める用途のエイリアス。 |
| `RectangleBitAnd` | 矩形内の bitwise and を求める用途のエイリアス。 |
| `RectangleBitOr` | 矩形内の bitwise or を求める用途のエイリアス。 |
| `RectangleAddMod998244353` | 矩形内の和を mod 998244353 で求める用途のエイリアス。 |
| `RectangleMulMod998244353` | 矩形内の積を mod 998244353 で求める用途のエイリアス。 |
