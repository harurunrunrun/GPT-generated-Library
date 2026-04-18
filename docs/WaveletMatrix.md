# Wavelet Matrix
[ライブラリのリンク](https://github.com/harurunrunrun/GPT-generated-Library/blob/main/library/WaveletMatrixMultifunction.hpp)

## 概要

- `BitVector`: 15 個の公開関数
- `WaveletMatrix<T>`: 239 個の公開関数
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
| `int range_freq(int l, int r, const T& upper) const` | O(log σ) | 区間 [l, r) において upper 未満の値の個数を返す。 つまり x < upper を数える。 |
| `int range_freq(int l, int r, const T& lower, const T& upper) const` | O(log σ) | 区間 [l, r) において lower <= x < upper を満たす個数を返す。 |
| `int count_less(int l, int r, const T& upper) const` | O(log σ) | 区間 [l, r) において x < upper の個数。 |
| `int count_less_equal(int l, int r, const T& upper) const` | O(log σ) | 区間 [l, r) において x <= upper の個数。 |
| `int count_greater(int l, int r, const T& lower) const` | O(log σ) | 区間 [l, r) において x > lower の個数。 |
| `int count_greater_equal(int l, int r, const T& lower) const` | O(log σ) | 区間 [l, r) において x >= lower の個数。 |

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
| `std::pair<int, sum_type> count_and_sum_less(int l, int r, const T& upper) const` | O(log σ) | 区間 [l, r) において x < upper を満たす要素の (個数, 総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_equal(int l, int r, const T& upper) const` | O(log σ) | 区間 [l, r) において x <= upper を満たす要素の (個数, 総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater(int l, int r, const T& lower) const` | O(log σ) | 区間 [l, r) において x > lower を満たす要素の (個数, 総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_equal(int l, int r, const T& lower) const` | O(log σ) | 区間 [l, r) において x >= lower を満たす要素の (個数, 総和) を返す。 |
| `sum_type sum_all(int l, int r) const` | O(1) | 区間 [l, r) の総和を返す。 |
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
| `int count_less_bitwise(int l, int r, unsigned long long mask, unsigned long long upper, BitwiseOperation op) const` | XOR のとき O(W) / OR / AND のとき O(M) / W は unsigned(T) のビット幅、M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において、各要素 value に bitwise 演算を施した 結果 op(value, mask) が upper 未満である個数を返す。 |
| `int range_freq_bitwise( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper, BitwiseOperation op ) const` | count_less_bitwise を 2 回呼ぶので、それに準ずる | 区間 [l, r) において、各要素 value に bitwise 演算を施した 結果 op(value, mask) が lower <= x < upper を満たす個数を返す。 |
| `int count_less_equal_bitwise(int l, int r, unsigned long long mask, unsigned long long upper, BitwiseOperation op) const` | count_less_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) において、各要素 value に bitwise 演算を施した 結果 op(value, mask) が upper 以下である個数を返す。 |
| `int count_equal_bitwise(int l, int r, unsigned long long mask, unsigned long long value, BitwiseOperation op) const` | count_less_bitwise を 2 回呼ぶので、それに準ずる | 区間 [l, r) において、各要素 value に bitwise 演算を施した 結果 op(value, mask) がちょうど value に等しい個数を返す。 |
| `int count_greater_bitwise(int l, int r, unsigned long long mask, unsigned long long lower, BitwiseOperation op) const` | count_less_equal_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) において、各要素 value に bitwise 演算を施した 結果 op(value, mask) が lower より大きい個数を返す。 |
| `int count_greater_equal_bitwise(int l, int r, unsigned long long mask, unsigned long long lower, BitwiseOperation op) const` | count_less_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) において、各要素 value に bitwise 演算を施した 結果 op(value, mask) が lower 以上である個数を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_bitwise( int l, int r, unsigned long long mask, unsigned long long upper, BitwiseOperation op ) const` | XOR のとき O(W) / OR / AND のとき O(M) / W は unsigned(T) のビット幅、M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において、各要素 value に bitwise 演算を施した 結果 op(value, mask) が upper 未満である要素について、 (個数, 元の value の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_equal_bitwise( int l, int r, unsigned long long mask, unsigned long long upper, BitwiseOperation op ) const` | count_and_sum_less_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) において、各要素 value に bitwise 演算を施した 結果 op(value, mask) が upper 以下である要素について、 (個数, 元の value の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_range_bitwise( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper, BitwiseOperation op ) const` | count_and_sum_less_bitwise を 2 回呼ぶので、それに準ずる | 区間 [l, r) において、各要素 value に bitwise 演算を施した 結果 op(value, mask) が lower <= x < upper を満たす要素について、 (個数, 元の value の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_bitwise( int l, int r, unsigned long long mask, unsigned long long lower, BitwiseOperation op ) const` | count_and_sum_less_equal_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) において、各要素 value に bitwise 演算を施した 結果 op(value, mask) が lower より大きい要素について、 (個数, 元の value の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_equal_bitwise( int l, int r, unsigned long long mask, unsigned long long lower, BitwiseOperation op ) const` | count_and_sum_less_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) において、各要素 value に bitwise 演算を施した 結果 op(value, mask) が lower 以上の要素について、 (個数, 元の value の総和) を返す。 |
| `sum_type sum_less_bitwise( int l, int r, unsigned long long mask, unsigned long long upper, BitwiseOperation op ) const` | count_and_sum_less_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) において、各要素 value に bitwise 演算を施した 結果 op(value, mask) が upper 未満である要素の、 元の value の総和を返す。 |
| `sum_type sum_less_equal_bitwise( int l, int r, unsigned long long mask, unsigned long long upper, BitwiseOperation op ) const` | count_and_sum_less_equal_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) において、各要素 value に bitwise 演算を施した 結果 op(value, mask) が upper 以下である要素の、 元の value の総和を返す。 |
| `sum_type sum_equal_bitwise( int l, int r, unsigned long long mask, unsigned long long value, BitwiseOperation op ) const` | count_and_sum_less_bitwise を 2 回呼ぶので、それに準ずる | 区間 [l, r) において、各要素 value に bitwise 演算を施した 結果 op(value, mask) がちょうど value に等しい要素の、 元の value の総和を返す。 |
| `sum_type sum_range_bitwise( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper, BitwiseOperation op ) const` | count_and_sum_less_bitwise を 2 回呼ぶので、それに準ずる | 区間 [l, r) において、各要素 value に bitwise 演算を施した 結果 op(value, mask) が lower <= x < upper を満たす要素の、 元の value の総和を返す。 |
| `sum_type sum_greater_bitwise( int l, int r, unsigned long long mask, unsigned long long lower, BitwiseOperation op ) const` | count_and_sum_greater_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) において、各要素 value に bitwise 演算を施した 結果 op(value, mask) が lower より大きい要素の、 元の value の総和を返す。 |
| `sum_type sum_greater_equal_bitwise( int l, int r, unsigned long long mask, unsigned long long lower, BitwiseOperation op ) const` | count_and_sum_greater_equal_bitwise を 1 回呼ぶので、それに準ずる | 区間 [l, r) において、各要素 value に bitwise 演算を施した 結果 op(value, mask) が lower 以上の要素の、 元の value の総和を返す。 |

### Bitwise 共通 API（変換後の値に対する検索・列挙）


| 関数   | 時間計算量 | 説明 |
|---|---|---|
| `unsigned long long access_bitwise(int i, unsigned long long mask, BitwiseOperation op) const` | O(1) | i 番目の要素 value に bitwise 演算を施した結果 op(value, mask) を返す。 |
| `int rank_bitwise(unsigned long long value, int r, unsigned long long mask, BitwiseOperation op) const` | count_equal_bitwise に準ずる | prefix 区間 [0, r) において、 op(element, mask) == value を満たす個数を返す。 |
| `int rank_bitwise(unsigned long long value, int l, int r, unsigned long long mask, BitwiseOperation op) const` | count_equal_bitwise に準ずる | 区間 [l, r) において、 op(element, mask) == value を満たす個数を返す。 |
| `int count_bitwise(unsigned long long value, unsigned long long mask, BitwiseOperation op) const` | count_equal_bitwise に準ずる | 配列全体において、 op(element, mask) == value を満たす個数を返す。 |
| `int count_bitwise(unsigned long long value, int l, int r, unsigned long long mask, BitwiseOperation op) const` | count_equal_bitwise に準ずる | 区間 [l, r) において、 op(element, mask) == value を満たす個数を返す。 |
| `bool contains_bitwise(unsigned long long value, unsigned long long mask, BitwiseOperation op) const` | count_bitwise に準ずる | 配列全体において、 op(element, mask) == value を満たす要素が存在するかを返す。 |
| `std::vector<unsigned long long> values_bitwise(unsigned long long mask, BitwiseOperation op) const` | list_frequencies_bitwise(0, n_, mask, op) に準ずる | 配列全体において現れる op(element, mask) の異なる値を昇順で返す。 |
| `int distinct_size_bitwise(unsigned long long mask, BitwiseOperation op) const` | values_bitwise(mask, op) に準ずる | 配列全体において現れる op(element, mask) の異なる値の個数を返す。 |
| `int index_of_bitwise(unsigned long long value, unsigned long long mask, BitwiseOperation op) const` | values_bitwise(mask, op) に準ずる | values_bitwise(mask, op) が返す昇順 distinct 配列の中で、 value が現れる位置を返す。存在しなければ -1 を返す。 |
| `int select_bitwise(unsigned long long value, int kth, unsigned long long mask, BitwiseOperation op) const` | O(log N * C) / C は rank_bitwise(value, r, mask, op) 1 回分 | 0-indexed で kth 番目に現れる op(element, mask) == value を満たす位置を返す。 存在しなければ -1 を返す。 |
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
| `int count_equal_xor(int l, int r, unsigned long long mask, unsigned long long value) const` | O(W) | 区間 [l, r) において (value xor mask) == value を満たす個数を返す。 |
| `int count_greater_xor(int l, int r, unsigned long long mask, unsigned long long lower) const` | O(W) | 区間 [l, r) において (value xor mask) > lower の個数を返す。 |
| `int count_greater_equal_xor(int l, int r, unsigned long long mask, unsigned long long lower) const` | O(W) | 区間 [l, r) において (value xor mask) >= lower の個数を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_xor( int l, int r, unsigned long long mask, unsigned long long upper ) const` | O(W) | 区間 [l, r) において (value xor mask) < upper を満たす要素について、 (個数, 元の value の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_equal_xor( int l, int r, unsigned long long mask, unsigned long long upper ) const` | O(W) | 区間 [l, r) において (value xor mask) <= upper を満たす要素について、 (個数, 元の value の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_range_xor( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper ) const` | O(W) | 区間 [l, r) において lower <= (value xor mask) < upper を満たす要素について、 (個数, 元の value の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_xor( int l, int r, unsigned long long mask, unsigned long long lower ) const` | O(W) | 区間 [l, r) において (value xor mask) > lower を満たす要素について、 (個数, 元の value の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_equal_xor( int l, int r, unsigned long long mask, unsigned long long lower ) const` | O(W) | 区間 [l, r) において (value xor mask) >= lower を満たす要素について、 (個数, 元の value の総和) を返す。 |
| `sum_type sum_less_xor(int l, int r, unsigned long long mask, unsigned long long upper) const` | O(W) | 区間 [l, r) において (value xor mask) < upper を満たす要素の、 元の value の総和を返す。 |
| `sum_type sum_less_equal_xor(int l, int r, unsigned long long mask, unsigned long long upper) const` | O(W) | 区間 [l, r) において (value xor mask) <= upper を満たす要素の、 元の value の総和を返す。 |
| `sum_type sum_equal_xor(int l, int r, unsigned long long mask, unsigned long long value) const` | O(W) | 区間 [l, r) において (value xor mask) == value を満たす要素の、 元の value の総和を返す。 |
| `sum_type sum_range_xor( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper ) const` | O(W) | 区間 [l, r) において lower <= (value xor mask) < upper を満たす要素の、 元の value の総和を返す。 |
| `sum_type sum_greater_xor(int l, int r, unsigned long long mask, unsigned long long lower) const` | O(W) | 区間 [l, r) において (value xor mask) > lower を満たす要素の、 元の value の総和を返す。 |
| `sum_type sum_greater_equal_xor(int l, int r, unsigned long long mask, unsigned long long lower) const` | O(W) | 区間 [l, r) において (value xor mask) >= lower を満たす要素の、 元の value の総和を返す。 |
| `unsigned long long access_xor(int i, unsigned long long mask) const` | O(1) | i 番目の要素に (value xor mask) を施した結果を返す。 |
| `int rank_xor(unsigned long long value, int r, unsigned long long mask) const` | O(log σ) | prefix 区間 [0, r) において、(element xor mask) == value を満たす個数を返す。 |
| `int rank_xor(unsigned long long value, int l, int r, unsigned long long mask) const` | O(log σ) | 区間 [l, r) において、(element xor mask) == value を満たす個数を返す。 |
| `int count_xor(unsigned long long value, unsigned long long mask) const` | O(log σ) | 配列全体において、(element xor mask) == value を満たす個数を返す。 |
| `int count_xor(unsigned long long value, int l, int r, unsigned long long mask) const` | O(log σ) | 区間 [l, r) において、(element xor mask) == value を満たす個数を返す。 |
| `bool contains_xor(unsigned long long value, unsigned long long mask) const` | O(log σ) | 配列全体において、(element xor mask) == value を満たす要素が存在するかを返す。 |
| `std::vector<unsigned long long> values_xor(unsigned long long mask) const` | おおよそ O((σ + 1)W) / W は unsigned(T) のビット幅 | 配列全体において現れる (element xor mask) の異なる値を昇順で返す。 |
| `int distinct_size_xor(unsigned long long mask) const` | O(1) | 配列全体において現れる (element xor mask) の異なる値の個数を返す。 XOR は全単射なので元の distinct 数と一致する。 |
| `int index_of_xor(unsigned long long value, unsigned long long mask) const` | おおよそ O((σ + 1)W) / W は unsigned(T) のビット幅 | values_xor(mask) が返す昇順 distinct 配列の中で、value が現れる位置を返す。存在しなければ -1 を返す。 |
| `int select_xor(unsigned long long value, int kth, unsigned long long mask) const` | O(log σ log N) | 0-indexed で kth 番目に現れる (element xor mask) == value を満たす位置を返す。存在しなければ -1 を返す。 |
| `unsigned long long kth_smallest_xor(int l, int r, int k, unsigned long long mask) const` | O(W) | 区間 [l, r) の (element xor mask) を昇順に見たときの k 番目の値を返す。k は 0-indexed。 |
| `unsigned long long kth_largest_xor(int l, int r, int k, unsigned long long mask) const` | O(W) | 区間 [l, r) の (element xor mask) を降順に見たときの k 番目の値を返す。k は 0-indexed。 |
| `unsigned long long quantile_xor(int l, int r, int k, unsigned long long mask) const` | O(W) | kth_smallest_xor の別名。 |
| `std::optional<unsigned long long> min_value_xor(int l, int r, unsigned long long mask) const` | O(W) | 区間 [l, r) の (element xor mask) の最小値を返す。空区間なら std::nullopt。 |
| `std::optional<unsigned long long> max_value_xor(int l, int r, unsigned long long mask) const` | O(W) | 区間 [l, r) の (element xor mask) の最大値を返す。空区間なら std::nullopt。 |
| `std::optional<unsigned long long> median_lower_xor(int l, int r, unsigned long long mask) const` | O(W) | 区間 [l, r) の (element xor mask) の下側中央値を返す。空区間なら std::nullopt。 |
| `std::optional<unsigned long long> median_upper_xor(int l, int r, unsigned long long mask) const` | O(W) | 区間 [l, r) の (element xor mask) の上側中央値を返す。空区間なら std::nullopt。 |
| `std::optional<unsigned long long> next_value_ge_xor(int l, int r, unsigned long long lower, unsigned long long mask) const` | O(W) | 区間 [l, r) の (element xor mask) の中で、lower 以上の最小値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> next_value_gt_xor(int l, int r, unsigned long long lower, unsigned long long mask) const` | O(W) | 区間 [l, r) の (element xor mask) の中で、lower より大きい最小値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> prev_value_lt_xor(int l, int r, unsigned long long upper, unsigned long long mask) const` | O(W) | 区間 [l, r) の (element xor mask) の中で、upper 未満の最大値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> prev_value_le_xor(int l, int r, unsigned long long upper, unsigned long long mask) const` | O(W) | 区間 [l, r) の (element xor mask) の中で、upper 以下の最大値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> next_value_xor(int l, int r, unsigned long long lower, unsigned long long mask) const` | O(W) | next_value_ge_xor の別名。 |
| `std::optional<unsigned long long> prev_value_xor(int l, int r, unsigned long long upper, unsigned long long mask) const` | O(W) | prev_value_lt_xor の別名。 |
| `std::vector<std::pair<unsigned long long, int>> list_frequencies_xor(int l, int r, unsigned long long mask) const` | おおよそ O((z + 1)W) / z は出力される異なる値の個数、W は unsigned(T) のビット幅 | 区間 [l, r) の (element xor mask) に現れる異なる値を (値, 個数) で昇順列挙する。 |
| `std::vector<std::pair<unsigned long long, int>> list_frequencies_xor(int l, int r, unsigned long long lower, unsigned long long upper, unsigned long long mask) const` | おおよそ O((z + 1)W) / z は出力される異なる値の個数、W は unsigned(T) のビット幅 | 区間 [l, r) の (element xor mask) について、値域 [lower, upper) に属する異なる値を (値, 個数) で昇順列挙する。 |
| `std::vector<unsigned long long> distinct_values_xor(int l, int r, unsigned long long mask) const` | list_frequencies_xor(l, r, mask) に準ずる | 区間 [l, r) の (element xor mask) に現れる異なる値を昇順で返す。 |
| `std::vector<std::pair<unsigned long long, int>> top_k_frequent_xor(int l, int r, int k, unsigned long long mask) const` | O((m + 1) log σ + m log k + k log k) / m は区間 [l, r) に現れる異なる元値の個数 | 区間 [l, r) の (element xor mask) に対し、頻度上位 k 個の (値, 出現回数) を返す。 |
| `std::optional<std::pair<unsigned long long, int>> mode_xor(int l, int r, unsigned long long mask) const` | O((m + 1) log σ) / m は区間 [l, r) に現れる異なる元値の個数 | 区間 [l, r) の (element xor mask) の最頻値を (値, 出現回数) で返す。空区間なら std::nullopt。 |
| `std::vector<std::tuple<unsigned long long, int, int>> intersect_xor(int l1, int r1, int l2, int r2, unsigned long long mask) const` | おおよそ O((z + 1)W) / z は共通して現れる異なる値の個数、W は unsigned(T) のビット幅 | 2 区間 [l1, r1), [l2, r2) について、共通して現れる (element xor mask) の値を列挙する。 |

### OR ラッパー

| 関数   | 時間計算量 | 説明 |
|---|---|---|
| `int count_less_or(int l, int r, unsigned long long mask, unsigned long long upper) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value \| mask) < upper の個数を返す。 |
| `int range_freq_or( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper ) const` | count_less_or を 2 回呼ぶので、それに準ずる | 区間 [l, r) において lower <= (value \| mask) < upper の個数を返す。 |
| `int count_less_equal_or(int l, int r, unsigned long long mask, unsigned long long upper) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value \| mask) <= upper の個数を返す。 |
| `int count_equal_or(int l, int r, unsigned long long mask, unsigned long long value) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value \| mask) == value を満たす個数を返す。 |
| `int count_greater_or(int l, int r, unsigned long long mask, unsigned long long lower) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value \| mask) > lower の個数を返す。 |
| `int count_greater_equal_or(int l, int r, unsigned long long mask, unsigned long long lower) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value \| mask) >= lower の個数を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_or( int l, int r, unsigned long long mask, unsigned long long upper ) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value \| mask) < upper を満たす要素について、 (個数, 元の value の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_equal_or( int l, int r, unsigned long long mask, unsigned long long upper ) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value \| mask) <= upper を満たす要素について、 (個数, 元の value の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_range_or( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper ) const` | count_and_sum_less_or を 2 回呼ぶので、それに準ずる | 区間 [l, r) において lower <= (value \| mask) < upper を満たす要素について、 (個数, 元の value の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_or( int l, int r, unsigned long long mask, unsigned long long lower ) const` | count_and_sum_less_equal_or を 1 回呼ぶので、それに準ずる | 区間 [l, r) において (value \| mask) > lower を満たす要素について、 (個数, 元の value の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_equal_or( int l, int r, unsigned long long mask, unsigned long long lower ) const` | count_and_sum_less_or を 1 回呼ぶので、それに準ずる | 区間 [l, r) において (value \| mask) >= lower を満たす要素について、 (個数, 元の value の総和) を返す。 |
| `sum_type sum_less_or(int l, int r, unsigned long long mask, unsigned long long upper) const` | count_and_sum_less_or を 1 回呼ぶので、それに準ずる | 区間 [l, r) において (value \| mask) < upper を満たす要素の、 元の value の総和を返す。 |
| `sum_type sum_less_equal_or(int l, int r, unsigned long long mask, unsigned long long upper) const` | count_and_sum_less_equal_or を 1 回呼ぶので、それに準ずる | 区間 [l, r) において (value \| mask) <= upper を満たす要素の、 元の value の総和を返す。 |
| `sum_type sum_equal_or(int l, int r, unsigned long long mask, unsigned long long value) const` | count_and_sum_less_or を 2 回呼ぶので、それに準ずる | 区間 [l, r) において (value \| mask) == value を満たす要素の、 元の value の総和を返す。 |
| `sum_type sum_range_or( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper ) const` | count_and_sum_less_or を 2 回呼ぶので、それに準ずる | 区間 [l, r) において lower <= (value \| mask) < upper を満たす要素の、 元の value の総和を返す。 |
| `sum_type sum_greater_or(int l, int r, unsigned long long mask, unsigned long long lower) const` | count_and_sum_greater_or を 1 回呼ぶので、それに準ずる | 区間 [l, r) において (value \| mask) > lower を満たす要素の、 元の value の総和を返す。 |
| `sum_type sum_greater_equal_or(int l, int r, unsigned long long mask, unsigned long long lower) const` | count_and_sum_greater_equal_or を 1 回呼ぶので、それに準ずる | 区間 [l, r) において (value \| mask) >= lower を満たす要素の、 元の value の総和を返す。 |
| `unsigned long long access_or(int i, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | i 番目の要素に (value \| mask) を施した結果を返す。 |
| `int rank_or(unsigned long long value, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | prefix 区間 [0, r) において、(value \| mask) == value を満たす個数を返す。 |
| `int rank_or(unsigned long long value, int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において、(value \| mask) == value を満たす個数を返す。 |
| `int count_or(unsigned long long value, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 配列全体において、(value \| mask) == value を満たす個数を返す。 |
| `int count_or(unsigned long long value, int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において、(value \| mask) == value を満たす個数を返す。 |
| `bool contains_or(unsigned long long value, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 配列全体において、(value \| mask) == value を満たす要素が存在するかを返す。 |
| `std::vector<unsigned long long> values_or(unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 配列全体において現れる (value \| mask) の異なる値を昇順で返す。 |
| `int distinct_size_or(unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 配列全体において現れる (value \| mask) の異なる値の個数を返す。 |
| `int index_of_or(unsigned long long value, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | values_or(mask) が返す昇順 distinct 配列の中で、value が現れる位置を返す。存在しなければ -1 を返す。 |
| `int select_or(unsigned long long value, int kth, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 0-indexed で kth 番目に現れる (value \| mask) == value を満たす位置を返す。存在しなければ -1 を返す。 |
| `unsigned long long kth_smallest_or(int l, int r, int k, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) を昇順に見たときの k 番目の値を返す。k は 0-indexed。 |
| `unsigned long long kth_largest_or(int l, int r, int k, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) を降順に見たときの k 番目の値を返す。k は 0-indexed。 |
| `unsigned long long quantile_or(int l, int r, int k, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | kth_smallest_or の別名。 |
| `std::optional<unsigned long long> min_value_or(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) の最小値を返す。空区間なら std::nullopt。 |
| `std::optional<unsigned long long> max_value_or(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) の最大値を返す。空区間なら std::nullopt。 |
| `std::optional<unsigned long long> median_lower_or(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) の下側中央値を返す。空区間なら std::nullopt。 |
| `std::optional<unsigned long long> median_upper_or(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) の上側中央値を返す。空区間なら std::nullopt。 |
| `std::optional<unsigned long long> next_value_ge_or(int l, int r, unsigned long long lower, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) の中で、lower 以上の最小値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> next_value_gt_or(int l, int r, unsigned long long lower, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) の中で、lower より大きい最小値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> prev_value_lt_or(int l, int r, unsigned long long upper, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) の中で、upper 未満の最大値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> prev_value_le_or(int l, int r, unsigned long long upper, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) の中で、upper 以下の最大値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> next_value_or(int l, int r, unsigned long long lower, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | next_value_ge_or の別名。 |
| `std::optional<unsigned long long> prev_value_or(int l, int r, unsigned long long upper, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | prev_value_lt_or の別名。 |
| `std::vector<std::pair<unsigned long long, int>> list_frequencies_or(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) に現れる異なる値を (値, 個数) で昇順列挙する。 |
| `std::vector<std::pair<unsigned long long, int>> list_frequencies_or(int l, int r, unsigned long long lower, unsigned long long upper, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) について、値域 [lower, upper) に属する異なる値を (値, 個数) で昇順列挙する。 |
| `std::vector<unsigned long long> distinct_values_or(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) に現れる異なる値を昇順で返す。 |
| `std::vector<std::pair<unsigned long long, int>> top_k_frequent_or(int l, int r, int k, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) に対し、頻度上位 k 個の (値, 出現回数) を返す。 |
| `std::optional<std::pair<unsigned long long, int>> mode_or(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value \| mask) の最頻値を (値, 出現回数) で返す。空区間なら std::nullopt。 |
| `std::vector<std::tuple<unsigned long long, int, int>> intersect_or(int l1, int r1, int l2, int r2, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 2 区間 [l1, r1), [l2, r2) について、共通して現れる (value \| mask) の値を列挙する。 |

### AND ラッパー


| 関数   | 時間計算量 | 説明 |
|---|---|---|
| `int count_less_and(int l, int r, unsigned long long mask, unsigned long long upper) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value & mask) < upper の個数を返す。 |
| `int range_freq_and( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper ) const` | count_less_and を 2 回呼ぶので、それに準ずる | 区間 [l, r) において lower <= (value & mask) < upper の個数を返す。 |
| `int count_less_equal_and(int l, int r, unsigned long long mask, unsigned long long upper) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value & mask) <= upper の個数を返す。 |
| `int count_equal_and(int l, int r, unsigned long long mask, unsigned long long value) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value & mask) == value を満たす個数を返す。 |
| `int count_greater_and(int l, int r, unsigned long long mask, unsigned long long lower) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value & mask) > lower の個数を返す。 |
| `int count_greater_equal_and(int l, int r, unsigned long long mask, unsigned long long lower) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value & mask) >= lower の個数を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_and( int l, int r, unsigned long long mask, unsigned long long upper ) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value & mask) < upper を満たす要素について、 (個数, 元の value の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_less_equal_and( int l, int r, unsigned long long mask, unsigned long long upper ) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において (value & mask) <= upper を満たす要素について、 (個数, 元の value の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_range_and( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper ) const` | count_and_sum_less_and を 2 回呼ぶので、それに準ずる | 区間 [l, r) において lower <= (value & mask) < upper を満たす要素について、 (個数, 元の value の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_and( int l, int r, unsigned long long mask, unsigned long long lower ) const` | count_and_sum_less_equal_and を 1 回呼ぶので、それに準ずる | 区間 [l, r) において (value & mask) > lower を満たす要素について、 (個数, 元の value の総和) を返す。 |
| `std::pair<int, sum_type> count_and_sum_greater_equal_and( int l, int r, unsigned long long mask, unsigned long long lower ) const` | count_and_sum_less_and を 1 回呼ぶので、それに準ずる | 区間 [l, r) において (value & mask) >= lower を満たす要素について、 (個数, 元の value の総和) を返す。 |
| `sum_type sum_less_and(int l, int r, unsigned long long mask, unsigned long long upper) const` | count_and_sum_less_and を 1 回呼ぶので、それに準ずる | 区間 [l, r) において (value & mask) < upper を満たす要素の、 元の value の総和を返す。 |
| `sum_type sum_less_equal_and(int l, int r, unsigned long long mask, unsigned long long upper) const` | count_and_sum_less_equal_and を 1 回呼ぶので、それに準ずる | 区間 [l, r) において (value & mask) <= upper を満たす要素の、 元の value の総和を返す。 |
| `sum_type sum_equal_and(int l, int r, unsigned long long mask, unsigned long long value) const` | count_and_sum_less_and を 2 回呼ぶので、それに準ずる | 区間 [l, r) において (value & mask) == value を満たす要素の、 元の value の総和を返す。 |
| `sum_type sum_range_and( int l, int r, unsigned long long mask, unsigned long long lower, unsigned long long upper ) const` | count_and_sum_less_and を 2 回呼ぶので、それに準ずる | 区間 [l, r) において lower <= (value & mask) < upper を満たす要素の、 元の value の総和を返す。 |
| `sum_type sum_greater_and(int l, int r, unsigned long long mask, unsigned long long lower) const` | count_and_sum_greater_and を 1 回呼ぶので、それに準ずる | 区間 [l, r) において (value & mask) > lower を満たす要素の、 元の value の総和を返す。 |
| `sum_type sum_greater_equal_and(int l, int r, unsigned long long mask, unsigned long long lower) const` | count_and_sum_greater_equal_and を 1 回呼ぶので、それに準ずる | 区間 [l, r) において (value & mask) >= lower を満たす要素の、 元の value の総和を返す。 |
| `unsigned long long access_and(int i, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | i 番目の要素に (value & mask) を施した結果を返す。 |
| `int rank_and(unsigned long long value, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | prefix 区間 [0, r) において、(value & mask) == value を満たす個数を返す。 |
| `int rank_and(unsigned long long value, int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において、(value & mask) == value を満たす個数を返す。 |
| `int count_and(unsigned long long value, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 配列全体において、(value & mask) == value を満たす個数を返す。 |
| `int count_and(unsigned long long value, int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) において、(value & mask) == value を満たす個数を返す。 |
| `bool contains_and(unsigned long long value, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 配列全体において、(value & mask) == value を満たす要素が存在するかを返す。 |
| `std::vector<unsigned long long> values_and(unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 配列全体において現れる (value & mask) の異なる値を昇順で返す。 |
| `int distinct_size_and(unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 配列全体において現れる (value & mask) の異なる値の個数を返す。 |
| `int index_of_and(unsigned long long value, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | values_and(mask) が返す昇順 distinct 配列の中で、value が現れる位置を返す。存在しなければ -1 を返す。 |
| `int select_and(unsigned long long value, int kth, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 0-indexed で kth 番目に現れる (value & mask) == value を満たす位置を返す。存在しなければ -1 を返す。 |
| `unsigned long long kth_smallest_and(int l, int r, int k, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) を昇順に見たときの k 番目の値を返す。k は 0-indexed。 |
| `unsigned long long kth_largest_and(int l, int r, int k, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) を降順に見たときの k 番目の値を返す。k は 0-indexed。 |
| `unsigned long long quantile_and(int l, int r, int k, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | kth_smallest_and の別名。 |
| `std::optional<unsigned long long> min_value_and(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) の最小値を返す。空区間なら std::nullopt。 |
| `std::optional<unsigned long long> max_value_and(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) の最大値を返す。空区間なら std::nullopt。 |
| `std::optional<unsigned long long> median_lower_and(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) の下側中央値を返す。空区間なら std::nullopt。 |
| `std::optional<unsigned long long> median_upper_and(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) の上側中央値を返す。空区間なら std::nullopt。 |
| `std::optional<unsigned long long> next_value_ge_and(int l, int r, unsigned long long lower, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) の中で、lower 以上の最小値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> next_value_gt_and(int l, int r, unsigned long long lower, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) の中で、lower より大きい最小値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> prev_value_lt_and(int l, int r, unsigned long long upper, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) の中で、upper 未満の最大値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> prev_value_le_and(int l, int r, unsigned long long upper, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) の中で、upper 以下の最大値を返す。存在しなければ std::nullopt。 |
| `std::optional<unsigned long long> next_value_and(int l, int r, unsigned long long lower, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | next_value_ge_and の別名。 |
| `std::optional<unsigned long long> prev_value_and(int l, int r, unsigned long long upper, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | prev_value_lt_and の別名。 |
| `std::vector<std::pair<unsigned long long, int>> list_frequencies_and(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) に現れる異なる値を (値, 個数) で昇順列挙する。 |
| `std::vector<std::pair<unsigned long long, int>> list_frequencies_and(int l, int r, unsigned long long lower, unsigned long long upper, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) について、値域 [lower, upper) に属する異なる値を (値, 個数) で昇順列挙する。 |
| `std::vector<unsigned long long> distinct_values_and(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) に現れる異なる値を昇順で返す。 |
| `std::vector<std::pair<unsigned long long, int>> top_k_frequent_and(int l, int r, int k, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) に対し、頻度上位 k 個の (値, 出現回数) を返す。 |
| `std::optional<std::pair<unsigned long long, int>> mode_and(int l, int r, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 区間 [l, r) の (value & mask) の最頻値を (値, 出現回数) で返す。空区間なら std::nullopt。 |
| `std::vector<std::tuple<unsigned long long, int, int>> intersect_and(int l1, int r1, int l2, int r2, unsigned long long mask) const` | O(M) / M は訪問ノード数（最悪 O(2^W)） | 2 区間 [l1, r1), [l2, r2) について、共通して現れる (value & mask) の値を列挙する。 |


# Weighted Wavelet Matrix
[ライブラリのリンク](https://github.com/harurunrunrun/GPT-generated-Library/blob/main/library/WeightedWaveletMatrixMonoid.hpp)


