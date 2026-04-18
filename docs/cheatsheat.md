# 区間更新チートシート

## 0. 最初に見る表

| 症状 | 切替先 |
|---|---|
| `mapping(f, op(x, y))` が書ける | 遅延セグ木 |
| 少数の補助情報を足せば閉じる | 拡張遅延セグ木 |
| `max1/max2/cntMax` や `min1/min2/cntMin` で処理できる | Segment Tree Beats |
| 値の分布そのものが要る | Wavelet Matrix / Binary Trie / 永続セグ木 / 平方分割 |
| 添字順が変わる | Implicit Treap / Splay |
| 差分にすると端点更新になる | 差分 + Fenwick / Segtree |
| オンラインが重い | Mo / CDQ / 永続化 / 並列二分探索 |
| 最悪計算量保証を捨ててもよい | ODT |

## 1. 遅延セグ木の判定

| 条件 | 判定 |
|---|---|
| `mapping(f, op(x, y)) = op(mapping(f, x), mapping(f, y))` | 区間をまとめた後に更新できる |
| `mapping(f, mapping(g, x)) = mapping(composition(f, g), x)` | タグ合成ができる |
| `id` がある | 何もしない更新がある |
| `S` が定数サイズ | ノード情報が有限個で済む |

## 2. `assign / set`

| query | lazy | `S` | `F` | `mapping` | `composition` | lazyで不可なら |
|---|---:|---|---|---|---|---|
| sum | 可 | `sum,len` | `set?` | `sum=c*len` | 新しい代入が優先 |  |
| min / max | 可 | `mn,mx` | `set?` | `mn=mx=c` | 新しい代入が優先 |  |
| xor | 可 | `xr,len` | `set?` | `xr=(len&1?c:0)` | 新しい代入が優先 |  |
| and / or | 可 | `agg` | `set?` | `agg=c` | 新しい代入が優先 |  |
| gcd | 可 | `g` | `set?` | `g=abs(c)` | 新しい代入が優先 |  |
| hash | 可 | `hash,len` | `set?` | `hash=const_hash(c,len)` | 新しい代入が優先 |  |
| max-subarray | 可 | `sum,pref,suf,best,len` | `set?` | 定数列として再計算 | 新しい代入が優先 |  |
| 0/1 longest run | 可 | `len,pref0/1,suf0/1,best0/1` | `set?` | 0列 / 1列として再計算 | 新しい代入が優先 |  |
| 小値域 freq | 可 | `cnt[σ],len` | `set?` | `cnt[c]=len` | 新しい代入が優先 |  |
| kth / median | 不可寄り |  |  |  |  | Wavelet Matrix / 永続セグ木 / 平方分割 |
| mex | 不可寄り |  |  |  |  | 値域セグ木 / 平方分割 / オフライン |
| mode | 不可寄り |  |  |  |  | 平方分割 / Mo / オフライン |
| distinct | 不可寄り |  |  |  |  | 平方分割 / Mo / オフライン |

## 3. `add`

| query | lazy | `S` | `F` | `mapping` | `composition` | lazyで不可なら |
|---|---:|---|---|---|---|---|
| sum | 可 | `sum,len` | `add` | `sum+=d*len` | `+` |  |
| min / max | 可 | `mn,mx` | `add` | `mn+=d,mx+=d` | `+` |  |
| min/max + count | 可 | `mn,mx,cnt...` | `add` | 値だけ平行移動 | `+` |  |
| moments | 可 | `m1,m2,...,len` | `add` | 二項展開 | `+` |  |
| hash(線形) | 可 | `hash,wSum` | `add` | `hash+=d*wSum` | `+` |  |
| gcd | 非推奨 |  |  |  |  | 差分 + Fenwick/Segtree |
| max-subarray | 不可寄り |  |  |  |  | 問題特化 / 別視点 |
| kth / median | 不可寄り |  |  |  |  | 平方分割 |
| mex | 不可寄り |  |  |  |  | 平方分割 / オフライン |
| mode | 不可寄り |  |  |  |  | 平方分割 / Mo / オフライン |
| distinct | 不可寄り |  |  |  |  | 平方分割 / Mo / オフライン |

### `add + gcd`

| 項目 | 形 |
|---|---|
| 差分 | `d[i]=a[i]-a[i-1]` |
| range add `[l,r]+=x` | `d[l]+=x`, `d[r+1]-=x` |
| `a[l]` | BIT / Segtree の prefix sum |
| `gcd(a[l..r])` | `gcd(a[l], d[l+1],...,d[r])` |

## 4. `mul` / `affine ax+b`

| query | lazy | `S` | `F` | `mapping` | `composition` | lazyで不可なら |
|---|---:|---|---|---|---|---|
| sum | 可 | `sum,len` | `{a,b}` | `sum=a*sum+b*len` | `{a1,b1}∘{a2,b2}={a1*a2,a1*b2+b1}` |  |
| min / max (`a>=0`) | 可 | `mn,mx` | `{a,b}` | `mn=a*mn+b, mx=a*mx+b` | 同上 |  |
| min / max (`a` 負あり) | 可 | `mn,mx` | `{a,b}` | 負なら swap 後に更新 | 同上 |  |
| moments | 可 | `m1,m2,...,len` | `{a,b}` | 二項展開 | 同上 |  |
| gcd (`b=0`) | 可 | `g` | `mul` | `g*=abs(a)` | `*` |  |
| max-subarray (`b=0`) | 条件付き | `sum,max/min pref/suf/sub` | `mul` | 負なら max/min 入替 | `*` |  |
| max-subarray (`b!=0`) | 不可寄り |  |  |  |  | 問題特化 |
| kth / median | 不可寄り |  |  |  |  | 平方分割 |
| mex / mode / distinct | 不可寄り |  |  |  |  | 平方分割 / Mo / オフライン |

## 5. `xor`

| query | lazy | `S` | `F` | `mapping` | `composition` | lazyで不可なら |
|---|---:|---|---|---|---|---|
| xor | 可 | `xr,len` | `mask` | `xr ^= (len&1?mask:0)` | `xor` |  |
| sum / popcount | 可 | `len,cnt[bit]` | `mask` | `mask` の立つ bit で `cnt[b]=len-cnt[b]` | `xor` |  |
| and / or | 可 | `len,cnt[bit]` | `mask` | 同上 | `xor` |  |
| 0/1 longest run | 可 | `len,pref0/1,suf0/1,best0/1` | `flip` | 0系と1系を swap | `xor` / toggle |  |
| 0/1 inversion | 可 | `cnt0,cnt1,inv01,inv10` | `flip` | 0/1, `inv01/inv10` を swap | `xor` / toggle |  |
| max / min | 不可寄り |  |  |  |  | Binary Trie / 平方分割 |
| kth / median | 不可寄り |  |  |  |  | Wavelet Matrix / Trie / 平方分割 |
| mex | 不可寄り |  |  |  |  | Trie / 平方分割 / 問題特化 |
| mode / distinct | 不可寄り |  |  |  |  | 平方分割 / Mo / オフライン |

## 6. `and / or`

| query | lazy | `S` | `F` | `mapping` | `composition` | lazyで不可なら |
|---|---:|---|---|---|---|---|
| sum / xor / and / or / popcount | 可 | `len,cnt[bit]` | `mask` または `{A,B}` | bit ごとに 0固定 / 1固定 / 維持 | `or`, `and`, あるいは `{A,B}` 合成 |  |
| max / min | 不可寄り |  |  |  |  | 平方分割 / 問題特化 |
| kth / median | 不可寄り |  |  |  |  | 平方分割 |
| mex / mode / distinct | 不可寄り |  |  |  |  | 平方分割 / Mo / オフライン |

### `and/or` 混在

| 項目 | 形 |
|---|---|
| タグ | `{A,B}` |
| 作用 | `x -> (x & A) | B` |
| 合成 | `{Af,Bf}∘{Ag,Bg} = {Af&Ag, (Bg&Af)|Bf}` |
| 単位元 | `{all_ones, 0}` |

## 7. `chmin / chmax / clamp`

| query | lazy | 切替先 | 持つもの | 条件 / 備考 |
|---|---:|---|---|---|
| `range chmin + sum/max` | 不可 | Segment Tree Beats | `sum,max1,max2,cntMax` | `max2 < x < max1` で一括更新 |
| `range chmax + sum/min` | 不可 | Segment Tree Beats | `sum,min1,min2,cntMin` | `min1 < x < min2` で一括更新 |
| `chmin/chmax/add + sum/min/max` | 不可 | Segment Tree Beats | 上の両側全部 |  |
| kth / median | 不可寄り | 平方分割 / 問題特化 |  |  |
| mex / mode / distinct | 不可寄り | 平方分割 / Mo / オフライン |  |  |

## 8. `mod`

| query | lazy | 切替先 | 持つもの | 条件 / 備考 |
|---|---:|---|---|---|
| `range mod + sum/max` | 不可 | Beats 系 | `sum,max` | `max < x` で打ち切り |
| `range mod + min/max` | 不可寄り | Beats 系 / 問題特化 | `min,max` | 実装依存 |
| kth / median / mex / mode / distinct | 不可寄り | 平方分割 / Mo / 問題特化 |  |  |

## 9. `reverse / rotate / cut-paste / insert / erase`

| 操作 | lazy | 切替先 | 持つもの |
|---|---:|---|---|
| reverse | 不可 | Implicit Treap / Splay | subtree aggregate + reverse flag |
| rotate | 不可 | Implicit Treap / Splay | split / merge |
| cut-paste | 不可 | Implicit Treap / Splay | split / merge |
| insert / erase | 不可 | Implicit Treap / Splay | size + aggregate |

| query | Treap / Splay で持つもの |
|---|---|
| sum / min / max / gcd | monoid |
| hash | forward / backward hash |
| max-subarray | `sum,pref,suf,best` |
| 非可換 fold | forward / backward |

## 10. 値分布系

| 欲しいもの | 切替先 |
|---|---|
| kth / median / quantile | Wavelet Matrix / 永続セグ木 / 平方分割 |
| mex | 値域セグ木 / 平方分割 / オフライン |
| mode | Mo / 平方分割 / 問題特化 |
| distinct count | Mo / オフライン / 平方分割 |
| predecessor / successor | ordered set / Wavelet Matrix / 平方分割 |
| max xor | Binary Trie / persistent trie / block trie |

## 11. 切替先ごとの手順

### 11.1 Binary Trie

| 手順 | 内容 |
|---|---|
| 1 | 値を bit 列で見る |
| 2 | 各ノードに `child[0/1]`, `cnt` を持つ |
| 3 | 挿入/削除を用意する |
| 4 | `xor mask` が全体一括なら root に lazy xor を持つ |
| 5 | `max xor(x)` は上位 bit から `bit^1` を優先して降りる |
| 6 | `min xor(x)` は上位 bit から `bit` を優先して降りる |
| 7 | `kth` が要るなら各枝の `cnt` で降りる |
| 8 | 区間 query にするなら `prefix persistent trie` にするか block ごとに trie を持つ |

| 向いている場面 | 形 |
|---|---|
| 全体集合に対する `xor`, `max xor`, `min xor`, `kth` | trie 単体 |
| `range [l,r]` に対する `max xor` | 永続 trie の差分 |
| `range xor update` が入る | 平方分割 + block trie |

### 11.2 Wavelet Matrix

| 手順 | 内容 |
|---|---|
| 1 | 更新が少ないか無いかを確認する |
| 2 | 値を座標圧縮する |
| 3 | bit ごとに配列を安定分割して構築する |
| 4 | `rank` を使って区間 `[l,r)` を各 level で写す |
| 5 | `kth`, `count <= x`, `range freq` を処理する |
| 6 | 更新が点更新なら別構造か再構築を考える |
| 7 | 区間更新があるならたいてい切る |

| 向いている場面 | 形 |
|---|---|
| static range kth | Wavelet Matrix |
| static range count | Wavelet Matrix |
| 区間更新あり | 平方分割 / 別方向 |

### 11.3 平方分割

| 手順 | 内容 |
|---|---|
| 1 | 配列を block に切る |
| 2 | 各 block に生配列と補助構造を持つ |
| 3 | 補助構造は query に応じて選ぶ |
| 4 | `range update` は丸ごと block に lazy を付ける |
| 5 | 端 block は展開して直接更新する |
| 6 | 端 block を再構築する |
| 7 | `kth` 系は block ごとの個数集計で二分探索する |
| 8 | `max xor` 系は block ごとに trie を持つ |

| query | block で持つもの |
|---|---|
| `kth`, `count <= x` | sort 済み配列 |
| `mex` | freq map / bucket |
| `max xor` | trie |
| `distinct` | hash map / freq map |

### 11.4 Segment Tree Beats

| 手順 | 内容 |
|---|---|
| 1 | 更新が `chmin`, `chmax`, `mod` などか確認する |
| 2 | `max1,max2,cntMax` / `min1,min2,cntMin` を持つ |
| 3 | ノードで一括更新できる条件を書く |
| 4 | 条件を満たさないときだけ子へ降りる |
| 5 | `sum` も同時に直す |
| 6 | push/pull の不変量を壊さないか確認する |
| 7 | `add` 混在なら両側の値を全部平行移動する |

| 更新 | 一括更新条件 |
|---|---|
| `chmin(x)` | `max2 < x < max1` |
| `chmax(x)` | `min1 < x < min2` |
| `mod(x)` | `max < x` なら打ち切り |

### 11.5 差分 + Fenwick / Segtree

| 手順 | 内容 |
|---|---|
| 1 | 更新を差分に落とせるか見る |
| 2 | `d[i]=a[i]-a[i-1]` を作る |
| 3 | 区間加算を端点 2 個の更新にする |
| 4 | 値そのものは prefix sum で戻す |
| 5 | `gcd` や符号変化などを `d` 上で取る |

| 典型 | 形 |
|---|---|
| `range add + point query` | BIT 1 本 |
| `range add + range gcd` | BIT + gcd segtree |
| `range add + first nonzero` | 差分 segtree |

### 11.6 Implicit Treap / Splay

| 手順 | 内容 |
|---|---|
| 1 | 添字順が動くか確認する |
| 2 | 各ノードに `size` と集約値を持つ |
| 3 | `split(root, k)` を作る |
| 4 | `merge(a, b)` を作る |
| 5 | 区間操作は `split` で切り出して処理する |
| 6 | `reverse` は lazy flag を反転する |
| 7 | 非可換なら forward/backward の両方を持つ |
| 8 | 処理後に `merge` で戻す |

| 操作 | 形 |
|---|---|
| reverse `[l,r)` | `A,B,C = split`; `B.rev ^= 1`; `merge` |
| cut-paste | `split` で抜く → 挿す |
| range fold | 区間ノードの aggregate を読む |

### 11.7 ODT

| 手順 | 内容 |
|---|---|
| 1 | 同値区間の集合で管理する |
| 2 | 更新位置で区間を split する |
| 3 | 更新範囲の区間を列挙する |
| 4 | `assign` なら区間を潰して 1 本にまとめる |
| 5 | `kth` なら列挙区間の `(値, 長さ)` を集めて処理する |
| 6 | `sum` なら列挙区間から計算する |
| 7 | 最悪計算量保証は捨てる |

| 向いている場面 | 形 |
|---|---|
| `assign + kth` | ODT |
| `assign + sum` | ODT |
| `add + kth` | ケースによる |

### 11.8 Mo / オフライン

| 手順 | 内容 |
|---|---|
| 1 | 更新をオンラインで処理しない |
| 2 | query を順序替えする |
| 3 | `add/remove` を書く |
| 4 | 更新付きなら time 次元を持つ |
| 5 | `mode`, `distinct`, 頻度系を処理する |

| 向いている場面 | 形 |
|---|---|
| `distinct` | Mo |
| `mode` | Mo |
| 更新付き頻度 query | Mo with updates |

## 12. 迷ったときの切替表

| 症状 | 切替先 |
|---|---|
| `mapping` が書けない | 遅延セグ木を切る |
| `composition` が壊れる | 遅延セグ木を切る |
| `S` が長さごとの最適値になりそう | 分布系 / 問題特化 |
| 2番目最大値/最小値が効きそう | Segment Tree Beats |
| 値の大小より値そのものの集合が欲しい | Wavelet Matrix / Trie / 平方分割 |
| 更新で添字順が変わる | Implicit Treap / Splay |
| 加算で gcd を取りたい | 差分 + Fenwick / Segtree |
| range assign を雑にさばきたい | ODT |
| オフラインでよい | Mo / CDQ / 永続化 / 並列二分探索 |

## 13. 典型対応表

| 問題型 | 切替先 |
|---|---|
| range add + range sum | 遅延セグ木 |
| range add + range min/max | 遅延セグ木 |
| range add + range gcd | 差分 + Fenwick / Segtree |
| range add + range kth | 平方分割 |
| range assign + range sum/min/max | 遅延セグ木 |
| range assign + range max-subarray | 遅延セグ木 |
| range assign + range kth | 平方分割 / ODT |
| range xor + range xor | 遅延セグ木 |
| range xor + range sum | bit-count 遅延セグ木 |
| range xor + range max | Binary Trie / 平方分割 |
| range xor + range kth | Wavelet Matrix / Trie / 平方分割 |
| range and/or + range sum | bit-count 遅延セグ木 |
| range chmin + range sum/max | Segment Tree Beats |
| range chmax + range sum/min | Segment Tree Beats |
| range mod + range sum | Beats 系 |
| reverse + range hash | Implicit Treap / Splay |
| cut-paste + range sum | Implicit Treap / Splay |
