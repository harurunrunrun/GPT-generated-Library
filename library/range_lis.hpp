#ifndef RANGE_LIS_HPP
#define RANGE_LIS_HPP

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

// Static range (weighted) LIS, following Korogi, arXiv:2610.20983v1.
//
// For any half-open interval [l, r), query(l, r) returns the maximum total
// weight of an increasing subsequence (strict or nondecreasing). With unit
// weights, this is its length. restore(l, r) returns original 0-based indices.
//
// Preprocessing: O(n log^2 n) time and space.
// query:        O(log n).
// restore:      O(log n + k), where k is the number of returned indices.
//
// Requirements:
//   * The input array is immutable after construction.
//   * All weights must be strictly positive.
//   * Value must be strictly ordered by operator<.
//   * Weight must support Weight{0}, Weight{1}, +, and >; additions and
//     comparisons must be exact for the stated mathematical guarantees.
//   * If Weight is an integer, the total weights must not overflow.
//
// Equal input values are handled according to LISMode. At equal optimum
// weights, a valid optimum is returned; lexicographically smallest indices
// are not guaranteed.
enum class LISMode { Strict, NonDecreasing };

template <class Value, class Weight = long long>
class RangeLIS {
public:
    explicit RangeLIS(const std::vector<Value>& a,
                      LISMode mode = LISMode::Strict)
        : RangeLIS(a, std::vector<Weight>(a.size(), Weight{1}), mode) {}

    RangeLIS(const std::vector<Value>& a, const std::vector<Weight>& w,
             LISMode mode = LISMode::Strict)
        : n_(static_cast<int>(a.size())) {
        if (a.size() != w.size()) {
            throw std::invalid_argument("RangeLIS: values and weights differ in size");
        }
        for (const auto& x : w) {
            if (!(x > Weight{0})) {
                throw std::invalid_argument("RangeLIS: weights must be positive");
            }
        }

        witnesses_.push_back(Witness{0, 0, -1, Weight{0}}); // Empty witness.
        functions_.push_back(Function{0, 0, 0, Weight{0}}); // Null sentinel.
        roots_.assign(n_ + 1, 0);
        if (n_ == 0) return;

        // Strict:        rank by (value, -index), forbidding equal successors.
        // NonDecreasing: rank by (value, +index), allowing equal successors.
        // Since selected indices must increase, either ranking makes the
        // requested LIS exactly a strict LIS on distinct ranks.
        std::vector<int> order(n_);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int i, int j) {
            if (a[i] < a[j]) return true;
            if (a[j] < a[i]) return false;
            return mode == LISMode::Strict ? i > j : i < j;
        });
        std::vector<int> rank(n_);
        for (int i = 0; i < n_; ++i) rank[order[i]] = i;

        zero_root_ = build_zero(0, n_);
        outer_.assign(4 * n_ + 1, zero_root_);
        roots_[0] = zero_root_;

        for (int i = 0; i < n_; ++i) {
            const int predecessor = prefix_envelope(1, 0, n_, rank[i], zero_root_);
            const int item = make_leaf(i, w[i]);
            const int endpoint = prefix_append(predecessor, 0, n_, i + 1, item);
            outer_set(1, 0, n_, rank[i], endpoint);
            roots_[i + 1] = outer_[1];
        }
    }

    int size() const noexcept { return n_; }

    // Maximum total weight on [l, r). For unit weights, the LIS length.
    Weight query(int l, int r) const {
        check_range(l, r);
        if (l == r) return Weight{0};
        Weight answer{0};
        int root = roots_[r];
        int b = 0, e = n_;
        while (true) {
            const Function& t = functions_[root];
            answer = answer + witnesses_[t.suffix].weight;
            if (e - b == 1) break;
            const int m = b + (e - b) / 2;
            if (l < m) {
                root = t.left;
                e = m;
            } else {
                root = t.right;
                b = m;
            }
        }
        return answer;
    }

    // Indices of an optimal subsequence (according to mode) on [l, r).
    // The output indices are in increasing order.
    std::vector<int> restore(int l, int r) const {
        return query_with_indices(l, r).second;
    }

    // Calculate the maximum weight and reconstruct an optimum in one query.
    std::pair<Weight, std::vector<int>> query_with_indices(int l, int r) const {
        check_range(l, r);
        if (l == r) return {Weight{0}, {}};

        Weight answer{0};
        std::vector<int> suffixes;
        int root = roots_[r];
        int b = 0, e = n_;
        while (true) {
            const Function& t = functions_[root];
            answer = answer + witnesses_[t.suffix].weight;
            if (t.suffix != 0) suffixes.push_back(t.suffix);
            if (e - b == 1) break;
            const int m = b + (e - b) / 2;
            if (l < m) {
                root = t.left;
                e = m;
            } else {
                root = t.right;
                b = m;
            }
        }

        // A leaf's suffix is appended before its ancestors' suffixes.
        // Expand the persistent concatenation DAG without recursive calls.
        std::vector<int> indices;
        std::vector<int> stack;
        for (auto it = suffixes.rbegin(); it != suffixes.rend(); ++it) {
            stack.push_back(*it);
            while (!stack.empty()) {
                const int id = stack.back();
                stack.pop_back();
                const Witness& t = witnesses_[id];
                if (t.index >= 0) {
                    indices.push_back(t.index);
                } else {
                    stack.push_back(t.right);
                    stack.push_back(t.left);
                }
            }
        }
        return {answer, std::move(indices)};
    }

private:
    struct Witness {
        int left, right; // Nonzero for a concatenation node.
        int index;       // >= 0 for an original array element.
        Weight weight;
    };
    struct Function {
        int left, right; // Children in the persistent tree over left endpoints.
        int suffix;      // Lazy suffix witness applied to the whole subtree.
        Weight last;     // Score at the rightmost left-endpoint coordinate.
    };

    int n_ = 0;
    int zero_root_ = 0;
    std::vector<Witness> witnesses_;
    std::vector<Function> functions_;
    std::vector<int> outer_; // Mutable segment tree indexed by value rank.
    std::vector<int> roots_; // Persistent root after each input prefix.

    void check_range(int l, int r) const {
        if (l < 0 || l > r || r > n_) {
            throw std::out_of_range("RangeLIS: expected 0 <= l <= r <= n");
        }
    }

    int make_leaf(int index, const Weight& weight) {
        const int id = static_cast<int>(witnesses_.size());
        witnesses_.push_back(Witness{0, 0, index, weight});
        return id;
    }

    int concat(int a, int b) {
        if (a == 0) return b;
        if (b == 0) return a;
        const Weight sum = witnesses_[a].weight + witnesses_[b].weight;
        const int id = static_cast<int>(witnesses_.size());
        witnesses_.push_back(Witness{a, b, -1, sum});
        return id;
    }

    int make_function(int left, int right, int suffix, const Weight& last) {
        const int id = static_cast<int>(functions_.size());
        functions_.push_back(Function{left, right, suffix, last});
        return id;
    }

    int build_zero(int b, int e) {
        if (e - b == 1) return make_function(0, 0, 0, Weight{0});
        const int m = b + (e - b) / 2;
        const int left = build_zero(b, m);
        const int right = build_zero(m, e);
        return make_function(left, right, 0, Weight{0});
    }

    int append(int root, int suffix) {
        if (suffix == 0) return root;
        const Function t = functions_[root]; // Copy before allocating nodes.
        return make_function(t.left, t.right, concat(t.suffix, suffix),
                             t.last + witnesses_[suffix].weight);
    }

    std::pair<int, int> children(int root) {
        const Function t = functions_[root];
        return {append(t.left, t.suffix), append(t.right, t.suffix)};
    }

    int join(int left, int right) {
        return make_function(left, right, 0, functions_[right].last);
    }

    int prefix_append(int root, int b, int e, int c, int suffix) {
        if (c <= b) return root;
        if (e <= c) return append(root, suffix);
        const auto [left, right] = children(root);
        const int m = b + (e - b) / 2;
        if (c <= m) return join(prefix_append(left, b, m, c, suffix), right);
        return join(append(left, suffix),
                    prefix_append(right, m, e, c, suffix));
    }

    // Pointwise max of two envelopes whose endpoint values are separated:
    // the second envelope wins strictly on a (possibly empty) prefix.
    int merge(int f, int g, int b, int e) {
        if (f == g || g == zero_root_) return f;
        if (f == zero_root_) return g;
        if (e - b == 1) {
            return (functions_[g].last > functions_[f].last) ? g : f;
        }
        const auto [fl, fr] = children(f);
        const auto [gl, gr] = children(g);
        const int m = b + (e - b) / 2;
        if (functions_[gl].last > functions_[fl].last) {
            return join(gl, merge(fr, gr, m, e));
        }
        return join(merge(fl, gl, b, m), fr);
    }

    // Fold full value-rank intervals from left to right.
    int prefix_envelope(int node, int b, int e, int c, int accum) {
        if (c <= b) return accum;
        if (e <= c) return merge(accum, outer_[node], 0, n_);
        const int m = b + (e - b) / 2;
        accum = prefix_envelope(node * 2, b, m, c, accum);
        return prefix_envelope(node * 2 + 1, m, e, c, accum);
    }

    void outer_set(int node, int b, int e, int pos, int value) {
        if (e - b == 1) {
            outer_[node] = value;
            return;
        }
        const int m = b + (e - b) / 2;
        if (pos < m) outer_set(node * 2, b, m, pos, value);
        else outer_set(node * 2 + 1, m, e, pos, value);
        outer_[node] = merge(outer_[node * 2], outer_[node * 2 + 1], 0, n_);
    }
};

#endif // RANGE_LIS_HPP
