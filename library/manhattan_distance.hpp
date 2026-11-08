#include <bits/stdc++.h>
using namespace std;

using i64 = long long;
using i128 = __int128_t;
using u128 = __uint128_t;

// Static weighted L1-ball sums over integer lattice points.
// Compute the answer centered at EVERY input point, in input order.
// The radius K is the same fixed, nonnegative integer for all points.
// Coordinates/weights are signed 64-bit. Weighted sums must fit into int64_t;
// intermediate coordinate arithmetic uses 128-bit integers.
//
// Input:
//   N M K
//   M lines: x_1 x_2 ... x_N weight
// Output: M weighted sums, one per input point (including its own weight).
//
// N=1: sorted prefix sums, O(M log M) total.
// N=2: transform L1 balls to squares; sweep + Fenwick, O(M log M) total.
// N>=3: aggregate duplicates. For small stencils, precompute nearby offsets
//       and use compact flat hashing. The padded mixed-radix encoding is
//       collision-free when it fits in 64 bits; otherwise matches are verified.
//       Expected time O(N*M + U*V) for U unique points and V lattice offsets,
//       memory O(N*M + N*V), when this method is selected.
//       For large stencils, partition dimensions into K+1 blocks (when
//       applicable), or use spatial-tree pruning/subtree sums. Their runtime
//       depends on point distribution; worst case can be O(N*U^2).

static u128 absdiff(i64 a, i64 b) {
    i128 z = (i128)a - b;
    return (u128)(z < 0 ? -z : z);
}

static vector<i64> computeOneDim(const vector<i64>& c, const vector<i64>& w, i64 k) {
    struct Point { i64 x, weight; int id; };
    const int m = (int)w.size();
    vector<Point> points(m);
    for (int i = 0; i < m; ++i) points[i] = {c[i], w[i], i};
    sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
        return a.x < b.x;
    });
    vector<i64> prefix(m + 1, 0), answer(m);
    for (int i = 0; i < m; ++i) prefix[i + 1] = prefix[i] + points[i].weight;
    int left = 0, right = 0;
    for (int i = 0; i < m; ++i) {
        while (right < m && (i128)points[right].x <= (i128)points[i].x + k) ++right;
        while (left < right && (i128)points[left].x < (i128)points[i].x - k) ++left;
        answer[points[i].id] = prefix[right] - prefix[left];
    }
    return answer;
}

class Fenwick {
    vector<i64> bit;

public:
    explicit Fenwick(int n) : bit(n + 1, 0) {}

    void add(int p, i128 delta) {
        for (++p; p < (int)bit.size(); p += p & -p)
            bit[p] = (i64)((i128)bit[p] + delta);
    }

    i64 prefix(int r) const {
        i64 s = 0;
        for (; r > 0; r -= r & -r) s += bit[r];
        return s;
    }
};

static vector<i64> computeTwoDim(const vector<i64>& c, const vector<i64>& w, i64 k) {
    // L1 in (x,y) becomes L-infinity in (u=x+y,v=x-y).
    // Sweep centers in u order and maintain active points in a v-Fenwick tree.
    struct Point { i128 u, v; i64 weight; int id, vrank; };
    const int m = (int)w.size();
    vector<Point> points(m);
    vector<i128> uniqueV(m);
    for (int i = 0; i < m; ++i) {
        i64 x = c[(size_t)i * 2], y = c[(size_t)i * 2 + 1];
        points[i] = {(i128)x + y, (i128)x - y, w[i], i, 0};
        uniqueV[i] = (i128)x - y;
    }
    sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
        return a.u < b.u;
    });
    sort(uniqueV.begin(), uniqueV.end());
    uniqueV.erase(unique(uniqueV.begin(), uniqueV.end()), uniqueV.end());
    for (Point& p : points)
        p.vrank = lower_bound(uniqueV.begin(), uniqueV.end(), p.v) - uniqueV.begin();

    Fenwick bit((int)uniqueV.size());
    vector<i64> answer(m);
    int added = 0, removed = 0;
    for (int i = 0; i < m; ++i) {
        i128 lowU = (i128)points[i].u - k, highU = (i128)points[i].u + k;
        while (added < m && (i128)points[added].u <= highU) {
            bit.add(points[added].vrank, points[added].weight);
            ++added;
        }
        while (removed < added && (i128)points[removed].u < lowU) {
            bit.add(points[removed].vrank, -(i128)points[removed].weight);
            ++removed;
        }
        i128 lowV = (i128)points[i].v - k, highV = (i128)points[i].v + k;
        int l = lower_bound(uniqueV.begin(), uniqueV.end(), lowV) - uniqueV.begin();
        int r = upper_bound(uniqueV.begin(), uniqueV.end(), highV) - uniqueV.begin();
        answer[points[i].id] = bit.prefix(r) - bit.prefix(l);
    }
    return answer;
}

// Computes min(limit, C(n,r)) without overflow.
static uint64_t binomLimited(i64 n, i64 r, uint64_t limit) {
    if (r < 0 || r > n) return 0;
    r = min(r, n - r);
    u128 value = 1;
    for (i64 i = 1; i <= r; ++i) {
        value = value * (n - r + i) / i;
        if (value >= limit) return limit;
    }
    return (uint64_t)value;
}

// Number of integer lattice sites within L1 distance K, capped at 'limit'.
static uint64_t ballSizeLimited(i64 n, i64 k, uint64_t limit) {
    if (limit == 0) return 0;
    uint64_t answer = 1;
    for (i64 j = 1; j <= min(n, k) && answer < limit; ++j) {
        if (j >= 63) return limit;
        uint64_t a = binomLimited(n, j, limit);
        uint64_t b = binomLimited(k, j, limit);
        u128 term = (u128)a * b * ((uint64_t)1 << j);
        if (term >= limit - answer) return limit;
        answer += (uint64_t)term;
    }
    return answer;
}

// Closed-address-free table with linear probing. Entries are never erased,
// and the capacity is at least twice the number of inserted distinct keys.
// A value of -1 means an empty slot; keys may take ANY uint64_t value.
class FlatHash {
    vector<uint64_t> keys;
    vector<int> heads;
    size_t mask;

    static uint64_t mix(uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    size_t position(uint64_t key) const {
        size_t p = mix(key) & mask;
        while (heads[p] >= 0 && keys[p] != key) p = (p + 1) & mask;
        return p;
    }

public:
    explicit FlatHash(size_t maxEntries) {
        size_t size = 2;
        while (size < maxEntries * 2) size *= 2;
        mask = size - 1;
        keys.resize(size);
        heads.assign(size, -1);
    }

    int get(uint64_t key) const {
        return heads[position(key)];
    }

    int& getOrInsert(uint64_t key) {
        size_t p = position(key);
        keys[p] = key;
        return heads[p];
    }

    void clear() {
        fill(heads.begin(), heads.end(), -1);
    }
};

// Bounding-box tree for few varying coordinates; metric-ball tree for many.
// Its worst case is still quadratic for all queries, but a whole subtree is
// skipped or summed in O(1) once distance bounds settle the question.
class StaticSpatialTree {
    struct Node {
        int center, left, right;
        i64 sum;
        u128 radius;
    };

    const vector<i64>& coords;
    const vector<int>& representatives;
    const vector<i64>& weights;
    const vector<int>& dims;
    int n, count, root = -1;
    bool boxes;
    vector<int> order;
    vector<Node> nodes;
    vector<i64> boxLow, boxHigh;
    vector<u128> pivotDistances;

    i64 coord(int group, int dim) const {
        return coords[(size_t)representatives[group] * n + dim];
    }

    u128 distance(int a, int b) const {
        u128 sum = 0;
        for (int dim : dims) sum += absdiff(coord(a, dim), coord(b, dim));
        return sum;
    }

    int buildBoxes(int l, int r) {
        if (l == r) return -1;
        const int dimensions = (int)dims.size();
        vector<i64> lo(dimensions), hi(dimensions);
        for (int j = 0; j < dimensions; ++j)
            lo[j] = hi[j] = coord(order[l], dims[j]);
        for (int p = l + 1; p < r; ++p)
            for (int j = 0; j < dimensions; ++j) {
                i64 value = coord(order[p], dims[j]);
                lo[j] = min(lo[j], value);
                hi[j] = max(hi[j], value);
            }

        int axis = 0;
        for (int j = 1; j < dimensions; ++j)
            if (absdiff(lo[j], hi[j]) > absdiff(lo[axis], hi[axis])) axis = j;

        int mid = l + (r - l) / 2;
        nth_element(order.begin() + l, order.begin() + mid, order.begin() + r,
                    [&](int a, int b) {
                        i64 x = coord(a, dims[axis]), y = coord(b, dims[axis]);
                        return x == y ? a < b : x < y;
                    });

        int index = (int)nodes.size();
        nodes.push_back({order[mid], -1, -1, 0, 0});
        for (int j = 0; j < dimensions; ++j) {
            boxLow[(size_t)index * dimensions + j] = lo[j];
            boxHigh[(size_t)index * dimensions + j] = hi[j];
        }
        int left = buildBoxes(l, mid);
        int right = buildBoxes(mid + 1, r);
        nodes[index].left = left;
        nodes[index].right = right;
        nodes[index].sum = weights[nodes[index].center]
                           + (left < 0 ? 0 : nodes[left].sum)
                           + (right < 0 ? 0 : nodes[right].sum);
        return index;
    }

    int buildBalls(int l, int r) {
        if (l == r) return -1;
        int center = order[l];
        int index = (int)nodes.size();
        nodes.push_back({center, -1, -1, 0, 0});

        u128 radius = 0;
        for (int p = l + 1; p < r; ++p) {
            int g = order[p];
            pivotDistances[g] = distance(center, g);
            radius = max(radius, pivotDistances[g]);
        }
        int left = -1, right = -1;
        if (r - l > 1) {
            int mid = l + 1 + (r - l - 1) / 2;
            nth_element(order.begin() + l + 1, order.begin() + mid, order.begin() + r,
                        [&](int a, int b) {
                            return pivotDistances[a] == pivotDistances[b]
                                   ? a < b : pivotDistances[a] < pivotDistances[b];
                        });
            left = buildBalls(l + 1, mid);
            right = buildBalls(mid, r);
        }
        nodes[index].left = left;
        nodes[index].right = right;
        nodes[index].radius = radius;
        nodes[index].sum = weights[center]
                           + (left < 0 ? 0 : nodes[left].sum)
                           + (right < 0 ? 0 : nodes[right].sum);
        return index;
    }

    i64 queryNode(int index, int center, u128 k) const {
        if (index < 0) return 0;
        const Node& node = nodes[index];
        u128 d;
        if (boxes) {
            u128 minimum = 0;
            for (int j = 0; j < (int)dims.size(); ++j) {
                i64 x = coord(center, dims[j]);
                i64 lo = boxLow[(size_t)index * dims.size() + j];
                i64 hi = boxHigh[(size_t)index * dims.size() + j];
                if (x < lo) minimum += absdiff(x, lo);
                else if (x > hi) minimum += absdiff(x, hi);
                if (minimum > k) return 0;
            }
            u128 maximum = 0;
            for (int j = 0; j < (int)dims.size(); ++j) {
                i64 x = coord(center, dims[j]);
                i64 lo = boxLow[(size_t)index * dims.size() + j];
                i64 hi = boxHigh[(size_t)index * dims.size() + j];
                maximum += max(absdiff(x, lo), absdiff(x, hi));
                if (maximum > k) break;
            }
            if (maximum <= k) return node.sum;
            d = distance(center, node.center);
        } else {
            d = distance(center, node.center);
            if (d > k + node.radius) return 0;
            if (d + node.radius <= k) return node.sum;
        }
        return (d <= k ? weights[node.center] : 0)
               + queryNode(node.left, center, k)
               + queryNode(node.right, center, k);
    }

public:
    StaticSpatialTree(const vector<i64>& c, const vector<int>& reps,
                      const vector<i64>& w, const vector<int>& active, int dimension)
        : coords(c), representatives(reps), weights(w), dims(active),
          n(dimension), count((int)reps.size()), boxes(active.size() <= 6) {
        order.resize(count);
        iota(order.begin(), order.end(), 0);
        nodes.reserve(count);
        if (boxes) {
            boxLow.resize((size_t)count * dims.size());
            boxHigh.resize((size_t)count * dims.size());
            root = buildBoxes(0, count);
        } else {
            pivotDistances.resize(count);
            root = buildBalls(0, count);
        }
    }

    vector<i64> computeAll(i64 k) const {
        vector<i64> result(count);
        for (int g = 0; g < count; ++g) result[g] = queryNode(root, g, (u128)k);
        return result;
    }
};

class SparseHighDim {
    struct Change { int dim; i64 delta; };
    struct Offset { i128 deltaCode; uint32_t begin, count; };

    int n, m;
    i64 k;
    vector<i64> coords;
    vector<int> representatives, originalGroup, nextCollision;
    vector<i64> groupWeights;
    vector<uint64_t> groupHashes;
    FlatHash hashHeads;
    vector<int> activeDims;
    vector<i64> low, high;
    vector<uint64_t> coefficient, randomCoefficient;
    bool packed = false;
    uint64_t salt;

    i64 coordinate(int id, int d) const { return coords[(size_t)id * n + d]; }

    static uint64_t mix(uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    uint64_t hashPoint(int id) const {
        uint64_t h = 0;
        for (int d = 0; d < n; ++d)
            h += coefficient[d] * (uint64_t)coordinate(id, d);
        return h;
    }

    bool equalPoints(int a, int b) const {
        for (int d = 0; d < n; ++d)
            if (coordinate(a, d) != coordinate(b, d)) return false;
        return true;
    }

    // For a packed key, coordinates are encoded injectively in a 64-bit mixed
    // radix integer. Otherwise use a linear randomized hash, verifying each
    // match exactly, so even an actual hash collision cannot change an answer.
    bool matchesOffset(int source, int target, const Offset& off,
                       const vector<Change>& changes) const {
        int nextChange = off.begin, end = off.begin + off.count;
        for (int d = 0; d < n; ++d) {
            i128 expected = coordinate(representatives[source], d);
            if (nextChange < end && changes[nextChange].dim == d)
                expected += changes[nextChange++].delta;
            if (expected != coordinate(representatives[target], d)) return false;
        }
        return true;
    }

    void computeWithOffsets(vector<i64>& answers) const {
        vector<Offset> offsets;
        vector<Change> changes, path;
        auto generate = [&](auto&& self, int depth, i64 remain, bool nonzero,
                            i128 deltaCode) -> void {
            if (depth == (int)activeDims.size()) {
                if (nonzero) {
                    offsets.push_back(Offset{deltaCode, (uint32_t)changes.size(),
                                             (uint32_t)path.size()});
                    changes.insert(changes.end(), path.begin(), path.end());
                }
                return;
            }
            int dim = activeDims[depth];
            for (i64 delta = -remain; delta <= remain; ++delta) {
                if (!nonzero && delta < 0) continue;
                if (delta) path.push_back(Change{dim, delta});
                self(self, depth + 1, remain - (delta < 0 ? -delta : delta),
                     nonzero || delta != 0,
                     deltaCode + (i128)coefficient[dim] * delta);
                if (delta) path.pop_back();
            }
        };
        generate(generate, 0, k, false, 0);

        for (int source = 0; source < (int)representatives.size(); ++source) {
            uint64_t base = groupHashes[source];
            for (const Offset& off : offsets) {
                // Packed keys have K empty slots on either side of every
                // coordinate range, so an offset cannot carry to another digit.
                uint64_t key = base + (uint64_t)off.deltaCode;
                int head = hashHeads.get(key);
                if (head < 0) continue;
                if (packed) {
                    int target = head;
                    answers[source] += groupWeights[target];
                    answers[target] += groupWeights[source];
                } else {
                    for (int target = head; target >= 0;
                         target = nextCollision[target]) {
                        if (matchesOffset(source, target, off, changes)) {
                            answers[source] += groupWeights[target];
                            answers[target] += groupWeights[source];
                        }
                    }
                }
            }
        }
    }

    // If K < the number of varying dimensions, partition those dimensions
    // into K+1 nonempty blocks. Any pair at L1 distance <= K differs in at
    // most K coordinates, so it is identical on at least one whole block.
    // Enumerate pairs sharing a block, counting each pair only at its first
    // identical block. Preflight rejects skewed buckets before quadratic work.
    bool tryProjectionPairs(vector<i64>& answers) const {
        if (k >= (i64)activeDims.size()) return false;
        const int u = (int)representatives.size();
        const int numBlocks = (int)k + 1;
        vector<int> dimensions = activeDims;
        sort(dimensions.begin(), dimensions.end(), [&](int a, int b) {
            return absdiff(high[a], low[a]) > absdiff(high[b], low[b]);
        });
        vector<vector<int>> blocks(numBlocks);
        for (int i = 0; i < (int)dimensions.size(); ++i)
            blocks[i % numBlocks].push_back(dimensions[i]);

        auto hashBlock = [&](int g, int block) -> uint64_t {
            uint64_t h = 0;
            int id = representatives[g];
            for (int d : blocks[block])
                h += randomCoefficient[d] * (uint64_t)coordinate(id, d);
            return h;
        };
        auto equalBlock = [&](int a, int b, int block) -> bool {
            for (int d : blocks[block])
                if (coordinate(representatives[a], d)
                    != coordinate(representatives[b], d)) return false;
            return true;
        };
        auto nearby = [&](int a, int b) -> bool {
            u128 distance = 0;
            for (int d : dimensions) {
                distance += absdiff(coordinate(representatives[a], d),
                                    coordinate(representatives[b], d));
                if (distance > (u128)k) return false;
            }
            return true;
        };

        FlatHash table(u);
        vector<int> previous(u), chainLength(u);
        // If candidate pairs are too numerous, geometric pruning is cheaper.
        const uint64_t budget = min((uint64_t)u * 8, (uint64_t)12000000);
        uint64_t candidates = 0;
        for (int block = 0; block < numBlocks; ++block) {
            table.clear();
            for (int g = 0; g < u; ++g) {
                int& head = table.getOrInsert(hashBlock(g, block));
                chainLength[g] = (head < 0 ? 1 : chainLength[head] + 1);
                candidates += (uint64_t)chainLength[g] - 1;
                if (candidates > budget) return false;
                head = g;
            }
        }

        for (int block = 0; block < numBlocks; ++block) {
            table.clear();
            for (int g = 0; g < u; ++g) {
                int& head = table.getOrInsert(hashBlock(g, block));
                previous[g] = head;
                head = g;
                for (int p = previous[g]; p >= 0; p = previous[p]) {
                    if (!nearby(g, p) || !equalBlock(g, p, block)) continue;
                    bool seen = false;
                    for (int earlier = 0; earlier < block; ++earlier)
                        if (equalBlock(g, p, earlier)) {
                            seen = true;
                            break;
                        }
                    if (!seen) {
                        answers[g] += groupWeights[p];
                        answers[p] += groupWeights[g];
                    }
                }
            }
        }
        return true;
    }

public:
    SparseHighDim(int dimension, i64 radius, vector<i64> c, vector<i64> w)
        : n(dimension), m((int)w.size()), k(radius), coords(move(c)),
          hashHeads(w.size()),
          salt((uint64_t)chrono::steady_clock::now().time_since_epoch().count()) {
        representatives.reserve(m);
        groupWeights.reserve(m);
        groupHashes.reserve(m);
        nextCollision.reserve(m);
        originalGroup.resize(m);
        coefficient.resize(n);
        randomCoefficient.resize(n);
        for (int d = 0; d < n; ++d)
            randomCoefficient[d] = coefficient[d] = mix(salt + (uint64_t)d) | 1ULL;
        for (int id = 0; id < m; ++id) {
            uint64_t h = hashPoint(id);
            int& head = hashHeads.getOrInsert(h);
            int group = -1;
            if (head >= 0) {
                for (int g = head; g >= 0; g = nextCollision[g]) {
                    if (equalPoints(representatives[g], id)) {
                        group = g;
                        break;
                    }
                }
            }
            if (group >= 0) {
                groupWeights[group] += w[id];
            } else {
                group = (int)representatives.size();
                representatives.push_back(id);
                groupWeights.push_back(w[id]);
                groupHashes.push_back(h);
                nextCollision.push_back(head);
                head = group;
            }
            originalGroup[id] = group;
        }

        if (!representatives.empty()) {
            low.resize(n);
            high.resize(n);
            for (int d = 0; d < n; ++d)
                low[d] = high[d] = coordinate(representatives[0], d);
            for (int id : representatives) {
                for (int d = 0; d < n; ++d) {
                    low[d] = min(low[d], coordinate(id, d));
                    high[d] = max(high[d], coordinate(id, d));
                }
            }
            for (int d = 0; d < n; ++d)
                if (low[d] != high[d]) activeDims.push_back(d);

            // Pad each active coordinate range by K on BOTH sides. This lets
            // us add any L1 offset of radius <=K directly to the packed key,
            // with no carry or borrow between digits and no per-point checks.
            const u128 capacity = (u128)ULLONG_MAX + 1;
            u128 product = 1;
            packed = true;
            for (int d : activeDims) {
                u128 span = absdiff(high[d], low[d]) + 1 + (u128)2 * (u128)k;
                if (span > capacity / product) {
                    packed = false;
                    break;
                }
                coefficient[d] = (uint64_t)product;
                product *= span;
            }
            if (packed) {
                hashHeads.clear();
                for (int g = 0; g < (int)representatives.size(); ++g) {
                    u128 key = 0;
                    for (int d : activeDims)
                        key += (u128)((i128)coordinate(representatives[g], d)
                                      - (i128)low[d] + (i128)k) * coefficient[d];
                    groupHashes[g] = (uint64_t)key;
                    nextCollision[g] = -1;
                    hashHeads.getOrInsert(groupHashes[g]) = g;
                }
            } else {
                // Restore the original random coefficients: they were modified
                // while attempting the mixed-radix encoding.
                coefficient = randomCoefficient;
            }
        }
    }

    vector<i64> computeAll() const {
        const int u = (int)representatives.size();
        vector<i64> answers = groupWeights;
        if (u > 1 && k > 0) {
            // If the L1 diameter upper bound is <= k, every point is nearby.
            u128 diameterBound = 0;
            for (int d = 0; d < n; ++d) diameterBound += absdiff(low[d], high[d]);
            // A half-stencil visits (V-1)/2 other positions per source.
            // Keep that work below roughly 80 million hash probes; also
            // avoid generating more offsets than there are unique points.
            uint64_t stencilLimit = min((uint64_t)u,
                max((uint64_t)128, (uint64_t)160000000 / (uint64_t)u + 1));
            if (diameterBound <= (u128)k) {
                i64 total = accumulate(groupWeights.begin(), groupWeights.end(), 0LL);
                fill(answers.begin(), answers.end(), total);
            } else if (ballSizeLimited((int)activeDims.size(), k, stencilLimit)
                       < stencilLimit) {
                // Compute all offset hashes just once, then reuse for every point.
                computeWithOffsets(answers);
            } else if (!tryProjectionPairs(answers)) {
                // Large stencil: use geometric bounds and subtree weight sums.
                StaticSpatialTree tree(coords, representatives, groupWeights, activeDims, n);
                answers = tree.computeAll(k);
            }
        }
        vector<i64> output(m);
        for (int id = 0; id < m; ++id) output[id] = answers[originalGroup[id]];
        return output;
    }
};

// Reusable batch API. The result has exactly M elements in the input order.
// The flattened coordinate array is [p0_dim0, ..., p0_dimN-1, p1_dim0, ...].
static vector<i64> computeAllWeightedManhattan(
    int n, i64 k, vector<i64> coords, vector<i64> weights) {
    if (n == 1) return computeOneDim(coords, weights, k);
    if (n == 2) return computeTwoDim(coords, weights, k);
    SparseHighDim index(n, k, move(coords), move(weights));
    return index.computeAll();
}

/*
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, m;
    i64 k;
    if (!(cin >> n >> m >> k)) return 0;
    if (n <= 0 || m < 0 || k < 0) return 1;
    vector<i64> coords((size_t)n * m), weights(m);
    for (int i = 0; i < m; ++i) {
        for (int d = 0; d < n; ++d) cin >> coords[(size_t)i * n + d];
        cin >> weights[i];
    }

    for (i64 answer : computeAllWeightedManhattan(n, k, move(coords), move(weights)))
        cout << answer << '\n';
}
*/
