#include <bits/stdc++.h>
using namespace std;

template <class T>
struct FastHashMap {
    struct Node {
        long long first = 0;
        T second = T();
    };

private:
    vector<Node> data;
    vector<unsigned char> used;
    size_t mask = 0;
    size_t sz = 0;
    size_t threshold = 0;

    static uint64_t splitmix64(uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    static uint64_t fixed_random() {
        static const uint64_t x =
            chrono::steady_clock::now().time_since_epoch().count();
        return x;
    }

    size_t hash_key(long long x) const {
        return splitmix64((uint64_t)x + fixed_random()) & mask;
    }

    void maybe_rehash() {
        if (sz >= threshold) rehash(data.size() << 1);
    }

    void insert_move(long long key, T&& val) {
        size_t i = hash_key(key);
        while (used[i]) {
            if (data[i].first == key) {
                data[i].second = move(val);
                return;
            }
            i = (i + 1) & mask;
        }
        used[i] = 1;
        data[i].first = key;
        data[i].second = move(val);
        ++sz;
    }

public:
    struct iterator {
        using iterator_category = forward_iterator_tag;
        using value_type = Node;
        using difference_type = ptrdiff_t;
        using pointer = Node*;
        using reference = Node&;

        FastHashMap* mp = nullptr;
        size_t idx = 0;

        iterator() = default;
        iterator(FastHashMap* mp_, size_t idx_) : mp(mp_), idx(idx_) {
            skip();
        }

        void skip() {
            while (idx < mp->data.size() && !mp->used[idx]) ++idx;
        }

        reference operator*() const { return mp->data[idx]; }
        pointer operator->() const { return &mp->data[idx]; }

        iterator& operator++() {
            ++idx;
            skip();
            return *this;
        }

        iterator operator++(int) {
            iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator==(const iterator& other) const {
            return mp == other.mp && idx == other.idx;
        }

        bool operator!=(const iterator& other) const {
            return !(*this == other);
        }
    };

    struct const_iterator {
        using iterator_category = forward_iterator_tag;
        using value_type = const Node;
        using difference_type = ptrdiff_t;
        using pointer = const Node*;
        using reference = const Node&;

        const FastHashMap* mp = nullptr;
        size_t idx = 0;

        const_iterator() = default;
        const_iterator(const FastHashMap* mp_, size_t idx_) : mp(mp_), idx(idx_) {
            skip();
        }

        void skip() {
            while (idx < mp->data.size() && !mp->used[idx]) ++idx;
        }

        reference operator*() const { return mp->data[idx]; }
        pointer operator->() const { return &mp->data[idx]; }

        const_iterator& operator++() {
            ++idx;
            skip();
            return *this;
        }

        const_iterator operator++(int) {
            const_iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator==(const const_iterator& other) const {
            return mp == other.mp && idx == other.idx;
        }

        bool operator!=(const const_iterator& other) const {
            return !(*this == other);
        }
    };

    FastHashMap(size_t n = 0) {
        reserve(max<size_t>(4, n));
    }

    void reserve(size_t n) {
        size_t cap = 1;
        while (cap < n * 2) cap <<= 1;
        rehash(cap);
    }

    void rehash(size_t new_cap) {
        vector<Node> old_data = move(data);
        vector<unsigned char> old_used = move(used);

        data.assign(new_cap, Node());
        used.assign(new_cap, 0);
        mask = new_cap - 1;
        sz = 0;
        threshold = new_cap * 7 / 10;

        for (size_t i = 0; i < old_data.size(); ++i) {
            if (old_used[i]) {
                insert_move(old_data[i].first, move(old_data[i].second));
            }
        }
    }

    T& operator[](long long key) {
        maybe_rehash();
        size_t i = hash_key(key);
        while (used[i]) {
            if (data[i].first == key) return data[i].second;
            i = (i + 1) & mask;
        }
        used[i] = 1;
        data[i].first = key;
        data[i].second = T();
        ++sz;
        return data[i].second;
    }

    bool contains(long long key) const {
        size_t i = hash_key(key);
        while (used[i]) {
            if (data[i].first == key) return true;
            i = (i + 1) & mask;
        }
        return false;
    }

    T* find_ptr(long long key) {
        size_t i = hash_key(key);
        while (used[i]) {
            if (data[i].first == key) return &data[i].second;
            i = (i + 1) & mask;
        }
        return nullptr;
    }

    const T* find_ptr(long long key) const {
        size_t i = hash_key(key);
        while (used[i]) {
            if (data[i].first == key) return &data[i].second;
            i = (i + 1) & mask;
        }
        return nullptr;
    }

    bool insert(long long key, const T& val) {
        maybe_rehash();
        size_t i = hash_key(key);
        while (used[i]) {
            if (data[i].first == key) return false;
            i = (i + 1) & mask;
        }
        used[i] = 1;
        data[i].first = key;
        data[i].second = val;
        ++sz;
        return true;
    }

    bool insert(long long key, T&& val) {
        maybe_rehash();
        size_t i = hash_key(key);
        while (used[i]) {
            if (data[i].first == key) return false;
            i = (i + 1) & mask;
        }
        used[i] = 1;
        data[i].first = key;
        data[i].second = move(val);
        ++sz;
        return true;
    }

    size_t size() const {
        return sz;
    }

    bool empty() const {
        return sz == 0;
    }

    void clear() {
        fill(used.begin(), used.end(), 0);
        sz = 0;
    }

    iterator begin() { return iterator(this, 0); }
    iterator end() { return iterator(this, data.size()); }

    const_iterator begin() const { return const_iterator(this, 0); }
    const_iterator end() const { return const_iterator(this, data.size()); }

    const_iterator cbegin() const { return const_iterator(this, 0); }
    const_iterator cend() const { return const_iterator(this, data.size()); }
};
