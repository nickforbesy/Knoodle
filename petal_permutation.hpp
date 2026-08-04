#pragma once

#include <algorithm>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * A petal permutation: a permutation of integers 1..N where N is odd.
 * Provides symmetry orbit, canonical form, and combinatorial statistics.
 */
class PetalPermutation {
public:
    explicit PetalPermutation(std::vector<int> values)
        : perm_(std::move(values)) {
        validate();
    }

    const std::vector<int>& values() const { return perm_; }
    int size() const { return static_cast<int>(perm_.size()); }

    std::string toString() const {
        std::ostringstream out;
        out << "{ ";
        for (std::size_t i = 0; i < perm_.size(); ++i) {
            out << perm_[i];
            if (i + 1 < perm_.size()) out << ", ";
        }
        out << " }";
        return out.str();
    }

    PetalPermutation reversed() const {
        std::vector<int> r = perm_;
        std::reverse(r.begin(), r.end());
        return PetalPermutation(r);
    }

    PetalPermutation inverse() const {
        const int n = size();
        std::vector<int> inv(n, 0);
        for (int i = 0; i < n; ++i)
            inv[perm_[i] - 1] = i + 1;
        return PetalPermutation(inv);
    }

    PetalPermutation rotatedLeft(int shift) const {
        const int n = size();
        shift %= n;
        if (shift < 0) shift += n;
        std::vector<int> r(n);
        for (int i = 0; i < n; ++i)
            r[i] = perm_[(i + shift) % n];
        return PetalPermutation(r);
    }

    std::vector<PetalPermutation> cyclicShifts() const {
        std::vector<PetalPermutation> result;
        result.reserve(size());
        for (int s = 0; s < size(); ++s)
            result.push_back(rotatedLeft(s));
        return result;
    }

    std::vector<PetalPermutation> symmetryOrbit() const {
        std::vector<PetalPermutation> seeds = {
            *this, reversed(), inverse(), inverse().reversed()
        };
        std::vector<PetalPermutation> all;
        for (const auto& p : seeds)
            for (const auto& q : p.cyclicShifts())
                all.push_back(q);

        std::set<std::string> seen;
        std::vector<PetalPermutation> unique;
        for (const auto& p : all)
            if (seen.insert(p.compactKey()).second)
                unique.push_back(p);
        return unique;
    }

    PetalPermutation canonicalRepresentative() const {
        auto orbit = symmetryOrbit();
        PetalPermutation best = orbit.front();
        for (const auto& p : orbit)
            if (p.compactKey() < best.compactKey())
                best = p;
        return best;
    }

    bool isIdentity() const {
        for (int i = 0; i < size(); ++i)
            if (perm_[i] != i + 1) return false;
        return true;
    }

    bool isReverseIdentity() const {
        const int n = size();
        for (int i = 0; i < n; ++i)
            if (perm_[i] != n - i) return false;
        return true;
    }

    std::vector<int> descents() const {
        std::vector<int> result;
        const int n = size();
        for (int i = 0; i < n; ++i)
            if (perm_[i] > perm_[(i + 1) % n])
                result.push_back(i);
        return result;
    }

    std::vector<int> ascents() const {
        std::vector<int> result;
        const int n = size();
        for (int i = 0; i < n; ++i)
            if (perm_[i] < perm_[(i + 1) % n])
                result.push_back(i);
        return result;
    }

    std::vector<int> peaks() const {
        std::vector<int> result;
        const int n = size();
        for (int i = 0; i < n; ++i) {
            int prev = perm_[(i - 1 + n) % n];
            int curr = perm_[i];
            int next = perm_[(i + 1) % n];
            if (curr > prev && curr > next)
                result.push_back(i);
        }
        return result;
    }

    std::vector<int> valleys() const {
        std::vector<int> result;
        const int n = size();
        for (int i = 0; i < n; ++i) {
            int prev = perm_[(i - 1 + n) % n];
            int curr = perm_[i];
            int next = perm_[(i + 1) % n];
            if (curr < prev && curr < next)
                result.push_back(i);
        }
        return result;
    }

    std::string compactKey() const {
        std::ostringstream out;
        for (std::size_t i = 0; i < perm_.size(); ++i) {
            if (i > 0) out << '-';
            out << perm_[i];
        }
        return out.str();
    }

private:
    std::vector<int> perm_;

    void validate() const {
        const int n = static_cast<int>(perm_.size());
        if (n == 0)
            throw std::invalid_argument("Permutation must not be empty.");
        if (n % 2 == 0)
            throw std::invalid_argument("Petal permutation length must be odd.");
        std::vector<int> expected(n);
        std::iota(expected.begin(), expected.end(), 1);
        std::vector<int> sorted = perm_;
        std::sort(sorted.begin(), sorted.end());
        if (sorted != expected)
            throw std::invalid_argument("Input must be a permutation of 1..N.");
    }
};
