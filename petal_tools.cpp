/**
 * @file petal_tools.cpp
 * @brief Command-line tool for analyzing petal permutations.
 *
 * This program provides functionality to work with petal permutations,
 * which are permutations of odd length. It can compute various properties
 * such as inverses, reverses, symmetry orbits, canonical representatives,
 * and cyclic descents.
 */

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * @class PetalPermutation
 * @brief Represents a petal permutation - a permutation of odd length.
 *
 * A petal permutation is a permutation of integers from 1 to N where N is odd.
 * This class provides methods to manipulate and analyze such permutations,
 * including symmetry operations and canonical forms.
 */
class PetalPermutation {
public:
    /**
     * @brief Constructs a PetalPermutation from a vector of integers.
     * @param values The permutation values (must be a valid odd-length permutation of 1..N).
     * @throws std::invalid_argument if the input is not a valid petal permutation.
     */
    explicit PetalPermutation(std::vector<int> values)
        : perm_(std::move(values)) {
        validate();
    }

    /**
     * @brief Returns the underlying permutation values.
     * @return Const reference to the vector of permutation values.
     */
    const std::vector<int>& values() const {
        return perm_;
    }

    /**
     * @brief Returns the size of the permutation.
     * @return The number of elements in the permutation.
     */
    int size() const {
        return static_cast<int>(perm_.size());
    }

    /**
     * @brief Converts the permutation to a string representation.
     * @return String representation in the format "{ a, b, c, ... }".
     */
    std::string toString() const {
        std::ostringstream out;
        out << "{ ";
        for (std::size_t i = 0; i < perm_.size(); ++i) {
            out << perm_[i];
            if (i + 1 < perm_.size()) {
                out << ", ";
            }
        }
        out << " }";
        return out.str();
    }

    /**
     * @brief Returns the reverse of this permutation.
     * @return A new PetalPermutation that is the reverse of this one.
     */
    PetalPermutation reversed() const {
        std::vector<int> r = perm_;
        std::reverse(r.begin(), r.end());
        return PetalPermutation(r);
    }

    /**
     * @brief Returns the inverse permutation.
     * @return A new PetalPermutation that is the inverse of this one.
     */
    PetalPermutation inverse() const {
        const int n = size();
        std::vector<int> inv(n, 0);

        for (int i = 0; i < n; ++i) {
            const int v = perm_[i];
            inv[v - 1] = i + 1;
        }

        return PetalPermutation(inv);
    }

    /**
     * @brief Rotates the permutation left by the specified number of positions.
     * @param shift The number of positions to rotate left (can be negative for right rotation).
     * @return A new PetalPermutation that is the rotated version.
     */
    PetalPermutation rotatedLeft(int shift) const {
        const int n = size();

        shift %= n;
        if (shift < 0) {
            shift += n;
        }

        std::vector<int> r(n);
        for (int i = 0; i < n; ++i) {
            r[i] = perm_[(i + shift) % n];
        }

        return PetalPermutation(r);
    }

    /**
     * @brief Generates all cyclic shifts of this permutation.
     * @return A vector containing all cyclic shifts of this permutation.
     */
    std::vector<PetalPermutation> cyclicShifts() const {
        std::vector<PetalPermutation> result;
        result.reserve(size());

        for (int s = 0; s < size(); ++s) {
            result.push_back(rotatedLeft(s));
        }

        return result;
    }

    /**
     * @brief Computes the symmetry orbit of this permutation.
     *
     * The symmetry orbit includes all permutations that can be obtained
     * by applying rotations, reversals, and inversions.
     * @return A vector of unique permutations in the symmetry orbit.
     */
    std::vector<PetalPermutation> symmetryOrbit() const {
        std::vector<PetalPermutation> seeds = {
            *this,
            reversed(),
            inverse(),
            inverse().reversed()
        };

        std::vector<PetalPermutation> all;
        for (const auto& p : seeds) {
            for (const auto& q : p.cyclicShifts()) {
                all.push_back(q);
            }
        }

        std::set<std::string> seen;
        std::vector<PetalPermutation> unique;

        for (const auto& p : all) {
            const std::string key = p.compactKey();
            if (seen.insert(key).second) {
                unique.push_back(p);
            }
        }

        return unique;
    }

    /**
     * @brief Returns the canonical representative of the symmetry orbit.
     *
     * The canonical representative is the lexicographically smallest
     * permutation in the symmetry orbit.
     * @return The canonical representative permutation.
     */
    PetalPermutation canonicalRepresentative() const {
        auto orbit = symmetryOrbit();
        PetalPermutation best = orbit.front();

        for (const auto& p : orbit) {
            if (p.compactKey() < best.compactKey()) {
                best = p;
            }
        }

        return best;
    }

    /**
     * @brief Checks if this permutation is the identity permutation.
     * @return True if this is the identity permutation (1, 2, 3, ..., N).
     */
    bool isIdentity() const {
        for (int i = 0; i < size(); ++i) {
            if (perm_[i] != i + 1) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Checks if this permutation is the reverse identity permutation.
     * @return True if this is the reverse identity (N, N-1, ..., 1).
     */
    bool isReverseIdentity() const {
        const int n = size();
        for (int i = 0; i < n; ++i) {
            if (perm_[i] != n - i) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Computes the cyclic descents of this permutation.
     *
     * A cyclic descent occurs at position i where perm[i] > perm[(i+1) % n].
     * @return A vector of positions where cyclic descents occur.
     */
    std::vector<int> descents() const {
        std::vector<int> result;
        const int n = size();

        for (int i = 0; i < n; ++i) {
            int current = perm_[i];
            int next = perm_[(i + 1) % n];

            if (current > next) {
                result.push_back(i);
            }
        }

        return result;
    }

    std::vector<int> ascents() const {
        std::vector<int> result;
        const int n = size();

        for (int i = 0; i < n; ++i) {
            int current = perm_[i];
            int next = perm_[(i + 1) % n];

            if (current < next) {
                result.push_back(i);
            }
        }

        return result;
    }

    std::vector<int> peaks() const {
    std::vector<int> result;
    const int n = size();

        for (int i = 0; i < n; ++i) {
            int previous = perm_[(i - 1 + n) % n];
            int current = perm_[i];
            int next = perm_[(i + 1) % n];

            if (current > previous && current > next) {
                result.push_back(i);
            }
        }

        return result;
    }

    std::vector<int> valleys() const {
        std::vector<int> result;
        const int n = size();

        for (int i = 0; i < n; ++i) {
            int previous = perm_[(i - 1 + n) % n];
            int current = perm_[i];
            int next = perm_[(i + 1) % n];

            if (current < previous && current < next) {
                result.push_back(i);
            }
        }

        return result;
    }

private:
    std::vector<int> perm_;

    /**
     * @brief Validates that the stored permutation is a valid petal permutation.
     * @throws std::invalid_argument if validation fails.
     */
    void validate() const {
        const int n = static_cast<int>(perm_.size());

        if (n == 0) {
            throw std::invalid_argument("Permutation must not be empty.");
        }

        if (n % 2 == 0) {
            throw std::invalid_argument("Petal permutation length must be odd.");
        }

        std::vector<int> expected(n);
        std::iota(expected.begin(), expected.end(), 1);

        std::vector<int> sorted = perm_;
        std::sort(sorted.begin(), sorted.end());

        if (sorted != expected) {
            throw std::invalid_argument("Input must be a permutation of 1..N.");
        }
    }

    /**
     * @brief Generates a compact string key for this permutation.
     * @return A string representation used for uniqueness checking.
     */
    std::string compactKey() const {
        std::ostringstream out;
        for (std::size_t i = 0; i < perm_.size(); ++i) {
            if (i > 0) {
                out << '-';
            }
            out << perm_[i];
        }
        return out.str();
    }
};

void printUsage(const std::string& programName) {
    std::cout
        << "Usage:\n"
        << "  " << programName << " [options] <odd-length permutation>\n\n"
        << "Examples:\n"
        << "  " << programName << " 1 5 3 7 2 4 6\n"
        << "  " << programName << " --orbit 1 5 3 7 2 4 6\n"
        << "  " << programName << " --canonical 1 5 3 7 2 4 6\n\n"
        << "Options:\n"
        << "  --orbit       Print the full symmetry orbit\n"
        << "  --canonical   Print only the canonical representative\n"
        << "  --help        Show this help message\n";
}

std::vector<int> parsePermutation(
    int argc,
    char* argv[],
    bool& printOrbit,
    bool& canonicalOnly,
    bool& orbitLines
) {
    std::vector<int> perm;
    printOrbit = false;
    canonicalOnly = false;
    orbitLines = false;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];

        if (arg == "--help") {
            printUsage(argv[0]);
            std::exit(0);
        }
        else if (arg == "--orbit") {
            printOrbit = true;
        }
        else if (arg == "--canonical") {
            canonicalOnly = true;
        }
        else if (arg == "--orbit-lines") {
    orbitLines = true;
        }
        else {
            try {
                std::size_t pos = 0;
                int value = std::stoi(arg, &pos);

                if (pos != arg.size()) {
                    throw std::invalid_argument("Non-integer argument.");
                }

                perm.push_back(value);
            }
            catch (const std::exception&) {
                throw std::invalid_argument("Could not parse argument as integer: " + arg);
            }
        }
    }

    if (perm.empty()) {
        throw std::invalid_argument("No permutation entries provided. Use --help for usage.");
    }

    return perm;
}

/**
 * @brief Prints a permutation as a single line of space-separated values.
 * @param p The permutation to print.
 */
void printPermutationLine(const PetalPermutation& p) {
    const auto& values = p.values();

    for (std::size_t i = 0; i < values.size(); ++i) {
        std::cout << values[i];
        if (i + 1 < values.size()) {
            std::cout << ' ';
        }
    }

    std::cout << '\n';
}

/**
 * @brief Main entry point for the petal_tools command-line application.
 *
 * Parses command-line arguments, processes the input permutation,
 * and outputs various properties based on the specified options.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @return 0 on success, 1 on error.
 */
int main(int argc, char* argv[]) {
    try {
        bool printOrbit = false;
        bool canonicalOnly = false;
        bool orbitLines = false;

        std::vector<int> values = parsePermutation(argc, argv, printOrbit, canonicalOnly, orbitLines);
        PetalPermutation p(values);

        if (orbitLines) {
            auto orbit = p.symmetryOrbit();

            for (const auto& q : orbit) {
                printPermutationLine(q);
            }

            return 0;
        }

        if (canonicalOnly) {
            std::cout << p.canonicalRepresentative().toString() << '\n';
            return 0;
        }

        std::cout << "Input permutation:      " << p.toString() << '\n';
        std::cout << "Reverse:                " << p.reversed().toString() << '\n';
        std::cout << "Inverse:                " << p.inverse().toString() << '\n';
        std::cout << "Canonical representative: "
                  << p.canonicalRepresentative().toString() << '\n';
        std::cout << "Identity permutation:    "
                  << (p.isIdentity() ? "yes" : "no") << '\n';

        std::cout << "Reverse identity:        "
                  << (p.isReverseIdentity() ? "yes" : "no") << '\n';


        auto orbit = p.symmetryOrbit();
        std::cout << "Symmetry orbit size:    " << orbit.size() << '\n';

        auto d = p.descents();
        std::cout << "Cyclic descents:         ";
        for (int x : d) {
            std::cout << x << ' ';
        }
        std::cout << '\n';

        auto a = p.ascents();
        std::cout << "Cyclic ascents:          ";
        for (int x : a) {
            std::cout << x << ' ';
        }
        std::cout << '\n';  

        auto pk = p.peaks();
        std::cout << "Peaks:                   ";
        for (int x : pk) {
            std::cout << x << ' ';
        }
        std::cout << '\n';

        auto vl = p.valleys();
        std::cout << "Valleys:                 ";
        for (int x : vl) {
            std::cout << x << ' ';
        }
        std::cout << '\n';

        if (printOrbit) {
            std::cout << "\nSymmetry orbit:\n";
            for (std::size_t i = 0; i < orbit.size(); ++i) {
                std::cout << "  [" << i << "] " << orbit[i].toString() << '\n';
            }
        }

        return 0;
        }
        catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << '\n';
            std::cerr << "Use --help for usage.\n";
            return 1;
        }

}