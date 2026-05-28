/**
 * @file petal_tools.cpp
 * @brief Command-line tool for analyzing petal permutations.
 */

#include "petal_permutation.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>


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