// test_grid_moves.cpp -- regression test for the Cromwell moves in
// grid_diagram.hpp.
//
// A move that quietly changed the knot type would corrupt every arc index
// this pipeline reports, and would not show up as a crash. So: for every
// rotation-orbit rep of length n, name the knot three ways --
//   (a) straight from the petal grid,
//   (b) after three random stabilizations,
//   (c) after a full simplify() run,
// -- and require all three names to agree, chirality included.
//
//   ./test_grid_moves 7     # 720 grids, ~1 minute
//
// Requires knoodleidentify on PATH with KNOODLE_KLUT_DIR set.
#include "grid_diagram.hpp"
#include "knot_naming.hpp"
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>

static std::string name(const GridDiagram& g, const char* path) {
    std::ofstream o(path); g.writePolygonTSV(o); o.close();
    std::vector<knotname::Summand> s;
    return knotname::identifyTSV(path, s);
}

int main(int argc, char** argv) {
    int n = argc > 1 ? std::atoi(argv[1]) : 7;
    int bad = 0, total = 0;
    std::vector<int> tail(n - 1); std::iota(tail.begin(), tail.end(), 2);
    std::mt19937_64 rng(7);
    do {
        std::vector<int> perm{1}; perm.insert(perm.end(), tail.begin(), tail.end());
        GridDiagram g = GridDiagram::fromPetalPermutation(perm);
        const std::string before = name(g, "/tmp/tm_before.tsv");

        // exercise stabilize explicitly, then the full search
        GridDiagram h = g;
        for (int k = 0; k < 3; ++k) h.stabilize(rng() % h.size(), (rng() & 1u) != 0);
        const std::string stab = name(h, "/tmp/tm_stab.tsv");

        GridDiagram s = g;
        s.simplify(rng(), 20000);
        const std::string after = name(s, "/tmp/tm_after.tsv");

        ++total;
        if (before != stab || before != after) {
            ++bad;
            std::cout << "MISMATCH  perm=";
            for (int v : perm) std::cout << v << ' ';
            std::cout << " before=" << before << " stabilized=" << stab
                      << " simplified=" << after << "  (" << g.size()
                      << " -> " << s.size() << ")\n";
        }
    } while (std::next_permutation(tail.begin(), tail.end()));
    std::cout << "checked " << total << " grids, " << bad << " mismatches\n";
    return bad == 0 ? 0 : 1;
}
