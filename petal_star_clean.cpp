#include "petal_permutation.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// -----------------------------------------------------------------------
// Petal-star polygon (deterministic, no random perturbation).
//
// Idea:
//   * The curve passes through the centre n times.  The t-th passage is a
//     straight chord carrying height h_t = heightScale*(perm[t]-(n+1)/2),
//     so the vertical order of the strands at the centre IS the petal
//     permutation -- exactly what a petal knot encodes.
//   * A petal knot's n strands must all mutually cross, which means n
//     near-diameters through the centre.  Exact diameters are concurrent
//     (a degenerate n-fold crossing).  We resolve this the way you
//     described: rotate every chord's far endpoint by a fixed small angle
//     delta (a uniform "shift left").  The chords then become tangent to a
//     small central circle, so all C(n,2) pairs cross simply, forming a
//     clean star of crossings.
//   * Over/under at every crossing is decided purely by the constant
//     strand heights, i.e. purely by the permutation.  Because the shift
//     is uniform (height-independent), a permutation and its reverse /
//     mirror give the same diagram up to symmetry -- no asymmetry bug.
//
// Traversal order uses step s = (n+1)/2 so consecutive passages exit and
// re-enter at adjacent boundary points; the connecting petals are then
// thin and add no crossings.
// -----------------------------------------------------------------------
struct Point3D {
    double x, y, z;
};

static double heightFromValue(int value, int n, double heightScale) {
    const double center = (n + 1) / 2.0;
    return heightScale * (value - center);
}

static double wrapTo2Pi(double a) {
    const double twoPi = 2.0 * M_PI;
    while (a < 0.0)     a += twoPi;
    while (a >= twoPi)  a -= twoPi;
    return a;
}

static void writeTSV(const std::vector<Point3D>& points, const std::string& filename) {
    std::ofstream out(filename);
    if (!out)
        throw std::runtime_error("Could not open output file: " + filename);
    out << std::fixed << std::setprecision(10);
    for (const auto& p : points)
        out << p.x << '\t' << p.y << '\t' << p.z << '\n';
}

std::vector<Point3D> buildPetalStar(
    const PetalPermutation& perm,
    double radius,
    double petalRadius,
    double heightScale
) {
    const int n = perm.size();
    const auto& vals = perm.values();
    const int step = (n + 1) / 2;            // traversal step; coprime to odd n
    const double twoPi = 2.0 * M_PI;

    // Uniform "shift left": rotate each chord's far endpoint by delta so the
    // n chords miss the centre. Must be < pi/n so the connecting petals stay
    // thin (petal gap = pi/n - delta > 0) and all chords still mutually cross.
    const double delta = M_PI / (2.0 * n);

    std::vector<Point3D> polygon;
    polygon.reserve(3 * n + 1);

    for (int t = 0; t < n; ++t) {
        const int    p        = (t * step) % n;
        const int    pNext    = ((t + 1) * step) % n;
        const double a        = twoPi * p / n;
        const double aNext    = twoPi * pNext / n;

        const double z     = heightFromValue(vals[t], n, heightScale);
        const double zNext = heightFromValue(vals[(t + 1) % n], n, heightScale);

        const double entryAng = a;
        const double exitAng  = a + M_PI + delta;

        const Point3D entry = { radius * std::cos(entryAng),
                                radius * std::sin(entryAng), z };
        const Point3D exit  = { radius * std::cos(exitAng),
                                radius * std::sin(exitAng),  z };

        // Thin petal connecting this passage's exit to the next entry.
        const double gap    = wrapTo2Pi(aNext - exitAng);   // ~ pi/n - delta
        const double tipAng = exitAng + 0.5 * gap;
        const Point3D tip   = { petalRadius * std::cos(tipAng),
                                petalRadius * std::sin(tipAng),
                                0.5 * (z + zNext) };

        polygon.push_back(entry);
        polygon.push_back(exit);
        polygon.push_back(tip);
    }

    polygon.push_back(polygon.front());      // close the loop

    return polygon;
}

// -----------------------------------------------------------------------
// CLI
// -----------------------------------------------------------------------
static void printUsage(const std::string& prog) {
    std::cout
        << "Usage:\n"
        << "  " << prog << " [options] <odd-length permutation>\n\n"
        << "Examples:\n"
        << "  " << prog << " 2 4 1 5 3\n"
        << "  " << prog << " --canonical 2 4 1 5 3\n"
        << "  " << prog << " --info 2 4 1 5 3\n\n"
        << "Permutation transforms (applied before building polygon):\n"
        << "  --canonical       Use canonical representative of symmetry orbit\n"
        << "  --inverse         Use inverse permutation\n"
        << "  --reverse         Use reversed permutation\n\n"
        << "Output options:\n"
        << "  --info            Print permutation analysis; skip polygon output\n"
        << "  --output <file>   Write TSV to <file>  (default: petal_star.tsv)\n"
        << "  --help            Show this help message\n";
}

static void printInfo(const PetalPermutation& p) {
    auto printVec = [](const std::string& label, const std::vector<int>& v) {
        std::cout << label;
        for (int x : v) std::cout << x << ' ';
        std::cout << '\n';
    };

    std::cout << "Input permutation:        " << p.toString() << '\n';
    std::cout << "Reverse:                  " << p.reversed().toString() << '\n';
    std::cout << "Inverse:                  " << p.inverse().toString() << '\n';
    std::cout << "Canonical representative: " << p.canonicalRepresentative().toString() << '\n';
    std::cout << "Identity:                 " << (p.isIdentity()        ? "yes" : "no") << '\n';
    std::cout << "Reverse identity:         " << (p.isReverseIdentity() ? "yes" : "no") << '\n';
    std::cout << "Symmetry orbit size:      " << p.symmetryOrbit().size() << '\n';
    printVec("Cyclic descents:          ", p.descents());
    printVec("Cyclic ascents:           ", p.ascents());
    printVec("Peaks:                    ", p.peaks());
    printVec("Valleys:                  ", p.valleys());
}

int main(int argc, char* argv[]) {
    try {
        bool doCanonical = false;
        bool doInverse   = false;
        bool doReverse   = false;
        bool infoOnly    = false;
        std::string outputFile = "petal_star.tsv";
        std::vector<int> vals;

        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];

            if (arg == "--help") {
                printUsage(argv[0]);
                return 0;
            } else if (arg == "--canonical") {
                doCanonical = true;
            } else if (arg == "--inverse") {
                doInverse = true;
            } else if (arg == "--reverse") {
                doReverse = true;
            } else if (arg == "--info") {
                infoOnly = true;
            } else if (arg == "--output") {
                if (i + 1 >= argc)
                    throw std::invalid_argument("--output requires a filename argument.");
                outputFile = argv[++i];
            } else {
                std::size_t pos = 0;
                int value = std::stoi(arg, &pos);
                if (pos != arg.size())
                    throw std::invalid_argument("Could not parse argument as integer: " + arg);
                vals.push_back(value);
            }
        }

        if (vals.empty()) {
            printUsage(argv[0]);
            return 1;
        }

        PetalPermutation p(vals);

        if (doCanonical) p = p.canonicalRepresentative();
        if (doInverse)   p = p.inverse();
        if (doReverse)   p = p.reversed();

        if (infoOnly) {
            printInfo(p);
            return 0;
        }

        const double radius      = 1.0;
        const double petalRadius = 1.1;
        const double heightScale = 1.0;

        std::vector<Point3D> polygon = buildPetalStar(p, radius, petalRadius, heightScale);

        writeTSV(polygon, outputFile);

        std::cout << "Permutation: " << p.toString() << '\n';
        std::cout << "Wrote " << outputFile << " with "
                  << polygon.size() << " vertices.\n";

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        std::cerr << "Use --help for usage.\n";
        return 1;
    }
}
