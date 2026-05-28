#include "petal_permutation.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// -----------------------------------------------------------------------
// 3D point and I/O
// -----------------------------------------------------------------------
struct Point3D {
    double x, y, z;
};

static Point3D polarPoint(double radius, double angle, double z) {
    return { radius * std::cos(angle), radius * std::sin(angle), z };
}

static double heightFromValue(int value, int n, double heightScale) {
    const double center = (n + 1) / 2.0;
    return heightScale * (value - center);
}

static void writeTSV(const std::vector<Point3D>& points, const std::string& filename) {
    std::ofstream out(filename);
    if (!out)
        throw std::runtime_error("Could not open output file: " + filename);
    out << std::fixed << std::setprecision(10);
    for (const auto& p : points)
        out << p.x << '\t' << p.y << '\t' << p.z << '\n';
}

// -----------------------------------------------------------------------
// Build star polygon from a PetalPermutation.
//
// For each strand i the polygon visits 5 key points in order:
//   Q_in  — outer entry point at unperturbed angle alpha, z=0
//   P_in  — inner entry point at perturbed angle alpha_in, z=z_i
//   P_out — inner exit point at perturbed angle beta_out, z=z_i
//   Q_out — outer exit point at perturbed angle beta_out, z=0
//   A     — petal apex between Q_out and the next strand's Q_in, z~0
//
// The implicit edge A(i) -> Q_in(i+1) completes each petal arc.
// The polygon is closed by appending polygon.front() at the end.
// -----------------------------------------------------------------------
std::vector<Point3D> buildStarPolygon(
    const PetalPermutation& perm,
    double innerRadius,
    double outerRadius,
    double petalRadius,
    double heightScale
) {
    const int n = perm.size();
    const auto& vals = perm.values();
    const double twoPi       = 2.0 * M_PI;
    const double epsAngle    = 0.03;
    const double innerJitter = 0.04;
    const double outerJitter = 0.03;
    const double petalJitter = 0.03;
    const double apexLift    = 0.05;

    std::vector<Point3D> polygon;
    polygon.reserve(5 * n + 1);

    for (int i = 0; i < n; ++i) {
        const int    next       = (i + 1) % n;
        const double alpha      = twoPi * i / n;
        const double alpha_next = twoPi * next / n;
        const double z          = heightFromValue(vals[i], n, heightScale);

        // Per-strand angle perturbation to break degeneracy at crossings.
        const double eps_i    = epsAngle * (i - 0.5 * (n - 1));
        const double alpha_in = alpha + eps_i;
        const double beta_out = alpha + M_PI - eps_i;

        // Midpoint angle for the petal apex between beta_out and alpha_next.
        double delta = alpha_next - beta_out;
        while (delta <= -M_PI) delta += twoPi;
        while (delta >   M_PI) delta -= twoPi;
        const double gamma = beta_out + 0.5 * delta;

        // Perturbed radii to avoid exact coincidences at arc crossings.
        const double r_i   = innerRadius * (1.0 + innerJitter * std::cos(3.0 * alpha));
        const double R_in  = outerRadius * (1.0 + outerJitter * std::sin(2.0 * alpha));
        const double R_out = outerRadius * (1.0 + outerJitter * std::cos(2.0 * alpha));
        const double Rp    = petalRadius * (1.0 + petalJitter * std::sin(3.0 * gamma));

        // Q_in uses the unperturbed alpha so outer-circle entry points sit at
        // the canonical strand positions (matching Perm_2_Poly.cpp).
        const Point3D Q_in  = polarPoint(R_in,  alpha,    0.0);
        const Point3D P_in  = polarPoint(r_i,   alpha_in, z);
        const Point3D P_out = polarPoint(r_i,   beta_out, z);
        const Point3D Q_out = polarPoint(R_out, beta_out, 0.0);
        const Point3D A     = polarPoint(Rp,    gamma,    apexLift * std::sin(gamma));

        polygon.push_back(Q_in);
        polygon.push_back(P_in);
        polygon.push_back(P_out);
        polygon.push_back(Q_out);
        polygon.push_back(A);
    }

    polygon.push_back(polygon.front());

    return polygon;
}

// -----------------------------------------------------------------------
// CLI helpers
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

// -----------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------
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

        const double innerRadius = 0.35;
        const double outerRadius = 2.0;
        const double petalRadius = 3.5;
        const double heightScale = 0.8;

        std::vector<Point3D> polygon = buildStarPolygon(
            p, innerRadius, outerRadius, petalRadius, heightScale
        );

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
