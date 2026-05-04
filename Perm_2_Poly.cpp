#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

// --------------------------------------------------
// Basic 3D point
// --------------------------------------------------
struct Point3D {
    double x;
    double y;
    double z;
};

// --------------------------------------------------
// Vector utilities
// --------------------------------------------------
Point3D lerp(const Point3D& a, const Point3D& b, double t) {
    return {
        (1.0 - t) * a.x + t * b.x,
        (1.0 - t) * a.y + t * b.y,
        (1.0 - t) * a.z + t * b.z
    };
}

Point3D polarPoint(double radius, double angle, double z) {
    return {
        radius * std::cos(angle),
        radius * std::sin(angle),
        z
    };
}

// --------------------------------------------------
// Angle wrapping into (-pi, pi]
// --------------------------------------------------
double wrapToPi(double angle) {
    const double twoPi = 2.0 * M_PI;

    while (angle <= -M_PI) {
        angle += twoPi;
    }
    while (angle > M_PI) {
        angle -= twoPi;
    }

    return angle;
}

// --------------------------------------------------
// Petal permutation validation:
// - odd length
// - values are exactly 1..N with no repeats
// --------------------------------------------------
bool isValidPetalPermutation(const std::vector<int>& perm) {
    const int n = static_cast<int>(perm.size());

    if (n == 0 || n % 2 == 0) {
        return false;
    }

    std::unordered_set<int> seen;
    for (int value : perm) {
        if (value < 1 || value > n) {
            return false;
        }
        if (seen.count(value) > 0) {
            return false;
        }
        seen.insert(value);
    }

    return true;
}

// --------------------------------------------------
// Convert permutation value to centered height
// z_i = h * (perm[i] - (N+1)/2)
// --------------------------------------------------
double permutationValueToHeight(int value, int n, double heightScale) {
    const double center = (n + 1) / 2.0;
    return heightScale * (value - center);
}

// --------------------------------------------------
// Append sampled points along a line segment.
// If skipFirst is true, do not repeat the first point.
// subdivisions = number of equal subsegments
// so subdivisions=1 adds only the endpoint
// --------------------------------------------------
void appendSegment(
    std::vector<Point3D>& polygon,
    const Point3D& a,
    const Point3D& b,
    int subdivisions,
    bool skipFirst
) {
    if (subdivisions < 1) {
        throw std::invalid_argument("subdivisions must be >= 1");
    }

    int start = skipFirst ? 1 : 0;
    for (int k = start; k <= subdivisions; ++k) {
        double t = static_cast<double>(k) / subdivisions;
        polygon.push_back(lerp(a, b, t));
    }
}

// --------------------------------------------------
// Write polygon to TSV
// --------------------------------------------------
void writePolygonTSV(const std::vector<Point3D>& polygon, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Could not open output file: " + filename);
    }

    out << std::fixed << std::setprecision(10);
    for (const auto& p : polygon) {
        out << p.x << '\t' << p.y << '\t' << p.z << '\n';
    }
}

// --------------------------------------------------
// Build a petal-style polygon embedding.
//
// Geometry:
// For each strand i:
//
// Q_in(i)   = outer entry point at angle alpha_i, z=0
// P_in(i)   = inner entry point at angle alpha_i, z=z_i
// P_out(i)  = inner exit point at angle beta_i = alpha_i + pi, z=z_i
// Q_out(i)  = outer exit point at angle beta_i, z=0
// A_i       = petal apex point at larger radius, between beta_i and alpha_{i+1}
//
// Then connect
// Q_in(i) -> P_in(i) -> P_out(i) -> Q_out(i) -> A_i -> Q_in(i+1)
//
// This gives one closed polygon.
// --------------------------------------------------
std::vector<Point3D> buildPetalPolygon(
    const std::vector<int>& perm,
    double innerRadius = 0.35,
    double outerRadius = 2.0,
    double petalRadius = 3.5,
    double heightScale = 0.8,
    int rampSubdivisions = 6,
    int chordSubdivisions = 16,
    int petalSubdivisions = 12
) {
    if (!isValidPetalPermutation(perm)) {
        throw std::invalid_argument("Invalid petal permutation.");
    }

    if (!(0.0 < innerRadius && innerRadius < outerRadius && outerRadius < petalRadius)) {
        throw std::invalid_argument("Require 0 < innerRadius < outerRadius < petalRadius.");
    }

    const int n = static_cast<int>(perm.size());
    const double twoPi = 2.0 * M_PI;

    // Small perturbation parameters to avoid exact geometric degeneracies
    const double epsAngle = 0.03;     // small angle offset in radians
    const double innerJitter = 0.04;  // relative inner radius variation
    const double outerJitter = 0.03;  // relative outer radius variation
    const double petalJitter = 0.03;  // relative petal radius variation
    const double apexLift = 0.05;     // tiny z-lift for outer apex points

    std::vector<Point3D> polygon;
    polygon.reserve(n * (2 * rampSubdivisions + chordSubdivisions + 2 * petalSubdivisions + 10));

    for (int i = 0; i < n; ++i) {
        const int next = (i + 1) % n;

        const double alpha_i = twoPi * i / n;
        const double alpha_next = twoPi * next / n;

        const double z_i = permutationValueToHeight(perm[i], n, heightScale);

        // Strand-dependent perturbation
        const double eps_i = epsAngle * (i - 0.5 * (n - 1));

        // Slightly perturbed entry/exit angles for the central strand
        const double alpha_in = alpha_i + eps_i;
        const double beta_out = alpha_i + M_PI - eps_i;

        // Slightly perturbed radii
        const double r_i = innerRadius * (1.0 + innerJitter * std::cos(3.0 * alpha_i));
        const double R_in_i = outerRadius * (1.0 + outerJitter * std::sin(2.0 * alpha_i));
        const double R_out_i = outerRadius * (1.0 + outerJitter * std::cos(2.0 * alpha_i));

        // Mid-angle for the petal arc from beta_out to alpha_next
        const double delta = wrapToPi(alpha_next - beta_out);
        const double gamma_i = beta_out + 0.5 * delta;

        const double Rp_i = petalRadius * (1.0 + petalJitter * std::sin(3.0 * gamma_i));
        const double z_apex = apexLift * std::sin(gamma_i);

        const Point3D Q_in   = polarPoint(R_in_i, alpha_i, 0.0);
        const Point3D P_in   = polarPoint(r_i, alpha_in, z_i);
        const Point3D P_out  = polarPoint(r_i, beta_out, z_i);
        const Point3D Q_out  = polarPoint(R_out_i, beta_out, 0.0);
        const Point3D A      = polarPoint(Rp_i, gamma_i, z_apex);

        // Next strand's entry point
        const double R_in_next = outerRadius * (1.0 + outerJitter * std::sin(2.0 * alpha_next));
        const Point3D Q_next   = polarPoint(R_in_next, alpha_next, 0.0);

        if (polygon.empty()) {
            appendSegment(polygon, Q_in,  P_in,   rampSubdivisions, false);
        } else {
            appendSegment(polygon, Q_in,  P_in,   rampSubdivisions, true);
        }

        appendSegment(polygon, P_in,  P_out,  chordSubdivisions, true);
        appendSegment(polygon, P_out, Q_out,  rampSubdivisions,  true);
        appendSegment(polygon, Q_out, A,      petalSubdivisions, true);
        appendSegment(polygon, A,     Q_next, petalSubdivisions, true);
    }

    return polygon;
}

// --------------------------------------------------
// Pretty-print permutation
// --------------------------------------------------
void printPermutation(const std::vector<int>& perm) {
    std::cout << "{ ";
    for (std::size_t i = 0; i < perm.size(); ++i) {
        std::cout << perm[i];
        if (i + 1 < perm.size()) {
            std::cout << ", ";
        }
    }
    std::cout << " }";
}

// --------------------------------------------------
// Main
// --------------------------------------------------
int main() {
    try {
        // Example: replace this with the petal permutation you want
        std::vector<int> petalPermutation{1,2,3,4,5};

        std::cout << "Petal permutation = ";
        printPermutation(petalPermutation);
        std::cout << '\n';

        if (!isValidPetalPermutation(petalPermutation)) {
            std::cerr << "Error: invalid petal permutation.\n";
            return 1;
        }

        std::vector<Point3D> polygon = buildPetalPolygon(
            petalPermutation,
            0.35, // innerRadius
            2.0,  // outerRadius
            3.5,  // petalRadius
            0.8,  // heightScale
            6,    // rampSubdivisions
            16,   // chordSubdivisions
            12    // petalSubdivisions
        );

        const std::string outputFile = "petal_knot.tsv";
        writePolygonTSV(polygon, outputFile);

        std::cout << "Generated polygon with " << polygon.size() << " vertices.\n";
        std::cout << "Wrote file: " << outputFile << '\n';
        std::cout << "Now run Knoodle on it.\n";

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << '\n';
        return 1;
    }
}