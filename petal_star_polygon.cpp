#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

struct Point3D {
    double x;
    double y;
    double z;
};

bool isValidPermutation(const std::vector<int>& perm) {
    const int n = static_cast<int>(perm.size());

    if (n == 0 || n % 2 == 0) {
        return false;
    }

    std::vector<int> sorted = perm;
    std::sort(sorted.begin(), sorted.end());

    for (int i = 0; i < n; ++i) {
        if (sorted[i] != i + 1) {
            return false;
        }
    }

    return true;
}

double heightFromValue(int value, int n, double heightScale) {
    const double center = (n + 1) / 2.0;
    return heightScale * (value - center);
}

Point3D polarPoint(double radius, double angle, double z) {
    return {
        radius * std::cos(angle),
        radius * std::sin(angle),
        z
    };
}

void writeTSV(const std::vector<Point3D>& points, const std::string& filename) {
    std::ofstream out(filename);

    if (!out) {
        throw std::runtime_error("Could not open output file: " + filename);
    }

    out << std::fixed << std::setprecision(10);

    for (const auto& p : points) {
        out << p.x << '\t' << p.y << '\t' << p.z << '\n';
    }
}

std::vector<int> parsePermutation(int argc, char* argv[]) {
    std::vector<int> perm;

    for (int i = 1; i < argc; ++i) {
        perm.push_back(std::stoi(argv[i]));
    }

    if (!isValidPermutation(perm)) {
        throw std::invalid_argument("Input must be an odd-length permutation of 1..N.");
    }

    return perm;
}

std::vector<Point3D> buildStarPolygon(
    const std::vector<int>& perm,
    double innerRadius,
    double outerRadius,
    double angleOffset,
    double heightScale
) {
    const int n = static_cast<int>(perm.size());
    const double twoPi = 2.0 * M_PI;

    std::vector<Point3D> polygon;
    polygon.reserve(3 * n + 1);

    for (int i = 0; i < n; ++i) {
        const double theta    = twoPi * i / n;
        const double zBefore  = heightFromValue(perm[(i - 1 + n) % n], n, heightScale);
        const double zAfter   = heightFromValue(perm[i], n, heightScale);

        Point3D centerBefore = polarPoint(innerRadius, theta - angleOffset, zBefore);
        Point3D outer        = polarPoint(outerRadius, theta, 0.0);
        Point3D centerAfter  = polarPoint(innerRadius, theta + angleOffset, zAfter);

        polygon.push_back(centerBefore);
        polygon.push_back(outer);
        polygon.push_back(centerAfter);
    }

    // Explicitly close the polygon for file formats/tools that like this.
    polygon.push_back(polygon.front());

    return polygon;
}

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            std::cerr << "Usage:\n";
            std::cerr << "  " << argv[0] << " <odd-length permutation>\n";
            std::cerr << "Example:\n";
            std::cerr << "  " << argv[0] << " 1 2 3 4 5\n";
            return 1;
        }

        std::vector<int> perm = parsePermutation(argc, argv);

        const double innerRadius = 0.08;
        const double outerRadius = 2.0;
        const double angleOffset = 0.08;
        const double heightScale = 0.25;

        std::vector<Point3D> polygon = buildStarPolygon(
            perm,
            innerRadius,
            outerRadius,
            angleOffset,
            heightScale
        );

        writeTSV(polygon, "petal_star.tsv");

        std::cout << "Wrote petal_star.tsv with "
                  << polygon.size() << " vertices.\n";

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}