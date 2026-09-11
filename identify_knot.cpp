// identify_knot.cpp -- read a PD code, compute the Alexander polynomial,
// look up the knot.
//
// Pipeline (kept separate from petal_star_clean.cpp on purpose):
//   ./petal_star_clean PERM           -> petal_star.tsv (3D polygon)
//   ./knoodletool petal_star.tsv      -> PD code on stdout
//   ./identify_knot < pd_code.txt     -> Alexander polynomial + knot name
//
// Input format: knoodletool's PD output.
//   Each crossing is a whitespace-separated line of 5 integers
//       e1 e2 e3 e4 sign
//   where e1..e4 are edge labels and sign is +1 (right-handed) or -1 (left).
//   Lines that don't parse as 5 integers (e.g. "k" or "s" separators)
//   are ignored, so we can pipe knoodletool's raw output through.
//
// What the program does, step by step, with the math spelled out:
//
//   1. Parse the PD code into a list of crossings (a, b, c, d, sign).
//      Convention: at a crossing X[a,b,c,d],
//          a = incoming under-edge, c = outgoing under-edge,
//          b and d are the two over-strand edges adjacent to the crossing.
//      The over-strand is unbroken at the crossing, so edges b and d
//      belong to the same ARC. (An "arc" is a maximal piece of the knot
//      bounded by under-crossings; for n crossings there are n arcs.)
//
//   2. Identify arcs by union-find: union(b, d) at each crossing.
//      The resulting connected components are the arcs.
//
//   3. Build the n x n Alexander matrix M(t) over Z[t]. One row per
//      crossing, one column per arc. For a positive crossing with
//      over-arc o, incoming-under u_in, outgoing-under u_out:
//          M[row, o]     += 1 - t
//          M[row, u_in]  += -1
//          M[row, u_out] += t
//      For a negative crossing:
//          M[row, o]     += 1 - t
//          M[row, u_in]  += t
//          M[row, u_out] += -1
//      Each row sums to 0, so M(t) is singular -- expected.
//
//   4. Delete the last row and last column to get the (n-1) x (n-1)
//      reduced Alexander matrix A(t). Its determinant is the Alexander
//      polynomial of the knot, up to multiplication by +/- t^k.
//
//   5. We compute det A(t) by evaluating at enough integer values of t
//      (numerically, using Bareiss-style integer Gaussian elimination)
//      and interpolating. Cheaper than polynomial Bareiss, and the
//      arithmetic stays in int64 for crossings <= 13.
//
//   6. Strip leading and trailing zero coefficients of det A(t). That
//      removes the t^k ambiguity, leaving a polynomial whose lowest
//      and highest coefficients are both nonzero.
//
//   7. Normalize the sign so that Delta(1) = +1 (standard convention;
//      for any knot, |Delta(1)| = 1).
//
//   8. Look up the coefficient list in a hard-coded table of prime
//      knots up to 7 crossings. Print Delta(t) and the knot name (or
//      "unknown" if the polynomial isn't in the table).

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using i64 = std::int64_t;

// ---------------------------------------------------------------------------
// Polynomial helpers. We represent a polynomial in t as a vector<i64>
// where coef[i] is the coefficient of t^i. Negative powers are handled
// by shifting (we strip trailing zeros at the end).
// ---------------------------------------------------------------------------

static std::string polyToString(const std::vector<i64>& p) {
    if (p.empty()) return "0";
    std::ostringstream os;
    bool first = true;
    for (int i = static_cast<int>(p.size()) - 1; i >= 0; --i) {
        i64 c = p[i];
        if (c == 0) continue;
        if (first) {
            if (c < 0) os << "-";
        } else {
            os << (c < 0 ? " - " : " + ");
        }
        i64 a = std::abs(c);
        if (a != 1 || i == 0) os << a;
        if (i > 0) {
            os << "t";
            if (i > 1) os << "^" << i;
        }
        first = false;
    }
    return os.str();
}

// Strip leading (highest-power) zero coefficients, then trailing
// (lowest-power) zero coefficients. The remaining polynomial has
// nonzero lowest and highest coefficients (or is empty if input was 0).
static std::vector<i64> stripZeros(std::vector<i64> p) {
    while (!p.empty() && p.back() == 0) p.pop_back();
    if (p.empty()) return p;
    size_t lead = 0;
    while (lead < p.size() && p[lead] == 0) ++lead;
    return std::vector<i64>(p.begin() + lead, p.end());
}

// Normalize: strip zeros, then make Delta(1) = +1 by negating if needed.
static std::vector<i64> normalize(std::vector<i64> p) {
    p = stripZeros(std::move(p));
    if (p.empty()) return p;
    i64 sum = 0;
    for (i64 c : p) sum += c;
    if (sum < 0) {
        for (i64& c : p) c = -c;
    }
    return p;
}

// ---------------------------------------------------------------------------
// Union-find for arc identification.
// ---------------------------------------------------------------------------

struct DSU {
    std::vector<int> p;
    DSU(int n) : p(n) {
        for (int i = 0; i < n; ++i) p[i] = i;
    }
    int find(int x) {
        while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
        return x;
    }
    void unite(int a, int b) {
        a = find(a); b = find(b);
        if (a != b) p[a] = b;
    }
};

// ---------------------------------------------------------------------------
// PD code parsing. Accepts whitespace-separated 5-int lines; any line
// that doesn't yield exactly 5 ints is silently skipped (so knoodletool
// separators like "k" / "s" pass through harmlessly).
// ---------------------------------------------------------------------------

struct Crossing {
    i64 a, b, c, d;
    int sign;  // +1 right-handed, -1 left-handed
};

static std::vector<Crossing> readPDCode(std::istream& in) {
    std::vector<Crossing> X;
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream ss(line);
        i64 a, b, c, d; int s;
        if (!(ss >> a >> b >> c >> d >> s)) continue;
        std::string extra;
        if (ss >> extra) continue;  // line has too many tokens; skip
        if (s != 1 && s != -1) {
            throw std::runtime_error("Crossing sign must be +1 or -1");
        }
        X.push_back({a, b, c, d, s});
    }
    return X;
}

// ---------------------------------------------------------------------------
// Arc identification. After this we have:
//   arc_of[edge] = arc index in 0..n-1
// for every edge label appearing in the PD code.
// ---------------------------------------------------------------------------

static std::unordered_map<i64, int> assignArcs(const std::vector<Crossing>& X) {
    // First, gather all distinct edge labels and assign each an internal id.
    std::unordered_map<i64, int> idx;
    auto getId = [&](i64 e) {
        auto it = idx.find(e);
        if (it != idx.end()) return it->second;
        int id = static_cast<int>(idx.size());
        idx[e] = id;
        return id;
    };
    for (const auto& c : X) {
        getId(c.a); getId(c.b); getId(c.c); getId(c.d);
    }
    // Union b and d at every crossing (over-strand is unbroken there).
    DSU dsu(static_cast<int>(idx.size()));
    for (const auto& c : X) {
        dsu.unite(idx[c.b], idx[c.d]);
    }
    // Re-number the connected components 0..n-1.
    std::unordered_map<int, int> rootToArc;
    std::unordered_map<i64, int> arc_of;
    for (auto& kv : idx) {
        int root = dsu.find(kv.second);
        auto it = rootToArc.find(root);
        int aid;
        if (it == rootToArc.end()) {
            aid = static_cast<int>(rootToArc.size());
            rootToArc[root] = aid;
        } else {
            aid = it->second;
        }
        arc_of[kv.first] = aid;
    }
    return arc_of;
}

// ---------------------------------------------------------------------------
// Build the (n-1) x (n-1) reduced Alexander matrix at a SPECIFIC integer
// value of t. Each entry is an int64. Last row and last column dropped.
// ---------------------------------------------------------------------------

static std::vector<std::vector<i64>>
buildReducedAlexanderAt(const std::vector<Crossing>& X,
                        const std::unordered_map<i64, int>& arc_of,
                        i64 t)
{
    const int n = static_cast<int>(X.size());
    std::vector<std::vector<i64>> M(n, std::vector<i64>(n, 0));
    for (int row = 0; row < n; ++row) {
        const Crossing& c = X[row];
        int o     = arc_of.at(c.b);   // over-arc (== arc_of[c.d])
        int u_in  = arc_of.at(c.a);
        int u_out = arc_of.at(c.c);
        i64 v_over   = 1 - t;
        i64 v_in     = (c.sign == 1) ? -1 :  t;
        i64 v_out    = (c.sign == 1) ?  t : -1;
        M[row][o]     += v_over;
        M[row][u_in]  += v_in;
        M[row][u_out] += v_out;
    }
    // Drop last row and last column.
    M.pop_back();
    for (auto& row : M) row.pop_back();
    return M;
}

// ---------------------------------------------------------------------------
// Integer determinant via Bareiss algorithm (fraction-free Gaussian
// elimination). Stays in int64 for our problem size.
// Returns 0 if matrix is singular at this t.
// ---------------------------------------------------------------------------

static i64 modPow(i64 b, i64 e, i64 p) {
    i64 r = 1; b %= p;
    while (e > 0) {
        if (e & 1) r = static_cast<i64>((__int128)r * b % p);
        b = static_cast<i64>((__int128)b * b % p);
        e >>= 1;
    }
    return r;
}

static i64 modInv(i64 a, i64 p) {
    a %= p; if (a < 0) a += p;
    return modPow(a, p - 2, p);           // p prime
}

// Determinant mod p by Gaussian elimination. The old fraction-free Bareiss
// version overflowed int64 once the diagram had more than ~13 crossings
// (the intermediate minors grow like t^n); working modulo a word-sized
// prime keeps every intermediate bounded by p^2, and the true integer
// coefficients are recovered afterwards by CRT.
static i64 detMod(std::vector<std::vector<i64>> M, i64 p) {
    const int n = static_cast<int>(M.size());
    if (n == 0) return 1 % p;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            M[i][j] %= p;
            if (M[i][j] < 0) M[i][j] += p;
        }
    i64 det = 1;
    for (int i = 0; i < n; ++i) {
        int piv = -1;
        for (int k = i; k < n; ++k) if (M[k][i] != 0) { piv = k; break; }
        if (piv == -1) return 0;
        if (piv != i) { std::swap(M[i], M[piv]); det = p - det; if (det == p) det = 0; }
        det = static_cast<i64>((__int128)det * M[i][i] % p);
        const i64 inv = modInv(M[i][i], p);
        for (int j = i + 1; j < n; ++j) {
            if (M[j][i] == 0) continue;
            const i64 f = static_cast<i64>((__int128)M[j][i] * inv % p);
            for (int k = i; k < n; ++k) {
                M[j][k] = static_cast<i64>((M[j][k] - (__int128)f * M[i][k]) % p);
                if (M[j][k] < 0) M[j][k] += p;
            }
        }
    }
    return det;
}

// ---------------------------------------------------------------------------
// Lagrange interpolation through (xs[i], ys[i]) pairs, evaluated to
// return the coefficient list of the unique polynomial of degree
// <= xs.size()-1. Uses rational arithmetic via __int128 to stay exact.
//
// We know a priori that det A(t) is an integer polynomial of degree
// at most (n-1) where n = crossing count, so we use (n+1) sample
// points around 0 and the answer is integer.
// ---------------------------------------------------------------------------

// Lagrange interpolation modulo p, returning coefficients (lowest power
// first) of the unique polynomial of degree <= xs.size()-1 through the
// given points.
static std::vector<i64> lagrangeInterpolateMod(const std::vector<i64>& xs,
                                               const std::vector<i64>& ys,
                                               i64 p) {
    const int m = static_cast<int>(xs.size());
    std::vector<i64> coefs(m, 0);
    for (int i = 0; i < m; ++i) {
        // num(x) = prod_{j != i} (x - x_j),  denom = prod_{j != i} (x_i - x_j)
        std::vector<i64> num(1, 1);
        i64 denom = 1;
        for (int j = 0; j < m; ++j) {
            if (j == i) continue;
            std::vector<i64> nx(num.size() + 1, 0);
            i64 xj = xs[j] % p; if (xj < 0) xj += p;
            for (std::size_t k = 0; k < num.size(); ++k) {
                nx[k + 1] = (nx[k + 1] + num[k]) % p;
                nx[k]     = (nx[k] - (__int128)num[k] * xj % p + p) % p;
            }
            num = std::move(nx);
            i64 diff = (xs[i] - xs[j]) % p; if (diff < 0) diff += p;
            denom = static_cast<i64>((__int128)denom * diff % p);
        }
        i64 yi = ys[i] % p; if (yi < 0) yi += p;
        const i64 scale = static_cast<i64>((__int128)yi * modInv(denom, p) % p);
        for (std::size_t k = 0; k < num.size(); ++k)
            coefs[k] = static_cast<i64>((coefs[k] + (__int128)num[k] * scale) % p);
    }
    return coefs;
}

// Primes just under 2^31, so products fit comfortably in __int128 and the
// combined modulus (~9.9e27) dwarfs any Alexander coefficient we meet.
static const i64 CRT_PRIMES[3] = { 2147483647LL, 2147483629LL, 2147483587LL };

// Combine per-prime residues into the true integer coefficients, taking
// symmetric representatives so negative coefficients come back correctly.
static std::vector<i64> crtLift(const std::vector<std::vector<i64>>& residues) {
    const std::size_t m = residues[0].size();
    __int128 M = 1;
    for (i64 pr : CRT_PRIMES) M *= pr;

    std::vector<i64> out(m);
    for (std::size_t k = 0; k < m; ++k) {
        __int128 acc = 0;
        for (int i = 0; i < 3; ++i) {
            const i64 p  = CRT_PRIMES[i];
            const __int128 Mi = M / p;
            const i64 Mi_mod = static_cast<i64>(Mi % p);
            const i64 inv = modInv(Mi_mod, p);
            const i64 term = static_cast<i64>((__int128)residues[i][k] * inv % p);
            acc = (acc + Mi * term) % M;
        }
        if (acc > M / 2) acc -= M;          // symmetric representative
        if (acc > static_cast<__int128>(std::numeric_limits<i64>::max()) ||
            acc < static_cast<__int128>(std::numeric_limits<i64>::min())) {
            throw std::runtime_error("Coefficient too large for int64 -- add another CRT prime.");
        }
        out[k] = static_cast<i64>(acc);
    }
    return out;
}

// ---------------------------------------------------------------------------
// Knot lookup table. Keyed on the normalized Alexander polynomial
// coefficient list (lowest power of t first), value is the knot name.
// Extended later for higher crossings.
// ---------------------------------------------------------------------------

static const std::map<std::vector<i64>, std::string>& knotTable() {
    static const std::map<std::vector<i64>, std::string> t = {
        { {1},                            "0_1 (unknot)" },
        { {1,-1,1},                       "3_1 (trefoil)" },
        { {-1,3,-1},                      "4_1 (figure-eight)" },
        { {1,-1,1,-1,1},                  "5_1" },
        { {2,-3,2},                       "5_2" },
        { {-2,5,-2},                      "6_1" },
        { {-1,3,-3,3,-1},                 "6_2" },
        { {1,-3,5,-3,1},                  "6_3" },
        { {1,-1,1,-1,1,-1,1},             "7_1" },
        { {3,-5,3},                       "7_2" },
        { {2,-3,3,-3,2},                  "7_3" },
        { {4,-7,4},                       "7_4" },
        { {2,-4,5,-4,2},                  "7_5" },
        { {-1,5,-7,5,-1},                 "7_6" },
        { {1,-5,9,-5,1},                  "7_7" },
    };
    return t;
}

// ---------------------------------------------------------------------------
// Main.
// ---------------------------------------------------------------------------

static void printUsage(const std::string& prog) {
    std::cout
        << "Usage:\n"
        << "  " << prog << "                  read PD code from stdin\n"
        << "  " << prog << " --input FILE     read PD code from FILE\n\n"
        << "Input is knoodletool's PD code format: each crossing on its own\n"
        << "line as 5 whitespace-separated integers (e1 e2 e3 e4 sign).\n"
        << "Non-matching lines (\"k\", \"s\", blanks, etc.) are skipped.\n";
}

int main(int argc, char* argv[]) {
    try {
        std::string inputFile;
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--help" || arg == "-h") { printUsage(argv[0]); return 0; }
            if (arg == "--input") {
                if (i + 1 >= argc) throw std::invalid_argument("--input requires a filename");
                inputFile = argv[++i];
            } else {
                throw std::invalid_argument("Unknown argument: " + arg);
            }
        }

        std::vector<Crossing> X;
        if (inputFile.empty()) {
            X = readPDCode(std::cin);
        } else {
            std::ifstream in(inputFile);
            if (!in) throw std::runtime_error("Cannot open input file: " + inputFile);
            X = readPDCode(in);
        }
        const int n = static_cast<int>(X.size());
        if (n == 0) {
            // No crossings parsed -- treat as unknot.
            std::cout << "Crossings : 0\n";
            std::cout << "Delta(t)  : 1\n";
            std::cout << "Knot      : 0_1 (unknot)\n";
            std::cout << "Coefs     : [1]\n";
            return 0;
        }

        auto arc_of = assignArcs(X);
        // Sanity check: at every crossing, b and d should map to the same arc.
        for (const auto& c : X) {
            if (arc_of.at(c.b) != arc_of.at(c.d)) {
                throw std::runtime_error("Inconsistent PD code: over-strand edges b and d landed in different arcs");
            }
        }

        // The reduced Alexander matrix is (n-1) x (n-1). Its determinant
        // is a polynomial of degree at most (n-1) in t. Sample at
        // 2*(n-1)+1 = 2n-1 integer points to be safe (the matrix can
        // actually have degree exactly n-1 only when every row contributes
        // a t-factor, but we don't need to be tight here).
        const int numPoints = std::max(2, 2 * n - 1);
        std::vector<i64> xs;
        for (int k = 0; k < numPoints; ++k)
            xs.push_back(static_cast<i64>(k - numPoints / 2));  // centered on 0

        // Evaluate det A(t) modulo each CRT prime, interpolate mod that
        // prime, then lift. This replaces the old int64 Bareiss + rational
        // Lagrange path, which overflowed past ~13 crossings.
        std::vector<std::vector<i64>> residues;
        for (i64 pr : CRT_PRIMES) {
            std::vector<i64> ys;
            ys.reserve(xs.size());
            for (i64 t : xs) {
                auto M = buildReducedAlexanderAt(X, arc_of, t);
                ys.push_back(detMod(std::move(M), pr));
            }
            residues.push_back(lagrangeInterpolateMod(xs, ys, pr));
        }

        auto poly = crtLift(residues);
        poly = normalize(std::move(poly));

        std::cout << "Crossings : " << n << "\n";
        std::cout << "Delta(t)  : " << polyToString(poly) << "\n";

        const auto& tbl = knotTable();
        auto it = tbl.find(poly);
        if (it != tbl.end()) {
            std::cout << "Knot      : " << it->second << "\n";
        } else {
            std::cout << "Knot      : unknown (polynomial not in table)\n";
        }

        // Always emit the coefficient list (lowest power first) so that
        // downstream tools -- notably torus_match -- can consume it without
        // re-parsing the pretty-printed polynomial.
        std::cout << "Coefs     : [";
        for (size_t i = 0; i < poly.size(); ++i) {
            if (i) std::cout << ", ";
            std::cout << poly[i];
        }
        std::cout << "]\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
