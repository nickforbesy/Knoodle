// torus_match.cpp -- decide whether an Alexander polynomial is that of a
// torus knot T(p,q), and if so report which one.
//
// Motivation: sample_const_diff.cpp recognises torus knots via a hand-typed
// TORUS_TABLE keyed on knot names, which cannot grow past the few entries
// someone has looked up by hand. The Alexander polynomial of a torus knot
// has a closed form,
//
//     Delta_{T(p,q)}(t) = (t^{pq} - 1)(t - 1) / ((t^p - 1)(t^q - 1)),
//
// so the recognition can be computed instead of looked up, for arbitrary
// p and q. That removes the table entirely and lets the constant-difference
// sweep run at any n.
//
// The search is cheap because the degree pins the candidates down:
//     deg Delta_{T(p,q)} = (p-1)(q-1) = 2 * genus,
// so p-1 must divide deg, leaving only a handful of (p,q) to test.
//
// Usage:
//   ./identify_knot < pd_code.txt | ./torus_match
//   ./torus_match --coefs "1,-1,0,1,0,-1,1"
//
// Reads identify_knot's "Coefs : [...]" line from stdin (other lines are
// ignored), or takes the list directly with --coefs.
//
// Exit status is 0 on a match, 2 when the polynomial is not a torus knot's.
// Alexander polynomials do not separate knots in general, so a match is
// strong evidence, not a proof.

#include <cstdint>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using i64 = std::int64_t;
using Poly = std::vector<i64>;   // coefficients, lowest power first

// ---------------------------------------------------------------------------
// Minimal exact polynomial arithmetic over Z.
// ---------------------------------------------------------------------------

static Poly polyTrim(Poly a) {
    while (!a.empty() && a.back() == 0) a.pop_back();
    return a;
}

static Poly polyMul(const Poly& a, const Poly& b) {
    if (a.empty() || b.empty()) return {};
    Poly r(a.size() + b.size() - 1, 0);
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j)
            r[i + j] += a[i] * b[j];
    return polyTrim(std::move(r));
}

// Exact division a / b. Returns false if the division leaves a remainder.
static bool polyDivExact(Poly a, const Poly& b, Poly& q) {
    a = polyTrim(std::move(a));
    if (b.empty()) return false;
    const std::size_t db = b.size() - 1;
    if (a.size() < b.size()) { q.clear(); return a.empty(); }
    q.assign(a.size() - db, 0);
    for (std::size_t k = q.size(); k-- > 0; ) {
        const i64 lead = a[k + db];
        if (lead % b[db] != 0) return false;
        const i64 c = lead / b[db];
        q[k] = c;
        if (c == 0) continue;
        for (std::size_t j = 0; j <= db; ++j)
            a[k + j] -= c * b[j];
    }
    return polyTrim(std::move(a)).empty();
}

// t^k - 1
static Poly tPowMinusOne(int k) {
    Poly p(static_cast<std::size_t>(k) + 1, 0);
    p[0] = -1;
    p[static_cast<std::size_t>(k)] = 1;
    return p;
}

// Normalise the way identify_knot does: strip zeros, then fix the sign so
// the lowest surviving coefficient is positive. Alexander polynomials are
// only defined up to +/- t^k, so this makes comparison well defined.
static Poly normalize(Poly p) {
    p = polyTrim(std::move(p));
    std::size_t lo = 0;
    while (lo < p.size() && p[lo] == 0) ++lo;
    if (lo == p.size()) return {};
    Poly r(p.begin() + static_cast<long>(lo), p.end());
    if (r.front() < 0)
        for (auto& c : r) c = -c;
    return r;
}

// Delta_{T(p,q)}(t), normalised.
static Poly torusAlexander(int p, int q) {
    Poly num = polyMul(tPowMinusOne(p * q), tPowMinusOne(1));
    Poly den = polyMul(tPowMinusOne(p), tPowMinusOne(q));
    Poly quo;
    if (!polyDivExact(std::move(num), den, quo)) return {};
    return normalize(std::move(quo));
}

static int gcdInt(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }

// Crossing number of T(p,q) for coprime 1 < p < q.
static int torusCrossingNumber(int p, int q) {
    const int a = p * (q - 1);
    const int b = q * (p - 1);
    return a < b ? a : b;
}

// ---------------------------------------------------------------------------
// Input parsing.
// ---------------------------------------------------------------------------

static Poly parseCoefList(const std::string& s) {
    Poly out;
    std::string cur;
    for (char ch : s) {
        if (ch == '-' || (ch >= '0' && ch <= '9')) cur += ch;
        else if (!cur.empty()) { out.push_back(std::stoll(cur)); cur.clear(); }
    }
    if (!cur.empty()) out.push_back(std::stoll(cur));
    return out;
}

// Pull the "Coefs     : [ ... ]" line out of identify_knot's report.
static bool readCoefsFromStream(std::istream& in, Poly& out) {
    std::string line;
    while (std::getline(in, line)) {
        const std::size_t c = line.find("Coefs");
        if (c == std::string::npos) continue;
        const std::size_t lb = line.find('[', c);
        if (lb == std::string::npos) continue;
        const std::size_t rb = line.find(']', lb);
        if (rb == std::string::npos) continue;
        out = parseCoefList(line.substr(lb + 1, rb - lb - 1));
        return true;
    }
    return false;
}

static std::string polyToString(const Poly& p) {
    if (p.empty()) return "0";
    std::ostringstream o;
    bool first = true;
    for (std::size_t i = p.size(); i-- > 0; ) {
        if (p[i] == 0) continue;
        if (!first) o << (p[i] < 0 ? " - " : " + ");
        else if (p[i] < 0) o << "-";
        const i64 a = p[i] < 0 ? -p[i] : p[i];
        if (a != 1 || i == 0) o << a;
        if (i >= 1) o << "t";
        if (i >= 2) o << "^" << i;
        first = false;
    }
    return o.str();
}

static void printUsage(const char* prog) {
    std::cout
        << "Usage:\n"
        << "  ./identify_knot < pd.txt | " << prog << "\n"
        << "  " << prog << " --coefs \"1,-1,0,1,0,-1,1\"\n\n"
        << "Matches an Alexander polynomial against the closed form\n"
        << "  Delta_T(p,q) = (t^pq - 1)(t - 1) / ((t^p - 1)(t^q - 1)).\n\n"
        << "Options:\n"
        << "  --coefs LIST   coefficients, lowest power first\n"
        << "  --max-p N      largest p to try (default 64)\n"
        << "  --quiet        print only the torus label\n";
}

int main(int argc, char* argv[]) {
    Poly delta;
    bool haveInput = false;
    bool quiet = false;
    int maxP = 64;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--help" || a == "-h") { printUsage(argv[0]); return 0; }
        else if (a == "--quiet") { quiet = true; }
        else if (a == "--coefs") {
            if (i + 1 >= argc) { std::cerr << "--coefs needs a list\n"; return 1; }
            delta = normalize(parseCoefList(argv[++i]));
            haveInput = true;
        } else if (a == "--max-p") {
            if (i + 1 >= argc) { std::cerr << "--max-p needs an integer\n"; return 1; }
            maxP = std::stoi(argv[++i]);
        } else {
            std::cerr << "Unknown argument: " << a << "\n";
            return 1;
        }
    }

    if (!haveInput) {
        Poly raw;
        if (!readCoefsFromStream(std::cin, raw)) {
            std::cerr << "Error: no \"Coefs : [...]\" line found on stdin.\n"
                      << "Pipe identify_knot's output in, or use --coefs.\n";
            return 1;
        }
        delta = normalize(std::move(raw));
    }

    if (delta.empty()) { std::cerr << "Error: empty polynomial.\n"; return 1; }

    const int deg = static_cast<int>(delta.size()) - 1;

    // The unknot: Delta = 1.
    if (deg == 0) {
        if (quiet) std::cout << "unknot\n";
        else {
            std::cout << "Delta(t)  : " << polyToString(delta) << "\n";
            std::cout << "Genus     : 0\n";
            std::cout << "Torus     : yes -- unknot\n";
        }
        return 0;
    }

    if (deg % 2 != 0) {
        if (quiet) std::cout << "no\n";
        else {
            std::cout << "Delta(t)  : " << polyToString(delta) << "\n";
            std::cout << "Torus     : no -- odd degree " << deg << "\n";
        }
        return 2;
    }

    // (p-1)(q-1) = deg, so p-1 runs over the divisors of deg.
    std::vector<std::string> hits;
    int bestP = 0, bestQ = 0;
    for (int p = 2; p <= maxP; ++p) {
        const int pm = p - 1;
        if (deg % pm != 0) continue;
        const int q = deg / pm + 1;
        if (q <= p) continue;                  // canonical order p < q
        if (gcdInt(p, q) != 1) continue;       // otherwise a link, not a knot
        if (torusAlexander(p, q) == delta) {
            hits.push_back("T(" + std::to_string(p) + "," + std::to_string(q) + ")");
            if (!bestP) { bestP = p; bestQ = q; }
        }
    }

    if (quiet) {
        if (hits.empty()) std::cout << "no\n";
        else {
            for (std::size_t i = 0; i < hits.size(); ++i) {
                if (i) std::cout << "|";
                std::cout << hits[i];
            }
            std::cout << "\n";
        }
        return hits.empty() ? 2 : 0;
    }

    std::cout << "Delta(t)  : " << polyToString(delta) << "\n";
    std::cout << "Genus     : " << deg / 2 << "  (if fibred)\n";

    if (hits.empty()) {
        std::cout << "Torus     : no -- no T(p,q) has this Alexander polynomial\n";
        return 2;
    }

    std::cout << "Torus     : yes -- ";
    for (std::size_t i = 0; i < hits.size(); ++i) {
        if (i) std::cout << " or ";
        std::cout << hits[i];
    }
    std::cout << "\n";
    if (hits.size() > 1)
        std::cout << "Note      : several candidates share this polynomial\n";
    std::cout << "Crossings : " << torusCrossingNumber(bestP, bestQ)
              << " expected for " << hits[0] << "\n";
    return 0;
}
