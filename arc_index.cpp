// arc_index.cpp -- grid diagrams and arc index for petal knots.
//
// Background, in one paragraph. An ARC PRESENTATION of a knot puts it into
// finitely many half-planes ("pages") around a binding axis, one arc per
// page, with the arc endpoints on the axis. The minimum number of pages is
// the ARC INDEX alpha(K). Flattening the pages turns an arc presentation
// into a GRID DIAGRAM: an n x n array with one X and one O per row and per
// column, joined by horizontal and vertical segments, verticals always over
// horizontals. So alpha(K) is also the minimal grid number, and the whole
// invariant is a statement about pairs of permutations -- no geometry.
//
// The point of this file for the thesis: a petal permutation ALREADY IS an
// arc presentation. A petal projection with n petals passes through the
// centre n times at heights pi_1..pi_n, and petal t joins the passage at
// height pi_t to the one at height pi_{t+1}. Read the central axis as the
// binding axis and each petal is an arc in its own page. Therefore
//
//      alpha(K) <= n                              (*)
//
// for every knot K represented by a length-n petal permutation, and in
// particular alpha(K) <= p(K) where p is the petal number. Bound (*) is the
// starting point; simplifying the grid by Cromwell moves drives it down
// towards the true arc index.
//
// Two theorems make the resulting table checkable:
//   * Bae-Park: for a prime ALTERNATING knot, alpha(K) = c(K) + 2 exactly.
//   * Jin-Park: a prime link is non-alternating IFF alpha(K) <= c(K).
// The alternating case is a hard prediction, so agreement between the search
// output and c+2 on every alternating knot is a real test of the code, and
// lends credibility to the non-alternating rows where no closed form exists.
//
// Modes:
//   --perm "1 2 3 4 5"    one permutation: grid, crossings, arc-index search
//   --sweep N             every rotation-orbit rep of length N (pi_1 = 1):
//                         identify the knot and search for its arc index
//   --tabulate FILE       aggregate a --sweep TSV into arc index by crossing
//   --verify N            check the petal->grid map: the grid polygon and the
//                         petal_star_clean polygon must name the same knot
//
// Identification shells out to knoodleidentify, exactly as the other
// sampling tools do, so the naming conventions stay consistent.

#include "grid_diagram.hpp"
#include "knot_naming.hpp"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<int> parsePerm(const std::string& s) {
    std::istringstream in(s);
    std::vector<int> v;
    int x;
    while (in >> x) v.push_back(x);
    return v;
}

std::string joinPerm(const std::vector<int>& v) {
    std::ostringstream out;
    for (std::size_t i = 0; i < v.size(); ++i) out << (i ? " " : "") << v[i];
    return out.str();
}

void writeGridTSV(const GridDiagram& g, const std::string& path) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open " + path);
    g.writePolygonTSV(out);
}

// Crossing number of the identified knot, summed over prime factors. Returns
// -1 when the knot could not be named.
int crossingNumberOf(const std::vector<knotname::Summand>& s, bool named) {
    if (!named) return -1;
    int c = 0;
    for (const auto& f : s) c += f.c * f.mult;
    return c;
}

bool allAlternating(const std::vector<knotname::Summand>& s) {
    for (const auto& f : s) if (!f.alternating) return false;
    return true;
}

// -----------------------------------------------------------------------
// --perm
// -----------------------------------------------------------------------
int runPerm(const std::string& permStr, std::uint64_t seed, int steps, int slack,
            const std::string& tsvOut, bool identify) {
    const std::vector<int> perm = parsePerm(permStr);
    GridDiagram g = GridDiagram::fromPetalPermutation(perm);
    const int n = g.size();

    std::cout << "permutation      " << joinPerm(perm) << "\n"
              << "grid size        " << n
              << "   (so arc index <= " << n << ")\n"
              << "components       " << g.componentCount() << "\n"
              << "grid crossings   " << g.crossingCount() << "\n"
              << g.toPermutationPair() << "\n\n"
              << g.toAscii() << "\n";

    if (!tsvOut.empty()) {
        writeGridTSV(g, tsvOut);
        std::cout << "grid polygon written to " << tsvOut << "\n";
    }

    std::string knot = "-";
    std::vector<knotname::Summand> summands;
    if (identify) {
        const std::string tmp = "/tmp/arc_index_grid.tsv";
        writeGridTSV(g, tmp);
        knot = knotname::identifyTSV(tmp, summands);
    }

    GridDiagram simplified = g;
    const int alpha = simplified.simplify(seed, steps, slack);

    std::cout << "\nafter Cromwell simplification\n"
              << "grid size        " << alpha << "\n"
              << "grid crossings   " << simplified.crossingCount() << "\n"
              << simplified.toPermutationPair() << "\n\n"
              << simplified.toAscii();

    if (identify) {
        const bool named = (knot != "parse-error" && knot.rfind("notfound", 0) != 0 &&
                            knot.rfind("unidentified", 0) != 0);
        const int c = crossingNumberOf(summands, named);
        std::cout << "\nknot             " << knot << "\n";
        if (named && c >= 0) {
            std::cout << "crossing number  " << c << "\n";
            if (!summands.empty() && summands.size() == 1 && summands[0].mult == 1 &&
                summands[0].alternating)
                std::cout << "Bae-Park alpha   " << c + 2
                          << (c + 2 == alpha ? "   (matches search)"
                                             : "   *** DISAGREES with search ***")
                          << "\n";
        }
    }
    std::cout << "\narc index <= " << alpha << "\n";
    return 0;
}

// -----------------------------------------------------------------------
// --sweep: one rep per rotation orbit, i.e. every permutation with pi_1 = 1.
// -----------------------------------------------------------------------
int runSweep(int n, std::uint64_t seed, int steps, int slack, std::ostream& out,
             const std::string& petalBin) {
    if (n < 3 || n % 2 == 0) throw std::invalid_argument("n must be odd and >= 3");

    const std::string petalTsv = "/tmp/arc_index_petal.tsv";

    out << "# arc index sweep, n=" << n << " (one rep per rotation orbit)\n";
    out << "permutation\tknot\tcrossings\talternating\tgrid_n\talpha_ub\tbae_park\n";

    std::vector<int> tail(n - 1);
    std::iota(tail.begin(), tail.end(), 2);
    int rows = 0, agree = 0, tested = 0;

    do {
        std::vector<int> perm;
        perm.reserve(n);
        perm.push_back(1);
        perm.insert(perm.end(), tail.begin(), tail.end());

        GridDiagram g = GridDiagram::fromPetalPermutation(perm);

        // Identify via the petal_star polygon -- the pipeline the rest of the
        // thesis uses -- so the arc index column is attached to a knot name
        // produced exactly the way every other table produces it.
        const std::string cmd = petalBin + " " + joinPerm(perm) +
                                " --output " + petalTsv + " > /dev/null 2>&1";
        if (std::system(cmd.c_str()) != 0)
            throw std::runtime_error("petal_star_clean failed on " + joinPerm(perm));

        std::vector<knotname::Summand> summands;
        const std::string knot = knotname::identifyTSV(petalTsv, summands);
        const bool named = (knot != "parse-error" && knot.rfind("notfound", 0) != 0 &&
                            knot.rfind("unidentified", 0) != 0);
        const int c = crossingNumberOf(summands, named);
        const bool alt = named && allAlternating(summands);

        GridDiagram simplified = g;
        const int alpha = simplified.simplify(seed + rows, steps, slack);

        std::string bp = "-";
        const bool primeAlt = named && summands.size() == 1 &&
                              summands[0].mult == 1 && summands[0].alternating;
        if (named && summands.empty()) bp = "2";            // unknot
        else if (primeAlt)             bp = std::to_string(c + 2);
        if (bp != "-") {
            ++tested;
            if (std::stoi(bp) == alpha) ++agree;
        }

        out << joinPerm(perm) << '\t' << knot << '\t'
            << (c >= 0 ? std::to_string(c) : "-") << '\t'
            << (named ? (alt ? "yes" : "no") : "-") << '\t'
            << g.size() << '\t' << alpha << '\t' << bp << '\n';
        ++rows;
    } while (std::next_permutation(tail.begin(), tail.end()));

    out << "# rows=" << rows << "  Bae-Park checkable=" << tested
        << "  agreeing=" << agree << '\n';
    return 0;
}

// -----------------------------------------------------------------------
// --tabulate: arc index by crossing number, from a sweep TSV.
// -----------------------------------------------------------------------
int runTabulate(const std::string& path, std::ostream& out) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);

    // crossing number -> arc index -> set of distinct knot names
    std::map<int, std::map<int, std::set<std::string>>> table;
    std::map<std::string, int> alphaOfKnot;
    std::map<std::string, int> cOfKnot;
    std::map<std::string, bool> altOfKnot;

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ls(line);
        std::string perm, knot, cs, alt, gridN, alphaS, bp;
        // permutation field contains spaces, so split on tabs.
        auto field = [&](std::string& dst) { return static_cast<bool>(std::getline(ls, dst, '\t')); };
        if (!(field(perm) && field(knot) && field(cs) && field(alt) &&
              field(gridN) && field(alphaS) && field(bp))) continue;
        if (knot == "knot" || cs == "-") continue;          // header / unnamed

        const int c = std::stoi(cs);
        const int a = std::stoi(alphaS);
        // Keep the smallest arc index found for each knot: different
        // permutations give different presentations of the same knot, and the
        // invariant is the minimum over all of them.
        auto it = alphaOfKnot.find(knot);
        if (it == alphaOfKnot.end() || a < it->second) alphaOfKnot[knot] = a;
        cOfKnot[knot] = c;
        altOfKnot[knot] = (alt == "yes");
    }

    for (const auto& kv : alphaOfKnot)
        table[cOfKnot[kv.first]][kv.second].insert(kv.first);

    out << "crossings\tarc_index\tcount\tknots\n";
    for (const auto& [c, byAlpha] : table)
        for (const auto& [alpha, knots] : byAlpha) {
            out << c << '\t' << alpha << '\t' << knots.size() << '\t';
            bool first = true;
            for (const auto& k : knots) { out << (first ? "" : ",") << k; first = false; }
            out << '\n';
        }

    out << "\n# Prediction check: prime alternating knots satisfy alpha = c + 2\n"
           "# (Bae-Park); connected sums satisfy alpha(K1 # K2) = alpha(K1) + alpha(K2) - 2;\n"
           "# prime non-alternating knots satisfy alpha <= c (Jin-Park).\n";
    out << "knot\tcrossings\talternating\talpha_found\tpredicted\tstatus\n";
    for (const auto& [knot, alpha] : alphaOfKnot) {
        const int c = cOfKnot[knot];
        const bool alt = altOfKnot[knot];
        const bool composite = knot.find('#') != std::string::npos;
        std::string status = "-";
        std::string predicted = (c == 0 ? "2" : std::to_string(c + 2));
        if (c == 0) {
            status = (alpha == 2 ? "ok" : "CHECK");
        } else if (composite) {
            // Arc index is additive on connected sums up to the shared pair
            // of binding points: alpha(K1 # K2) = alpha(K1) + alpha(K2) - 2.
            // Check it whenever every factor was also seen on its own.
            int sum = 0, factors = 0;
            bool haveAll = true;
            std::size_t pos = 0;
            while (pos <= knot.size()) {
                const std::size_t hit = knot.find('#', pos);
                const std::string factor =
                    knot.substr(pos, hit == std::string::npos ? std::string::npos : hit - pos);
                auto f = alphaOfKnot.find(factor);
                if (f == alphaOfKnot.end()) { haveAll = false; break; }
                sum += f->second; ++factors;
                if (hit == std::string::npos) break;
                pos = hit + 1;
            }
            if (!haveAll) { status = "composite"; predicted = "-"; }
            else {
                const int expect = sum - 2 * (factors - 1);
                predicted = std::to_string(expect);
                status = (alpha == expect ? "ok (additivity)" : "CHECK");
            }
        } else if (alt) {
            status = (alpha == c + 2 ? "ok" : "CHECK");
        } else {
            status = (alpha <= c ? "ok (Jin-Park)" : "CHECK");
        }
        out << knot << '\t' << c << '\t' << (alt ? "yes" : "no") << '\t'
            << alpha << '\t' << predicted << '\t' << status << '\n';
    }
    return 0;
}

// -----------------------------------------------------------------------
// --verify: the petal->grid map must preserve the knot type.
// -----------------------------------------------------------------------
int runVerify(int n, const std::string& petalBin, std::ostream& out) {
    if (n < 3 || n % 2 == 0) throw std::invalid_argument("n must be odd and >= 3");
    const std::string gridTsv  = "/tmp/arc_index_verify_grid.tsv";
    const std::string petalTsv = "/tmp/arc_index_verify_petal.tsv";

    std::vector<int> tail(n - 1);
    std::iota(tail.begin(), tail.end(), 2);
    int rows = 0, mismatches = 0;

    out << "permutation\tvia_petal_star\tvia_grid\tstatus\n";
    do {
        std::vector<int> perm{1};
        perm.insert(perm.end(), tail.begin(), tail.end());

        const std::string cmd = petalBin + " " + joinPerm(perm) +
                                " --output " + petalTsv + " > /dev/null 2>&1";
        if (std::system(cmd.c_str()) != 0)
            throw std::runtime_error("petal_star_clean failed");

        std::vector<knotname::Summand> s1, s2;
        const std::string viaPetal = knotname::identifyTSV(petalTsv, s1);

        writeGridTSV(GridDiagram::fromPetalPermutation(perm), gridTsv);
        const std::string viaGrid = knotname::identifyTSV(gridTsv, s2);

        const bool ok = (viaPetal == viaGrid);
        if (!ok) ++mismatches;
        out << joinPerm(perm) << '\t' << viaPetal << '\t' << viaGrid << '\t'
            << (ok ? "ok" : "MISMATCH") << '\n';
        ++rows;
    } while (std::next_permutation(tail.begin(), tail.end()));

    out << "# rows=" << rows << "  mismatches=" << mismatches << '\n';
    return mismatches == 0 ? 0 : 1;
}

void printUsage(const char* prog) {
    std::cout <<
        "Usage:\n"
        "  " << prog << " --perm \"1 2 3 4 5\" [--tsv FILE] [--no-identify]\n"
        "  " << prog << " --sweep N [--output FILE]\n"
        "  " << prog << " --tabulate SWEEP.tsv [--output FILE]\n"
        "  " << prog << " --verify N [--output FILE]\n"
        "\n"
        "Common options:\n"
        "  --seed S     RNG seed for the Cromwell-move search (default 12345)\n"
        "  --steps K    search steps per grid (default 40000)\n"
        "  --slack S    how far above the best size the search may wander (default 2)\n"
        "  --petal-star-clean PATH   default ./petal_star_clean\n"
        "\n"
        "--perm       build the grid diagram of one petal permutation, print it,\n"
        "             and search for its arc index.\n"
        "--sweep      do that for every rotation-orbit rep of length N, naming each\n"
        "             knot with knoodleidentify. Emits TSV.\n"
        "--tabulate   turn a sweep TSV into 'arc index by crossing number', and\n"
        "             check every prime alternating knot against Bae-Park (c + 2).\n"
        "--verify     confirm the petal->grid map preserves knot type by naming the\n"
        "             same permutation through both polygon constructions.\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        std::string mode, permStr, tsvOut, outputFile, tabFile;
        std::string petalBin = "./petal_star_clean";
        int n = -1, steps = 40000, slack = 2;
        std::uint64_t seed = 12345u;
        bool identify = true;

        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            auto need = [&](const char* what) -> std::string {
                if (i + 1 >= argc) throw std::invalid_argument(std::string(what) + " requires a value");
                return argv[++i];
            };
            if (a == "--help" || a == "-h") { printUsage(argv[0]); return 0; }
            else if (a == "--perm")     { mode = "perm";     permStr = need("--perm"); }
            else if (a == "--sweep")    { mode = "sweep";    n = std::stoi(need("--sweep")); }
            else if (a == "--verify")   { mode = "verify";   n = std::stoi(need("--verify")); }
            else if (a == "--tabulate") { mode = "tabulate"; tabFile = need("--tabulate"); }
            else if (a == "--tsv")      { tsvOut = need("--tsv"); }
            else if (a == "--output")   { outputFile = need("--output"); }
            else if (a == "--seed")     { seed = std::stoull(need("--seed")); }
            else if (a == "--steps")    { steps = std::stoi(need("--steps")); }
            else if (a == "--slack")    { slack = std::stoi(need("--slack")); }
            else if (a == "--petal-star-clean") { petalBin = need("--petal-star-clean"); }
            else if (a == "--no-identify") { identify = false; }
            else throw std::invalid_argument("Unknown argument: " + a);
        }

        if (mode.empty()) { printUsage(argv[0]); return 1; }

        std::ofstream ofs;
        std::ostream* out = &std::cout;
        if (!outputFile.empty()) {
            ofs.open(outputFile);
            if (!ofs) throw std::runtime_error("cannot open output: " + outputFile);
            out = &ofs;
        }

        if (mode == "perm")     return runPerm(permStr, seed, steps, slack, tsvOut, identify);
        if (mode == "sweep")    return runSweep(n, seed, steps, slack, *out, petalBin);
        if (mode == "tabulate") return runTabulate(tabFile, *out);
        if (mode == "verify")   return runVerify(n, petalBin, *out);
        printUsage(argv[0]);
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
