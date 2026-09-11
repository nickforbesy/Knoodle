// sample_const_diff.cpp -- test the conjecture that every constant-difference
// petal permutation produces a torus knot.
//
// A constant-difference permutation of length n (odd) with step d, where
// gcd(d, n) = 1, is
//     pi_i = ((i - 1) * d) mod n + 1,   i = 1, ..., n.
// These form an arithmetic progression mod n. There are phi(n) of them
// (up to cyclic shift of the start, which doesn't change the knot).
//
// For each d coprime to n this program:
//   1. builds the permutation,
//   2. shells out to petal_star_clean and knoodleidentify,
//   3. converts KnotSymbol to compact form (Parsley chirality convention),
//   4. looks the base knot up in a small torus-knot table and annotates,
//   5. prints one TSV row: n, d, permutation, knot, torus_type.
//
// Torus-knot recognition uses a hardcoded table keyed by compact base
// name "<crossing>.<index>" (chirality prefix stripped). Extend
// TORUS_TABLE below to recognize more. Anything not in the table is
// annotated "?" so the conjecture-tester can spot unknowns at a glance.
//
// Requires petal_star_clean at ./petal_star_clean and knoodleidentify on
// PATH with KNOODLE_KLUT_DIR set in the environment.

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// ---------- knot-name conversion (Parsley chirality convention) ----------

static std::string chiralityPrefix(const std::string& sym) {
    bool hasE = false, hasM = false;
    std::size_t start = 0;
    while (start < sym.size()) {
        std::size_t end = sym.find('/', start);
        if (end == std::string::npos) end = sym.size();
        std::string tok = sym.substr(start, end - start);
        if (tok == "e")      hasE = true;
        else if (tok == "m") hasM = true;
        start = end + 1;
    }
    if (hasE && hasM) return "a";  // amphichiral
    // Parsley uses the opposite mirror convention from knoodleidentify, so
    // the two branches below are swapped to make our output match Parsley's.
    // Original (knoodleidentify convention):
    //     if (hasE) return "p";
    //     if (hasM) return "m";
    if (hasE)         return "m";
    if (hasM)         return "p";
    return "h";
}

// Parsed summand of the knoodleidentify result: crossing count, table
// index, alternating flag, chirality prefix, and multiplicity.
struct KnotSummand {
    int c = 0;
    int idx = 0;
    bool alternating = true;
    std::string prefix;
    int mult = 1;
};

static std::vector<KnotSummand> parseSummands(const std::string& out) {
    static const std::regex knotRe(
        R"regex(KnotSymbol\[(\d+),(\d+),(True|False),"([^"]+)"\]\s*->\s*(\d+))regex"
    );
    std::vector<KnotSummand> res;
    for (auto it = std::sregex_iterator(out.begin(), out.end(), knotRe),
              end = std::sregex_iterator();
         it != end; ++it) {
        KnotSummand s;
        s.c   = std::stoi((*it)[1].str());
        s.idx = std::stoi((*it)[2].str());
        s.alternating = ((*it)[3].str() == "True");
        s.prefix = chiralityPrefix((*it)[4].str());
        s.mult = std::stoi((*it)[5].str());
        res.push_back(s);
    }
    return res;
}

static std::string summandName(const KnotSummand& s) {
    std::string suffix = (s.c >= 11 && !s.alternating) ? "n" : "";
    return s.prefix + std::to_string(s.c) + "." +
           std::to_string(s.idx) + suffix;
}

// Compact name for the full (possibly composite) knot.
static std::string knotCompact(const std::vector<KnotSummand>& summands) {
    if (summands.empty()) return "a0.1";  // knoodleidentify's <||> case
    std::string out;
    for (const auto& s : summands) {
        for (int k = 0; k < s.mult; ++k) {
            if (!out.empty()) out += "#";
            out += summandName(s);
        }
    }
    return out;
}

// ---------- torus-knot recognition ----------

// Keyed by compact base name "<crossing>.<index>" (no chirality prefix,
// no "n" suffix). Extend as more constant-diff results become known.
static const std::unordered_map<std::string, std::string> TORUS_TABLE = {
    { "0.1",     "unknot"   },
    { "3.1",     "T(2,3)"   },
    { "5.1",     "T(2,5)"   },
    { "7.1",     "T(2,7)"   },
    { "8.19",    "T(3,4)"   },
    { "9.1",     "T(2,9)"   },
    { "10.124",  "T(3,5)"   },
};

// Look up a torus-knot label for a prime summand. Returns "?" if unknown.
static std::string torusOfSummand(const KnotSummand& s) {
    std::string base = std::to_string(s.c) + "." + std::to_string(s.idx);
    auto it = TORUS_TABLE.find(base);
    return (it == TORUS_TABLE.end()) ? "?" : it->second;
}

// Annotate a whole (possibly composite) knot with torus info. A knot is
// "torus" iff every prime summand is torus. Composites are usually NOT
// torus knots -- flag them but pass through the T(p,q) labels for
// context.
static std::string torusAnnotation(const std::vector<KnotSummand>& summands) {
    if (summands.empty()) return "unknot";
    std::vector<std::string> labels;
    bool anyUnknown = false;
    for (const auto& s : summands) {
        std::string t = torusOfSummand(s);
        if (t == "?") anyUnknown = true;
        for (int k = 0; k < s.mult; ++k) labels.push_back(t);
    }
    if (anyUnknown && labels.size() == 1) return "?";
    std::string joined;
    for (std::size_t i = 0; i < labels.size(); ++i) {
        if (i) joined += " # ";
        joined += labels[i];
    }
    if (labels.size() > 1) joined = "composite: " + joined;
    return joined;
}

// ---------- unresolved results ----------

// knoodleidentify reports knots it cannot name as NotFound[c, {...}] or
// Unidentified[c, {...}], where c is the crossing count of the diagram it
// got stuck on. Neither matches the KnotSymbol regex, so before this fix
// they fell through parseSummands as an empty summand list -- which
// knotCompact renders as "a0.1", the unknot. Every knot past the KLUT's
// range was therefore silently recorded as trivial, and then counted as
// confirming the torus conjecture. Detect them explicitly instead.
static bool isUnresolved(const std::string& out, std::string& label) {
    static const std::regex re(R"regex((NotFound|Unidentified)\[(\d+),)regex");
    std::smatch m;
    if (!std::regex_search(out, m, re)) return false;
    label = (m[1].str() == "NotFound" ? "notfound[" : "unidentified[")
          + m[2].str() + "]";
    return true;
}

// The genuine unknot is knoodleidentify's empty result "<||>".
static bool isTrueUnknot(const std::string& out) {
    return out.find("<||>") != std::string::npos;
}

// ---------- combinatorics + shell helpers ----------

static int gcd(int a, int b) {
    while (b) { a %= b; std::swap(a, b); }
    return a;
}

static std::vector<int> coprimeTo(int n) {
    std::vector<int> out;
    for (int d = 1; d < n; ++d)
        if (gcd(d, n) == 1) out.push_back(d);
    return out;
}

static std::vector<int> constDiffPerm(int n, int d) {
    std::vector<int> p(n);
    for (int i = 0; i < n; ++i) p[i] = (i * d) % n + 1;
    return p;
}

static std::string join(const std::vector<int>& v, char sep = ' ') {
    std::ostringstream ss;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (i) ss << sep;
        ss << v[i];
    }
    return ss.str();
}

static std::string runAndCapture(const std::string& cmd) {
    FILE* p = popen(cmd.c_str(), "r");
    if (!p) throw std::runtime_error("popen failed: " + cmd);
    std::string out;
    char buf[4096];
    while (std::fgets(buf, sizeof(buf), p)) out += buf;
    int status = pclose(p);
    if (status != 0) {
        throw std::runtime_error(
            "Command exited with status " + std::to_string(status) + ": " + cmd);
    }
    return out;
}

// Decide the torus type from the Alexander polynomial rather than from a
// name lookup. This works at any crossing count, including the cases the
// KLUT cannot name, and a negative answer is rigorous: deg Delta =
// (p-1)(q-1) leaves only finitely many candidate T(p,q), and torus_match
// checks them all.
static std::string torusViaAlexander(const std::string& tsv) {
    const std::string cmd =
        "knoodlesimplify --streaming-mode < " + tsv + " 2>/dev/null"
        " | ./identify_knot 2>/dev/null"
        " | ./torus_match --quiet 2>/dev/null";
    FILE* p = popen(cmd.c_str(), "r");
    if (!p) return "?";
    std::string out;
    char buf[512];
    while (std::fgets(buf, sizeof(buf), p)) out += buf;
    pclose(p);
    while (!out.empty() && (out.back() == '\n' || out.back() == '\r')) out.pop_back();
    return out.empty() ? "?" : out;
}

static void printUsage(const std::string& prog) {
    std::cout
        << "Usage:\n"
        << "  " << prog << " --n N [--d D] [--output FILE] [--petal-star-clean PATH]\n\n"
        << "Sweeps constant-difference permutations of length N (odd) with step d\n"
        << "coprime to N, identifies each knot, and annotates it with a torus-knot\n"
        << "label when the compact name matches TORUS_TABLE in the source.\n\n"
        << "Output columns (TSV): n, d, permutation, knot, torus_type.\n"
        << "torus_type is one of: 'unknot', 'T(p,q)', 'composite: ...', or '?'\n"
        << "when the compact name isn't in the table. A trailing summary line\n"
        << "reports how many rows were recognized as torus knots.\n";
}

int main(int argc, char* argv[]) {
    try {
        int n = -1;
        int single_d = -1;
        std::string outputFile;
        std::string petalBin = "./petal_star_clean";

        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "--help" || a == "-h") { printUsage(argv[0]); return 0; }
            else if (a == "--n") {
                if (i + 1 >= argc) throw std::invalid_argument("--n requires an integer");
                n = std::stoi(argv[++i]);
            }
            else if (a == "--d") {
                if (i + 1 >= argc) throw std::invalid_argument("--d requires an integer");
                single_d = std::stoi(argv[++i]);
            }
            else if (a == "--output") {
                if (i + 1 >= argc) throw std::invalid_argument("--output requires a filename");
                outputFile = argv[++i];
            }
            else if (a == "--petal-star-clean") {
                if (i + 1 >= argc) throw std::invalid_argument("--petal-star-clean requires a path");
                petalBin = argv[++i];
            }
            else throw std::invalid_argument("Unknown argument: " + a);
        }

        if (n <= 0) { printUsage(argv[0]); return 1; }
        if (n % 2 == 0) throw std::invalid_argument("n must be odd");
        if (n < 3) throw std::invalid_argument("n must be at least 3");

        std::vector<int> ds;
        if (single_d > 0) {
            if (gcd(single_d, n) != 1)
                throw std::invalid_argument("--d must be coprime to n");
            ds.push_back(single_d);
        } else {
            ds = coprimeTo(n);
        }

        std::ostream* out = &std::cout;
        std::ofstream ofs;
        if (!outputFile.empty()) {
            ofs.open(outputFile);
            if (!ofs) throw std::runtime_error("cannot open output: " + outputFile);
            out = &ofs;
        }

        *out << "# n=" << n << ", " << ds.size() << " permutation(s)\n";
        *out << "n\td\tpermutation\tknot\ttorus_type\n";

        const std::string tmp_tsv = "/tmp/sample_const_diff.tsv";
        int torusHits = 0, unknowns = 0, nonTorus = 0;

        for (int d : ds) {
            auto perm = constDiffPerm(n, d);
            std::string perm_str = join(perm);

            std::string psc_cmd =
                petalBin + " " + perm_str +
                " --output " + tmp_tsv + " > /dev/null";
            (void)runAndCapture(psc_cmd);

            std::string ki_cmd = "knoodleidentify " + tmp_tsv + " 2>/dev/null";
            std::string ki_out = runAndCapture(ki_cmd);

            auto summands = parseSummands(ki_out);

            std::string knot;
            std::string unresolvedLabel;
            if (!summands.empty()) {
                knot = knotCompact(summands);
            } else if (isTrueUnknot(ki_out)) {
                knot = "a0.1";
            } else if (isUnresolved(ki_out, unresolvedLabel)) {
                knot = unresolvedLabel;      // NOT the unknot -- just unnamed
            } else {
                knot = "parse-error";
            }

            // Torus type always comes from the Alexander polynomial, which
            // is defined regardless of whether the KLUT can name the knot.
            std::string torus = torusViaAlexander(tmp_tsv);
            if (torus == "?" || torus == "parse-error") ++unknowns;
            else if (torus == "no")                     ++nonTorus;
            else                                        ++torusHits;

            *out << n << '\t' << d << '\t' << perm_str << '\t'
                 << knot << '\t' << torus << '\n';
            out->flush();
        }

        *out << "# summary: " << torusHits << "/" << ds.size()
             << " torus or unknot, " << nonTorus
             << " PROVABLY NOT torus, " << unknowns << " undetermined\n";

        std::remove(tmp_tsv.c_str());
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
