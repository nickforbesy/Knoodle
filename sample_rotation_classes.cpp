// sample_rotation_classes.cpp -- enumerate canonical representatives of the
// equivalence classes of petal permutations of odd length n under the
// combined action of cyclic rotation of position AND cyclic shift of value:
//
//     (pi_1, pi_2, ..., pi_n) ~ (pi_{i+a} + b, ...), for any a, b in Z_n.
//
// Two permutations are equivalent iff their cyclic-difference tuples
//
//     pd_i = (pi_{i+1} - pi_i) mod n,   pi_{n+1} := pi_1,
//
// are cyclic rotations of one another.
//
// Canonical representative of each class:
//   - pi_1 = 1
//   - pd tuple equals the lex-smallest cyclic rotation of the class's
//     pd cyclic orbit.
//
// A class is "robust" iff every pd_i is in {2, ..., n-2} -- equivalently,
// no cyclically adjacent pair of petals has consecutive height values.
// Robust rows are emitted first, then non-robust. Within each block rows
// are sorted by (winding_number, pd tuple), where winding_number = sum(pd)/n.
//
// Output columns (TSV): counter, pd, associated_permutation, knot_type.
//
// Requires petal_star_clean built at ./petal_star_clean (override with
// --petal-star-clean PATH) and knoodleidentify on PATH with
// KNOODLE_KLUT_DIR set in the environment.

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Perm = std::vector<int>;
using Pd   = std::vector<int>;

static std::string tupleFormat(const std::vector<int>& v) {
    std::ostringstream ss;
    ss << "(";
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (i) ss << ",";
        ss << v[i];
    }
    ss << ")";
    return ss.str();
}

static std::string joinPerm(const Perm& v) {
    std::ostringstream ss;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (i) ss << ' ';
        ss << v[i];
    }
    return ss.str();
}

static Pd computePd(const Perm& perm) {
    int n = (int)perm.size();
    Pd pd(n);
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        int d = (perm[j] - perm[i]) % n;
        if (d < 0) d += n;
        pd[i] = d;
    }
    return pd;
}

// True iff `pd` is lex-<= every cyclic rotation of itself, i.e. it is the
// canonical (lex-min) representative of its cyclic-rotation class.
static bool isLexMinRotation(const Pd& pd) {
    int n = (int)pd.size();
    for (int k = 1; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            int a = pd[i];
            int b = pd[(i + k) % n];
            if (a < b) break;
            if (a > b) return false;
        }
    }
    return true;
}

static std::string runAndCapture(const std::string& cmd) {
    FILE* p = popen(cmd.c_str(), "r");
    if (!p) throw std::runtime_error("popen failed: " + cmd);
    std::string out;
    char buf[4096];
    while (std::fgets(buf, sizeof(buf), p)) out += buf;
    int status = pclose(p);
    if (status != 0)
        throw std::runtime_error("Command exited with status " +
                                 std::to_string(status) + ": " + cmd);
    return out;
}

// Convert a symmetry-coset string like "e/r" or "m/mr" or "e/r/m/mr" to
// the single-letter chirality prefix used in the compact notation.
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
    if (hasE)         return "m";  // e-coset -> Parsley "m"
    if (hasM)         return "p";  // m-coset -> Parsley "p"
    return "h";                    // rare (chiral, non-invertible cosets)
}

// Convert knoodleidentify's default output (a Wolfram-language association)
// into the compact "a0.1" / "p3.1" / "m5.2" / "p3.1#m3.1" style.
static std::string convertKnot(const std::string& out) {
    if (out.find("<||>") != std::string::npos) return "a0.1";

    static const std::regex knotRe(
        R"regex(KnotSymbol\[(\d+),(\d+),(True|False),"([^"]+)"\]\s*->\s*(\d+))regex"
    );

    std::vector<std::string> summands;
    auto begin = std::sregex_iterator(out.begin(), out.end(), knotRe);
    auto end   = std::sregex_iterator();
    for (auto it = begin; it != end; ++it) {
        int c   = std::stoi((*it)[1].str());
        int idx = std::stoi((*it)[2].str());
        std::string alt = (*it)[3].str();
        std::string sym = (*it)[4].str();
        int mult = std::stoi((*it)[5].str());

        std::string prefix = chiralityPrefix(sym);
        std::string suffix = (c >= 11 && alt == "False") ? "n" : "";
        std::string name = prefix + std::to_string(c) + "."
                           + std::to_string(idx) + suffix;
        for (int k = 0; k < mult; ++k) summands.push_back(name);
    }

    if (summands.empty()) return "?";

    std::string result;
    for (std::size_t k = 0; k < summands.size(); ++k) {
        if (k) result += "#";
        result += summands[k];
    }
    return result;
}

struct Row {
    bool robust;
    int  winding;
    Pd   pd;
    Perm perm;
};

static void printUsage(const std::string& prog) {
    std::cout
        << "Usage:\n"
        << "  " << prog << " --n N [--output FILE] [--petal-star-clean PATH]\n\n"
        << "Enumerates canonical representatives of the equivalence classes of\n"
        << "petal permutations of length N (odd) under the combined action of\n"
        << "cyclic rotation of position AND cyclic shift of value. Each class is\n"
        << "represented by the permutation starting with 1 whose cyclic-difference\n"
        << "tuple is the lex-smallest rotation of its cyclic-difference orbit.\n\n"
        << "Robust rows (all pd_i in {2, ..., N-2}) are emitted first; non-robust\n"
        << "follow. Within each block rows are sorted by (winding_number, pd tuple).\n\n"
        << "Output columns (TSV): counter, pd, associated_permutation, knot_type.\n\n"
        << "Requires KNOODLE_KLUT_DIR to be set in the environment. Warning:\n"
        << "the enumeration walks all (N-1)! perms with pi_1=1 to find canonical\n"
        << "reps; for N=11 this is 3.6M candidates.\n";
}

int main(int argc, char* argv[]) {
    try {
        int n = -1;
        std::string outputFile;
        std::string petalBin = "./petal_star_clean";

        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "--help" || a == "-h") { printUsage(argv[0]); return 0; }
            else if (a == "--n") {
                if (i + 1 >= argc) throw std::invalid_argument("--n requires an integer");
                n = std::stoi(argv[++i]);
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
        if (n < 3)      throw std::invalid_argument("n must be at least 3");

        std::ostream* out = &std::cout;
        std::ofstream ofs;
        if (!outputFile.empty()) {
            ofs.open(outputFile);
            if (!ofs) throw std::runtime_error("cannot open output: " + outputFile);
            out = &ofs;
        }

        // Walk all perms with pi_1 = 1 (i.e., all perms of {2..n} in
        // positions 2..n). Keep only those whose pd tuple is the lex-min
        // rotation of its cyclic orbit -- those are the canonical class reps.
        std::vector<Row> rows;
        Perm perm(n);
        perm[0] = 1;
        std::vector<int> tail(n - 1);
        std::iota(tail.begin(), tail.end(), 2);   // {2, 3, ..., n}

        do {
            for (int i = 0; i < n - 1; ++i) perm[i + 1] = tail[i];
            Pd pd = computePd(perm);
            if (!isLexMinRotation(pd)) continue;

            int wsum = 0;
            for (int d : pd) wsum += d;
            int winding = wsum / n;

            bool robust = true;
            for (int d : pd) if (d == 1 || d == n - 1) { robust = false; break; }

            rows.push_back({robust, winding, pd, perm});
        } while (std::next_permutation(tail.begin(), tail.end()));

        // Robust first, then by winding, then by pd tuple lex.
        std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) {
            if (a.robust != b.robust) return a.robust > b.robust;
            if (a.winding != b.winding) return a.winding < b.winding;
            return a.pd < b.pd;
        });

        *out << "counter\tpd\tassociated_permutation\tknot_type\n";

        const std::string tmp_tsv = "/tmp/sample_rotation_classes.tsv";
        int counter = 0;

        for (const Row& r : rows) {
            ++counter;
            std::string perm_str = joinPerm(r.perm);
            std::string psc_cmd =
                petalBin + " " + perm_str +
                " --output " + tmp_tsv + " > /dev/null";
            (void)runAndCapture(psc_cmd);

            std::string ki_cmd = "knoodleidentify " + tmp_tsv + " 2>/dev/null";
            std::string ki_out = runAndCapture(ki_cmd);
            std::string knot   = convertKnot(ki_out);

            *out << counter << '\t'
                 << tupleFormat(r.pd) << '\t'
                 << tupleFormat(r.perm) << '\t'
                 << knot << '\n';
            out->flush();
        }

        std::remove(tmp_tsv.c_str());
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
