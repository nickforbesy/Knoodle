// sample_all_perms.cpp -- enumerate one representative per rotation orbit of
// petal permutations of length n (odd) -- i.e. all perms with pi_1 = 1 --
// in lexicographic order, identify each as a petal knot, and print
//
//     # 1 2 3 4 5 6 7 8 9 10 11
//     a0.1
//
//     # 1 2 3 4 5 6 7 8 9 11 10
//     a0.1
//
// Cyclic shifts of a petal permutation produce the same knot, so we only
// walk the (n-1)! reps with pi_1 = 1 instead of all n!.
//
// Warning: still combinatorial. For n=11, 10! = 3,628,800 reps at ~10 ms
// each is on the order of hours-to-a-day. Use small n (3, 5, 7, or 9)
// unless you have a plan for parallel / batched execution.
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

static std::string joinPerm(const std::vector<int>& v) {
    std::ostringstream ss;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (i) ss << ' ';
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

    if (summands.empty()) {
        // Unidentified / NotFound / Link / Invalid -- return a marker so
        // the row remains diagnosable rather than silently blank.
        return "?";
    }

    std::string result;
    for (std::size_t k = 0; k < summands.size(); ++k) {
        if (k) result += "#";
        result += summands[k];
    }
    return result;
}

static void printUsage(const std::string& prog) {
    std::cout
        << "Usage:\n"
        << "  " << prog << " --n N [--output FILE] [--petal-star-clean PATH]\n\n"
        << "Enumerates one representative per rotation orbit of petal\n"
        << "permutations of length N (odd) -- i.e. every perm with pi_1 = 1 --\n"
        << "in lex order, and identifies each as a petal knot. Output:\n\n"
        << "    # 1 2 3 4 5 6 7 8 9 10 11\n"
        << "    a0.1\n\n"
        << "    # 1 2 3 4 5 6 7 8 9 11 10\n"
        << "    a0.1\n\n"
        << "Composite knots print as e.g. 'p3.1#m3.1'. Requires KNOODLE_KLUT_DIR\n"
        << "to be set in the environment. Warning: (N-1)! still grows fast\n"
        << "(10! = 3.6M).\n";
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
        if (n < 3) throw std::invalid_argument("n must be at least 3");

        std::ostream* out = &std::cout;
        std::ofstream ofs;
        if (!outputFile.empty()) {
            ofs.open(outputFile);
            if (!ofs) throw std::runtime_error("cannot open output: " + outputFile);
            out = &ofs;
        }

        const std::string tmp_tsv = "/tmp/sample_all_perms.tsv";

        // Enumerate only rotation-orbit representatives: perms with pi_1 = 1.
        // Cyclic shifts of a petal permutation give the same knot, so one rep
        // per orbit is enough. That's (n-1)! perms instead of n!.
        std::vector<int> perm(n);
        perm[0] = 1;
        std::vector<int> tail(n - 1);
        std::iota(tail.begin(), tail.end(), 2);   // [2, 3, ..., n]

        bool first = true;
        do {
            for (int i = 0; i < n - 1; ++i) perm[i + 1] = tail[i];
            std::string perm_str = joinPerm(perm);

            std::string psc_cmd =
                petalBin + " " + perm_str +
                " --output " + tmp_tsv + " > /dev/null";
            (void)runAndCapture(psc_cmd);

            std::string ki_cmd = "knoodleidentify " + tmp_tsv + " 2>/dev/null";
            std::string ki_out = runAndCapture(ki_cmd);

            std::string knot = convertKnot(ki_out);

            if (!first) *out << "\n";  // blank line separator between entries only
            first = false;
            *out << "# " << perm_str << "\n";
            *out << knot << "\n";
            out->flush();
        } while (std::next_permutation(tail.begin(), tail.end()));

        std::remove(tmp_tsv.c_str());
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
