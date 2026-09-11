#pragma once

// Shared helpers for turning knoodleidentify's Wolfram-language output into
// the compact names described in thesis Section "Compact knot notation".
// Extracted so tools other than sample_const_diff.cpp can reuse them; the
// conventions (including the Parsley chirality swap) are identical.

#include <cstdio>
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>

namespace knotname {

inline std::string chiralityPrefix(const std::string& sym) {
    bool hasE = false, hasM = false;
    std::size_t start = 0;
    while (start < sym.size()) {
        std::size_t end = sym.find('/', start);
        if (end == std::string::npos) end = sym.size();
        const std::string tok = sym.substr(start, end - start);
        if (tok == "e")      hasE = true;
        else if (tok == "m") hasM = true;
        start = end + 1;
    }
    if (hasE && hasM) return "a";
    // Parsley uses the opposite mirror convention from knoodleidentify; the
    // two branches are swapped to match. See sample_const_diff.cpp.
    if (hasE) return "m";
    if (hasM) return "p";
    return "h";
}

struct Summand {
    int c = 0;
    int idx = 0;
    bool alternating = true;
    std::string prefix;
    int mult = 1;
};

inline std::vector<Summand> parseSummands(const std::string& out) {
    static const std::regex knotRe(
        R"regex(KnotSymbol\[(\d+),(\d+),(True|False),"([^"]+)"\]\s*->\s*(\d+))regex");
    std::vector<Summand> res;
    for (auto it = std::sregex_iterator(out.begin(), out.end(), knotRe),
              end = std::sregex_iterator(); it != end; ++it) {
        Summand s;
        s.c           = std::stoi((*it)[1].str());
        s.idx         = std::stoi((*it)[2].str());
        s.alternating = ((*it)[3].str() == "True");
        s.prefix      = chiralityPrefix((*it)[4].str());
        s.mult        = std::stoi((*it)[5].str());
        res.push_back(s);
    }
    return res;
}

inline std::string summandName(const Summand& s) {
    const std::string suffix = (s.c >= 11 && !s.alternating) ? "n" : "";
    return s.prefix + std::to_string(s.c) + "." + std::to_string(s.idx) + suffix;
}

inline std::string compactName(const std::vector<Summand>& summands) {
    if (summands.empty()) return "a0.1";
    std::string out;
    for (const auto& s : summands)
        for (int k = 0; k < s.mult; ++k) {
            if (!out.empty()) out += "#";
            out += summandName(s);
        }
    return out;
}

inline bool isUnresolved(const std::string& out, std::string& label) {
    static const std::regex re(R"regex((NotFound|Unidentified)\[(\d+),)regex");
    std::smatch m;
    if (!std::regex_search(out, m, re)) return false;
    label = (m[1].str() == "NotFound" ? "notfound[" : "unidentified[") + m[2].str() + "]";
    return true;
}

inline bool isTrueUnknot(const std::string& out) {
    return out.find("<||>") != std::string::npos;
}

inline std::string runAndCapture(const std::string& cmd) {
    FILE* p = popen(cmd.c_str(), "r");
    if (!p) throw std::runtime_error("popen failed: " + cmd);
    std::string out;
    char buf[4096];
    while (std::fgets(buf, sizeof(buf), p)) out += buf;
    pclose(p);
    return out;
}

// Identify the knot described by a polygon TSV. Returns the compact name and
// fills `summands` with the parsed prime factors (empty for the unknot).
inline std::string identifyTSV(const std::string& tsv,
                               std::vector<Summand>& summands) {
    const std::string out = runAndCapture("knoodleidentify " + tsv + " 2>/dev/null");
    summands = parseSummands(out);
    if (!summands.empty())      return compactName(summands);
    if (isTrueUnknot(out))      return "a0.1";
    std::string label;
    if (isUnresolved(out, label)) return label;
    return "parse-error";
}

}  // namespace knotname
