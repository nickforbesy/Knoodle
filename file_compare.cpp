// file_compare -- compare two text files after stripping all whitespace.
//
// Two files are considered identical iff their character streams agree once
// every whitespace character (spaces, tabs, newlines, CR, etc.) is removed.
// This ignores blank-line differences, trailing newlines, indentation, and
// the fact that one file may be a strict prefix of the other only when it
// really is a prefix -- we still require the whitespace-free streams to
// have equal length.

#include <cctype>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

// Advance `it` to the next non-whitespace byte (or end).
static void skipWhitespace(std::istreambuf_iterator<char>& it,
                           const std::istreambuf_iterator<char>& end) {
    while (it != end && std::isspace(static_cast<unsigned char>(*it))) {
        ++it;
    }
}

bool areFilesIdentical(const std::string& file1, const std::string& file2,
                       std::size_t& firstMismatchByte) {
    std::ifstream f1(file1, std::ios::binary);
    std::ifstream f2(file2, std::ios::binary);

    if (!f1.is_open()) throw std::runtime_error("Failed to open " + file1);
    if (!f2.is_open()) throw std::runtime_error("Failed to open " + file2);

    std::istreambuf_iterator<char> it1(f1), it2(f2), end;
    std::size_t pos = 0;

    while (true) {
        skipWhitespace(it1, end);
        skipWhitespace(it2, end);

        if (it1 == end && it2 == end) return true;
        if (it1 == end || it2 == end) {
            firstMismatchByte = pos;
            return false;
        }
        if (*it1 != *it2) {
            firstMismatchByte = pos;
            return false;
        }
        ++it1;
        ++it2;
        ++pos;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " file1 file2\n";
        return 1;
    }

    try {
        std::size_t mismatch = 0;
        if (areFilesIdentical(argv[1], argv[2], mismatch)) {
            std::cout << "The files are identical (ignoring whitespace).\n";
            return 0;
        } else {
            std::cout << "The files differ at non-whitespace character #"
                      << mismatch << ".\n";
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
