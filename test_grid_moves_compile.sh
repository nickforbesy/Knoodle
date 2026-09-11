#!/bin/bash
# Build script for test_grid_moves -- single-file, no external dependencies.

set -e

g++ -Wall -Wextra -std=c++17 -O2 -I. -o test_grid_moves test_grid_moves.cpp

echo "Build complete: ./test_grid_moves  (run as ./test_grid_moves 7)"
