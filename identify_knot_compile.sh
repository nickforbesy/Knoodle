#!/bin/bash
# Build script for identify_knot -- single-file, no external dependencies.

set -e

if [ "$1" = "debug" ]; then
    OPT_FLAGS="-g -O0 -DDEBUG"
    echo "Building identify_knot (debug)..."
else
    OPT_FLAGS="-O2"
    echo "Building identify_knot (release)..."
fi

g++ \
    -Wall \
    -Wextra \
    -std=c++17 \
    $OPT_FLAGS \
    -o identify_knot \
    identify_knot.cpp

echo "Build complete: ./identify_knot"
