#!/bin/bash
# Build script for arc_index -- single-file, no external dependencies.

set -e

if [ "$1" = "debug" ]; then
    OPT_FLAGS="-g -O0 -DDEBUG"
    echo "Building arc_index (debug)..."
else
    OPT_FLAGS="-O2"
    echo "Building arc_index (release)..."
fi

g++ \
    -Wall \
    -Wextra \
    -std=c++17 \
    $OPT_FLAGS \
    -o arc_index \
    arc_index.cpp

echo "Build complete: ./arc_index"
