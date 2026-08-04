#!/bin/bash
# Build script for sample_const_diff -- single-file, no external dependencies.

set -e

if [ "$1" = "debug" ]; then
    OPT_FLAGS="-g -O0 -DDEBUG"
    echo "Building sample_const_diff (debug)..."
else
    OPT_FLAGS="-O2"
    echo "Building sample_const_diff (release)..."
fi

g++ \
    -Wall \
    -Wextra \
    -std=c++17 \
    $OPT_FLAGS \
    -o sample_const_diff \
    sample_const_diff.cpp

echo "Build complete: ./sample_const_diff"
