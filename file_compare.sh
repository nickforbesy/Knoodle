#!/bin/bash
# Build script for sample_rotation_classes.

set -e

if [ "$1" = "debug" ]; then
    OPT_FLAGS="-g -O0 -DDEBUG"
    echo "Building sample_rotation_classes (debug)..."
else
    OPT_FLAGS="-O2"
    echo "Building sample_rotation_classes (release)..."
fi

g++ \
    -Wall \
    -Wextra \
    -std=c++17 \
    $OPT_FLAGS \
    -o file_compare \
    file_compare.cpp

echo "Build complete: ./file_compare"