#!/bin/bash
# Build script for sample_all_perms.

set -e

if [ "$1" = "debug" ]; then
    OPT_FLAGS="-g -O0 -DDEBUG"
    echo "Building sample_all_perms (debug)..."
else
    OPT_FLAGS="-O2"
    echo "Building sample_all_perms (release)..."
fi

g++ \
    -Wall \
    -Wextra \
    -std=c++17 \
    $OPT_FLAGS \
    -o sample_all_perms \
    sample_all_perms.cpp

echo "Build complete: ./sample_all_perms"
