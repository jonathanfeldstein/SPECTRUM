#!/usr/bin/env bash

set -e

DEBUG=OFF

for arg in "$@"; do
    case "$arg" in
        --debug)
            DEBUG=ON
            ;;
    esac
done

# Apple Silicon:
# Make OpenMP installed by Homebrew visible to CMake.
if [[ "$OSTYPE" == "darwin"* ]]; then
    export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
    export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"
    export CMAKE_PREFIX_PATH="/opt/homebrew/opt/libomp:/opt/homebrew/opt/jemalloc"
fi

# Conda / Micromamba:
# Make packages installed by conda visible to CMake.
if [[ -n "$CONDA_PREFIX" ]]; then
    export CMAKE_PREFIX_PATH="$CONDA_PREFIX:$CMAKE_PREFIX_PATH"
fi

rm -rf build

cmake \
    -S . \
    -B build \
    -DENABLE_DEBUG="$DEBUG"

cmake --build build --parallel