#!/bin/bash
# configure_and_build.sh
# Install dependencies and build KratosMultiphysics with MPMApplication and
# LinearSolversApplication on Ubuntu/Debian Linux (Build from Source – developer
# approach, as described in applications/MPMApplication/README.md and INSTALL.md).
#
# Quick alternative for users (no compilation required):
#   pip3 install KratosMPMApplication matplotlib meshio
#
# Usage (developer build):
#   cd /path/to/Kratos3.3
#   bash examples/sediment_dumping/configure_and_build.sh

set -e

# --------------------------------------------------------------------------
# 1. Install system dependencies
# --------------------------------------------------------------------------
sudo apt-get update -y
sudo apt-get install -y \
    python3-dev \
    gcc \
    g++ \
    cmake \
    libboost-all-dev \
    python3-pip \
    python3-numpy \
    python3-matplotlib

# Post-processing dependencies (matplotlib, meshio)
pip3 install --break-system-packages matplotlib meshio

# --------------------------------------------------------------------------
# 2. Environment variables
# --------------------------------------------------------------------------
add_app () {
    export KRATOS_APPLICATIONS="${KRATOS_APPLICATIONS}$1;"
}

export CC=gcc
export CXX=g++
export KRATOS_SOURCE="${KRATOS_SOURCE:-"$( cd "$(dirname "$0")/../.." ; pwd -P )"}"
export KRATOS_BUILD="${KRATOS_SOURCE}/build"
export KRATOS_APP_DIR="${KRATOS_SOURCE}/applications"
export KRATOS_BUILD_TYPE="Release"
export PYTHON_EXECUTABLE="/usr/bin/python3"

export KRATOS_APPLICATIONS=
add_app "${KRATOS_APP_DIR}/MPMApplication"
add_app "${KRATOS_APP_DIR}/LinearSolversApplication"

echo "KRATOS_SOURCE       = ${KRATOS_SOURCE}"
echo "KRATOS_BUILD        = ${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}"
echo "KRATOS_APPLICATIONS = ${KRATOS_APPLICATIONS}"

# --------------------------------------------------------------------------
# 3. Clean previous CMake cache
# --------------------------------------------------------------------------
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/cmake_install.cmake"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeCache.txt"
rm -rf "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}/CMakeFiles"

# --------------------------------------------------------------------------
# 4. Configure
# --------------------------------------------------------------------------
cmake \
    -H"${KRATOS_SOURCE}" \
    -B"${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" \
    -DUSE_MPI=OFF \
    -DUSE_EIGEN_MKL=OFF \
    -DPYTHON_EXECUTABLE="${PYTHON_EXECUTABLE}" \
    -DCMAKE_BUILD_TYPE="${KRATOS_BUILD_TYPE}"

# --------------------------------------------------------------------------
# 5. Build and install
# --------------------------------------------------------------------------
cmake \
    --build "${KRATOS_BUILD}/${KRATOS_BUILD_TYPE}" \
    --target install \
    -- -j"$(nproc)"

echo ""
echo "Build complete. Kratos is installed under ${KRATOS_SOURCE}/bin."
echo ""
echo "To run the sediment dumping simulation:"
echo "  cd ${KRATOS_SOURCE}/examples/sediment_dumping"
echo "  python3 run_simulation.py"
