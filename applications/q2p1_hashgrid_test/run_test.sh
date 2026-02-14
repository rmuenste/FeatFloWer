#!/bin/bash
# HashGrid Test Runner
# This script helps automate the testing process

set -e  # Exit on error

REPO_ROOT="/data/warehouse17/rmuenste/code/FF-ATC-NEW/ff-accel"
NP=4  # Number of MPI processes

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=========================================="
echo "HashGrid Acceleration Test Runner"
echo "=========================================="
echo ""

# Check if we're in the right directory
if [ ! -f "q2p1_hashgrid_test.f90" ]; then
    echo -e "${RED}Error: Not in test application directory${NC}"
    echo "Run this script from: applications/q2p1_hashgrid_test/"
    exit 1
fi

# Function to build version
build_version() {
    local version=$1
    local accel_flag=$2
    local build_dir="${REPO_ROOT}/build_test_${version}"

    echo -e "${YELLOW}Building ${version} version...${NC}"

    # Create build directory
    mkdir -p "${build_dir}"
    cd "${build_dir}"

    # Configure
    echo "  Configuring with CMake..."
    cmake -DCMAKE_BUILD_TYPE=Release \
          -DUSE_PE=ON \
          -DBUILD_APPLICATIONS=ON \
          .. > /dev/null 2>&1

    cmake -DUSE_PE_SERIAL_MODE=ON \
          -DUSE_ACCELERATED_POINT_QUERY=${accel_flag} \
          .. > /dev/null 2>&1

    # Build
    echo "  Building q2p1_hashgrid_test..."
    make -j8 q2p1_hashgrid_test > /dev/null 2>&1

    if [ $? -eq 0 ]; then
        echo -e "  ${GREEN}✓ Build successful${NC}"
    else
        echo -e "  ${RED}✗ Build failed${NC}"
        exit 1
    fi
}

# Function to run test
run_test() {
    local version=$1
    local build_dir="${REPO_ROOT}/build_test_${version}"
    local output_file="${REPO_ROOT}/output_${version}.txt"

    echo -e "${YELLOW}Running ${version} test...${NC}"

    cd "${build_dir}/applications/q2p1_hashgrid_test"

    if [ ! -f "./q2p1_hashgrid_test" ]; then
        echo -e "  ${RED}✗ Executable not found${NC}"
        exit 1
    fi

    # Run test
    mpirun -np ${NP} ./q2p1_hashgrid_test > "${output_file}" 2>&1

    if [ $? -eq 0 ]; then
        echo -e "  ${GREEN}✓ Test completed${NC}"
    else
        echo -e "  ${RED}✗ Test failed${NC}"
        echo "  See output in: ${output_file}"
        exit 1
    fi
}

# Function to extract count
extract_count() {
    local file=$1
    grep "Total dofs inside" "${file}" | head -1 | awk '{print $5}'
}

# Main test sequence
main() {
    echo "Step 1: Building baseline version (no acceleration)"
    build_version "baseline" "OFF"
    echo ""

    echo "Step 2: Building accelerated version"
    build_version "accelerated" "ON"
    echo ""

    echo "Step 3: Running baseline test"
    run_test "baseline"
    echo ""

    echo "Step 4: Running accelerated test"
    run_test "accelerated"
    echo ""

    echo "=========================================="
    echo "Results Comparison"
    echo "=========================================="

    baseline_count=$(extract_count "${REPO_ROOT}/output_baseline.txt")
    accelerated_count=$(extract_count "${REPO_ROOT}/output_accelerated.txt")

    echo "Baseline count:    ${baseline_count}"
    echo "Accelerated count: ${accelerated_count}"
    echo ""

    if [ "${baseline_count}" == "${accelerated_count}" ]; then
        echo -e "${GREEN}✓ TEST PASSED: Counts match!${NC}"
        echo ""
        echo "HashGrid acceleration produces identical results to baseline."
        exit 0
    else
        echo -e "${RED}✗ TEST FAILED: Counts differ!${NC}"
        echo ""
        echo "There is a discrepancy in the point-inside query results."
        echo "Check the output files for details:"
        echo "  - ${REPO_ROOT}/output_baseline.txt"
        echo "  - ${REPO_ROOT}/output_accelerated.txt"
        exit 1
    fi
}

# Parse command line options
case "${1:-}" in
    build-only)
        echo "Building both versions only..."
        build_version "baseline" "OFF"
        build_version "accelerated" "ON"
        echo -e "${GREEN}Build complete${NC}"
        ;;
    baseline-only)
        echo "Running baseline test only..."
        run_test "baseline"
        echo -e "${GREEN}Baseline test complete${NC}"
        echo "Count: $(extract_count ${REPO_ROOT}/output_baseline.txt)"
        ;;
    accelerated-only)
        echo "Running accelerated test only..."
        run_test "accelerated"
        echo -e "${GREEN}Accelerated test complete${NC}"
        echo "Count: $(extract_count ${REPO_ROOT}/output_accelerated.txt)"
        ;;
    compare-only)
        echo "Comparing existing results..."
        baseline_count=$(extract_count "${REPO_ROOT}/output_baseline.txt")
        accelerated_count=$(extract_count "${REPO_ROOT}/output_accelerated.txt")
        echo "Baseline:    ${baseline_count}"
        echo "Accelerated: ${accelerated_count}"
        [ "${baseline_count}" == "${accelerated_count}" ] && echo -e "${GREEN}✓ Match${NC}" || echo -e "${RED}✗ Differ${NC}"
        ;;
    help|--help|-h)
        echo "Usage: $0 [option]"
        echo ""
        echo "Options:"
        echo "  (none)           - Run complete test sequence"
        echo "  build-only       - Build both versions only"
        echo "  baseline-only    - Run baseline test only"
        echo "  accelerated-only - Run accelerated test only"
        echo "  compare-only     - Compare existing results"
        echo "  help             - Show this help"
        ;;
    *)
        main
        ;;
esac
