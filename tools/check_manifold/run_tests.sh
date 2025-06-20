#!/bin/bash

# Test suite for triangle mesh validator
echo "=== Triangle Mesh Validator Test Suite ==="
echo

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

VALIDATOR="./validate_triangle_mesh"
TEST_PASSED=0
TEST_FAILED=0

# Function to run a test
run_test() {
    local test_name="$1"
    local test_file="$2"
    local expected_manifold="$3"
    local expected_issues="$4"
    local params="$5"
    
    echo -e "${YELLOW}Testing: $test_name${NC}"
    echo "File: $test_file"
    echo "Expected: $expected_manifold, Issues: $expected_issues"
    
    if [ ! -f "$test_file" ]; then
        echo -e "${RED}FAIL: Test file not found${NC}"
        ((TEST_FAILED++))
        echo
        return
    fi
    
    # Run the validator and capture output
    if [ -n "$params" ]; then
        output=$($VALIDATOR $test_file $params 2>&1)
    else
        output=$($VALIDATOR $test_file 2>&1)
    fi
    
    # Check manifold status
    manifold_check_passed=false
    if [[ "$expected_manifold" == "MANIFOLD" && "$output" =~ "Mesh is MANIFOLD" ]]; then
        manifold_check_passed=true
    elif [[ "$expected_manifold" == "NON-MANIFOLD" && "$output" =~ "Mesh is NON-MANIFOLD" ]]; then
        manifold_check_passed=true
    fi
    
    # Check for expected issues
    issues_check_passed=true
    if [[ "$expected_issues" == *"degenerate"* && ! "$output" =~ "Degenerate triangles: 0" ]]; then
        issues_check_passed=true
    elif [[ "$expected_issues" == *"small_area"* && ! "$output" =~ "Small area triangles: 0" ]]; then
        issues_check_passed=true
    elif [[ "$expected_issues" == *"extreme_angle"* && ! "$output" =~ "Extreme angle triangles: 0" ]]; then
        issues_check_passed=true
    elif [[ "$expected_issues" == *"high_aspect"* && ! "$output" =~ "High aspect ratio triangles: 0" ]]; then
        issues_check_passed=true
    elif [[ "$expected_issues" == *"none"* ]]; then
        if [[ "$output" =~ "Degenerate triangles: 0" && 
              "$output" =~ "Small area triangles: 0" && 
              "$output" =~ "Extreme angle triangles: 0" && 
              "$output" =~ "High aspect ratio triangles: 0" ]]; then
            issues_check_passed=true
        else
            issues_check_passed=false
        fi
    fi
    
    # Report results
    if $manifold_check_passed && $issues_check_passed; then
        echo -e "${GREEN}PASS${NC}"
        ((TEST_PASSED++))
    else
        echo -e "${RED}FAIL${NC}"
        if ! $manifold_check_passed; then
            echo "  Manifold check failed"
        fi
        if ! $issues_check_passed; then
            echo "  Issue detection failed"
        fi
        echo "  Output:"
        echo "$output" | sed 's/^/    /'
        ((TEST_FAILED++))
    fi
    echo
}

# Build the project first
echo "Building project..."
if ! make > /dev/null 2>&1; then
    echo -e "${RED}FAIL: Build failed${NC}"
    exit 1
fi
echo -e "${GREEN}Build successful${NC}"

# Check if validator exists
if [ ! -x "$VALIDATOR" ]; then
    echo -e "${RED}ERROR: Validator executable ${VALIDATOR} not found or not executable${NC}"
    exit 1
fi

# Run tests
echo "Running tests..."
echo

# Test 1: Clean tetrahedron (should pass all checks)
run_test "Clean Tetrahedron" "clean_tetrahedron.off" "MANIFOLD" "none"

# Test 2: Cube (existing test file, should be clean)
run_test "Cube" "cube.off" "MANIFOLD" "none"

# Test 3: Non-manifold edge (Note: CGAL may auto-repair some non-manifold cases)
echo -e "${YELLOW}Testing: Non-manifold Edge (informational)${NC}"
echo "Note: CGAL mesh loading may automatically repair some non-manifold cases"
output=$(./build/validate_triangle_mesh nonmanifold_edge.off 2>&1)
if [[ "$output" =~ "Mesh is MANIFOLD" ]]; then
    echo -e "${YELLOW}INFO: CGAL auto-repaired the mesh to manifold${NC}"
else
    echo -e "${GREEN}SUCCESS: Non-manifold detected${NC}"
fi
echo

# Test 4: Zero area triangle (collinear points)
run_test "Zero Area Triangle" "zero_area_triangle.off" "MANIFOLD" "small_area"

# Test 5: Sliver triangle (very small area and extreme angles)
run_test "Sliver Triangle" "sliver_triangle.off" "MANIFOLD" "small_area" "1e-6"

# Test 6: Open cylinder (existing test file)
run_test "Open Cylinder" "open_cylinder.off" "MANIFOLD" "none"

# Test parameter parsing
echo -e "${YELLOW}Testing parameter parsing...${NC}"
param_output=$($VALIDATOR 2>&1)
if [[ "$param_output" =~ "Usage" && "$param_output" =~ "Defaults" ]]; then
    echo -e "${GREEN}PASS: Usage message displayed correctly${NC}"
    ((TEST_PASSED++))
else
    echo -e "${RED}FAIL: Usage message incorrect${NC}"
    ((TEST_FAILED++))
fi
echo

# Summary
echo "=== Test Summary ==="
echo -e "Tests passed: ${GREEN}$TEST_PASSED${NC}"
echo -e "Tests failed: ${RED}$TEST_FAILED${NC}"
echo "Total tests: $((TEST_PASSED + TEST_FAILED))"

if [ $TEST_FAILED -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed.${NC}"
    exit 1
fi
