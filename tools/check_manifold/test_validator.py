#!/usr/bin/env python3
"""
Test suite for the triangle mesh validator.
This script runs various test cases and verifies the output.
"""

import subprocess
import sys
import os
from dataclasses import dataclass
from typing import List, Optional

@dataclass
class TestCase:
    name: str
    filename: str
    expected_manifold: bool
    expected_degenerate: int
    expected_small_area: int
    expected_extreme_angles: int
    params: Optional[List[str]] = None

class TestRunner:
    def __init__(self, validator_path: str):
        self.validator_path = validator_path
        self.passed = 0
        self.failed = 0
        
    def run_validator(self, filename: str, params: List[str] = None) -> str:
        """Run the validator and return output."""
        cmd = [self.validator_path, filename]
        if params:
            cmd.extend(params)
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            return result.stdout + result.stderr
        except subprocess.TimeoutExpired:
            return "ERROR: Validator timed out"
        except Exception as e:
            return f"ERROR: {str(e)}"
    
    def parse_output(self, output: str) -> dict:
        """Parse validator output to extract key metrics."""
        metrics = {
            'manifold': None,
            'degenerate': 0,
            'small_area': 0,
            'extreme_angles': 0,
            'high_aspect': 0,
            'high_circumradius': 0
        }
        
        lines = output.split('\n')
        for line in lines:
            line = line.strip()
            if 'Final verdict: Mesh is MANIFOLD' in line:
                metrics['manifold'] = True
            elif 'Final verdict: Mesh is NON-MANIFOLD' in line:
                metrics['manifold'] = False
            elif 'Degenerate triangles:' in line:
                metrics['degenerate'] = int(line.split(':')[-1].strip())
            elif 'Small area triangles:' in line:
                metrics['small_area'] = int(line.split(':')[-1].strip())
            elif 'Extreme angle triangles:' in line:
                metrics['extreme_angles'] = int(line.split(':')[-1].strip())
            elif 'High aspect ratio triangles:' in line:
                metrics['high_aspect'] = int(line.split(':')[-1].strip())
            elif 'High circumradius/inradius triangles:' in line:
                metrics['high_circumradius'] = int(line.split(':')[-1].strip())
        
        return metrics
    
    def run_test(self, test_case: TestCase) -> bool:
        """Run a single test case."""
        print(f"Testing: {test_case.name}")
        print(f"File: {test_case.filename}")
        
        if not os.path.exists(test_case.filename):
            print(f"❌ FAIL: Test file {test_case.filename} not found")
            return False
        
        output = self.run_validator(test_case.filename, test_case.params)
        if output.startswith("ERROR:"):
            print(f"❌ FAIL: {output}")
            return False
        
        metrics = self.parse_output(output)
        
        # Check results
        success = True
        
        if metrics['manifold'] != test_case.expected_manifold:
            print(f"❌ Manifold check failed: expected {test_case.expected_manifold}, got {metrics['manifold']}")
            success = False
        
        if metrics['degenerate'] != test_case.expected_degenerate:
            print(f"❌ Degenerate count failed: expected {test_case.expected_degenerate}, got {metrics['degenerate']}")
            success = False
            
        if metrics['small_area'] != test_case.expected_small_area:
            print(f"❌ Small area count failed: expected {test_case.expected_small_area}, got {metrics['small_area']}")
            success = False
            
        if metrics['extreme_angles'] != test_case.expected_extreme_angles:
            print(f"❌ Extreme angles count failed: expected {test_case.expected_extreme_angles}, got {metrics['extreme_angles']}")
            success = False
        
        if success:
            print("✅ PASS")
            self.passed += 1
        else:
            print("❌ FAIL")
            self.failed += 1
            print("Output:")
            print(output.replace('\n', '\n  '))
        
        print()
        return success
    
    def run_all_tests(self):
        """Run all test cases."""
        print("=== Triangle Mesh Validator Test Suite ===\n")
        
        # Build if needed
        if not os.path.exists(self.validator_path):
            print("Building project...")
            try:
                subprocess.run(["cmake", ".", "-B", "build"], check=True, capture_output=True)
                subprocess.run(["make", "-C", "build"], check=True, capture_output=True)
                print("✅ Build successful\n")
            except subprocess.CalledProcessError as e:
                print(f"❌ Build failed: {e}")
                return False
        
        # Define test cases
        test_cases = [
            TestCase(
                name="Clean Tetrahedron",
                filename="test_files/clean_tetrahedron.off",
                expected_manifold=True,
                expected_degenerate=0,
                expected_small_area=0,
                expected_extreme_angles=0
            ),
            TestCase(
                name="Cube",
                filename="test_files/cube.off",
                expected_manifold=True,
                expected_degenerate=0,
                expected_small_area=0,
                expected_extreme_angles=0
            ),
            TestCase(
                name="Zero Area Triangle",
                filename="test_files/zero_area_triangle.off",
                expected_manifold=True,
                expected_degenerate=0,
                expected_small_area=1,
                expected_extreme_angles=1
            ),
            TestCase(
                name="Sliver Triangle (strict params)",
                filename="test_files/sliver_triangle.off",
                expected_manifold=True,
                expected_degenerate=0,
                expected_small_area=1,
                expected_extreme_angles=1,
                params=["1e-6", "1", "179", "30", "10"]
            ),
            TestCase(
                name="Non-manifold Edge (should be detected)",
                filename="test_files/nonmanifold_edge.off",
                expected_manifold=False,
                expected_degenerate=0,
                expected_small_area=0,
                expected_extreme_angles=0
            ),
            TestCase(
                name="Open Cylinder",
                filename="test_files/open_cylinder.off",
                expected_manifold=True,
                expected_degenerate=0,
                expected_small_area=0,
                expected_extreme_angles=0
            )
        ]
        
        # Run tests
        for test_case in test_cases:
            self.run_test(test_case)
        
        # Test parameter parsing
        print("Testing parameter parsing...")
        output = subprocess.run([self.validator_path], capture_output=True, text=True).stderr
        if "Usage" in output and "Defaults" in output:
            print("✅ PASS: Usage message displayed correctly")
            self.passed += 1
        else:
            print("❌ FAIL: Usage message incorrect")
            self.failed += 1
        print()
        
        # Summary
        print("=== Test Summary ===")
        print(f"Tests passed: {self.passed}")
        print(f"Tests failed: {self.failed}")
        print(f"Total tests: {self.passed + self.failed}")
        
        if self.failed == 0:
            print("🎉 All tests passed!")
            return True
        else:
            print("❌ Some tests failed.")
            return False

def main():
    validator_path = "validate_triangle_mesh"
    
    runner = TestRunner(validator_path)
    success = runner.run_all_tests()
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
