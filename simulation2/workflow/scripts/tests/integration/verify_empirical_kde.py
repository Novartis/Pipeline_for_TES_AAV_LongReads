#!/usr/bin/env python3
"""
Copyright 2025 Novartis Institutes for BioMedical Research Inc.
 
Licensed under the MIT License (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
 
https://www.mit.edu/~amini/LICENSE.md
 
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

"""
Verification script for empirical_kde method implementation.
This script tests the complete workflow of the new method.
"""

import sys
import os
import tempfile

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from aav_fragment_generator import AAVFragmentGenerator
import numpy as np


def test_empirical_kde_workflow():
    """Test the complete empirical_kde workflow."""
    
    print("=" * 70)
    print("EMPIRICAL_KDE METHOD - COMPLETE WORKFLOW VERIFICATION")
    print("=" * 70)
    
    # Step 1: Create empirical data file
    print("\n[1] Creating empirical size data file...")
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
        empirical_file = f.name
        f.write("# Example empirical AAV fragment sizes\n")
        # Create a bimodal distribution
        sizes = list(range(400, 800, 50)) * 3 + list(range(1000, 1400, 50)) * 5
        for size in sizes:
            f.write(f"{size}\n")
    print(f"   Created: {empirical_file}")
    print(f"   Data points: {len(sizes)}")
    
    # Step 2: Create test AAV sequence
    print("\n[2] Creating test AAV sequence...")
    test_aav_seq = "ACGT" * 1000  # 4000 bp sequence
    print(f"   AAV length: {len(test_aav_seq)} bp")
    
    # Step 3: Initialize AAVFragmentGenerator with empirical_kde method
    print("\n[3] Initializing AAVFragmentGenerator with empirical_kde...")
    method_str = "empirical_kde,300,1500,50,10000"
    print(f"   Method: {method_str}")
    
    try:
        gen = AAVFragmentGenerator(
            aav_seq=test_aav_seq,
            method_str=method_str,
            empirical_sizes=empirical_file,
            seed=12345
        )
        print("   ✓ Generator initialized successfully")
    except Exception as e:
        print(f"   ✗ Failed to initialize: {e}")
        os.unlink(empirical_file)
        return False
    
    # Step 4: Generate fragments
    print("\n[4] Generating AAV fragments...")
    num_fragments = 50
    fragments_data = []
    
    try:
        for i in range(num_fragments):
            frag_seq, start, end, length = gen.generate()
            fragments_data.append({
                'seq': frag_seq,
                'start': start,
                'end': end,
                'length': length
            })
        print(f"   ✓ Generated {num_fragments} fragments successfully")
    except Exception as e:
        print(f"   ✗ Failed to generate fragments: {e}")
        os.unlink(empirical_file)
        return False
    
    # Step 5: Validate results
    print("\n[5] Validating results...")
    lengths = [f['length'] for f in fragments_data]
    
    # Check all fragments are within bounds
    min_constraint = 300
    max_constraint = 1500
    within_bounds = all(min_constraint <= l <= max_constraint for l in lengths)
    print(f"   All fragments within [{min_constraint}, {max_constraint}] bp: {within_bounds}")
    
    # Check all fragments are valid substrings
    all_valid = all(f['seq'] in test_aav_seq for f in fragments_data)
    print(f"   All fragments are valid AAV substrings: {all_valid}")
    
    # Check position consistency
    consistent = all(
        f['end'] - f['start'] + 1 == f['length']
        for f in fragments_data
    )
    print(f"   Position calculations consistent: {consistent}")
    
    # Step 6: Statistical analysis
    print("\n[6] Statistical analysis of generated fragments:")
    print(f"   Mean length:   {np.mean(lengths):.1f} bp")
    print(f"   Median length: {np.median(lengths):.1f} bp")
    print(f"   Std deviation: {np.std(lengths):.1f} bp")
    print(f"   Min length:    {np.min(lengths)} bp")
    print(f"   Max length:    {np.max(lengths)} bp")
    
    # Step 7: Check distribution follows empirical data
    print("\n[7] Distribution characteristics:")
    # The empirical data was bimodal (400-800 and 1000-1400)
    mode1_count = sum(1 for l in lengths if 300 <= l < 900)
    mode2_count = sum(1 for l in lengths if 900 <= l <= 1500)
    print(f"   Fragments in mode 1 [300-900): {mode1_count}")
    print(f"   Fragments in mode 2 [900-1500]: {mode2_count}")
    print(f"   Expected: More fragments in mode 2 (higher density in empirical data)")
    
    # Cleanup
    os.unlink(empirical_file)
    
    # Final verdict
    print("\n" + "=" * 70)
    if within_bounds and all_valid and consistent:
        print("✓ VERIFICATION PASSED - empirical_kde method working correctly!")
        print("=" * 70)
        return True
    else:
        print("✗ VERIFICATION FAILED - Some tests did not pass")
        print("=" * 70)
        return False


def test_error_handling():
    """Test error handling for empirical_kde method."""
    
    print("\n" + "=" * 70)
    print("ERROR HANDLING VERIFICATION")
    print("=" * 70)
    
    test_aav_seq = "ACGT" * 500
    
    # Test 1: Missing empirical_sizes parameter
    print("\n[1] Testing missing empirical_sizes parameter...")
    try:
        gen = AAVFragmentGenerator(test_aav_seq, "empirical_kde,100,1000,50,10000")
        print("   ✗ Should have raised ValueError")
        return False
    except ValueError as e:
        if "requires empirical_sizes" in str(e):
            print("   ✓ Correctly raises ValueError for missing empirical_sizes")
        else:
            print(f"   ✗ Wrong error message: {e}")
            return False
    
    # Test 2: Non-existent file
    print("\n[2] Testing non-existent empirical sizes file...")
    try:
        gen = AAVFragmentGenerator(
            test_aav_seq, 
            "empirical_kde,100,1000,50,10000",
            empirical_sizes="/nonexistent/file.txt"
        )
        print("   ✗ Should have raised FileNotFoundError")
        return False
    except FileNotFoundError:
        print("   ✓ Correctly raises FileNotFoundError for missing file")
    
    # Test 3: Wrong parameter count
    print("\n[3] Testing wrong parameter count...")
    try:
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write("500\n1000\n")
            temp_file = f.name
        
        gen = AAVFragmentGenerator(
            test_aav_seq,
            "empirical_kde,100,1000",  # Only 3 params instead of 5
            empirical_sizes=temp_file
        )
        print("   ✗ Should have raised ValueError")
        os.unlink(temp_file)
        return False
    except ValueError as e:
        os.unlink(temp_file)
        if "requires 5 parameters" in str(e):
            print("   ✓ Correctly raises ValueError for wrong parameter count")
        else:
            print(f"   ✗ Wrong error message: {e}")
            return False
    
    print("\n" + "=" * 70)
    print("✓ ERROR HANDLING VERIFICATION PASSED")
    print("=" * 70)
    return True


if __name__ == "__main__":
    print("\n" + "#" * 70)
    print("# EMPIRICAL_KDE METHOD - COMPREHENSIVE VERIFICATION")
    print("#" * 70)
    
    # Run workflow test
    workflow_passed = test_empirical_kde_workflow()
    
    # Run error handling test
    error_handling_passed = test_error_handling()
    
    # Final summary
    print("\n" + "#" * 70)
    print("# FINAL SUMMARY")
    print("#" * 70)
    print(f"\nWorkflow Test: {'PASSED ✓' if workflow_passed else 'FAILED ✗'}")
    print(f"Error Handling Test: {'PASSED ✓' if error_handling_passed else 'FAILED ✗'}")
    
    if workflow_passed and error_handling_passed:
        print("\n🎉 ALL VERIFICATIONS PASSED!")
        print("\nThe empirical_kde method is fully functional and ready to use.")
        sys.exit(0)
    else:
        print("\n❌ SOME VERIFICATIONS FAILED")
        print("\nPlease review the output above for details.")
        sys.exit(1)
