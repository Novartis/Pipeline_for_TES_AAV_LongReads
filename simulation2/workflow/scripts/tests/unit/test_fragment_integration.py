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
Unit tests for fragment_integration function.
Tests fragmentation of sequences containing AAV insertions.
"""

import sys
import os

# Add parent directory to path to import modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

# Try to import the actual function
try:
    from add_integration import fragment_integration
    print("✓ Successfully imported fragment_integration from add_integration.py\n")
except ImportError as e:
    print(f"✗ Could not import from add_integration.py: {e}")
    print("  Make sure numpy and biopython are installed:")
    print("  pip install numpy biopython")
    sys.exit(1)


class MockFragmentLengths:
    """Mock object to simulate FragmentLengths generator with predetermined lengths."""
    
    def __init__(self, lengths):
        self.lengths = lengths
        self.index = 0
    
    def get_fragment_length(self):
        if self.index < len(self.lengths):
            length = self.lengths[self.index]
            self.index += 1
            return length
        return self.lengths[-1] if self.lengths else 100


def test_case1_spans_entire_insertion():
    """Test Case 1: Single fragment spans the entire insertion"""
    print("Test 1: Fragment spans entire insertion...")
    
    left_flank = "A" * 35
    aav = "T" * 30
    right_flank = "G" * 35
    intgr_seq = left_flank + aav + right_flank
    intgr_start = 35
    intgr_end = 65
    
    frag_lengths = MockFragmentLengths([100])
    result = fragment_integration(intgr_seq, intgr_start, intgr_end, frag_lengths)
    frag_seq_lst, frag_side_lst, flank_left, read_split, flank_right = result
    
    assert len(frag_seq_lst) == 1, f"Expected 1 fragment, got {len(frag_seq_lst)}"
    assert frag_side_lst == ['M'], f"Expected ['M'], got {frag_side_lst}"
    assert flank_left == 35, f"Expected flank_left=35, got {flank_left}"
    assert read_split == 0, f"Expected read_split=0, got {read_split}"
    assert flank_right == 35, f"Expected flank_right=35, got {flank_right}"
    print("  ✓ PASSED")


def test_case2_entirely_within_insertion():
    """Test Case 2: Fragments entirely within insertion are discarded"""
    print("Test 2: Fragments entirely within insertion...")
    
    # Setup: 100bp sequence with 60bp AAV inserted at position 20-80
    # Fragments: [0-30] (L: 20bp flank + 10bp AAV), 
    #            [30-60] (middle, discarded), 
    #            [60-90] (R: 20bp AAV + 10bp flank), 
    #            [90-100] (after insertion)
    left_flank = "A" * 20
    aav = "T" * 60
    right_flank = "G" * 20
    intgr_seq = left_flank + aav + right_flank
    intgr_start = 20
    intgr_end = 80
    
    frag_lengths = MockFragmentLengths([30, 30, 30, 30])
    result = fragment_integration(intgr_seq, intgr_start, intgr_end, frag_lengths)
    frag_seq_lst, frag_side_lst, flank_left, read_split, flank_right = result
    
    assert len(frag_seq_lst) == 2, f"Expected 2 fragments, got {len(frag_seq_lst)}"
    assert frag_side_lst == ['L', 'R'], f"Expected ['L', 'R'], got {frag_side_lst}"
    assert flank_left == 20, f"Expected flank_left=20, got {flank_left}"
    assert isinstance(read_split, tuple), f"Expected tuple for read_split, got {type(read_split)}"
    assert read_split[0] == 10, f"Expected read_split[0]=10, got {read_split[0]}"
    # Right fragment starts at position 60, which is 40bp into the AAV (60-20=40)
    assert read_split[1] == 40, f"Expected read_split[1]=40, got {read_split[1]}"
    assert flank_right == 10, f"Expected flank_right=10, got {flank_right}"
    print("  ✓ PASSED")


def test_case3_left_flank_only():
    """Test Case 3: Fragment with left flank only"""
    print("Test 3: Fragment with left flank only...")
    
    left_flank = "A" * 20
    aav = "T" * 30
    intgr_seq = left_flank + aav
    intgr_start = 20
    intgr_end = 50
    
    frag_lengths = MockFragmentLengths([50])
    result = fragment_integration(intgr_seq, intgr_start, intgr_end, frag_lengths)
    frag_seq_lst, frag_side_lst, flank_left, read_split, flank_right = result
    
    assert len(frag_seq_lst) == 1, f"Expected 1 fragment, got {len(frag_seq_lst)}"
    assert frag_side_lst == ['L'], f"Expected ['L'], got {frag_side_lst}"
    assert flank_left == 20, f"Expected flank_left=20, got {flank_left}"
    assert read_split == 30, f"Expected read_split=30, got {read_split}"
    assert flank_right == 0, f"Expected flank_right=0, got {flank_right}"
    print("  ✓ PASSED")


def test_case4_right_flank_only():
    """Test Case 4: Fragment with right flank only"""
    print("Test 4: Fragment with right flank only...")
    
    aav = "T" * 30
    right_flank = "G" * 20
    intgr_seq = aav + right_flank
    intgr_start = 0
    intgr_end = 30
    
    frag_lengths = MockFragmentLengths([50])
    result = fragment_integration(intgr_seq, intgr_start, intgr_end, frag_lengths)
    frag_seq_lst, frag_side_lst, flank_left, read_split, flank_right = result
    
    assert len(frag_seq_lst) == 1, f"Expected 1 fragment, got {len(frag_seq_lst)}"
    assert frag_side_lst == ['R'], f"Expected ['R'], got {frag_side_lst}"
    assert flank_left == 0, f"Expected flank_left=0, got {flank_left}"
    assert read_split == 0, f"Expected read_split=0, got {read_split}"
    assert flank_right == 20, f"Expected flank_right=20, got {flank_right}"
    print("  ✓ PASSED")


def test_multiple_small_fragments():
    """Test: Multiple small fragments spanning insertion"""
    print("Test 5: Multiple small fragments...")
    
    # Setup: 100bp sequence with 50bp AAV inserted at position 25-75
    # With 20bp fragments:
    # [0-20] (before insertion), [20-40] (L: 5bp flank + 15bp AAV),
    # [40-60] (middle, discarded), [60-80] (R: 15bp AAV + 5bp flank),
    # [80-100] (after insertion)
    left_flank = "A" * 25
    aav = "T" * 50
    right_flank = "G" * 25
    intgr_seq = left_flank + aav + right_flank
    intgr_start = 25
    intgr_end = 75
    
    frag_lengths = MockFragmentLengths([20] * 10)
    result = fragment_integration(intgr_seq, intgr_start, intgr_end, frag_lengths)
    frag_seq_lst, frag_side_lst, flank_left, read_split, flank_right = result
    
    assert len(frag_seq_lst) == 2, f"Expected 2 fragments, got {len(frag_seq_lst)}"
    assert frag_side_lst == ['L', 'R'], f"Expected ['L', 'R'], got {frag_side_lst}"
    # Left fragment starts at 20, insertion at 25, so flank is 5bp
    assert flank_left == 5, f"Expected flank_left=5, got {flank_left}"
    assert isinstance(read_split, tuple), f"Expected tuple for read_split, got {type(read_split)}"
    # Right fragment ends at 80, insertion ends at 75, so flank is 5bp
    assert flank_right == 5, f"Expected flank_right=5, got {flank_right}"
    print("  ✓ PASSED")


def test_very_small_insertion():
    """Test: Very small insertion spanned by one fragment"""
    print("Test 6: Very small insertion...")
    
    left_flank = "A" * 47
    aav = "T" * 5
    right_flank = "G" * 48
    intgr_seq = left_flank + aav + right_flank
    intgr_start = 47
    intgr_end = 52
    
    frag_lengths = MockFragmentLengths([100])
    result = fragment_integration(intgr_seq, intgr_start, intgr_end, frag_lengths)
    frag_seq_lst, frag_side_lst, flank_left, read_split, flank_right = result
    
    assert len(frag_seq_lst) == 1, f"Expected 1 fragment, got {len(frag_seq_lst)}"
    assert frag_side_lst == ['M'], f"Expected ['M'], got {frag_side_lst}"
    assert flank_left == 47, f"Expected flank_left=47, got {flank_left}"
    assert read_split == 0, f"Expected read_split=0, got {read_split}"
    assert flank_right == 48, f"Expected flank_right=48, got {flank_right}"
    print("  ✓ PASSED")


if __name__ == '__main__':
    print("\n" + "="*60)
    print("Running fragment_integration unit tests")
    print("="*60 + "\n")
    
    tests = [
        test_case1_spans_entire_insertion,
        test_case2_entirely_within_insertion,
        test_case3_left_flank_only,
        test_case4_right_flank_only,
        test_multiple_small_fragments,
        test_very_small_insertion,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            passed += 1
        except AssertionError as e:
            failed += 1
            print(f"  ✗ FAILED: {e}")
        except Exception as e:
            failed += 1
            print(f"  ✗ ERROR: {e}")
    
    print("\n" + "="*60)
    print(f"Results: {passed} passed, {failed} failed out of {len(tests)} tests")
    print("="*60 + "\n")
    
    if failed == 0:
        print("✓ All tests passed!")
        exit(0)
    else:
        print("✗ Some tests failed")
        exit(1)
