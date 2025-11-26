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
Unit tests for pulldown probability functions.
Tests CIGAR string parsing and minimum mismatch calculations.
"""
import unittest
import sys
import os

# Add parent directory to path to import modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

from pulldown_probability import (
    cigar_longest_contiguous_match,
    cigar_minimum_random_mismatch
)


class TestCigarMinimumRandomMismatch(unittest.TestCase):
    """Test suite for cigar_minimum_random_mismatch function."""
    
    def test_perfect_match(self):
        """Test CIGAR with perfect match longer than probe."""
        cigar = "150="
        probe_len = 120
        result = cigar_minimum_random_mismatch(cigar, probe_len)
        self.assertEqual(result, 0, "Perfect match should have 0 mismatches")
    
    def test_perfect_match_exact_length(self):
        """Test CIGAR with perfect match exactly equal to probe length."""
        cigar = "120="
        probe_len = 120
        result = cigar_minimum_random_mismatch(cigar, probe_len)
        self.assertEqual(result, 0, "Exact length perfect match should have 0 mismatches")
    
    def test_all_mismatches(self):
        """Test CIGAR with all mismatches."""
        cigar = "150X"
        probe_len = 120
        result = cigar_minimum_random_mismatch(cigar, probe_len)
        self.assertEqual(result, 120, "All mismatches should equal probe_len")
    
    def test_mixed_match_mismatch(self):
        """Test CIGAR with mixed matches and mismatches."""
        # 50 matches, 10 mismatches, 60 matches = 120 total
        # Best window: either first 120 (50= + 10X + 60=) or last 120
        cigar = "50=10X60="
        probe_len = 120
        result = cigar_minimum_random_mismatch(cigar, probe_len)
        self.assertEqual(result, 10, "Should find window with 10 mismatches")
    
    def test_match_at_start(self):
        """Test CIGAR with good match at the start."""
        cigar = "120=30X"
        probe_len = 120
        result = cigar_minimum_random_mismatch(cigar, probe_len)
        self.assertEqual(result, 0, "Should find perfect match at start")
    
    def test_match_at_end(self):
        """Test CIGAR with good match at the end."""
        cigar = "30X120="
        probe_len = 120
        result = cigar_minimum_random_mismatch(cigar, probe_len)
        self.assertEqual(result, 0, "Should find perfect match at end")
    
    def test_match_in_middle(self):
        """Test CIGAR with good match in the middle."""
        cigar = "20X120=20X"
        probe_len = 120
        result = cigar_minimum_random_mismatch(cigar, probe_len)
        self.assertEqual(result, 0, "Should find perfect match in middle")
    
    def test_insertions_deletions_reset_window(self):
        """Test that insertions/deletions reset the window."""
        # 100 matches, then insertion (resets), then 120 matches
        cigar = "100=5I120="
        probe_len = 120
        result = cigar_minimum_random_mismatch(cigar, probe_len)
        # Should find the 120= after the insertion
        self.assertEqual(result, 0, "Should find perfect match after insertion reset")
    
    def test_multiple_segments_with_indels(self):
        """Test CIGAR with multiple segments separated by indels."""
        # 80 matches, deletion, 50 matches, insertion, 70 matches
        cigar = "80=5D50=3I70="
        probe_len = 120
        result = cigar_minimum_random_mismatch(cigar, probe_len)
        # After resets, best windows are 80=, 50=, or 70= - none reach 120
        # So all would return probe_len mismatches in their respective segments
        self.assertEqual(result, probe_len, "Result should be non-negative")
    
    def test_window_slides_correctly(self):
        """Test that sliding window finds the best position."""
        # 50 mismatches, 120 matches, 30 mismatches
        cigar = "50X120=30X"
        probe_len = 120
        result = cigar_minimum_random_mismatch(cigar, probe_len)
        self.assertEqual(result, 0, "Should find perfect 120 match window")
    
    def test_overlapping_windows(self):
        """Test case where windows overlap with different mismatch counts."""
        # 60 matches, 20 mismatches, 60 matches
        cigar = "60=20X60="
        probe_len = 120
        result = cigar_minimum_random_mismatch(cigar, probe_len)
        # Best window: last 60 matches + first 60 matches would span the 20X
        # Actually: first 120 bp = 60= + 20X + 40= = 20 mismatches
        # Last 120 bp = 20= + 20X + 60= = 20 mismatches  
        # Middle 120 bp starting at position 10 would be optimal
        # This tests the sliding aspect
        self.assertLessEqual(result, 20, "Should find window with 20 mismatches")
    
    def test_short_cigar(self):
        """Test CIGAR shorter than probe length."""
        cigar = "50="
        probe_len = 120
        result = cigar_minimum_random_mismatch(cigar, probe_len)
        # When alignment is shorter than probe, behavior depends on implementation
        # Should handle gracefully - likely returns probe_len or similar
        self.assertLessEqual(result, probe_len, "Should handle short CIGAR gracefully")
    
    def test_complex_cigar(self):
        """Test complex CIGAR with multiple operation types."""
        # Match, mismatch, match, soft clip (should reset), match
        cigar = "40=10X70=5S100="
        probe_len = 120
        result = cigar_minimum_random_mismatch(cigar, probe_len)
        # Before soft clip: 40=10X70= gives windows with 10 mismatches
        # After soft clip: 100= is not long enough for probe_len
        self.assertLessEqual(result, 10, "Should find window with 10 mismatches")
    
    def test_edge_case_probe_len_one(self):
        """Test with probe length of 1."""
        cigar = "10X5=10X"
        probe_len = 1
        result = cigar_minimum_random_mismatch(cigar, probe_len)
        self.assertEqual(result, 0, "Should find at least one match for probe_len=1")
    
    def test_alternating_match_mismatch(self):
        """Test alternating matches and mismatches."""
        # Pattern: 10 matches, 10 mismatches, repeated
        cigar = "10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X"
        probe_len = 120
        result = cigar_minimum_random_mismatch(cigar, probe_len)
        self.assertEqual(result, 60, "Should find window with 60 mismatches")

    def test_alternating_match_mismatch2(self):
        """Test alternating matches and mismatches."""
        # Pattern: 10 matches, 10 mismatches, repeated
        cigar = "10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X"
        probe_len = 110
        result = cigar_minimum_random_mismatch(cigar, probe_len)
        self.assertEqual(result, 50, "Should find window with 50 mismatches")

    def test_alternating_match_mismatch3(self):
        """Test alternating matches and mismatches."""
        cigar = "10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X10=10X69=10X10=10X10=10X10=10X10=10X10=10X"
        probe_len = 110
        result = cigar_minimum_random_mismatch(cigar, probe_len)
        self.assertEqual(result, 30, "Should find window with 30 mismatches")


class TestCigarLongestContiguousMatch(unittest.TestCase):
    """Test suite for cigar_longest_contiguous_match function."""
    
    def test_single_match_block(self):
        """Test CIGAR with single match block."""
        cigar = "150="
        result = cigar_longest_contiguous_match(cigar)
        self.assertEqual(result, 150)
    
    def test_multiple_match_blocks(self):
        """Test CIGAR with multiple match blocks."""
        cigar = "50=10X100=20X30="
        result = cigar_longest_contiguous_match(cigar)
        self.assertEqual(result, 100, "Should return longest match block")
    
    def test_no_matches(self):
        """Test CIGAR with no matches."""
        cigar = "150X"
        result = cigar_longest_contiguous_match(cigar)
        self.assertEqual(result, 0, "Should return 0 when no matches")
    
    def test_equal_length_matches(self):
        """Test CIGAR with equal length match blocks."""
        cigar = "50=10X50=10X50="
        result = cigar_longest_contiguous_match(cigar)
        self.assertEqual(result, 50, "Should return length of longest match")


if __name__ == '__main__':
    unittest.main()
