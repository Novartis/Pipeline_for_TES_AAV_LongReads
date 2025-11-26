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
Unit tests for BED probability file validation.
Tests format validation, overlap detection, and error handling.
"""
import unittest
import tempfile
import os
import sys

# Add parent directory to path to import modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

from bed_validator import (
    validate_bed_probability_file,
    validate_integration_probability_file,
    validate_aav_breakpoint_probability_file,
    BEDValidationError
)


class TestBEDValidator(unittest.TestCase):
    """Test suite for BED probability file validator."""
    
    def setUp(self):
        """Set up temporary file for testing."""
        self.temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.bed')
        self.temp_path = self.temp_file.name
        self.temp_file.close()
    
    def tearDown(self):
        """Clean up temporary file."""
        if os.path.exists(self.temp_path):
            os.remove(self.temp_path)
    
    def test_valid_bed_file(self):
        """Test that a valid BED file passes validation."""
        with open(self.temp_path, 'w') as f:
            f.write("# Comment line\n")
            f.write("chr1\t100\t200\t2.0\n")
            f.write("chr1\t300\t400\t1.5\n")
            f.write("chr2\t50\t150\t3.0\n")
        
        result = validate_bed_probability_file(self.temp_path)
        
        self.assertIn("chr1", result)
        self.assertIn("chr2", result)
        self.assertEqual(len(result["chr1"]), 2)
        self.assertEqual(len(result["chr2"]), 1)
        self.assertEqual(result["chr1"][0], (100, 200, 2.0))
        self.assertEqual(result["chr1"][1], (300, 400, 1.5))
        self.assertEqual(result["chr2"][0], (50, 150, 3.0))
    
    def test_empty_file_path(self):
        """Test that empty file path is allowed."""
        result = validate_bed_probability_file("")
        self.assertEqual(result, {})
        
        result = validate_bed_probability_file(None)
        self.assertEqual(result, {})
    
    def test_file_not_found(self):
        """Test that missing file raises BEDValidationError."""
        with self.assertRaises(BEDValidationError) as ctx:
            validate_bed_probability_file("/nonexistent/file.bed")
        self.assertIn("not found", str(ctx.exception))
    
    def test_invalid_column_count(self):
        """Test that files with wrong number of columns raise error."""
        with open(self.temp_path, 'w') as f:
            f.write("chr1\t100\t200\n")  # Missing probability column
        
        with self.assertRaises(BEDValidationError) as ctx:
            validate_bed_probability_file(self.temp_path)
        self.assertIn("Expected 4 tab-separated columns", str(ctx.exception))
    
    def test_invalid_start_position(self):
        """Test that non-integer start position raises error."""
        with open(self.temp_path, 'w') as f:
            f.write("chr1\tabc\t200\t1.0\n")
        
        with self.assertRaises(BEDValidationError) as ctx:
            validate_bed_probability_file(self.temp_path)
        self.assertIn("Start position must be an integer", str(ctx.exception))
    
    def test_invalid_end_position(self):
        """Test that non-integer end position raises error."""
        with open(self.temp_path, 'w') as f:
            f.write("chr1\t100\txyz\t1.0\n")
        
        with self.assertRaises(BEDValidationError) as ctx:
            validate_bed_probability_file(self.temp_path)
        self.assertIn("End position must be an integer", str(ctx.exception))
    
    def test_invalid_probability(self):
        """Test that non-numeric probability raises error."""
        with open(self.temp_path, 'w') as f:
            f.write("chr1\t100\t200\tinvalid\n")
        
        with self.assertRaises(BEDValidationError) as ctx:
            validate_bed_probability_file(self.temp_path)
        self.assertIn("Probability must be a number", str(ctx.exception))
    
    def test_start_greater_than_end(self):
        """Test that start >= end raises error."""
        with open(self.temp_path, 'w') as f:
            f.write("chr1\t200\t100\t1.0\n")
        
        with self.assertRaises(BEDValidationError) as ctx:
            validate_bed_probability_file(self.temp_path)
        self.assertIn("Start position", str(ctx.exception))
        self.assertIn("must be less than end position", str(ctx.exception))
    
    def test_start_equal_to_end(self):
        """Test that start == end raises error."""
        with open(self.temp_path, 'w') as f:
            f.write("chr1\t100\t100\t1.0\n")
        
        with self.assertRaises(BEDValidationError) as ctx:
            validate_bed_probability_file(self.temp_path)
        self.assertIn("Start position", str(ctx.exception))
        self.assertIn("must be less than end position", str(ctx.exception))
    
    def test_negative_probability(self):
        """Test that negative probability raises error."""
        with open(self.temp_path, 'w') as f:
            f.write("chr1\t100\t200\t-1.0\n")
        
        with self.assertRaises(BEDValidationError) as ctx:
            validate_bed_probability_file(self.temp_path)
        self.assertIn("Probability must be positive", str(ctx.exception))
    
    def test_zero_probability(self):
        """Test that zero probability raises error."""
        with open(self.temp_path, 'w') as f:
            f.write("chr1\t100\t200\t0.0\n")
        
        with self.assertRaises(BEDValidationError) as ctx:
            validate_bed_probability_file(self.temp_path)
        self.assertIn("Probability must be positive", str(ctx.exception))
    
    def test_overlapping_regions(self):
        """Test that overlapping regions are detected."""
        with open(self.temp_path, 'w') as f:
            f.write("chr1\t100\t300\t1.0\n")
            f.write("chr1\t250\t400\t2.0\n")  # Overlaps with previous region
        
        with self.assertRaises(BEDValidationError) as ctx:
            validate_bed_probability_file(self.temp_path)
        self.assertIn("Overlapping regions detected", str(ctx.exception))
        self.assertIn("100-300", str(ctx.exception))
        self.assertIn("250-400", str(ctx.exception))
    
    def test_adjacent_regions_allowed(self):
        """Test that adjacent (non-overlapping) regions are allowed."""
        with open(self.temp_path, 'w') as f:
            f.write("chr1\t100\t200\t1.0\n")
            f.write("chr1\t200\t300\t2.0\n")  # Starts where previous ends
        
        # Should not raise an error
        result = validate_bed_probability_file(self.temp_path)
        self.assertEqual(len(result["chr1"]), 2)
    
    def test_overlapping_regions_different_sequences_allowed(self):
        """Test that overlapping coordinates on different sequences are allowed."""
        with open(self.temp_path, 'w') as f:
            f.write("chr1\t100\t300\t1.0\n")
            f.write("chr2\t100\t300\t2.0\n")  # Same coordinates, different sequence
        
        # Should not raise an error
        result = validate_bed_probability_file(self.temp_path)
        self.assertEqual(len(result["chr1"]), 1)
        self.assertEqual(len(result["chr2"]), 1)
    
    def test_multiple_overlapping_regions(self):
        """Test detection of overlaps with multiple regions."""
        with open(self.temp_path, 'w') as f:
            f.write("chr1\t100\t200\t1.0\n")
            f.write("chr1\t300\t400\t1.0\n")
            f.write("chr1\t350\t450\t1.0\n")  # Overlaps with second region
        
        with self.assertRaises(BEDValidationError) as ctx:
            validate_bed_probability_file(self.temp_path)
        self.assertIn("Overlapping regions detected", str(ctx.exception))
        self.assertIn("300-400", str(ctx.exception))
        self.assertIn("350-450", str(ctx.exception))
    
    def test_unsorted_regions_detected(self):
        """Test that overlaps are detected even if regions are unsorted in file."""
        with open(self.temp_path, 'w') as f:
            f.write("chr1\t300\t400\t1.0\n")  # Later region first
            f.write("chr1\t100\t350\t2.0\n")  # Earlier region that overlaps
        
        with self.assertRaises(BEDValidationError) as ctx:
            validate_bed_probability_file(self.temp_path)
        self.assertIn("Overlapping regions detected", str(ctx.exception))
    
    def test_comments_and_empty_lines_ignored(self):
        """Test that comments and empty lines are properly ignored."""
        with open(self.temp_path, 'w') as f:
            f.write("# This is a comment\n")
            f.write("\n")
            f.write("chr1\t100\t200\t1.0\n")
            f.write("\n")
            f.write("# Another comment\n")
            f.write("chr1\t300\t400\t2.0\n")
        
        result = validate_bed_probability_file(self.temp_path)
        self.assertEqual(len(result["chr1"]), 2)
    
    def test_integration_probability_file_wrapper(self):
        """Test integration site probability file validator wrapper."""
        with open(self.temp_path, 'w') as f:
            f.write("chr1\t100\t200\t1.0\n")
        
        result = validate_integration_probability_file(self.temp_path)
        self.assertIn("chr1", result)
        
        # Test error message includes correct description
        with open(self.temp_path, 'w') as f:
            f.write("chr1\tabc\t200\t1.0\n")
        
        with self.assertRaises(BEDValidationError) as ctx:
            validate_integration_probability_file(self.temp_path)
        self.assertIn("Integration site probability file", str(ctx.exception))
    
    def test_aav_breakpoint_probability_file_wrapper(self):
        """Test AAV breakpoint probability file validator wrapper."""
        with open(self.temp_path, 'w') as f:
            f.write("AAV\t100\t200\t3.0\n")
        
        result = validate_aav_breakpoint_probability_file(self.temp_path)
        self.assertIn("AAV", result)
        
        # Test error message includes correct description
        with open(self.temp_path, 'w') as f:
            f.write("AAV\t100\t50\t1.0\n")
        
        with self.assertRaises(BEDValidationError) as ctx:
            validate_aav_breakpoint_probability_file(self.temp_path)
        self.assertIn("AAV breakpoint probability file", str(ctx.exception))
    
    def test_nested_overlaps(self):
        """Test detection of nested/contained regions."""
        with open(self.temp_path, 'w') as f:
            f.write("chr1\t100\t500\t1.0\n")  # Large region
            f.write("chr1\t200\t300\t2.0\n")  # Nested inside first region
        
        with self.assertRaises(BEDValidationError) as ctx:
            validate_bed_probability_file(self.temp_path)
        self.assertIn("Overlapping regions detected", str(ctx.exception))


if __name__ == '__main__':
    unittest.main()
