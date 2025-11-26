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
Unit tests for DeletionGenerator and InsertionGenerator classes.
Tests boundary conditions, distribution formats, and error handling.
"""
import unittest
import random
import sys
import os

# Add parent directory to path to import modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

from indel_generator import DeletionGenerator, InsertionGenerator


class TestDeletionGenerator(unittest.TestCase):
    def test_uniform_within_bounds(self):
        dg = DeletionGenerator(max_size=10, size_dist_str='uniform')
        # draw many times to ensure values are within [0, max_size]
        for _ in range(1000):
            val = dg.draw()
            self.assertIsInstance(val, int)
            self.assertGreaterEqual(val, 0)
            self.assertLessEqual(val, 10)

    def test_triuniform_structure_and_bounds(self):
        # triuniform requiring 7 parts: name,max1,max2,max3,prob1,prob2,prob3
        s = 'triuniform,3,6,9,0.2,0.3,0.5'
        dg = DeletionGenerator(max_size=9, size_dist_str=s)
        # draw many times to ensure values fall within 0..max3 and obey int type
        vals = [dg.draw() for _ in range(1000)]
        self.assertTrue(all(isinstance(v, int) for v in vals))
        self.assertTrue(all(0 <= v <= 9 for v in vals))
        # ensure at least some values hit within each range
        self.assertTrue(any(0 <= v <= 3 for v in vals))
        self.assertTrue(any(3 <= v <= 6 for v in vals))
        self.assertTrue(any(6 <= v <= 9 for v in vals))

    def test_empty_string_returns_zero(self):
        dg = DeletionGenerator(max_size=5, size_dist_str='')
        for _ in range(10):
            self.assertEqual(dg.draw(), 0)

    def test_invalid_format_raises(self):
        with self.assertRaises(ValueError):
            DeletionGenerator(max_size=5, size_dist_str='notasupported,1,2')


class TestInsertionGenerator(unittest.TestCase):
    def test_dualside_diuniform_invalid(self):
        # wrong number of parts should raise ValueError
        with self.assertRaises(ValueError):
            InsertionGenerator(size_dist_str='dualside_diuniform,1,2')

    def test_draw_dualside_returns_tuple_of_strings(self):
        # For the current API, InsertionGenerator.draw() returns a tuple of two strings
        # The implementation expects 7 comma-separated parts; put placeholders for unused fields
        ig = InsertionGenerator(size_dist_str='dualside_diuniform,0,5,0.5,0.5')
        random.seed(123)
        seq1, seq2 = ig.draw()
        self.assertIsInstance(seq1, str)
        self.assertIsInstance(seq2, str)
        # lengths should be non-negative
        self.assertGreaterEqual(len(seq1), 0)
        self.assertGreaterEqual(len(seq2), 0)

    def test_draw_lengths_respect_distribution(self):
        ig = InsertionGenerator(size_dist_str='dualside_diuniform,0,6,0.3,0.7')
        random.seed(42)
        vals = [len(s[0]) + len(s[1]) for s in (ig.draw() for _ in range(200))]
        # ensure some small and some larger combined lengths exist
        self.assertTrue(any(v <= 4 for v in vals))
        self.assertTrue(any(v > 4 for v in vals))

    def test_uniform_power_invalid_parameters(self):
        # Test wrong number of parameters
        with self.assertRaises(ValueError) as context:
            InsertionGenerator(size_dist_str='uniform_power,0,10000,28.4,-0.7488,0.6')
        self.assertIn('7 parameters', str(context.exception))

    def test_uniform_power_max1_greater_than_max2(self):
        # Test max1 >= max2 raises error
        with self.assertRaises(ValueError) as context:
            InsertionGenerator(size_dist_str='uniform_power,100,50,28.4,-0.7488,0.6,0.4')
        self.assertIn('max1 < max2', str(context.exception))

    def test_uniform_power_probabilities_not_sum_to_one(self):
        # Test prob1 + prob2 != 1.0 raises error
        with self.assertRaises(ValueError) as context:
            InsertionGenerator(size_dist_str='uniform_power,0,10000,28.4,-0.7488,0.5,0.6')
        self.assertIn('prob1 + prob2 = 1.0', str(context.exception))

    def test_uniform_power_returns_tuple_of_strings(self):
        # Test that draw() returns tuple of two DNA sequences
        ig = InsertionGenerator(size_dist_str='uniform_power,0,10000,28.4,-0.7488,0.6,0.4')
        random.seed(123)
        seq1, seq2 = ig.draw()
        self.assertIsInstance(seq1, str)
        self.assertIsInstance(seq2, str)
        # Ensure sequences only contain valid DNA bases
        valid_bases = set('ATCG')
        self.assertTrue(all(c in valid_bases for c in seq1))
        self.assertTrue(all(c in valid_bases for c in seq2))

    def test_uniform_power_lengths_within_bounds(self):
        # Test that generated lengths are within expected bounds
        ig = InsertionGenerator(size_dist_str='uniform_power,0,1000,28.4,-0.7488,0.6,0.4')
        random.seed(42)
        lengths = [ig.draw_length() for _ in range(500)]
        # All lengths should be within [0, 1000]
        self.assertTrue(all(0 <= l <= 1000 for l in lengths))
        # Should have some zero-length insertions (from uniform component)
        self.assertTrue(any(l == 0 for l in lengths))
        # Should have some non-zero insertions (from power-law component)
        self.assertTrue(any(l > 0 for l in lengths))

    def test_uniform_power_probability_mixture(self):
        # Test that the probability mixture is approximately correct
        ig = InsertionGenerator(size_dist_str='uniform_power,0,100,1.0,-0.5,0.7,0.3')
        random.seed(999)
        # With max1=0, uniform component always gives 0
        # So we can test the mixture ratio by counting zeros vs non-zeros
        lengths = [ig.draw_length() for _ in range(2000)]
        zero_count = sum(1 for l in lengths if l == 0)
        nonzero_count = sum(1 for l in lengths if l > 0)
        # Should be approximately 70% zeros, 30% non-zeros (allow 10% tolerance)
        zero_ratio = zero_count / len(lengths)
        self.assertGreater(zero_ratio, 0.6)  # At least 60%
        self.assertLess(zero_ratio, 0.8)     # At most 80%

    def test_uniform_power_power_law_weighting(self):
        # Test that power-law weighting favors smaller values with negative exponent
        ig = InsertionGenerator(size_dist_str='uniform_power,10,1000,1.0,-1.0,0.0,1.0')
        random.seed(555)
        # With prob1=0.0, all samples come from power-law component
        lengths = [ig.draw_length() for _ in range(1000)]
        # Filter out any accidental zeros (shouldn't happen with prob1=0)
        nonzero_lengths = [l for l in lengths if l > 10]
        # With negative exponent, median should be closer to min than max
        import numpy as np
        median = np.median(nonzero_lengths)
        # Median should be much closer to 11 than to 1000
        self.assertLess(median, 300)  # Should be heavily weighted toward small values

    def test_uniform_power_precomputed_arrays(self):
        # Test that x and x_prob arrays are correctly precomputed
        ig = InsertionGenerator(size_dist_str='uniform_power,5,15,2.0,-0.5,0.5,0.5')
        # Check that x array has correct range
        import numpy as np
        self.assertEqual(len(ig.x), 10)  # [6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        self.assertEqual(ig.x[0], 6)
        self.assertEqual(ig.x[-1], 15)
        # Check that probabilities sum to 1
        self.assertAlmostEqual(ig.x_prob.sum(), 1.0, places=10)
        # Check that probabilities are all non-negative
        self.assertTrue(all(p >= 0 for p in ig.x_prob))


if __name__ == '__main__':
    unittest.main()
