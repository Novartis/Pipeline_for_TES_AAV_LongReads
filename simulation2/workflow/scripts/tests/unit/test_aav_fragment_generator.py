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
Unit tests for AAVFragmentGenerator class.
Tests API, error handling, and basic functionality for both gamma_stick and gamma_biased_stick methods.
"""
import unittest
import random
import os
import sys

# Add parent directory to path to import modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

from aav_fragment_generator import AAVFragmentGenerator


class TestAAVFragmentGenerator(unittest.TestCase):
    """Test suite for AAVFragmentGenerator class."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test AAV sequence and probability BED file."""
        # Create a test AAV sequence (simplified version)
        cls.test_aav_seq = "ACGT" * 250  # 1000 bp test sequence
        
        # Create a temporary BED file for testing
        cls.test_bed_file = "/tmp/test_aav_probs.bed"
        with open(cls.test_bed_file, 'w') as f:
            f.write("# Test BED file for AAV breakpoint probabilities\n")
            f.write("test_aav\t0\t250\t0.5\n")      # Low probability region
            f.write("test_aav\t250\t750\t3.0\n")    # High probability region (hotspot)
            f.write("test_aav\t750\t1000\t1.0\n")   # Baseline probability region
        
        # Create a temporary empirical sizes file for testing
        cls.test_empirical_file = "/tmp/test_empirical_sizes.txt"
        with open(cls.test_empirical_file, 'w') as f:
            f.write("# Test empirical fragment sizes\n")
            # Create a distribution centered around 400-600bp
            for size in [300, 350, 400, 450, 500, 550, 600, 650, 700]:
                for _ in range(5):  # Repeat each size 5 times
                    f.write(f"{size}\n")
            # Add some more variation
            for size in [400, 500, 600]:
                for _ in range(10):
                    f.write(f"{size}\n")
    
    @classmethod
    def tearDownClass(cls):
        """Clean up temporary BED file."""
        if os.path.exists(cls.test_bed_file):
            os.remove(cls.test_bed_file)
        if os.path.exists(cls.test_empirical_file):
            os.remove(cls.test_empirical_file)
    
    def test_gamma_stick_basic_functionality(self):
        """Test that gamma_stick generates valid fragments."""
        method_str = "gamma_stick,500,100,50"
        gen = AAVFragmentGenerator(self.test_aav_seq, method_str)
        
        frag_seq, start, end, length = gen.generate()
        
        # Check return types
        self.assertIsInstance(frag_seq, str)
        self.assertIsInstance(start, int)
        self.assertIsInstance(end, int)
        self.assertIsInstance(length, int)
        
        # Check fragment is from AAV sequence
        self.assertIn(frag_seq, self.test_aav_seq)
        
        # Check positions are 1-based and consistent
        self.assertGreaterEqual(start, 1)
        self.assertLessEqual(end, len(self.test_aav_seq))
        self.assertEqual(length, len(frag_seq))
        self.assertEqual(end - start + 1, length)
        
        # Check minimum length is respected
        self.assertGreaterEqual(length, 50)
    
    def test_gamma_stick_respects_min_length(self):
        """Test that gamma_stick enforces minimum fragment length."""
        method_str = "gamma_stick,100,50,80"  # min_length = 80
        gen = AAVFragmentGenerator(self.test_aav_seq, method_str)
        
        # Generate multiple fragments
        for _ in range(10):
            frag_seq, start, end, length = gen.generate()
            self.assertGreaterEqual(length, 80, "Fragment length should be >= min_length")
    
    def test_gamma_biased_stick_basic_functionality(self):
        """Test that gamma_biased_stick generates valid fragments."""
        method_str = f"gamma_biased_stick,500,100,50,100,10000"
        gen = AAVFragmentGenerator(self.test_aav_seq, method_str, prob_file=self.test_bed_file)
        
        frag_seq, start, end, length = gen.generate()
        
        # Check return types
        self.assertIsInstance(frag_seq, str)
        self.assertIsInstance(start, int)
        self.assertIsInstance(end, int)
        self.assertIsInstance(length, int)
        
        # Check fragment is from AAV sequence
        self.assertIn(frag_seq, self.test_aav_seq)
        
        # Check positions are 1-based and consistent
        self.assertGreaterEqual(start, 1)
        self.assertLessEqual(end, len(self.test_aav_seq))
        self.assertEqual(length, len(frag_seq))
        
        # Check minimum length is respected
        self.assertGreaterEqual(length, 50)
    
    def test_gamma_biased_stick_respects_tolerance(self):
        """Test that gamma_biased_stick generates fragments within tolerance."""
        random.seed(42)
        target_mean = 400
        tolerance = 50
        method_str = f"gamma_biased_stick,{target_mean},100,50,{tolerance},10000"
        gen = AAVFragmentGenerator(self.test_aav_seq, method_str, prob_file=self.test_bed_file)
        
        # Generate multiple fragments and check tolerance
        fragment_lengths = []
        for _ in range(20):
            frag_seq, start, end, length = gen.generate()
            fragment_lengths.append(length)
        
        # Most fragments should be reasonably close to the target mean
        # (not all will be exactly within tolerance due to sampling, but average should be close)
        avg_length = sum(fragment_lengths) / len(fragment_lengths)
        self.assertGreater(avg_length, target_mean - 150, 
                          f"Average length {avg_length} too far below target {target_mean}")
        self.assertLess(avg_length, target_mean + 150,
                       f"Average length {avg_length} too far above target {target_mean}")
    
    def test_gamma_biased_stick_bias_distribution(self):
        """Test that biased sampling favors high-probability regions."""
        random.seed(123)
        method_str = f"gamma_biased_stick,300,50,50,100,10000"
        gen = AAVFragmentGenerator(self.test_aav_seq, method_str, prob_file=self.test_bed_file)
        
        # Sample many breakpoints and count which regions they fall in
        # BED file has: 0-250 (prob=0.5), 250-750 (prob=3.0), 750-1000 (prob=1.0)
        region_counts = {"low": 0, "high": 0, "baseline": 0}
        
        for _ in range(100):
            frag_seq, start, end, length = gen.generate()
            # Check where fragment center is located (convert to 0-based)
            center = (start + end) // 2 - 1
            
            if center < 250:
                region_counts["low"] += 1
            elif center < 750:
                region_counts["high"] += 1
            else:
                region_counts["baseline"] += 1
        
        # High probability region (3.0x) should have most fragments
        self.assertGreater(region_counts["high"], region_counts["low"],
                          "High-probability region should have more fragments than low-probability")
        self.assertGreater(region_counts["high"], region_counts["baseline"],
                          "High-probability region should have more fragments than baseline")
    
    def test_invalid_method_name_raises(self):
        """Test that invalid method name raises ValueError."""
        with self.assertRaises(ValueError) as ctx:
            AAVFragmentGenerator(self.test_aav_seq, "invalid_method,100,50,25")
        self.assertIn("Unsupported", str(ctx.exception))
    
    def test_gamma_stick_wrong_param_count_raises(self):
        """Test that wrong parameter count for gamma_stick raises ValueError."""
        with self.assertRaises(ValueError) as ctx:
            AAVFragmentGenerator(self.test_aav_seq, "gamma_stick,100,50")  # Missing min param
        self.assertIn("requires 4 parameters", str(ctx.exception))
    
    def test_gamma_biased_stick_wrong_param_count_raises(self):
        """Test that wrong parameter count for gamma_biased_stick raises ValueError."""
        with self.assertRaises(ValueError) as ctx:
            AAVFragmentGenerator(self.test_aav_seq, "gamma_biased_stick,100,50,25", prob_file=self.test_bed_file)
        self.assertIn("requires 6 parameters", str(ctx.exception))
    
    def test_gamma_biased_stick_missing_bed_file_raises(self):
        """Test that missing BED file raises FileNotFoundError."""
        with self.assertRaises(FileNotFoundError):
            AAVFragmentGenerator(self.test_aav_seq, "gamma_biased_stick,100,50,25,50,10000", prob_file="/nonexistent/file.bed")
    
    def test_gamma_biased_stick_no_prob_file_raises(self):
        """Test that gamma_biased_stick without prob_file raises ValueError."""
        with self.assertRaises(ValueError) as ctx:
            AAVFragmentGenerator(self.test_aav_seq, "gamma_biased_stick,100,50,25,50,10000")
        self.assertIn("requires prob_file parameter", str(ctx.exception))
    
    def test_empty_method_string_raises(self):
        """Test that empty method string raises ValueError."""
        with self.assertRaises(ValueError) as ctx:
            AAVFragmentGenerator(self.test_aav_seq, "")
        self.assertIn("Empty", str(ctx.exception))
    
    def test_gamma_biased_stick_max_resets_raises(self):
        """Test that exceeding 100 length resamples raises RuntimeError."""
        # Create parameters that make it nearly impossible to find valid fragments
        # Very tight tolerance + very small max_iter + mismatch with sequence length
        method_str = "gamma_biased_stick,5000,100,4000,10,5"  # target 5000bp, tolerance 10bp, max_iter=5
        gen = AAVFragmentGenerator(self.test_aav_seq, method_str, prob_file=self.test_bed_file)
        
        # Should fail after 100 resets and raise RuntimeError
        with self.assertRaises(RuntimeError) as ctx:
            gen.generate()
        
        self.assertIn("100 length resamples", str(ctx.exception))
        self.assertIn("mismatch", str(ctx.exception))
    
    def test_empirical_kde_basic_functionality(self):
        """Test that empirical_kde generates valid fragments."""
        method_str = "empirical_kde,50,800,50,10000"
        gen = AAVFragmentGenerator(self.test_aav_seq, method_str, empirical_sizes=self.test_empirical_file)
        
        frag_seq, start, end, length = gen.generate()
        
        # Check return types
        self.assertIsInstance(frag_seq, str)
        self.assertIsInstance(start, int)
        self.assertIsInstance(end, int)
        self.assertIsInstance(length, int)
        
        # Check fragment is from AAV sequence
        self.assertIn(frag_seq, self.test_aav_seq)
        
        # Check positions are 1-based and consistent
        self.assertGreaterEqual(start, 1)
        self.assertLessEqual(end, len(self.test_aav_seq))
        self.assertEqual(length, len(frag_seq))
        self.assertEqual(end - start + 1, length)
        
        # Check min/max constraints are respected
        self.assertGreaterEqual(length, 50)
        self.assertLessEqual(length, 800)
    
    def test_empirical_kde_respects_min_max(self):
        """Test that empirical_kde enforces min/max constraints."""
        min_length = 100
        max_length = 700
        method_str = f"empirical_kde,{min_length},{max_length},50,10000"
        gen = AAVFragmentGenerator(self.test_aav_seq, method_str, empirical_sizes=self.test_empirical_file)
        
        # Generate multiple fragments
        for _ in range(20):
            frag_seq, start, end, length = gen.generate()
            self.assertGreaterEqual(length, min_length, 
                                  f"Fragment length {length} should be >= min_length {min_length}")
            self.assertLessEqual(length, max_length,
                               f"Fragment length {length} should be <= max_length {max_length}")
    
    def test_empirical_kde_distribution(self):
        """Test that empirical_kde generates sizes following the empirical distribution."""
        random.seed(789)
        method_str = "empirical_kde,50,800,50,10000"
        gen = AAVFragmentGenerator(self.test_aav_seq, method_str, empirical_sizes=self.test_empirical_file, seed=789)
        
        # Generate many fragments and check distribution
        fragment_lengths = []
        for _ in range(50):
            frag_seq, start, end, length = gen.generate()
            fragment_lengths.append(length)
        
        # Average should be close to the empirical data mean (around 400-600)
        avg_length = sum(fragment_lengths) / len(fragment_lengths)
        self.assertGreater(avg_length, 300, 
                          f"Average length {avg_length} too far below expected range")
        self.assertLess(avg_length, 750,
                       f"Average length {avg_length} too far above expected range")
    
    def test_empirical_kde_wrong_param_count_raises(self):
        """Test that wrong parameter count for empirical_kde raises ValueError."""
        with self.assertRaises(ValueError) as ctx:
            AAVFragmentGenerator(self.test_aav_seq, "empirical_kde,50,800", 
                               empirical_sizes=self.test_empirical_file)
        self.assertIn("requires 5 parameters", str(ctx.exception))
    
    def test_empirical_kde_missing_file_raises(self):
        """Test that missing empirical sizes file raises FileNotFoundError."""
        with self.assertRaises(FileNotFoundError):
            AAVFragmentGenerator(self.test_aav_seq, "empirical_kde,50,800,50,10000", 
                               empirical_sizes="/nonexistent/sizes.txt")
    
    def test_empirical_kde_no_empirical_sizes_raises(self):
        """Test that empirical_kde without empirical_sizes raises ValueError."""
        with self.assertRaises(ValueError) as ctx:
            AAVFragmentGenerator(self.test_aav_seq, "empirical_kde,50,800,50,10000")
        self.assertIn("requires empirical_sizes parameter", str(ctx.exception))
    
    def test_empirical_kde_impossible_constraints_raises(self):
        """Test that impossible constraints raise RuntimeError."""
        # Create empirical data with sizes 300-700, but set min=800, max=900
        # This makes it impossible to sample valid sizes
        method_str = "empirical_kde,800,900,50,100"  # Very few iterations to fail fast
        gen = AAVFragmentGenerator(self.test_aav_seq, method_str, empirical_sizes=self.test_empirical_file)
        
        # Should fail after max_iter attempts
        with self.assertRaises(RuntimeError) as ctx:
            gen.generate()
        
        self.assertIn("Failed to sample valid fragment length", str(ctx.exception))
    
    # Tests for empirical_kde_biased method
    def test_empirical_kde_biased_basic(self):
        """Test basic functionality of empirical_kde_biased method."""
        random.seed(999)
        method_str = "empirical_kde_biased,50,800,50,10000"
        gen = AAVFragmentGenerator(self.test_aav_seq, method_str, 
                                  prob_file=self.test_bed_file,
                                  empirical_sizes=self.test_empirical_file,
                                  seed=999)
        
        frag_seq, start, end, length = gen.generate()
        
        # Verify output format
        self.assertIsInstance(frag_seq, str)
        self.assertIsInstance(start, int)
        self.assertIsInstance(end, int)
        self.assertIsInstance(length, int)
        
        # Verify fragment is within AAV sequence
        self.assertGreater(start, 0)
        self.assertLessEqual(end, len(self.test_aav_seq))
        self.assertEqual(length, end - start + 1)
        self.assertEqual(length, len(frag_seq))
    
    def test_empirical_kde_biased_respects_constraints(self):
        """Test that empirical_kde_biased respects min/max length constraints."""
        random.seed(1000)
        min_length = 100
        max_length = 700
        method_str = f"empirical_kde_biased,{min_length},{max_length},50,10000"
        gen = AAVFragmentGenerator(self.test_aav_seq, method_str,
                                  prob_file=self.test_bed_file,
                                  empirical_sizes=self.test_empirical_file,
                                  seed=1000)
        
        # Generate multiple fragments and check constraints
        for _ in range(20):
            frag_seq, start, end, length = gen.generate()
            self.assertGreaterEqual(length, min_length,
                                  f"Fragment length {length} should be >= min_length {min_length}")
            self.assertLessEqual(length, max_length,
                               f"Fragment length {length} should be <= max_length {max_length}")
    
    def test_empirical_kde_biased_breakpoint_bias(self):
        """Test that empirical_kde_biased uses biased breakpoints from prob_file."""
        random.seed(1001)
        method_str = "empirical_kde_biased,200,600,50,10000"
        gen = AAVFragmentGenerator(self.test_aav_seq, method_str,
                                  prob_file=self.test_bed_file,
                                  empirical_sizes=self.test_empirical_file,
                                  seed=1001)
        
        # Generate many fragments and check if breakpoints are biased toward high-prob region (250-750)
        breakpoints = []
        for _ in range(100):
            frag_seq, start, end, length = gen.generate()
            breakpoints.append(start - 1)  # Convert to 0-based
            breakpoints.append(end)
        
        # Count breakpoints in high-probability region (250-750) vs other regions
        high_prob_count = sum(1 for bp in breakpoints if 250 <= bp <= 750)
        
        # With 3x probability in region 250-750 (500bp) vs baseline in other regions (500bp total),
        # we expect roughly 75% of breakpoints in high-prob region
        # Use a loose threshold to account for randomness
        high_prob_ratio = high_prob_count / len(breakpoints)
        self.assertGreater(high_prob_ratio, 0.60,
                          f"Expected >60% breakpoints in high-prob region, got {high_prob_ratio:.2%}")
    
    def test_empirical_kde_biased_missing_prob_file_raises(self):
        """Test that empirical_kde_biased without prob_file raises ValueError."""
        with self.assertRaises(ValueError) as ctx:
            AAVFragmentGenerator(self.test_aav_seq, "empirical_kde_biased,50,800,50,10000",
                               empirical_sizes=self.test_empirical_file)
        self.assertIn("requires prob_file parameter", str(ctx.exception))
    
    def test_empirical_kde_biased_missing_empirical_sizes_raises(self):
        """Test that empirical_kde_biased without empirical_sizes raises ValueError."""
        with self.assertRaises(ValueError) as ctx:
            AAVFragmentGenerator(self.test_aav_seq, "empirical_kde_biased,50,800,50,10000",
                               prob_file=self.test_bed_file)
        self.assertIn("requires empirical_sizes parameter", str(ctx.exception))
    
    def test_empirical_kde_biased_wrong_param_count_raises(self):
        """Test that wrong parameter count for empirical_kde_biased raises ValueError."""
        with self.assertRaises(ValueError) as ctx:
            AAVFragmentGenerator(self.test_aav_seq, "empirical_kde_biased,50,800",
                               prob_file=self.test_bed_file,
                               empirical_sizes=self.test_empirical_file)
        self.assertIn("requires 5 parameters", str(ctx.exception))


if __name__ == '__main__':
    unittest.main()
