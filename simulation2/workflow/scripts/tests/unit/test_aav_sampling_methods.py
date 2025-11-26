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
Unit tests for AAV fragment generation sampling methods.
Tests biased sampling behavior, weight calculations, breakpoint distributions,
and edge cases for both gamma_stick and gamma_biased_stick methods.
"""
import unittest
import random
import os
import sys
import tempfile
from collections import Counter

# Add parent directory to path to import modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

from aav_fragment_generator import AAVFragmentGenerator


class TestBiasedSampling(unittest.TestCase):
    """Test suite for biased breakpoint sampling in gamma_biased_stick method."""
    
    def setUp(self):
        """Set up test sequence and create temporary BED file for each test."""
        # Create a simple test sequence
        self.test_seq = "ACGT" * 250  # 1000 bp sequence
        
        # Create a temporary file for BED data
        self.temp_bed = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.bed')
        self.temp_bed_path = self.temp_bed.name
    
    def tearDown(self):
        """Clean up temporary BED file after each test."""
        if os.path.exists(self.temp_bed_path):
            os.remove(self.temp_bed_path)
    
    def _write_bed_file(self, regions):
        """
        Helper to write regions to temporary BED file.
        
        Args:
            regions: List of tuples (start, end, probability)
        """
        with open(self.temp_bed_path, 'w') as f:
            f.write("# Test BED file\n")
            for start, end, prob in regions:
                f.write(f"test_seq\t{start}\t{end}\t{prob}\n")
    
    def test_uniform_sampling_no_bed_regions(self):
        """Test that breakpoints are uniform when no regions defined in BED file."""
        # Empty BED file (no regions)
        self._write_bed_file([])
        
        random.seed(42)
        gen = AAVFragmentGenerator(
            self.test_seq, 
            "gamma_biased_stick,200,50,50,50,10000",
            prob_file=self.temp_bed_path
        )
        
        # Sample many breakpoints and check distribution
        breakpoints = []
        for _ in range(100):
            frag_seq, start, end, length = gen.generate()
            breakpoints.append(start)
            breakpoints.append(end)
        
        # Check that breakpoints span the full sequence
        self.assertGreater(max(breakpoints), 800, "Breakpoints should reach end of sequence")
        self.assertLess(min(breakpoints), 200, "Breakpoints should start near beginning")
    
    def test_single_hotspot_enrichment(self):
        """Test that breakpoints are enriched in high-probability region."""
        # Define a single hotspot region (3x probability)
        self._write_bed_file([
            (400, 600, 3.0)  # Middle 200bp with 3x probability
        ])
        
        random.seed(123)
        gen = AAVFragmentGenerator(
            self.test_seq,
            "gamma_biased_stick,200,50,50,50,10000",
            prob_file=self.temp_bed_path
        )
        
        # Sample many fragments and count breakpoints in hotspot vs outside
        hotspot_count = 0
        outside_count = 0
        
        for _ in range(100):
            frag_seq, start, end, length = gen.generate()
            # Check both breakpoints
            for bp in [start, end]:
                if 400 <= bp <= 600:
                    hotspot_count += 1
                else:
                    outside_count += 1
        
        # Hotspot is 20% of sequence (200/1000) but has 3x probability
        # Expected proportion in hotspot: (200*3) / (800*1 + 200*3) = 600/1400 = 42.9%
        # We should see significantly more than 20% in hotspot
        hotspot_proportion = hotspot_count / (hotspot_count + outside_count)
        
        self.assertGreater(hotspot_proportion, 0.30, 
                          f"Hotspot should have >30% of breakpoints, got {hotspot_proportion:.2%}")
        self.assertLess(hotspot_proportion, 0.60,
                       f"Hotspot proportion should be <60%, got {hotspot_proportion:.2%}")
    
    def test_coldspot_depletion(self):
        """Test that breakpoints are depleted in low-probability region."""
        # Define a coldspot region (0.2x probability)
        self._write_bed_file([
            (300, 700, 0.2)  # Middle 400bp with 0.2x probability
        ])
        
        random.seed(456)
        gen = AAVFragmentGenerator(
            self.test_seq,
            "gamma_biased_stick,200,50,50,50,10000",
            prob_file=self.temp_bed_path
        )
        
        # Sample many fragments
        coldspot_count = 0
        outside_count = 0
        
        for _ in range(100):
            frag_seq, start, end, length = gen.generate()
            for bp in [start, end]:
                if 300 <= bp <= 700:
                    coldspot_count += 1
                else:
                    outside_count += 1
        
        # Coldspot is 40% of sequence but has 0.2x probability
        # Expected: (400*0.2) / (600*1 + 400*0.2) = 80/680 = 11.8%
        # Should be much less than 40%
        coldspot_proportion = coldspot_count / (coldspot_count + outside_count)
        
        self.assertLess(coldspot_proportion, 0.25,
                       f"Coldspot should have <25% of breakpoints, got {coldspot_proportion:.2%}")
    
    def test_multiple_regions_weighted_sampling(self):
        """Test weighted sampling across multiple regions with different probabilities."""
        # Define three regions with different probabilities
        self._write_bed_file([
            (0, 200, 0.5),      # Low probability (10% weight vs 20% size)
            (400, 600, 5.0),    # Very high probability (50% weight vs 20% size)
            (800, 1000, 1.0)    # Baseline probability (20% weight vs 20% size)
        ])
        
        random.seed(789)
        gen = AAVFragmentGenerator(
            self.test_seq,
            "gamma_biased_stick,150,30,50,40,10000",
            prob_file=self.temp_bed_path
        )
        
        # Count breakpoints in each region
        region_counts = {"low": 0, "high": 0, "baseline": 0, "undefined": 0}
        
        for _ in range(150):
            frag_seq, start, end, length = gen.generate()
            for bp in [start, end]:
                if bp < 200:
                    region_counts["low"] += 1
                elif 400 <= bp < 600:
                    region_counts["high"] += 1
                elif 800 <= bp <= 1000:
                    region_counts["baseline"] += 1
                else:
                    region_counts["undefined"] += 1
        
        # High probability region should have most breakpoints
        self.assertGreater(region_counts["high"], region_counts["low"],
                          "High-prob region should have more breakpoints than low-prob")
        self.assertGreater(region_counts["high"], region_counts["baseline"],
                          "High-prob region should have more breakpoints than baseline")
    
    def test_tolerance_constraint(self):
        """Test that tolerance parameter correctly constrains fragment length."""
        self._write_bed_file([
            (0, 500, 2.0),
            (500, 1000, 1.0)
        ])
        
        random.seed(111)
        target_mean = 300
        tolerance = 30
        
        gen = AAVFragmentGenerator(
            self.test_seq,
            f"gamma_biased_stick,{target_mean},50,50,{tolerance},10000",
            prob_file=self.temp_bed_path
        )
        
        # Generate fragments and check they respect tolerance
        lengths_outside_tolerance = 0
        
        for _ in range(50):
            frag_seq, start, end, length = gen.generate()
            
            # Most fragments should be near target_mean
            # (allowing some variance due to gamma distribution sampling)
            if abs(length - target_mean) > tolerance * 3:  # Allow 3x tolerance for statistical variation
                lengths_outside_tolerance += 1
        
        # Most should be reasonably close to target
        self.assertLess(lengths_outside_tolerance, 10,
                       f"Too many fragments far from target: {lengths_outside_tolerance}/50")
    
    def test_min_length_constraint_with_biased_sampling(self):
        """Test that minimum length is generally respected with biased sampling."""
        self._write_bed_file([
            (200, 800, 3.0)
        ])
        
        random.seed(222)
        min_length = 100
        
        gen = AAVFragmentGenerator(
            self.test_seq,
            f"gamma_biased_stick,200,50,{min_length},50,10000",
            prob_file=self.temp_bed_path
        )
        
        # Most fragments should be at least min_length, but tolerance allows some variation
        # The tolerance parameter means we accept fragments within +/- tolerance of sampled length
        below_min_count = 0
        
        for _ in range(30):
            frag_seq, start, end, length = gen.generate()
            if length < min_length:
                below_min_count += 1
        
        # Allow some fragments to be below minimum due to tolerance parameter
        # but most should still respect it
        self.assertLess(below_min_count, 10,
                       f"Too many fragments ({below_min_count}/30) below minimum {min_length}")
    
    def test_edge_case_tiny_hotspot(self):
        """Test behavior with very small high-probability region."""
        # Tiny 10bp hotspot with very high probability
        self._write_bed_file([
            (495, 505, 10.0)  # 10bp with 10x probability
        ])
        
        random.seed(333)
        gen = AAVFragmentGenerator(
            self.test_seq,
            "gamma_biased_stick,200,50,50,50,10000",
            prob_file=self.temp_bed_path
        )
        
        # Should still be able to generate fragments
        # Some should have breakpoints near the hotspot
        near_hotspot = 0
        
        for _ in range(50):
            frag_seq, start, end, length = gen.generate()
            for bp in [start, end]:
                if 490 <= bp <= 510:  # Within 5bp of hotspot
                    near_hotspot += 1
        
        # With 10x probability on 1% of sequence, should see some enrichment
        self.assertGreater(near_hotspot, 2, "Should see some breakpoints near tiny hotspot")
    
    def test_adjacent_regions_different_probs(self):
        """Test behavior with adjacent regions of different probabilities."""
        # Two adjacent regions with different probabilities
        self._write_bed_file([
            (0, 500, 0.5),      # First half: low probability
            (500, 1000, 3.0)    # Second half: high probability
        ])
        
        random.seed(444)
        gen = AAVFragmentGenerator(
            self.test_seq,
            "gamma_biased_stick,200,50,50,50,10000",
            prob_file=self.temp_bed_path
        )
        
        first_half_count = 0
        second_half_count = 0
        
        for _ in range(100):
            frag_seq, start, end, length = gen.generate()
            for bp in [start, end]:
                if bp < 500:
                    first_half_count += 1
                else:
                    second_half_count += 1
        
        # Second half should have significantly more breakpoints
        self.assertGreater(second_half_count, first_half_count * 1.5,
                          "Second half (high prob) should have >1.5x more breakpoints than first half (low prob)")
    
    def test_weight_calculation_accuracy(self):
        """Test that internal weight calculations are correct."""
        # Test that the generator correctly calculates total weights
        regions = [
            (0, 100, 0.5),      # Weight: 50
            (200, 300, 2.0),    # Weight: 200
            # 100-200: baseline 1.0, weight: 100
            # 300-1000: baseline 1.0, weight: 700
        ]
        self._write_bed_file(regions)
        
        gen = AAVFragmentGenerator(
            self.test_seq,
            "gamma_biased_stick,200,50,50,50,10000",
            prob_file=self.temp_bed_path
        )
        
        # Total weight should be: 50 + 100 + 200 + 700 = 1050
        expected_weight = 50 + 100 + 200 + 700
        self.assertEqual(gen.total_weight, expected_weight,
                        f"Total weight should be {expected_weight}, got {gen.total_weight}")
    
    def test_fragments_contain_valid_sequence(self):
        """Test that generated fragments are valid subsequences with biased sampling."""
        self._write_bed_file([
            (250, 750, 2.5)
        ])
        
        random.seed(555)
        gen = AAVFragmentGenerator(
            self.test_seq,
            "gamma_biased_stick,200,50,50,50,10000",
            prob_file=self.temp_bed_path
        )
        
        for _ in range(20):
            frag_seq, start, end, length = gen.generate()
            
            # Fragment should be in the original sequence
            self.assertIn(frag_seq, self.test_seq,
                         "Fragment sequence should be in original sequence")
            
            # Verify the extracted fragment matches the positions
            expected_frag = self.test_seq[start-1:end]  # Convert 1-based to 0-based
            self.assertEqual(frag_seq, expected_frag,
                           f"Fragment at {start}-{end} doesn't match expected sequence")


class TestGammaStickSampling(unittest.TestCase):
    """Test suite specifically for gamma_stick method sampling behavior."""
    
    def setUp(self):
        """Set up test sequence."""
        self.test_seq = "ATCG" * 500  # 2000 bp sequence
    
    def test_sequence_coverage(self):
        """Test that gamma_stick generates valid fragments from the sequence."""
        random.seed(12345)
        gen = AAVFragmentGenerator(self.test_seq, "gamma_stick,400,100,100")
        
        # Collect fragment information
        fragments = []
        
        for _ in range(50):
            frag_seq, start, end, length = gen.generate()
            fragments.append((frag_seq, start, end, length))
        
        # All fragments should be valid
        for frag_seq, start, end, length in fragments:
            self.assertIn(frag_seq, self.test_seq,
                         "Fragment should be in original sequence")
            self.assertGreater(length, 0,
                             "Fragment length should be positive")
            self.assertEqual(length, len(frag_seq),
                           "Length should match fragment sequence length")
    
    def test_stick_breaking_produces_valid_fragments(self):
        """Test that stick-breaking always produces valid fragments."""
        random.seed(67890)
        gen = AAVFragmentGenerator(self.test_seq, "gamma_stick,300,80,50")
        
        for _ in range(50):
            frag_seq, start, end, length = gen.generate()
            
            # Fragment must be valid
            self.assertIsInstance(frag_seq, str)
            self.assertGreater(len(frag_seq), 0)
            
            # Must be in original sequence
            self.assertIn(frag_seq, self.test_seq)
            
            # Length consistency
            self.assertEqual(len(frag_seq), length)
            self.assertEqual(end - start + 1, length)
    
    def test_gamma_distribution_length_sampling(self):
        """Test that fragment lengths follow gamma distribution parameters."""
        random.seed(11111)
        target_mean = 500
        target_std = 150
        min_length = 100
        
        gen = AAVFragmentGenerator(
            self.test_seq,
            f"gamma_stick,{target_mean},{target_std},{min_length}"
        )
        
        lengths = []
        for _ in range(100):
            frag_seq, start, end, length = gen.generate()
            lengths.append(length)
        
        # Check mean is roughly correct (within 20%)
        actual_mean = sum(lengths) / len(lengths)
        self.assertGreater(actual_mean, target_mean * 0.8,
                          f"Mean length {actual_mean} too far below target {target_mean}")
        self.assertLess(actual_mean, target_mean * 1.2,
                       f"Mean length {actual_mean} too far above target {target_mean}")
        
        # All should be above minimum
        self.assertTrue(all(l >= min_length for l in lengths),
                       "Some fragments below minimum length")


class TestEdgeCasesAndErrors(unittest.TestCase):
    """Test edge cases and error handling."""
    
    def setUp(self):
        """Set up test sequence and temporary file."""
        self.test_seq = "GCTA" * 100  # 400 bp
        self.temp_bed = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.bed')
        self.temp_bed_path = self.temp_bed.name
    
    def tearDown(self):
        """Clean up temporary file."""
        if os.path.exists(self.temp_bed_path):
            os.remove(self.temp_bed_path)
    
    def test_empty_bed_file_uses_uniform(self):
        """Test that empty BED file falls back to uniform sampling."""
        # Write empty BED file (only comments)
        with open(self.temp_bed_path, 'w') as f:
            f.write("# Empty BED file\n")
            f.write("# No regions defined\n")
        
        # Should work without error
        gen = AAVFragmentGenerator(
            self.test_seq,
            "gamma_biased_stick,150,50,50,30,10000",
            prob_file=self.temp_bed_path
        )
        
        # Should be able to generate fragments
        frag_seq, start, end, length = gen.generate()
        self.assertIsInstance(frag_seq, str)
        self.assertGreater(len(frag_seq), 0)
    
    def test_very_tight_tolerance(self):
        """Test behavior with tight tolerance."""
        with open(self.temp_bed_path, 'w') as f:
            f.write("test_seq\t100\t300\t2.0\n")
        
        random.seed(99999)
        # Tight tolerance
        gen = AAVFragmentGenerator(
            self.test_seq,
            "gamma_biased_stick,200,30,50,10,10000",  # tolerance = 10bp
            prob_file=self.temp_bed_path
        )
        
        # Should still be able to generate fragments
        frag_seq, start, end, length = gen.generate()
        
        # Should be a valid fragment
        self.assertIsInstance(frag_seq, str)
        self.assertGreater(length, 0)
        self.assertIn(frag_seq, self.test_seq)
    
    def test_sequence_shorter_than_target_length(self):
        """Test behavior when sequence is shorter than target fragment length."""
        short_seq = "ATCG" * 10  # Only 40bp
        
        with open(self.temp_bed_path, 'w') as f:
            f.write("test_seq\t0\t40\t1.0\n")
        
        # Target length (mean=500bp, min=200bp) much longer than sequence (40bp)
        # With min_length=200, all sampled lengths will be >= 200bp, which is impossible
        # for a 40bp sequence. This should raise RuntimeError after 100 resamples.
        gen = AAVFragmentGenerator(
            short_seq,
            "gamma_biased_stick,500,50,200,50,10000",
            prob_file=self.temp_bed_path
        )
        
        # Should raise RuntimeError because parameters don't match sequence length
        with self.assertRaises(RuntimeError) as ctx:
            gen.generate()
        self.assertIn("100 length resamples", str(ctx.exception))


if __name__ == '__main__':
    unittest.main()
