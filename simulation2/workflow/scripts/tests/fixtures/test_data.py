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
Shared test fixtures and configuration for pytest/unittest.

This module provides common test data and utilities that can be
imported by any test module.
"""
import os
import tempfile


class TestSequences:
    """Standard test sequences for consistent testing."""
    
    # Simple repeating pattern sequences
    SIMPLE_100BP = "ACGT" * 25
    SIMPLE_400BP = "GCTA" * 100
    SIMPLE_1000BP = "ATCG" * 250
    SIMPLE_2000BP = "ATCG" * 500
    
    # More complex sequences
    GC_RICH = "GCGCGCGC" * 50  # 400bp, 100% GC
    AT_RICH = "ATATAT" * 50  # 300bp, 100% AT
    MIXED = "ACGTACGT" * 50  # 400bp, 50% GC


class TempBEDFile:
    """Context manager for creating temporary BED files."""
    
    def __init__(self, regions=None):
        """
        Create a temporary BED file.
        
        Args:
            regions: List of tuples (start, end, probability)
                    If None, creates an empty BED file.
        """
        self.regions = regions or []
        self.temp_file = None
        self.path = None
    
    def __enter__(self):
        """Create and write the temporary BED file."""
        self.temp_file = tempfile.NamedTemporaryFile(
            mode='w', 
            delete=False, 
            suffix='.bed'
        )
        self.path = self.temp_file.name
        
        # Write header
        self.temp_file.write("# Test BED file\n")
        
        # Write regions
        for start, end, prob in self.regions:
            self.temp_file.write(f"test_seq\t{start}\t{end}\t{prob}\n")
        
        self.temp_file.close()
        return self.path
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Clean up the temporary file."""
        if self.path and os.path.exists(self.path):
            os.remove(self.path)


def create_mock_fragment_lengths(lengths):
    """
    Create a mock FragmentLengths object with predetermined lengths.
    
    Args:
        lengths: List of predetermined fragment lengths
    
    Returns:
        Mock object with get_fragment_length() method
    """
    class MockFragmentLengths:
        def __init__(self, lengths):
            self.lengths = lengths
            self.index = 0
        
        def get_fragment_length(self):
            if self.index < len(self.lengths):
                length = self.lengths[self.index]
                self.index += 1
                return length
            return self.lengths[-1] if self.lengths else 100
    
    return MockFragmentLengths(lengths)


# Standard test parameters for AAV fragment generation
DEFAULT_AAV_PARAMS = {
    "gamma_stick": "gamma_stick,1000,500,25",
    "gamma_biased_stick": "gamma_biased_stick,1000,500,25,100,10000",
    "empirical_kde": "empirical_kde,25,2000,50,10000",
}

# Standard BED file configurations for testing
STANDARD_BED_CONFIGS = {
    "single_hotspot": [(400, 600, 3.0)],
    "single_coldspot": [(300, 700, 0.2)],
    "multiple_regions": [
        (0, 200, 0.5),
        (400, 600, 5.0),
        (800, 1000, 1.0)
    ],
    "adjacent_different": [
        (0, 500, 0.5),
        (500, 1000, 3.0)
    ],
    "tiny_hotspot": [(495, 505, 10.0)],
}
