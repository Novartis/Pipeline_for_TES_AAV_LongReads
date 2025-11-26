#!/usr/bin/env python
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

import re
import numpy as np


def cigar_longest_contiguous_match(cigar):
    """
    Find the longest contiguous match in a CIGAR string.
    
    Args:
        cigar: CIGAR string with = and X operators for match/mismatch
        
    Returns:
        int: Length of the longest contiguous match segment
    """
    matches = re.findall(r'(\d+)=', cigar)
    if not matches:
        return 0
    return max(int(match) for match in matches)


def cigar_minimum_random_mismatch(cigar, probe_len):
    """
    Find the minimum number of mismatches in any sliding window of length probe_len.
    
    Uses a sliding window approach over the alignment, tracking contiguous match/mismatch
    segments. Insertions and deletions reset the window tracking since they create
    discontinuities in the alignment.
    
    Args:
        cigar: CIGAR string with =/X operators for match/mismatch
        probe_len: Length of the sliding window (probe length)
        
    Returns:
        int: Minimum number of mismatches found in any probe_len window
    """
    min_mismatches = probe_len
    window_lengths = []
    mismatches = []
    offset = 0
    
    # Parse CIGAR string in reverse (right to left)
    cigar_split = re.split(r'([MIDNSH=X])', cigar)[:-1]
    
    while cigar_split:
        op = cigar_split.pop()
        num = int(cigar_split.pop())
        
        if op in '=X':
            # Match or mismatch - add to current window
            window_lengths.append(num)
            mismatches.append(0 if op == '=' else num)
        else:
            # Indel or other operation - reset window
            window_lengths = []
            mismatches = []
            offset = 0
            continue

        # Slide the window by removing segments from the left that are no longer needed
        while offset < len(window_lengths) - 1 and sum(window_lengths[offset + 1:]) >= probe_len:
            offset += 1

        # Check if we have a valid window of at least probe_len
        if sum(window_lengths[offset:]) >= probe_len:
            total_mismatches = sum(mismatches[offset:])
            min_mismatches = min(min_mismatches, total_mismatches)
    
    return min_mismatches



class PulldownProbability:
    """Calculate target enrichment pulldown probabilities from alignment CIGAR strings."""
    
    SUPPORTED_MODELS = ["twist", "twist_relaxed"]
    
    def __init__(self, model_name):
        """
        Initialize pulldown probability calculator.
        
        Args:
            model_name: Name of the enrichment model to use
            
        Raises:
            ValueError: If model_name is not supported
        """
        self.model_name = model_name
        self.parse_model_name()

    def parse_model_name(self):
        """Validate model name and set up the appropriate probability function."""
        if self.model_name.lower() not in self.SUPPORTED_MODELS:
            raise ValueError(
                f"Model '{self.model_name}' not recognized. "
                f"Supported models: {', '.join(self.SUPPORTED_MODELS)}"
            )

        if self.model_name.lower() == "twist":
            self.probability_function = self.twist_model
        elif self.model_name.lower() == "twist_relaxed":
            self.probability_function = self.twist_relaxed_model

    def twist_model(self, cigar, probe_len=120):
        """
        Calculate pulldown probability using the Twist Bioscience enrichment model.
        
        This model considers both contiguous matches and random mismatches within
        a probe length window to determine enrichment probability.
        
        Args:
            cigar: CIGAR string from alignment (must use = and X operators)
            probe_len: Length of enrichment probe (default: 120bp)
            
        Returns:
            float: Probability of successful pulldown [0, 1]
        """
        # Find longest contiguous match
        max_conti_len = cigar_longest_contiguous_match(cigar)
        min_conti_mm = max(probe_len - max_conti_len, 0)
        
        # Early return if perfect contiguous match
        p_conti = self.p_contiguous_mm(min_conti_mm)
        if p_conti >= 1.0:
            return 1.0

        # Find minimum mismatches in any sliding window
        min_random_mm = cigar_minimum_random_mismatch(cigar, probe_len)
        
        # Validate mismatch counts
        assert 0 <= min_conti_mm <= probe_len, f"Invalid contiguous mismatches: {min_conti_mm}"
        assert 0 <= min_random_mm <= probe_len, f"Invalid random mismatches: {min_random_mm}"
        
        # Return maximum probability from both models
        return max(p_conti, self.p_random_mm(min_random_mm))

    def p_random_mm(self, x, k=-0.293237, x0=10.198399, c=0.949740):
        """
        Calculate pulldown probability based on random mismatches.
        
        Uses a logistic function fitted to empirical Twist enrichment data.
        
        Args:
            x: Number of mismatches in probe_len window
            k: Logistic curve slope parameter
            x0: Logistic curve midpoint parameter
            c: Logistic curve offset parameter
            
        Returns:
            float: Probability [0, 1]
        """
        prob = 1.0 / (c + np.exp(-k * (x - x0)))
        return np.clip(prob, 0.0, 1.0)

    def p_contiguous_mm(self, x, a=-0.013703, b=1.022119):
        """
        Calculate pulldown probability based on contiguous match length.
        
        Uses a linear function fitted to empirical Twist enrichment data.
        
        Args:
            x: Number of mismatches (probe_len - longest_contiguous_match)
            a: Linear slope parameter
            b: Linear intercept parameter
            
        Returns:
            float: Probability [0, 1]
        """
        prob = a * x + b
        return np.clip(prob, 0.0, 1.0)
    
    def twist_relaxed_model(self, cigar, probe_len=120):
        """
        Calculate pulldown probability using the relaxed Twist enrichment model.
        
        This model uses more permissive parameters than the standard Twist model,
        allowing higher mismatch tolerance and capping maximum probability at 0.8.
        
        Args:
            cigar: CIGAR string from alignment (must use = and X operators)
            probe_len: Length of enrichment probe (default: 120bp)
            
        Returns:
            float: Probability of successful pulldown [0, 0.8]
        """
        # Find longest contiguous match
        max_conti_len = cigar_longest_contiguous_match(cigar)
        min_conti_mm = max(probe_len - max_conti_len, 0)
        
        # Calculate contiguous probability
        p_conti = self.p_contiguous_mm_relaxed(min_conti_mm)
        if p_conti >= 0.8:  # Max probability for relaxed model
            return 0.8

        # Find minimum mismatches in any sliding window
        min_random_mm = cigar_minimum_random_mismatch(cigar, probe_len)
        
        # Validate mismatch counts
        assert 0 <= min_conti_mm <= probe_len, f"Invalid contiguous mismatches: {min_conti_mm}"
        assert 0 <= min_random_mm <= probe_len, f"Invalid random mismatches: {min_random_mm}"
        
        # Return maximum probability from both models
        return max(p_conti, self.p_random_mm_relaxed(min_random_mm))
    
    def p_random_mm_relaxed(self, x, k=-0.08, x0=20, c=0.7, max_prob=0.8):
        """
        Calculate pulldown probability based on random mismatches (relaxed model).
        
        Uses a logistic function with more permissive parameters than the standard
        Twist model, allowing higher mismatch tolerance.
        
        Args:
            x: Number of mismatches in probe_len window
            k: Logistic curve slope parameter
            x0: Logistic curve midpoint parameter
            c: Logistic curve offset parameter
            max_prob: Maximum probability cap
            
        Returns:
            float: Probability [0, max_prob]
        """
        prob = 1.0 / (c + np.exp(-k * (x - x0)))
        return np.clip(prob, 0.0, max_prob)
    
    def p_contiguous_mm_relaxed(self, x, a=-0.01, b=1, max_prob=0.8):
        """
        Calculate pulldown probability based on contiguous match length (relaxed model).
        
        Uses a linear function with more permissive parameters than the standard
        Twist model.
        
        Args:
            x: Number of mismatches (probe_len - longest_contiguous_match)
            a: Linear slope parameter
            b: Linear intercept parameter
            max_prob: Maximum probability cap
            
        Returns:
            float: Probability [0, max_prob]
        """
        prob = a * x + b
        return np.clip(prob, 0.0, max_prob)

