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

import sys
import random
import os
import numpy as np
from scipy.stats import gaussian_kde
from badread_dependencies import FragmentLengths
try:
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend
    import matplotlib.pyplot as plt
    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False


class AAVFragmentGenerator:
    """
    Generate AAV fragments using different methods specified by a configuration string.
    
    Currently supports:
    - gamma_stick,mean,stdev,min: Generate fragment length from gamma distribution, 
                                   then use stick-breaking to randomly extract a fragment 
                                   of that length from the AAV sequence.
                                   Example: "gamma_stick,1000,500,25"
    
    - gamma_biased_stick,mean,stdev,min,tolerance,max_iter: Generate fragment length from 
                                   gamma distribution, then sample two breakpoints using 
                                   biased sampling from BED file. Repeat until fragment 
                                   length is within tolerance of target. If max_iter attempts 
                                   fail, resample a new target length. Raises exception after 
                                   100 length resamples (indicates parameter mismatch).
                                   Example: "gamma_biased_stick,1000,500,25,50,10000"
                                   Note: Requires prob_file parameter to be set.
    
    - empirical_kde,min,max,tolerance,max_iter: Generate fragment length by sampling from
                                   a Kernel Density Estimation (KDE) fitted to empirical
                                   size data, then use stick-breaking to extract a fragment.
                                   Breakpoints are uniform random (does NOT use prob_file).
                                   Example: "empirical_kde,25,5000,50,10000"
                                   Note: Requires empirical_sizes parameter to be set.
    
    - empirical_kde_biased,min,max,tolerance,max_iter: Generate fragment length by sampling from
                                   KDE fitted to empirical sizes, then sample two breakpoints using
                                   biased sampling from BED file. Repeat until fragment length is
                                   within tolerance of target. Similar to gamma_biased_stick but
                                   uses empirical KDE instead of gamma distribution.
                                   Example: "empirical_kde_biased,25,5000,50,10000"
                                   Note: Requires both empirical_sizes and prob_file parameters.
    """
    
    def __init__(self, aav_seq, method_str, verbose_output=None, prob_file=None, empirical_sizes=None, seed=None):
        """
        Initialize the AAV fragment generator.
        
        Args:
            aav_seq: The full AAV sequence string
            method_str: Method specification string with parameters 
                       (e.g., "gamma_stick,1000,500,25" for mean=1000, stdev=500, min=25)
            verbose_output: File handle for verbose output (default: None)
            prob_file: Path to BED file with breakpoint probabilities for gamma_biased_stick (default: None)
            empirical_sizes: Path to file with empirical fragment sizes for empirical_kde (default: None)
            seed: Random seed for reproducibility (default: None)
        """
        self.aav_seq = aav_seq
        self.method_str = method_str
        self.verbose_output = verbose_output if verbose_output else open('/dev/null', 'w')
        self.prob_file = prob_file
        self.empirical_sizes = empirical_sizes
        self.seed = seed
        
        # Set random seed if provided
        if self.seed is not None:
            random.seed(self.seed)
            np.random.seed(self.seed)
        
        self.parse_method()
    
    def parse_method(self):
        """Parse the method string and set up the appropriate generation method."""
        if not self.method_str or not self.method_str.strip():
            raise ValueError("Empty AAV fragment generation method string")
        
        parts = self.method_str.split(',')
        
        if len(parts) == 0:
            raise ValueError("Empty AAV fragment generation method string")
        
        method_name = parts[0]
        
        if method_name == "gamma_stick":
            if len(parts) != 4:
                raise ValueError(
                    f"gamma_stick requires 4 parameters: method,mean,stdev,min. "
                    f"Got {len(parts)} parts: {self.method_str}"
                )
            self.mean_length = int(parts[1])
            self.std_length = int(parts[2])
            self.min_length = int(parts[3])
            self.generate = self.gamma_stick
            # Initialize the fragment length generator for gamma distribution
            self.length_generator = FragmentLengths(self.mean_length, self.std_length, self.verbose_output)
        
        elif method_name == "gamma_biased_stick":
            if len(parts) != 6:
                raise ValueError(
                    f"gamma_biased_stick requires 6 parameters: method,mean,stdev,min,tolerance,max_iter. "
                    f"Got {len(parts)} parts: {self.method_str}"
                )
            self.mean_length = int(parts[1])
            self.std_length = int(parts[2])
            self.min_length = int(parts[3])
            self.tolerance = int(parts[4])
            self.max_iter = int(parts[5])
            
            # prob_file must be provided separately
            if not self.prob_file:
                raise ValueError(
                    "gamma_biased_stick requires prob_file parameter to be set. "
                    "Pass prob_file argument to AAVFragmentGenerator constructor."
                )
            
            self.prob_bed_file = self.prob_file
            self.generate = self.gamma_biased_stick
            # Initialize the fragment length generator for gamma distribution
            self.length_generator = FragmentLengths(self.mean_length, self.std_length, self.verbose_output)
            # Load probability regions from BED file
            self.prob_regions = self._read_probability_bed(self.prob_bed_file)
            self.baseline_prob = 1.0
            # Calculate total weight for the AAV sequence
            self._setup_biased_sampling()
        
        elif method_name == "empirical_kde":
            if len(parts) != 5:
                raise ValueError(
                    f"empirical_kde requires 5 parameters: method,min,max,tolerance,max_iter. "
                    f"Got {len(parts)} parts: {self.method_str}"
                )
            self.min_length = int(parts[1])
            self.max_length = int(parts[2])
            self.tolerance = int(parts[3])
            self.max_iter = int(parts[4])
            
            # empirical_sizes file must be provided separately
            if not self.empirical_sizes:
                raise ValueError(
                    "empirical_kde requires empirical_sizes parameter to be set. "
                    "Pass empirical_sizes argument to AAVFragmentGenerator constructor."
                )
            
            self.generate = self.empirical_kde
            # Load empirical size data and fit KDE
            self._load_empirical_sizes()
            self._fit_kde()
        
        elif method_name == "empirical_kde_biased":
            if len(parts) != 5:
                raise ValueError(
                    f"empirical_kde_biased requires 5 parameters: method,min,max,tolerance,max_iter. "
                    f"Got {len(parts)} parts: {self.method_str}"
                )
            self.min_length = int(parts[1])
            self.max_length = int(parts[2])
            self.tolerance = int(parts[3])
            self.max_iter = int(parts[4])
            
            # Both empirical_sizes and prob_file must be provided
            if not self.empirical_sizes:
                raise ValueError(
                    "empirical_kde_biased requires empirical_sizes parameter to be set. "
                    "Pass empirical_sizes argument to AAVFragmentGenerator constructor."
                )
            if not self.prob_file:
                raise ValueError(
                    "empirical_kde_biased requires prob_file parameter to be set. "
                    "Pass prob_file argument to AAVFragmentGenerator constructor."
                )
            
            self.prob_bed_file = self.prob_file
            self.generate = self.empirical_kde_biased
            # Load empirical size data and fit KDE
            self._load_empirical_sizes()
            self._fit_kde()
            # Load probability regions from BED file
            self.prob_regions = self._read_probability_bed(self.prob_bed_file)
            self.baseline_prob = 1.0
            # Calculate total weight for the AAV sequence
            self._setup_biased_sampling()
        
        else:
            raise ValueError(f"Unsupported AAV fragment generation method: {method_name}")
    
    def gamma_stick(self):
        """
        Generate AAV fragment using gamma-distributed length and stick-breaking.
        
        Returns:
            tuple: (fragment_sequence, start_position, end_position, fragment_length)
                - fragment_sequence: The AAV fragment string
                - start_position: 1-based start position in original AAV sequence
                - end_position: 1-based end position in original AAV sequence
                - fragment_length: Length of the fragment
        """
        # Sample fragment length from gamma distribution, enforce minimum
        while True:
            frag_len = self.length_generator.get_fragment_length()
            if frag_len >= self.min_length:
                break
        
        # Generate fragment using stick-breaking
        frag_seq = self._break_seq_to_length(self.aav_seq, frag_len)
        
        # Find position in original sequence
        start_pos = self.aav_seq.index(frag_seq) + 1  # 1-based
        end_pos = start_pos + len(frag_seq) - 1
        
        return (frag_seq, start_pos, end_pos, len(frag_seq))
    
    def _break_seq_to_length(self, seq, target_length):
        """
        Randomly break a sequence using stick-breaking until a fragment of target_length is found.
        
        Args:
            seq: Sequence to break
            target_length: Desired fragment length
        
        Returns:
            A fragment of exactly target_length from seq
        """
        if target_length >= len(seq):
            return seq
        
        while True:
            fragments = self._break_seq(seq, target_length)
            if fragments:
                return random.choice(fragments)
    
    def _break_seq(self, seq, k):
        """
        Recursively break a sequence into pieces, returning all pieces of length k.
        
        Args:
            seq: Sequence to break
            k: Target length for fragments
        
        Returns:
            List of all fragments with length exactly k
        """
        # Base case: if the sequence is too short, return an empty list
        if len(seq) < k:
            return []
        
        # Randomly select a breaking point
        break_point = random.randint(1, len(seq) - 1)
        
        # Define the two pieces
        piece1 = seq[0:break_point]
        piece2 = seq[break_point:]
        
        # Initialize a list to store fragments of length k
        fragments = []
        
        # Check if any piece has length k and add to fragments list,
        # else start recursive process of further breaking
        if len(piece1) == k:
            fragments.append(piece1)
        else:
            fragments.extend(self._break_seq(piece1, k))
        
        if len(piece2) == k:
            fragments.append(piece2)
        else:
            fragments.extend(self._break_seq(piece2, k))
        
        return fragments
    
    def _read_probability_bed(self, bed_file):
        """
        Read a BED file with sampling probabilities for AAV breakpoints.
        Format: sequence_name\tstart\tend\tprobability
        Returns a list of (start, end, prob) tuples sorted by start position.
        """
        if not os.path.exists(bed_file):
            raise FileNotFoundError(f"Probability BED file not found: {bed_file}")
        
        prob_regions = []
        with open(bed_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 4:
                    # We expect format: sequence_name, start, end, probability
                    # For AAV, we only care about start, end, prob (ignore sequence name)
                    start = int(parts[1])
                    end = int(parts[2])
                    prob = float(parts[3])
                    prob_regions.append((start, end, prob))
        
        # Sort regions by start position for efficient lookup
        prob_regions.sort(key=lambda x: x[0])

        # Make sure the max end does not exceed the AAV sequence length (if regions exist):
        if prob_regions and prob_regions[-1][1] > len(self.aav_seq):
            raise ValueError(
                f"Probability BED file contains region end {prob_regions[-1][1]} "
                f"that exceeds AAV sequence length {len(self.aav_seq)}."
            )

        return prob_regions
    
    def _setup_biased_sampling(self):
        """
        Calculate cumulative weights for biased sampling across the AAV sequence.
        This allows for efficient weighted sampling of breakpoints.
        """
        aav_length = len(self.aav_seq)
        total_weight = 0.0
        
        if self.prob_regions:
            last_end = 0
            
            for start, end, prob in self.prob_regions:
                # Add baseline weight for gap before this region
                if start > last_end:
                    total_weight += (start - last_end) * self.baseline_prob
                
                # Add custom weight for this region
                total_weight += (end - start) * prob
                last_end = end
            
            # Add baseline weight for remaining part of sequence
            if last_end < aav_length:
                total_weight += (aav_length - last_end) * self.baseline_prob
        else:
            # Entire sequence has baseline probability
            total_weight = aav_length * self.baseline_prob
        
        self.total_weight = total_weight
    
    def _sample_breakpoint(self):
        """
        Sample a single breakpoint position in the AAV sequence using biased sampling.
        Returns a 0-based position in the AAV sequence.
        """
        aav_length = len(self.aav_seq)
        target_weight = random.uniform(0, self.total_weight)
        cumulative_weight = 0.0
        
        if self.prob_regions:
            last_end = 0
            
            for start, end, prob in self.prob_regions:
                # Check baseline region before this custom region
                if start > last_end:
                    gap_weight = (start - last_end) * self.baseline_prob
                    if cumulative_weight + gap_weight >= target_weight:
                        # Position is in this baseline gap
                        offset = (target_weight - cumulative_weight) / self.baseline_prob
                        return int(last_end + offset)
                    cumulative_weight += gap_weight
                
                # Check custom probability region
                region_weight = (end - start) * prob
                if cumulative_weight + region_weight >= target_weight:
                    # Position is in this custom region
                    offset = (target_weight - cumulative_weight) / prob
                    return int(start + offset)
                cumulative_weight += region_weight
                last_end = end
            
            # Check remaining baseline region
            if last_end < aav_length:
                gap_weight = (aav_length - last_end) * self.baseline_prob
                if cumulative_weight + gap_weight >= target_weight:
                    offset = (target_weight - cumulative_weight) / self.baseline_prob
                    return int(last_end + offset)
        else:
            # Entire sequence is baseline - uniform sampling
            offset = target_weight / self.baseline_prob
            return int(offset)
        
        # Fallback (shouldn't happen)
        return int(aav_length * random.random())
    
    def gamma_biased_stick(self):
        """
        Generate AAV fragment using gamma-distributed length and biased breakpoint sampling.
        
        Process:
        1. Sample target fragment length from gamma distribution
        2. Sample two breakpoints using biased sampling from BED file
        3. Check if resulting fragment length is within tolerance of target
        4. Repeat steps 2-3 up to max_iter times
        5. If max_iter reached without success, go back to step 1 (resample length)
        6. Raise exception if more than 100 length resamples occur (indicates parameter mismatch)
        
        Returns:
            tuple: (fragment_sequence, start_position, end_position, fragment_length)
                - fragment_sequence: The AAV fragment string
                - start_position: 1-based start position in original AAV sequence
                - end_position: 1-based end position in original AAV sequence
                - fragment_length: Length of the fragment
        
        Raises:
            RuntimeError: If more than 100 length resamples are needed, indicating a mismatch
                         between length sampling parameters and AAV breakpoint probabilities
        """
        max_length_resamples = 100
        length_resample_count = 0
        
        while length_resample_count < max_length_resamples:
            # Sample target fragment length from gamma distribution, enforce minimum
            target_length = self._sample_valid_target_length()
            if target_length is None:
                length_resample_count += 1
                continue
            
            # Try to find breakpoints that yield a fragment within tolerance
            result = self._find_matching_breakpoints(target_length)
            if result is not None:
                return result
            
            # No suitable breakpoints found - resample length
            length_resample_count += 1
            if self.verbose_output and self.verbose_output != open('/dev/null', 'w'):
                print(f"Warning: Could not find fragment within tolerance after {self.max_iter} attempts "
                    f"for target length {target_length}. Resampling length (reset #{length_resample_count}).",
                    file=sys.stderr)
        
        # If we've exceeded max_length_resamples, there's likely a parameter mismatch
        raise RuntimeError(
            f"Failed to generate AAV fragment after {max_length_resamples} length resamples. "
            f"This indicates a mismatch between the length sampling parameters "
            f"(mean={self.mean_length}, stdev={self.std_length}, min={self.min_length}, "
            f"tolerance={self.tolerance}) and the AAV sequence length ({len(self.aav_seq)}bp) "
            f"and/or breakpoint probability distribution. Consider adjusting these parameters."
        )
    
    def _sample_valid_target_length(self):
        """
        Sample a target fragment length from gamma distribution that meets minimum requirements.
        
        Returns:
            int: Valid target length, or None if couldn't find one after max_iter attempts
        """
        for _ in range(self.max_iter):
            target_length = self.length_generator.get_fragment_length()
            if target_length >= self.min_length:
                return target_length
        
        if self.verbose_output and self.verbose_output != open('/dev/null', 'w'):
            print(f"Warning: Could not find fragment larger than {self.min_length} "
                f"after {self.max_iter} attempts.", file=sys.stderr)
        return None
    
    def _find_matching_breakpoints(self, target_length):
        """
        Find two breakpoints that yield a fragment within tolerance of target length.
        
        Args:
            target_length: Desired fragment length
        
        Returns:
            tuple: (fragment_seq, start_pos, end_pos, actual_length) or None if no match found
        """
        for _ in range(self.max_iter):
            # Sample two breakpoints
            bp1 = self._sample_breakpoint()
            bp2 = self._sample_breakpoint()
            
            # Ensure bp1 < bp2
            if bp1 > bp2:
                bp1, bp2 = bp2, bp1
            
            # Skip if breakpoints are too close (would give empty fragment)
            if bp1 == bp2:
                continue
            
            # Calculate actual fragment length
            actual_length = bp2 - bp1
            
            # Check if within tolerance of target length
            if abs(actual_length - target_length) <= self.tolerance:
                # Extract fragment (0-based slicing)
                frag_seq = self.aav_seq[bp1:bp2]
                
                # Return 1-based positions for consistency with gamma_stick
                return (frag_seq, bp1 + 1, bp2, actual_length)
        
        return None
    
    def _load_empirical_sizes(self):
        """
        Load empirical fragment sizes from file.
        Expected format: one integer per line representing fragment size.
        """
        if self.empirical_sizes is None:
            raise ValueError("No empirical sizes file specified")

        if not os.path.exists(self.empirical_sizes):
            raise FileNotFoundError(f"Empirical sizes file not found: {self.empirical_sizes}")
        
        sizes = []
        with open(self.empirical_sizes, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                try:
                    size = int(line)
                    if size > 0:  # Only keep positive sizes
                        sizes.append(size)
                except ValueError:
                    print(f"Warning: Skipping invalid line in empirical sizes file: {line}", 
                          file=sys.stderr)
        
        if not sizes:
            raise ValueError(f"No valid sizes found in empirical sizes file: {self.empirical_sizes}")
        
        self.empirical_data = np.array(sizes)
        print(f"Loaded {len(self.empirical_data)} empirical fragment sizes", file=self.verbose_output)
        print(f"  min={np.min(self.empirical_data)}, max={np.max(self.empirical_data)}, "
              f"median={np.median(self.empirical_data):.0f}, mean={np.mean(self.empirical_data):.0f}",
              file=self.verbose_output)
    
    def _fit_kde(self):
        """
        Fit a Kernel Density Estimation to the empirical size data.
        Uses Scott's rule for bandwidth selection by default.
        """
        print("Fitting KDE to empirical fragment size distribution...", file=self.verbose_output)
        self.kde = gaussian_kde(self.empirical_data, bw_method='scott')
        print(f"  KDE bandwidth: {self.kde.factor * self.empirical_data.std():.2f}", 
              file=self.verbose_output)
    
    def empirical_kde(self):
        """
        Generate AAV fragment by sampling length from empirical KDE distribution.
        
        Process:
        1. Sample target fragment length from KDE
        2. Enforce min/max length constraints
        3. Use stick-breaking to randomly extract fragment of that length
        
        Returns:
            tuple: (fragment_sequence, start_position, end_position, fragment_length)
                - fragment_sequence: The AAV fragment string
                - start_position: 1-based start position in original AAV sequence
                - end_position: 1-based end position in original AAV sequence
                - fragment_length: Length of the fragment
        
        Raises:
            RuntimeError: If cannot find valid fragment length after max_iter attempts
        """
        # Sample fragment length from KDE, enforce min/max constraints
        for attempt in range(self.max_iter):
            # Sample from KDE
            frag_len = int(round(self.kde.resample(1)[0, 0]))
            
            # Check constraints
            if self.min_length <= frag_len <= self.max_length:
                # Also ensure it doesn't exceed AAV sequence length
                if frag_len <= len(self.aav_seq):
                    break
        else:
            raise RuntimeError(
                f"Failed to sample valid fragment length from KDE after {self.max_iter} attempts. "
                f"Constraints: min={self.min_length}, max={self.max_length}, "
                f"AAV_seq_len={len(self.aav_seq)}. "
                f"Consider adjusting min/max parameters or checking empirical size distribution."
            )
        
        # Generate fragment using stick-breaking
        frag_seq = self._break_seq_to_length(self.aav_seq, frag_len)
        
        # Find position in original sequence
        start_pos = self.aav_seq.index(frag_seq) + 1  # 1-based
        end_pos = start_pos + len(frag_seq) - 1
        
        return (frag_seq, start_pos, end_pos, len(frag_seq))
    
    def empirical_kde_biased(self):
        """
        Generate AAV fragment using KDE-sampled length and biased breakpoint sampling.
        
        Process:
        1. Sample target fragment length from KDE fitted to empirical data
        2. Sample two breakpoints using biased sampling from BED file
        3. Check if resulting fragment length is within tolerance of target
        4. Repeat steps 2-3 up to max_iter times
        5. If max_iter reached without success, go back to step 1 (resample length)
        6. Raise exception if more than 100 length resamples occur (indicates parameter mismatch)
        
        Returns:
            tuple: (fragment_sequence, start_position, end_position, fragment_length)
                - fragment_sequence: The AAV fragment string
                - start_position: 1-based start position in original AAV sequence
                - end_position: 1-based end position in original AAV sequence
                - fragment_length: Length of the fragment
        
        Raises:
            RuntimeError: If more than 100 length resamples are needed, indicating a mismatch
                         between length sampling parameters and AAV breakpoint probabilities
        """
        max_length_resamples = 100
        length_resample_count = 0
        
        while length_resample_count < max_length_resamples:
            # Sample target fragment length from KDE, enforce constraints
            target_length = self._sample_valid_kde_length()
            if target_length is None:
                length_resample_count += 1
                continue
            
            # Try to find breakpoints that yield a fragment within tolerance
            result = self._find_matching_breakpoints(target_length)
            if result is not None:
                return result
            
            # No suitable breakpoints found - resample length
            length_resample_count += 1
            if self.verbose_output and self.verbose_output != open('/dev/null', 'w'):
                print(f"Warning: Could not find fragment within tolerance after {self.max_iter} attempts "
                    f"for target length {target_length}. Resampling length (reset #{length_resample_count}).",
                    file=sys.stderr)
        
        # If we've exceeded max_length_resamples, there's likely a parameter mismatch
        raise RuntimeError(
            f"Failed to generate AAV fragment after {max_length_resamples} length resamples. "
            f"This indicates a mismatch between the KDE-sampled lengths "
            f"(min={self.min_length}, max={self.max_length}, tolerance={self.tolerance}) "
            f"and the AAV sequence length ({len(self.aav_seq)}bp) "
            f"and/or breakpoint probability distribution. Consider adjusting these parameters."
        )
    
    def _sample_valid_kde_length(self):
        """
        Sample a target fragment length from KDE that meets min/max requirements.
        
        Returns:
            int: Valid target length, or None if couldn't find one after max_iter attempts
        """
        for _ in range(self.max_iter):
            # Sample from KDE
            target_length = int(round(self.kde.resample(1)[0, 0]))
            
            # Check constraints
            if self.min_length <= target_length <= self.max_length:
                # Also ensure it doesn't exceed AAV sequence length
                if target_length <= len(self.aav_seq):
                    return target_length
        
        if self.verbose_output and self.verbose_output != open('/dev/null', 'w'):
            print(f"Warning: Could not sample valid fragment length from KDE "
                f"(min={self.min_length}, max={self.max_length}) "
                f"after {self.max_iter} attempts.", file=sys.stderr)
        return None
    
    def plot_empirical_kde(self, output_file, title="AAV Fragment Size Distribution"):
        """
        Plot the empirical data histogram and fitted KDE curve.
        
        Args:
            output_file: Path to save the plot (e.g., 'kde_plot.png')
            title: Title for the plot
        
        Raises:
            RuntimeError: If matplotlib is not available or method is not empirical_kde/empirical_kde_biased
        """
        if not PLOTTING_AVAILABLE:
            raise RuntimeError(
                "Matplotlib is not available. Install matplotlib to use plotting functionality."
            )
        
        if not hasattr(self, 'empirical_data') or not hasattr(self, 'kde'):
            raise RuntimeError(
                "plot_empirical_kde can only be called for empirical_kde or empirical_kde_biased methods. "
                "Make sure the generator was initialized with one of these methods."
            )
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Plot histogram of empirical data
        ax.hist(self.empirical_data, bins=50, density=True, alpha=0.6, 
                color='skyblue', edgecolor='black', label='Empirical Data')
        
        # Plot KDE curve
        x_range = np.linspace(
            max(0, self.empirical_data.min() - 200),
            self.empirical_data.max() + 200,
            1000
        )
        kde_values = self.kde(x_range)
        ax.plot(x_range, kde_values, 'r-', linewidth=2, label='KDE Fit')
        
        # Add vertical lines for min/max constraints if they exist
        if hasattr(self, 'min_length'):
            ax.axvline(self.min_length, color='green', linestyle='--', 
                      linewidth=1.5, alpha=0.7, label=f'Min: {self.min_length} bp')
        if hasattr(self, 'max_length'):
            ax.axvline(self.max_length, color='orange', linestyle='--',
                      linewidth=1.5, alpha=0.7, label=f'Max: {self.max_length} bp')
        
        # Formatting
        ax.set_xlabel('Fragment Size (bp)', fontsize=12)
        ax.set_ylabel('Density', fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.legend(loc='best', fontsize=10)
        ax.grid(True, alpha=0.3, linestyle=':')
        
        # Add statistics text box
        stats_text = (
            f'n = {len(self.empirical_data)}\n'
            f'Mean = {np.mean(self.empirical_data):.1f} bp\n'
            f'Median = {np.median(self.empirical_data):.1f} bp\n'
            f'Std = {np.std(self.empirical_data):.1f} bp'
        )
        ax.text(0.98, 0.97, stats_text, transform=ax.transAxes,
               fontsize=10, verticalalignment='top', horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Save figure
        plt.tight_layout()
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"Plot saved to: {output_file}", file=self.verbose_output)
