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
import argparse
import random
import numpy as np



class DeletionGenerator:
    """Generate deletion sizes at integration junctions using configurable distributions."""
    
    def __init__(self, max_size, size_dist_str):
        """
        Initialize deletion size generator.
        
        Args:
            max_size: Maximum deletion size
            size_dist_str: Distribution specification string
        """
        self.max_size = max_size
        self.size_dist_str = size_dist_str
        self.parse_size_distribution()

    def parse_size_distribution(self):
        """Parse the size distribution string and set up the appropriate sampling method."""
        VALID_DISTRIBUTIONS = ["uniform", "triuniform"]
        
        if not self.size_dist_str:
            self.draw = self.zero
            return
        
        parts = self.size_dist_str.split(',')
        dist_name = parts[0]
        
        if dist_name not in VALID_DISTRIBUTIONS:
            raise ValueError(
                f"Unsupported deletion size distribution: '{dist_name}'. "
                f"Must be one of: {', '.join(VALID_DISTRIBUTIONS)}"
            )
        
        if dist_name == "uniform":
            if len(parts) != 1:
                raise ValueError(
                    f"uniform distribution requires 1 parameter (just the name), got {len(parts)}: {self.size_dist_str}"
                )
            self.draw = self.uniform
            
        elif dist_name == "triuniform":
            if len(parts) != 7:
                raise ValueError(
                    f"triuniform distribution requires 7 parameters (name,max1,max2,max3,prob1,prob2,prob3), "
                    f"got {len(parts)}: {self.size_dist_str}"
                )
            self.max1 = int(parts[1])
            self.max2 = int(parts[2])
            self.max3 = int(parts[3])
            if not (0 <= self.max1 <= self.max2 <= self.max3 <= self.max_size):
                raise ValueError(
                    f"triuniform distribution requires 0 <= max1 <= max2 <= max3 <= max_size, "
                    f"got max1={self.max1}, max2={self.max2}, max3={self.max3}, max_size={self.max_size}"
                )
            self.prob1 = float(parts[4])
            self.prob2 = float(parts[5])
            self.prob3 = float(parts[6])
            if not np.isclose(self.prob1 + self.prob2 + self.prob3, 1.0):
                raise ValueError(
                    f"triuniform distribution requires prob1 + prob2 + prob3 = 1.0, "
                    f"got prob1={self.prob1}, prob2={self.prob2}, prob3={self.prob3}"
                )
            self.draw = self.triuniform
    
    def zero(self):
        """Return zero deletion size."""
        return 0

    def uniform(self):
        """Sample deletion size from uniform distribution [0, max_size]."""
        return random.randint(0, self.max_size)

    def triuniform(self):
        """Sample deletion size from three-piece uniform distribution."""
        mode = random.choices([1, 2, 3], weights=[self.prob1, self.prob2, self.prob3], k=1)[0]
        if mode == 1:
            return random.randint(0, self.max1)
        elif mode == 2:
            return random.randint(self.max1, self.max2)
        else:
            return random.randint(self.max2, self.max3)



class InsertionGenerator:
    """Generate insertion sequences at integration junctions using configurable distributions."""
    
    def __init__(self, size_dist_str):
        """
        Initialize insertion sequence generator.
        
        Args:
            size_dist_str: Distribution specification string
        """
        self.size_dist_str = size_dist_str
        self.parse_size_distribution()

    def parse_size_distribution(self):
        """Parse the size distribution string and set up the appropriate sampling method."""
        VALID_DISTRIBUTIONS = ["uniform", "dualside_diuniform", "uniform_power"]
        
        if not self.size_dist_str:
            self.draw_length = self.zero_length
            self.draw = self.draw_dualside
            return
        
        parts = self.size_dist_str.split(',')
        dist_name = parts[0]
        
        if dist_name not in VALID_DISTRIBUTIONS:
            raise ValueError(
                f"Unsupported insertion size distribution: '{dist_name}'. "
                f"Must be one of: {', '.join(VALID_DISTRIBUTIONS)}"
            )

        if dist_name == "uniform":
            if len(parts) != 2:
                raise ValueError(
                    f"uniform distribution requires 2 parameters (name,max), got {len(parts)}: {self.size_dist_str}"
                )
            self.max1 = int(parts[1])
            self.draw_length = self.uniform_length
            self.draw = self.draw_dualside
            
        elif dist_name == "dualside_diuniform":
            if len(parts) != 5:
                raise ValueError(
                    f"dualside_diuniform distribution requires 5 parameters (name,max1,max2,prob1,prob2), "
                    f"got {len(parts)}: {self.size_dist_str}"
                )
            self.max1 = int(parts[1])
            self.max2 = int(parts[2])
            self.prob1 = float(parts[3])
            self.prob2 = float(parts[4])
            self.draw_length = self.dualside_diuniform_length
            self.draw = self.draw_dualside
            
        elif dist_name == "uniform_power":
            if len(parts) != 7:
                raise ValueError(
                    f"uniform_power distribution requires 7 parameters (name,max1,max2,a,b,prob1,prob2), "
                    f"got {len(parts)}: {self.size_dist_str}"
                )
            self.max1 = int(parts[1])
            self.max2 = int(parts[2])
            if self.max1 >= self.max2:
                raise ValueError(
                    f"uniform_power distribution requires max1 < max2, got max1={self.max1}, max2={self.max2}"
                )
            self.a = float(parts[3])
            self.b = float(parts[4])
            self.prob1 = float(parts[5])
            self.prob2 = float(parts[6])
            if self.prob1 + self.prob2 != 1.0:
                raise ValueError(
                    f"uniform_power distribution requires prob1 + prob2 = 1.0, got prob1={self.prob1}, prob2={self.prob2}"
                )
            # Precompute x and x_prob to avoid repeat calculations
            self.x = np.array(range(self.max1 + 1, self.max2 + 1))
            x_prob = self.a * self.x**self.b
            self.x_prob = x_prob / x_prob.sum()
            self.draw_length = self.uniform_power_length
            self.draw = self.draw_dualside

    def zero_length(self):
        """Return zero insertion length."""
        return 0

    def uniform_length(self):
        """Sample insertion length from uniform distribution [0, max1]."""
        return random.randint(0, self.max1)

    def dualside_diuniform_length(self):
        """Sample insertion length from two-piece uniform distribution."""
        mode = random.choices([1, 2], weights=[self.prob1, self.prob2], k=1)[0]
        if mode == 1:
            return random.randint(0, self.max1)
        else:
            return random.randint(self.max1, self.max2)

    def uniform_power_length(self):
        """Sample insertion length from uniform-power distribution.
        
        Uses a two-component mixture:
        - Component 1 (prob1): Uniform [0, max1]
        - Component 2 (prob2): Power-law weighted [max1+1, max2] with weights a*x^b
        """
        mode = random.choices([1, 2], weights=[self.prob1, self.prob2], k=1)[0]
        if mode == 1:
            return random.randint(0, self.max1)
        else:
            return np.random.choice(self.x, size=1, p=self.x_prob)[0]

    def draw_dualside(self):
        """
        Draw two independent random DNA insertion sequences.
        
        Returns:
            tuple: (left_insertion_seq, right_insertion_seq)
        """
        left_length = self.draw_length()
        right_length = self.draw_length()
        left_seq = ''.join(random.choices(['A', 'T', 'C', 'G'], k=left_length))
        right_seq = ''.join(random.choices(['A', 'T', 'C', 'G'], k=right_length))
        return (left_seq, right_seq)





