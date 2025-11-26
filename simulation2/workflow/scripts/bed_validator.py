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

import os


class BEDValidationError(Exception):
    """Exception raised for BED file validation errors."""
    pass


def validate_bed_probability_file(bed_file, file_description="BED probability file"):
    """
    Validate a BED-like file with probability values.
    
    Checks:
    1. File exists and is readable
    2. Format is correct (4 columns: name, start, end, probability)
    3. Start < End for all regions
    4. No overlapping regions within the same sequence
    5. Probability values are numeric and positive
    
    Args:
        bed_file: Path to the BED file
        file_description: Description of the file for error messages
    
    Raises:
        BEDValidationError: If validation fails
    
    Returns:
        dict: Dictionary mapping sequence names to list of (start, end, prob) tuples
    """
    if not bed_file or bed_file == "":
        # Empty file path is allowed (means no probability file specified)
        return {}
    
    # Check file exists
    if not os.path.exists(bed_file):
        raise BEDValidationError(
            f"{file_description} not found: {bed_file}"
        )
    
    # Check file is readable
    try:
        with open(bed_file, 'r') as f:
            pass
    except PermissionError:
        raise BEDValidationError(
            f"{file_description} is not readable: {bed_file}"
        )
    
    # Parse and validate file
    regions_by_seq = {}
    line_num = 0
    
    with open(bed_file, 'r') as f:
        for line in f:
            line_num += 1
            line = line.strip()
            
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue
            
            # Parse line
            parts = line.split('\t')
            if len(parts) < 4:
                raise BEDValidationError(
                    f"{file_description} line {line_num}: Expected 4 tab-separated columns "
                    f"(sequence_name, start, end, probability), got {len(parts)} columns"
                )
            
            seq_name = parts[0]
            
            # Validate start position
            try:
                start = int(parts[1])
            except ValueError:
                raise BEDValidationError(
                    f"{file_description} line {line_num}: Start position must be an integer, got '{parts[1]}'"
                )
            
            # Validate end position
            try:
                end = int(parts[2])
            except ValueError:
                raise BEDValidationError(
                    f"{file_description} line {line_num}: End position must be an integer, got '{parts[2]}'"
                )
            
            # Validate probability
            try:
                prob = float(parts[3])
            except ValueError:
                raise BEDValidationError(
                    f"{file_description} line {line_num}: Probability must be a number, got '{parts[3]}'"
                )
            
            # Check start < end
            if start >= end:
                raise BEDValidationError(
                    f"{file_description} line {line_num}: Start position ({start}) must be less than end position ({end})"
                )
            
            # Check probability is positive
            if prob <= 0:
                raise BEDValidationError(
                    f"{file_description} line {line_num}: Probability must be positive, got {prob}"
                )
            
            # Store region
            if seq_name not in regions_by_seq:
                regions_by_seq[seq_name] = []
            regions_by_seq[seq_name].append((start, end, prob, line_num))
    
    # Check for overlapping regions within each sequence
    for seq_name, regions in regions_by_seq.items():
        # Sort by start position
        sorted_regions = sorted(regions, key=lambda x: x[0])
        
        # Check each pair of consecutive regions
        for i in range(len(sorted_regions) - 1):
            curr_start, curr_end, curr_prob, curr_line = sorted_regions[i]
            next_start, next_end, next_prob, next_line = sorted_regions[i + 1]
            
            # Check if regions overlap
            if curr_end > next_start:
                raise BEDValidationError(
                    f"{file_description}: Overlapping regions detected for sequence '{seq_name}':\n"
                    f"  Region 1 (line {curr_line}): {curr_start}-{curr_end}\n"
                    f"  Region 2 (line {next_line}): {next_start}-{next_end}\n"
                    f"  Overlap: {next_start}-{curr_end}"
                )
    
    # Return regions without line numbers
    result = {}
    for seq_name, regions in regions_by_seq.items():
        result[seq_name] = [(start, end, prob) for start, end, prob, _ in regions]
    
    return result


def validate_integration_probability_file(bed_file):
    """
    Validate integration site probability BED file.
    
    Args:
        bed_file: Path to the BED file
    
    Raises:
        BEDValidationError: If validation fails
    
    Returns:
        dict: Dictionary mapping sequence names to list of (start, end, prob) tuples
    """
    return validate_bed_probability_file(bed_file, "Integration site probability file")


def validate_aav_breakpoint_probability_file(bed_file):
    """
    Validate AAV breakpoint probability BED file.
    
    Args:
        bed_file: Path to the BED file
    
    Raises:
        BEDValidationError: If validation fails
    
    Returns:
        dict: Dictionary mapping sequence names to list of (start, end, prob) tuples
    """
    return validate_bed_probability_file(bed_file, "AAV breakpoint probability file")
