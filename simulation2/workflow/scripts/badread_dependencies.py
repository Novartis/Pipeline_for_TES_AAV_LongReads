"""
This script contains only the minimal Badread dependencies required for fragment length distribution simulation.
Extracted from Badread by Ryan Wick (https://github.com/rrwick/Badread), GPLv3.

This module provides gamma distribution-based fragment length sampling for sequencing simulation.
"""

import numpy as np
import scipy.special
import scipy.stats
import sys


class FragmentLengths:
    """
    Generate fragment lengths from a gamma distribution or constant value.
    
    This class handles fragment length sampling for sequencing simulations,
    supporting both constant lengths and gamma-distributed lengths.
    
    Attributes:
        mean: Mean fragment length in base pairs
        stdev: Standard deviation of fragment length distribution
        gamma_k: Gamma distribution shape parameter (k)
        gamma_t: Gamma distribution scale parameter (theta)
    """
    
    def __init__(self, mean, stdev, output=sys.stderr):
        """
        Initialize fragment length generator.
        
        Args:
            mean: Mean fragment length in base pairs
            stdev: Standard deviation of fragment length (0 for constant length)
            output: File handle for verbose output (default: sys.stderr)
        """
        self.mean = mean
        self.stdev = stdev
        
        print('', file=output)
        if self.stdev == 0:
            self.gamma_k, self.gamma_t = None, None
            print(f'Using a constant fragment length of {mean} bp', file=output)
        else:
            print('Generating fragment lengths from a gamma distribution:', file=output)
            gamma_a, gamma_b, self.gamma_k, self.gamma_t = gamma_parameters(mean, stdev)
            n50 = int(round(find_n_value(gamma_a, gamma_b, 50)))
            print_in_two_columns(f'  mean  = {float_to_str(mean):>6} bp',
                                 f'  stdev = {float_to_str(stdev):>6} bp',
                                 f'  N50   = {n50:>6} bp',
                                 'parameters:',
                                 f'  k (shape)     = {self.gamma_k:.4e}',
                                 f'  theta (scale) = {self.gamma_t:.4e}',
                                 output=output)
            quickhist_gamma(gamma_a, gamma_b, n50, 8, output=output)
    
    def get_fragment_length(self):
        """
        Sample a fragment length.
        
        Returns:
            int: Fragment length in base pairs (minimum 1bp)
        """
        if self.stdev == 0:
            return int(round(self.mean))
        else:
            fragment_length = int(round(np.random.gamma(self.gamma_k, self.gamma_t)))
            return max(fragment_length, 1)

def gamma_parameters(gamma_mean, gamma_stdev):
    """
    Calculate gamma distribution parameters from mean and standard deviation.
    
    Args:
        gamma_mean: Desired mean of the distribution
        gamma_stdev: Desired standard deviation of the distribution
    
    Returns:
        tuple: (gamma_a, gamma_b, gamma_k, gamma_t) where:
            - gamma_a: Shape parameter (alpha)
            - gamma_b: Rate parameter (beta)
            - gamma_k: Shape parameter (k, same as gamma_a)
            - gamma_t: Scale parameter (theta)
    """
    gamma_a = (gamma_mean ** 2) / (gamma_stdev ** 2)
    gamma_b = gamma_mean / (gamma_stdev ** 2)
    gamma_k = (gamma_mean ** 2) / (gamma_stdev ** 2)
    gamma_t = (gamma_stdev ** 2) / gamma_mean
    return gamma_a, gamma_b, gamma_k, gamma_t


def find_n_value(a, b, n):
    """
    Find the N-value (e.g., N50) for a gamma distribution.
    
    The N-value is the length threshold where n% of the total bases
    are in fragments of that length or longer.
    
    Args:
        a: Gamma shape parameter (alpha)
        b: Gamma rate parameter (beta)
        n: Percentage for N-value (e.g., 50 for N50)
    
    Returns:
        float: The N-value length threshold
    """
    target = 1.0 - (n / 100.0)
    bottom_range = 0.0
    top_range = 1.0
    
    # Find upper bound
    while base_distribution_integral(a, b, top_range) < target:
        bottom_range = top_range
        top_range *= 2
    
    # Binary search for precise value
    guess = (bottom_range + top_range) / 2.0
    while True:
        integral = base_distribution_integral(a, b, guess)
        if top_range - bottom_range < 0.01:
            return guess
        elif integral < target:
            bottom_range = guess
            guess = (bottom_range + top_range) / 2.0
        else:
            top_range = guess
            guess = (bottom_range + top_range) / 2.0


def base_distribution_integral(a, b, x):
    """
    Calculate the cumulative distribution for the base-weighted gamma distribution.
    
    Args:
        a: Gamma shape parameter (alpha)
        b: Gamma rate parameter (beta)
        x: Upper limit for integration
    
    Returns:
        float: Cumulative probability up to x
    """
    integral = 1.0 - np.exp(inc_gamma_ln(a+1, b*x) - scipy.special.gammaln(a+1))
    return integral


def inc_gamma_ln(a, b):
    """
    Calculate the natural log of the incomplete gamma function.
    
    Args:
        a: Shape parameter
        b: Upper limit
    
    Returns:
        float: ln(incomplete_gamma(a, b))
    """
    return scipy.special.gammaln(a) + np.log(1-scipy.stats.gamma.cdf(b, a))


def float_to_str(v, decimals=1, trim_zeros=False):
    """
    Convert a float to a string with specified decimal places.
    
    Args:
        v: Value to convert
        decimals: Number of decimal places (default: 1)
        trim_zeros: Whether to trim trailing zeros (default: False)
    
    Returns:
        str: Formatted string representation
    """
    if float(int(v)) == v:
        return str(int(v))
    else:
        formatter = '%.' + str(decimals) + 'f'
        result = formatter % v
        if trim_zeros:
            while result.endswith('0'):
                result = result[:-1]
        return result


def print_in_two_columns(l1p1, l2p1, l3p1, l1p2, l2p2, l3p2, output, space_between=6):
    """
    Print three lines in two-column format.
    
    Args:
        l1p1, l2p1, l3p1: Left column content for lines 1, 2, 3
        l1p2, l2p2, l3p2: Right column content for lines 1, 2, 3
        output: File handle for output
        space_between: Space between columns (default: 6)
    """
    part_1_len = max(len(l1p1), len(l2p1), len(l3p1)) + space_between
    format_str = '{:<' + str(part_1_len) + '}'
    l1p1 = format_str.format(l1p1)
    l2p1 = format_str.format(l2p1)
    l3p1 = format_str.format(l3p1)
    print(l1p1 + l1p2, file=output)
    print(l2p1 + l2p2, file=output)
    print(l3p1 + l3p2, file=output)


def quickhist_gamma(a, b, n50, height, output=sys.stderr):
    """
    Display a quick histogram of fragment and base distributions.
    
    Args:
        a: Gamma shape parameter (alpha)
        b: Gamma rate parameter (beta)
        n50: N50 value for scaling
        height: Height of histogram in terminal lines
        output: File handle for output (default: sys.stderr)
    """
    import math
    
    def get_max_width():
        """Get terminal width for histogram display."""
        try:
            import shutil
            return shutil.get_terminal_size((80, 20)).columns
        except Exception:
            return 80
    
    def draw_hist(y, shape, bins, height, x_tick_interval, y_label='', y_label_space=0,
                  print_labels=True, output=sys.stderr):
        """Draw a histogram to the terminal."""
        max_count = max(y)
        normed_hist_list = [float(x) * height / max_count for x in y]
        
        # Draw bars
        for depth in range(height-1, -1, -1):
            print(' ', end='', file=output)
            print(' |', end='', file=output)
            for item in normed_hist_list:
                floored_item = math.floor(item)
                if floored_item >= depth:
                    print('█', end='', file=output)
                else:
                    print(' ', end='', file=output)
            print('', file=output)
        
        # Draw X axis
        print('  └' + '─' * bins, file=output)
    
    # Calculate histogram parameters
    hist_max = int(np.ceil(n50 * 3 / 2000) * 2000)
    if get_max_width() > 120:
        bin_size = int(hist_max / 100)
    else:
        bin_size = int(hist_max / 50)
    
    bins = np.asarray([bin_size * (i + 1) for i in range(int(hist_max / bin_size))])
    frags_y, bases_y = [], []
    
    # Calculate distribution values for each bin
    for b_val in bins:
        x = b_val - (bin_size / 2)
        frag_y = np.exp((-x*b) + ((a-1)*np.log(x)) + (a*np.log(b)) - scipy.special.gammaln(a))
        base_y = np.exp((-x*b) + (a*np.log(x)) + ((a+1)*np.log(b)) - scipy.special.gammaln(a+1))
        frags_y.append(frag_y)
        bases_y.append(base_y)
    
    shape = (0, hist_max)
    draw_hist(frags_y, shape, len(bins), height, 10, 'frags', 2, print_labels=False, output=output)
    draw_hist(bases_y, shape, len(bins), height, 10, 'bases', 2, output=output)
