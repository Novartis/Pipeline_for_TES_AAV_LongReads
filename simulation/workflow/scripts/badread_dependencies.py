"""
This script contains only the minimal Badread dependencies required for fragment length distribution simulation.
Extracted from Badread by Ryan Wick (https://github.com/rrwick/Badread), GPLv3.
"""

import numpy as np
import scipy.special
import scipy.stats
import sys

# --- quickhist_gamma and helpers (minimal) ---
def quickhist_gamma(a, b, n50, height, output=sys.stderr):
    import math
    def get_max_width():
        try:
            import shutil
            return shutil.get_terminal_size((80, 20)).columns
        except Exception:
            return 80
    def draw_hist(y, shape, bins, height, x_tick_interval, y_label='', y_label_space=0,
                  print_labels=True, output=sys.stderr):
        max_count = max(y)
        normed_hist_list = [float(x) * height / max_count for x in y]
        i = 0
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
            i += 1
        # Draw X axis (labels omitted for brevity)
        print('  └' + '─' * bins, file=output)
    hist_max = int(np.ceil(n50 * 3 / 2000) * 2000)
    if get_max_width() > 120:
        bin_size = int(hist_max / 100)
    else:
        bin_size = int(hist_max / 50)
    bins = np.asarray([bin_size * (i + 1) for i in range(int(hist_max / bin_size))])
    frags_y, bases_y = [], []
    for b_val in bins:
        x = b_val - (bin_size / 2)
        frag_y = np.exp((-x*b) + ((a-1)*np.log(x)) + (a*np.log(b)) - scipy.special.gammaln(a))
        base_y = np.exp((-x*b) + (a*np.log(x)) + ((a+1)*np.log(b)) - scipy.special.gammaln(a+1))
        frags_y.append(frag_y)
        bases_y.append(base_y)
    shape = (0, hist_max)
    draw_hist(frags_y, shape, len(bins), height, 10, 'frags', 2, print_labels=False, output=output)
    draw_hist(bases_y, shape, len(bins), height, 10, 'bases', 2, output=output)

# --- float_to_str and print_in_two_columns (minimal) ---
def float_to_str(v, decimals=1, trim_zeros=False):
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
    part_1_len = max(len(l1p1), len(l2p1), len(l3p1)) + space_between
    format_str = '{:<' + str(part_1_len) + '}'
    l1p1 = format_str.format(l1p1)
    l2p1 = format_str.format(l2p1)
    l3p1 = format_str.format(l3p1)
    print(l1p1 + l1p2, file=output)
    print(l2p1 + l2p2, file=output)
    print(l3p1 + l3p2, file=output)

# --- FragmentLengths and helpers ---
class FragmentLengths(object):
    def __init__(self, mean, stdev, output=sys.stderr):
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
        if self.stdev == 0:
            return int(round(self.mean))
        else:
            fragment_length = int(round(np.random.gamma(self.gamma_k, self.gamma_t)))
            return max(fragment_length, 1)

def gamma_parameters(gamma_mean, gamma_stdev):
    gamma_a = (gamma_mean ** 2) / (gamma_stdev ** 2)
    gamma_b = gamma_mean / (gamma_stdev ** 2)
    gamma_k = (gamma_mean ** 2) / (gamma_stdev ** 2)
    gamma_t = (gamma_stdev ** 2) / gamma_mean
    return gamma_a, gamma_b, gamma_k, gamma_t

def find_n_value(a, b, n):
    target = 1.0 - (n / 100.0)
    bottom_range = 0.0
    top_range = 1.0
    while base_distribution_integral(a, b, top_range) < target:
        bottom_range = top_range
        top_range *= 2
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
    integral = 1.0 - np.exp(inc_gamma_ln(a+1, b*x) - scipy.special.gammaln(a+1))
    return integral

def inc_gamma_ln(a, b):
    return scipy.special.gammaln(a) + np.log(1-scipy.stats.gamma.cdf(b, a))
