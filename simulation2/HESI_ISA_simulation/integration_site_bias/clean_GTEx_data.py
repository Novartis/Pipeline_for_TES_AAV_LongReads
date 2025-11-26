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

import sys
import pandas as pd

# Read all lines from stdin
lines = sys.stdin.readlines()

# Find the header line starting with 'Name'
header_index = None
for i, line in enumerate(lines):
    if line.startswith('Name'):
        header_index = i
        break

if header_index is None:
    sys.exit("Error: No header line starting with 'Name' found.")

# Extract header and data lines
trimmed_lines = lines[header_index:]

# Load into pandas using a temporary buffer
from io import StringIO
data = StringIO(''.join(trimmed_lines))
df = pd.read_csv(data, sep='\t')

# Keep only required columns
required_cols = ['Name', 'Liver', 'Muscle_Skeletal']
missing_cols = [col for col in required_cols if col not in df.columns]
if missing_cols:
    sys.exit(f"Error: Missing columns: {missing_cols}")

df_filtered = df[required_cols].copy()
df_filtered["Name"] = [n.split('.')[0] for n in df_filtered["Name"]]
df_filtered = df_filtered.groupby('Name', as_index=False)[df_filtered.select_dtypes(include='number').columns].mean().copy()



# Print to stdout as TSV
print(df_filtered.to_csv(sep='\t', index=False))

