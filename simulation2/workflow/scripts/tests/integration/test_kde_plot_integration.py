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
Test script to verify KDE plotting integration with add_integration.py
"""

import subprocess
import tempfile
import os
import sys

def test_kde_plot_integration():
    """Test that add_integration.py generates KDE plot correctly."""
    
    print("=" * 70)
    print("TESTING KDE PLOT INTEGRATION WITH add_integration.py")
    print("=" * 70)
    
    # Create temporary files
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create empirical sizes file
        empirical_file = os.path.join(tmpdir, "empirical_sizes.txt")
        with open(empirical_file, 'w') as f:
            f.write("# Test empirical sizes\n")
            for size in range(400, 1200, 50):
                for _ in range(5):
                    f.write(f"{size}\n")
        
        # Create test genome fragment (simulating stdin input)
        test_input = ">chr1:1000-5000\nACGTACGTACGT" * 100 + "\n"
        
        # Output files
        events_tsv = os.path.join(tmpdir, "events.tsv")
        kde_plot = os.path.join(tmpdir, "kde_plot.png")
        
        # Build command
        cmd = [
            "python", "workflow/scripts/add_integration.py",
            "--fq_out",
            "--events", "1",
            "--mfl", "2000",
            "--sfl", "500",
            "--aav_method", "empirical_kde,300,1500,50,10000",
            "--aav_empirical_sizes", empirical_file,
            "--aav_kde_plot", kde_plot,
            "--aav", "test_run/test_AAV.fa",
            "--int_out", events_tsv,
            "--del_max", "100",
            "--del_flag", "triuniform,0,50,100,0.5,0.4,0.1",
            "--ins_flag", "dualside_diuniform,0,100,0.5,0.5",
            "--p_trans", "0.0",
            "--exp", "1",
            "--seed", "42",
            "--verbose"
        ]
        
        print("\n[1] Running add_integration.py with empirical_kde and KDE plot...")
        print(f"    Command: {' '.join(cmd[:10])}...")
        
        try:
            # Run the command
            result = subprocess.run(
                cmd,
                input=test_input,
                capture_output=True,
                text=True,
                cwd="/Users/davidkr2/bitbucket/davidkr2/development_AAV_integration/read_simulation",
                timeout=30
            )
            
            if result.returncode != 0:
                print(f"    ✗ Command failed with return code {result.returncode}")
                print(f"    stderr: {result.stderr[:500]}")
                return False
            
            print("    ✓ Command executed successfully")
            
        except subprocess.TimeoutExpired:
            print("    ✗ Command timed out")
            return False
        except Exception as e:
            print(f"    ✗ Error running command: {e}")
            return False
        
        # Check if KDE plot was created
        print("\n[2] Checking if KDE plot was created...")
        if os.path.exists(kde_plot):
            file_size = os.path.getsize(kde_plot)
            print(f"    ✓ Plot file exists: {kde_plot}")
            print(f"    ✓ File size: {file_size} bytes")
            
            if file_size > 1000:  # Should be at least 1KB for a valid PNG
                print("    ✓ File size looks reasonable")
            else:
                print("    ✗ File size too small, may be corrupted")
                return False
        else:
            print(f"    ✗ Plot file not found: {kde_plot}")
            return False
        
        # Check if events TSV was created
        print("\n[3] Checking if integration events TSV was created...")
        if os.path.exists(events_tsv):
            print(f"    ✓ Events TSV exists: {events_tsv}")
        else:
            print(f"    ✗ Events TSV not found: {events_tsv}")
            return False
        
        # Check stderr for KDE plot message
        print("\n[4] Checking for KDE plot message in output...")
        if "KDE plot saved to:" in result.stderr or "Plot saved to:" in result.stderr:
            print("    ✓ Found plot generation message in stderr")
        else:
            print("    ⚠ Plot generation message not found (may have gone to stdout)")
        
        print("\n" + "=" * 70)
        print("✓ ALL TESTS PASSED!")
        print("=" * 70)
        print("\nKDE plot integration is working correctly.")
        print(f"Test plot was generated at: {kde_plot}")
        return True

if __name__ == "__main__":
    success = test_kde_plot_integration()
    sys.exit(0 if success else 1)
