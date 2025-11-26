#!/usr/bin/env python3
"""
Test runner for AAV integration simulation pipeline.

This script provides convenient ways to run tests with different configurations.
"""
import sys
import os
import unittest
import argparse


def run_all_tests(verbosity=2):
    """Run all tests in the test suite."""
    print("Running all tests...\n")
    loader = unittest.TestLoader()
    start_dir = os.path.join(os.path.dirname(__file__), 'tests')
    suite = loader.discover(start_dir, pattern='test_*.py')
    
    runner = unittest.TextTestRunner(verbosity=verbosity)
    result = runner.run(suite)
    
    return 0 if result.wasSuccessful() else 1


def run_unit_tests(verbosity=2):
    """Run only unit tests."""
    print("Running unit tests...\n")
    loader = unittest.TestLoader()
    start_dir = os.path.join(os.path.dirname(__file__), 'tests/unit')
    suite = loader.discover(start_dir, pattern='test_*.py')
    
    runner = unittest.TextTestRunner(verbosity=verbosity)
    result = runner.run(suite)
    
    return 0 if result.wasSuccessful() else 1


def run_integration_tests(verbosity=2):
    """Run only integration tests."""
    print("Running integration tests...\n")
    loader = unittest.TestLoader()
    start_dir = os.path.join(os.path.dirname(__file__), 'tests/integration')
    suite = loader.discover(start_dir, pattern='test_*.py')
    
    runner = unittest.TextTestRunner(verbosity=verbosity)
    result = runner.run(suite)
    
    return 0 if result.wasSuccessful() else 1


def run_specific_test(test_path, verbosity=2):
    """
    Run a specific test module, class, or method.
    
    Args:
        test_path: Path in format 'module' or 'module.Class' or 'module.Class.method'
    """
    print(f"Running test: {test_path}\n")
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromName(test_path)
    
    runner = unittest.TextTestRunner(verbosity=verbosity)
    result = runner.run(suite)
    
    return 0 if result.wasSuccessful() else 1


def list_tests():
    """List all available tests."""
    print("Available tests:\n")
    
    test_dir = os.path.join(os.path.dirname(__file__), 'tests')
    
    for root, dirs, files in os.walk(test_dir):
        for file in files:
            if file.startswith('test_') and file.endswith('.py'):
                rel_path = os.path.relpath(os.path.join(root, file), test_dir)
                print(f"  - {rel_path}")
    
    return 0


def main():
    """Main entry point for test runner."""
    parser = argparse.ArgumentParser(
        description='Test runner for AAV integration simulation pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run all tests
  python run_tests.py

  # Run only unit tests
  python run_tests.py --unit

  # Run only integration tests
  python run_tests.py --integration

  # Run a specific test module
  python run_tests.py --test tests.unit.test_indel_generators

  # Run a specific test class
  python run_tests.py --test tests.unit.test_indel_generators.TestDeletionGenerator

  # Run a specific test method
  python run_tests.py --test tests.unit.test_indel_generators.TestDeletionGenerator.test_uniform_within_bounds

  # List all available tests
  python run_tests.py --list

  # Run with minimal output
  python run_tests.py --quiet
        """
    )
    
    parser.add_argument(
        '--unit', '-u',
        action='store_true',
        help='Run only unit tests'
    )
    
    parser.add_argument(
        '--integration', '-i',
        action='store_true',
        help='Run only integration tests'
    )
    
    parser.add_argument(
        '--test', '-t',
        type=str,
        help='Run a specific test (module, class, or method)'
    )
    
    parser.add_argument(
        '--list', '-l',
        action='store_true',
        help='List all available tests'
    )
    
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Minimal output'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose output'
    )
    
    args = parser.parse_args()
    
    # Determine verbosity
    if args.quiet:
        verbosity = 0
    elif args.verbose:
        verbosity = 2
    else:
        verbosity = 1
    
    # Add scripts directory to path
    sys.path.insert(0, os.path.dirname(__file__))
    
    # Execute requested test mode
    if args.list:
        return list_tests()
    elif args.test:
        return run_specific_test(args.test, verbosity)
    elif args.unit:
        return run_unit_tests(verbosity)
    elif args.integration:
        return run_integration_tests(verbosity)
    else:
        return run_all_tests(verbosity)


if __name__ == '__main__':
    sys.exit(main())
