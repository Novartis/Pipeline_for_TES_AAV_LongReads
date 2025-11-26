# Test Structure and Guidelines

## Directory Structure

```
workflow/scripts/
├── tests/
│   ├── __init__.py
│   ├── unit/                    # Unit tests for individual components
│   │   ├── __init__.py
│   │   ├── test_indel_generators.py       # DeletionGenerator, InsertionGenerator
│   │   ├── test_aav_fragment_generator.py # AAVFragmentGenerator API & error handling
│   │   ├── test_aav_sampling_methods.py   # Sampling behavior & distributions
│   │   └── test_fragment_integration.py   # fragment_integration function
│   ├── integration/             # Integration tests for full workflows
│   │   └── __init__.py
│   └── fixtures/                # Shared test data and utilities
│       ├── __init__.py
│       └── test_data.py         # Common test sequences, BED configs, mocks
├── run_tests.py                 # Main test runner script
└── [source files...]
```

## Test Organization

### Unit Tests (`tests/unit/`)
Each file tests a single module or class:

- **test_indel_generators.py** (7 tests)
  - `TestDeletionGenerator`: Tests deletion size sampling with various distributions
  - `TestInsertionGenerator`: Tests insertion sequence generation

- **test_aav_fragment_generator.py** (12 tests)
  - `TestAAVFragmentGenerator`: API tests, error handling, parameter validation, max_iter functionality

- **test_aav_sampling_methods.py** (16 tests)
  - `TestBiasedSampling`: Biased breakpoint sampling, hotspots, coldspots, weight calculations
  - `TestGammaStickSampling`: Stick-breaking algorithm, gamma distribution
  - `TestEdgeCasesAndErrors`: Edge cases, empty files, tight tolerances

- **test_fragment_integration.py** (tests for fragment_integration function)
  - Tests fragmentation of sequences containing AAV insertions

- **test_bed_validator.py** (20 tests)
  - `TestBEDValidator`: BED probability file validation, format checking, overlap detection
  - Tests file format, start/end validation, probability validation, overlap detection

### Integration Tests (`tests/integration/`)
Currently empty. Future integration tests will go here:
- End-to-end workflow tests
- Snakemake pipeline tests
- Multi-component integration tests

### Fixtures (`tests/fixtures/`)
Shared test utilities:

- **test_data.py**: Common test data and helper classes
  - `TestSequences`: Standard test sequences (100bp, 400bp, 1000bp, etc.)
  - `TempBEDFile`: Context manager for temporary BED files
  - `create_mock_fragment_lengths()`: Mock FragmentLengths objects
  - `DEFAULT_AAV_PARAMS`: Standard parameter strings
  - `STANDARD_BED_CONFIGS`: Pre-defined BED configurations

---

## Running Tests

### Using the Test Runner Script

The `run_tests.py` script provides convenient ways to run tests:

```bash
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

# Run with verbose output
python run_tests.py --verbose
```

### Using unittest Directly

```bash
# From the scripts directory:
cd workflow/scripts

# Run all tests
python -m unittest discover -s tests -p 'test_*.py'

# Run only unit tests
python -m unittest discover -s tests/unit -p 'test_*.py'

# Run a specific test file
python -m unittest tests.unit.test_indel_generators

# Run a specific test class
python -m unittest tests.unit.test_indel_generators.TestDeletionGenerator

# Run a specific test method
python -m unittest tests.unit.test_indel_generators.TestDeletionGenerator.test_uniform_within_bounds

# Run with verbose output
python -m unittest discover -s tests -p 'test_*.py' -v
```

### Quick Test Commands

```bash
# Run all tests (short form)
python -m unittest discover -s tests

# Run unit tests only
python -m unittest discover -s tests/unit

# Run integration tests only
python -m unittest discover -s tests/integration
```

---

## Adding New Tests

### Step 1: Choose Test Category

**Unit Test**: Tests a single function, class, or component in isolation
- Add to `tests/unit/`
- Name: `test_<module_name>.py` or `test_<component_name>.py`

**Integration Test**: Tests multiple components working together
- Add to `tests/integration/`
- Name: `test_<workflow_name>.py` or `test_<integration_scenario>.py`

### Step 2: Create Test File

Template for a new unit test file:

```python
#!/usr/bin/env python3
"""
Unit tests for [ComponentName].
Brief description of what is being tested.
"""
import unittest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

# Import what you're testing
from module_name import ComponentName

# Optionally import shared fixtures
from tests.fixtures.test_data import TestSequences, TempBEDFile


class TestComponentName(unittest.TestCase):
    """Test suite for ComponentName."""
    
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Initialize test data
        pass
    
    def tearDown(self):
        """Clean up after each test method."""
        # Clean up resources
        pass
    
    def test_basic_functionality(self):
        """Test basic usage of ComponentName."""
        # Arrange
        component = ComponentName()
        
        # Act
        result = component.do_something()
        
        # Assert
        self.assertEqual(result, expected_value)
    
    def test_error_handling(self):
        """Test that ComponentName handles errors correctly."""
        with self.assertRaises(ValueError):
            ComponentName(invalid_input)


if __name__ == '__main__':
    unittest.main()
```

### Step 3: Follow Naming Conventions

**Test Files**: `test_<what_is_tested>.py`
- Example: `test_aav_fragment_generator.py`

**Test Classes**: `Test<ClassName>` or `Test<Functionality>`
- Example: `TestAAVFragmentGenerator`, `TestBiasedSampling`

**Test Methods**: `test_<what_is_tested>`
- Use descriptive names that explain what is being tested
- Examples:
  - `test_uniform_within_bounds`
  - `test_single_hotspot_enrichment`
  - `test_empty_string_returns_zero`

### Step 4: Write Good Tests

**Structure**: Arrange-Act-Assert
```python
def test_something(self):
    """Describe what this test validates."""
    # Arrange: Set up test data and conditions
    generator = SomeGenerator(params)
    
    # Act: Execute the code being tested
    result = generator.generate()
    
    # Assert: Verify the result is correct
    self.assertEqual(result, expected)
```

**Use Descriptive Assertions**:
```python
# Good - clear message on failure
self.assertGreater(count, 10, 
                  f"Expected >10 fragments in hotspot, got {count}")

# Less helpful
self.assertGreater(count, 10)
```

**Test One Thing Per Test**:
```python
# Good - focused test
def test_returns_valid_length(self):
    result = generator.generate()
    self.assertGreater(len(result), 0)

def test_returns_string_type(self):
    result = generator.generate()
    self.assertIsInstance(result, str)

# Less ideal - testing multiple things
def test_generate(self):
    result = generator.generate()
    self.assertGreater(len(result), 0)
    self.assertIsInstance(result, str)
    self.assertIn('A', result)
```

### Step 5: Use Shared Fixtures

Leverage `tests/fixtures/test_data.py` for common test data:

```python
from tests.fixtures.test_data import TestSequences, TempBEDFile, DEFAULT_AAV_PARAMS

# Use standard test sequences
seq = TestSequences.SIMPLE_1000BP

# Use temporary BED files
with TempBEDFile([(400, 600, 3.0)]) as bed_path:
    generator = AAVFragmentGenerator(seq, method_str, prob_file=bed_path)
    result = generator.generate()

# Use standard parameters
method_str = DEFAULT_AAV_PARAMS["gamma_stick"]
```

### Step 6: Handle Temporary Files

Always clean up temporary files:

```python
import tempfile
import os

class TestWithTempFiles(unittest.TestCase):
    def setUp(self):
        """Create temporary file."""
        self.temp_file = tempfile.NamedTemporaryFile(delete=False)
        self.temp_path = self.temp_file.name
        self.temp_file.close()
    
    def tearDown(self):
        """Clean up temporary file."""
        if os.path.exists(self.temp_path):
            os.remove(self.temp_path)
```

Or use context managers:
```python
with TempBEDFile(regions) as bed_path:
    # Use bed_path
    pass
# File automatically cleaned up
```

---

## Test Guidelines

### Best Practices

1. **Independence**: Tests should not depend on each other
   - Each test should be runnable in isolation
   - Use `setUp()` and `tearDown()` for test-specific setup

2. **Reproducibility**: Use fixed random seeds for deterministic tests
   ```python
   def test_something(self):
       random.seed(12345)
       # Now randomized operations are reproducible
   ```

3. **Speed**: Keep tests fast
   - Unit tests should complete in milliseconds
   - Use small test data when possible
   - Mock expensive operations

4. **Clarity**: Tests are documentation
   - Use descriptive names
   - Add docstrings explaining what is tested
   - Include comments for complex test logic

5. **Coverage**: Test both success and failure paths
   - Normal operation
   - Edge cases (empty input, maximum values)
   - Error conditions (invalid input, missing files)

### What to Test

**For Algorithms**:
- Correct output for typical inputs
- Boundary conditions (min/max values)
- Edge cases (empty, single element)
- Statistical properties (distributions, means)

**For APIs**:
- Valid parameter combinations
- Invalid parameter combinations (error handling)
- Return types and structures
- Side effects

**For Error Handling**:
- Appropriate exception types
- Error messages are helpful
- Graceful degradation

### What NOT to Test

- External libraries (assume they work)
- Language features (Python built-ins)
- Obvious getters/setters without logic
- Auto-generated code

---

## Continuous Testing

### During Development

Run tests frequently while developing:

```bash
# Quick check - run related tests
python run_tests.py --test tests.unit.test_aav_fragment_generator

# Full check before committing
python run_tests.py --unit
```

### Before Committing

Always run the full test suite:

```bash
python run_tests.py
```

All tests should pass before committing code.

### Test-Driven Development (TDD)

Consider writing tests first:

1. Write a failing test that describes desired behavior
2. Implement the minimum code to make it pass
3. Refactor while keeping tests passing
4. Repeat

---

## Troubleshooting

### Import Errors

If you see `ModuleNotFoundError`:
- Check that `sys.path` is set correctly in your test file
- Verify you're running tests from the `scripts/` directory
- Ensure `__init__.py` files exist in all test directories

### Tests Pass Individually But Fail Together

- Tests may have hidden dependencies
- Check for shared state (module-level variables)
- Ensure proper `setUp()` and `tearDown()` in each test

### Flaky Tests

If tests sometimes pass and sometimes fail:
- Missing random seed
- Timing-dependent behavior
- Insufficiently strict assertions for stochastic operations
- External dependencies (network, files)

---

## Summary

**Test Structure**:
- `tests/unit/` - Unit tests (55 tests currently)
- `tests/integration/` - Integration tests (empty, for future)
- `tests/fixtures/` - Shared test data and utilities

**Running Tests**:
- Use `run_tests.py` for convenient test execution
- Use `python -m unittest` for standard unittest discovery
- Run specific tests by specifying module/class/method path

**Adding Tests**:
- Choose appropriate category (unit vs integration)
- Follow naming conventions
- Use shared fixtures from `tests/fixtures/test_data.py`
- Write clear, focused tests with good assertions
- Clean up resources in `tearDown()`

**Test Quality**:
- Independent, reproducible, fast
- Test success and failure paths
- Use descriptive names and docstrings
- Run tests frequently during development
