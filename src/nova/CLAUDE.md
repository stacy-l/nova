# CLAUDE.md - `nova` Source Code Guidelines

This file provides guidance for identifying and implementing reusable utilities within the `nova` codebase.

## Utility Identification Philosophy

When working with `nova`'s source code, look for opportunities to extract reusable functionality into dedicated utility modules. The goal is to create clean, testable, and extensible code that serves both current needs and anticipates future use cases.

## When to Create Utility Modules

### ✅ **Create a utility module when you find:**

1. **Cross-functional operations** that span multiple modules
   - Example: Sequence comparison logic used in testing, validation, and analysis

2. **Standalone algorithms** that could be useful beyond their original context
   - Example: FASTA content hashing for output validation → useful for variant calling validation

3. **Data processing functions** that operate on standard bioinformatics formats
   - Example: BAM file statistics, FASTA manipulation, VCF parsing

4. **Validation and comparison logic** that could apply to different data types
   - Example: Statistical comparisons, file format validation, data integrity checks

5. **Common analysis patterns** that appear in multiple workflows
   - Example: Coverage analysis, sequence similarity, performance benchmarking

### ❌ **Don't create utilities for:**

1. **Highly specific business logic** tightly coupled to `nova`'s simulation workflow
2. **One-off helper functions** used in a single location
3. **Configuration or initialization code** specific to particular modules
4. **UI/CLI presentation logic** that's tied to specific output formats

## Utility Module Design Principles

### 1. **Single Responsibility**
Each utility module should focus on one domain:
- `sequence_utils.py` - FASTA/sequence operations
- `file_utils.py` - General file I/O operations  
- `stats_utils.py` - Statistical analysis functions

### 2. **Framework Agnostic**
Utilities should work with standard data types and avoid deep dependencies:
```python
# Good: Works with standard types
def extract_sequence_content_hashes(fasta_path: str) -> Set[str]:

# Avoid: Tightly coupled to internal classes
def extract_hashes_from_registry(registry: VariantRegistry) -> Set[str]:
```

### 3. **Comprehensive Documentation**
Each utility function should include:
- Clear docstring with purpose and use cases
- Type hints for all parameters and returns
- Example usage in docstring
- Raised exceptions documented

### 4. **Future-Oriented Design**
Consider how the utility might be extended:
```python
# Current need: Compare two FASTA files
def compare_fasta_content(fasta1_path: str, fasta2_path: str) -> Dict[str, Any]:

# Future-oriented: Support multiple files
def compare_multiple_fasta_files(fasta_paths: List[str]) -> Dict[str, Any]:
```

## Current Utility Modules

### `sequence_utils.py`
**Purpose**: FASTA file analysis and sequence comparison operations

**Key Functions**:
- `extract_sequence_content_hashes()` - Hash sequence content for comparison
- `compare_fasta_content()` - Compare sequences between two FASTA files
- `analyze_sequence_recovery()` - Variant calling validation metrics
- `get_sequence_statistics()` - Basic sequence statistics

**Use Cases**:
- Output validation between different simulation methods
- Variant calling performance assessment  
- Sequence set analysis and comparison
- Quality control and testing

**Future Extensions**:
- Support for different hash algorithms
- Sequence clustering and similarity analysis
- Integration with variant calling pipelines
- Performance optimization for large files

## Implementation Checklist

When creating a new utility module:

- [ ] **Identify the core domain** (sequences, files, statistics, etc.)
- [ ] **Design clean function signatures** with standard types
- [ ] **Include comprehensive docstrings** with examples
- [ ] **Add type hints** for all parameters and returns
- [ ] **Handle edge cases** and provide meaningful error messages
- [ ] **Write unit tests** covering normal and edge cases
- [ ] **Consider future extensions** in the API design
- [ ] **Update this CLAUDE.md** with the new utility

## Integration Guidelines

### Importing Utilities
```python
# In source code modules
from nova.sequence_utils import compare_fasta_content

# In test files  
from nova.sequence_utils import extract_sequence_content_hashes

# In external scripts
sys.path.insert(0, 'path/to/nova/src')
from nova.sequence_utils import analyze_sequence_recovery
```

### Testing Utilities
- Unit tests go in `tests/test_<utility_name>.py`
- Integration tests can import utilities as needed
- Memory/performance scripts should use utilities for consistency

### CLI Integration
Consider whether utilities should be exposed via CLI commands:
```bash
# Future possibility
nova compare-sequences file1.fasta file2.fasta
nova analyze-recovery ground_truth.fasta detected.fasta
```

## Examples of Good Utility Identification

### ✅ **sequence_utils.py Creation**
**Original context**: Testing memory optimization equivalence
**Utility insight**: Sequence comparison is fundamental to bioinformatics validation
**Future applications**: 
- Variant calling validation
- Simulation comparison
- Quality control pipelines
- Cross-method benchmarking

### ✅ **Future Opportunities**

**BAM Processing Utilities** (`bam_utils.py`):
- Read quality statistics
- Coverage analysis
- Region sampling utilities
- Index validation

**Statistical Analysis** (`stats_utils.py`):
- Distribution comparisons
- Performance metrics calculation
- Significance testing
- Benchmarking utilities

**File Management** (`file_utils.py`):
- Safe file operations with cleanup
- Temporary directory management
- Path validation and normalization
- Format detection and validation

## Conclusion

The goal is to build a robust, extensible utility ecosystem that makes `nova` more powerful both as a simulation tool and as a foundation for genomics analysis workflows. When in doubt, err on the side of creating well-designed utilities that anticipate future needs.

**Remember**: Good utilities turn one-off solutions into reusable building blocks that accelerate future development and enable new use cases we haven't thought of yet.