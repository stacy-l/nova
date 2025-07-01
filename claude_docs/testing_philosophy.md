# Nova Testing Philosophy & Best Practices

## Core Principle

**Tests must validate actual algorithm behavior, not just pass.** A passing test that doesn't verify the intended functionality is worse than a failing test that exposes real issues.

## Testing Strategy

### 1. Extensive Mock Testing

Nova uses comprehensive mocking to avoid requiring large BAM files during development. However, mocks must be realistic and thorough.

#### Good Mock Example
```python
# Creates realistic genomic data distribution
chrom_read_pools = {
    "chr1": [(f"read_chr1_{i}", i * 5000) for i in range(800)],  # 800 reads
    "chr2": [(f"read_chr2_{i}", i * 5000) for i in range(400)],  # 400 reads
    "chr3": [(f"read_chr3_{i}", i * 5000) for i in range(300)],
}
```

#### Bad Mock Example
```python
# Too minimal - doesn't test real behavior
mock_reads = [("read1", 0), ("read2", 100)]  # Won't expose distribution issues
```

### 2. Coordinate-Aware Mocking

Use the `_create_coordinate_aware_mock_fetch()` pattern to simulate realistic BAM file behavior:

```python
def _create_coordinate_aware_mock_fetch(read_pools):
    def mock_fetch(chrom=None, start=None, end=None):
        if chrom in read_pools:
            for read_name, pos in read_pools[chrom]:
                if start <= pos < end:
                    yield _create_mock_read(read_name, chrom, pos)
    return mock_fetch
```

### 3. Algorithm Validation

Tests should verify the algorithm's logic, not just code execution:

```python
# Good: Verifies window limiting actually works
def test_window_limit_enforcement():
    selector = ReadSelector(max_reads_per_window=1)
    # Create overlapping reads in same window
    reads = create_clustered_reads()
    selected = selector.select_reads(reads, n=10)
    
    # Verify no two reads from same 1MB window
    windows = defaultdict(int)
    for read in selected:
        window = read.pos // 1_000_000
        windows[window] += 1
    
    assert all(count <= 1 for count in windows.values())
```

## Test Categories

### 1. Unit Tests
- Test individual methods in isolation
- Use mocks for external dependencies
- Fast execution (<1s per test)
- Location: `tests/test_*.py`

### 2. Integration Tests
- Test module interactions
- May use real BAM files if available
- Verify end-to-end workflows
- Location: `tests/test_integration.py`

### 3. Performance Tests
- Memory usage validation
- Scaling behavior
- Processing speed benchmarks
- Location: `tests/test_memory_optimization.py`

## Mock Testing Guidelines

### Do's
✅ Create hundreds of mock reads for realistic coverage
✅ Distribute reads across multiple chromosomes
✅ Include edge cases (chromosome boundaries, quality scores)
✅ Verify statistical properties of results
✅ Test both success and failure paths

### Don'ts
❌ Use minimal data that just makes tests pass
❌ Ignore coordinate space in mocks
❌ Test implementation details instead of behavior
❌ Assume mocks perfectly represent real data

## Real Data Testing

When BAM files are available (`tests/test_data/test_reads.bam`):

### Requirements
- Long reads (10-20kb)
- MAPQ ≥ 20
- Minimal soft clipping (<10%)
- Multiple chromosomes
- At least 100 mapped reads

### Test Scenarios
1. **Small scale**: 20 variants/reads
2. **Medium scale**: 100 variants with window limiting
3. **Large scale**: 500 reads with proportional sampling

## Debugging Failed Tests

### Step 1: Verify Mock Data
```python
# Check if mock has sufficient reads
print(f"Total mock reads: {sum(len(pool) for pool in chrom_read_pools.values())}")
# Should be >> n_reads requested
```

### Step 2: Check Coordinate Distribution
```python
# Verify reads span realistic genomic space
for chrom, reads in chrom_read_pools.items():
    positions = [pos for _, pos in reads]
    print(f"{chrom}: {min(positions)}-{max(positions)}")
```

### Step 3: Validate Algorithm Logic
- Add debug logging to trace selection process
- Verify filters are applied correctly
- Check randomization is working

## Test-Driven Development

### Process
1. Write test for new feature/bug fix
2. Verify test fails for the right reason
3. Implement minimal code to pass
4. Refactor while keeping tests green
5. Add edge case tests

### Example: Adding New Variant Type
```python
# 1. Write test first
def test_complex_variant_generation():
    generator = VariantGenerator()
    variant = generator.generate_complex_variant(
        type="inversion",
        length=500
    )
    assert len(variant.sequence) == 500
    assert variant.type == "inversion"
    # This will fail - feature doesn't exist yet

# 2. Implement feature
# 3. Verify test passes
# 4. Add edge cases
```

## Continuous Integration

### Pre-commit Checks
```bash
# Run before committing
python -m pytest tests/ -v
python -m pytest tests/ --cov=nova --cov-report=term-missing
```

### Coverage Requirements
- Aim for >80% code coverage
- 100% coverage for critical paths
- Document why uncovered code is acceptable

## Common Pitfalls

### 1. Reward Hacking
**Problem**: Test passes without testing intended behavior
**Solution**: Verify test fails when algorithm is broken

### 2. Brittle Tests
**Problem**: Tests break with minor refactoring
**Solution**: Test behavior, not implementation

### 3. Slow Tests
**Problem**: Test suite takes too long
**Solution**: Use mocks, parallelize, profile bottlenecks

### 4. Flaky Tests
**Problem**: Tests pass/fail randomly
**Solution**: Control randomness with seeds, avoid time dependencies

## Best Practices Summary

1. **Realistic mocks** > minimal mocks
2. **Behavior validation** > code coverage
3. **Fast feedback** > comprehensive slowness
4. **Clear failures** > cryptic errors
5. **Isolated tests** > interdependent suites

## Test Maintenance

### Regular Reviews
- Remove obsolete tests
- Update mocks when algorithms change
- Refactor duplicate test code
- Keep tests readable

### Documentation
- Document why complex tests exist
- Explain non-obvious assertions
- Link tests to issues/features

### Performance
- Monitor test suite execution time
- Profile slow tests
- Balance thoroughness with speed

## Conclusion

Good tests are an investment in code quality and development speed. They should give confidence that the code works correctly, catch regressions early, and serve as living documentation of intended behavior. When in doubt, err on the side of more realistic, thorough testing rather than superficial coverage.