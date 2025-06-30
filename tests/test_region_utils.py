"""
Tests for region filtering utilities.
"""

import pytest
import tempfile
import os
from pathlib import Path
from nova.region_utils import RegionFilter, validate_bed_file, load_multiple_bed_files


@pytest.fixture
def sample_bed_content():
    """Sample BED file content for testing."""
    return """chr1\t1000\t2000\tregion1
chr1\t5000\t6000\tregion2
chr2\t3000\t4000\tregion3
chr2\t7000\t8000\tregion4"""


@pytest.fixture
def sample_bed_file(sample_bed_content):
    """Create a temporary BED file for testing."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write(sample_bed_content)
        f.flush()
        yield f.name
    os.unlink(f.name)


@pytest.fixture
def overlapping_bed_content():
    """BED content with overlapping regions for merge testing."""
    return """chr1\t1000\t2000\tregion1
chr1\t1500\t2500\tregion2
chr1\t5000\t6000\tregion3"""


@pytest.fixture
def overlapping_bed_file(overlapping_bed_content):
    """Create a temporary BED file with overlapping regions."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write(overlapping_bed_content)
        f.flush()
        yield f.name
    os.unlink(f.name)


class TestRegionFilter:
    """Test RegionFilter functionality."""
    
    def test_empty_filter(self):
        """Test empty filter behavior."""
        filter_obj = RegionFilter()
        assert len(filter_obj) == 0
        assert not filter_obj
        assert not filter_obj.overlaps_read('chr1', 1000, 2000)
        assert not filter_obj.contains_point('chr1', 1500)
    
    def test_load_from_bed(self, sample_bed_file):
        """Test loading regions from BED file."""
        filter_obj = RegionFilter(bed_path=sample_bed_file)
        assert len(filter_obj) == 4
        assert bool(filter_obj)
        
        stats = filter_obj.get_statistics()
        assert stats['total_regions'] == 4
        assert stats['chromosomes'] == 2
        assert 'chr1' in stats['regions_per_chromosome']
        assert 'chr2' in stats['regions_per_chromosome']
        assert stats['regions_per_chromosome']['chr1'] == 2
        assert stats['regions_per_chromosome']['chr2'] == 2
    
    def test_overlaps_read(self, sample_bed_file):
        """Test read overlap detection."""
        filter_obj = RegionFilter(bed_path=sample_bed_file)
        
        # Test overlapping reads
        assert filter_obj.overlaps_read('chr1', 1500, 1600)  # Inside region1
        assert filter_obj.overlaps_read('chr1', 900, 1100)   # Overlaps start of region1
        assert filter_obj.overlaps_read('chr1', 1900, 2100)  # Overlaps end of region1
        assert filter_obj.overlaps_read('chr2', 3500, 4500)  # Spans region3
        
        # Test non-overlapping reads
        assert not filter_obj.overlaps_read('chr1', 2500, 3000)  # Between regions
        assert not filter_obj.overlaps_read('chr2', 500, 1000)   # Before regions
        assert not filter_obj.overlaps_read('chr3', 1000, 2000)  # Different chromosome
        assert not filter_obj.overlaps_read('chr1', 8500, 9000)  # After regions
    
    def test_contains_point(self, sample_bed_file):
        """Test point containment."""
        filter_obj = RegionFilter(bed_path=sample_bed_file)
        
        # Test points inside regions
        assert filter_obj.contains_point('chr1', 1500)  # Inside region1
        assert filter_obj.contains_point('chr2', 3500)  # Inside region3
        
        # Test boundary conditions (BED is 0-based, half-open)
        assert filter_obj.contains_point('chr1', 1000)  # Start of region1 (inclusive)
        assert not filter_obj.contains_point('chr1', 2000)  # End of region1 (exclusive)
        
        # Test points outside regions
        assert not filter_obj.contains_point('chr1', 2500)  # Between regions
        assert not filter_obj.contains_point('chr3', 1500)  # Different chromosome
    
    def test_get_chromosome_names(self, sample_bed_file):
        """Test chromosome name extraction."""
        filter_obj = RegionFilter(bed_path=sample_bed_file)
        chromosomes = filter_obj.get_chromosome_names()
        assert sorted(chromosomes) == ['chr1', 'chr2']
    
    def test_filter_regions_by_chromosome(self, sample_bed_file):
        """Test chromosome-specific region filtering."""
        filter_obj = RegionFilter(bed_path=sample_bed_file)
        
        chr1_regions = filter_obj.filter_regions_by_chromosome('chr1')
        assert len(chr1_regions) == 2
        
        chr2_regions = filter_obj.filter_regions_by_chromosome('chr2')
        assert len(chr2_regions) == 2
        
        chr3_regions = filter_obj.filter_regions_by_chromosome('chr3')
        assert len(chr3_regions) == 0
    
    def test_merge_overlapping(self, overlapping_bed_file):
        """Test merging of overlapping regions."""
        filter_obj = RegionFilter(bed_path=overlapping_bed_file)
        assert len(filter_obj) == 3  # Before merge
        
        filter_obj.merge_overlapping()
        assert len(filter_obj) == 2  # After merge (first two should be merged)
        
        # Test that merged region covers the original span
        assert filter_obj.overlaps_read('chr1', 1000, 1100)  # Original region1 start
        assert filter_obj.overlaps_read('chr1', 2400, 2500)  # Original region2 end
        assert filter_obj.overlaps_read('chr1', 5500, 5600)  # Original region3 (separate)


class TestValidateBedFile:
    """Test BED file validation functionality."""
    
    def test_valid_bed_file(self, sample_bed_file):
        """Test validation of a valid BED file."""
        is_valid, errors = validate_bed_file(sample_bed_file)
        assert is_valid
        assert len(errors) == 0
    
    def test_nonexistent_file(self):
        """Test validation of non-existent file."""
        is_valid, errors = validate_bed_file('/nonexistent/file.bed')
        assert not is_valid
        assert len(errors) == 1
        assert 'not found' in errors[0]
    
    def test_invalid_bed_format(self):
        """Test validation of invalid BED format."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            f.write("chr1\t-100\t2000\tregion1\n")  # Negative start
            f.write("chr2\t5000\t4000\tregion2\n")   # End before start
            f.flush()
            
            try:
                is_valid, errors = validate_bed_file(f.name)
                assert not is_valid
                assert len(errors) >= 1
            finally:
                os.unlink(f.name)
    
    def test_empty_bed_file(self):
        """Test validation of empty BED file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            f.write("")  # Empty file
            f.flush()
            
            try:
                is_valid, errors = validate_bed_file(f.name)
                assert not is_valid
                assert any('pyranges error' in error for error in errors)
            finally:
                os.unlink(f.name)


class TestLoadMultipleBedFiles:
    """Test loading and combining multiple BED files."""
    
    def test_load_multiple_files(self, sample_bed_file):
        """Test loading multiple BED files."""
        # Create a second BED file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            f.write("chr3\t1000\t2000\tregion5\nchr3\t3000\t4000\tregion6\n")
            f.flush()
            second_bed_file = f.name
        
        try:
            combined_filter = load_multiple_bed_files([sample_bed_file, second_bed_file])
            assert len(combined_filter) == 6  # 4 from first file + 2 from second
            
            chromosomes = combined_filter.get_chromosome_names()
            assert sorted(chromosomes) == ['chr1', 'chr2', 'chr3']
            
        finally:
            os.unlink(second_bed_file)
    
    def test_load_empty_list(self):
        """Test loading empty list of BED files."""
        combined_filter = load_multiple_bed_files([])
        assert len(combined_filter) == 0
        assert not combined_filter
    
    def test_load_with_invalid_file(self, sample_bed_file):
        """Test loading with one invalid file."""
        with pytest.raises(Exception):
            load_multiple_bed_files([sample_bed_file, '/nonexistent/file.bed'])


class TestIntegration:
    """Integration tests for region filtering."""
    
    def test_read_selection_workflow(self, sample_bed_file):
        """Test typical workflow for read selection."""
        # Load target regions
        target_filter = RegionFilter(bed_path=sample_bed_file)
        
        # Create exclusion regions
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            f.write("chr1\t1800\t1900\texclude1\n")  # Overlaps with region1
            f.flush()
            exclusion_file = f.name
        
        try:
            exclusion_filter = RegionFilter(bed_path=exclusion_file)
            
            # Test read that overlaps target but not exclusion
            read_chr, read_start, read_end = 'chr1', 1200, 1300
            assert target_filter.overlaps_read(read_chr, read_start, read_end)
            assert not exclusion_filter.overlaps_read(read_chr, read_start, read_end)
            
            # Test read that overlaps both target and exclusion
            read_chr, read_start, read_end = 'chr1', 1850, 1950
            assert target_filter.overlaps_read(read_chr, read_start, read_end)
            assert exclusion_filter.overlaps_read(read_chr, read_start, read_end)
            
            # Test read that overlaps neither
            read_chr, read_start, read_end = 'chr1', 3000, 3100
            assert not target_filter.overlaps_read(read_chr, read_start, read_end)
            assert not exclusion_filter.overlaps_read(read_chr, read_start, read_end)
            
        finally:
            os.unlink(exclusion_file)