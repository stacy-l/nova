"""
Tests for RegionAwareWindowGenerator functionality.
"""

import pytest
import os
import pyranges as pr
from nova.region_utils import RegionFilter, RegionAwareWindowGenerator


@pytest.fixture
def test_bam_path():
    """Path to test BAM file."""
    return "tests/test_data/test_reads.bam"


@pytest.fixture
def test_centromeres_bed():
    """Path to test centromeres BED file."""
    return "tests/test_data/test_centromeres.bed"


@pytest.fixture
def test_superdups_bed():
    """Path to test genomic super duplications BED file."""
    return "tests/test_data/test_genomic_superdups.bed"


class TestRegionAwareWindowGenerator:
    """Test RegionAwareWindowGenerator class."""
    
    def test_init(self, test_bam_path):
        """Test initialization."""
        generator = RegionAwareWindowGenerator(test_bam_path)
        assert generator.bam_path == test_bam_path
        assert generator._valid_regions_cache == {}
        assert generator._genome_regions is None
    
    def test_load_genome_regions(self, test_bam_path):
        """Test loading genome regions from real BAM file."""
        generator = RegionAwareWindowGenerator(test_bam_path)
        regions = generator._load_genome_regions()
        
        # Should have multiple chromosomes/contigs
        assert len(regions) > 0
        
        # Check that regions start at 0
        df = regions.df
        assert all(df['Start'] == 0)
        
        # Should include chr1 and chr2 at minimum
        chromosomes = set(df['Chromosome'])
        assert 'chr1' in chromosomes
        assert 'chr2' in chromosomes
        
        # Main chromosomes should be reasonably large
        chr1_length = df[df['Chromosome'] == 'chr1']['End'].iloc[0]
        chr2_length = df[df['Chromosome'] == 'chr2']['End'].iloc[0]
        assert chr1_length > 100_000_000  # chr1 should be > 100MB
        assert chr2_length > 100_000_000  # chr2 should be > 100MB
    
    def test_make_cache_key(self, test_centromeres_bed, test_superdups_bed):
        """Test cache key generation."""
        generator = RegionAwareWindowGenerator('test.bam')
        
        # Test with different filter combinations
        target_filter = RegionFilter(bed_path=test_superdups_bed)
        exclusion_filter = RegionFilter(bed_path=test_centromeres_bed)
        empty_filter = RegionFilter()
        
        key1 = generator._make_cache_key(target_filter, exclusion_filter)
        key2 = generator._make_cache_key(target_filter, exclusion_filter)
        key3 = generator._make_cache_key(None, exclusion_filter)
        key4 = generator._make_cache_key(empty_filter, exclusion_filter)
        
        # Same objects should generate same key
        assert key1 == key2
        # Different combinations should generate different keys
        assert key1 != key3
        # None and empty filter should generate same key
        assert key3 == key4
    
    def test_get_valid_regions_no_filters(self, test_bam_path):
        """Test getting valid regions with no filters."""
        generator = RegionAwareWindowGenerator(test_bam_path)
        regions = generator.get_valid_regions(None, None)
        
        # Should return entire genome
        assert len(regions) > 0
        
        # Should include chr1 and chr2
        chromosomes = set(regions.Chromosome)
        assert 'chr1' in chromosomes
        assert 'chr2' in chromosomes
        
        # Total length should be substantial (multiple GBs)
        total_length = (regions.End - regions.Start).sum()
        assert total_length > 1_000_000_000  # > 1GB
    
    def test_get_valid_regions_target_only(self, test_bam_path, test_superdups_bed):
        """Test getting valid regions with target regions only."""
        target_filter = RegionFilter(bed_path=test_superdups_bed)
        generator = RegionAwareWindowGenerator(test_bam_path)
        
        regions = generator.get_valid_regions(target_filter, None)
        
        # Should return only target regions
        assert len(regions) > 0
        
        # All regions should be from chr1 (based on our test file)
        chromosomes = set(regions.Chromosome)
        assert 'chr1' in chromosomes
        
        # Total length should be much smaller than whole genome
        total_length = (regions.End - regions.Start).sum()
        assert total_length < 100_000_000  # < 100MB
    
    def test_get_valid_regions_exclusion_only(self, test_bam_path, test_centromeres_bed):
        """Test getting valid regions with exclusion regions only."""
        exclusion_filter = RegionFilter(bed_path=test_centromeres_bed)
        generator = RegionAwareWindowGenerator(test_bam_path)
        
        regions = generator.get_valid_regions(None, exclusion_filter)
        
        # Should be genome minus exclusions
        assert len(regions) > 0
        
        # Should still have chr1 and chr2, but with gaps
        chromosomes = set(regions.Chromosome)
        assert 'chr1' in chromosomes
        assert 'chr2' in chromosomes
        
        # Total length should be less than whole genome
        whole_genome = generator.get_valid_regions(None, None)
        whole_length = (whole_genome.End - whole_genome.Start).sum()
        excluded_length = (regions.End - regions.Start).sum()
        assert excluded_length < whole_length
    
    def test_get_valid_regions_both_filters(self, test_bam_path, test_superdups_bed, test_centromeres_bed):
        """Test getting valid regions with both target and exclusion filters."""
        target_filter = RegionFilter(bed_path=test_superdups_bed)
        exclusion_filter = RegionFilter(bed_path=test_centromeres_bed)
        generator = RegionAwareWindowGenerator(test_bam_path)
        
        regions = generator.get_valid_regions(target_filter, exclusion_filter)
        
        # Should be target regions minus exclusions
        assert len(regions) >= 0  # Could be empty if all targets are excluded
        
        # Check caching works
        regions2 = generator.get_valid_regions(target_filter, exclusion_filter)
        assert len(generator._valid_regions_cache) == 1
    
    def test_calculate_required_windows(self):
        """Test window calculation."""
        generator = RegionAwareWindowGenerator('test.bam')
        
        assert generator.calculate_required_windows(100, 1) == 100
        assert generator.calculate_required_windows(100, 10) == 10
        assert generator.calculate_required_windows(105, 10) == 11  # ceil
    
    def test_generate_windows_empty_regions(self):
        """Test window generation with empty regions."""
        generator = RegionAwareWindowGenerator('test.bam')
        empty_regions = pr.PyRanges()
        
        windows = generator.generate_windows_for_regions(empty_regions, 10, 1000)
        assert windows == []
    
    def test_generate_windows_basic(self, test_bam_path, test_superdups_bed):
        """Test basic window generation with real regions."""
        generator = RegionAwareWindowGenerator(test_bam_path)
        target_filter = RegionFilter(bed_path=test_superdups_bed)
        regions = target_filter.regions
        
        windows = generator.generate_windows_for_regions(regions, 4, 10000, seed=42)
        
        # Should generate some windows (might be less than 4 if regions are small)
        assert len(windows) > 0
        
        # Check windows are within regions and have correct size
        for chrom, start, end in windows:
            assert chrom.startswith('chr')
            assert start >= 0
            assert end > start
            # Window size should be 10000 or the region size if smaller
            assert end - start <= 10000
    
    def test_generate_windows_small_regions(self):
        """Test window generation with regions smaller than window size."""
        generator = RegionAwareWindowGenerator('test.bam')
        
        # Create small regions
        regions_data = {
            'Chromosome': ['chr1'],
            'Start': [100],
            'End': [200]  # Only 100bp region
        }
        regions = pr.from_dict(regions_data)
        
        windows = generator.generate_windows_for_regions(regions, 2, 1000, seed=42)
        
        assert len(windows) == 2
        # Both windows should be the entire small region
        for chrom, start, end in windows:
            assert chrom == 'chr1'
            assert start == 100
            assert end == 200
    
    def test_generate_windows_proportional_allocation(self, test_bam_path):
        """Test that windows are allocated proportionally to chromosome sizes."""
        generator = RegionAwareWindowGenerator(test_bam_path)
        
        # Get real genome regions for chr1 and chr2
        genome_regions = generator._load_genome_regions()
        chr1_chr2 = genome_regions[genome_regions.Chromosome.isin(['chr1', 'chr2'])]
        
        windows = generator.generate_windows_for_regions(chr1_chr2, 10, 1000000, seed=42)
        
        # Count windows per chromosome
        chr1_windows = sum(1 for chrom, _, _ in windows if chrom == 'chr1')
        chr2_windows = sum(1 for chrom, _, _ in windows if chrom == 'chr2')
        
        # chr1 is larger than chr2, so should get more windows
        assert chr1_windows + chr2_windows == len(windows)
        # chr1 should get at least as many windows as chr2 (usually more)
        assert chr1_windows >= chr2_windows


class TestRegionAwareIntegration:
    """Integration tests for region-aware functionality."""
    
    def test_no_regions_equivalent_to_genome_wide(self, test_bam_path):
        """Test that no regions specified gives genome-wide sampling."""
        generator = RegionAwareWindowGenerator(test_bam_path)
        
        regions = generator.get_valid_regions(None, None)
        windows = generator.generate_windows_for_regions(regions, 6, 1000000, seed=42)
        
        # Should get windows from multiple chromosomes
        chromosomes = set(chrom for chrom, _, _ in windows)
        assert len(chromosomes) > 1
        assert 'chr1' in chromosomes
        
        # Should generate requested number of windows
        assert len(windows) == 6
    
    def test_end_to_end_with_real_data(self, test_bam_path, test_centromeres_bed, test_superdups_bed):
        """Test complete workflow with real BAM and BED files."""
        # Test exclusion only
        exclusion_filter = RegionFilter(bed_path=test_centromeres_bed)
        generator = RegionAwareWindowGenerator(test_bam_path)
        
        excluded_regions = generator.get_valid_regions(None, exclusion_filter)
        windows1 = generator.generate_windows_for_regions(excluded_regions, 5, 100000, seed=42)
        
        assert len(windows1) > 0
        
        # Test target only  
        target_filter = RegionFilter(bed_path=test_superdups_bed)
        target_regions = generator.get_valid_regions(target_filter, None)
        windows2 = generator.generate_windows_for_regions(target_regions, 5, 100000, seed=42)
        
        assert len(windows2) > 0
        
        # Test both together
        both_regions = generator.get_valid_regions(target_filter, exclusion_filter)
        windows3 = generator.generate_windows_for_regions(both_regions, 5, 100000, seed=42)
        
        # Should have fewer or equal windows than target-only
        assert len(windows3) <= len(windows2)