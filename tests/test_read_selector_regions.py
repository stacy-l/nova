"""
Tests for ReadSelector region filtering functionality.
"""

import pytest
import tempfile
import os
from unittest.mock import Mock, patch
import pysam
from nova.read_selector import ReadSelector
from nova.region_utils import RegionFilter


@pytest.fixture
def sample_target_bed():
    """Create a temporary BED file for target regions."""
    content = """chr1\t1000\t2000\ttarget1
chr1\t5000\t6000\ttarget2"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write(content)
        f.flush()
        yield f.name
    os.unlink(f.name)


@pytest.fixture
def sample_exclusion_bed():
    """Create a temporary BED file for exclusion regions."""
    content = """chr1\t1500\t1600\texclude1
chr2\t3000\t4000\texclude2"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write(content)
        f.flush()
        yield f.name
    os.unlink(f.name)


@pytest.fixture
def mock_read():
    """Create a mock aligned read for testing."""
    read = Mock(spec=pysam.AlignedSegment)
    read.is_secondary = False
    read.is_supplementary = False
    read.is_unmapped = False
    read.mapping_quality = 30
    read.query_length = 15000
    read.cigartuples = [(0, 15000)]  # All matches
    read.query_sequence = 'A' * 15000
    read.reference_name = 'chr1'
    read.reference_start = 1200
    read.reference_end = 1800
    return read


class TestReadSelectorRegionFiltering:
    """Test ReadSelector region filtering functionality."""
    
    def test_init_with_region_filters(self, sample_target_bed, sample_exclusion_bed):
        """Test ReadSelector initialization with region filters."""
        target_filter = RegionFilter(bed_path=sample_target_bed)
        exclusion_filter = RegionFilter(bed_path=sample_exclusion_bed)
        
        selector = ReadSelector(
            bam_path='dummy.bam',
            target_regions=target_filter,
            exclusion_regions=exclusion_filter
        )
        
        assert selector.target_regions is target_filter
        assert selector.exclusion_regions is exclusion_filter
    
    def test_init_without_region_filters(self):
        """Test ReadSelector initialization without region filters."""
        selector = ReadSelector(bam_path='dummy.bam')
        
        assert selector.target_regions is None
        assert selector.exclusion_regions is None
    
    def test_passes_region_filters_no_filters(self, mock_read):
        """Test _passes_region_filters with no filters set."""
        selector = ReadSelector(bam_path='dummy.bam')
        assert selector._passes_region_filters(mock_read)
    
    def test_passes_region_filters_target_only(self, sample_target_bed, mock_read):
        """Test _passes_region_filters with target regions only."""
        target_filter = RegionFilter(bed_path=sample_target_bed)
        selector = ReadSelector(
            bam_path='dummy.bam',
            target_regions=target_filter
        )
        
        # Mock read overlaps target region (chr1:1200-1800 overlaps chr1:1000-2000)
        assert selector._passes_region_filters(mock_read)
        
        # Mock read outside target region
        mock_read.reference_start = 3000
        mock_read.reference_end = 3500
        assert not selector._passes_region_filters(mock_read)
    
    def test_passes_region_filters_exclusion_only(self, sample_exclusion_bed, mock_read):
        """Test _passes_region_filters with exclusion regions only."""
        exclusion_filter = RegionFilter(bed_path=sample_exclusion_bed)
        selector = ReadSelector(
            bam_path='dummy.bam',
            exclusion_regions=exclusion_filter
        )
        
        # Mock read overlaps exclusion region (chr1:1200-1800 overlaps chr1:1500-1600)
        assert not selector._passes_region_filters(mock_read)
        
        # Mock read outside exclusion region
        mock_read.reference_start = 2000
        mock_read.reference_end = 2500
        assert selector._passes_region_filters(mock_read)
    
    def test_passes_region_filters_both(self, sample_target_bed, sample_exclusion_bed, mock_read):
        """Test _passes_region_filters with both target and exclusion regions."""
        target_filter = RegionFilter(bed_path=sample_target_bed)
        exclusion_filter = RegionFilter(bed_path=sample_exclusion_bed)
        selector = ReadSelector(
            bam_path='dummy.bam',
            target_regions=target_filter,
            exclusion_regions=exclusion_filter
        )
        
        # Mock read overlaps target but also exclusion (should be excluded)
        mock_read.reference_start = 1200
        mock_read.reference_end = 1800
        assert not selector._passes_region_filters(mock_read)
        
        # Mock read overlaps target but not exclusion (should pass)
        mock_read.reference_start = 1000
        mock_read.reference_end = 1200
        assert selector._passes_region_filters(mock_read)
        
        # Mock read doesn't overlap target (should fail regardless of exclusion)
        mock_read.reference_start = 3000
        mock_read.reference_end = 3500
        assert not selector._passes_region_filters(mock_read)
    
    def test_passes_region_filters_invalid_coordinates(self, sample_target_bed):
        """Test _passes_region_filters with invalid read coordinates."""
        target_filter = RegionFilter(bed_path=sample_target_bed)
        selector = ReadSelector(
            bam_path='dummy.bam',
            target_regions=target_filter
        )
        
        read = Mock(spec=pysam.AlignedSegment)
        read.reference_name = None
        read.reference_start = None
        read.reference_end = None
        
        assert not selector._passes_region_filters(read)
    
    def test_passes_filters_integration(self, sample_target_bed, mock_read):
        """Test that _passes_filters calls _passes_region_filters."""
        target_filter = RegionFilter(bed_path=sample_target_bed)
        selector = ReadSelector(
            bam_path='dummy.bam',
            target_regions=target_filter
        )
        
        # Read should pass all filters including region filters
        assert selector._passes_filters(mock_read)
        
        # Move read outside target region
        mock_read.reference_start = 3000
        mock_read.reference_end = 3500
        
        # Should fail region filter even though other filters pass
        assert not selector._passes_filters(mock_read)


class TestReadSelectorRegionStatistics:
    """Test ReadSelector region filter statistics and logging."""
    
    def test_region_statistics_logging(self, sample_target_bed, sample_exclusion_bed):
        """Test that region filter statistics are logged during initialization."""
        target_filter = RegionFilter(bed_path=sample_target_bed)
        exclusion_filter = RegionFilter(bed_path=sample_exclusion_bed)
        
        with patch('nova.read_selector.logging.getLogger') as mock_logger:
            mock_log = Mock()
            mock_logger.return_value = mock_log
            
            ReadSelector(
                bam_path='dummy.bam',
                target_regions=target_filter,
                exclusion_regions=exclusion_filter
            )
            
            # Check that info messages were logged
            info_calls = mock_log.info.call_args_list
            assert len(info_calls) == 2
            
            # Check target region logging
            target_call = info_calls[0][0][0]
            assert 'Targeting' in target_call
            assert '2 regions' in target_call
            assert '1 chromosomes' in target_call
            
            # Check exclusion region logging
            exclusion_call = info_calls[1][0][0]
            assert 'Excluding' in exclusion_call
            assert '2 regions' in exclusion_call
            assert '2 chromosomes' in exclusion_call