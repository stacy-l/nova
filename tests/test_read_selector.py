"""
Unit tests for read_selector module.
"""

import unittest
import tempfile
import os
import pysam
from unittest.mock import Mock, patch, MagicMock
from pathlib import Path
from typing import Dict, List

from nova.read_selector import ReadSelector, ReadMetadata, LazyReadReference
from nova.region_utils import RegionFilter


class TestReadSelector(unittest.TestCase):
    """Test cases for ReadSelector class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_bam = "test.bam"
        self.selector = ReadSelector(
            bam_path=self.test_bam,
            min_mapq=20,
            max_soft_clip_ratio=0.1,
            min_read_length=10000,
            max_read_length=20000,
            reads_per_window=1  # Default setting - prevents genomic clustering
        )
    
    def test_calculate_gc_content(self):
        """Test GC content calculation."""
        test_cases = [
            ("ATCG", 0.5),
            ("AAAA", 0.0),
            ("GGGG", 1.0),
            ("ATCGATCG", 0.5),
            ("", 0.0),
            ("GCGCGCGC", 1.0)
        ]
        
        for sequence, expected_gc in test_cases:
            with self.subTest(sequence=sequence):
                gc_content = self.selector._calculate_gc_content(sequence)
                self.assertAlmostEqual(gc_content, expected_gc, places=3)
    
    def test_calculate_soft_clip_ratio(self):
        """Test soft clipping ratio calculation."""
        # Mock read with soft clipping
        mock_read = Mock()
        mock_read.query_length = 100
        
        # Test no soft clipping
        mock_read.cigartuples = [(0, 100)]  # 100 match operations
        ratio = self.selector._calculate_soft_clip_ratio(mock_read)
        self.assertEqual(ratio, 0.0)
        
        # Test with soft clipping (10 bases soft clipped out of 100)
        mock_read.cigartuples = [(4, 5), (0, 90), (4, 5)]  # 5+5=10 soft clipped
        ratio = self.selector._calculate_soft_clip_ratio(mock_read)
        self.assertEqual(ratio, 0.1)
        
        # Test with no CIGAR
        mock_read.cigartuples = None
        ratio = self.selector._calculate_soft_clip_ratio(mock_read)
        self.assertEqual(ratio, 0.0)
    
    def test_passes_filters(self):
        """Test read filtering criteria."""
        mock_read = Mock()
        mock_read.is_secondary = False
        mock_read.is_supplementary = False
        mock_read.is_unmapped = False
        mock_read.mapping_quality = 30
        mock_read.query_length = 15000
        mock_read.cigartuples = [(0, 15000)]  # All matches, no soft clipping
        
        # Should pass all filters
        self.assertTrue(self.selector._passes_filters(mock_read))
        
        # Test secondary read rejection
        mock_read.is_secondary = True
        self.assertFalse(self.selector._passes_filters(mock_read))
        mock_read.is_secondary = False
        
        # Test supplementary read rejection
        mock_read.is_supplementary = True
        self.assertFalse(self.selector._passes_filters(mock_read))
        mock_read.is_supplementary = False
        
        # Test unmapped read rejection
        mock_read.is_unmapped = True
        self.assertFalse(self.selector._passes_filters(mock_read))
        mock_read.is_unmapped = False
        
        # Test MAPQ filtering
        mock_read.mapping_quality = 10
        self.assertFalse(self.selector._passes_filters(mock_read))
        mock_read.mapping_quality = 30
        
        # Test read too short
        mock_read.query_length = 5000
        self.assertFalse(self.selector._passes_filters(mock_read))
        
        # Test read too long
        mock_read.query_length = 25000
        self.assertFalse(self.selector._passes_filters(mock_read))
        mock_read.query_length = 15000
        
        # Test too much soft clipping
        mock_read.cigartuples = [(4, 2000), (0, 11000), (4, 2000)]  # 4000/15000 = 0.27 > 0.1
        self.assertFalse(self.selector._passes_filters(mock_read))


class TestReadMetadata(unittest.TestCase):
    """Test ReadMetadata dataclass."""
    
    def test_creation(self):
        """Test ReadMetadata creation and attribute access."""
        metadata = ReadMetadata(
            read_name="test_read",
            original_chr="chr1",
            original_pos=1000,
            read_length=15000,
            mapq=30,
            gc_content=0.45,
            soft_clip_ratio=0.05
        )
        
        self.assertEqual(metadata.read_name, "test_read")
        self.assertEqual(metadata.original_chr, "chr1")
        self.assertEqual(metadata.original_pos, 1000)
        self.assertEqual(metadata.read_length, 15000)
        self.assertEqual(metadata.mapq, 30)
        self.assertEqual(metadata.gc_content, 0.45)
        self.assertEqual(metadata.soft_clip_ratio, 0.05)


class TestReadSelectorRegionInitialization(unittest.TestCase):
    """Test ReadSelector initialization with region filters."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create temporary BED files for testing
        self.target_bed_content = """chr1\t1000\t2000\ttarget1
chr1\t5000\t6000\ttarget2"""
        
        self.exclusion_bed_content = """chr1\t1500\t1600\texclude1
chr2\t3000\t4000\texclude2"""
        
        # Create target BED file
        self.target_bed_file = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
        self.target_bed_file.write(self.target_bed_content)
        self.target_bed_file.flush()
        self.target_bed_path = self.target_bed_file.name
        
        # Create exclusion BED file
        self.exclusion_bed_file = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
        self.exclusion_bed_file.write(self.exclusion_bed_content)
        self.exclusion_bed_file.flush()
        self.exclusion_bed_path = self.exclusion_bed_file.name
        
        # Test data paths
        self.test_bam_path = "tests/test_data/test_reads.bam"
        self.test_centromeres_bed = "tests/test_data/test_centromeres.bed"
        self.test_superdups_bed = "tests/test_data/test_genomic_superdups.bed"
    
    def tearDown(self):
        """Clean up temporary files."""
        os.unlink(self.target_bed_path)
        os.unlink(self.exclusion_bed_path)
    
    def test_init_with_region_filters(self):
        """Test ReadSelector initialization with region filters."""
        target_filter = RegionFilter(bed_path=self.target_bed_path)
        exclusion_filter = RegionFilter(bed_path=self.exclusion_bed_path)
        
        selector = ReadSelector(
            bam_path='dummy.bam',
            target_regions=target_filter,
            exclusion_regions=exclusion_filter
        )
        
        self.assertIs(selector.target_regions, target_filter)
        self.assertIs(selector.exclusion_regions, exclusion_filter)
    
    def test_init_without_region_filters(self):
        """Test ReadSelector initialization without region filters."""
        selector = ReadSelector(bam_path='dummy.bam')
        
        self.assertIsNone(selector.target_regions)
        self.assertIsNone(selector.exclusion_regions)
    
    def test_region_statistics_logging(self):
        """Test that region filter statistics are logged during initialization."""
        target_filter = RegionFilter(bed_path=self.target_bed_path)
        exclusion_filter = RegionFilter(bed_path=self.exclusion_bed_path)
        
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
            self.assertEqual(len(info_calls), 2)
            
            # Check target region logging
            target_call = info_calls[0][0][0]
            self.assertIn('Targeting', target_call)
            self.assertIn('2 regions', target_call)
            self.assertIn('1 chromosomes', target_call)
            
            # Check exclusion region logging
            exclusion_call = info_calls[1][0][0]
            self.assertIn('Excluding', exclusion_call)
            self.assertIn('2 regions', exclusion_call)
            self.assertIn('2 chromosomes', exclusion_call)


class TestReadSelectorFilters(unittest.TestCase):
    """Test ReadSelector filtering functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.selector = ReadSelector(bam_path='dummy.bam')
        
        # Create mock read
        self.mock_read = Mock()
        self.mock_read.is_secondary = False
        self.mock_read.is_supplementary = False
        self.mock_read.is_unmapped = False
        self.mock_read.mapping_quality = 30
        self.mock_read.query_length = 15000
        self.mock_read.cigartuples = [(0, 15000)]  # All matches
        self.mock_read.query_sequence = 'A' * 15000
        self.mock_read.reference_name = 'chr1'
        self.mock_read.reference_start = 1200
        self.mock_read.reference_end = 1800
    
    def test_passes_filters_basic(self):
        """Test that _passes_filters checks basic criteria (no region checking)."""
        # Should pass all basic filters
        self.assertTrue(self.selector._passes_filters(self.mock_read))
        
        # Test MAPQ filter
        self.mock_read.mapping_quality = 10
        self.assertFalse(self.selector._passes_filters(self.mock_read))
        self.mock_read.mapping_quality = 30  # Reset
        
        # Test read length filter
        self.mock_read.query_length = 5000  # Too short
        self.assertFalse(self.selector._passes_filters(self.mock_read))
        self.mock_read.query_length = 15000  # Reset
        
        # Test secondary/supplementary
        self.mock_read.is_secondary = True
        self.assertFalse(self.selector._passes_filters(self.mock_read))
        self.mock_read.is_secondary = False  # Reset


class TestReadSelectorRegionAwareIntegration(unittest.TestCase):
    """Test region-aware window generation integration."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_bam_path = "tests/test_data/test_reads.bam"
        self.test_centromeres_bed = "tests/test_data/test_centromeres.bed"
        self.test_superdups_bed = "tests/test_data/test_genomic_superdups.bed"
        
        # Create temporary target BED file
        target_bed_content = """chr1\t1000\t2000\ttarget1
chr1\t5000\t6000\ttarget2"""
        
        self.target_bed_file = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
        self.target_bed_file.write(target_bed_content)
        self.target_bed_file.flush()
        self.target_bed_path = self.target_bed_file.name
    
    def tearDown(self):
        """Clean up temporary files."""
        os.unlink(self.target_bed_path)
    
    @patch('pysam.AlignmentFile')
    def test_lazy_select_no_regions(self, mock_pysam_file):
        """Test that _lazy_select works without regions (genome-wide)."""
        # Setup mock BAM file
        mock_bam_file = Mock()
        mock_bam_file.references = ['chr1', 'chr2']
        mock_bam_file.lengths = [10000, 8000]
        
        def get_ref_length(chrom):
            length_map = {'chr1': 10000, 'chr2': 8000}
            return length_map.get(chrom)
        mock_bam_file.get_reference_length = get_ref_length
        
        # Mock the pysam.AlignmentFile call in RegionAwareWindowGenerator
        mock_pysam_file.return_value.__enter__.return_value = mock_bam_file
        
        # Create reads to return from fetch
        mock_reads = []
        for i in range(5):
            read = Mock()
            read.is_secondary = False
            read.is_supplementary = False
            read.is_unmapped = False
            read.mapping_quality = 30
            read.query_length = 15000
            read.cigartuples = [(0, 15000)]
            read.query_sequence = 'A' * 15000
            read.reference_name = 'chr1'
            read.reference_start = 1000 + i * 100
            read.query_name = f'read_{i}'
            mock_reads.append(read)
        
        mock_bam_file.fetch.return_value = mock_reads
        
        selector = ReadSelector(bam_path='test.bam', reads_per_window=2)
        selected = selector._lazy_select(mock_bam_file, 3)
        
        self.assertLessEqual(len(selected), 3)
        self.assertTrue(all(hasattr(read, 'read_name') for read in selected))
    
    def test_lazy_select_with_target_regions(self):
        """Test that _lazy_select works with target regions."""
        # Setup mock BAM file
        mock_bam_file = Mock()
        mock_bam_file.references = ['chr1', 'chr2']
        mock_bam_file.lengths = [10000, 8000]
        
        def get_ref_length(chrom):
            length_map = {'chr1': 10000, 'chr2': 8000}
            return length_map.get(chrom)
        mock_bam_file.get_reference_length = get_ref_length
        
        target_filter = RegionFilter(bed_path=self.target_bed_path)
        
        # Create reads that would be in target regions
        mock_reads = []
        for i in range(3):
            read = Mock()
            read.is_secondary = False
            read.is_supplementary = False
            read.is_unmapped = False
            read.mapping_quality = 30
            read.query_length = 15000
            read.cigartuples = [(0, 15000)]
            read.query_sequence = 'A' * 15000
            read.reference_name = 'chr1'
            read.reference_start = 1100 + i * 100  # Within target region
            read.query_name = f'target_read_{i}'
            mock_reads.append(read)
        
        mock_bam_file.fetch.return_value = mock_reads
        
        selector = ReadSelector(
            bam_path='test.bam', 
            target_regions=target_filter,
            reads_per_window=1
        )
        
        selected = selector._lazy_select(mock_bam_file, 2)
        self.assertLessEqual(len(selected), 2)
    
    @patch('pysam.AlignmentFile')
    def test_lazy_select_no_valid_windows(self, mock_pysam_file):
        """Test _lazy_select behavior when no valid windows can be generated."""
        # Setup mock BAM file
        mock_bam_file = Mock()
        mock_bam_file.references = ['chr1', 'chr2']
        mock_bam_file.lengths = [10000, 8000]
        
        def get_ref_length(chrom):
            length_map = {'chr1': 10000, 'chr2': 8000}
            return length_map.get(chrom)
        mock_bam_file.get_reference_length = get_ref_length
        
        # Mock the pysam.AlignmentFile call in RegionAwareWindowGenerator
        mock_pysam_file.return_value.__enter__.return_value = mock_bam_file
        
        # Create an empty target filter
        empty_target = RegionFilter()
        
        selector = ReadSelector(
            bam_path='test.bam',
            target_regions=empty_target
        )
        
        selected = selector._lazy_select(mock_bam_file, 5)
        
        # Should return empty list when no valid regions
        self.assertEqual(selected, [])
    
    def test_select_reads_integration(self):
        """Test select_reads method with region-aware approach."""
        selector = ReadSelector(bam_path='test.bam')
        
        # Mock the _lazy_select method to avoid file I/O
        mock_reads = []
        for i in range(3):
            mock_read = Mock()
            mock_read.original_chr = f'chr{i+1}'
            mock_reads.append(mock_read)
            
        with patch.object(selector, '_lazy_select', return_value=mock_reads):
            with patch('nova.read_selector.pysam.AlignmentFile'):
                selected = selector.select_reads(3)
        
        self.assertEqual(len(selected), 3)
        self.assertEqual(selected, mock_reads)


class TestRegionAwarePerformance(unittest.TestCase):
    """Test performance characteristics of region-aware approach."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create temporary target BED file
        target_bed_content = """chr1\t1000\t2000\ttarget1
chr1\t5000\t6000\ttarget2"""
        
        self.target_bed_file = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False)
        self.target_bed_file.write(target_bed_content)
        self.target_bed_file.flush()
        self.target_bed_path = self.target_bed_file.name
    
    def tearDown(self):
        """Clean up temporary files."""
        os.unlink(self.target_bed_path)
    
    def test_no_per_read_region_checking(self):
        """Test that region checking is not called per-read anymore."""
        target_filter = RegionFilter(bed_path=self.target_bed_path)
        selector = ReadSelector(bam_path='test.bam', target_regions=target_filter)
        
        # Verify the old per-read region checking method doesn't exist
        self.assertFalse(hasattr(selector, '_passes_region_filters'))
        
        # Verify _passes_filters doesn't do region checking
        mock_read = Mock()
        mock_read.is_secondary = False
        mock_read.is_supplementary = False  
        mock_read.is_unmapped = False
        mock_read.mapping_quality = 30
        mock_read.query_length = 15000
        mock_read.cigartuples = [(0, 15000)]
        
        # Should pass regardless of coordinates (no region checking)
        mock_read.reference_name = 'chr999'  # Non-existent chromosome
        mock_read.reference_start = 999999   # Outside any regions
        mock_read.reference_end = 999999 + 1000
        
        # This should still pass basic filters
        self.assertTrue(selector._passes_filters(mock_read))


if __name__ == '__main__':
    unittest.main()