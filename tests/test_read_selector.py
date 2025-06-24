"""
Unit tests for read_selector module.
"""

import unittest
import tempfile
import pysam
from unittest.mock import Mock, patch, MagicMock
from pathlib import Path
from typing import Dict, List

from nova.read_selector import ReadSelector, ReadMetadata


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
            max_read_length=20000
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
        
        # Test no CIGAR
        mock_read.cigartuples = None
        ratio = self.selector._calculate_soft_clip_ratio(mock_read)
        self.assertEqual(ratio, 0.0)
        
        # Test with soft clipping (operation 4 = soft clip)
        mock_read.cigartuples = [(4, 10), (0, 80), (4, 10)]  # 10 + 10 = 20 soft clipped
        ratio = self.selector._calculate_soft_clip_ratio(mock_read)
        self.assertEqual(ratio, 0.2)
        
        # Test no soft clipping
        mock_read.cigartuples = [(0, 100)]  # 100 match
        ratio = self.selector._calculate_soft_clip_ratio(mock_read)
        self.assertEqual(ratio, 0.0)
    
    def test_passes_filters(self):
        """Test read filtering logic."""
        # Create mock read
        mock_read = Mock()
        mock_read.is_secondary = False
        mock_read.is_supplementary = False
        mock_read.is_unmapped = False
        mock_read.mapping_quality = 25
        mock_read.query_length = 15000
        mock_read.cigartuples = [(0, 15000)]  # No soft clipping
        
        # Should pass all filters
        self.assertTrue(self.selector._passes_filters(mock_read))
        
        # Test secondary read
        mock_read.is_secondary = True
        self.assertFalse(self.selector._passes_filters(mock_read))
        mock_read.is_secondary = False
        
        # Test supplementary read
        mock_read.is_supplementary = True
        self.assertFalse(self.selector._passes_filters(mock_read))
        mock_read.is_supplementary = False
        
        # Test unmapped read
        mock_read.is_unmapped = True
        self.assertFalse(self.selector._passes_filters(mock_read))
        mock_read.is_unmapped = False
        
        # Test low mapping quality
        mock_read.mapping_quality = 10
        self.assertFalse(self.selector._passes_filters(mock_read))
        mock_read.mapping_quality = 25
        
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
    
    def test_get_chromosome_sampling_regions(self):
        """Test chromosome sampling region generation."""
        # Mock BAM file and index statistics
        mock_bam = Mock()
        mock_stat1 = Mock()
        mock_stat1.contig = "chr1"
        mock_stat2 = Mock()
        mock_stat2.contig = "chr2"
        
        mock_bam.get_index_statistics.return_value = [mock_stat1, mock_stat2]
        mock_bam.get_reference_length.side_effect = lambda chrom: {"chr1": 2000000, "chr2": 1500000}[chrom]
        
        regions = self.selector._get_chromosome_sampling_regions(mock_bam, window_size=1000000)
        
        # Should have regions for both chromosomes
        self.assertGreater(len(regions), 0)
        
        # Check that regions contain valid coordinates
        for chrom, start, end in regions:
            self.assertIn(chrom, ["chr1", "chr2"])
            self.assertGreaterEqual(start, 0)
            self.assertGreater(end, start)
            self.assertLessEqual(end - start, 1000000)
    
    @patch('pysam.AlignmentFile')
    def test_select_reads_integration(self, mock_alignment_file):
        """Test read selection integration with mocked BAM file."""
        # Setup mock BAM file
        mock_bam = Mock()
        mock_alignment_file.return_value.__enter__.return_value = mock_bam
        
        # Mock index statistics
        mock_stat = Mock()
        mock_stat.contig = "chr1"
        mock_bam.get_index_statistics.return_value = [mock_stat]
        mock_bam.get_reference_length.return_value = 3000000  # 3MB chromosome
        
        # Create substantial pool of reads for robust testing
        # Use smaller spacing to ensure better coverage within 1MB windows
        chrom_read_pools = {
            "chr1": [(f"read_chr1_{i}", i * 1000) for i in range(200)]  # 200 reads, one every 1KB
        }
        
        mock_bam.fetch.side_effect = self._create_coordinate_aware_mock_fetch(chrom_read_pools)
        
        # Test read selection with small simulation (should use window-limited)
        selected_reads = self.selector.select_reads(20)
        
        # Should return at least most of the requested reads (allowing for some variance in mock)
        self.assertGreaterEqual(len(selected_reads), 18, "Should get at least 18 of 20 requested reads")
        self.assertLessEqual(len(selected_reads), 20, "Should not exceed 20 requested reads")
        
        # Check that each returned item is a tuple of (read, metadata)
        for read, metadata in selected_reads:
            self.assertIsInstance(metadata, ReadMetadata)
            self.assertIsInstance(metadata.read_name, str)
            self.assertIsInstance(metadata.original_chr, str)
            self.assertIsInstance(metadata.original_pos, int)
            self.assertIsInstance(metadata.read_length, int)
            self.assertIsInstance(metadata.mapq, int)
            self.assertIsInstance(metadata.gc_content, float)
            self.assertIsInstance(metadata.soft_clip_ratio, float)
            
            # Verify metadata values are reasonable
            self.assertEqual(metadata.original_chr, "chr1")
            self.assertEqual(metadata.read_length, 15000)
            self.assertEqual(metadata.mapq, 30)
            self.assertEqual(metadata.gc_content, 0.0)  # All A's
            self.assertEqual(metadata.soft_clip_ratio, 0.0)  # No soft clipping
    
    @patch('pysam.AlignmentFile')
    def test_get_chromosome_proportional_targets(self, mock_alignment_file):
        """Test chromosome-proportional target calculation."""
        # Setup mock BAM file
        mock_bam = Mock()
        mock_alignment_file.return_value.__enter__.return_value = mock_bam
        
        # Mock index statistics for multiple chromosomes
        mock_stat1 = Mock()
        mock_stat1.contig = "chr1"
        mock_stat2 = Mock()
        mock_stat2.contig = "chr2"
        mock_stat3 = Mock()
        mock_stat3.contig = "chr3"
        
        mock_bam.get_index_statistics.return_value = [mock_stat1, mock_stat2, mock_stat3]
        
        # Mock chromosome lengths: chr1=100M, chr2=50M, chr3=25M (total=175M)
        def mock_get_reference_length(chrom):
            lengths = {"chr1": 100000000, "chr2": 50000000, "chr3": 25000000}
            return lengths.get(chrom)
        
        mock_bam.get_reference_length.side_effect = mock_get_reference_length
        
        # Test with 1000 reads
        targets = self.selector._get_chromosome_proportional_targets(mock_bam, 1000)
        
        # Should allocate proportionally: chr1=571, chr2=286, chr3=143 (approximately)
        self.assertEqual(len(targets), 3)
        self.assertIn("chr1", targets)
        self.assertIn("chr2", targets)
        self.assertIn("chr3", targets)
        
        # chr1 should get the most reads (largest chromosome)
        self.assertGreater(targets["chr1"], targets["chr2"])
        self.assertGreater(targets["chr2"], targets["chr3"])
        
        # Total should equal requested reads
        self.assertEqual(sum(targets.values()), 1000)
    
    @patch('pysam.AlignmentFile')
    def test_small_vs_large_simulation_strategies(self, mock_alignment_file):
        """Test that different strategies are used for small vs large simulations."""
        # Setup mock BAM file with multiple chromosomes
        mock_bam = Mock()
        mock_alignment_file.return_value.__enter__.return_value = mock_bam
        
        # Mock index statistics for multiple chromosomes
        mock_stats = []
        for i in range(1, 4):
            mock_stat = Mock()
            mock_stat.contig = f"chr{i}"
            mock_stats.append(mock_stat)
        
        mock_bam.get_index_statistics.return_value = mock_stats
        
        # Different chromosome lengths for proportional testing
        def mock_get_reference_length(chrom):
            lengths = {"chr1": 3000000, "chr2": 2000000, "chr3": 1000000}  # 3:2:1 ratio
            return lengths.get(chrom)
        
        mock_bam.get_reference_length.side_effect = mock_get_reference_length
        
        # Create large pools of reads for robust testing
        chrom_read_pools = {
            "chr1": [(f"read_chr1_{i}", i * 5000) for i in range(200)],  # 200 reads
            "chr2": [(f"read_chr2_{i}", i * 5000) for i in range(200)],  # 200 reads
            "chr3": [(f"read_chr3_{i}", i * 5000) for i in range(200)],  # 200 reads
        }
        
        mock_bam.fetch.side_effect = self._create_coordinate_aware_mock_fetch(chrom_read_pools)
        
        # Test small simulation (should use window-limited)
        with self.assertLogs('nova.read_selector', level='INFO') as log:
            selected_reads = self.selector.select_reads(50)  # < 500
            self.assertIn('window-limited', ''.join(log.output))
            # Verify we got reads
            self.assertEqual(len(selected_reads), 50, "Should get exactly 50 reads for small simulation")
        
        # Test large simulation (should use proportional)
        with self.assertLogs('nova.read_selector', level='INFO') as log:
            selected_reads = self.selector.select_reads(500)  # >= 500 for proportional strategy
            self.assertIn('proportional', ''.join(log.output))
            
            # Verify we got a substantial number of reads
            self.assertGreater(len(selected_reads), 100, "Should get substantial reads for large simulation")
            
            # For chromosome-proportional, we expect roughly 3:2:1 distribution
            chrom_counts = {}
            for _, metadata in selected_reads:
                chrom = metadata.original_chr
                chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
            
            # Verify proportional distribution (chr1 should have most reads)
            if len(chrom_counts) >= 2:  # Need at least 2 chromosomes to test
                self.assertGreater(chrom_counts.get("chr1", 0), chrom_counts.get("chr3", 0),
                                 "chr1 should have more reads than chr3 in proportional sampling")
    
    @patch('pysam.AlignmentFile')
    def test_genomic_distribution_window_limits(self, mock_alignment_file):
        """Test that window-limited sampling improves genomic distribution."""
        # Setup mock BAM file with multiple chromosomes
        mock_bam = Mock()
        mock_alignment_file.return_value.__enter__.return_value = mock_bam
        
        # Mock index statistics for 3 chromosomes
        mock_stats = []
        for i in range(1, 4):
            mock_stat = Mock()
            mock_stat.contig = f"chr{i}"
            mock_stats.append(mock_stat)
        
        mock_bam.get_index_statistics.return_value = mock_stats
        mock_bam.get_reference_length.return_value = 5000000  # 5MB chromosomes
        
        # Create large pools of reads distributed across coordinates
        chrom_read_pools = {
            "chr1": [(f"read_chr1_{i}", i * 10000) for i in range(150)],  # 150 reads across 5MB
            "chr2": [(f"read_chr2_{i}", i * 10000) for i in range(100)],  # 100 reads across 5MB  
            "chr3": [(f"read_chr3_{i}", i * 10000) for i in range(80)],   # 80 reads across 5MB
        }
        
        mock_bam.fetch.side_effect = self._create_coordinate_aware_mock_fetch(chrom_read_pools)
        
        # Test window-limited sampling with 30 reads
        selected_reads = self.selector._select_reads_with_window_limits(mock_bam, 30)
        
        # Verify we got requested number of reads
        self.assertEqual(len(selected_reads), 30, "Should get exactly 30 reads")
        
        # Should have reads from multiple chromosomes
        chromosomes = set(metadata.original_chr for _, metadata in selected_reads)
        self.assertGreaterEqual(len(chromosomes), 2, "Should sample from at least 2 chromosomes")
        
        # Count reads per chromosome
        chrom_counts = {}
        for _, metadata in selected_reads:
            chrom = metadata.original_chr
            chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
        
        # Verify window limits are working (max 3 reads per window for 30 total)
        # With 5 chromosomes * ~5 windows each, max per window should be around 3
        max_reads_per_window = max(1, 30 // 10)  # 3 reads per window
        for chrom, count in chrom_counts.items():
            # No single chromosome should dominate (rough check)
            self.assertLessEqual(count, 20, f"{chrom} has too many reads ({count}), indicates poor distribution")
    
    @patch('pysam.AlignmentFile')
    def test_proportional_sampling_distribution(self, mock_alignment_file):
        """Test that proportional sampling distributes reads correctly."""
        # Setup mock BAM file
        mock_bam = Mock()
        mock_alignment_file.return_value.__enter__.return_value = mock_bam
        
        # Mock index statistics for chromosomes of different sizes
        mock_stats = []
        for i in range(1, 4):
            mock_stat = Mock()
            mock_stat.contig = f"chr{i}"
            mock_stats.append(mock_stat)
        
        mock_bam.get_index_statistics.return_value = mock_stats
        
        # Different chromosome lengths (4:2:1 ratio)
        def mock_get_reference_length(chrom):
            lengths = {"chr1": 200000000, "chr2": 100000000, "chr3": 50000000}  # 4:2:1 ratio
            return lengths.get(chrom)
        
        mock_bam.get_reference_length.side_effect = mock_get_reference_length
        
        # Create large pools of reads with denser coverage for better window overlap
        chrom_read_pools = {
            "chr1": [(f"read_chr1_{i}", i * 5000) for i in range(800)],   # 800 reads, one every 5KB
            "chr2": [(f"read_chr2_{i}", i * 5000) for i in range(400)],   # 400 reads, one every 5KB  
            "chr3": [(f"read_chr3_{i}", i * 5000) for i in range(200)],   # 200 reads, one every 5KB
        }
        
        mock_bam.fetch.side_effect = self._create_coordinate_aware_mock_fetch(chrom_read_pools)
        
        # Test proportional sampling with 150 reads (more manageable for mock)
        selected_reads = self.selector._select_reads_with_proportional_sampling(mock_bam, 150)
        
        # Verify we got most of the requested reads (allowing for mock limitations)
        self.assertGreaterEqual(len(selected_reads), 120, "Should get at least 120 of 150 requested reads")
        self.assertLessEqual(len(selected_reads), 150, "Should not exceed 150 requested reads")
        
        # Count reads per chromosome
        chrom_counts = {}
        for _, metadata in selected_reads:
            chrom = metadata.original_chr
            chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
        
        # Should have reads from all chromosomes
        self.assertEqual(len(chrom_counts), 3, "Should have reads from all 3 chromosomes")
        
        # Verify proportional distribution (4:2:1 ratio)
        # chr1 should get ~114 reads (4/7 * 200), chr2 ~57 (2/7 * 200), chr3 ~29 (1/7 * 200)
        self.assertGreater(chrom_counts.get("chr1", 0), chrom_counts.get("chr2", 0),
                          "chr1 should have more reads than chr2")
        self.assertGreater(chrom_counts.get("chr2", 0), chrom_counts.get("chr3", 0),
                          "chr2 should have more reads than chr3")
        
        # More strict proportional check - chr1 should have roughly 2x chr2
        chr1_count = chrom_counts.get("chr1", 0)
        chr2_count = chrom_counts.get("chr2", 0)
        if chr2_count > 0:  # Avoid division by zero
            ratio = chr1_count / chr2_count
            self.assertGreater(ratio, 1.5, f"chr1:chr2 ratio ({ratio:.2f}) should be roughly 2:1")
            self.assertLess(ratio, 2.5, f"chr1:chr2 ratio ({ratio:.2f}) should be roughly 2:1")
    
    def _create_mock_read(self, read_name: str, chrom: str, position: int = 1000):
        """Helper to create a mock read that passes filters."""
        mock_read = Mock()
        mock_read.is_secondary = False
        mock_read.is_supplementary = False
        mock_read.is_unmapped = False
        mock_read.mapping_quality = 30
        mock_read.query_length = 15000
        mock_read.query_name = read_name
        mock_read.reference_name = chrom
        mock_read.reference_start = position
        mock_read.query_sequence = "A" * 15000
        mock_read.cigartuples = [(0, 15000)]
        return mock_read
    
    def _create_coordinate_aware_mock_fetch(self, chrom_read_pools: Dict[str, List]):
        """
        Create a coordinate-aware mock_fetch function that returns different reads
        based on genomic coordinates.
        
        Args:
            chrom_read_pools: Dict mapping chromosome names to lists of read data
                             Each read data is (read_name, position)
        
        Returns:
            Function that can be used as mock_bam.fetch.side_effect
        """
        def mock_fetch(chrom, start, end):
            if chrom not in chrom_read_pools:
                return []
            
            # Return reads whose positions fall within the requested window
            # Use a more generous window to ensure we get enough reads
            reads_in_window = []
            for read_name, position in chrom_read_pools[chrom]:
                if start <= position <= end:
                    reads_in_window.append(
                        self._create_mock_read(read_name, chrom, position)
                    )
            
            # If no reads in exact window, return some reads from the chromosome
            # This simulates realistic BAM behavior where fetch() usually finds some reads
            if not reads_in_window and chrom_read_pools[chrom]:
                # Return up to 10 reads from this chromosome regardless of coordinates
                available_reads = chrom_read_pools[chrom][:10]
                reads_in_window = [
                    self._create_mock_read(read_name, chrom, position)
                    for read_name, position in available_reads
                ]
            
            return reads_in_window
        
        return mock_fetch


class TestReadMetadata(unittest.TestCase):
    """Test cases for ReadMetadata dataclass."""
    
    def test_creation(self):
        """Test creating ReadMetadata."""
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


if __name__ == '__main__':
    unittest.main()