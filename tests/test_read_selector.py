"""
Unit tests for read_selector module.
"""

import unittest
import tempfile
import pysam
from unittest.mock import Mock, patch, MagicMock
from pathlib import Path

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
        mock_bam.get_reference_length.return_value = 1000000
        
        # Create mock reads
        mock_reads = []
        for i in range(10):
            mock_read = Mock()
            mock_read.is_secondary = False
            mock_read.is_supplementary = False
            mock_read.is_unmapped = False
            mock_read.mapping_quality = 30
            mock_read.query_length = 15000
            mock_read.query_name = f"read_{i}"
            mock_read.reference_name = "chr1"
            mock_read.reference_start = i * 1000
            mock_read.query_sequence = "A" * 15000
            mock_read.cigartuples = [(0, 15000)]
            mock_reads.append(mock_read)
        
        # Mock fetch to return our reads
        mock_bam.fetch.return_value = iter(mock_reads)
        
        # Test read selection
        selected_reads = self.selector.select_reads(5)
        
        # Should return requested number of reads (or available reads)
        self.assertLessEqual(len(selected_reads), 5)
        self.assertGreater(len(selected_reads), 0)
        
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
        # Setup mock BAM file
        mock_bam = Mock()
        mock_alignment_file.return_value.__enter__.return_value = mock_bam
        
        # Mock index statistics
        mock_stat = Mock()
        mock_stat.contig = "chr1"
        mock_bam.get_index_statistics.return_value = [mock_stat]
        mock_bam.get_reference_length.return_value = 1000000
        
        # Create mock reads
        mock_reads = []
        for i in range(20):
            mock_read = Mock()
            mock_read.is_secondary = False
            mock_read.is_supplementary = False
            mock_read.is_unmapped = False
            mock_read.mapping_quality = 30
            mock_read.query_length = 15000
            mock_read.query_name = f"read_{i}"
            mock_read.reference_name = "chr1"
            mock_read.reference_start = i * 1000
            mock_read.query_sequence = "A" * 15000
            mock_read.cigartuples = [(0, 15000)]
            mock_reads.append(mock_read)
        
        mock_bam.fetch.return_value = iter(mock_reads)
        
        # Test small simulation (should use window-limited)
        with self.assertLogs('nova.read_selector', level='INFO') as log:
            selected_reads = self.selector.select_reads(50)  # < 500
            self.assertIn('window-limited', ''.join(log.output))
        
        # Test large simulation (should use proportional)
        with self.assertLogs('nova.read_selector', level='INFO') as log:
            selected_reads = self.selector.select_reads(500)  # >= 500
            self.assertIn('proportional', ''.join(log.output))
    
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
        mock_bam.get_reference_length.return_value = 1000000
        
        # Create reads distributed across chromosomes
        def mock_fetch(chrom, start, end):
            # Return different numbers of reads per chromosome/region
            if chrom == "chr1":
                return [self._create_mock_read(f"read_{chrom}_{i}", chrom) for i in range(5)]
            elif chrom == "chr2":
                return [self._create_mock_read(f"read_{chrom}_{i}", chrom) for i in range(3)]
            else:  # chr3
                return [self._create_mock_read(f"read_{chrom}_{i}", chrom) for i in range(2)]
        
        mock_bam.fetch.side_effect = mock_fetch
        
        # Test window-limited sampling
        selected_reads = self.selector._select_reads_with_window_limits(mock_bam, 20)
        
        # Should have reads from multiple chromosomes
        chromosomes = set(metadata.original_chr for _, metadata in selected_reads)
        self.assertGreater(len(chromosomes), 1, "Should sample from multiple chromosomes")
    
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
        
        # Different chromosome lengths
        def mock_get_reference_length(chrom):
            lengths = {"chr1": 200000000, "chr2": 100000000, "chr3": 50000000}  # 4:2:1 ratio
            return lengths.get(chrom)
        
        mock_bam.get_reference_length.side_effect = mock_get_reference_length
        
        # Create reads for each chromosome
        def mock_fetch(chrom, start, end):
            # Return plenty of reads for each chromosome
            return [self._create_mock_read(f"read_{chrom}_{i}", chrom) for i in range(20)]
        
        mock_bam.fetch.side_effect = mock_fetch
        
        # Test proportional sampling with 700 reads
        selected_reads = self.selector._select_reads_with_proportional_sampling(mock_bam, 700)
        
        # Count reads per chromosome
        chrom_counts = {}
        for _, metadata in selected_reads:
            chrom = metadata.original_chr
            chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
        
        # Should have reads from all chromosomes
        self.assertEqual(len(chrom_counts), 3)
        
        # chr1 should have most reads (largest), chr3 should have least
        self.assertGreater(chrom_counts.get("chr1", 0), chrom_counts.get("chr2", 0))
        self.assertGreater(chrom_counts.get("chr2", 0), chrom_counts.get("chr3", 0))
    
    def _create_mock_read(self, read_name: str, chrom: str):
        """Helper to create a mock read that passes filters."""
        mock_read = Mock()
        mock_read.is_secondary = False
        mock_read.is_supplementary = False
        mock_read.is_unmapped = False
        mock_read.mapping_quality = 30
        mock_read.query_length = 15000
        mock_read.query_name = read_name
        mock_read.reference_name = chrom
        mock_read.reference_start = 1000
        mock_read.query_sequence = "A" * 15000
        mock_read.cigartuples = [(0, 15000)]
        return mock_read


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