"""
Tests for memory optimization features: lazy read processing and streaming insertion.
Uses integration test scaffolding with lazy/streaming methods substituted.
"""

import unittest
import tempfile
import json
import os
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from nova.variant_registry import VariantRegistry
from nova.variant_generator import VariantGenerator
from nova.read_inserter import ReadInserter
from nova.read_selector import ReadSelector, LazyReadReference, ReadMetadata


class TestLazyReadProcessing(unittest.TestCase):
    """Test lazy read processing functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_bam_path = "test_file.bam"
        self.selector = ReadSelector(
            self.test_bam_path, 
            min_mapq=20, 
            max_soft_clip_ratio=0.1
        )
    
    def test_lazy_read_reference_creation(self):
        """Test LazyReadReference object creation and metadata."""
        lazy_ref = LazyReadReference(
            bam_path="test.bam",
            read_name="test_read_1",
            original_chr="chr1",
            original_pos=1000,
            read_length=15000,
            mapq=30,
            gc_content=0.45,
            soft_clip_ratio=0.02
        )
        
        # Test basic attributes
        self.assertEqual(lazy_ref.read_name, "test_read_1")
        self.assertEqual(lazy_ref.original_chr, "chr1")
        self.assertEqual(lazy_ref.original_pos, 1000)
        self.assertEqual(lazy_ref.read_length, 15000)
        self.assertEqual(lazy_ref.mapq, 30)
        self.assertEqual(lazy_ref.gc_content, 0.45)
        self.assertEqual(lazy_ref.soft_clip_ratio, 0.02)
        
        # Test metadata conversion
        metadata = lazy_ref.get_metadata()
        self.assertIsInstance(metadata, ReadMetadata)
        self.assertEqual(metadata.read_name, "test_read_1")
        self.assertEqual(metadata.original_chr, "chr1")
        self.assertEqual(metadata.original_pos, 1000)
        self.assertEqual(metadata.read_length, 15000)
        self.assertEqual(metadata.mapq, 30)
        self.assertEqual(metadata.gc_content, 0.45)
        self.assertEqual(metadata.soft_clip_ratio, 0.02)
    
    @patch('nova.read_selector.pysam.AlignmentFile')
    def test_lazy_sequence_fetching(self, mock_alignment_file):
        """Test on-demand sequence fetching from LazyReadReference."""
        # Mock BAM file and read
        mock_bam = Mock()
        mock_read = Mock()
        mock_read.query_name = "test_read_1"
        mock_read.reference_start = 1000
        mock_read.query_sequence = "ATCGATCGATCGATCG"
        
        mock_bam.fetch.return_value = [mock_read]
        mock_alignment_file.return_value.__enter__.return_value = mock_bam
        
        lazy_ref = LazyReadReference(
            bam_path="test.bam",
            read_name="test_read_1",
            original_chr="chr1",
            original_pos=1000,
            read_length=16,
            mapq=30,
            gc_content=0.5,
            soft_clip_ratio=0.0
        )
        
        # Test sequence fetching
        sequence = lazy_ref.get_sequence()
        self.assertEqual(sequence, "ATCGATCGATCGATCG")
        
        # Verify correct BAM fetch call
        mock_bam.fetch.assert_called_once_with("chr1", 1000, 1001)


class TestMemoryOptimizedIntegration(unittest.TestCase):
    """Memory-optimized integration tests based on test_integration.py scaffold."""
    
    def _create_mock_lazy_reads(self, count, bam_path="test.bam"):
        """Create mock LazyReadReference objects for testing."""
        lazy_reads = []
        chromosomes = ["chr1", "chr2", "chr3"]
        
        for i in range(count):
            chrom = chromosomes[i % len(chromosomes)]
            lazy_ref = LazyReadReference(
                bam_path=bam_path,
                read_name=f"read_{i}",
                original_chr=chrom,
                original_pos=i * 50000,
                read_length=15000,
                mapq=30,
                gc_content=0.45,
                soft_clip_ratio=0.02
            )
            lazy_reads.append(lazy_ref)
        
        return lazy_reads
    
    @patch.object(LazyReadReference, 'get_sequence')
    def test_memory_optimized_end_to_end_simulation(self, mock_get_sequence):
        """Test complete simulation workflow using memory-optimized methods."""
        # Mock sequence fetching (mirroring test_end_to_end_simulation)
        mock_get_sequence.return_value = "A" * 2000  # 2kb reads
        
        # Create mock data (same as original test)
        registry = VariantRegistry()
        generator = VariantGenerator(registry, random_seed=42)
        
        # Create test FASTA file (same as original test)
        sequences = [
            SeqRecord(Seq("ATCGATCGATCGATCG"), id="TestAlu", description=""),
        ]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            temp_fasta = f.name
        
        try:
            SeqIO.write(sequences, temp_fasta, "fasta")
            
            # Generate variants from config (same as original test)
            config = {
                'random': {'n': 2, 'length': 10},
                'simple': {'n': 1, 'repeat': 'CAG', 'units': 5},
                'predefined': {
                    'Test': {
                        'fasta': temp_fasta,
                        'spec': {'TestAlu': 1}
                    }
                }
            }
            
            insertion_ids = generator.generate_from_config(config)
            
            # Verify generation (same as original test)
            self.assertEqual(len(insertion_ids), 4)  # 2 + 1 + 1
            self.assertEqual(len(registry), 4)
            
            # Check registry statistics (same as original test)
            stats = registry.get_statistics()
            self.assertEqual(stats['total_sequences'], 4)
            self.assertEqual(stats['type_counts']['random'], 2)
            self.assertEqual(stats['type_counts']['simple'], 1)
            self.assertEqual(stats['type_counts']['TestAlu'], 1)
            
            # Create mock lazy reads instead of mock pysam reads
            lazy_reads = self._create_mock_lazy_reads(4)
            
            # Update read lengths to match expected 2000 length
            for lazy_ref in lazy_reads:
                lazy_ref.read_length = 2000
            
            # Test streaming insertion (replacing insert_random_mode)
            inserter = ReadInserter(registry, min_distance_from_ends=500, random_seed=42)
            
            # Use streaming mode instead of batch mode
            results = list(inserter.insert_streaming_mode(lazy_reads, insertion_ids))
            insertion_records = [result[0] for result in results]
            modified_sequences = [result[1] for result in results]
            
            # Verify insertions (same checks as original test)
            self.assertEqual(len(insertion_records), 4)
            self.assertEqual(len(modified_sequences), 4)
            
            # Check insertion records (same checks as original test)
            for record in insertion_records:
                self.assertIsInstance(record.base_read_name, str)
                self.assertIsInstance(record.modified_read_name, str)
                self.assertIsInstance(record.insertion_id, str)
                self.assertGreater(record.insertion_pos, 500)
                self.assertLess(record.insertion_pos, 1500)  # 2000 - 500
                self.assertIn(record.original_chr, ["chr1", "chr2", "chr3"])
            
            # Check modified sequences (same checks as original test)
            for seq_record in modified_sequences:
                self.assertGreater(len(seq_record.seq), 2000)  # Original + insertion
            
            # Test statistics (same as original test)
            insertion_stats = inserter.get_insertion_statistics(insertion_records)
            self.assertEqual(insertion_stats['total_insertions'], 4)
            
            # Test saving/loading using streaming methods
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                
                # Save registry (same as original test)
                registry_file = temp_path / "registry.json"
                registry.save_to_json(str(registry_file))
                self.assertTrue(registry_file.exists())
                
                # Save using streaming methods
                output_prefix = temp_path / "streaming_output"
                stats = inserter.save_streaming_results(lazy_reads, insertion_ids, str(output_prefix))
                
                # Verify files exist (adapting original test checks)
                self.assertTrue(Path(stats['records_file']).exists())
                self.assertTrue(Path(stats['sequences_file']).exists())
                
                # Verify file contents (same checks as original test)
                with open(registry_file) as f:
                    registry_data = json.load(f)
                    self.assertEqual(len(registry_data['sequences']), 4)
                
                with open(stats['records_file']) as f:
                    records_data = json.load(f)
                    self.assertEqual(len(records_data), 4)
                
                fasta_records = list(SeqIO.parse(stats['sequences_file'], "fasta"))
                self.assertEqual(len(fasta_records), 4)
        
        finally:
            Path(temp_fasta).unlink()
    
    @unittest.skipUnless(
        os.path.exists('tests/test_data/test_reads.bam'), 
        "Real BAM file not available - skipping memory-optimized small-scale BAM test"
    )
    def test_memory_optimized_small_scale_with_real_bam(self):
        """Test memory-optimized small-scale pipeline with real BAM file."""
        bam_file = 'tests/test_data/test_reads.bam'
        
        # Setup registry and generator (same as integration test)
        registry = VariantRegistry()
        generator = VariantGenerator(registry, random_seed=42)
        
        # Use the test config and sequences (same as integration test)
        config_file = 'tests/test_data/test_config_small.json'
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        # Generate variants (same as integration test)
        insertion_ids = generator.generate_from_config(config)
        self.assertGreater(len(insertion_ids), 0)
        
        # Test LAZY read selection instead of regular selection
        selector = ReadSelector(bam_file, min_mapq=20, max_soft_clip_ratio=0.1)
        lazy_reads = selector.select_lazy_reads(len(insertion_ids))
        self.assertGreater(len(lazy_reads), 0)
        
        # Test STREAMING insertion instead of batch insertion
        inserter = ReadInserter(registry, random_seed=42)
        num_to_insert = min(len(lazy_reads), len(insertion_ids))
        
        results = list(inserter.insert_streaming_mode(
            lazy_reads[:num_to_insert], insertion_ids[:num_to_insert]
        ))
        
        # Verify results (same checks as integration test)
        insertion_records = [result[0] for result in results]
        modified_sequences = [result[1] for result in results]
        
        self.assertGreater(len(insertion_records), 0)
        self.assertGreater(len(modified_sequences), 0)
        
        # Test saving results using streaming methods
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Save using both old method (for registry) and new streaming methods
            registry_file = temp_path / "registry.json"
            registry.save_to_json(str(registry_file))
            
            output_prefix = temp_path / "streaming_output"
            stats = inserter.save_streaming_results(
                lazy_reads[:num_to_insert], insertion_ids[:num_to_insert], str(output_prefix)
            )
            
            # Verify files exist and have content (same checks as integration test)
            self.assertTrue(registry_file.exists())
            self.assertTrue(Path(stats['records_file']).exists())
            self.assertTrue(Path(stats['sequences_file']).exists())
            
            # Quick content verification (same as integration test)
            with open(registry_file) as f:
                registry_data = json.load(f)
                self.assertGreater(len(registry_data['sequences']), 0)
    
    @unittest.skipUnless(
        os.path.exists('tests/test_data/test_reads.bam'), 
        "Real BAM file not available - skipping memory-optimized medium-scale BAM test"
    )
    def test_memory_optimized_medium_scale_with_real_bam(self):
        """Test memory-optimized medium-scale pipeline (100 reads) with real BAM data."""
        bam_file = 'tests/test_data/test_reads.bam'
        
        # Setup registry and generator (same as integration test)
        registry = VariantRegistry()
        generator = VariantGenerator(registry, random_seed=42)
        
        # Use the updated test config with 100 total variants (same as integration test)
        config_file = 'tests/test_data/test_config_medium.json'
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        # Generate variants (same as integration test)
        insertion_ids = generator.generate_from_config(config)
        expected_count = 100  # 20 random + 10 simple + 70 predefined
        self.assertEqual(len(insertion_ids), expected_count, 
                        f"Expected {expected_count} variants, got {len(insertion_ids)}")
        
        # Test LAZY read selection with 100 reads (should trigger window-limited sampling)
        selector = ReadSelector(bam_file, min_mapq=20, max_soft_clip_ratio=0.1)
        
        # Capture logs to verify window-limited sampling strategy is used (same as integration test)
        with self.assertLogs('nova.read_selector', level='INFO') as log:
            lazy_reads = selector.select_lazy_reads(100)  # Using lazy version!
            
            # Verify window-limited sampling strategy was used (same check as integration test)
            log_output = ''.join(log.output)
            self.assertIn('window-limited', log_output.lower(), 
                         "Expected window-limited sampling strategy for 100 reads")
        
        self.assertLessEqual(len(lazy_reads), 100, 
                            "Should not exceed 100 requested reads")
        
        # Test STREAMING insertion with matching numbers
        inserter = ReadInserter(registry, random_seed=42)
        num_to_insert = min(len(lazy_reads), len(insertion_ids))
        
        results = list(inserter.insert_streaming_mode(
            lazy_reads[:num_to_insert], insertion_ids[:num_to_insert]
        ))
        
        insertion_records = [result[0] for result in results]
        modified_sequences = [result[1] for result in results]
        
        # Verify medium-scale results (same checks as integration test)
        min_successful_expected = int(100 * 0.95)  # 95
        self.assertGreaterEqual(len(insertion_records), min_successful_expected, 
                               f"Should have at least {min_successful_expected} insertion records (95%+ of 100)")
        self.assertGreaterEqual(len(modified_sequences), min_successful_expected, 
                               f"Should have at least {min_successful_expected} modified sequences (95%+ of 100)")
        
        # Verify chromosome distribution in selected reads (same check as integration test)
        chrom_counts = {}
        for lazy_ref in lazy_reads:
            chrom = lazy_ref.original_chr
            chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
        
        # Should have reads from multiple chromosomes due to window-limited sampling
        self.assertGreaterEqual(len(chrom_counts), 2, 
                               "Should sample from at least 2 chromosomes")
        
        # Test saving results at medium scale using streaming methods
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            registry_file = temp_path / "registry.json"
            registry.save_to_json(str(registry_file))
            
            output_prefix = temp_path / "streaming_output"
            stats = inserter.save_streaming_results(
                lazy_reads[:num_to_insert], insertion_ids[:num_to_insert], str(output_prefix)
            )
            
            # Verify files exist and have substantial content (same checks as integration test)
            self.assertTrue(registry_file.exists())
            self.assertTrue(Path(stats['records_file']).exists())
            self.assertTrue(Path(stats['sequences_file']).exists())
            
            # Verify content scale (same checks as integration test)
            with open(registry_file) as f:
                registry_data = json.load(f)
                self.assertEqual(len(registry_data['sequences']), 100, 
                               "Registry should contain 100 sequences")
            
            with open(stats['records_file']) as f:
                records_data = json.load(f)
                self.assertGreaterEqual(len(records_data), min_successful_expected, 
                                      f"Records should contain at least {min_successful_expected} insertions (95% success)")
            
            fasta_records = list(SeqIO.parse(stats['sequences_file'], "fasta"))
            self.assertGreaterEqual(len(fasta_records), min_successful_expected, 
                                  f"FASTA should contain at least {min_successful_expected} sequences (95% success)")
    
    @unittest.skipUnless(
        os.path.exists('tests/test_data/test_reads.bam'), 
        "Real BAM file not available - skipping memory-optimized large-scale BAM test"
    )
    def test_memory_optimized_large_scale_with_real_bam(self):
        """Test memory-optimized large-scale pipeline (500 reads) with real BAM data."""
        bam_file = 'tests/test_data/test_reads.bam'
        
        # Setup registry and generator (same as integration test)
        registry = VariantRegistry()
        generator = VariantGenerator(registry, random_seed=42)
        
        # Use the large config with 500 total variants (same as integration test)
        config_file = 'tests/test_data/test_config_large.json'
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        # Generate variants (same as integration test)
        insertion_ids = generator.generate_from_config(config)
        expected_count = 500  # 100 random + 50 simple + 350 predefined
        # Allow for >95% success rate due to real BAM file limitations (same as integration test)
        min_expected = int(expected_count * 0.95)  # 475
        self.assertGreaterEqual(len(insertion_ids), min_expected, 
                        f"Expected at least {min_expected} variants (95% of {expected_count}), got {len(insertion_ids)}")
        
        actual_count = len(insertion_ids)
        
        # Test LAZY read selection with 500 reads (should trigger proportional sampling)
        selector = ReadSelector(bam_file, min_mapq=20, max_soft_clip_ratio=0.1)
        
        # Capture logs to verify proportional sampling strategy is used (same as integration test)
        with self.assertLogs('nova.read_selector', level='INFO') as log:
            lazy_reads = selector.select_lazy_reads(500)  # Using lazy version!
            
            # Verify proportional sampling strategy was used (same check as integration test)
            log_output = ''.join(log.output)
            self.assertIn('proportional', log_output.lower(), 
                         "Expected proportional sampling strategy for 500 reads")
        
        # Should get close to 500 reads (same check as integration test)
        self.assertLessEqual(len(lazy_reads), 500, 
                            "Should not exceed 500 requested reads")
        
        # Test STREAMING insertion with matching numbers
        inserter = ReadInserter(registry, random_seed=42)
        num_to_insert = min(len(lazy_reads), len(insertion_ids))
        
        # Calculate expected minimum successful insertions (same as integration test)
        min_insertions_expected = int(num_to_insert * 0.95)
        
        results = list(inserter.insert_streaming_mode(
            lazy_reads[:num_to_insert], insertion_ids[:num_to_insert]
        ))
        
        insertion_records = [result[0] for result in results]
        modified_sequences = [result[1] for result in results]
        
        # Verify large-scale results (same checks as integration test)
        self.assertGreaterEqual(len(insertion_records), min_insertions_expected, 
                               f"Should have at least {min_insertions_expected} insertion records")
        self.assertGreaterEqual(len(modified_sequences), min_insertions_expected, 
                               f"Should have at least {min_insertions_expected} modified sequences")
        
        # Verify chromosome distribution shows proportional sampling (same check as integration test)
        chrom_counts = {}
        for lazy_ref in lazy_reads:
            chrom = lazy_ref.original_chr
            chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
        
        # Should have reads from multiple chromosomes (same check as integration test)
        self.assertGreaterEqual(len(chrom_counts), 3, 
                               "Should sample from at least 3 chromosomes")
        
        # Test saving results at large scale using streaming methods
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            registry_file = temp_path / "registry.json"
            registry.save_to_json(str(registry_file))
            
            output_prefix = temp_path / "streaming_output"
            stats = inserter.save_streaming_results(
                lazy_reads[:num_to_insert], insertion_ids[:num_to_insert], str(output_prefix)
            )
            
            # Verify files exist and have substantial content (same checks as integration test)
            self.assertTrue(registry_file.exists())
            self.assertTrue(Path(stats['records_file']).exists())
            self.assertTrue(Path(stats['sequences_file']).exists())
            
            # Verify content scale (same checks as integration test)
            with open(registry_file) as f:
                registry_data = json.load(f)
                saved_count = len(registry_data['sequences'])
                
                self.assertLessEqual(saved_count, actual_count, 
                                    f"Registry should contain at most {actual_count} sequences")
                # Verify we got >95% of target (475+ variants)
                self.assertGreaterEqual(saved_count, int(500 * 0.95), 
                                      "Should generate >95% of target variants")
            
            fasta_records = list(SeqIO.parse(stats['sequences_file'], "fasta"))
            self.assertGreaterEqual(len(fasta_records), min_insertions_expected, 
                                  f"FASTA should contain at least {min_insertions_expected} sequences (95% success)")


if __name__ == '__main__':
    unittest.main()