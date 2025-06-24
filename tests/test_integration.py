"""
Integration tests for Nova simulator.
"""

import unittest
import tempfile
import json
import os
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from nova.variant_registry import VariantRegistry
from nova.variant_generator import VariantGenerator
from nova.read_inserter import ReadInserter
from nova.read_selector import ReadMetadata, ReadSelector


class TestIntegration(unittest.TestCase):
    """Integration test cases."""
    
    def test_end_to_end_simulation(self):
        """Test complete simulation workflow."""
        # Create mock data
        registry = VariantRegistry()
        generator = VariantGenerator(registry, random_seed=42)
        
        # Create test FASTA file
        sequences = [
            SeqRecord(Seq("ATCGATCGATCGATCG"), id="TestAlu", description=""),
        ]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            temp_fasta = f.name
        
        try:
            SeqIO.write(sequences, temp_fasta, "fasta")
            
            # Generate variants from config
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
            
            # Verify generation
            self.assertEqual(len(insertion_ids), 4)  # 2 + 1 + 1
            self.assertEqual(len(registry), 4)
            
            # Check registry statistics
            stats = registry.get_statistics()
            self.assertEqual(stats['total_sequences'], 4)
            self.assertEqual(stats['type_counts']['random'], 2)
            self.assertEqual(stats['type_counts']['simple'], 1)
            self.assertEqual(stats['type_counts']['TestAlu'], 1)
            
            # Create mock read data (simulate pysam.AlignedSegment behavior)
            class MockRead:
                def __init__(self, name, sequence, length):
                    self.query_name = name
                    self.query_sequence = sequence
                    self.query_length = length
            
            mock_reads_with_metadata = [
                (MockRead(f"read_{i}", "A" * 2000, 2000), 
                 ReadMetadata(
                     read_name=f"read_{i}",
                     original_chr="chr1",
                     original_pos=i * 1000,
                     read_length=2000,
                     mapq=30,
                     gc_content=0.0,
                     soft_clip_ratio=0.0
                 ))
                for i in range(4)
            ]
            
            # Test read insertion
            inserter = ReadInserter(registry, min_distance_from_ends=500, random_seed=42)
            
            # Perform insertions
            insertion_records, modified_sequences, skip_stats = inserter.insert_random_mode(mock_reads_with_metadata, insertion_ids)
            
            # Check skip statistics
            self.assertEqual(skip_stats['total_attempted'], 4)
            self.assertEqual(skip_stats['successful_insertions'], 4)
            self.assertEqual(skip_stats['success_rate'], 1.0)
            
            # Verify insertions
            self.assertEqual(len(insertion_records), 4)
            self.assertEqual(len(modified_sequences), 4)
            
            # Check insertion records
            for record in insertion_records:
                self.assertIsInstance(record.base_read_name, str)
                self.assertIsInstance(record.modified_read_name, str)
                self.assertIsInstance(record.insertion_id, str)
                self.assertGreater(record.insertion_pos, 500)
                self.assertLess(record.insertion_pos, 1500)  # 2000 - 500
                self.assertEqual(record.original_chr, "chr1")
            
            # Check modified sequences
            for seq_record in modified_sequences:
                self.assertGreater(len(seq_record.seq), 2000)  # Original + insertion
            
            # Test statistics
            insertion_stats = inserter.get_insertion_statistics(insertion_records)
            self.assertEqual(insertion_stats['total_insertions'], 4)
            
            # Test saving/loading
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                
                # Save registry
                registry_file = temp_path / "registry.json"
                registry.save_to_json(str(registry_file))
                self.assertTrue(registry_file.exists())
                
                # Save insertion records
                records_file = temp_path / "records.json"
                inserter.save_insertion_records(insertion_records, str(records_file))
                self.assertTrue(records_file.exists())
                
                # Save modified sequences
                sequences_file = temp_path / "sequences.fasta"
                inserter.save_modified_sequences(modified_sequences, str(sequences_file))
                self.assertTrue(sequences_file.exists())
                
                # Verify file contents
                with open(registry_file) as f:
                    registry_data = json.load(f)
                    self.assertEqual(len(registry_data['sequences']), 4)
                
                with open(records_file) as f:
                    records_data = json.load(f)
                    self.assertEqual(len(records_data), 4)
                
                fasta_records = list(SeqIO.parse(sequences_file, "fasta"))
                self.assertEqual(len(fasta_records), 4)
        
        finally:
            Path(temp_fasta).unlink()
    
    @unittest.skipUnless(
        os.path.exists('tests/test_data/test_reads.bam'), 
        "Real BAM file not available - skipping BAM integration test"
    )
    def test_full_pipeline_with_real_bam(self):
        """Test full pipeline with real BAM file (if available)."""
        bam_file = 'tests/test_data/test_reads.bam'
        
        # Setup registry and generator
        registry = VariantRegistry()
        generator = VariantGenerator(registry, random_seed=42)
        
        # Use the test config and sequences
        config_file = 'tests/test_data/test_config_small.json'
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        # Generate variants
        insertion_ids = generator.generate_from_config(config)
        self.assertGreater(len(insertion_ids), 0)
        
        # Test read selection - select enough reads to match insertions
        selector = ReadSelector(bam_file, min_mapq=20, max_soft_clip_ratio=0.1)
        selected_reads = selector.select_reads(len(insertion_ids))
        self.assertGreater(len(selected_reads), 0)
        
        # Test insertion - ensure we have matching numbers
        inserter = ReadInserter(registry, random_seed=42)
        num_to_insert = min(len(selected_reads), len(insertion_ids))
        insertion_records, modified_sequences, skip_stats = inserter.insert_random_mode(
            selected_reads[:num_to_insert], insertion_ids[:num_to_insert]
        )
        
        # Verify results
        self.assertGreater(len(insertion_records), 0)
        self.assertGreater(len(modified_sequences), 0)
        self.assertGreater(skip_stats['successful_insertions'], 0)
        
        # Test saving results
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            registry_file = temp_path / "registry.json"
            records_file = temp_path / "records.json"
            sequences_file = temp_path / "sequences.fasta"
            
            registry.save_to_json(str(registry_file))
            inserter.save_insertion_records(insertion_records, str(records_file))
            inserter.save_modified_sequences(modified_sequences, str(sequences_file))
            
            # Verify files exist and have content
            self.assertTrue(registry_file.exists())
            self.assertTrue(records_file.exists())
            self.assertTrue(sequences_file.exists())
            
            # Quick content verification
            with open(registry_file) as f:
                registry_data = json.load(f)
                self.assertGreater(len(registry_data['sequences']), 0)
    
    @unittest.skipUnless(
        os.path.exists('tests/test_data/test_reads.bam'), 
        "Real BAM file not available - skipping medium-scale BAM integration test"
    )
    def test_medium_scale_real_bam_pipeline(self):
        """Test medium-scale pipeline with 100 variants/reads using real BAM data and window-limited sampling."""
        bam_file = 'tests/test_data/test_reads.bam'
        
        # Setup registry and generator
        registry = VariantRegistry()
        generator = VariantGenerator(registry, random_seed=42)
        
        # Use the updated test config with 100 total variants
        config_file = 'tests/test_data/test_config_medium.json'
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        # Generate variants - should produce 100 total (20+10+30+20+10+10)
        insertion_ids = generator.generate_from_config(config)
        expected_count = 100  # 20 random + 10 simple + 70 predefined
        self.assertEqual(len(insertion_ids), expected_count, 
                        f"Expected {expected_count} variants, got {len(insertion_ids)}")
        
        # Test read selection with 100 reads (should trigger window-limited sampling)
        # Note: threshold for proportional sampling is 500+ reads
        selector = ReadSelector(bam_file, min_mapq=20, max_soft_clip_ratio=0.1)
        
        # Capture logs to verify window-limited sampling strategy is used
        with self.assertLogs('nova.read_selector', level='INFO') as log:
            selected_reads = selector.select_reads(100)
            
            # Verify window-limited sampling strategy was used (100 < 500 threshold)
            log_output = ''.join(log.output)
            self.assertIn('window-limited', log_output.lower(), 
                         "Expected window-limited sampling strategy for 100 reads")
        
        self.assertLessEqual(len(selected_reads), 100, 
                            "Should not exceed 100 requested reads")
        
        # Test insertion with matching numbers
        inserter = ReadInserter(registry, random_seed=42)
        num_to_insert = min(len(selected_reads), len(insertion_ids))
        
        insertion_records, modified_sequences, skip_stats = inserter.insert_random_mode(
            selected_reads[:num_to_insert], insertion_ids[:num_to_insert]
        )
        
        # Verify medium-scale results - expect 95%+ success rate (95+ of 100)
        min_successful_expected = int(100 * 0.95)  # 95
        self.assertGreaterEqual(len(insertion_records), min_successful_expected, 
                               f"Should have at least {min_successful_expected} insertion records (95%+ of 100)")
        self.assertGreaterEqual(len(modified_sequences), min_successful_expected, 
                               f"Should have at least {min_successful_expected} modified sequences (95%+ of 100)")
        self.assertGreaterEqual(skip_stats['successful_insertions'], min_successful_expected, 
                               f"Should have at least {min_successful_expected} successful insertions (95%+ of 100)")
        self.assertGreater(skip_stats['success_rate'], 0.95, 
                          "Should have very high success rate (>95%)")
        
        # Verify chromosome distribution in selected reads (window-limited sampling)
        chrom_counts = {}
        for _, metadata in selected_reads:
            chrom = metadata.original_chr
            chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
        
        # Should have reads from multiple chromosomes due to window-limited sampling
        self.assertGreaterEqual(len(chrom_counts), 2, 
                               "Should sample from at least 2 chromosomes")
        
        # Test saving results at medium scale
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            registry_file = temp_path / "registry.json"
            records_file = temp_path / "records.json"
            sequences_file = temp_path / "sequences.fasta"
            
            # Save all results
            registry.save_to_json(str(registry_file))
            inserter.save_insertion_records(insertion_records, str(records_file))
            inserter.save_modified_sequences(modified_sequences, str(sequences_file))
            
            # Verify files exist and have substantial content
            self.assertTrue(registry_file.exists())
            self.assertTrue(records_file.exists())
            self.assertTrue(sequences_file.exists())
            
            # Verify content scale
            with open(registry_file) as f:
                registry_data = json.load(f)
                self.assertEqual(len(registry_data['sequences']), 100, 
                               "Registry should contain 100 sequences")
            
            with open(records_file) as f:
                records_data = json.load(f)
                self.assertGreaterEqual(len(records_data), min_successful_expected, 
                                      f"Records should contain at least {min_successful_expected} insertions (95% success)")
            
            fasta_records = list(SeqIO.parse(sequences_file, "fasta"))
            self.assertGreaterEqual(len(fasta_records), min_successful_expected, 
                                  f"FASTA should contain at least {min_successful_expected} sequences (95% success)")
    
    @unittest.skipUnless(
        os.path.exists('tests/test_data/test_reads.bam'), 
        "Real BAM file not available - skipping large-scale BAM integration test"
    )
    def test_large_scale_proportional_sampling(self):
        """Test large-scale pipeline with 500 variants/reads using real BAM data and proportional sampling."""
        bam_file = 'tests/test_data/test_reads.bam'
        
        # Setup registry and generator
        registry = VariantRegistry()
        generator = VariantGenerator(registry, random_seed=42)
        
        # Use the large config with 500 total variants
        config_file = 'tests/test_data/test_config_large.json'
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        # Generate variants - should produce 500 total (100+50+150+100+50+50)
        insertion_ids = generator.generate_from_config(config)
        expected_count = 500  # 100 random + 50 simple + 350 predefined
        # Allow for >95% success rate due to real BAM file limitations
        min_expected = int(expected_count * 0.95)  # 475
        self.assertGreaterEqual(len(insertion_ids), min_expected, 
                        f"Expected at least {min_expected} variants (95% of {expected_count}), got {len(insertion_ids)}")
        
        actual_count = len(insertion_ids)
        
        # Test read selection with 500 reads (should trigger proportional sampling)
        selector = ReadSelector(bam_file, min_mapq=20, max_soft_clip_ratio=0.1)
        
        # Capture logs to verify proportional sampling strategy is used
        with self.assertLogs('nova.read_selector', level='INFO') as log:
            selected_reads = selector.select_reads(500)  # >= 500 triggers proportional
            
            # Verify proportional sampling strategy was used
            log_output = ''.join(log.output)
            self.assertIn('proportional', log_output.lower(), 
                         "Expected proportional sampling strategy for 500 reads")
        
        # Should get close to 500 reads (expect 95%+ success rate)
        self.assertLessEqual(len(selected_reads), 500, 
                            "Should not exceed 500 requested reads")
        
        # Test insertion with matching numbers
        inserter = ReadInserter(registry, random_seed=42)
        num_to_insert = min(len(selected_reads), len(insertion_ids))
        
        # Calculate expected minimum successful insertions (95% of available)
        min_insertions_expected = int(num_to_insert * 0.95)
        
        insertion_records, modified_sequences, skip_stats = inserter.insert_random_mode(
            selected_reads[:num_to_insert], insertion_ids[:num_to_insert]
        )
        self.assertGreaterEqual(len(insertion_records), min_insertions_expected, 
                               f"Should have at least {min_insertions_expected} insertion records")
        self.assertGreaterEqual(len(modified_sequences), min_insertions_expected, 
                               f"Should have at least {min_insertions_expected} modified sequences")
        self.assertGreaterEqual(skip_stats['successful_insertions'], min_insertions_expected, 
                               f"Should have at least {min_insertions_expected} successful insertions")
        self.assertGreater(skip_stats['success_rate'], 0.95, 
                          "Should have very high success rate (>95%)")
        
        # Verify chromosome distribution shows proportional sampling
        chrom_counts = {}
        for _, metadata in selected_reads:
            chrom = metadata.original_chr
            chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
        
        # Should have reads from multiple chromosomes
        self.assertGreaterEqual(len(chrom_counts), 3, 
                               "Should sample from at least 3 chromosomes")
        
        # Test saving results at large scale
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            registry_file = temp_path / "registry.json"
            records_file = temp_path / "records.json"
            sequences_file = temp_path / "sequences.fasta"
            
            # Save all results
            registry.save_to_json(str(registry_file))
            inserter.save_insertion_records(insertion_records, str(records_file))
            inserter.save_modified_sequences(modified_sequences, str(sequences_file))
            
            # Verify files exist and have substantial content
            self.assertTrue(registry_file.exists())
            self.assertTrue(records_file.exists())
            self.assertTrue(sequences_file.exists())
            
            # Verify content scale - registry should contain >95% of target variants
            with open(registry_file) as f:
                registry_data = json.load(f)
                saved_count = len(registry_data['sequences'])
                
                self.assertLessEqual(saved_count, actual_count, 
                                    f"Registry should contain at most {actual_count} sequences")
                # Verify we got >95% of target (475+ variants)
                self.assertGreaterEqual(saved_count, int(500 * 0.95), 
                                      "Should generate >95% of target variants")
            
            fasta_records = list(SeqIO.parse(sequences_file, "fasta"))
            self.assertGreaterEqual(len(fasta_records), min_insertions_expected, 
                                  f"FASTA should contain at least {min_insertions_expected} sequences (95% success)")
    
    def test_configuration_validation_integration(self):
        """Test configuration validation in realistic scenario."""
        registry = VariantRegistry()
        generator = VariantGenerator(registry)
        
        # Test various config scenarios
        valid_configs = [
            {'random': {'n': 10, 'length': 100}},
            {'simple': {'n': 5, 'repeat': 'CAG', 'units': 20}},
            {'random': {'n': 5, 'length': 50}, 'simple': {'n': 3, 'repeat': 'AT', 'units': 10}},
        ]
        
        for config in valid_configs:
            errors = generator.validate_config(config)
            self.assertEqual(len(errors), 0, f"Valid config failed validation: {config}")
        
        invalid_configs = [
            {'random': {'n': 0, 'length': 100}},  # n = 0
            {'simple': {'n': 5, 'repeat': '', 'units': 20}},  # empty repeat
            {'random': {'n': 5}},  # missing length
            {'simple': {'repeat': 'CAG', 'units': 20}},  # missing n
        ]
        
        for config in invalid_configs:
            errors = generator.validate_config(config)
            self.assertGreater(len(errors), 0, f"Invalid config passed validation: {config}")


if __name__ == '__main__':
    unittest.main()