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
        config_file = 'tests/test_data/test_config.json'
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        # Generate variants
        insertion_ids = generator.generate_from_config(config)
        self.assertGreater(len(insertion_ids), 0)
        
        # Test read selection
        selector = ReadSelector(bam_file, min_mapq=20, max_soft_clip_ratio=0.1)
        selected_reads = selector.select_reads(min(len(insertion_ids), 10))
        self.assertGreater(len(selected_reads), 0)
        
        # Test insertion
        inserter = ReadInserter(registry, random_seed=42)
        insertion_records, modified_sequences, skip_stats = inserter.insert_random_mode(
            selected_reads[:len(insertion_ids)], insertion_ids
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