"""
Tests for verifying output equivalence between traditional and memory-optimized methods.
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
from nova.read_selector import ReadSelector
from nova.read_inserter import ReadInserter
from nova.sequence_utils import extract_sequence_content_hashes, compare_fasta_content


class TestOutputEquivalence(unittest.TestCase):
    """Test equivalence between traditional and memory-optimized simulation outputs."""
    
    @unittest.skipUnless(
        os.path.exists('tests/test_data/test_reads.bam'), 
        "Real BAM file not available - skipping output equivalence test"
    )
    def test_small_scale_output_equivalence(self):
        """Test that traditional and memory-optimized methods produce identical sequence content."""
        bam_file = 'tests/test_data/test_reads.bam'
        config_file = 'tests/test_data/test_config_small.json'
        
        # Load configuration
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # === Run Traditional Method ===
            registry_trad = VariantRegistry()
            generator_trad = VariantGenerator(registry_trad, random_seed=42)
            insertion_ids_trad = generator_trad.generate_from_config(config)
            
            selector_trad = ReadSelector(bam_file, min_mapq=20, max_soft_clip_ratio=0.1)
            selected_reads_trad = selector_trad.select_reads(len(insertion_ids_trad))
            
            inserter_trad = ReadInserter(registry_trad, random_seed=42)
            num_to_insert = min(len(selected_reads_trad), len(insertion_ids_trad))
            
            insertion_records_trad, modified_sequences_trad, skip_stats_trad = inserter_trad.insert_random_mode(
                selected_reads_trad[:num_to_insert], 
                insertion_ids_trad[:num_to_insert]
            )
            
            # Save traditional results
            trad_fasta = temp_path / "traditional_sequences.fasta"
            inserter_trad.save_modified_sequences(modified_sequences_trad, str(trad_fasta))
            
            # === Run Memory-Optimized Method ===
            registry_mem = VariantRegistry()
            generator_mem = VariantGenerator(registry_mem, random_seed=42)
            insertion_ids_mem = generator_mem.generate_from_config(config)
            
            selector_mem = ReadSelector(bam_file, min_mapq=20, max_soft_clip_ratio=0.1)
            lazy_reads_mem = selector_mem.select_lazy_reads(len(insertion_ids_mem))
            
            inserter_mem = ReadInserter(registry_mem, random_seed=42)
            num_to_insert_mem = min(len(lazy_reads_mem), len(insertion_ids_mem))
            
            # Use streaming mode and collect results
            results_mem = list(inserter_mem.insert_streaming_mode(
                lazy_reads_mem[:num_to_insert_mem], 
                insertion_ids_mem[:num_to_insert_mem]
            ))
            
            insertion_records_mem = [result[0] for result in results_mem]
            modified_sequences_mem = [result[1] for result in results_mem]
            
            # Save memory-optimized results
            mem_fasta = temp_path / "memory_optimized_sequences.fasta"
            inserter_mem.save_modified_sequences(modified_sequences_mem, str(mem_fasta))
            
            # === Compare Results ===
            
            # Compare sequence counts
            self.assertEqual(len(modified_sequences_trad), len(modified_sequences_mem),
                           "Traditional and memory-optimized methods produced different numbers of sequences")
            
            # Compare sequence content
            comparison_info = compare_fasta_content(str(trad_fasta), str(mem_fasta))
            
            # Assert sequence content is identical
            self.assertTrue(comparison_info['sequences_identical'], 
                          f"Sequence content differs between methods. Comparison: {comparison_info}")
            
            # Verify no sequences are unique to either method
            self.assertEqual(comparison_info['unique_to_file1'], 0,
                           "Traditional method produced sequences not found in memory-optimized output")
            self.assertEqual(comparison_info['unique_to_file2'], 0,
                           "Memory-optimized method produced sequences not found in traditional output")
            
            # Log successful comparison
            print(f"✓ Output equivalence verified: {comparison_info['common_sequences']} identical sequences")
    
    def test_sequence_hash_function(self):
        """Test the sequence hashing function with known inputs."""
        # Create test FASTA file
        test_sequences = [
            SeqRecord(Seq("ATCGATCGATCG"), id="seq1", description=""),
            SeqRecord(Seq("GGCCTTAAGGCC"), id="seq2", description=""),
            SeqRecord(Seq("ATCGATCGATCG"), id="seq3", description=""),  # Duplicate sequence
        ]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            temp_fasta = f.name
        
        try:
            SeqIO.write(test_sequences, temp_fasta, "fasta")
            
            # Extract hashes
            hashes = extract_sequence_content_hashes(temp_fasta)
            
            # Should have 2 unique sequences (seq1 and seq3 are identical)
            self.assertEqual(len(hashes), 2, "Should detect 2 unique sequences from 3 records")
            
            # Test comparison with itself
            info = compare_fasta_content(temp_fasta, temp_fasta)
            self.assertTrue(info['sequences_identical'], "File should be identical to itself")
            self.assertEqual(info['file1_sequences'], 2, "Should count 2 unique sequences")
            self.assertEqual(info['file2_sequences'], 2, "Should count 2 unique sequences")
            self.assertEqual(info['common_sequences'], 2, "All sequences should be common")
            
        finally:
            Path(temp_fasta).unlink()
    
    @unittest.skipUnless(
        os.path.exists('tests/test_data/test_reads.bam'), 
        "Real BAM file not available - skipping medium scale equivalence test"
    )
    def test_medium_scale_output_equivalence(self):
        """Test output equivalence at medium scale (100 variants)."""
        bam_file = 'tests/test_data/test_reads.bam'
        config_file = 'tests/test_data/test_config_medium.json'
        
        # Load configuration
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # === Run Both Methods ===
            
            # Traditional
            registry_trad = VariantRegistry()
            generator_trad = VariantGenerator(registry_trad, random_seed=12345)
            insertion_ids_trad = generator_trad.generate_from_config(config)
            
            selector_trad = ReadSelector(bam_file, min_mapq=20, max_soft_clip_ratio=0.1)
            selected_reads_trad = selector_trad.select_reads(len(insertion_ids_trad))
            
            inserter_trad = ReadInserter(registry_trad, random_seed=12345)
            num_to_insert = min(len(selected_reads_trad), len(insertion_ids_trad))
            
            insertion_records_trad, modified_sequences_trad, _ = inserter_trad.insert_random_mode(
                selected_reads_trad[:num_to_insert], 
                insertion_ids_trad[:num_to_insert]
            )
            
            trad_fasta = temp_path / "traditional_medium.fasta"
            inserter_trad.save_modified_sequences(modified_sequences_trad, str(trad_fasta))
            
            # Memory-optimized
            registry_mem = VariantRegistry()
            generator_mem = VariantGenerator(registry_mem, random_seed=12345)
            insertion_ids_mem = generator_mem.generate_from_config(config)
            
            selector_mem = ReadSelector(bam_file, min_mapq=20, max_soft_clip_ratio=0.1)
            lazy_reads_mem = selector_mem.select_lazy_reads(len(insertion_ids_mem))
            
            inserter_mem = ReadInserter(registry_mem, random_seed=12345)
            num_to_insert_mem = min(len(lazy_reads_mem), len(insertion_ids_mem))
            
            results_mem = list(inserter_mem.insert_streaming_mode(
                lazy_reads_mem[:num_to_insert_mem], 
                insertion_ids_mem[:num_to_insert_mem]
            ))
            
            modified_sequences_mem = [result[1] for result in results_mem]
            
            mem_fasta = temp_path / "memory_optimized_medium.fasta"
            inserter_mem.save_modified_sequences(modified_sequences_mem, str(mem_fasta))
            
            # === Compare Results ===
            comparison_info = compare_fasta_content(str(trad_fasta), str(mem_fasta))
            
            # At medium scale, we expect high success rates (>95%)
            min_expected_sequences = int(100 * 0.95)  # 95 sequences
            
            self.assertGreaterEqual(len(modified_sequences_trad), min_expected_sequences,
                                  f"Traditional method should produce at least {min_expected_sequences} sequences")
            self.assertGreaterEqual(len(modified_sequences_mem), min_expected_sequences,
                                  f"Memory-optimized method should produce at least {min_expected_sequences} sequences")
            
            # The sequence content should be identical
            self.assertTrue(comparison_info['sequences_identical'], 
                          f"Medium-scale sequence content differs. Comparison: {comparison_info}")
            
            print(f"✓ Medium-scale equivalence verified: {comparison_info['common_sequences']} identical sequences")


if __name__ == '__main__':
    unittest.main()