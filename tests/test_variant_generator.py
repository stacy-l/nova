"""
Unit tests for variant_generator module.
"""

import unittest
import tempfile
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from nova.variant_registry import VariantRegistry
from nova.variant_generator import VariantGenerator


class TestVariantGenerator(unittest.TestCase):
    """Test cases for VariantGenerator class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.registry = VariantRegistry()
        self.generator = VariantGenerator(self.registry, random_seed=42)
    
    def test_generate_random_insertions(self):
        """Test generating random insertions."""
        n = 5
        length = 100
        
        insertion_ids = self.generator.generate_random_insertions(n, length)
        
        self.assertEqual(len(insertion_ids), n)
        self.assertEqual(len(self.registry), n)
        
        # Check each insertion
        for insertion_id in insertion_ids:
            seq = self.registry.get_sequence(insertion_id)
            self.assertIsNotNone(seq)
            self.assertEqual(seq.insertion_length, length)
            self.assertEqual(seq.insertion_type, "random")
            self.assertEqual(len(seq.sequence), length)
            # Check sequence only contains valid nucleotides
            self.assertTrue(all(c in 'ATGC' for c in seq.sequence))
    
    def test_generate_random_insertions_with_gc_content(self):
        """Test generating random insertions with specified GC content."""
        n = 10
        length = 1000
        gc_content = 0.6
        
        insertion_ids = self.generator.generate_random_insertions(n, length, gc_content)
        
        # Check GC content is approximately correct
        for insertion_id in insertion_ids:
            seq = self.registry.get_sequence(insertion_id)
            actual_gc = (seq.sequence.count('G') + seq.sequence.count('C')) / len(seq.sequence)
            self.assertAlmostEqual(actual_gc, gc_content, delta=0.1)
            self.assertEqual(seq.metadata['target_gc_content'], gc_content)
    
    def test_generate_simple_repeat_insertions(self):
        """Test generating simple repeat insertions."""
        n = 3
        repeat_unit = "CAG"
        units = 40
        expected_sequence = repeat_unit * units
        
        insertion_ids = self.generator.generate_simple_repeat_insertions(n, repeat_unit, units)
        
        self.assertEqual(len(insertion_ids), n)
        self.assertEqual(len(self.registry), n)
        
        # Check each insertion
        for insertion_id in insertion_ids:
            seq = self.registry.get_sequence(insertion_id)
            self.assertIsNotNone(seq)
            self.assertEqual(seq.sequence, expected_sequence)
            self.assertEqual(seq.insertion_type, "simple")
            self.assertEqual(seq.metadata['repeat_unit'], repeat_unit)
            self.assertEqual(seq.metadata['repeat_units'], units)
    
    def test_generate_predefined_insertions(self):
        """Test generating predefined insertions from FASTA."""
        # Create temporary FASTA file
        sequences = [
            SeqRecord(Seq("ATCGATCGATCG"), id="AluYa5", description=""),
            SeqRecord(Seq("GCTAGCTAGCTA"), id="AluYb8", description=""),
        ]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            temp_fasta = f.name
        
        try:
            SeqIO.write(sequences, temp_fasta, "fasta")
            
            sequence_counts = {"AluYa5": 2, "AluYb8": 3}
            insertion_ids = self.generator.generate_predefined_insertions(temp_fasta, sequence_counts)
            
            self.assertEqual(len(insertion_ids), 5)  # 2 + 3
            self.assertEqual(len(self.registry), 5)
            
            # Check AluYa5 sequences
            aluya5_seqs = self.registry.get_sequences_by_type("AluYa5")
            self.assertEqual(len(aluya5_seqs), 2)
            for seq in aluya5_seqs:
                self.assertEqual(seq.sequence, "ATCGATCGATCG")
                self.assertEqual(seq.metadata['source_sequence_name'], "AluYa5")
            
            # Check AluYb8 sequences
            aluyb8_seqs = self.registry.get_sequences_by_type("AluYb8")
            self.assertEqual(len(aluyb8_seqs), 3)
            for seq in aluyb8_seqs:
                self.assertEqual(seq.sequence, "GCTAGCTAGCTA")
                self.assertEqual(seq.metadata['source_sequence_name'], "AluYb8")
        
        finally:
            Path(temp_fasta).unlink()
    
    def test_generate_from_config(self):
        """Test generating from configuration dictionary."""
        # Create temporary FASTA file
        sequences = [
            SeqRecord(Seq("ATCGATCGATCG"), id="AluYa5", description=""),
        ]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            temp_fasta = f.name
        
        try:
            SeqIO.write(sequences, temp_fasta, "fasta")
            
            config = {
                'random': {'n': 2, 'length': 50},
                'simple': {'n': 3, 'repeat': 'CAG', 'units': 10},
                'predefined': {
                    'Alu': {
                        'fasta': temp_fasta,
                        'spec': {'AluYa5': 2}
                    }
                }
            }
            
            insertion_ids = self.generator.generate_from_config(config)
            
            self.assertEqual(len(insertion_ids), 7)  # 2 + 3 + 2
            self.assertEqual(len(self.registry), 7)
            
            # Check types
            random_seqs = self.registry.get_sequences_by_type("random")
            simple_seqs = self.registry.get_sequences_by_type("simple")
            alu_seqs = self.registry.get_sequences_by_type("AluYa5")
            
            self.assertEqual(len(random_seqs), 2)
            self.assertEqual(len(simple_seqs), 3)
            self.assertEqual(len(alu_seqs), 2)
        
        finally:
            Path(temp_fasta).unlink()
    
    def test_validate_config_valid(self):
        """Test validating valid configuration."""
        config = {
            'random': {'n': 5, 'length': 100, 'gc_content': 0.5},
            'simple': {'n': 3, 'repeat': 'CAG', 'units': 40},
            'predefined': {
                'Alu': {
                    'fasta': 'test.fasta',
                    'spec': {'AluYa5': 2}
                }
            }
        }
        
        errors = self.generator.validate_config(config)
        self.assertEqual(len(errors), 0)
    
    def test_validate_config_invalid(self):
        """Test validating invalid configuration."""
        config = {
            'random': {'n': -1, 'length': 0, 'gc_content': 1.5},
            'simple': {'n': 0, 'repeat': '', 'units': -1},
            'predefined': {
                'Alu': {
                    'fasta': 'test.fasta'
                    # missing 'spec'
                }
            }
        }
        
        errors = self.generator.validate_config(config)
        self.assertGreater(len(errors), 0)
        
        # Check for specific error types
        error_text = ' '.join(errors)
        self.assertIn("positive", error_text)
        self.assertIn("empty", error_text)
        self.assertIn("spec", error_text)
    
    def test_validate_config_not_dict(self):
        """Test validating non-dictionary configuration."""
        errors = self.generator.validate_config("not a dict")
        self.assertEqual(len(errors), 1)
    
    def test_generate_from_config_multi_entry(self):
        """Test generating from configuration with multiple entries per type."""
        config = {
            'random': [
                {'n': 5, 'length': 100, 'gc_content': 0.4},
                {'n': 3, 'length': 200}
            ],
            'simple': [
                {'n': 4, 'repeat': 'CAG', 'units': 20},
                {'n': 2, 'repeat': 'GC', 'units': 10}
            ]
        }
        
        insertion_ids = self.generator.generate_from_config(config)
        
        # Should generate 5+3+4+2 = 14 total insertions
        self.assertEqual(len(insertion_ids), 14)
        self.assertEqual(len(self.registry), 14)
        
        # Check that we have the right types
        random_seqs = self.registry.get_sequences_by_type('random')
        simple_seqs = self.registry.get_sequences_by_type('simple')
        
        self.assertEqual(len(random_seqs), 8)  # 5 + 3
        self.assertEqual(len(simple_seqs), 6)  # 4 + 2
        
        # Check that different lengths were generated for random sequences
        random_lengths = [seq.insertion_length for seq in random_seqs]
        self.assertIn(100, random_lengths)
        self.assertIn(200, random_lengths)
        
        # Check that different repeat units were generated for simple sequences
        simple_metadata = [seq.metadata for seq in simple_seqs]
        repeat_units = [meta.get('repeat_unit') for meta in simple_metadata]
        self.assertIn('CAG', repeat_units)
        self.assertIn('GC', repeat_units)
    
    def test_validate_config_multi_entry_valid(self):
        """Test validating valid multi-entry configuration."""
        config = {
            'random': [
                {'n': 5, 'length': 100, 'gc_content': 0.5},
                {'n': 3, 'length': 200}
            ],
            'simple': [
                {'n': 4, 'repeat': 'CAG', 'units': 40},
                {'n': 2, 'repeat': 'GC', 'units': 10}
            ]
        }
        
        errors = self.generator.validate_config(config)
        self.assertEqual(len(errors), 0)
    
    def test_validate_config_multi_entry_invalid(self):
        """Test validating invalid multi-entry configuration."""
        config = {
            'random': [
                {'n': -1, 'length': 100},  # Invalid n
                {'n': 3, 'length': 0}      # Invalid length
            ],
            'simple': [
                {'n': 4, 'repeat': '', 'units': 40},  # Invalid repeat
                {'n': 0, 'repeat': 'GC', 'units': -1}  # Invalid n and units
            ]
        }
        
        errors = self.generator.validate_config(config)
        self.assertGreater(len(errors), 0)
        
        # Check that multiple configs are identified in error messages
        error_text = ' '.join(errors)
        self.assertIn("Random config 1:", error_text)
        self.assertIn("Random config 2:", error_text)
        self.assertIn("Simple config 1:", error_text)
        self.assertIn("Simple config 2:", error_text)
    
    def test_backwards_compatibility_single_entry(self):
        """Test that single-entry (old format) configs still work."""
        config = {
            'random': {'n': 5, 'length': 100, 'gc_content': 0.5},
            'simple': {'n': 3, 'repeat': 'CAG', 'units': 40}
        }
        
        insertion_ids = self.generator.generate_from_config(config)
        
        # Should generate 5+3 = 8 total insertions
        self.assertEqual(len(insertion_ids), 8)
        self.assertEqual(len(self.registry), 8)
        
        # Validation should also work
        errors = self.generator.validate_config(config)
        self.assertEqual(len(errors), 0)
    
    def test_normalize_to_list_helper(self):
        """Test the _normalize_to_list helper method."""
        # Test with single object
        single_obj = {'n': 5, 'length': 100}
        result = self.generator._normalize_to_list(single_obj)
        self.assertEqual(result, [single_obj])
        
        # Test with list
        list_obj = [{'n': 5, 'length': 100}, {'n': 3, 'length': 200}]
        result = self.generator._normalize_to_list(list_obj)
        self.assertEqual(result, list_obj)


if __name__ == '__main__':
    unittest.main()