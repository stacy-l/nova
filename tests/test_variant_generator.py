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
    
    def test_mutation_seed_generation(self):
        """Test mutation seed generation is deterministic and unique."""
        # Same parameters should produce same seed
        seed1 = self.generator._generate_mutation_seed(42, 0)
        seed2 = self.generator._generate_mutation_seed(42, 0)
        self.assertEqual(seed1, seed2)
        
        # Different sequence indices should produce different seeds
        seed3 = self.generator._generate_mutation_seed(42, 1)
        self.assertNotEqual(seed1, seed3)
        
        # Different base seeds should produce different seeds
        seed4 = self.generator._generate_mutation_seed(99, 0)
        self.assertNotEqual(seed1, seed4)
    
    def test_apply_mutations_basic(self):
        """Test basic mutation application."""
        sequence = "ATGCATGCATGC"  # 12 bases
        mutation_config = {"substitution_rate": 0.25}  # 3 mutations expected
        
        mutated_seq, metadata = self.generator._apply_mutations(sequence, mutation_config, 0)
        
        # Check metadata
        self.assertEqual(metadata.substitution_rate, 0.25)
        self.assertGreater(metadata.total_mutations, 0)
        self.assertLessEqual(metadata.total_mutations, len(sequence))
        self.assertIsNotNone(metadata.mutation_seed)
        
        # Check that sequence changed (with high probability)
        self.assertNotEqual(sequence, mutated_seq)
        self.assertEqual(len(sequence), len(mutated_seq))
        
        # Check that mutations are recorded
        self.assertEqual(len(metadata.mutations_applied), metadata.total_mutations)
        
        # Verify each mutation record
        for mutation in metadata.mutations_applied:
            self.assertIn(mutation.original_base, 'ATGC')
            self.assertIn(mutation.mutated_base, 'ATGC')
            self.assertNotEqual(mutation.original_base, mutation.mutated_base)
            self.assertLess(mutation.position, len(sequence))
    
    def test_apply_mutations_no_config(self):
        """Test mutation application with no configuration."""
        sequence = "ATGCATGC"
        
        # Test with None config
        mutated_seq, metadata = self.generator._apply_mutations(sequence, None, 0)
        self.assertEqual(sequence, mutated_seq)
        self.assertEqual(metadata.total_mutations, 0)
        
        # Test with empty config
        mutated_seq, metadata = self.generator._apply_mutations(sequence, {}, 0)
        self.assertEqual(sequence, mutated_seq)
        self.assertEqual(metadata.total_mutations, 0)
        
        # Test with zero rate
        mutation_config = {"substitution_rate": 0.0}
        mutated_seq, metadata = self.generator._apply_mutations(sequence, mutation_config, 0)
        self.assertEqual(sequence, mutated_seq)
        self.assertEqual(metadata.total_mutations, 0)
    
    def test_apply_mutations_deterministic(self):
        """Test that mutations are deterministic with same seed."""
        sequence = "ATGCATGCATGCATGC"
        mutation_config = {"substitution_rate": 0.2}
        
        # Same sequence index should produce same mutations
        mutated1, metadata1 = self.generator._apply_mutations(sequence, mutation_config, 5)
        mutated2, metadata2 = self.generator._apply_mutations(sequence, mutation_config, 5)
        
        self.assertEqual(mutated1, mutated2)
        self.assertEqual(metadata1.mutation_seed, metadata2.mutation_seed)
        self.assertEqual(len(metadata1.mutations_applied), len(metadata2.mutations_applied))
        
        # Different sequence index should produce different mutations
        mutated3, metadata3 = self.generator._apply_mutations(sequence, mutation_config, 10)
        self.assertNotEqual(mutated1, mutated3)
        self.assertNotEqual(metadata1.mutation_seed, metadata3.mutation_seed)
    
    def test_mutation_validation(self):
        """Test mutation configuration validation."""
        # Valid mutation config
        valid_config = {"substitution_rate": 0.02}
        errors = self.generator._validate_mutation_config(valid_config, "Test: ")
        self.assertEqual(len(errors), 0)
        
        # Invalid substitution rate - too high
        invalid_config = {"substitution_rate": 1.5}
        errors = self.generator._validate_mutation_config(invalid_config, "Test: ")
        self.assertGreater(len(errors), 0)
        self.assertIn("between 0 and 1", errors[0])
        
        # Invalid substitution rate - negative
        invalid_config = {"substitution_rate": -0.1}
        errors = self.generator._validate_mutation_config(invalid_config, "Test: ")
        self.assertGreater(len(errors), 0)
        
        # Invalid substitution rate - not a number
        invalid_config = {"substitution_rate": "not_a_number"}
        errors = self.generator._validate_mutation_config(invalid_config, "Test: ")
        self.assertGreater(len(errors), 0)
    
    def test_random_insertions_with_mutations(self):
        """Test random insertion generation with mutations."""
        n = 5
        length = 50
        mutation_config = {"substitution_rate": 0.1}
        
        insertion_ids = self.generator.generate_random_insertions(
            n, length, gc_content=0.5, mutation_config=mutation_config
        )
        
        self.assertEqual(len(insertion_ids), n)
        
        # Check that mutations were applied
        for insertion_id in insertion_ids:
            seq = self.registry.get_sequence(insertion_id)
            self.assertIsNotNone(seq)
            self.assertEqual(seq.insertion_length, length)
            
            # Should have mutation metadata
            self.assertIn('mutations', seq.metadata)
            mut_meta = seq.metadata['mutations']
            self.assertEqual(mut_meta['substitution_rate'], 0.1)
            self.assertGreater(mut_meta['total_mutations'], 0)
    
    def test_simple_insertions_with_mutations(self):
        """Test simple repeat insertion generation with mutations."""
        n = 3
        repeat_unit = "CAG"
        units = 10
        mutation_config = {"substitution_rate": 0.05}
        
        insertion_ids = self.generator.generate_simple_repeat_insertions(
            n, repeat_unit, units, mutation_config
        )
        
        self.assertEqual(len(insertion_ids), n)
        
        # Check that mutations were applied to each sequence
        for insertion_id in insertion_ids:
            seq = self.registry.get_sequence(insertion_id)
            self.assertIsNotNone(seq)
            self.assertEqual(seq.insertion_length, len(repeat_unit) * units)
            
            # Should have mutation metadata
            self.assertIn('mutations', seq.metadata)
            mut_meta = seq.metadata['mutations']
            self.assertEqual(mut_meta['substitution_rate'], 0.05)
    
    def test_config_with_mutations_integration(self):
        """Test end-to-end configuration with mutations."""
        config = {
            'random': [
                {
                    'n': 3,
                    'length': 50,
                    'mutations': {'substitution_rate': 0.02}
                },
                {
                    'n': 2,
                    'length': 100
                    # No mutations for this one
                }
            ],
            'simple': [
                {
                    'n': 2,
                    'repeat': 'GC',
                    'units': 5,
                    'mutations': {'substitution_rate': 0.1}
                }
            ]
        }
        
        # Validate config
        errors = self.generator.validate_config(config)
        self.assertEqual(len(errors), 0)
        
        # Generate sequences
        insertion_ids = self.generator.generate_from_config(config)
        self.assertEqual(len(insertion_ids), 7)  # 3 + 2 + 2
        
        # Check mutation metadata presence
        mutated_count = 0
        unmutated_count = 0
        
        for insertion_id in insertion_ids:
            seq = self.registry.get_sequence(insertion_id)
            if 'mutations' in seq.metadata:
                mutated_count += 1
            else:
                unmutated_count += 1
        
        # Should have 5 mutated (3 + 2) and 2 unmutated
        self.assertEqual(mutated_count, 5)
        self.assertEqual(unmutated_count, 2)
    
    def test_mutation_config_file_integration(self):
        """Test loading and using actual mutation configuration file."""
        import json
        from pathlib import Path
        
        # Load the test mutation config file
        config_path = Path(__file__).parent / "test_data" / "test_simple_mutation_config.json"
        
        with open(config_path, 'r') as f:
            config = json.load(f)
        
        # Validate the configuration
        errors = self.generator.validate_config(config)
        self.assertEqual(len(errors), 0, f"Config validation failed: {errors}")
        
        # Generate sequences from the config
        insertion_ids = self.generator.generate_from_config(config)
        
        # Verify expected counts: 10 + 10 + 5 = 25 total
        self.assertEqual(len(insertion_ids), 25)
        
        # Verify that all sequences have mutation metadata
        for insertion_id in insertion_ids:
            seq = self.registry.get_sequence(insertion_id)
            self.assertIn('mutations', seq.metadata, 
                         f"Sequence {insertion_id} missing mutation metadata")
            
            # Check mutation metadata structure
            mut_meta = seq.metadata['mutations']
            self.assertIn('substitution_rate', mut_meta)
            self.assertIn('total_mutations', mut_meta)
            self.assertIn('mutation_seed', mut_meta)
            self.assertIn('mutation_records', mut_meta)
            
            # Verify substitution rate is as expected
            if seq.insertion_type == 'random':
                if seq.insertion_length == 50:
                    self.assertEqual(mut_meta['substitution_rate'], 0.02)
                else:  # length 100
                    self.assertEqual(mut_meta['substitution_rate'], 0.05)
            elif seq.insertion_type == 'simple':
                self.assertEqual(mut_meta['substitution_rate'], 0.03)
            
            # Verify mutations were actually applied
            self.assertGreater(mut_meta['total_mutations'], 0)


if __name__ == '__main__':
    unittest.main()