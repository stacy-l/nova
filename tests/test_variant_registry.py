"""
Unit tests for variant_registry module.
"""

import unittest
import tempfile
import json
from pathlib import Path

from nova.variant_registry import VariantRegistry, InsertionSequence


class TestVariantRegistry(unittest.TestCase):
    """Test cases for VariantRegistry class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.registry = VariantRegistry()
    
    def test_add_sequence(self):
        """Test adding sequences to registry."""
        sequence = "ATCGATCG"
        insertion_type = "random"
        metadata = {"test": "value"}
        
        insertion_id = self.registry.add_sequence(sequence, insertion_type, metadata)
        
        self.assertIsInstance(insertion_id, str)
        self.assertTrue(insertion_id.startswith(insertion_type))
        self.assertEqual(len(self.registry), 1)
        self.assertIn(insertion_id, self.registry)
    
    def test_get_sequence(self):
        """Test retrieving sequences from registry."""
        sequence = "ATCGATCG"
        insertion_type = "random"
        
        insertion_id = self.registry.add_sequence(sequence, insertion_type)
        retrieved = self.registry.get_sequence(insertion_id)
        
        self.assertIsNotNone(retrieved)
        self.assertEqual(retrieved.sequence, sequence)
        self.assertEqual(retrieved.insertion_type, insertion_type)
        self.assertEqual(retrieved.insertion_length, len(sequence))
    
    def test_get_nonexistent_sequence(self):
        """Test retrieving non-existent sequence."""
        result = self.registry.get_sequence("nonexistent_id")
        self.assertIsNone(result)
    
    def test_get_sequences_by_type(self):
        """Test filtering sequences by type."""
        self.registry.add_sequence("ATCG", "random")
        self.registry.add_sequence("GCTA", "random")
        self.registry.add_sequence("AAAA", "simple")
        
        random_seqs = self.registry.get_sequences_by_type("random")
        simple_seqs = self.registry.get_sequences_by_type("simple")
        
        self.assertEqual(len(random_seqs), 2)
        self.assertEqual(len(simple_seqs), 1)
        self.assertEqual(simple_seqs[0].sequence, "AAAA")
    
    def test_statistics(self):
        """Test registry statistics."""
        self.registry.add_sequence("ATCG", "random")
        self.registry.add_sequence("ATCGATCG", "random")
        self.registry.add_sequence("AAAA", "simple")
        
        stats = self.registry.get_statistics()
        
        self.assertEqual(stats['total_sequences'], 3)
        self.assertEqual(stats['type_counts']['random'], 2)
        self.assertEqual(stats['type_counts']['simple'], 1)
        self.assertEqual(stats['length_statistics']['random']['count'], 2)
        self.assertEqual(stats['length_statistics']['simple']['mean'], 4)
    
    def test_save_load_json(self):
        """Test saving and loading registry to/from JSON."""
        # Add some sequences
        id1 = self.registry.add_sequence("ATCG", "random", {"test": "value1"})
        id2 = self.registry.add_sequence("GCTA", "simple", {"test": "value2"})
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            temp_path = f.name
        
        try:
            # Save to JSON
            self.registry.save_to_json(temp_path)
            
            # Create new registry and load
            new_registry = VariantRegistry()
            new_registry.load_from_json(temp_path)
            
            # Verify loaded data
            self.assertEqual(len(new_registry), 2)
            self.assertIn(id1, new_registry)
            self.assertIn(id2, new_registry)
            
            seq1 = new_registry.get_sequence(id1)
            self.assertEqual(seq1.sequence, "ATCG")
            self.assertEqual(seq1.metadata["test"], "value1")
            
        finally:
            Path(temp_path).unlink()
    
    def test_remove_sequence(self):
        """Test removing sequences from registry."""
        insertion_id = self.registry.add_sequence("ATCG", "random")
        
        self.assertTrue(self.registry.remove_sequence(insertion_id))
        self.assertEqual(len(self.registry), 0)
        self.assertNotIn(insertion_id, self.registry)
        
        # Try to remove non-existent sequence
        self.assertFalse(self.registry.remove_sequence("nonexistent"))
    
    def test_clear(self):
        """Test clearing registry."""
        self.registry.add_sequence("ATCG", "random")
        self.registry.add_sequence("GCTA", "simple")
        
        self.assertEqual(len(self.registry), 2)
        
        self.registry.clear()
        self.assertEqual(len(self.registry), 0)


class TestInsertionSequence(unittest.TestCase):
    """Test cases for InsertionSequence class."""
    
    def test_creation(self):
        """Test creating InsertionSequence."""
        seq = InsertionSequence(
            insertion_id="test_id",
            sequence="ATCG",
            insertion_type="random",
            insertion_length=4,
            metadata={"test": "value"}
        )
        
        self.assertEqual(seq.insertion_id, "test_id")
        self.assertEqual(seq.sequence, "ATCG")
        self.assertEqual(seq.insertion_type, "random")
        self.assertEqual(seq.insertion_length, 4)
        self.assertEqual(seq.metadata, {"test": "value"})
    
    def test_to_dict(self):
        """Test converting to dictionary."""
        seq = InsertionSequence(
            insertion_id="test_id",
            sequence="ATCG",
            insertion_type="random",
            insertion_length=4,
            metadata={"test": "value"}
        )
        
        result = seq.to_dict()
        expected = {
            'insertion_id': 'test_id',
            'sequence': 'ATCG',
            'insertion_type': 'random',
            'insertion_length': 4,
            'metadata': {'test': 'value'}
        }
        
        self.assertEqual(result, expected)
    
    def test_from_dict(self):
        """Test creating from dictionary."""
        data = {
            'insertion_id': 'test_id',
            'sequence': 'ATCG',
            'insertion_type': 'random',
            'insertion_length': 4,
            'metadata': {'test': 'value'}
        }
        
        seq = InsertionSequence.from_dict(data)
        
        self.assertEqual(seq.insertion_id, "test_id")
        self.assertEqual(seq.sequence, "ATCG")
        self.assertEqual(seq.insertion_type, "random")
        self.assertEqual(seq.insertion_length, 4)
        self.assertEqual(seq.metadata, {"test": "value"})


if __name__ == '__main__':
    unittest.main()