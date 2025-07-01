"""
Unit tests for read_inserter module.
"""

import unittest
import tempfile
import json
from pathlib import Path
from unittest.mock import Mock
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from nova.variant_registry import VariantRegistry
from nova.read_inserter import ReadInserter, InsertionRecord
from nova.read_selector import ReadMetadata, LazyReadReference


class TestReadInserter(unittest.TestCase):
    """Test cases for ReadInserter class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.registry = VariantRegistry()
        self.inserter = ReadInserter(
            registry=self.registry,
            min_distance_from_ends=500,
            random_seed=42
        )
        
        # Add some test sequences to registry
        self.random_id = self.registry.add_sequence("ATCGATCG", "random", {"test": "data"})
        self.simple_id = self.registry.add_sequence("CAGCAGCAG", "simple", {"repeat": "CAG"})
    
    def test_get_valid_insertion_position(self):
        """Test getting valid insertion positions."""
        # Test feasible insertion
        read_length = 2000
        insertion_length = 100
        
        position = self.inserter._get_valid_insertion_position(read_length, insertion_length)
        self.assertIsNotNone(position)
        self.assertGreaterEqual(position, 500)  # min_distance_from_ends
        self.assertLessEqual(position, 1500)  # read_length - min_distance_from_ends
        
        # Test infeasible insertion (read too short)
        read_length = 800  # Too short for min_distance_from_ends = 500
        position = self.inserter._get_valid_insertion_position(read_length, insertion_length)
        self.assertIsNone(position)
    
    def test_generate_modified_read_name(self):
        """Test modified read name generation."""
        insertion_id = "test_insertion_123"
        base_name = "original_read"
        
        modified_name = self.inserter._generate_modified_read_name(insertion_id, base_name)
        expected = f"{insertion_id}.{base_name}"
        self.assertEqual(modified_name, expected)
    
    def test_insert_sequence_into_read(self):
        """Test sequence insertion into read."""
        read_sequence = "AAAAAAAAAATTTTTTTTTT"  # 20 bases
        insertion_sequence = "GGGG"
        insertion_pos = 10
        
        result = self.inserter._insert_sequence_into_read(
            read_sequence, insertion_sequence, insertion_pos
        )
        
        expected = "AAAAAAAAAAGGGGTTTTTTTTTT"  # GGGG inserted at position 10
        self.assertEqual(result, expected)
        self.assertEqual(len(result), len(read_sequence) + len(insertion_sequence))
    
    
    
    
    
    
    def test_insert_streaming_success(self):
        """Test successful streaming insertion."""
        # Create mock lazy read reference
        lazy_read = LazyReadReference(
            bam_path="test.bam",
            read_name="test_read",
            original_chr="chr1",
            original_pos=1000,
            read_length=2000,
            mapq=30,
            gc_content=0.0,
            soft_clip_ratio=0.05
        )
        
        # Mock the get_sequence method to return a test sequence
        lazy_read.get_sequence = Mock(return_value="A" * 2000)
        
        lazy_reads = [lazy_read]
        insertion_ids = [self.random_id]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as records_f, \
             tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as sequences_f:
            records_path = records_f.name
            sequences_path = sequences_f.name
        
        try:
            # Perform streaming insertion
            stats = self.inserter.insert_streaming(lazy_reads, insertion_ids, records_path, sequences_path)
            
            # Check statistics
            self.assertEqual(stats['total_insertions'], 1)
            
            # Verify files were created and contain correct data
            self.assertTrue(Path(records_path).exists())
            self.assertTrue(Path(sequences_path).exists())
            
            # Check records file
            with open(records_path, 'r') as f:
                records_data = json.load(f)
            self.assertEqual(len(records_data), 1)
            self.assertEqual(records_data[0]['base_read_name'], "test_read")
            self.assertEqual(records_data[0]['insertion_id'], self.random_id)
            
            # Check sequences file
            loaded_seqs = list(SeqIO.parse(sequences_path, "fasta"))
            self.assertEqual(len(loaded_seqs), 1)
            self.assertGreater(len(loaded_seqs[0].seq), 2000)  # Original + insertion
            self.assertIn("ATCGATCG", str(loaded_seqs[0].seq))  # Insertion sequence present
            
        finally:
            Path(records_path).unlink()
            Path(sequences_path).unlink()
    
    def test_stream_insertion_infeasible_read(self):
        """Test streaming insertion with infeasible read."""
        # Create lazy read that's too short
        lazy_read = LazyReadReference(
            bam_path="test.bam", 
            read_name="short_read",
            original_chr="chr1",
            original_pos=1000,
            read_length=800,  # Too short for min_distance_from_ends = 500
            mapq=30,
            gc_content=0.0,
            soft_clip_ratio=0.05
        )
        
        lazy_reads = [lazy_read]
        insertion_ids = [self.random_id]
        
        # Test _stream_insertion generator directly
        results = list(self.inserter._stream_insertion(lazy_reads, insertion_ids))
        
        # Should have no successful insertions
        self.assertEqual(len(results), 0)
    
    def test_stream_insertion_missing_sequence(self):
        """Test streaming insertion with missing sequence."""
        lazy_read = LazyReadReference(
            bam_path="test.bam",
            read_name="test_read", 
            original_chr="chr1",
            original_pos=1000,
            read_length=2000,
            mapq=30,
            gc_content=0.0,
            soft_clip_ratio=0.05
        )
        
        lazy_reads = [lazy_read]
        insertion_ids = ["nonexistent_id"]
        
        # Test _stream_insertion generator directly
        results = list(self.inserter._stream_insertion(lazy_reads, insertion_ids))
        
        # Should have no successful insertions
        self.assertEqual(len(results), 0)
    
    def test_get_insertion_statistics(self):
        """Test calculating insertion statistics."""
        # Create sample insertion records
        records = [
            InsertionRecord("read1", "mod_read1", "chr1", 1000, "id1", "random", 10, 500),
            InsertionRecord("read2", "mod_read2", "chr1", 2000, "id2", "random", 15, 600),
            InsertionRecord("read3", "mod_read3", "chr2", 3000, "id3", "simple", 20, 700),
        ]
        
        stats = self.inserter.get_insertion_statistics(records)
        
        # Check overall statistics
        self.assertEqual(stats['total_insertions'], 3)
        self.assertEqual(stats['type_counts']['random'], 2)
        self.assertEqual(stats['type_counts']['simple'], 1)
        
        # Check length statistics
        self.assertEqual(stats['length_statistics']['random']['count'], 2)
        self.assertEqual(stats['length_statistics']['random']['min'], 10)
        self.assertEqual(stats['length_statistics']['random']['max'], 15)
        self.assertEqual(stats['length_statistics']['random']['mean'], 12.5)
        
        self.assertEqual(stats['length_statistics']['simple']['count'], 1)
        self.assertEqual(stats['length_statistics']['simple']['min'], 20)
        self.assertEqual(stats['length_statistics']['simple']['max'], 20)
        self.assertEqual(stats['length_statistics']['simple']['mean'], 20)
    
    def test_get_insertion_statistics_empty(self):
        """Test statistics with empty records."""
        stats = self.inserter.get_insertion_statistics([])
        self.assertEqual(stats['total_insertions'], 0)


class TestInsertionRecord(unittest.TestCase):
    """Test cases for InsertionRecord dataclass."""
    
    def test_creation(self):
        """Test creating InsertionRecord."""
        record = InsertionRecord(
            base_read_name="test_read",
            modified_read_name="mod_test_read",
            original_chr="chr1",
            original_pos=1000,
            insertion_id="test_id",
            insertion_type="random",
            insertion_length=10,
            insertion_pos=500
        )
        
        self.assertEqual(record.base_read_name, "test_read")
        self.assertEqual(record.modified_read_name, "mod_test_read")
        self.assertEqual(record.original_chr, "chr1")
        self.assertEqual(record.original_pos, 1000)
        self.assertEqual(record.insertion_id, "test_id")
        self.assertEqual(record.insertion_type, "random")
        self.assertEqual(record.insertion_length, 10)
        self.assertEqual(record.insertion_pos, 500)
    
    def test_to_dict(self):
        """Test converting to dictionary."""
        record = InsertionRecord(
            base_read_name="test_read",
            modified_read_name="mod_test_read",
            original_chr="chr1",
            original_pos=1000,
            insertion_id="test_id",
            insertion_type="random",
            insertion_length=10,
            insertion_pos=500
        )
        
        result = record.to_dict()
        expected = {
            'base_read_name': 'test_read',
            'modified_read_name': 'mod_test_read',
            'original_chr': 'chr1',
            'original_pos': 1000,
            'insertion_id': 'test_id',
            'insertion_type': 'random',
            'insertion_length': 10,
            'insertion_pos': 500
        }
        
        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()