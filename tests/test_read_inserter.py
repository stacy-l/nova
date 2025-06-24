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
from nova.read_selector import ReadMetadata


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
        expected = f"{insertion_id}_{base_name}"
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
    
    def test_insert_random_mode_success(self):
        """Test successful random insertion."""
        # Create mock read and metadata
        mock_read = Mock()
        mock_read.query_name = "test_read"
        mock_read.query_length = 2000
        mock_read.query_sequence = "A" * 2000
        
        metadata = ReadMetadata(
            read_name="test_read",
            original_chr="chr1",
            original_pos=1000,
            read_length=2000,
            mapq=30,
            gc_content=0.0,
            soft_clip_ratio=0.05
        )
        
        reads_with_metadata = [(mock_read, metadata)]
        insertion_ids = [self.random_id]
        
        # Perform insertion
        insertion_records, modified_sequences, skip_stats = self.inserter.insert_random_mode(
            reads_with_metadata, insertion_ids
        )
        
        # Check results
        self.assertEqual(len(insertion_records), 1)
        self.assertEqual(len(modified_sequences), 1)
        self.assertEqual(skip_stats['successful_insertions'], 1)
        self.assertEqual(skip_stats['success_rate'], 1.0)
        
        # Check insertion record
        record = insertion_records[0]
        self.assertEqual(record.base_read_name, "test_read")
        self.assertEqual(record.insertion_id, self.random_id)
        self.assertEqual(record.original_chr, "chr1")
        self.assertEqual(record.original_pos, 1000)
        self.assertGreater(record.insertion_pos, 500)
        self.assertLess(record.insertion_pos, 1500)
        
        # Check modified sequence
        seq_record = modified_sequences[0]
        self.assertGreater(len(seq_record.seq), 2000)  # Original + insertion
        self.assertIn("ATCGATCG", str(seq_record.seq))  # Insertion sequence present
    
    def test_insert_random_mode_infeasible_read(self):
        """Test insertion with infeasible read (too short)."""
        # Create mock read that's too short
        mock_read = Mock()
        mock_read.query_name = "short_read"
        mock_read.query_length = 800  # Too short for min_distance_from_ends = 500
        mock_read.query_sequence = "A" * 800
        
        metadata = ReadMetadata(
            read_name="short_read",
            original_chr="chr1",
            original_pos=1000,
            read_length=800,
            mapq=30,
            gc_content=0.0,
            soft_clip_ratio=0.05
        )
        
        reads_with_metadata = [(mock_read, metadata)]
        insertion_ids = [self.random_id]
        
        # Perform insertion
        insertion_records, modified_sequences, skip_stats = self.inserter.insert_random_mode(
            reads_with_metadata, insertion_ids
        )
        
        # Should have no successful insertions
        self.assertEqual(len(insertion_records), 0)
        self.assertEqual(len(modified_sequences), 0)
        self.assertEqual(skip_stats['successful_insertions'], 0)
        self.assertEqual(skip_stats['success_rate'], 0.0)
        self.assertEqual(skip_stats['skipped_infeasible_reads'], 1)
    
    def test_insert_random_mode_missing_sequence(self):
        """Test insertion with missing sequence in registry."""
        mock_read = Mock()
        mock_read.query_name = "test_read"
        mock_read.query_length = 2000
        mock_read.query_sequence = "A" * 2000
        
        metadata = ReadMetadata(
            read_name="test_read",
            original_chr="chr1",
            original_pos=1000,
            read_length=2000,
            mapq=30,
            gc_content=0.0,
            soft_clip_ratio=0.05
        )
        
        reads_with_metadata = [(mock_read, metadata)]
        insertion_ids = ["nonexistent_id"]
        
        # Perform insertion
        insertion_records, modified_sequences, skip_stats = self.inserter.insert_random_mode(
            reads_with_metadata, insertion_ids
        )
        
        # Should have no successful insertions
        self.assertEqual(len(insertion_records), 0)
        self.assertEqual(len(modified_sequences), 0)
        self.assertEqual(skip_stats['successful_insertions'], 0)
        self.assertEqual(skip_stats['skipped_missing_sequences'], 1)
    
    def test_save_insertion_records(self):
        """Test saving insertion records to JSON."""
        # Create sample insertion record
        record = InsertionRecord(
            base_read_name="test_read",
            modified_read_name="insertion_123_test_read",
            original_chr="chr1",
            original_pos=1000,
            insertion_id="insertion_123",
            insertion_type="random",
            insertion_length=8,
            insertion_pos=750
        )
        
        records = [record]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            temp_path = f.name
        
        try:
            # Save records
            self.inserter.save_insertion_records(records, temp_path)
            
            # Verify file was created and contains correct data
            self.assertTrue(Path(temp_path).exists())
            
            with open(temp_path, 'r') as f:
                loaded_data = json.load(f)
            
            self.assertEqual(len(loaded_data), 1)
            self.assertEqual(loaded_data[0]['base_read_name'], "test_read")
            self.assertEqual(loaded_data[0]['insertion_id'], "insertion_123")
            
        finally:
            Path(temp_path).unlink()
    
    def test_save_modified_sequences(self):
        """Test saving modified sequences to FASTA."""
        # Create sample sequence records
        sequences = [
            SeqRecord(Seq("ATCGATCGATCG"), id="seq1", description="Test sequence 1"),
            SeqRecord(Seq("GCTAGCTAGCTA"), id="seq2", description="Test sequence 2")
        ]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            temp_path = f.name
        
        try:
            # Save sequences
            self.inserter.save_modified_sequences(sequences, temp_path)
            
            # Verify file was created and contains correct data
            self.assertTrue(Path(temp_path).exists())
            
            loaded_seqs = list(SeqIO.parse(temp_path, "fasta"))
            self.assertEqual(len(loaded_seqs), 2)
            self.assertEqual(str(loaded_seqs[0].seq), "ATCGATCGATCG")
            self.assertEqual(str(loaded_seqs[1].seq), "GCTAGCTAGCTA")
            
        finally:
            Path(temp_path).unlink()
    
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
        
        # Check position statistics
        self.assertEqual(stats['position_statistics']['random']['count'], 2)
        self.assertEqual(stats['position_statistics']['random']['min'], 500)
        self.assertEqual(stats['position_statistics']['random']['max'], 600)
        self.assertEqual(stats['position_statistics']['random']['mean'], 550)
        
        # Check length statistics
        self.assertEqual(stats['length_statistics']['random']['count'], 2)
        self.assertEqual(stats['length_statistics']['random']['min'], 10)
        self.assertEqual(stats['length_statistics']['random']['max'], 15)
        self.assertEqual(stats['length_statistics']['random']['mean'], 12.5)
    
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