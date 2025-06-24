"""
Read insertion module for inserting generated sequences into reads.
"""

import random
import json
from typing import List, Dict, Any, Tuple, Optional
from dataclasses import dataclass, asdict
import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import logging

from .variant_registry import VariantRegistry, InsertionSequence
from .read_selector import ReadMetadata


@dataclass
class InsertionRecord:
    """Record of an insertion made into a read."""
    base_read_name: str
    modified_read_name: str
    original_chr: str
    original_pos: int
    insertion_id: str
    insertion_type: str
    insertion_length: int
    insertion_pos: int
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class ReadInserter:
    """
    Insert generated sequences into reads.
    """
    
    def __init__(self, registry: VariantRegistry, 
                 min_distance_from_ends: int = 1000,
                 random_seed: Optional[int] = None):
        """
        Initialize ReadInserter.
        
        Args:
            registry: VariantRegistry containing insertion sequences
            min_distance_from_ends: Minimum distance from read ends for insertion
            random_seed: Optional seed for reproducible random positioning
        """
        self.registry = registry
        self.min_distance_from_ends = min_distance_from_ends
        self.logger = logging.getLogger(__name__)
        
        if random_seed is not None:
            random.seed(random_seed)
    
    def _get_valid_insertion_position(self, read_length: int, 
                                    insertion_length: int) -> Optional[int]:
        """
        Get a valid random insertion position within the read.
        
        Args:
            read_length: Length of the target read
            insertion_length: Length of sequence to insert
            
        Returns:
            Random insertion position or None if not feasible
        """
        start_pos = self.min_distance_from_ends
        end_pos = read_length - self.min_distance_from_ends
        
        if end_pos <= start_pos:
            self.logger.warning(f"Read too short ({read_length}) for insertion with "
                              f"minimum distance {self.min_distance_from_ends}")
            return None
        
        return random.randint(start_pos, end_pos)
    
    def _generate_modified_read_name(self, insertion_id: str, base_read_name: str) -> str:
        """
        Generate a modified read name incorporating the insertion ID.
        
        Args:
            insertion_id: Unique identifier for the insertion
            base_read_name: Original read name
            
        Returns:
            Modified read name
        """
        return f"{insertion_id}_{base_read_name}"
    
    def _insert_sequence_into_read(self, read_sequence: str, 
                                 insertion_sequence: str, 
                                 insertion_pos: int) -> str:
        """
        Insert sequence into read at specified position.
        
        Args:
            read_sequence: Original read sequence
            insertion_sequence: Sequence to insert
            insertion_pos: Position for insertion
            
        Returns:
            Modified read sequence
        """
        return (read_sequence[:insertion_pos] + 
                insertion_sequence + 
                read_sequence[insertion_pos:])
    
    def insert_random_mode(self, reads_with_metadata: List[Tuple[pysam.AlignedSegment, ReadMetadata]],
                          insertion_ids: List[str]) -> Tuple[List[InsertionRecord], List[SeqRecord]]:
        """
        Insert sequences into reads using random positioning.
        
        Args:
            reads_with_metadata: List of (read, metadata) tuples
            insertion_ids: List of insertion IDs to use
            
        Returns:
            Tuple of (insertion_records, modified_sequences, skip_statistics)
        """
        if len(reads_with_metadata) != len(insertion_ids):
            raise ValueError(f"Number of reads ({len(reads_with_metadata)}) must match "
                           f"number of insertions ({len(insertion_ids)})")
        
        insertion_records = []
        modified_sequences = []
        skipped_missing_sequences = []
        skipped_infeasible_reads = []
        
        for (read, metadata), insertion_id in zip(reads_with_metadata, insertion_ids):
            insertion_seq = self.registry.get_sequence(insertion_id)
            if insertion_seq is None:
                self.logger.error(f"Insertion sequence {insertion_id} not found in registry")
                skipped_missing_sequences.append(insertion_id)
                continue
            
            insertion_pos = self._get_valid_insertion_position(
                read.query_length, insertion_seq.insertion_length
            )
            
            if insertion_pos is None:
                self.logger.warning(f"Skipping read {read.query_name}: no valid insertion position "
                                  f"(read_length={read.query_length}, insertion_length={insertion_seq.insertion_length}, "
                                  f"min_distance_required={self.min_distance_from_ends * 2})")
                skipped_infeasible_reads.append({
                    'base_read_name': read.query_name,
                    'read_length': read.query_length,
                    'insertion_length': insertion_seq.insertion_length,
                    'insertion_id': insertion_id
                })
                continue
            
            base_read_name = read.query_name
            modified_read_name = self._generate_modified_read_name(insertion_id, base_read_name)
            
            original_sequence = read.query_sequence
            modified_sequence = self._insert_sequence_into_read(
                original_sequence, insertion_seq.sequence, insertion_pos
            )
            
            insertion_record = InsertionRecord(
                base_read_name=base_read_name,
                modified_read_name=modified_read_name,
                original_chr=metadata.original_chr,
                original_pos=metadata.original_pos,
                insertion_id=insertion_id,
                insertion_type=insertion_seq.insertion_type,
                insertion_length=insertion_seq.insertion_length,
                insertion_pos=insertion_pos
            )
            
            insertion_records.append(insertion_record)
            
            seq_record = SeqRecord(
                Seq(modified_sequence),
                id=modified_read_name,
                description=f"Modified with {insertion_seq.insertion_type} insertion at pos {insertion_pos}"
            )
            modified_sequences.append(seq_record)
        
        total_attempted = len(reads_with_metadata)
        successful_insertions = len(insertion_records)
        success_rate = successful_insertions / total_attempted if total_attempted > 0 else 0
        
        self.logger.info(f"Insertion summary: {successful_insertions}/{total_attempted} successful "
                        f"({success_rate:.2%} success rate)")
        
        if skipped_missing_sequences:
            self.logger.warning(f"Skipped {len(skipped_missing_sequences)} pairs due to missing sequences: "
                              f"{', '.join(set(skipped_missing_sequences))}")
        
        if skipped_infeasible_reads:
            self.logger.warning(f"Skipped {len(skipped_infeasible_reads)} pairs due to insufficient space")
        
        skip_stats = {
            'total_attempted': total_attempted,
            'successful_insertions': successful_insertions,
            'success_rate': success_rate,
            'skipped_missing_sequences': len(skipped_missing_sequences),
            'skipped_infeasible_reads': len(skipped_infeasible_reads)
        }
        
        return insertion_records, modified_sequences, skip_stats
    
    def save_insertion_records(self, insertion_records: List[InsertionRecord], 
                             output_path: str) -> None:
        """
        Save insertion records to JSON file.
        
        Args:
            insertion_records: List of insertion records
            output_path: Path to output JSON file
        """
        records_data = [record.to_dict() for record in insertion_records]
        
        with open(output_path, 'w') as f:
            json.dump(records_data, f, indent=2)
        
        self.logger.info(f"Saved {len(insertion_records)} insertion records to {output_path}")
    
    def save_modified_sequences(self, modified_sequences: List[SeqRecord], 
                              output_path: str) -> None:
        """
        Save modified sequences to FASTA file.
        
        Args:
            modified_sequences: List of modified sequence records
            output_path: Path to output FASTA file
        """
        with open(output_path, 'w') as f:
            SeqIO.write(modified_sequences, f, "fasta")
        
        self.logger.info(f"Saved {len(modified_sequences)} modified sequences to {output_path}")
    
    def get_insertion_statistics(self, insertion_records: List[InsertionRecord]) -> Dict[str, Any]:
        """
        Calculate statistics from insertion records.
        
        Args:
            insertion_records: List of insertion records
            
        Returns:
            Dictionary containing insertion statistics
        """
        if not insertion_records:
            return {'total_insertions': 0}
        
        type_counts = {}
        position_stats = {}
        length_stats = {}
        
        for record in insertion_records:
            type_counts[record.insertion_type] = type_counts.get(record.insertion_type, 0) + 1
            
            if record.insertion_type not in position_stats:
                position_stats[record.insertion_type] = []
            position_stats[record.insertion_type].append(record.insertion_pos)
            
            if record.insertion_type not in length_stats:
                length_stats[record.insertion_type] = []
            length_stats[record.insertion_type].append(record.insertion_length)
        for insertion_type, positions in position_stats.items():
            position_stats[insertion_type] = {
                'count': len(positions),
                'min': min(positions),
                'max': max(positions),
                'mean': sum(positions) / len(positions)
            }
        
        for insertion_type, lengths in length_stats.items():
            length_stats[insertion_type] = {
                'count': len(lengths),
                'min': min(lengths),
                'max': max(lengths),
                'mean': sum(lengths) / len(lengths)
            }
        
        return {
            'total_insertions': len(insertion_records),
            'type_counts': type_counts,
            'position_statistics': position_stats,
            'length_statistics': length_stats
        }
    
