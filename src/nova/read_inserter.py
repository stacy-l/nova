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
from .read_selector import ReadMetadata, LazyReadReference


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
        return f"{insertion_id}.{base_read_name}"
    
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
    
    def _stream_insertion(self, lazy_reads: List[LazyReadReference], 
                            insertion_ids: List[str]):
        """
        Insert sequences into reads using streaming to minimize memory usage.
        
        Args:
            lazy_reads: List of LazyReadReference objects
            insertion_ids: List of insertion IDs to use
            
        Yields:
            Tuple of (InsertionRecord, SeqRecord) for each successful insertion
        """
        if len(lazy_reads) != len(insertion_ids):
            raise ValueError(f"Number of reads ({len(lazy_reads)}) must match "
                           f"number of insertions ({len(insertion_ids)})")
        
        skipped_missing_sequences = []
        skipped_infeasible_reads = []
        successful_insertions = 0
        total_attempted = len(lazy_reads)
        
        for lazy_read, insertion_id in zip(lazy_reads, insertion_ids):
            insertion_seq = self.registry.get_sequence(insertion_id)
            if insertion_seq is None:
                self.logger.error(f"Insertion sequence {insertion_id} not found in registry")
                skipped_missing_sequences.append(insertion_id)
                continue
            
            insertion_pos = self._get_valid_insertion_position(
                lazy_read.read_length, insertion_seq.insertion_length
            )
            
            if insertion_pos is None:
                self.logger.warning(f"Skipping read {lazy_read.read_name}: no valid insertion position "
                                  f"(read_length={lazy_read.read_length}, insertion_length={insertion_seq.insertion_length}, "
                                  f"min_distance_required={self.min_distance_from_ends * 2})")
                skipped_infeasible_reads.append({
                    'base_read_name': lazy_read.read_name,
                    'read_length': lazy_read.read_length,
                    'insertion_length': insertion_seq.insertion_length,
                    'insertion_id': insertion_id
                })
                continue
            
            # Fetch the sequence on-demand
            original_sequence = lazy_read.get_sequence()
            if original_sequence is None:
                self.logger.warning(f"Could not fetch sequence for read {lazy_read.read_name}")
                skipped_infeasible_reads.append({
                    'base_read_name': lazy_read.read_name,
                    'read_length': lazy_read.read_length,
                    'insertion_length': insertion_seq.insertion_length,
                    'insertion_id': insertion_id,
                    'reason': 'sequence_fetch_failed'
                })
                continue
            
            base_read_name = lazy_read.read_name
            modified_read_name = self._generate_modified_read_name(insertion_id, base_read_name)
            
            modified_sequence = self._insert_sequence_into_read(
                original_sequence, insertion_seq.sequence, insertion_pos
            )
            
            insertion_record = InsertionRecord(
                base_read_name=base_read_name,
                modified_read_name=modified_read_name,
                original_chr=lazy_read.original_chr,
                original_pos=lazy_read.original_pos,
                insertion_id=insertion_id,
                insertion_type=insertion_seq.insertion_type,
                insertion_length=insertion_seq.insertion_length,
                insertion_pos=insertion_pos
            )
            
            seq_record = SeqRecord(
                Seq(modified_sequence),
                id=modified_read_name,
                description=""
            )
            
            successful_insertions += 1
            yield insertion_record, seq_record
        
        # Log final statistics
        success_rate = successful_insertions / total_attempted if total_attempted > 0 else 0
        
        self.logger.info(f"Insertion summary: {successful_insertions}/{total_attempted} successful "
                        f"({success_rate:.2%} success rate)")
        
        if skipped_missing_sequences:
            self.logger.warning(f"Skipped {len(skipped_missing_sequences)} pairs due to missing sequences: "
                              f"{', '.join(set(skipped_missing_sequences))}")
        
        if skipped_infeasible_reads:
            self.logger.warning(f"Skipped {len(skipped_infeasible_reads)} pairs due to insufficient space or fetch failures")
    
    def insert_streaming(self, lazy_reads: List[LazyReadReference], 
                             insertion_ids: List[str],
                             records_path: str,
                             sequences_path: str) -> Dict[str, Any]:
        """
        Process insertions in streaming mode and save results directly to files.
        
        Args:
            lazy_reads: List of LazyReadReference objects
            insertion_ids: List of insertion IDs to use
            records_path: Path to output JSON file
            sequences_path: Path to output FASTA file
            
        Returns:
            Dictionary with statistics and file paths
        """
        insertion_records = []
        sequences_written = 0
        
        # Open files for streaming write
        with open(records_path, 'w') as records_f, open(sequences_path, 'w') as sequences_f:
            records_f.write('[\n')  # Start JSON array
            first_record = True
            
            for insertion_record, seq_record in self._stream_insertion(lazy_reads, insertion_ids):
                # Write insertion record to JSON (streaming)
                if not first_record:
                    records_f.write(',\n')
                else:
                    first_record = False
                
                json.dump(insertion_record.to_dict(), records_f, indent=2)
                insertion_records.append(insertion_record)  # Keep for statistics
                
                # Write sequence record to FASTA (streaming)
                SeqIO.write([seq_record], sequences_f, "fasta")
                sequences_written += 1
            
            records_f.write('\n]')  # Close JSON array
        
        if insertion_records:
            insertion_stats = self.get_insertion_statistics(insertion_records)
        
        self.logger.info(f"Results saved: {len(insertion_records)} records to {records_path}, "
                        f"{sequences_written} sequences to {sequences_path}")
        
        return insertion_stats
    
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
        length_stats = {}
        
        for record in insertion_records:
            type_counts[record.insertion_type] = type_counts.get(record.insertion_type, 0) + 1
            
            if record.insertion_type not in length_stats:
                length_stats[record.insertion_type] = []
            length_stats[record.insertion_type].append(record.insertion_length)
        
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
            'length_statistics': length_stats
        }
    
