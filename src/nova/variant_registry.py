"""
Variant registry module for managing insertion sequences and their metadata.
"""

import json
import hashlib
import time
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, asdict
import logging


@dataclass
class InsertionSequence:
    """Represents an insertion sequence with metadata."""
    insertion_id: str
    sequence: str
    insertion_type: str
    insertion_length: int
    metadata: Dict[str, Any]
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'InsertionSequence':
        """Create from dictionary representation."""
        return cls(**data)


class VariantRegistry:
    """
    Registry for managing insertion sequences with unique IDs and metadata.
    """
    
    def __init__(self):
        """Initialize empty registry."""
        self.sequences: Dict[str, InsertionSequence] = {}
        self.logger = logging.getLogger(__name__)
        
    def _generate_unique_id(self, sequence: str, insertion_type: str) -> str:
        """
        Generate unique ID for insertion sequence.
        Should be robust to collisions, even if the sequence is the same.
        
        Args:
            sequence: The insertion sequence
            insertion_type: Type of insertion
            
        Returns:
            Unique identifier string
        """
        content = f"{sequence}_{insertion_type}_{time.time_ns()}"
        hash_obj = hashlib.sha256(content.encode())
        return f"{insertion_type}_{hash_obj.hexdigest()[:12]}"
    
    def add_sequence(self, sequence: str, insertion_type: str, 
                    metadata: Optional[Dict[str, Any]] = None) -> str:
        """
        Add a new insertion sequence to the registry.
        
        Args:
            sequence: The insertion sequence
            insertion_type: Type of insertion (e.g., 'random', 'simple', 'Alu')
            metadata: Additional metadata
            
        Returns:
            Unique insertion ID
        """
        if metadata is None:
            metadata = {}
        
        insertion_id = self._generate_unique_id(sequence, insertion_type)
        
        insertion_seq = InsertionSequence(
            insertion_id=insertion_id,
            sequence=sequence,
            insertion_type=insertion_type,
            insertion_length=len(sequence),
            metadata=metadata
        )
        
        self.sequences[insertion_id] = insertion_seq
        self.logger.debug(f"Added sequence {insertion_id} of type {insertion_type}")
        
        return insertion_id
    
    def get_sequence(self, insertion_id: str) -> Optional[InsertionSequence]:
        """
        Get insertion sequence by ID.
        
        Args:
            insertion_id: Unique identifier
            
        Returns:
            InsertionSequence object or None if not found
        """
        return self.sequences.get(insertion_id)
    
    def get_sequences_by_type(self, insertion_type: str) -> List[InsertionSequence]:
        """
        Get all sequences of a specific type.
        
        Args:
            insertion_type: Type of insertion
            
        Returns:
            List of InsertionSequence objects
        """
        return [seq for seq in self.sequences.values() 
                if seq.insertion_type == insertion_type]
    
    def list_sequences(self) -> List[InsertionSequence]:
        """
        Get all sequences in the registry.
        
        Returns:
            List of all InsertionSequence objects
        """
        return list(self.sequences.values())
    
    def get_statistics(self) -> Dict[str, Any]:
        """
        Get registry statistics.
        
        Returns:
            Dictionary with statistics
        """
        type_counts = {}
        length_stats = {}
        
        for seq in self.sequences.values():
            type_counts[seq.insertion_type] = type_counts.get(seq.insertion_type, 0) + 1
            
            if seq.insertion_type not in length_stats:
                length_stats[seq.insertion_type] = []
            length_stats[seq.insertion_type].append(seq.insertion_length)
        
        for insertion_type, lengths in length_stats.items():
            length_stats[insertion_type] = {
                'count': len(lengths),
                'min': min(lengths),
                'max': max(lengths),
                'mean': sum(lengths) / len(lengths)
            }
        
        return {
            'total_sequences': len(self.sequences),
            'type_counts': type_counts,
            'length_statistics': length_stats
        }
    
    def save_to_json(self, filepath: str) -> None:
        """
        Save registry to JSON file.
        
        Args:
            filepath: Path to output JSON file
        """
        registry_data = {
            'sequences': {seq_id: seq.to_dict() 
                         for seq_id, seq in self.sequences.items()},
            'statistics': self.get_statistics()
        }
        
        with open(filepath, 'w') as f:
            json.dump(registry_data, f, indent=2)
        
        self.logger.info(f"Saved registry to {filepath}")
    
    def load_from_json(self, filepath: str) -> None:
        """
        Load registry from JSON file.
        
        Args:
            filepath: Path to input JSON file
        """
        with open(filepath, 'r') as f:
            registry_data = json.load(f)
        
        self.sequences = {}
        for seq_id, seq_data in registry_data['sequences'].items():
            self.sequences[seq_id] = InsertionSequence.from_dict(seq_data)
        
        self.logger.info(f"Loaded {len(self.sequences)} sequences from {filepath}")
    
    def clear(self) -> None:
        """Clear all sequences from registry."""
        self.sequences.clear()
        self.logger.info("Cleared registry")
    
    def remove_sequence(self, insertion_id: str) -> bool:
        """
        Remove sequence from registry.
        
        Args:
            insertion_id: Unique identifier
            
        Returns:
            True if removed, False if not found
        """
        if insertion_id in self.sequences:
            del self.sequences[insertion_id]
            self.logger.debug(f"Removed sequence {insertion_id}")
            return True
        return False
    
    def __len__(self) -> int:
        """Return number of sequences in registry."""
        return len(self.sequences)
    
    def __contains__(self, insertion_id: str) -> bool:
        """Check if insertion ID exists in registry."""
        return insertion_id in self.sequences