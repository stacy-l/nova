"""
Variant generator module for creating insertion sequences.
"""

import random
import string
from typing import List, Dict, Any, Optional
from Bio import SeqIO
from Bio.Seq import Seq
import logging

from .variant_registry import VariantRegistry


class VariantGenerator:
    """
    Generate insertion sequences for simulation.
    """
    
    def __init__(self, registry: VariantRegistry, random_seed: Optional[int] = None):
        """
        Initialize VariantGenerator.
        
        Args:
            registry: VariantRegistry instance to store generated sequences
            random_seed: Optional seed for reproducible random generation
        """
        self.registry = registry
        self.logger = logging.getLogger(__name__)
        
        if random_seed is not None:
            random.seed(random_seed)
    
    def _normalize_to_list(self, config_item: Any) -> List[Any]:
        """
        Normalize configuration item to a list format.
        
        Args:
            config_item: Either a single config object or a list of config objects
            
        Returns:
            List of configuration objects
        """
        if isinstance(config_item, list):
            return config_item
        else:
            return [config_item]
    
    def generate_random_insertions(self, n: int, length: int, 
                                 gc_content: Optional[float] = None) -> List[str]:
        """
        Generate random insertion sequences.
        
        Args:
            n: Number of variants to generate
            length: Length of each insertion
            gc_content: Target GC content (0.0-1.0). If None, uses equal probabilities.
            
        Returns:
            List of insertion IDs
        """
        nucleotides = ['A', 'T', 'G', 'C']
        insertion_ids = []
        
        if gc_content is not None:
            g_prob = c_prob = gc_content / 2
            a_prob = t_prob = (1 - gc_content) / 2
            weights = [a_prob, t_prob, g_prob, c_prob]
        else:
            weights = [0.25, 0.25, 0.25, 0.25]
        
        for i in range(n):
            sequence = ''.join(random.choices(nucleotides, weights=weights, k=length))
            
            metadata = {
                'generation_method': 'random',
                'target_gc_content': gc_content,
                'actual_gc_content': (sequence.count('G') + sequence.count('C')) / len(sequence)
            }
            
            insertion_id = self.registry.add_sequence(sequence, 'random', metadata)
            insertion_ids.append(insertion_id)
        
        self.logger.info(f"Generated {n} random insertions of length {length}")
        return insertion_ids
    
    def generate_simple_repeat_insertions(self, n: int, repeat_unit: str, 
                                        units: int) -> List[str]:
        """
        Generate simple repeat insertion sequences.
        
        Args:
            n: Number of variants to generate
            repeat_unit: Repeat unit sequence (e.g., 'CAG')
            units: Number of repeat units
            
        Returns:
            List of insertion IDs
        """
        insertion_ids = []
        sequence = repeat_unit * units
        
        for i in range(n):
            metadata = {
                'generation_method': 'simple_repeat',
                'repeat_unit': repeat_unit,
                'repeat_units': units,
                'unit_length': len(repeat_unit)
            }
            
            insertion_id = self.registry.add_sequence(sequence, 'simple', metadata)
            insertion_ids.append(insertion_id)
        
        self.logger.info(f"Generated {n} simple repeat insertions: {repeat_unit} x {units}")
        return insertion_ids
    
    def generate_predefined_insertions(self, fasta_path: str, 
                                     sequence_counts: Dict[str, int]) -> List[str]:
        """
        Generate predefined insertion sequences from FASTA file.
        
        Args:
            fasta_path: Path to FASTA file containing predefined sequences
            sequence_counts: Dictionary mapping sequence names to counts
                           e.g., {'AluYa5': 30, 'AluYb8': 20}
        
        Returns:
            List of insertion IDs
        """
        insertion_ids = []
        
        sequences = {}
        try:
            for record in SeqIO.parse(fasta_path, "fasta"):
                sequences[record.id] = str(record.seq)
        except Exception as e:
            self.logger.error(f"Failed to read FASTA file {fasta_path}: {e}")
            raise
        
        for seq_name, count in sequence_counts.items():
            if seq_name not in sequences:
                self.logger.warning(f"Sequence '{seq_name}' not found in {fasta_path}")
                continue
            
            sequence = sequences[seq_name]
            
            for i in range(count):
                metadata = {
                    'generation_method': 'predefined',
                    'source_file': fasta_path,
                    'source_sequence_name': seq_name,
                    'copy_number': i + 1
                }
                
                insertion_id = self.registry.add_sequence(sequence, seq_name, metadata)
                insertion_ids.append(insertion_id)
        
        total_generated = sum(sequence_counts.values())
        self.logger.info(f"Generated {total_generated} predefined insertions from {fasta_path}")
        return insertion_ids
    
    def generate_from_config(self, config: Dict[str, Any]) -> List[str]:
        """
        Generate insertion sequences from configuration dictionary.
        
        Args:
            config: Variant configuration dictionary specifying the number of variants to generate for each type
                   Example (single entry format):
                   {
                       'random': {'n': 50, 'length': 100},
                       'simple': {'n': 50, 'repeat': 'CAG', 'units': 40}
                   }
                   
                   Example (multiple entries format):
                   {
                       'random': [
                           {'n': 50, 'length': 100, 'gc_content': 0.4},
                           {'n': 50, 'length': 200}
                       ],
                       'simple': [
                           {'n': 100, 'repeat': 'CAG', 'units': 20},
                           {'n': 100, 'repeat': 'GC', 'units': 10}
                       ]
                   }
        
        Returns:
            List of all generated insertion IDs
        """
        all_insertion_ids = []
        
        if 'random' in config:
            random_configs = self._normalize_to_list(config['random'])
            for random_config in random_configs:
                n = random_config['n']
                length = random_config['length']
                gc_content = random_config.get('gc_content')
                
                ids = self.generate_random_insertions(n, length, gc_content)
                all_insertion_ids.extend(ids)
        
        if 'simple' in config:
            simple_configs = self._normalize_to_list(config['simple'])
            for simple_config in simple_configs:
                n = simple_config['n']
                repeat_unit = simple_config['repeat']
                units = simple_config['units']
                
                ids = self.generate_simple_repeat_insertions(n, repeat_unit, units)
                all_insertion_ids.extend(ids)
        
        if 'predefined' in config:
            predefined_configs = self._normalize_to_list(config['predefined'])
            for pred_config in predefined_configs:
                for pred_type, type_config in pred_config.items():
                    fasta_path = type_config['fasta']
                    sequence_counts = type_config['spec']
                    
                    ids = self.generate_predefined_insertions(fasta_path, sequence_counts)
                    all_insertion_ids.extend(ids)
        
        self.logger.info(f"Generated total of {len(all_insertion_ids)} insertion sequences")
        return all_insertion_ids
    
    def validate_config(self, config: Dict[str, Any]) -> List[str]:
        """
        Validate configuration dictionary.
        
        Args:
            config: Configuration dictionary
            
        Returns:
            List of validation error messages (empty if valid)
        """
        errors = []
        
        if not isinstance(config, dict):
            errors.append("Config must be a dictionary")
            return errors
        
        if 'random' in config:
            random_configs = self._normalize_to_list(config['random'])
            for i, random_config in enumerate(random_configs):
                prefix = f"Random config {i+1}: " if len(random_configs) > 1 else "Random config: "
                
                if not isinstance(random_config, dict):
                    errors.append(f"{prefix}must be a dictionary")
                    continue
                
                if 'n' not in random_config or 'length' not in random_config:
                    errors.append(f"{prefix}missing required 'n' or 'length' fields")
                if random_config.get('n', 0) <= 0:
                    errors.append(f"{prefix}'n' must be positive")
                if random_config.get('length', 0) <= 0:
                    errors.append(f"{prefix}'length' must be positive")
                if 'gc_content' in random_config:
                    gc = random_config['gc_content']
                    if not 0 <= gc <= 1:
                        errors.append(f"{prefix}'gc_content' must be between 0 and 1")
        
        if 'simple' in config:
            simple_configs = self._normalize_to_list(config['simple'])
            for i, simple_config in enumerate(simple_configs):
                prefix = f"Simple config {i+1}: " if len(simple_configs) > 1 else "Simple config: "
                
                if not isinstance(simple_config, dict):
                    errors.append(f"{prefix}must be a dictionary")
                    continue
                
                if 'n' not in simple_config or 'repeat' not in simple_config or 'units' not in simple_config:
                    errors.append(f"{prefix}missing required 'n', 'repeat', or 'units' fields")
                if simple_config.get('n', 0) <= 0:
                    errors.append(f"{prefix}'n' must be positive")
                if simple_config.get('units', 0) <= 0:
                    errors.append(f"{prefix}'units' must be positive")
                if not simple_config.get('repeat', ''):
                    errors.append(f"{prefix}'repeat' cannot be empty")
        
        if 'predefined' in config:
            predefined_configs = self._normalize_to_list(config['predefined'])
            for i, pred_config in enumerate(predefined_configs):
                prefix = f"Predefined config {i+1}: " if len(predefined_configs) > 1 else "Predefined config: "
                
                if not isinstance(pred_config, dict):
                    errors.append(f"{prefix}must be a dictionary")
                    continue
                
                for pred_type, type_config in pred_config.items():
                    if 'fasta' not in type_config or 'spec' not in type_config:
                        errors.append(f"{prefix}'{pred_type}' missing 'fasta' or 'spec' fields")
                    if not isinstance(type_config.get('spec', {}), dict):
                        errors.append(f"{prefix}'{pred_type}' 'spec' must be a dictionary")
        
        return errors