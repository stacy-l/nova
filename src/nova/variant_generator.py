"""
Variant generator module for creating insertion sequences.
"""

import random
import string
import hashlib
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
from Bio import SeqIO
from Bio.Seq import Seq
import logging

from .variant_registry import VariantRegistry


@dataclass
class MutationRecord:
    """Record of a single mutation applied to a sequence."""
    position: int
    original_base: str
    mutated_base: str


@dataclass
class MutationMetadata:
    """Metadata about mutations applied to a sequence."""
    substitution_rate: float
    mutations_applied: List[MutationRecord]
    mutation_seed: int
    total_mutations: int


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
        self.random_seed = random_seed
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
    
    def _generate_mutation_seed(self, base_seed: Optional[int], sequence_index: int) -> int:
        """
        Generate a deterministic mutation seed for a specific sequence.
        
        Args:
            base_seed: Base random seed from generator
            sequence_index: Index of the sequence being mutated
            
        Returns:
            Deterministic mutation seed for this sequence
        """
        if base_seed is None:
            base_seed = 42  # Default fallback
        
        # Create deterministic but unique seed for each sequence
        seed_string = f"{base_seed}_{sequence_index}_mutation"
        seed_hash = hashlib.md5(seed_string.encode()).hexdigest()
        return int(seed_hash[:8], 16)  # Use first 8 hex chars as seed
    
    def _apply_mutations(self, sequence: str, mutation_config: Dict[str, Any], 
                        sequence_index: int) -> Tuple[str, MutationMetadata]:
        """
        Apply mutations to a sequence based on configuration.
        
        Args:
            sequence: Original sequence to mutate
            mutation_config: Dictionary containing mutation parameters
            sequence_index: Index of this sequence for deterministic seeding
            
        Returns:
            Tuple of (mutated_sequence, mutation_metadata)
        """
        if not mutation_config or 'substitution_rate' not in mutation_config:
            # No mutations configured
            empty_metadata = MutationMetadata(
                substitution_rate=0.0,
                mutations_applied=[],
                mutation_seed=0,
                total_mutations=0
            )
            return sequence, empty_metadata
        
        substitution_rate = mutation_config['substitution_rate']
        
        if substitution_rate <= 0:
            # No mutations to apply
            empty_metadata = MutationMetadata(
                substitution_rate=substitution_rate,
                mutations_applied=[],
                mutation_seed=0,
                total_mutations=0
            )
            return sequence, empty_metadata
        
        # Generate mutation seed and create local random state
        mutation_seed = self._generate_mutation_seed(self.random_seed, sequence_index)
        mutation_random = random.Random(mutation_seed)
        
        # Calculate number of positions to mutate
        num_mutations = max(1, int(len(sequence) * substitution_rate))
        
        if num_mutations >= len(sequence):
            num_mutations = len(sequence) - 1  # Don't mutate entire sequence
        
        # Select random positions to mutate (without replacement)
        positions_to_mutate = mutation_random.sample(range(len(sequence)), num_mutations)
        
        # Convert sequence to list for easy mutation
        mutated_sequence = list(sequence)
        mutations_applied = []
        
        # Define nucleotide alternatives
        nucleotide_alternatives = {
            'A': ['T', 'G', 'C'],
            'T': ['A', 'G', 'C'],
            'G': ['A', 'T', 'C'],
            'C': ['A', 'T', 'G']
        }
        
        # Apply mutations
        for pos in positions_to_mutate:
            original_base = sequence[pos].upper()
            
            if original_base in nucleotide_alternatives:
                # Choose random alternative base
                alternatives = nucleotide_alternatives[original_base]
                mutated_base = mutation_random.choice(alternatives)
                
                # Apply mutation
                mutated_sequence[pos] = mutated_base
                
                # Record mutation
                mutation_record = MutationRecord(
                    position=pos,
                    original_base=original_base,
                    mutated_base=mutated_base
                )
                mutations_applied.append(mutation_record)
        
        # Create metadata
        metadata = MutationMetadata(
            substitution_rate=substitution_rate,
            mutations_applied=mutations_applied,
            mutation_seed=mutation_seed,
            total_mutations=len(mutations_applied)
        )
        
        return ''.join(mutated_sequence), metadata
    
    def _validate_mutation_config(self, mutation_config: Dict[str, Any], prefix: str) -> List[str]:
        """
        Validate mutation configuration parameters.
        
        Args:
            mutation_config: Dictionary containing mutation parameters
            prefix: Prefix for error messages
            
        Returns:
            List of validation error messages
        """
        errors = []
        
        if not isinstance(mutation_config, dict):
            errors.append(f"{prefix}mutations config must be a dictionary")
            return errors
        
        if 'substitution_rate' in mutation_config:
            sub_rate = mutation_config['substitution_rate']
            if not isinstance(sub_rate, (int, float)):
                errors.append(f"{prefix}substitution_rate must be a number")
            elif not 0 <= sub_rate <= 1:
                errors.append(f"{prefix}substitution_rate must be between 0 and 1")
        
        # Add validation for future mutation types here
        
        return errors
    
    def generate_random_insertions(self, n: int, length: int, 
                                 gc_content: Optional[float] = None,
                                 mutation_config: Optional[Dict[str, Any]] = None) -> List[str]:
        """
        Generate random insertion sequences.
        
        Args:
            n: Number of variants to generate
            length: Length of each insertion
            gc_content: Target GC content (0.0-1.0). If None, uses equal probabilities.
            mutation_config: Optional mutation configuration
            
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
            
            # Apply mutations if configured
            if mutation_config:
                sequence, mutation_metadata = self._apply_mutations(sequence, mutation_config, i)
            else:
                mutation_metadata = None
            
            metadata = {
                'generation_method': 'random',
                'target_gc_content': gc_content,
                'actual_gc_content': (sequence.count('G') + sequence.count('C')) / len(sequence)
            }
            
            # Add mutation metadata if mutations were applied
            if mutation_metadata:
                metadata['mutations'] = {
                    'substitution_rate': mutation_metadata.substitution_rate,
                    'total_mutations': mutation_metadata.total_mutations,
                    'mutation_seed': mutation_metadata.mutation_seed,
                    'mutation_records': [
                        {
                            'position': record.position,
                            'original_base': record.original_base,
                            'mutated_base': record.mutated_base
                        }
                        for record in mutation_metadata.mutations_applied
                    ]
                }
            
            insertion_id = self.registry.add_sequence(sequence, 'random', metadata)
            insertion_ids.append(insertion_id)
        
        mutation_text = f" with mutations" if mutation_config else ""
        self.logger.info(f"Generated {n} random insertions of length {length}{mutation_text}")
        return insertion_ids
    
    def generate_simple_repeat_insertions(self, n: int, repeat_unit: str, 
                                        units: int,
                                        mutation_config: Optional[Dict[str, Any]] = None) -> List[str]:
        """
        Generate simple repeat insertion sequences.
        
        Args:
            n: Number of variants to generate
            repeat_unit: Repeat unit sequence (e.g., 'CAG')
            units: Number of repeat units
            mutation_config: Optional mutation configuration
            
        Returns:
            List of insertion IDs
        """
        insertion_ids = []
        base_sequence = repeat_unit * units
        
        for i in range(n):
            sequence = base_sequence
            
            # Apply mutations if configured
            if mutation_config:
                sequence, mutation_metadata = self._apply_mutations(sequence, mutation_config, i)
            else:
                mutation_metadata = None
            
            metadata = {
                'generation_method': 'simple_repeat',
                'repeat_unit': repeat_unit,
                'repeat_units': units,
                'unit_length': len(repeat_unit)
            }
            
            # Add mutation metadata if mutations were applied
            if mutation_metadata:
                metadata['mutations'] = {
                    'substitution_rate': mutation_metadata.substitution_rate,
                    'total_mutations': mutation_metadata.total_mutations,
                    'mutation_seed': mutation_metadata.mutation_seed,
                    'mutation_records': [
                        {
                            'position': record.position,
                            'original_base': record.original_base,
                            'mutated_base': record.mutated_base
                        }
                        for record in mutation_metadata.mutations_applied
                    ]
                }
            
            insertion_id = self.registry.add_sequence(sequence, 'simple', metadata)
            insertion_ids.append(insertion_id)
        
        mutation_text = f" with mutations" if mutation_config else ""
        self.logger.info(f"Generated {n} simple repeat insertions: {repeat_unit} x {units}{mutation_text}")
        return insertion_ids
    
    def generate_predefined_insertions(self, fasta_path: str, 
                                     sequence_counts: Dict[str, int],
                                     mutation_config: Optional[Dict[str, Any]] = None) -> List[str]:
        """
        Generate predefined insertion sequences from FASTA file.
        
        Args:
            fasta_path: Path to FASTA file containing predefined sequences
            sequence_counts: Dictionary mapping sequence names to counts
                           e.g., {'AluYa5': 30, 'AluYb8': 20}
            mutation_config: Optional mutation configuration
        
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
        
        sequence_index = 0  # Global index for mutation seeding
        
        for seq_name, count in sequence_counts.items():
            if seq_name not in sequences:
                self.logger.warning(f"Sequence '{seq_name}' not found in {fasta_path}")
                continue
            
            base_sequence = sequences[seq_name]
            
            for i in range(count):
                sequence = base_sequence
                
                # Apply mutations if configured
                if mutation_config:
                    sequence, mutation_metadata = self._apply_mutations(sequence, mutation_config, sequence_index)
                else:
                    mutation_metadata = None
                
                metadata = {
                    'generation_method': 'predefined',
                    'source_file': fasta_path,
                    'source_sequence_name': seq_name,
                    'copy_number': i + 1
                }
                
                # Add mutation metadata if mutations were applied
                if mutation_metadata:
                    metadata['mutations'] = {
                        'substitution_rate': mutation_metadata.substitution_rate,
                        'total_mutations': mutation_metadata.total_mutations,
                        'mutation_seed': mutation_metadata.mutation_seed,
                        'mutation_records': [
                            {
                                'position': record.position,
                                'original_base': record.original_base,
                                'mutated_base': record.mutated_base
                            }
                            for record in mutation_metadata.mutations_applied
                        ]
                    }
                
                insertion_id = self.registry.add_sequence(sequence, seq_name, metadata)
                insertion_ids.append(insertion_id)
                sequence_index += 1  # Increment for next sequence
        
        total_generated = sum(sequence_counts.values())
        mutation_text = f" with mutations" if mutation_config else ""
        self.logger.info(f"Generated {total_generated} predefined insertions from {fasta_path}{mutation_text}")
        return insertion_ids
    
    def generate_from_config_with_regions(self, config: Dict[str, Any]) -> Dict[str, List[str]]:
        """
        Generate insertion sequences grouped by their target region requirements.
        
        Args:
            config: Variant configuration dictionary
            
        Returns:
            Dictionary mapping target BED file paths (or 'global' for no targeting) to lists of insertion IDs
        """
        region_groups = {}
        
        if 'random' in config:
            random_configs = self._normalize_to_list(config['random'])
            for random_config in random_configs:
                n = random_config['n']
                length = random_config['length']
                gc_content = random_config.get('gc_content')
                mutation_config = random_config.get('mutations')
                target_regions = random_config.get('target_regions', 'global')
                
                ids = self.generate_random_insertions(n, length, gc_content, mutation_config)
                
                if target_regions not in region_groups:
                    region_groups[target_regions] = []
                region_groups[target_regions].extend(ids)
        
        if 'simple' in config:
            simple_configs = self._normalize_to_list(config['simple'])
            for simple_config in simple_configs:
                n = simple_config['n']
                repeat_unit = simple_config['repeat']
                units = simple_config['units']
                mutation_config = simple_config.get('mutations')
                target_regions = simple_config.get('target_regions', 'global')
                
                ids = self.generate_simple_repeat_insertions(n, repeat_unit, units, mutation_config)
                
                if target_regions not in region_groups:
                    region_groups[target_regions] = []
                region_groups[target_regions].extend(ids)
        
        if 'predefined' in config:
            predefined_configs = self._normalize_to_list(config['predefined'])
            for pred_config in predefined_configs:
                for pred_type, type_config in pred_config.items():
                    fasta_path = type_config['fasta']
                    sequence_counts = type_config['spec']
                    mutation_config = type_config.get('mutations')
                    target_regions = type_config.get('target_regions', 'global')
                    
                    ids = self.generate_predefined_insertions(fasta_path, sequence_counts, mutation_config)
                    
                    if target_regions not in region_groups:
                        region_groups[target_regions] = []
                    region_groups[target_regions].extend(ids)
        
        total_generated = sum(len(ids) for ids in region_groups.values())
        self.logger.info(f"Generated total of {total_generated} insertion sequences across {len(region_groups)} region groups")
        
        return region_groups

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
                mutation_config = random_config.get('mutations')
                
                ids = self.generate_random_insertions(n, length, gc_content, mutation_config)
                all_insertion_ids.extend(ids)
        
        if 'simple' in config:
            simple_configs = self._normalize_to_list(config['simple'])
            for simple_config in simple_configs:
                n = simple_config['n']
                repeat_unit = simple_config['repeat']
                units = simple_config['units']
                mutation_config = simple_config.get('mutations')
                
                ids = self.generate_simple_repeat_insertions(n, repeat_unit, units, mutation_config)
                all_insertion_ids.extend(ids)
        
        if 'predefined' in config:
            predefined_configs = self._normalize_to_list(config['predefined'])
            for pred_config in predefined_configs:
                for pred_type, type_config in pred_config.items():
                    fasta_path = type_config['fasta']
                    sequence_counts = type_config['spec']
                    mutation_config = type_config.get('mutations')
                    
                    ids = self.generate_predefined_insertions(fasta_path, sequence_counts, mutation_config)
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
                
                # Validate mutations config if present
                if 'mutations' in random_config:
                    mutation_errors = self._validate_mutation_config(random_config['mutations'], prefix)
                    errors.extend(mutation_errors)
                
                # Validate target_regions if present
                if 'target_regions' in random_config:
                    target_bed = random_config['target_regions']
                    if not isinstance(target_bed, str) or not target_bed.strip():
                        errors.append(f"{prefix}'target_regions' must be a non-empty string path to BED file")
        
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
                
                # Validate mutations config if present
                if 'mutations' in simple_config:
                    mutation_errors = self._validate_mutation_config(simple_config['mutations'], prefix)
                    errors.extend(mutation_errors)
                
                # Validate target_regions if present
                if 'target_regions' in simple_config:
                    target_bed = simple_config['target_regions']
                    if not isinstance(target_bed, str) or not target_bed.strip():
                        errors.append(f"{prefix}'target_regions' must be a non-empty string path to BED file")
        
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
                    for seq_name, count in type_config['spec'].items():
                        if not isinstance(count, int) or count <= 0:
                            errors.append(f"{prefix}'{pred_type}' 'spec' must contain positive integer counts")
                    
                    # Validate mutations config if present
                    if 'mutations' in type_config:
                        mutation_errors = self._validate_mutation_config(type_config['mutations'], f"{prefix}'{pred_type}' ")
                        errors.extend(mutation_errors)
                    
                    # Validate target_regions if present
                    if 'target_regions' in type_config:
                        target_bed = type_config['target_regions']
                        if not isinstance(target_bed, str) or not target_bed.strip():
                            errors.append(f"{prefix}'{pred_type}' 'target_regions' must be a non-empty string path to BED file")
        
        return errors