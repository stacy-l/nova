"""
Sequence comparison and analysis utilities for Nova variant simulation.

This module provides utilities for comparing FASTA files and analyzing sequence content,
useful for validating simulation outputs, comparing methods, and analyzing variant 
calling results.
"""

import hashlib
from typing import Set, Dict, List, Tuple, Any, Optional
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import logging


def extract_sequence_content_hashes(fasta_path: str) -> Set[str]:
    """
    Extract SHA256 hashes of sequence content from FASTA file.
    
    This function reads a FASTA file and creates a set of hashes based solely on
    the sequence content, ignoring sequence IDs and descriptions. Useful for
    comparing sequence sets regardless of naming differences.
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        Set of SHA256 hash strings representing unique sequences
        
    Raises:
        FileNotFoundError: If FASTA file doesn't exist
        ValueError: If file cannot be parsed as FASTA
    """
    if not Path(fasta_path).exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    
    sequence_hashes = set()
    
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            # Hash just the sequence content (uppercase for consistency)
            sequence_str = str(record.seq).upper()
            sequence_hash = hashlib.sha256(sequence_str.encode()).hexdigest()
            sequence_hashes.add(sequence_hash)
    except Exception as e:
        raise ValueError(f"Error parsing FASTA file {fasta_path}: {e}")
    
    return sequence_hashes


def compare_fasta_content(fasta1_path: str, fasta2_path: str) -> Dict[str, Any]:
    """
    Compare sequence content between two FASTA files.
    
    Performs comprehensive comparison of sequence content between two FASTA files,
    identifying common sequences, unique sequences, and providing summary statistics.
    
    Args:
        fasta1_path: Path to first FASTA file
        fasta2_path: Path to second FASTA file
        
    Returns:
        Dictionary containing comparison results:
        - sequences_identical: bool, True if both files contain identical sequence sets
        - file1_sequences: int, number of unique sequences in file 1
        - file2_sequences: int, number of unique sequences in file 2
        - common_sequences: int, number of sequences present in both files
        - unique_to_file1: int, number of sequences only in file 1
        - unique_to_file2: int, number of sequences only in file 2
        - overlap_percentage: float, percentage of sequences that are common
    """
    # Extract sequence hashes from both files
    hashes1 = extract_sequence_content_hashes(fasta1_path)
    hashes2 = extract_sequence_content_hashes(fasta2_path)
    
    # Calculate set relationships
    common = hashes1 & hashes2
    unique_to_1 = hashes1 - hashes2
    unique_to_2 = hashes2 - hashes1
    
    # Calculate overlap percentage
    total_unique_sequences = len(hashes1 | hashes2)
    overlap_percentage = (len(common) / total_unique_sequences * 100) if total_unique_sequences > 0 else 0.0
    
    return {
        'sequences_identical': hashes1 == hashes2,
        'file1_sequences': len(hashes1),
        'file2_sequences': len(hashes2),
        'common_sequences': len(common),
        'unique_to_file1': len(unique_to_1),
        'unique_to_file2': len(unique_to_2),
        'overlap_percentage': overlap_percentage,
        'file1_path': fasta1_path,
        'file2_path': fasta2_path
    }


def analyze_sequence_recovery(ground_truth_fasta: str, detected_fasta: str) -> Dict[str, Any]:
    """
    Analyze sequence recovery between ground truth and detected sequences.
    
    Designed for variant calling validation: compares a ground truth set of 
    sequences (e.g., Nova simulation output) against sequences detected by
    variant calling pipelines.
    
    Args:
        ground_truth_fasta: Path to FASTA file containing true sequences
        detected_fasta: Path to FASTA file containing detected sequences
        
    Returns:
        Dictionary containing recovery analysis:
        - total_ground_truth: int, number of ground truth sequences
        - total_detected: int, number of detected sequences  
        - true_positives: int, ground truth sequences that were detected
        - false_negatives: int, ground truth sequences not detected
        - false_positives: int, detected sequences not in ground truth
        - sensitivity: float, recall rate (TP / (TP + FN))
        - precision: float, precision rate (TP / (TP + FP))
        - f1_score: float, harmonic mean of precision and recall
    """
    ground_truth_hashes = extract_sequence_content_hashes(ground_truth_fasta)
    detected_hashes = extract_sequence_content_hashes(detected_fasta)
    
    true_positives = len(ground_truth_hashes & detected_hashes)
    false_negatives = len(ground_truth_hashes - detected_hashes)
    false_positives = len(detected_hashes - ground_truth_hashes)
    
    total_ground_truth = len(ground_truth_hashes)
    total_detected = len(detected_hashes)
    
    sensitivity = true_positives / total_ground_truth if total_ground_truth > 0 else 0.0
    precision = true_positives / total_detected if total_detected > 0 else 0.0
    f1_score = (2 * precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) > 0 else 0.0
    
    return {
        'total_ground_truth': total_ground_truth,
        'total_detected': total_detected,
        'true_positives': true_positives,
        'false_negatives': false_negatives,
        'false_positives': false_positives,
        'sensitivity': sensitivity,
        'precision': precision,
        'f1_score': f1_score,
        'ground_truth_path': ground_truth_fasta,
        'detected_path': detected_fasta
    }


def compare_multiple_fasta_files(fasta_paths: List[str], labels: Optional[List[str]] = None) -> Dict[str, Any]:
    """
    Compare sequence content across multiple FASTA files.
    
    Useful for comparing multiple variant calling methods or simulation runs
    against ground truth data.
    
    Args:
        fasta_paths: List of paths to FASTA files to compare
        labels: Optional list of labels for each file (defaults to filenames)
        
    Returns:
        Dictionary containing multi-way comparison results:
        - file_info: List of dicts with file paths, labels, and sequence counts
        - pairwise_comparisons: Dict of pairwise comparison results
        - common_to_all: int, sequences present in all files
        - unique_sequences: Dict mapping each file to its unique sequence count
        - overlap_matrix: 2D list showing pairwise overlap percentages
    """
    if not fasta_paths:
        raise ValueError("At least one FASTA file path must be provided")
    
    if labels is None:
        labels = [Path(path).stem for path in fasta_paths]
    
    if len(labels) != len(fasta_paths):
        raise ValueError("Number of labels must match number of FASTA files")
    
    all_hashes = []
    file_info = []
    
    for i, fasta_path in enumerate(fasta_paths):
        hashes = extract_sequence_content_hashes(fasta_path)
        all_hashes.append(hashes)
        file_info.append({
            'path': fasta_path,
            'label': labels[i],
            'sequence_count': len(hashes)
        })
    
    common_to_all = set.intersection(*all_hashes) if all_hashes else set()
    
    unique_sequences = {}
    for i, (label, hashes) in enumerate(zip(labels, all_hashes)):
        other_hashes = set.union(*(all_hashes[:i] + all_hashes[i+1:]))
        unique_sequences[label] = len(hashes - other_hashes)
    
    pairwise_comparisons = {}
    overlap_matrix = [[0.0 for _ in range(len(fasta_paths))] for _ in range(len(fasta_paths))]
    
    for i in range(len(fasta_paths)):
        for j in range(i, len(fasta_paths)):
            if i == j:
                overlap_matrix[i][j] = 100.0  # File overlaps 100% with itself
            else:
                comparison = compare_fasta_content(fasta_paths[i], fasta_paths[j])
                key = f"{labels[i]}_vs_{labels[j]}"
                pairwise_comparisons[key] = comparison
                overlap_matrix[i][j] = comparison['overlap_percentage']
                overlap_matrix[j][i] = comparison['overlap_percentage']  # Symmetric
    
    return {
        'file_info': file_info,
        'pairwise_comparisons': pairwise_comparisons,
        'common_to_all': len(common_to_all),
        'unique_sequences': unique_sequences,
        'overlap_matrix': overlap_matrix,
        'labels': labels
    }


def get_sequence_statistics(fasta_path: str) -> Dict[str, Any]:
    """
    Get basic statistics about sequences in a FASTA file.
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        Dictionary containing sequence statistics:
        - total_sequences: int, total number of sequence records
        - unique_sequences: int, number of unique sequences (by content)
        - duplicate_sequences: int, number of duplicate sequences
        - min_length: int, minimum sequence length
        - max_length: int, maximum sequence length
        - mean_length: float, average sequence length
        - total_bases: int, total number of bases across all sequences
    """
    if not Path(fasta_path).exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    
    lengths = []
    seen_sequences = set()
    total_sequences = 0
    duplicates = 0
    
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            total_sequences += 1
            sequence_str = str(record.seq).upper()
            lengths.append(len(sequence_str))
            
            sequence_hash = hashlib.sha256(sequence_str.encode()).hexdigest()
            if sequence_hash in seen_sequences:
                duplicates += 1
            else:
                seen_sequences.add(sequence_hash)
    
    except Exception as e:
        raise ValueError(f"Error parsing FASTA file {fasta_path}: {e}")
    
    if not lengths:
        return {
            'total_sequences': 0,
            'unique_sequences': 0,
            'duplicate_sequences': 0,
            'min_length': 0,
            'max_length': 0,
            'mean_length': 0.0,
            'total_bases': 0
        }
    
    return {
        'total_sequences': total_sequences,
        'unique_sequences': len(seen_sequences),
        'duplicate_sequences': duplicates,
        'min_length': min(lengths),
        'max_length': max(lengths),
        'mean_length': sum(lengths) / len(lengths),
        'total_bases': sum(lengths)
    }