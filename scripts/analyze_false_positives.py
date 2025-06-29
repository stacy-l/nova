#!/usr/bin/env python3
"""
Analyze false positive variants in nova simulation results.

This script identifies and categorizes variants that contain nova reads but are not
single-read calls, which represent false positives in the de novo simulation.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
from collections import defaultdict, Counter
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys

def load_data(output_dir, output_prefix="nova"):
    """Load variant data and insertion metadata."""
    output_path = Path(output_dir)
    
    # Load tabular variant data
    tabular_file = output_path / f"{output_prefix}_variant_analysis_tabular.csv"
    if not tabular_file.exists():
        raise FileNotFoundError(f"Tabular data file not found: {tabular_file}")
    
    variants_df = pd.read_csv(tabular_file)
    print(f"Loaded {len(variants_df)} variant records from tabular data")
    
    # Load insertion metadata
    insertions_file = output_path / f"{output_prefix}_insertions.json"
    if not insertions_file.exists():
        raise FileNotFoundError(f"Insertions file not found: {insertions_file}")
        
    with open(insertions_file, 'r') as f:
        insertions_data = json.load(f)
    
    print(f"Loaded {len(insertions_data)} insertion records")
    
    return variants_df, insertions_data

def create_insertion_lookup(insertions_data):
    """Create lookup dictionaries for insertion metadata."""
    read_to_insertion = {}
    read_to_sequence = {}
    
    for insertion in insertions_data:
        read_name = insertion.get('modified_read_name', '')
        if read_name:
            read_to_insertion[read_name] = {
                'insertion_type': insertion.get('insertion_type', 'unknown'),
                'insertion_id': insertion.get('insertion_id', ''),
                'insertion_length': insertion.get('insertion_length', 0),
                'original_chr': insertion.get('original_chr', ''),
                'original_pos': insertion.get('original_pos', 0),
                'base_read_name': insertion.get('base_read_name', '')
            }
            
            # Try to get sequence if available (may be in different format)
            sequence = None
            if 'insertion_sequence' in insertion:
                seq_data = insertion['insertion_sequence']
                if isinstance(seq_data, dict):
                    sequence = seq_data.get('sequence', '')
                elif isinstance(seq_data, str):
                    sequence = seq_data
            
            if sequence:
                read_to_sequence[read_name] = sequence
    
    print(f"Created lookup for {len(read_to_insertion)} insertions")
    return read_to_insertion, read_to_sequence

def categorize_false_positives(variants_df):
    """Categorize variants into false positive types."""
    
    # Group by variant_index to get unique variants
    variant_groups = variants_df.groupby('variant_index').agg({
        'support_reads': 'first',
        'nova_reads': 'first', 
        'is_single_read_call': 'first',
        'nova_read_name': lambda x: list(x.dropna()),
        'insertion_type': lambda x: list(x.dropna()),
        'svtype': 'first',
        'precision_category': 'first',
        'chrom': 'first',
        'pos': 'first',
        'svlen': 'first'
    }).reset_index()
    
    # Categorize variants
    successful_calls = variant_groups[variant_groups['is_single_read_call'] == True]
    false_positives = variant_groups[variant_groups['is_single_read_call'] == False]
    
    print(f"\nVariant Classification:")
    print(f"  Successful calls (single nova read): {len(successful_calls)}")
    print(f"  False positives: {len(false_positives)}")
    
    # Categorize false positives
    fp_categories = {}
    
    for _, variant in false_positives.iterrows():
        support_reads = variant['support_reads']
        nova_reads = variant['nova_reads']
        non_nova_reads = support_reads - nova_reads
        
        variant_idx = variant['variant_index']
        
        if nova_reads > 1 and non_nova_reads == 0:
            category = "multiple_nova_only"
        elif nova_reads >= 1 and non_nova_reads >= 1:
            if nova_reads >= non_nova_reads:
                category = "mixed_nova_majority"
            else:
                category = "mixed_nova_minority"
        else:
            category = "other"  # Should not happen given our filtering
        
        fp_categories[variant_idx] = {
            'category': category,
            'support_reads': support_reads,
            'nova_reads': nova_reads,
            'non_nova_reads': non_nova_reads,
            'nova_read_names': variant['nova_read_name'],
            'insertion_types': variant['insertion_type'],
            'svtype': variant['svtype'],
            'precision': variant['precision_category'],
            'chrom': variant['chrom'],
            'pos': variant['pos'],
            'svlen': variant['svlen']
        }
    
    # Print category counts
    category_counts = Counter([fp['category'] for fp in fp_categories.values()])
    print(f"\nFalse Positive Categories:")
    for category, count in category_counts.items():
        print(f"  {category}: {count}")
    
    return fp_categories, successful_calls

def analyze_sequence_patterns(fp_categories, read_to_insertion, read_to_sequence):
    """Analyze sequence patterns in false positive variants."""
    
    pattern_analysis = {}
    
    for variant_idx, fp_data in fp_categories.items():
        nova_reads = fp_data['nova_read_names']
        
        # Get insertion details for all nova reads
        insertion_details = []
        for read_name in nova_reads:
            if read_name in read_to_insertion:
                details = read_to_insertion[read_name].copy()
                details['read_name'] = read_name
                if read_name in read_to_sequence:
                    details['sequence'] = read_to_sequence[read_name]
                insertion_details.append(details)
        
        # Analyze patterns
        insertion_types = [d['insertion_type'] for d in insertion_details]
        insertion_ids = [d['insertion_id'] for d in insertion_details]
        original_positions = [(d['original_chr'], d['original_pos']) for d in insertion_details]
        
        # Check for duplicate insertion types
        type_counts = Counter(insertion_types)
        has_duplicate_types = any(count > 1 for count in type_counts.values())
        
        # Check for identical insertion IDs (same exact sequence)
        id_counts = Counter(insertion_ids)
        has_identical_sequences = any(count > 1 for count in id_counts.values())
        
        # Check for genomic proximity (reads from same region)
        position_clusters = defaultdict(list)
        for i, (chrom, pos) in enumerate(original_positions):
            position_clusters[chrom].append((pos, i))
        
        # Find reads within 5kb of each other
        nearby_reads = []
        for chrom, positions in position_clusters.items():
            positions.sort()  # Sort by position
            for i in range(len(positions)):
                for j in range(i + 1, len(positions)):
                    pos1, idx1 = positions[i]
                    pos2, idx2 = positions[j]
                    if abs(pos1 - pos2) <= 5000:  # Within 5kb
                        nearby_reads.append((idx1, idx2, abs(pos1 - pos2)))
        
        pattern_analysis[variant_idx] = {
            'insertion_details': insertion_details,
            'type_counts': dict(type_counts),
            'id_counts': dict(id_counts),
            'has_duplicate_types': has_duplicate_types,
            'has_identical_sequences': has_identical_sequences,
            'nearby_read_pairs': nearby_reads,
            'has_genomic_clustering': len(nearby_reads) > 0,
            'num_unique_types': len(type_counts),
            'num_unique_sequences': len(id_counts)
        }
    
    return pattern_analysis

def analyze_root_causes(fp_categories, pattern_analysis):
    """Identify potential root causes of false positives."""
    
    root_cause_summary = {
        'identical_sequences': 0,
        'similar_types_nearby': 0,
        'genomic_clustering': 0,
        'multiple_diverse_types': 0,
        'precision_related': 0
    }
    
    detailed_causes = {}
    
    for variant_idx, fp_data in fp_categories.items():
        patterns = pattern_analysis.get(variant_idx, {})
        
        causes = []
        
        # Check for identical sequences
        if patterns.get('has_identical_sequences', False):
            causes.append('identical_sequences')
            root_cause_summary['identical_sequences'] += 1
        
        # Check for similar types in nearby regions
        if patterns.get('has_duplicate_types', False) and patterns.get('has_genomic_clustering', False):
            causes.append('similar_types_nearby')
            root_cause_summary['similar_types_nearby'] += 1
        
        # Check for general genomic clustering
        if patterns.get('has_genomic_clustering', False):
            causes.append('genomic_clustering')
            root_cause_summary['genomic_clustering'] += 1
        
        # Check for multiple diverse types (potential mapping artifacts)
        if patterns.get('num_unique_types', 0) >= 3:
            causes.append('multiple_diverse_types')
            root_cause_summary['multiple_diverse_types'] += 1
        
        # Check if precision might be related
        if fp_data.get('precision', '') == 'IMPRECISE':
            causes.append('precision_related')
            root_cause_summary['precision_related'] += 1
        
        detailed_causes[variant_idx] = causes
    
    return root_cause_summary, detailed_causes

def generate_summary_statistics(fp_categories, pattern_analysis, root_cause_summary):
    """Generate summary statistics for the analysis."""
    
    total_fps = len(fp_categories)
    
    # Category breakdown
    category_stats = Counter([fp['category'] for fp in fp_categories.values()])
    
    # Insertion type involvement in FPs
    type_involvement = defaultdict(int)
    for fp_data in fp_categories.values():
        for ins_type in fp_data['insertion_types']:
            type_involvement[ins_type] += 1
    
    # Pattern statistics
    pattern_stats = {
        'variants_with_identical_sequences': sum(1 for p in pattern_analysis.values() if p.get('has_identical_sequences', False)),
        'variants_with_genomic_clustering': sum(1 for p in pattern_analysis.values() if p.get('has_genomic_clustering', False)),
        'variants_with_duplicate_types': sum(1 for p in pattern_analysis.values() if p.get('has_duplicate_types', False))
    }
    
    summary = {
        'total_false_positives': total_fps,
        'category_breakdown': dict(category_stats),
        'insertion_type_involvement': dict(type_involvement),
        'pattern_statistics': pattern_stats,
        'root_cause_summary': root_cause_summary
    }
    
    return summary

def save_results(fp_categories, pattern_analysis, detailed_causes, summary_stats, output_prefix="output/false_positives"):
    """Save analysis results to JSON and CSV formats."""
    
    # Prepare JSON output
    json_output = {
        'metadata': {
            'analysis_type': 'false_positive_analysis',
            'timestamp': pd.Timestamp.now().isoformat(),
            'total_false_positives': len(fp_categories)
        },
        'summary_statistics': summary_stats,
        'false_positive_details': fp_categories,
        'pattern_analysis': pattern_analysis,
        'root_causes': detailed_causes
    }
    
    # Save JSON
    json_file = f"{output_prefix}_analysis.json"
    with open(json_file, 'w') as f:
        json.dump(json_output, f, indent=2, default=str)
    print(f"JSON results saved to: {json_file}")
    
    # Prepare CSV output
    csv_records = []
    for variant_idx, fp_data in fp_categories.items():
        patterns = pattern_analysis.get(variant_idx, {})
        causes = detailed_causes.get(variant_idx, [])
        
        record = {
            'variant_index': variant_idx,
            'category': fp_data['category'],
            'support_reads': fp_data['support_reads'],
            'nova_reads': fp_data['nova_reads'],
            'non_nova_reads': fp_data['non_nova_reads'],
            'svtype': fp_data['svtype'],
            'precision': fp_data['precision'],
            'chrom': fp_data['chrom'],
            'pos': fp_data['pos'],
            'svlen': fp_data['svlen'],
            'insertion_types': ','.join(fp_data['insertion_types']),
            'num_unique_types': patterns.get('num_unique_types', 0),
            'num_unique_sequences': patterns.get('num_unique_sequences', 0),
            'has_identical_sequences': patterns.get('has_identical_sequences', False),
            'has_genomic_clustering': patterns.get('has_genomic_clustering', False),
            'has_duplicate_types': patterns.get('has_duplicate_types', False),
            'root_causes': ','.join(causes),
            'nova_read_names': ','.join(fp_data['nova_read_names'])
        }
        csv_records.append(record)
    
    # Save CSV
    csv_file = f"{output_prefix}_tabular.csv"
    csv_df = pd.DataFrame(csv_records)
    csv_df.to_csv(csv_file, index=False)
    print(f"CSV results saved to: {csv_file}")
    print(f"CSV shape: {csv_df.shape[0]} rows × {csv_df.shape[1]} columns")
    
    return json_file, csv_file

def print_analysis_summary(summary_stats):
    """Print a comprehensive analysis summary."""
    
    print("\n" + "="*60)
    print("FALSE POSITIVE ANALYSIS SUMMARY")
    print("="*60)
    
    print(f"\nTotal False Positives: {summary_stats['total_false_positives']}")
    
    print("\n1. FALSE POSITIVE CATEGORIES:")
    for category, count in summary_stats['category_breakdown'].items():
        percentage = (count / summary_stats['total_false_positives']) * 100
        print(f"   {category}: {count} ({percentage:.1f}%)")
    
    print("\n2. INSERTION TYPE INVOLVEMENT:")
    type_involvement = summary_stats['insertion_type_involvement']
    for ins_type, count in sorted(type_involvement.items(), key=lambda x: x[1], reverse=True):
        print(f"   {ins_type}: {count} false positives")
    
    print("\n3. PATTERN ANALYSIS:")
    patterns = summary_stats['pattern_statistics']
    total_fps = summary_stats['total_false_positives']
    for pattern, count in patterns.items():
        percentage = (count / total_fps) * 100 if total_fps > 0 else 0
        print(f"   {pattern}: {count} ({percentage:.1f}%)")
    
    print("\n4. ROOT CAUSE ANALYSIS:")
    root_causes = summary_stats['root_cause_summary']
    for cause, count in root_causes.items():
        percentage = (count / total_fps) * 100 if total_fps > 0 else 0
        print(f"   {cause}: {count} ({percentage:.1f}%)")

def main():
    """Main analysis function."""
    parser = argparse.ArgumentParser(description='Analyze false positive variants in nova simulation results')
    parser.add_argument('output_dir', help='Output directory containing nova simulation results')
    parser.add_argument('--output-prefix', default='nova', help='Output file prefix (default: nova)')
    
    args = parser.parse_args()
    
    # Check if output directory exists
    if not Path(args.output_dir).exists():
        print(f"Error: Output directory {args.output_dir} does not exist")
        return
    
    print("Starting false positive analysis...")
    
    # Load data
    variants_df, insertions_data = load_data(args.output_dir, args.output_prefix)
    read_to_insertion, read_to_sequence = create_insertion_lookup(insertions_data)
    
    # Categorize false positives
    print("\nCategorizing false positives...")
    fp_categories, successful_calls = categorize_false_positives(variants_df)
    
    # Analyze sequence patterns
    print("\nAnalyzing sequence patterns...")
    pattern_analysis = analyze_sequence_patterns(fp_categories, read_to_insertion, read_to_sequence)
    
    # Analyze root causes
    print("\nAnalyzing root causes...")
    root_cause_summary, detailed_causes = analyze_root_causes(fp_categories, pattern_analysis)
    
    # Generate summary statistics
    summary_stats = generate_summary_statistics(fp_categories, pattern_analysis, root_cause_summary)
    
    # Print summary
    print_analysis_summary(summary_stats)
    
    # Save results
    print("\nSaving results...")
    output_prefix = f"{args.output_dir}/false_positives"
    json_file, csv_file = save_results(fp_categories, pattern_analysis, detailed_causes, summary_stats, output_prefix)
    
    print(f"\nAnalysis complete!")
    print(f"Results saved to:")
    print(f"  - {json_file}")
    print(f"  - {csv_file}")

if __name__ == "__main__":
    main()