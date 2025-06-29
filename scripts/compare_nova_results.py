#!/usr/bin/env python3
"""
Compare nova simulation results between two output directories.

This script analyzes differences between two nova result directories by matching
variants based on coordinates, SVTYPE, and size, then categorizing them as
identical, dir1_only, or dir2_only for flexible downstream analysis.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
import argparse
import sys
from collections import defaultdict, Counter


def load_variants(directory, output_prefix="nova"):
    """Load variant data from tabular CSV file."""
    directory = Path(directory)
    
    # Load tabular variant data
    tabular_file = directory / f"{output_prefix}_variant_analysis_tabular.csv"
    if not tabular_file.exists():
        raise FileNotFoundError(f"Tabular data file not found: {tabular_file}")
    
    variants_df = pd.read_csv(tabular_file)
    print(f"Loaded {len(variants_df)} variants from {directory}")
    
    return variants_df


def match_variants(df1, df2, pos_tolerance=50, size_tolerance=0.2):
    """
    Match variants between two dataframes using coordinate, type, and size criteria.
    
    Args:
        df1, df2: DataFrames with variant data
        pos_tolerance: Maximum position difference (bp)
        size_tolerance: Maximum size difference as fraction (0.2 = 20%)
    
    Returns:
        matched_pairs: List of (idx1, idx2, match_type) tuples
        unmatched_df1: DataFrame of unmatched variants from df1
        unmatched_df2: DataFrame of unmatched variants from df2
    """
    
    # Group variants by unique variant_index to avoid duplicates
    df1_unique = df1.groupby('variant_index').first().reset_index()
    df2_unique = df2.groupby('variant_index').first().reset_index()
    
    matched_pairs = []
    matched_indices_1 = set()
    matched_indices_2 = set()
    
    # For each variant in df1, find best match in df2
    for idx1, var1 in df1_unique.iterrows():
        best_match = None
        best_score = float('inf')
        match_type = None
        
        # Look for matches in df2
        for idx2, var2 in df2_unique.iterrows():
            if idx2 in matched_indices_2:
                continue
                
            # Must match chromosome and SVTYPE exactly
            if var1['chrom'] != var2['chrom'] or var1['svtype'] != var2['svtype']:
                continue
            
            # Check position tolerance
            pos_diff = abs(var1['pos'] - var2['pos'])
            if pos_diff > pos_tolerance:
                continue
            
            # Check size tolerance (handle NaN values)
            size1 = var1.get('svlen', np.nan)
            size2 = var2.get('svlen', np.nan)
            
            # If both have sizes, check size tolerance
            if pd.notna(size1) and pd.notna(size2):
                if abs(size1) > 0:  # Avoid division by zero
                    size_diff_frac = abs(size1 - size2) / abs(size1)
                    if size_diff_frac > size_tolerance:
                        continue
                elif abs(size2) > 0:
                    size_diff_frac = abs(size1 - size2) / abs(size2)
                    if size_diff_frac > size_tolerance:
                        continue
                else:
                    size_diff_frac = 0  # Both are zero
            elif pd.isna(size1) and pd.isna(size2):
                size_diff_frac = 0  # Both missing, consider match
            else:
                continue  # One has size, other doesn't - no match
            
            # Calculate match score (lower is better)
            score = pos_diff + size_diff_frac * 100
            
            if score < best_score:
                best_match = idx2
                best_score = score
                # Determine match type
                if pos_diff == 0 and size_diff_frac < 0.01:
                    match_type = "exact"
                else:
                    match_type = "approximate"
        
        # Record best match if found
        if best_match is not None:
            matched_pairs.append((idx1, best_match, match_type))
            matched_indices_1.add(idx1)
            matched_indices_2.add(best_match)
    
    # Get unmatched variants
    unmatched_df1 = df1_unique[~df1_unique.index.isin(matched_indices_1)].copy()
    unmatched_df2 = df2_unique[~df2_unique.index.isin(matched_indices_2)].copy()
    
    print(f"Found {len(matched_pairs)} matched variant pairs")
    print(f"  - Exact matches: {sum(1 for _, _, t in matched_pairs if t == 'exact')}")
    print(f"  - Approximate matches: {sum(1 for _, _, t in matched_pairs if t == 'approximate')}")
    print(f"Unmatched variants: {len(unmatched_df1)} from dir1, {len(unmatched_df2)} from dir2")
    
    return matched_pairs, unmatched_df1, unmatched_df2


def categorize_variants(df1, df2, matched_pairs, unmatched_df1, unmatched_df2):
    """Categorize variants into identical, dir1_only, dir2_only."""
    
    # Group by unique variant_index
    df1_unique = df1.groupby('variant_index').first().reset_index()
    df2_unique = df2.groupby('variant_index').first().reset_index()
    
    results = []
    
    # Process matched pairs
    for idx1, idx2, match_type in matched_pairs:
        var1 = df1_unique.iloc[idx1]
        var2 = df2_unique.iloc[idx2]
        
        # Create comparison record
        record = {
            'category': 'identical',
            'match_type': match_type,
            'variant_index_dir1': var1['variant_index'],
            'variant_index_dir2': var2['variant_index'],
            'chrom': var1['chrom'],
            'pos_dir1': var1['pos'],
            'pos_dir2': var2['pos'],
            'pos_diff': abs(var1['pos'] - var2['pos']),
            'svtype': var1['svtype'],
            'svlen_dir1': var1.get('svlen', np.nan),
            'svlen_dir2': var2.get('svlen', np.nan),
            'support_reads_dir1': var1['support_reads'],
            'support_reads_dir2': var2['support_reads'],
            'nova_reads_dir1': var1['nova_reads'],
            'nova_reads_dir2': var2['nova_reads'],
            'is_single_read_call_dir1': var1['is_single_read_call'],
            'is_single_read_call_dir2': var2['is_single_read_call'],
            'nova_read_name_dir1': var1.get('nova_read_name', ''),
            'nova_read_name_dir2': var2.get('nova_read_name', ''),
            'insertion_type_dir1': var1.get('insertion_type', ''),
            'insertion_type_dir2': var2.get('insertion_type', ''),
            'precision_dir1': var1.get('precision_category', ''),
            'precision_dir2': var2.get('precision_category', '')
        }
        
        # Calculate size difference if both have sizes
        if pd.notna(record['svlen_dir1']) and pd.notna(record['svlen_dir2']):
            record['size_diff'] = abs(record['svlen_dir1'] - record['svlen_dir2'])
            if abs(record['svlen_dir1']) > 0:
                record['size_diff_frac'] = record['size_diff'] / abs(record['svlen_dir1'])
            else:
                record['size_diff_frac'] = 0
        else:
            record['size_diff'] = np.nan
            record['size_diff_frac'] = np.nan
        
        results.append(record)
    
    # Process unmatched variants from dir1
    for _, var in unmatched_df1.iterrows():
        record = {
            'category': 'dir1_only',
            'match_type': 'none',
            'variant_index_dir1': var['variant_index'],
            'variant_index_dir2': '',
            'chrom': var['chrom'],
            'pos_dir1': var['pos'],
            'pos_dir2': np.nan,
            'pos_diff': np.nan,
            'svtype': var['svtype'],
            'svlen_dir1': var.get('svlen', np.nan),
            'svlen_dir2': np.nan,
            'support_reads_dir1': var['support_reads'],
            'support_reads_dir2': np.nan,
            'nova_reads_dir1': var['nova_reads'],
            'nova_reads_dir2': np.nan,
            'is_single_read_call_dir1': var['is_single_read_call'],
            'is_single_read_call_dir2': np.nan,
            'nova_read_name_dir1': var.get('nova_read_name', ''),
            'nova_read_name_dir2': '',
            'insertion_type_dir1': var.get('insertion_type', ''),
            'insertion_type_dir2': '',
            'precision_dir1': var.get('precision_category', ''),
            'precision_dir2': '',
            'size_diff': np.nan,
            'size_diff_frac': np.nan
        }
        results.append(record)
    
    # Process unmatched variants from dir2
    for _, var in unmatched_df2.iterrows():
        record = {
            'category': 'dir2_only',
            'match_type': 'none',
            'variant_index_dir1': '',
            'variant_index_dir2': var['variant_index'],
            'chrom': var['chrom'],
            'pos_dir1': np.nan,
            'pos_dir2': var['pos'],
            'pos_diff': np.nan,
            'svtype': var['svtype'],
            'svlen_dir1': np.nan,
            'svlen_dir2': var.get('svlen', np.nan),
            'support_reads_dir1': np.nan,
            'support_reads_dir2': var['support_reads'],
            'nova_reads_dir1': np.nan,
            'nova_reads_dir2': var['nova_reads'],
            'is_single_read_call_dir1': np.nan,
            'is_single_read_call_dir2': var['is_single_read_call'],
            'nova_read_name_dir1': '',
            'nova_read_name_dir2': var.get('nova_read_name', ''),
            'insertion_type_dir1': '',
            'insertion_type_dir2': var.get('insertion_type', ''),
            'precision_dir1': '',
            'precision_dir2': var.get('precision_category', ''),
            'size_diff': np.nan,
            'size_diff_frac': np.nan
        }
        results.append(record)
    
    return pd.DataFrame(results)


def generate_summary_statistics(comparison_df):
    """Generate summary statistics for the comparison."""
    
    category_counts = comparison_df['category'].value_counts()
    
    # Match type breakdown for identical variants
    identical_variants = comparison_df[comparison_df['category'] == 'identical']
    match_type_counts = identical_variants['match_type'].value_counts() if len(identical_variants) > 0 else pd.Series()
    
    # SVTYPE breakdown by category
    svtype_breakdown = {}
    for category in ['identical', 'dir1_only', 'dir2_only']:
        cat_data = comparison_df[comparison_df['category'] == category]
        if len(cat_data) > 0:
            svtype_counts = cat_data['svtype'].value_counts()
            svtype_breakdown[category] = svtype_counts.to_dict()
        else:
            svtype_breakdown[category] = {}
    
    # Single read call analysis for identical variants
    single_read_analysis = {}
    if len(identical_variants) > 0:
        for col in ['is_single_read_call_dir1', 'is_single_read_call_dir2']:
            if col in identical_variants.columns:
                single_read_analysis[col] = identical_variants[col].value_counts().to_dict()
    
    summary = {
        'total_variants': len(comparison_df),
        'category_counts': category_counts.to_dict(),
        'match_type_breakdown': match_type_counts.to_dict(),
        'svtype_breakdown': svtype_breakdown,
        'single_read_call_analysis': single_read_analysis
    }
    
    return summary


def save_results(comparison_df, summary_stats, output_dir, dir1_name, dir2_name):
    """Save comparison results to files."""
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Save detailed comparison CSV
    csv_file = output_path / f"variant_comparison_{dir1_name}_vs_{dir2_name}.csv"
    comparison_df.to_csv(csv_file, index=False)
    print(f"Detailed comparison saved to: {csv_file}")
    print(f"CSV shape: {comparison_df.shape[0]} rows × {comparison_df.shape[1]} columns")
    
    # Save summary statistics JSON
    json_file = output_path / f"comparison_summary_{dir1_name}_vs_{dir2_name}.json"
    summary_output = {
        'metadata': {
            'analysis_type': 'nova_results_comparison',
            'timestamp': pd.Timestamp.now().isoformat(),
            'dir1_name': dir1_name,
            'dir2_name': dir2_name
        },
        'summary_statistics': summary_stats
    }
    
    with open(json_file, 'w') as f:
        json.dump(summary_output, f, indent=2, default=str)
    print(f"Summary statistics saved to: {json_file}")
    
    return csv_file, json_file


def print_summary(summary_stats, dir1_name, dir2_name):
    """Print a comprehensive comparison summary."""
    
    print("\n" + "="*60)
    print(f"NOVA RESULTS COMPARISON: {dir1_name} vs {dir2_name}")
    print("="*60)
    
    total_variants = summary_stats['total_variants']
    print(f"\nTotal variant records: {total_variants}")
    
    print("\n1. VARIANT CATEGORIES:")
    category_counts = summary_stats['category_counts']
    for category, count in category_counts.items():
        percentage = (count / total_variants) * 100 if total_variants > 0 else 0
        print(f"   {category}: {count} ({percentage:.1f}%)")
    
    print("\n2. MATCH TYPE BREAKDOWN (for identical variants):")
    match_types = summary_stats['match_type_breakdown']
    for match_type, count in match_types.items():
        print(f"   {match_type}: {count}")
    
    print("\n3. SVTYPE BREAKDOWN BY CATEGORY:")
    svtype_breakdown = summary_stats['svtype_breakdown']
    for category, svtypes in svtype_breakdown.items():
        if svtypes:
            print(f"   {category}:")
            for svtype, count in sorted(svtypes.items(), key=lambda x: x[1], reverse=True):
                print(f"     {svtype}: {count}")
    
    print("\n4. SINGLE READ CALL ANALYSIS (for identical variants):")
    single_read_analysis = summary_stats['single_read_call_analysis']
    for col, counts in single_read_analysis.items():
        dir_name = "dir1" if "dir1" in col else "dir2"
        print(f"   {dir_name}:")
        for value, count in counts.items():
            print(f"     {value}: {count}")


def main():
    """Main comparison function."""
    parser = argparse.ArgumentParser(description='Compare nova simulation results between two directories')
    parser.add_argument('dir1', help='First directory path')
    parser.add_argument('dir2', help='Second directory path')
    parser.add_argument('--output-dir', default='comparison_results', help='Output directory (default: comparison_results)')
    parser.add_argument('--output-prefix', default='nova', help='Output file prefix (default: nova)')
    parser.add_argument('--pos-tolerance', type=int, default=50, help='Position tolerance in bp (default: 50)')
    parser.add_argument('--size-tolerance', type=float, default=0.2, help='Size tolerance as fraction (default: 0.2)')
    
    args = parser.parse_args()
    
    # Check if directories exist
    if not Path(args.dir1).exists():
        print(f"Error: Directory {args.dir1} does not exist")
        return
    if not Path(args.dir2).exists():
        print(f"Error: Directory {args.dir2} does not exist")
        return
    
    # Extract directory names for output files
    dir1_name = Path(args.dir1).name
    dir2_name = Path(args.dir2).name
    
    print(f"Comparing nova results: {args.dir1} vs {args.dir2}")
    
    # Load variant data
    print("\nLoading variant data...")
    df1 = load_variants(args.dir1, args.output_prefix)
    df2 = load_variants(args.dir2, args.output_prefix)
    
    # Match variants
    print("\nMatching variants...")
    matched_pairs, unmatched_df1, unmatched_df2 = match_variants(
        df1, df2, 
        pos_tolerance=args.pos_tolerance, 
        size_tolerance=args.size_tolerance
    )
    
    # Categorize variants
    print("\nCategorizing variants...")
    comparison_df = categorize_variants(df1, df2, matched_pairs, unmatched_df1, unmatched_df2)
    
    # Generate summary statistics
    summary_stats = generate_summary_statistics(comparison_df)
    
    # Print summary
    print_summary(summary_stats, dir1_name, dir2_name)
    
    # Save results
    print("\nSaving results...")
    csv_file, json_file = save_results(comparison_df, summary_stats, args.output_dir, dir1_name, dir2_name)
    
    print(f"\nComparison complete!")
    print(f"Results saved to:")
    print(f"  - {csv_file}")
    print(f"  - {json_file}")


if __name__ == "__main__":
    main()