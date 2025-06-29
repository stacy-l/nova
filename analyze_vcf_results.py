#!/usr/bin/env python3
"""
Analyze Sniffles2 VCF results to identify variants supported by nova simulated reads.

This script loads the simulation VCF data and categorizes variants based on:
- SVTYPE (insertion, deletion, etc.)
- Precision (PRECISE/IMPRECISE)
- Read support (nova reads vs total reads)
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
import joblib
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess

# Try to import pysam, but handle gracefully if it fails
try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError as e:
    print(f"Warning: pysam not available ({e}). BAM analysis will be skipped.")
    PYSAM_AVAILABLE = False

def load_vcf_data(jl_file):
    """Load VCF data from joblib file."""
    try:
        df = joblib.load(jl_file)
        print(f"Loaded {len(df)} variants from {jl_file}")
        return df
    except Exception as e:
        print(f"Error loading {jl_file}: {e}")
        return None

def parse_rnames(rnames_data):
    """Parse RNAMES field to extract read names."""
    if pd.isna(rnames_data) or rnames_data == '':
        return []
    
    # Handle different data types
    if isinstance(rnames_data, tuple):
        return list(rnames_data)
    elif isinstance(rnames_data, list):
        return rnames_data
    elif isinstance(rnames_data, str):
        # Handle comma-separated strings
        return [name.strip() for name in rnames_data.split(',')]
    else:
        # Convert to string and try to parse
        return [str(rnames_data)]

def identify_nova_variants(df):
    """Identify variants that contain nova simulated reads."""
    nova_variants = []
    
    print(f"Scanning {len(df)} variants for nova reads...")
    
    for i, (idx, row) in enumerate(df.iterrows()):
        if i % 10000 == 0:
            print(f"  Processed {i:,} variants...")
            
        # Handle both uppercase and lowercase column names
        rnames_col = 'RNAMES' if 'RNAMES' in row else 'rnames'
        if rnames_col not in row or pd.isna(row[rnames_col]):
            continue
            
        rnames = parse_rnames(row[rnames_col])
        nova_reads = [r for r in rnames if 'nova' in str(r)]
        
        if nova_reads:
            # Use flexible column name lookup
            variant_info = {
                'index': idx,
                'CHROM': row.get('CHROM', row.get('chrom', '')),
                'POS': row.get('POS', row.get('start', 0)),
                'SVTYPE': row.get('SVTYPE', row.get('svtype', 'UNKNOWN')),
                'SVLEN': row.get('SVLEN', row.get('svlen', 0)),
                'PRECISE': row.get('PRECISE', row.get('precise', False)),
                'SUPPORT': row.get('SUPPORT', row.get('support', 0)),
                'support_reads': len(rnames), # backup for case where SUPPORT != len(rnames)
                'nova_reads': len(nova_reads),
                'nova_read_names': nova_reads,
                'all_read_names': rnames
            }
            nova_variants.append(variant_info)
    
    print(f"Found {len(nova_variants)} variants with nova reads!")
    return nova_variants

def categorize_variants(nova_variants):
    """Categorize variants by type, precision, and support."""
    categories = {
        'by_svtype': defaultdict(int),
        'by_precision': defaultdict(int),
        'by_support_level': defaultdict(int),
        'by_nova_fraction': defaultdict(int),
        'by_read_composition': defaultdict(int),
        'single_read_calls': {'single_nova_only': 0, 'total_variants': len(nova_variants)}
    }
    
    # Detailed read support breakdown
    support_breakdown = defaultdict(lambda: {'total': 0, 'all_nova': 0, 'mixed': 0})
    
    for variant in nova_variants:
        # By SVTYPE
        svtype = variant['SVTYPE']
        categories['by_svtype'][svtype] += 1
        
        # By precision
        precision = 'PRECISE' if variant['PRECISE'] else 'IMPRECISE'
        categories['by_precision'][precision] += 1
        
        # By support level
        support = variant['SUPPORT']
        if support <= 2:
            support_level = 'Low (1-2)'
        elif support <= 5:
            support_level = 'Medium (3-5)'
        else:
            support_level = 'High (6+)'
        categories['by_support_level'][support_level] += 1
        
        # By nova read fraction
        nova_fraction = variant['nova_reads'] / variant['support_reads']
        if nova_fraction == 1.0:
            fraction_level = 'All nova'
        elif nova_fraction >= 0.5:
            fraction_level = 'Majority nova'
        else:
            fraction_level = 'Minority nova'
        categories['by_nova_fraction'][fraction_level] += 1
        
        # Detailed read composition analysis
        support_reads = variant['support_reads']
        nova_reads = variant['nova_reads']
        non_nova_reads = support_reads - nova_reads
        
        # Track successful calls (single nova read, no other reads)
        if support_reads == 1 and nova_reads == 1:
            categories['single_read_calls']['single_nova_only'] += 1
        
        # Read composition categories
        if support_reads == 1:
            if nova_reads == 1:
                categories['by_read_composition']['1_nova_only'] += 1
            else:
                categories['by_read_composition']['1_non_nova_only'] += 1
        elif support_reads == 2:
            if nova_reads == 2:
                categories['by_read_composition']['2_nova_only'] += 1
            elif nova_reads == 1:
                categories['by_read_composition']['2_mixed_1nova_1other'] += 1
            else:
                categories['by_read_composition']['2_non_nova_only'] += 1
        elif support_reads == 3:
            if nova_reads == 3:
                categories['by_read_composition']['3_nova_only'] += 1
            elif nova_reads == 2:
                categories['by_read_composition']['3_mixed_2nova_1other'] += 1
            elif nova_reads == 1:
                categories['by_read_composition']['3_mixed_1nova_2other'] += 1
            else:
                categories['by_read_composition']['3_non_nova_only'] += 1
        else:
            # 4+ reads
            if nova_reads == support_reads:
                categories['by_read_composition'][f'{support_reads}_all_nova'] += 1
            elif nova_reads >= support_reads / 2:
                categories['by_read_composition'][f'{support_reads}_majority_nova'] += 1
            else:
                categories['by_read_composition'][f'{support_reads}_minority_nova'] += 1
        
        # Support breakdown for detailed analysis
        support_key = str(support)
        support_breakdown[support_key]['total'] += 1
        if nova_reads == support_reads:
            support_breakdown[support_key]['all_nova'] += 1
        else:
            support_breakdown[support_key]['mixed'] += 1
    
    categories['support_breakdown'] = dict(support_breakdown)
    return categories

def analyze_insertion_types(nova_variants, insertions_file):
    """Analyze which insertion types are being detected."""
    try:
        # Load insertion metadata
        with open(insertions_file, 'r') as f:
            insertions = json.load(f)
        
        print(f"Loaded {len(insertions)} insertion records")
        
        # Create mapping from read name to insertion metadata
        read_to_insertion = {}
        for insertion in insertions:
            read_name = insertion.get('modified_read_name', '')
            if read_name:
                read_to_insertion[read_name] = {
                    'insertion_type': insertion.get('insertion_type', 'unknown'),
                    'insertion_id': insertion.get('insertion_id', ''),
                    'original_read_name': insertion.get('original_read_name', ''),
                    'insertion_length': insertion.get('insertion_sequence', {}).get('length', 0)
                }
        
        # Categorize variants by insertion type
        type_detection = defaultdict(list)
        matched_reads = 0
        
        for variant in nova_variants:
            for read_name in variant['nova_read_names']:
                if read_name in read_to_insertion:
                    matched_reads += 1
                    insertion_info = read_to_insertion[read_name]
                    insertion_type = insertion_info['insertion_type']
                    # Add insertion metadata to variant info
                    variant_with_insertion = variant.copy()
                    variant_with_insertion['insertion_info'] = insertion_info
                    type_detection[insertion_type].append(variant_with_insertion)
        
        print(f"Matched {matched_reads} nova reads to insertion metadata")
        return type_detection
        
    except Exception as e:
        print(f"Error analyzing insertion types: {e}")
        return defaultdict(list)

def compare_alignment_positions(nova_variants, insertions_file, original_bam_file, modified_bam_file):
    """Compare alignment positions between original and modified nova reads."""
    try:
        # Load insertion metadata to get original positions
        with open(insertions_file, 'r') as f:
            insertions = json.load(f)
        
        # Create mapping from modified read name to original position
        read_to_original_pos = {}
        for insertion in insertions:
            modified_name = insertion.get('modified_read_name', '')
            original_name = insertion.get('base_read_name', '')
            original_chrom = insertion.get('original_chr', '')
            original_pos = insertion.get('original_pos', 0)
            
            if modified_name:
                read_to_original_pos[modified_name] = {
                    'original_name': original_name,
                    'chrom': original_chrom,
                    'pos': original_pos
                }
        
        print(f"Loaded original positions for {len(read_to_original_pos)} nova reads")
        
        # Read modified BAM to get new alignment positions
        alignment_comparisons = []
        
        if not Path(modified_bam_file).exists():
            print(f"Warning: Modified BAM file {modified_bam_file} not found")
            return []
        
        with pysam.AlignmentFile(modified_bam_file, 'rb') as modified_bam:
            for variant in nova_variants:
                for nova_read_name in variant['nova_read_names']:
                    if nova_read_name in read_to_original_pos:
                        original_info = read_to_original_pos[nova_read_name]
                        
                        # Find the read in the modified BAM
                        found_read = False
                        for read in modified_bam.fetch():
                            if read.query_name == nova_read_name:
                                found_read = True
                                
                                # Compare positions (binary match/no-match)
                                original_chrom = original_info['chrom']
                                original_pos = original_info['pos']
                                modified_chrom = read.reference_name
                                modified_pos = read.reference_start
                                
                                position_match = (original_chrom == modified_chrom and 
                                                abs(original_pos - modified_pos) <= 1000)  # Allow 1kb tolerance
                                
                                alignment_comparisons.append({
                                    'variant_index': variant['index'],
                                    'read_name': nova_read_name,
                                    'original_chrom': original_chrom,
                                    'original_pos': original_pos,
                                    'modified_chrom': modified_chrom,
                                    'modified_pos': modified_pos,
                                    'position_match': position_match,
                                    'position_diff': abs(original_pos - modified_pos) if original_chrom == modified_chrom else None
                                })
                                break
                        
                        if not found_read:
                            alignment_comparisons.append({
                                'variant_index': variant['index'],
                                'read_name': nova_read_name,
                                'original_chrom': original_info['chrom'],
                                'original_pos': original_info['pos'],
                                'modified_chrom': None,
                                'modified_pos': None,
                                'position_match': False,
                                'position_diff': None
                            })
        
        print(f"Compared alignment positions for {len(alignment_comparisons)} nova reads")
        return alignment_comparisons
        
    except Exception as e:
        print(f"Error comparing alignment positions: {e}")
        return []

def compare_insertion_sizes(nova_variants, insertions_file):
    """Compare generated insertion sizes with detected variant SVLEN values."""
    try:
        # Load insertion metadata to get generated sizes
        with open(insertions_file, 'r') as f:
            insertions = json.load(f)
        
        # Create mapping from modified read name to insertion size and type
        read_to_insertion_size = {}
        read_to_insertion_type = {}
        for insertion in insertions:
            modified_name = insertion.get('modified_read_name', '')
            insertion_size = insertion.get('insertion_length', 0)
            insertion_type = insertion.get('insertion_type', 'unknown')
            
            if modified_name:
                if insertion_size > 0:
                    read_to_insertion_size[modified_name] = insertion_size
                read_to_insertion_type[modified_name] = insertion_type
        
        print(f"Loaded insertion sizes for {len(read_to_insertion_size)} nova reads")
        
        # Compare with variant SVLEN values
        size_comparisons = []
        
        for variant in nova_variants:
            variant_svlen = abs(variant.get('SVLEN', 0))  # Use absolute value
            
            for nova_read_name in variant['nova_read_names']:
                if nova_read_name in read_to_insertion_size:
                    generated_size = read_to_insertion_size[nova_read_name]
                    
                    # Calculate size accuracy
                    size_diff = abs(generated_size - variant_svlen)
                    size_ratio = min(generated_size, variant_svlen) / max(generated_size, variant_svlen) if max(generated_size, variant_svlen) > 0 else 0
                    
                    # Categorize accuracy
                    exact_match = (size_diff == 0)
                    close_match = (size_diff <= 10)  # Within 10bp
                    reasonable_match = (size_ratio >= 0.8)  # Within 20% of each other
                    
                    size_comparisons.append({
                        'variant_index': variant['index'],
                        'read_name': nova_read_name,
                        'generated_size': generated_size,
                        'variant_svlen': variant_svlen,
                        'size_diff': size_diff,
                        'size_ratio': size_ratio,
                        'exact_match': exact_match,
                        'close_match': close_match,
                        'reasonable_match': reasonable_match,
                        'insertion_type': read_to_insertion_type.get(nova_read_name, 'unknown')
                    })
        
        print(f"Compared insertion sizes for {len(size_comparisons)} nova reads")
        return size_comparisons
        
    except Exception as e:
        print(f"Error comparing insertion sizes: {e}")
        return []

def load_expected_counts(statistics_file):
    """Load expected counts from statistics file."""
    try:
        with open(statistics_file, 'r') as f:
            stats = json.load(f)
        return stats['insertion_statistics']['type_counts']
    except Exception as e:
        print(f"Warning: Could not load statistics file {statistics_file}: {e}")
        return {}

def generate_summary_report(nova_variants, categories, type_detection, alignment_comparisons=None, size_comparisons=None):
    """Generate a comprehensive summary report."""
    print("\n" + "="*60)
    print("NOVA VARIANT DETECTION ANALYSIS")
    print("="*60)
    
    print(f"\nTotal variants with nova reads: {len(nova_variants)}")
    
    # Load actual expected counts from statistics
    statistics_file = "output/nova_statistics.json"
    expected_counts = load_expected_counts(statistics_file)
    
    print("\n1. EXPECTED vs DETECTED INSERTION COUNTS:")
    total_expected = 0
    total_detected = 0
    for ins_type in sorted(expected_counts.keys()):
        expected = expected_counts[ins_type]
        detected = len(type_detection.get(ins_type, []))
        detection_rate = (detected / expected * 100) if expected > 0 else 0
        print(f"   {ins_type}: {detected}/{expected} ({detection_rate:.1f}%)")
        total_expected += expected
        total_detected += detected
    
    overall_rate = (total_detected / total_expected * 100) if total_expected > 0 else 0
    print(f"   OVERALL: {total_detected}/{total_expected} ({overall_rate:.1f}%)")
    
    print("\n2. SUCCESSFUL CALLS (Single Nova Read Only):")
    successful = categories['single_read_calls']['single_nova_only']
    success_rate = (successful / len(nova_variants) * 100) if nova_variants else 0
    print(f"   Single nova read only: {successful}/{len(nova_variants)} ({success_rate:.1f}%)")
    
    print("\n3. READ COMPOSITION BREAKDOWN:")
    for comp_type, count in sorted(categories['by_read_composition'].items()):
        print(f"   {comp_type}: {count}")
    
    print("\n4. SUPPORT BREAKDOWN (by number of reads):")
    for support_level, breakdown in sorted(categories['support_breakdown'].items(), key=lambda x: int(x[0])):
        total = breakdown['total']
        all_nova = breakdown['all_nova']
        mixed = breakdown['mixed']
        print(f"   {support_level} reads: {total} variants ({all_nova} all-nova, {mixed} mixed)")
    
    print("\n5. VARIANT TYPES (SVTYPE):")
    for svtype, count in sorted(categories['by_svtype'].items()):
        print(f"   {svtype}: {count}")
    
    print("\n6. PRECISION:")
    for precision, count in sorted(categories['by_precision'].items()):
        print(f"   {precision}: {count}")
    
    # Detailed statistics
    print("\n7. DETAILED STATISTICS:")
    supports = [v['SUPPORT'] for v in nova_variants]
    nova_fractions = [v['nova_reads']/v['support_reads'] for v in nova_variants]
    
    print(f"   Support reads - Mean: {np.mean(supports):.1f}, Median: {np.median(supports):.1f}")
    print(f"   Nova fraction - Mean: {np.mean(nova_fractions):.2f}, Median: {np.median(nova_fractions):.2f}")
    
    # Alignment position analysis
    if alignment_comparisons:
        print("\n8. ALIGNMENT POSITION ANALYSIS:")
        total_comparisons = len(alignment_comparisons)
        position_matches = sum(1 for c in alignment_comparisons if c['position_match'])
        match_rate = (position_matches / total_comparisons * 100) if total_comparisons > 0 else 0
        
        print(f"   Total alignments compared: {total_comparisons}")
        print(f"   Position matches (±1kb): {position_matches}/{total_comparisons} ({match_rate:.1f}%)")
        
        # Position differences for matches on same chromosome
        same_chrom_diffs = [c['position_diff'] for c in alignment_comparisons 
                           if c['position_diff'] is not None]
        if same_chrom_diffs:
            print(f"   Position difference - Mean: {np.mean(same_chrom_diffs):.1f}bp, Median: {np.median(same_chrom_diffs):.1f}bp")
    
    # Size comparison analysis
    if size_comparisons:
        print("\n9. INSERTION SIZE ANALYSIS:")
        total_size_comparisons = len(size_comparisons)
        exact_matches = sum(1 for c in size_comparisons if c['exact_match'])
        close_matches = sum(1 for c in size_comparisons if c['close_match'])
        reasonable_matches = sum(1 for c in size_comparisons if c['reasonable_match'])
        
        exact_rate = (exact_matches / total_size_comparisons * 100) if total_size_comparisons > 0 else 0
        close_rate = (close_matches / total_size_comparisons * 100) if total_size_comparisons > 0 else 0
        reasonable_rate = (reasonable_matches / total_size_comparisons * 100) if total_size_comparisons > 0 else 0
        
        print(f"   Total size comparisons: {total_size_comparisons}")
        print(f"   Exact matches (0bp diff): {exact_matches} ({exact_rate:.1f}%)")
        print(f"   Close matches (≤10bp diff): {close_matches} ({close_rate:.1f}%)")
        print(f"   Reasonable matches (≥80% ratio): {reasonable_matches} ({reasonable_rate:.1f}%)")
        
        # Size accuracy statistics
        size_diffs = [c['size_diff'] for c in size_comparisons]
        size_ratios = [c['size_ratio'] for c in size_comparisons]
        
        if size_diffs:
            print(f"   Size difference - Mean: {np.mean(size_diffs):.1f}bp, Median: {np.median(size_diffs):.1f}bp")
        if size_ratios:
            print(f"   Size ratio - Mean: {np.mean(size_ratios):.2f}, Median: {np.median(size_ratios):.2f}")

def save_detailed_results(nova_variants, categories, type_detection, alignment_comparisons, size_comparisons, output_file):
    """Save comprehensive analysis results to structured JSON."""
    
    # Calculate summary statistics
    expected_counts = load_expected_counts("output/nova_statistics.json")
    
    # Detection rates by type
    detection_rates = {}
    total_expected = 0
    total_detected = 0
    for ins_type in sorted(expected_counts.keys()):
        expected = expected_counts[ins_type]
        detected = len(type_detection.get(ins_type, []))
        detection_rates[ins_type] = {
            'expected': expected,
            'detected': detected,
            'rate': (detected / expected * 100) if expected > 0 else 0
        }
        total_expected += expected
        total_detected += detected
    
    # Alignment statistics
    alignment_stats = {}
    if alignment_comparisons:
        total_alignments = len(alignment_comparisons)
        position_matches = sum(1 for c in alignment_comparisons if c['position_match'])
        same_chrom_diffs = [c['position_diff'] for c in alignment_comparisons if c['position_diff'] is not None]
        
        alignment_stats = {
            'total_comparisons': total_alignments,
            'position_matches': position_matches,
            'match_rate': (position_matches / total_alignments * 100) if total_alignments > 0 else 0,
            'position_differences': {
                'mean': float(np.mean(same_chrom_diffs)) if same_chrom_diffs else 0,
                'median': float(np.median(same_chrom_diffs)) if same_chrom_diffs else 0,
                'values': same_chrom_diffs[:100]  # Save first 100 for plotting
            }
        }
    
    # Size comparison statistics
    size_stats = {}
    if size_comparisons:
        total_size_comparisons = len(size_comparisons)
        exact_matches = sum(1 for c in size_comparisons if c['exact_match'])
        close_matches = sum(1 for c in size_comparisons if c['close_match'])
        reasonable_matches = sum(1 for c in size_comparisons if c['reasonable_match'])
        
        size_diffs = [c['size_diff'] for c in size_comparisons]
        size_ratios = [c['size_ratio'] for c in size_comparisons]
        
        size_stats = {
            'total_comparisons': total_size_comparisons,
            'exact_matches': exact_matches,
            'close_matches': close_matches,
            'reasonable_matches': reasonable_matches,
            'exact_rate': (exact_matches / total_size_comparisons * 100) if total_size_comparisons > 0 else 0,
            'close_rate': (close_matches / total_size_comparisons * 100) if total_size_comparisons > 0 else 0,
            'reasonable_rate': (reasonable_matches / total_size_comparisons * 100) if total_size_comparisons > 0 else 0,
            'size_differences': {
                'mean': float(np.mean(size_diffs)) if size_diffs else 0,
                'median': float(np.median(size_diffs)) if size_diffs else 0,
                'values': size_diffs[:100]  # Save first 100 for plotting
            },
            'size_ratios': {
                'mean': float(np.mean(size_ratios)) if size_ratios else 0,
                'median': float(np.median(size_ratios)) if size_ratios else 0,
                'values': size_ratios[:100]  # Save first 100 for plotting
            }
        }
        
        # Size accuracy by insertion type
        size_by_type = {}
        for comp in size_comparisons:
            ins_type = comp['insertion_type']
            if ins_type not in size_by_type:
                size_by_type[ins_type] = {'exact': 0, 'close': 0, 'reasonable': 0, 'total': 0}
            
            size_by_type[ins_type]['total'] += 1
            if comp['exact_match']:
                size_by_type[ins_type]['exact'] += 1
            if comp['close_match']:
                size_by_type[ins_type]['close'] += 1
            if comp['reasonable_match']:
                size_by_type[ins_type]['reasonable'] += 1
        
        size_stats['by_insertion_type'] = size_by_type
    
    # Support statistics
    supports = [v['SUPPORT'] for v in nova_variants]
    nova_fractions = [v['nova_reads']/v['support_reads'] for v in nova_variants]
    
    output_data = {
        'metadata': {
            'analysis_timestamp': pd.Timestamp.now().isoformat(),
            'total_nova_variants': len(nova_variants),
            'total_expected_insertions': total_expected,
            'total_detected_insertions': total_detected,
            'overall_detection_rate': (total_detected / total_expected * 100) if total_expected > 0 else 0
        },
        'detection_rates': detection_rates,
        'variant_categories': {
            'by_svtype': dict(categories['by_svtype']),
            'by_precision': dict(categories['by_precision']),
            'by_support_level': dict(categories['by_support_level']),
            'by_nova_fraction': dict(categories['by_nova_fraction']),
            'by_read_composition': dict(categories['by_read_composition']),
            'support_breakdown': categories['support_breakdown']
        },
        'single_read_calls': categories['single_read_calls'],
        'alignment_analysis': alignment_stats,
        'size_analysis': size_stats,
        'summary_statistics': {
            'support_reads': {
                'mean': float(np.mean(supports)) if supports else 0,
                'median': float(np.median(supports)) if supports else 0,
                'values': supports[:100]  # First 100 for plotting
            },
            'nova_fractions': {
                'mean': float(np.mean(nova_fractions)) if nova_fractions else 0,
                'median': float(np.median(nova_fractions)) if nova_fractions else 0,
                'values': nova_fractions[:100]  # First 100 for plotting
            }
        }
    }
    
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\nStructured analysis results saved to: {output_file}")

def save_tabular_data(nova_variants, categories, type_detection, alignment_comparisons, size_comparisons, output_file):
    """Save variant data in tabular format for advanced visualizations."""
    
    # Create variant-level records with all categorical information
    variant_records = []
    
    # Create lookup dictionaries for quick access
    alignment_lookup = {c['variant_index']: c for c in alignment_comparisons} if alignment_comparisons else {}
    size_lookup = {}
    if size_comparisons:
        for comp in size_comparisons:
            variant_idx = comp['variant_index']
            if variant_idx not in size_lookup:
                size_lookup[variant_idx] = []
            size_lookup[variant_idx].append(comp)
    
    for variant in nova_variants:
        variant_idx = variant['index']
        
        # Basic variant information
        base_record = {
            'variant_index': variant_idx,
            'chrom': variant['CHROM'],
            'pos': variant['POS'],
            'svtype': variant['SVTYPE'],
            'svlen': variant['SVLEN'],
            'precise': variant['PRECISE'],
            'support_reads': variant['SUPPORT'],
            'support_reads': variant['support_reads'],
            'nova_reads': variant['nova_reads'],
            'nova_fraction': variant['nova_reads'] / variant['support_reads'],
            
            # Categorizations
            'is_single_read_call': variant['support_reads'] == 1 and variant['nova_reads'] == 1,
            'support_level': ('Low (1-2)' if variant['SUPPORT'] <= 2 else 
                            'Medium (3-5)' if variant['SUPPORT'] <= 5 else 'High (6+)'),
            'nova_fraction_category': ('All nova' if variant['nova_reads'] / variant['support_reads'] == 1.0 else
                                     'Majority nova' if variant['nova_reads'] / variant['support_reads'] >= 0.5 else
                                     'Minority nova'),
            'precision_category': 'PRECISE' if variant['PRECISE'] else 'IMPRECISE'
        }
        
        # Add alignment information if available
        if variant_idx in alignment_lookup:
            alignment = alignment_lookup[variant_idx]
            base_record.update({
                'alignment_position_match': alignment['position_match'],
                'position_difference': alignment['position_diff'],
                'original_chrom': alignment['original_chrom'],
                'original_pos': alignment['original_pos'],
                'modified_chrom': alignment['modified_chrom'],
                'modified_pos': alignment['modified_pos']
            })
        else:
            base_record.update({
                'alignment_position_match': None,
                'position_difference': None,
                'original_chrom': None,
                'original_pos': None,
                'modified_chrom': None,
                'modified_pos': None
            })
        
        # Handle multiple nova reads per variant (create one record per nova read)
        if variant_idx in size_lookup:
            for size_comp in size_lookup[variant_idx]:
                record = base_record.copy()
                record.update({
                    'nova_read_name': size_comp['read_name'],
                    'insertion_type': size_comp['insertion_type'],
                    'generated_size': size_comp['generated_size'],
                    'size_difference': size_comp['size_diff'],
                    'size_ratio': size_comp['size_ratio'],
                    'exact_size_match': size_comp['exact_match'],
                    'close_size_match': size_comp['close_match'],
                    'reasonable_size_match': size_comp['reasonable_match']
                })
                variant_records.append(record)
        else:
            # No size data available, but try to get insertion type from variant metadata
            insertion_type = None
            if 'insertion_info' in variant and variant['insertion_info']:
                insertion_type = variant['insertion_info'].get('insertion_type', None)
            
            record = base_record.copy()
            record.update({
                'nova_read_name': variant['nova_read_names'][0] if variant['nova_read_names'] else None,
                'insertion_type': insertion_type,
                'generated_size': None,
                'size_difference': None,
                'size_ratio': None,
                'exact_size_match': None,
                'close_size_match': None,
                'reasonable_size_match': None
            })
            variant_records.append(record)
    
    # Convert to DataFrame and save as CSV
    df = pd.DataFrame(variant_records)
    
    # Save as CSV
    csv_file = output_file.replace('.json', '_tabular.csv')
    df.to_csv(csv_file, index=False)
    
    print(f"Tabular data saved to: {csv_file}")
    print(f"Shape: {df.shape[0]} rows × {df.shape[1]} columns")
    
    # Print sample of available columns for reference
    print(f"Columns available for visualization: {', '.join(df.columns)}")
    
    return df

def main():
    """Main analysis function."""
    # File paths
    simulation_jl = "output/nova_simulation.jl"
    insertions_json = "output/nova_insertions.json"
    output_json = "output/nova_variant_analysis.json"
    original_bam = "tests/test_data/test_reads.bam"  # From statistics file
    modified_bam = "output/nova_modified_reads.bam"
    
    # Check if files exist
    if not Path(simulation_jl).exists():
        print(f"Error: {simulation_jl} not found")
        return
    
    if not Path(insertions_json).exists():
        print(f"Error: {insertions_json} not found")
        return
    
    # Load and analyze data
    print("Loading VCF data...")
    df = load_vcf_data(simulation_jl)
    if df is None:
        return
    
    print("Identifying variants with nova reads...")
    nova_variants = identify_nova_variants(df)
    
    if not nova_variants:
        print("No variants found with nova reads!")
        return
    
    print("Categorizing variants...")
    categories = categorize_variants(nova_variants)
    
    print("Analyzing insertion types...")
    type_detection = analyze_insertion_types(nova_variants, insertions_json)
    
    print("Comparing alignment positions...")
    alignment_comparisons = compare_alignment_positions(nova_variants, insertions_json, original_bam, modified_bam)
    
    print("Comparing insertion sizes...")
    size_comparisons = compare_insertion_sizes(nova_variants, insertions_json)
    
    # Generate reports
    generate_summary_report(nova_variants, categories, type_detection, alignment_comparisons, size_comparisons)
    save_detailed_results(nova_variants, categories, type_detection, alignment_comparisons, size_comparisons, output_json)
    
    # Save tabular data for advanced visualizations
    print("Generating tabular data for visualization...")
    tabular_df = save_tabular_data(nova_variants, categories, type_detection, alignment_comparisons, size_comparisons, output_json)

if __name__ == "__main__":
    main()