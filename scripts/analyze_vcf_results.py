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
from collections import defaultdict, Counter
import argparse
import sys

try:
    import pysam
except ImportError as e:
    print(f"Error: pysam is required to analyze BAM files. ({e})")
    sys.exit(1)

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

def track_nova_read_usage(nova_variants):
    """Track which nova reads appear in multiple variant calls.
    
    Returns:
        Dict mapping nova read names to list of variant indices where they appear
    """
    read_to_variants = defaultdict(list)
    
    for variant in nova_variants:
        for read_name in variant['nova_read_names']:
            read_to_variants[read_name].append(variant['index'])
    
    return read_to_variants

def categorize_variant(variant):
    """Categorize a single variant and return all its categories."""
    support_reads = variant['support_reads']
    nova_reads = variant['nova_reads']
    nova_fraction = nova_reads / support_reads
    
    # Simplified read composition - only 3 categories
    if support_reads == 1 and nova_reads == 1:
        composition = 'single_nova_only'
    elif nova_fraction >= 0.5:
        composition = 'majority_nova'
    else:
        composition = 'minority_nova'
    
    # Support level
    support = variant['SUPPORT']
    if support <= 2:
        support_level = 'Low (1-2)'
    elif support <= 5:
        support_level = 'Medium (3-5)'
    else:
        support_level = 'High (6+)'
    
    return {
        'svtype': variant['SVTYPE'],
        'precision': 'PRECISE' if variant['PRECISE'] else 'IMPRECISE',
        'support_level': support_level,
        'composition': composition,
        'nova_fraction': nova_fraction,
        'is_single_nova_only': composition == 'single_nova_only'
    }

def categorize_variants(nova_variants, mapping_verification=None):
    """Categorize variants by type, precision, and support.
    
    Args:
        nova_variants: List of variants to categorize
        mapping_verification: Optional dict mapping read names to verification status
    
    Returns:
        Dict with categorization results
    """
    categories = {
        'by_svtype': defaultdict(int),
        'by_precision': defaultdict(int),
        'by_support_level': defaultdict(int),
        'by_composition': defaultdict(int),
        'single_read_calls': {'single_nova_only': 0, 'total_variants': len(nova_variants)}
    }
    
    # If mapping verification is available, add detailed single call tracking
    if mapping_verification:
        categories['single_read_calls'].update({
            'single_nova_only_correct_mapping': 0,
            'single_nova_only_incorrect_mapping': 0,
            'single_nova_only_unknown_mapping': 0,
            'single_nova_only_total': 0
        })
    
    for variant in nova_variants:
        cat = categorize_variant(variant)
        
        # Increment counters
        categories['by_svtype'][cat['svtype']] += 1
        categories['by_precision'][cat['precision']] += 1
        categories['by_support_level'][cat['support_level']] += 1
        categories['by_composition'][cat['composition']] += 1
        
        if cat['is_single_nova_only']:
            categories['single_read_calls']['single_nova_only'] += 1
            
            # Add mapping verification tracking if available
            if mapping_verification and len(variant['nova_read_names']) == 1:
                nova_read = variant['nova_read_names'][0]
                mapping_status = mapping_verification.get(nova_read, 'unknown')
                
                categories['single_read_calls']['single_nova_only_total'] += 1
                
                if mapping_status == 'correct_mapping':
                    categories['single_read_calls']['single_nova_only_correct_mapping'] += 1
                elif mapping_status == 'incorrect_mapping':
                    categories['single_read_calls']['single_nova_only_incorrect_mapping'] += 1
                else:
                    categories['single_read_calls']['single_nova_only_unknown_mapping'] += 1
    
    return categories

def categorize_variants_detailed(nova_variants, mapping_verification, read_to_variants=None):
    """Categorize variants using detailed true/false positive definitions.
    
    True Positives: Single nova read variants with correct mapping
      - Exclusive: Nova read only appears in this variant
      - Shared: Nova read also appears in other variants
      
    False Positives:
      - Mapping Error: Single nova read with incorrect mapping
      - Multi-read: Multiple reads including nova reads
    
    Args:
        nova_variants: List of variants containing nova reads
        mapping_verification: Dict mapping read names to verification status
        read_to_variants: Optional pre-computed read usage tracking
        
    Returns:
        Dict with refined categorization
    """
    # Track nova read usage if not provided
    if read_to_variants is None:
        read_to_variants = track_nova_read_usage(nova_variants)
    
    categories = {
        'true_positives': {
            'exclusive': [],      # Single nova, correct mapping, not used elsewhere
            'shared': [],         # Single nova, correct mapping, used elsewhere
            'total': 0
        },
        'false_positives': {
            'mapping_error': [],  # Single nova, incorrect mapping
            'multi_read': [],     # Multiple reads including nova
            'total': 0
        },
        'summary': {
            'total_variants': len(nova_variants),
            'by_svtype': defaultdict(int),
            'by_precision': defaultdict(int)
        }
    }
    
    # Categorize each variant
    for variant in nova_variants:
        support_reads = variant['support_reads']
        nova_reads = variant['nova_reads']
        nova_read_names = variant['nova_read_names']
        
        # Track basic stats
        categories['summary']['by_svtype'][variant['SVTYPE']] += 1
        categories['summary']['by_precision']['PRECISE' if variant['PRECISE'] else 'IMPRECISE'] += 1
        
        # Single nova read variant
        if support_reads == 1 and nova_reads == 1:
            nova_read_name = nova_read_names[0]
            mapping_status = mapping_verification.get(nova_read_name, 'unknown')
            
            # Check if read is used in multiple variants
            read_usage_count = len(read_to_variants.get(nova_read_name, []))
            is_shared = read_usage_count > 1
            
            if mapping_status == 'correct_mapping':
                # True positive
                variant_info = variant.copy()
                variant_info['mapping_status'] = mapping_status
                variant_info['is_shared'] = is_shared
                variant_info['usage_count'] = read_usage_count
                variant_info['used_in_variants'] = read_to_variants.get(nova_read_name, [])
                
                if is_shared:
                    categories['true_positives']['shared'].append(variant_info)
                else:
                    categories['true_positives']['exclusive'].append(variant_info)
                    
                categories['true_positives']['total'] += 1
                
            else:
                # False positive - mapping error
                variant_info = variant.copy()
                variant_info['mapping_status'] = mapping_status
                variant_info['category'] = 'mapping_error'
                
                categories['false_positives']['mapping_error'].append(variant_info)
                categories['false_positives']['total'] += 1
                
        else:
            # Multi-read variant (false positive)
            variant_info = variant.copy()
            variant_info['category'] = 'multi_read'
            variant_info['nova_fraction'] = nova_reads / support_reads
            
            # Track nova read sharing for multi-read variants too
            shared_nova_reads = []
            for read_name in nova_read_names:
                if len(read_to_variants.get(read_name, [])) > 1:
                    shared_nova_reads.append(read_name)
            
            variant_info['shared_nova_reads'] = shared_nova_reads
            variant_info['has_shared_nova_reads'] = len(shared_nova_reads) > 0
            
            categories['false_positives']['multi_read'].append(variant_info)
            categories['false_positives']['total'] += 1
    
    return categories

def analyze_insertion_types(nova_variants, insertions_file):
    """Analyze which insertion types are being detected."""
    try:
        # Load insertion metadata
        with open(insertions_file, 'r') as f:
            insertions = json.load(f)
        
        print(f"Loaded {len(insertions)} insertion records")
        
        # Use unified lookup creation
        read_to_insertion = create_insertion_lookup(insertions)
        
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

def get_alignment_comparisons(nova_variants, mapping_verification, read_to_insertion):
    """Convert mapping verification results to alignment comparison format for compatibility.
    
    This function bridges the gap between the new unified mapping verification
    and the existing report generation that expects alignment comparison data.
    """
    alignment_comparisons = []
    
    for variant in nova_variants:
        for nova_read_name in variant['nova_read_names']:
            if nova_read_name in mapping_verification and nova_read_name in read_to_insertion:
                original_info = read_to_insertion[nova_read_name]
                mapping_status = mapping_verification[nova_read_name]
                
                # Convert mapping status to position match
                position_match = mapping_status == 'correct_mapping'
                
                alignment_comparisons.append({
                    'variant_index': variant['index'],
                    'read_name': nova_read_name,
                    'original_chrom': original_info['original_chr'],
                    'original_pos': original_info['original_pos'],
                    'modified_chrom': original_info['original_chr'] if position_match else None,
                    'modified_pos': original_info['original_pos'] if position_match else None,
                    'position_match': position_match,
                    'position_diff': 0 if position_match else None
                })
    
    return alignment_comparisons

def verify_mapping_locations(nova_variants, insertions_file, modified_bam_file, position_tolerance=1000):
    """
    Verify that nova reads in variants map to the same locations as their original base reads.
    
    Args:
        nova_variants: List of variants containing nova reads
        insertions_file: Path to JSON file with insertion records
        modified_bam_file: Path to BAM file with modified reads
        position_tolerance: Maximum allowed position difference (default: 1000bp)
        
    Returns:
        Dictionary mapping read names to mapping verification status
    """
    try:
        # Load insertion metadata to get original positions
        with open(insertions_file, 'r') as f:
            insertions = json.load(f)
        
        # Use unified lookup creation (without size info for efficiency)
        read_to_insertion = create_insertion_lookup(insertions, include_size=False)
        
        # Extract just the position info we need
        read_to_original_pos = {}
        for read_name, info in read_to_insertion.items():
            if info['original_chr'] and info['original_pos'] > 0:
                read_to_original_pos[read_name] = {
                    'chrom': info['original_chr'],
                    'pos': info['original_pos']
                }
        
        print(f"Loaded original positions for {len(read_to_original_pos)} nova reads")
        
        # Verify mapping locations for nova reads in variants
        mapping_verification = {}
        
        if not Path(modified_bam_file).exists():
            print(f"Warning: Modified BAM file {modified_bam_file} not found")
            # Return all reads as "unknown" status
            for variant in nova_variants:
                for nova_read_name in variant.get('nova_read_names', []):
                    mapping_verification[nova_read_name] = 'unknown'
            return mapping_verification
        
        # Create index of reads we need to check
        target_reads = set()
        for variant in nova_variants:
            target_reads.update(variant.get('nova_read_names', []))
        
        with pysam.AlignmentFile(modified_bam_file, 'rb') as modified_bam:
            # Iterate through BAM file once to find all target reads
            for read in modified_bam.fetch():
                if read.query_name in target_reads:
                    read_name = read.query_name
                    
                    if read_name in read_to_original_pos:
                        original_info = read_to_original_pos[read_name]
                        
                        # Compare positions
                        original_chrom = original_info['chrom']
                        original_pos = original_info['pos']
                        modified_chrom = read.reference_name
                        modified_pos = read.reference_start
                        
                        # Determine mapping status
                        if (original_chrom == modified_chrom and 
                            modified_pos is not None and
                            abs(original_pos - modified_pos) <= position_tolerance):
                            mapping_verification[read_name] = 'correct_mapping'
                        else:
                            mapping_verification[read_name] = 'incorrect_mapping'
                    else:
                        # Read not found in insertion records
                        mapping_verification[read_name] = 'no_original_position'
        
        # Mark any remaining target reads as not found in BAM
        for read_name in target_reads:
            if read_name not in mapping_verification:
                if read_name in read_to_original_pos:
                    mapping_verification[read_name] = 'not_found_in_bam'
                else:
                    mapping_verification[read_name] = 'no_original_position'
        
        print(f"Verified mapping locations for {len(mapping_verification)} nova reads")
        return mapping_verification
        
    except Exception as e:
        print(f"Error verifying mapping locations: {e}")
        # Return all reads with unknown status
        mapping_verification = {}
        for variant in nova_variants:
            for nova_read_name in variant.get('nova_read_names', []):
                mapping_verification[nova_read_name] = 'error'
        return mapping_verification

def compare_insertion_sizes(nova_variants, insertions_file):
    """Compare generated insertion sizes with detected variant SVLEN values."""
    try:
        # Load insertion metadata to get generated sizes
        with open(insertions_file, 'r') as f:
            insertions = json.load(f)
        
        # Use unified lookup creation with size info
        read_to_insertion = create_insertion_lookup(insertions, include_size=True)
        
        # Filter to only reads with valid sizes
        reads_with_sizes = {name: info for name, info in read_to_insertion.items() 
                           if info.get('insertion_length', 0) > 0}
        
        print(f"Loaded insertion sizes for {len(reads_with_sizes)} nova reads")
        
        # Compare with variant SVLEN values
        size_comparisons = []
        
        for variant in nova_variants:
            variant_svlen = abs(variant.get('SVLEN', 0))  # Use absolute value
            
            for nova_read_name in variant['nova_read_names']:
                if nova_read_name in reads_with_sizes:
                    insertion_info = reads_with_sizes[nova_read_name]
                    generated_size = insertion_info['insertion_length']
                    
                    # Calculate size accuracy
                    size_diff = abs(generated_size - variant_svlen)
                    size_ratio = min(generated_size, variant_svlen) / max(generated_size, variant_svlen) if max(generated_size, variant_svlen) > 0 else 0
                    size_accuracy = None
                    
                    # Categorize accuracy
                    if size_diff == 0:
                        size_accuracy = 'exact'
                    elif size_diff <= 10:
                        size_accuracy = 'close'
                    elif size_ratio >= 0.8:
                        size_accuracy = 'reasonable'
                    
                    size_comparisons.append({
                        'variant_index': variant['index'],
                        'read_name': nova_read_name,
                        'generated_size': generated_size,
                        'variant_svlen': variant_svlen,
                        'size_diff': size_diff,
                        'size_ratio': size_ratio,
                        'size_accuracy': size_accuracy,
                        'insertion_type': insertion_info['insertion_type']
                    })
        
        print(f"Compared insertion sizes for {len(size_comparisons)} nova reads")
        return size_comparisons
        
    except Exception as e:
        print(f"Error comparing insertion sizes: {e}")
        return []

def create_insertion_lookup(insertions_data, include_size=True):
    """Create lookup dictionary for insertion metadata.
    
    Args:
        insertions_data: List of insertion records
        include_size: Whether to include insertion size info (default: True)
    
    Returns:
        Dict mapping modified read names to insertion metadata
    """
    read_to_insertion = {}
    
    for insertion in insertions_data:
        read_name = insertion.get('modified_read_name', '')
        if read_name:
            metadata = {
                'insertion_type': insertion.get('insertion_type', 'unknown'),
                'insertion_id': insertion.get('insertion_id', ''),
                'original_chr': insertion.get('original_chr', ''),
                'original_pos': insertion.get('original_pos', 0),
                'base_read_name': insertion.get('base_read_name', ''),
                'original_read_name': insertion.get('original_read_name', '')
            }
            
            if include_size:
                # Handle both direct insertion_length and nested structure
                insertion_length = insertion.get('insertion_length', 0)
                if insertion_length == 0 and 'insertion_sequence' in insertion:
                    insertion_length = insertion.get('insertion_sequence', {}).get('length', 0)
                metadata['insertion_length'] = insertion_length
            
            read_to_insertion[read_name] = metadata
    
    return read_to_insertion

def analyze_false_positives(nova_variants, read_to_insertion):
    """Analyze false positive variants and their patterns."""
    
    # Separate successful calls from false positives
    successful_calls = []
    false_positives = []
    
    for variant in nova_variants:
        cat = categorize_variant(variant)
        if cat['is_single_nova_only']:
            successful_calls.append(variant)
        else:
            false_positives.append(variant)
    
    print(f"\nFalse Positive Analysis:")
    print(f"  Successful calls (single nova only): {len(successful_calls)}")
    print(f"  False positives: {len(false_positives)}")
    
    if not false_positives:
        # Return complete structure even with no false positives
        return {
            'false_positives': [],
            'patterns': {},
            'summary': {
                'total_false_positives': 0,
                'by_composition': {},
                'with_genomic_clustering': 0,
                'with_identical_sequences': 0,
                'clustering_rate': 0.0,
                'identical_sequence_rate': 0.0
            }
        }
    
    # Analyze patterns in false positives
    fp_patterns = {}
    category_counts = Counter()
    
    for variant in false_positives:
        variant_idx = variant['index']
        cat = categorize_variant(variant)
        category_counts[cat['composition']] += 1
        
        # Analyze sequence patterns for this variant
        nova_reads = variant['nova_read_names']
        insertion_details = []
        
        for read_name in nova_reads:
            if read_name in read_to_insertion:
                details = read_to_insertion[read_name].copy()
                details['read_name'] = read_name
                insertion_details.append(details)
        
        # Pattern analysis
        insertion_types = [d['insertion_type'] for d in insertion_details]
        insertion_ids = [d['insertion_id'] for d in insertion_details]
        original_positions = [(d['original_chr'], d['original_pos']) for d in insertion_details]
        
        # Check for patterns
        type_counts = Counter(insertion_types)
        id_counts = Counter(insertion_ids)
        has_duplicate_types = any(count > 1 for count in type_counts.values())
        has_identical_sequences = any(count > 1 for count in id_counts.values())
        
        # Check for genomic clustering (within 5kb)
        has_genomic_clustering = False
        position_clusters = defaultdict(list)
        for i, (chrom, pos) in enumerate(original_positions):
            position_clusters[chrom].append(pos)
        
        for chrom, positions in position_clusters.items():
            if len(positions) > 1:
                positions.sort()
                for i in range(len(positions) - 1):
                    if abs(positions[i+1] - positions[i]) <= 5000:
                        has_genomic_clustering = True
                        break
        
        fp_patterns[variant_idx] = {
            'composition': cat['composition'],
            'nova_reads': variant['nova_reads'],
            'support_reads': variant['support_reads'],
            'insertion_types': insertion_types,
            'has_duplicate_types': has_duplicate_types,
            'has_identical_sequences': has_identical_sequences,
            'has_genomic_clustering': has_genomic_clustering,
            'num_unique_types': len(type_counts),
            'num_unique_sequences': len(id_counts)
        }
    
    # Summary statistics
    total_fps = len(false_positives)
    clustering_fps = sum(1 for p in fp_patterns.values() if p['has_genomic_clustering'])
    identical_seq_fps = sum(1 for p in fp_patterns.values() if p['has_identical_sequences'])
    
    summary = {
        'total_false_positives': total_fps,
        'by_composition': dict(category_counts),
        'with_genomic_clustering': clustering_fps,
        'with_identical_sequences': identical_seq_fps,
        'clustering_rate': (clustering_fps / total_fps * 100) if total_fps > 0 else 0,
        'identical_sequence_rate': (identical_seq_fps / total_fps * 100) if total_fps > 0 else 0
    }
    
    print(f"  False positive composition: {dict(category_counts)}")
    print(f"  With genomic clustering: {clustering_fps} ({summary['clustering_rate']:.1f}%)")
    print(f"  With identical sequences: {identical_seq_fps} ({summary['identical_sequence_rate']:.1f}%)")
    
    return {
        'false_positives': false_positives,
        'patterns': fp_patterns,
        'summary': summary
    }

def calculate_expected_counts(insertions_data):
    """Calculate expected counts directly from insertions data."""
    type_counts = Counter()
    
    for insertion in insertions_data:
        insertion_type = insertion.get('insertion_type', 'unknown')
        type_counts[insertion_type] += 1
    
    return dict(type_counts)

def generate_detailed_summary_report(nova_variants, categories, type_detection, insertions_data, read_to_variants, size_comparisons=None):
    """Generate a comprehensive summary report using detailed categorization."""
    print("\n" + "="*60)
    print("NOVA VARIANT DETECTION ANALYSIS")
    print("="*60)
    
    # Calculate expected counts from insertions data
    expected_counts = calculate_expected_counts(insertions_data)
    total_expected = sum(expected_counts.values())

    print(f"\nExpected variants with nova reads: {total_expected}")
    print(f"\nTotal variants with nova reads: {categories['summary']['total_variants']}")
    
    print("\n1. TRUE POSITIVES (Single nova with correct mapping):")
    tp_exclusive = len(categories['true_positives']['exclusive'])
    tp_shared = len(categories['true_positives']['shared'])
    tp_total = categories['true_positives']['total']
    
    print(f"   Exclusive (nova read unique to variant): {tp_exclusive}")
    print(f"   Shared (nova read in multiple variants): {tp_shared}")
    print(f"   TOTAL TRUE POSITIVES: {tp_total}")
    
    # Analyze shared reads
    if tp_shared > 0:
        shared_variants = categories['true_positives']['shared']
        avg_usage = np.mean([v['usage_count'] for v in shared_variants])
        max_usage = max(v['usage_count'] for v in shared_variants)
        print(f"   - Average usage per shared read: {avg_usage:.1f}")
        print(f"   - Maximum usage count: {max_usage}")
    
    print("\n2. FALSE POSITIVES:")
    fp_mapping = len(categories['false_positives']['mapping_error'])
    fp_multi = len(categories['false_positives']['multi_read'])
    fp_total = categories['false_positives']['total']
    
    print(f"   Mapping errors (single nova, wrong position): {fp_mapping}")
    print(f"   Multi-read variants: {fp_multi}")
    print(f"   TOTAL FALSE POSITIVES: {fp_total}")
    
    # Analyze multi-read variants
    if fp_multi > 0:
        multi_variants = categories['false_positives']['multi_read']
        with_shared = sum(1 for v in multi_variants if v['has_shared_nova_reads'])
        avg_nova_fraction = np.mean([v['nova_fraction'] for v in multi_variants])
        print(f"   - Multi-read variants with shared nova reads: {with_shared}")
        print(f"   - Average nova fraction in multi-read variants: {avg_nova_fraction:.2f}")

    print("\n3. PERFORMANCE METRICS:")
    print(f"   RECALL: {tp_total}/{total_expected} ({tp_total/total_expected*100:.1f}%)")
    print(f"   PRECISION: {tp_total}/{tp_total+fp_total} ({tp_total/(tp_total+fp_total)*100:.1f}%)")
    print(f"   F1 SCORE (Harmonic mean of recall and precision): {2*tp_total/(2*tp_total+fp_total+(total_expected-tp_total)):.2f}")

    print("\n4. NOVA READ REUSE ANALYSIS:")
    # Count reads that appear in multiple variants
    reused_reads = sum(1 for read, variants in read_to_variants.items() if len(variants) > 1)
    total_unique_nova_reads = len(read_to_variants)
    
    print(f"   Total unique nova reads: {total_unique_nova_reads}")
    print(f"   Reads appearing in multiple variants: {reused_reads} ({reused_reads/total_unique_nova_reads*100:.1f}%)")
    
    # Distribution of read usage
    usage_counts = Counter(len(variants) for variants in read_to_variants.values())
    print("   Read usage distribution:")
    for count, num_reads in sorted(usage_counts.items()):
        print(f"     {count} variant(s): {num_reads} reads")
    
    print("\n5. DETECTION BY INSERTION TYPE:")
    for ins_type in sorted(expected_counts.keys()):
        expected = expected_counts[ins_type]
        type_variants = type_detection.get(ins_type, [])
        
        # Count true positives for this type
        tp_count = 0
        for variant in type_variants:
            if variant['support_reads'] == 1 and variant['nova_reads'] == 1:
                nova_read = variant['nova_read_names'][0]
                # Check if this variant is in our true positives
                for tp_var in (categories['true_positives']['exclusive'] + 
                             categories['true_positives']['shared']):
                    if tp_var['index'] == variant['index']:
                        tp_count += 1
                        break
        
        detection_rate = (tp_count / expected * 100) if expected > 0 else 0
        print(f"   {ins_type}: {tp_count}/{expected} ({detection_rate:.1f}%)")
    
    print("\n6. VARIANT CHARACTERISTICS:")
    print("   By SV Type:")
    for svtype, count in sorted(categories['summary']['by_svtype'].items()):
        print(f"     {svtype}: {count}")
    
    print("   By Precision Flag:")
    for precision, count in sorted(categories['summary']['by_precision'].items()):
        print(f"     {precision}: {count}")
    
    # Size comparison analysis
    if size_comparisons:
        print("\n7. INSERTION SIZE ACCURACY:")
        total_size_comparisons = len(size_comparisons)
        exact_matches = sum(1 for c in size_comparisons if c['size_accuracy'] == 'exact')
        close_matches = sum(1 for c in size_comparisons if c['size_accuracy'] == 'close')
        reasonable_matches = sum(1 for c in size_comparisons if c['size_accuracy'] == 'reasonable')
        
        print(f"   Total size comparisons: {total_size_comparisons}")
        print(f"   Exact matches (0bp diff): {exact_matches} ({exact_matches/total_size_comparisons*100:.1f}%)")
        print(f"   Close matches (≤10bp diff): {close_matches} ({close_matches/total_size_comparisons*100:.1f}%)")
        print(f"   Reasonable matches (≥80% ratio): {reasonable_matches} ({reasonable_matches/total_size_comparisons*100:.1f}%)")

def generate_summary_report(nova_variants, categories, type_detection, insertions_data, alignment_comparisons=None, size_comparisons=None, fp_analysis=None):
    """Generate a comprehensive summary report."""
    print("\n" + "="*60)
    print("NOVA VARIANT DETECTION ANALYSIS")
    print("="*60)
    
    print(f"\nTotal variants with nova reads: {len(nova_variants)}")
    
    # Calculate expected counts from insertions data
    expected_counts = calculate_expected_counts(insertions_data)
    
    print("\n1. TRUE POSITIVE ANALYSIS (Single nova-Only Calls):")
    total_expected = 0
    total_single_nova_calls = 0
    unique_nova_reads = set()
    
    # Collect all unique nova reads
    for variant in nova_variants:
        for read_name in variant.get('nova_read_names', []):
            unique_nova_reads.add(read_name)
    
    for ins_type in sorted(expected_counts.keys()):
        expected = expected_counts[ins_type]
        type_variants = type_detection.get(ins_type, [])
        single_nova_calls = sum(1 for v in type_variants if v.get('support_reads', 0) == 1 and v.get('nova_reads', 0) == 1)
        
        # Count unique nova reads for this type
        type_unique_reads = set()
        for variant in type_variants:
            for read_name in variant.get('nova_read_names', []):
                type_unique_reads.add(read_name)
        
        true_positive_rate = (single_nova_calls / expected * 100) if expected > 0 else 0
        read_utilization = (len(type_unique_reads) / expected * 100) if expected > 0 else 0
        
        print(f"   {ins_type}: {single_nova_calls}/{expected} true positives ({true_positive_rate:.1f}%), {len(type_unique_reads)} reads utilized ({read_utilization:.1f}%)")
        total_expected += expected
        total_single_nova_calls += single_nova_calls
    
    overall_true_positive_rate = (total_single_nova_calls / total_expected * 100) if total_expected > 0 else 0
    overall_read_utilization = (len(unique_nova_reads) / total_expected * 100) if total_expected > 0 else 0
    print(f"   OVERALL: {total_single_nova_calls}/{total_expected} true positives ({overall_true_positive_rate:.1f}%)")
    print(f"   READ UTILIZATION: {len(unique_nova_reads)}/{total_expected} nova reads detected ({overall_read_utilization:.1f}%)")
    
    print("\n2. DETECTION QUALITY ANALYSIS:")
    total_variants = len(nova_variants)
    
    # Check if mapping verification was performed
    single_calls = categories['single_read_calls']
    if 'single_nova_only_correct_mapping' in single_calls:
        # Mapping-aware analysis
        correct_mapping = single_calls['single_nova_only_correct_mapping']
        incorrect_mapping = single_calls['single_nova_only_incorrect_mapping']
        unknown_mapping = single_calls['single_nova_only_unknown_mapping']
        total_single_calls = single_calls['single_nova_only_total']
        
        true_positives = correct_mapping
        mapping_errors = incorrect_mapping
        unknown_mapping_status = unknown_mapping
        false_positives = total_variants - total_single_calls
        
        true_positive_rate = (true_positives / total_expected * 100) if total_expected > 0 else 0
        mapping_error_rate = (mapping_errors / total_single_calls * 100) if total_single_calls > 0 else 0
        
        print(f"   Total variants detected: {total_variants}")
        print(f"   Single nova-only calls: {total_single_calls}")
        print(f"     - Correct mapping (true positives): {correct_mapping} ({true_positive_rate:.1f}%)")
        print(f"     - Incorrect mapping (mapping errors): {incorrect_mapping} ({mapping_error_rate:.1f}%)")
        print(f"     - Unknown mapping status: {unknown_mapping}")
        print(f"   False positives (multi-read/mixed): {false_positives}")
    else:
        # Original analysis (mapping verification disabled or failed)
        false_positives = total_variants - total_single_nova_calls
        detection_quality = (total_single_nova_calls / total_variants * 100) if total_variants > 0 else 0
        false_positive_rate = (false_positives / total_variants * 100) if total_variants > 0 else 0
        
        print(f"   Total variants detected: {total_variants}")
        print(f"   True positives (single nova-only): {total_single_nova_calls} ({detection_quality:.1f}%)")
        print(f"   False positives (multi-read/mixed): {false_positives} ({false_positive_rate:.1f}%)")
        print("   Note: Mapping verification not performed - use refined analysis for better accuracy")
    
    print("\n3. READ COMPOSITION BREAKDOWN:")
    for comp_type, count in sorted(categories['by_composition'].items()):
        print(f"   {comp_type}: {count}")
    
    print("\n4. VARIANT TYPES (SVTYPE):")
    for svtype, count in sorted(categories['by_svtype'].items()):
        print(f"   {svtype}: {count}")
    
    print("\n5. PRECISION:")
    for precision, count in sorted(categories['by_precision'].items()):
        print(f"   {precision}: {count}")
    
    # Detailed statistics
    print("\n6. DETAILED STATISTICS:")
    supports = [v['SUPPORT'] for v in nova_variants]
    nova_fractions = [v['nova_reads']/v['support_reads'] for v in nova_variants]
    
    print(f"   Support reads - Mean: {np.mean(supports):.1f}, Median: {np.median(supports):.1f}")
    print(f"   nova fraction - Mean: {np.mean(nova_fractions):.2f}, Median: {np.median(nova_fractions):.2f}")
    
    # Alignment position analysis
    if alignment_comparisons:
        print("\n7. ALIGNMENT POSITION ANALYSIS:")
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
        print("\n8. INSERTION SIZE ANALYSIS:")
        total_size_comparisons = len(size_comparisons)
        exact_matches = sum(1 for c in size_comparisons if c['size_accuracy'] == 'exact')
        close_matches = sum(1 for c in size_comparisons if c['size_accuracy'] == 'close')
        reasonable_matches = sum(1 for c in size_comparisons if c['size_accuracy'] == 'reasonable')
        
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
    
    # False positive analysis summary
    if fp_analysis:
        print("\n9. FALSE POSITIVE ANALYSIS:")
        fp_summary = fp_analysis['summary']
        print(f"   Total false positives: {fp_summary['total_false_positives']}")
        print(f"   Composition breakdown: {fp_summary['by_composition']}")
        print(f"   With genomic clustering: {fp_summary['with_genomic_clustering']} ({fp_summary['clustering_rate']:.1f}%)")
        print(f"   With identical sequences: {fp_summary['with_identical_sequences']} ({fp_summary['identical_sequence_rate']:.1f}%)")

def save_detailed_tabular_data(categories, type_detection, size_comparisons, output_file, read_to_variants):
    """Save detailed variant data in tabular format for advanced visualizations."""
    
    variant_records = []
    
    # Process true positives - exclusive
    for variant in categories['true_positives']['exclusive']:
        base_record = {
            'variant_index': variant['index'],
            'chrom': variant['CHROM'],
            'pos': variant['POS'],
            'svtype': variant['SVTYPE'],
            'svlen': variant['SVLEN'],
            'precise': variant['PRECISE'],
            'support_reads': variant['support_reads'],
            'nova_reads': variant['nova_reads'],
            'variant_category': 'true_positive_exclusive',
            'mapping_status': variant['mapping_status'],
            'nova_read_reused': False,
            'nova_read_usage_count': 1,
            'reused_in_variants': []
        }
        
        # Add size comparison if available
        nova_read_name = variant['nova_read_names'][0]
        size_comp = next((c for c in size_comparisons if c['read_name'] == nova_read_name), None) if size_comparisons else None
        
        if size_comp:
            base_record.update({
                'nova_read_name': nova_read_name,
                'insertion_type': size_comp['insertion_type'],
                'generated_size': size_comp['generated_size'],
                'size_difference': size_comp['size_diff'],
                'size_ratio': size_comp['size_ratio'],
                'size_accuracy': size_comp['size_accuracy']
            })
        else:
            base_record.update({
                'nova_read_name': nova_read_name,
                'insertion_type': None,
                'generated_size': None,
                'size_difference': None,
                'size_ratio': None,
                'size_accuracy': None
            })
        
        variant_records.append(base_record)
    
    # Process true positives - shared
    for variant in categories['true_positives']['shared']:
        nova_read_name = variant['nova_read_names'][0]
        other_variants = [v for v in variant['used_in_variants'] if v != variant['index']]
        
        base_record = {
            'variant_index': variant['index'],
            'chrom': variant['CHROM'],
            'pos': variant['POS'],
            'svtype': variant['SVTYPE'],
            'svlen': variant['SVLEN'],
            'precise': variant['PRECISE'],
            'support_reads': variant['support_reads'],
            'nova_reads': variant['nova_reads'],
            'variant_category': 'true_positive_shared',
            'mapping_status': variant['mapping_status'],
            'nova_read_reused': True,
            'nova_read_usage_count': variant['usage_count'],
            'reused_in_variants': other_variants
        }
        
        # Add size comparison if available
        size_comp = next((c for c in size_comparisons if c['read_name'] == nova_read_name), None) if size_comparisons else None
        
        if size_comp:
            base_record.update({
                'nova_read_name': nova_read_name,
                'insertion_type': size_comp['insertion_type'],
                'generated_size': size_comp['generated_size'],
                'size_difference': size_comp['size_diff'],
                'size_ratio': size_comp['size_ratio'],
                'size_accuracy': size_comp['size_accuracy']
            })
        else:
            base_record.update({
                'nova_read_name': nova_read_name,
                'insertion_type': None,
                'generated_size': None,
                'size_difference': None,
                'size_ratio': None,
                'size_accuracy': None
            })
        
        variant_records.append(base_record)
    
    # Process false positives - mapping errors
    for variant in categories['false_positives']['mapping_error']:
        nova_read_name = variant['nova_read_names'][0]
        usage_count = len(read_to_variants.get(nova_read_name, [variant['index']]))
        
        base_record = {
            'variant_index': variant['index'],
            'chrom': variant['CHROM'],
            'pos': variant['POS'],
            'svtype': variant['SVTYPE'],
            'svlen': variant['SVLEN'],
            'precise': variant['PRECISE'],
            'support_reads': variant['support_reads'],
            'nova_reads': variant['nova_reads'],
            'variant_category': 'false_positive_mapping',
            'mapping_status': variant['mapping_status'],
            'nova_read_reused': usage_count > 1,
            'nova_read_usage_count': usage_count,
            'reused_in_variants': read_to_variants.get(nova_read_name, [])
        }
        
        # Add placeholder fields for consistency
        base_record.update({
            'nova_read_name': nova_read_name,
            'insertion_type': None,
            'generated_size': None,
            'size_difference': None,
            'size_ratio': None,
            'size_accuracy': None
        })
        
        variant_records.append(base_record)
    
    # Process false positives - multi-read
    for variant in categories['false_positives']['multi_read']:
        # For multi-read variants, create one record per nova read
        if variant['nova_reads'] > 0:
            for nova_read_name in variant['nova_read_names']:
                usage_count = len(read_to_variants.get(nova_read_name, [variant['index']]))
                
                base_record = {
                    'variant_index': variant['index'],
                    'chrom': variant['CHROM'],
                    'pos': variant['POS'],
                    'svtype': variant['SVTYPE'],
                    'svlen': variant['SVLEN'],
                    'precise': variant['PRECISE'],
                    'support_reads': variant['support_reads'],
                    'nova_reads': variant['nova_reads'],
                    'variant_category': 'false_positive_multi',
                    'mapping_status': 'multi_read',
                    'nova_read_reused': usage_count > 1,
                    'nova_read_usage_count': usage_count,
                    'reused_in_variants': read_to_variants.get(nova_read_name, []),
                    'nova_fraction': variant['nova_fraction']
                }
                
                # Add size comparison if available
                size_comp = next((c for c in size_comparisons if c['read_name'] == nova_read_name), None) if size_comparisons else None
                
                if size_comp:
                    base_record.update({
                        'nova_read_name': nova_read_name,
                        'insertion_type': size_comp['insertion_type'],
                        'generated_size': size_comp['generated_size'],
                        'size_difference': size_comp['size_diff'],
                        'size_ratio': size_comp['size_ratio'],
                        'size_accuracy': size_comp['size_accuracy']
                    })
                else:
                    base_record.update({
                        'nova_read_name': nova_read_name,
                        'insertion_type': None,
                        'generated_size': None,
                        'size_difference': None,
                        'size_ratio': None,
                        'size_accuracy': None
                    })
                
                variant_records.append(base_record)
    
    # Convert to DataFrame and save
    df = pd.DataFrame(variant_records)
    
    # Convert lists to strings for CSV compatibility
    df['reused_in_variants'] = df['reused_in_variants'].apply(lambda x: ','.join(map(str, x)) if x else '')
    
    # Save as CSV
    df.to_csv(output_file, index=False)
    
    print(f"CSV data saved to: {output_file}")
    print(f"Shape: {df.shape[0]} rows × {df.shape[1]} columns")
    
    return df

def save_tabular_data(nova_variants, categories, type_detection, alignment_comparisons, size_comparisons, fp_analysis, output_file, mapping_verification):
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
        
        # Get unified categorization
        cat = categorize_variant(variant)
        
        # Get false positive pattern data if available
        fp_pattern = fp_analysis['patterns'].get(variant_idx, {}) if fp_analysis else {}
        
        # Basic variant information with unified categorization
        base_record = {
            'variant_index': variant_idx,
            'chrom': variant['CHROM'],
            'pos': variant['POS'],
            'svtype': variant['SVTYPE'],
            'svlen': variant['SVLEN'],
            'precise': variant['PRECISE'],
            'support_reads': variant['support_reads'],
            'nova_reads': variant['nova_reads'],
            'nova_fraction': cat['nova_fraction'],
            'composition': cat['composition'],
            'support_level': cat['support_level'],
            'precision_category': cat['precision'],
            'is_single_nova_only': cat['is_single_nova_only'],
            
            # False positive pattern data
            'has_genomic_clustering': fp_pattern.get('has_genomic_clustering', False),
            'has_identical_sequences': fp_pattern.get('has_identical_sequences', False),
            'has_duplicate_types': fp_pattern.get('has_duplicate_types', False),
            'num_unique_types': fp_pattern.get('num_unique_types', 0),
            'num_unique_sequences': fp_pattern.get('num_unique_sequences', 0)
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
        
        # Add mapping verification data if available (for single nova reads)
        if (mapping_verification is not None and 
            variant['support_reads'] == 1 and variant['nova_reads'] == 1):
            nova_read_name = variant['nova_read_names'][0]
            mapping_status = mapping_verification.get(nova_read_name, 'unknown')
            
            base_record.update({
                'mapping_verification_status': mapping_status,
                'is_refined_true_positive': mapping_status == 'correct_mapping',
                'is_mapping_error': mapping_status == 'incorrect_mapping'
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
                    'size_accuracy': size_comp['size_accuracy']
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
    df.to_csv(output_file, index=False)
    
    print(f"CSV data saved to: {output_file}")
    print(f"Shape: {df.shape[0]} rows × {df.shape[1]} columns")
    
    return df

def generate_detailed_summary_json(categories, expected_counts, read_to_variants, output_file):
    """Generate a JSON summary using detailed categorization."""
    
    total_expected = sum(expected_counts.values())
    tp_total = categories['true_positives']['total']
    fp_total = categories['false_positives']['total']
    
    # Nova read reuse statistics
    reused_reads = sum(1 for read, variants in read_to_variants.items() if len(variants) > 1)
    total_unique_nova_reads = len(read_to_variants)
    
    # Usage distribution
    usage_counts = Counter(len(variants) for variants in read_to_variants.values())
    
    summary = {
        'analysis_timestamp': pd.Timestamp.now().isoformat(),
        'overall_metrics': {
            'total_variants_detected': categories['summary']['total_variants'],
            'total_expected': total_expected,
            'total_unique_nova_reads': total_unique_nova_reads,
            'reads_appearing_in_multiple_variants': reused_reads,
            'read_reuse_rate': (reused_reads / total_unique_nova_reads * 100) if total_unique_nova_reads > 0 else 0
        },
        'true_positives': {
            'total': tp_total,
            'exclusive': len(categories['true_positives']['exclusive']),
            'shared': len(categories['true_positives']['shared']),
            'detection_rate': (tp_total / total_expected * 100) if total_expected > 0 else 0
        },
        'false_positives': {
            'total': fp_total,
            'mapping_errors': len(categories['false_positives']['mapping_error']),
            'multi_read_variants': len(categories['false_positives']['multi_read'])
        },
        'variant_characteristics': {
            'by_svtype': dict(categories['summary']['by_svtype']),
            'by_precision': dict(categories['summary']['by_precision'])
        },
        'read_usage_distribution': dict(usage_counts),
        'detection_by_insertion_type': {}
    }
    
    # Add per-type detection rates if we have the data
    if expected_counts:
        for ins_type, expected in expected_counts.items():
            # Count true positives for this type
            tp_count = 0
            for variant in (categories['true_positives']['exclusive'] + 
                          categories['true_positives']['shared']):
                # We'd need insertion type info added to variants for this
                # For now, just use placeholder
                pass
            
            summary['detection_by_insertion_type'][ins_type] = {
                'expected': expected,
                'true_positives': tp_count,
                'detection_rate': (tp_count / expected * 100) if expected > 0 else 0
            }
    
    with open(output_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Summary JSON saved to: {output_file}")
    return summary

def generate_summary_json(df, expected_counts, fp_analysis, output_file):
    """Generate a lightweight JSON summary from CSV data."""
    total_variants = len(df)
    single_nova_only = len(df[df['is_single_nova_only'] == True])
    
    # Basic metrics
    composition_counts = df['composition'].value_counts().to_dict()
    svtype_counts = df['svtype'].value_counts().to_dict()
    precision_counts = df['precision_category'].value_counts().to_dict()
    
    # Calculate true positive rates by insertion type if available
    type_performance = {}
    if 'insertion_type' in df.columns:
        for ins_type in expected_counts.keys():
            type_df = df[df['insertion_type'] == ins_type]
            if len(type_df) > 0:
                tp_count = len(type_df[type_df['is_single_nova_only'] == True])
                expected = expected_counts[ins_type]
                type_performance[ins_type] = {
                    'expected': expected,
                    'total_detected': len(type_df),
                    'true_positives': tp_count,
                    'true_positive_rate': (tp_count / expected * 100) if expected > 0 else 0,
                }
    
    # Overall metrics
    total_expected = sum(expected_counts.values()) if expected_counts else 0
    overall_tp_rate = (single_nova_only / total_expected * 100) if total_expected > 0 else 0
    detection_quality = (single_nova_only / total_variants * 100) if total_variants > 0 else 0
    
    summary = {
        'analysis_timestamp': pd.Timestamp.now().isoformat(),
        'overall_metrics': {
            'total_variants_detected': total_variants,
            'true_positives': single_nova_only,
            'false_positives': total_variants - single_nova_only,
            'total_expected': total_expected,
            'true_positive_rate': overall_tp_rate,
            'detection_quality': detection_quality
        },
        'composition_breakdown': composition_counts,
        'variant_types': svtype_counts,
        'precision_breakdown': precision_counts,
        'by_insertion_type': type_performance,
        'false_positive_analysis': fp_analysis['summary'] if fp_analysis else {}
    }
    
    with open(output_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Summary JSON saved to: {output_file}")
    return summary

def main():
    """Main analysis function."""
    parser = argparse.ArgumentParser(description='Analyze nova VCF results for variant detection performance')
    parser.add_argument('output_dir', help='Output directory containing nova simulation results')
    parser.add_argument('--output-prefix', default='nova', help='Output file prefix (default: nova)')
    parser.add_argument('--original-bam', default='tests/test_data/test_reads.bam', 
                       help='Path to original BAM file (default: tests/test_data/test_reads.bam)')

    args = parser.parse_args()
    
    # Construct file paths using provided output directory
    output_dir = Path(args.output_dir)
    if not output_dir.exists():
        print(f"Error: Output directory {output_dir} does not exist")
        return
    
    simulation_jl = output_dir / f"{args.output_prefix}_simulation.jl"
    insertions_json = output_dir / f"{args.output_prefix}_insertions.json"
    original_bam = args.original_bam
    modified_bam = output_dir / f"{args.output_prefix}_modified_reads.bam"
    
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
    
    print("Verifying mapping locations...")
    mapping_verification = verify_mapping_locations(nova_variants, str(insertions_json), str(modified_bam))
    
    print("Analyzing insertion types...")
    type_detection = analyze_insertion_types(nova_variants, insertions_json)
    
    # Create insertion lookup
    print("Creating insertion lookup...")
    with open(insertions_json, 'r') as f:
        insertions_data = json.load(f)
    read_to_insertion = create_insertion_lookup(insertions_data)
    
    print("Comparing insertion sizes...")
    size_comparisons = compare_insertion_sizes(nova_variants, insertions_json)
    
    # Track nova read usage
    print("Tracking nova read usage across variants...")
    read_to_variants = track_nova_read_usage(nova_variants)
    
    # Perform categorization
    print("Categorizing variants...")
    categories = categorize_variants_detailed(nova_variants, mapping_verification, read_to_variants)
    
    # Generate reports
    generate_detailed_summary_report(nova_variants, categories, type_detection, 
                                  insertions_data, read_to_variants, size_comparisons)
    
    # Save tabular data
    print("\nGenerating tabular data...")
    csv_file = output_dir / f"{args.output_prefix}_analysis.csv"
    tabular_df = save_detailed_tabular_data(categories, type_detection, size_comparisons, 
                                         str(csv_file), read_to_variants)
    
    # Generate summary JSON
    print("Generating summary JSON...")
    expected_counts = calculate_expected_counts(insertions_data)
    summary_json = output_dir / f"{args.output_prefix}_analysis_summary.json"
    generate_detailed_summary_json(categories, expected_counts, read_to_variants, str(summary_json))

if __name__ == "__main__":
    main()