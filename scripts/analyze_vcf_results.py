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
from typing import List, Dict, Tuple, Set

try:
    import pysam
except ImportError as e:
    print(f"Error: pysam is required to analyze BAM files. ({e})")
    sys.exit(1)

def load_vcf_data(jl_file) -> pd.DataFrame:
    """Load VCF data from joblib file."""
    try:
        df = joblib.load(jl_file)
        print(f"Loaded {len(df)} variants from {jl_file}")
        return df
    except Exception as e:
        print(f"ERROR: Error loading {jl_file}: {e}")
        return None

def parse_rnames(rnames_data) -> List[str]:
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

def identify_nova_variants(df) -> List[Dict]:
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
            # Normalize FILTER field to handle both strings and lists
            filter_raw = row.get('FILTER', row.get('filter', 'UNKNOWN'))
            if isinstance(filter_raw, list):
                filter_status = ';'.join(str(f) for f in filter_raw) if filter_raw else 'UNKNOWN'
            elif filter_raw is None or filter_raw == '':
                filter_status = 'UNKNOWN'
            else:
                filter_status = str(filter_raw)
            
            # Use flexible column name lookup
            variant_info = {
                'original_index': idx,  # Track original dataframe index
                'ID': row.get('ID', row.get('id', '')),
                'CHROM': row.get('CHROM', row.get('chrom', '')),
                'POS': row.get('POS', row.get('start', 0)),
                'SVTYPE': row.get('SVTYPE', row.get('svtype', 'UNKNOWN')),
                'SVLEN': row.get('SVLEN', row.get('svlen', 0)),
                'PRECISE': row.get('PRECISE', row.get('precise', False)),
                'SUPPORT': row.get('SUPPORT', row.get('support', 0)),
                'FILTER': filter_status,
                'support_reads': len(rnames), # backup for case where SUPPORT != len(rnames)
                'nova_reads': len(nova_reads),
                'nova_read_names': nova_reads,
                'all_read_names': rnames
            }
            nova_variants.append(variant_info)
    
    print(f"Found {len(nova_variants)} variants with nova reads!")
    return nova_variants

def categorize_variants(nova_variants, mapping_verification, read_to_variants) -> Dict[str, Dict]:
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
    
    categories = {
        'true_positives': {
            'exclusive': [],      # Single nova, correct mapping, not used elsewhere
            'shared': [],         # Single nova, correct mapping, used elsewhere
            'total': 0,
            'by_filter': defaultdict(int)  # Track by FILTER status
        },
        'false_positives': {
            'mapping_error': [],  # Single nova, incorrect mapping
            'multi_read_majority_nova': [],     # Multiple reads (majority nova)
            'multi_read_minority_nova': [], # Multiple reads (minority nova)
            'total': 0,
            'by_filter': defaultdict(int)  # Track by FILTER status
        },
        'variant_filters': []  # Track FILTER for each variant
    }
    
    for variant in nova_variants:
        support_reads = variant['support_reads'] # count
        nova_reads = variant['nova_reads'] # count
        nova_read_names = variant['nova_read_names'] # list strs
        filter_status = variant.get('FILTER', 'UNKNOWN')
        
        # Track filter for this variant
        categories['variant_filters'].append(filter_status)
        
        if support_reads == 1 and nova_reads == 1:
            nova_read_name = nova_read_names[0]
            mapping_status = mapping_verification.get(nova_read_name, 'unknown')
            
            # Check if read is used in multiple variants
            read_usage_count = len(read_to_variants.get(nova_read_name, []))
            is_shared = read_usage_count > 1
            
            if mapping_status == 'correct_mapping':
                if is_shared:
                    categories['true_positives']['shared'].append(nova_read_name)
                else:
                    categories['true_positives']['exclusive'].append(nova_read_name)
                categories['true_positives']['total'] += 1
                categories['true_positives']['by_filter'][filter_status] += 1
                
            else:
                categories['false_positives']['mapping_error'].append(nova_read_name)
                categories['false_positives']['total'] += 1
                categories['false_positives']['by_filter'][filter_status] += 1
                
        else:
            # Multi-read variant (false positive)
            nova_fraction = nova_reads / support_reads
            
            if nova_fraction > 0.5:
                categories['false_positives']['multi_read_majority_nova'].extend(nova_read_names)
            else:
                categories['false_positives']['multi_read_minority_nova'].extend(nova_read_names)

            categories['false_positives']['total'] += 1
            categories['false_positives']['by_filter'][filter_status] += 1
    
    return categories

def find_mapping_locations(all_inserted_nova_read_names, modified_bam_file):
    """Find mapping locations for inserted nova reads in modified BAM file."""
    read_to_mapped_pos = defaultdict(dict)
    try:
        with pysam.AlignmentFile(modified_bam_file, 'rb') as modified_bam:
            # Iterate through BAM file once to find all target reads
            for read in modified_bam.fetch():
                if read.query_name in all_inserted_nova_read_names:
                    read_name = read.query_name
                    read_to_mapped_pos[read_name] = {
                        'mapped_chr': read.reference_name,
                        'mapped_pos': read.reference_start
                    }
        
        print(f"Identified mapping locations for {len(read_to_mapped_pos)} nova reads")
        return read_to_mapped_pos
        
    except Exception as e:
        print(f"ERROR: Error identifying mapping locations: {e}")
        for nova_read_name in all_inserted_nova_read_names:
            read_to_mapped_pos[nova_read_name] = {
                'mapped_chr': 'unknown',
                'mapped_pos': 'unknown'
            }
        return read_to_mapped_pos


def verify_mapping_locations(all_inserted_nova_read_names, read_to_insertion, read_to_mapped_pos, position_tolerance=1000) -> Tuple[Dict[str, str], Dict[str, int]]:
    """
    Verify that nova reads in variants map to the same locations as their original base reads.
    
    Args:
        all_inserted_nova_read_names: List of nova read names
        read_to_insertion: Dict mapping read names to insertion metadata
        read_to_mapped_pos: Dict mapping read names to mapped position metadata
        position_tolerance: Maximum allowed position difference (default: 1000bp)
        
    Returns:
        Dictionary mapping read names to mapping verification status
    """
    try:
        mapping_verification = {}
            
        for read_name in all_inserted_nova_read_names:
            mapped_chr = read_to_mapped_pos.get(read_name, {}).get('mapped_chr', 'unknown')  
            mapped_pos = read_to_mapped_pos.get(read_name, {}).get('mapped_pos', 'unknown')

            # should never be the case, otherwise we wouldn't have been able to pull the read from the bam
            original_chr = read_to_insertion.get(read_name, {}).get('original_chr', 'unknown')
            original_pos = read_to_insertion.get(read_name, {}).get('original_pos', 'unknown')

            if mapped_chr == 'unknown' or mapped_pos == 'unknown':
                mapping_verification[read_name] = 'unknown_mapped_position'
            elif original_chr == 'unknown' or original_pos == 'unknown':
                mapping_verification[read_name] = 'unknown_original_position'
            elif mapped_chr == original_chr and abs(mapped_pos - original_pos) <= position_tolerance:
                mapping_verification[read_name] = 'correct_mapping'
            else:
                mapping_verification[read_name] = 'incorrect_mapping'
        
        n_unknown_mapped = sum(status == 'unknown_mapped_position' for status in mapping_verification.values())
        n_unknown_original = sum(status == 'unknown_original_position' for status in mapping_verification.values())
        n_correct = sum(status == 'correct_mapping' for status in mapping_verification.values())
        n_incorrect = sum(status == 'incorrect_mapping' for status in mapping_verification.values())

        mapping_report = {'unknown_mapping': n_unknown_mapped, 
                          'unknown_original_mapping': n_unknown_original, 
                          'correct_mapping': n_correct, 
                          'incorrect_mapping': n_incorrect}
        
        print(f"Verified mapping locations for {n_correct + n_incorrect} nova reads")
        return mapping_verification, mapping_report
        
    except Exception as e:
        print(f"ERROR: Error verifying mapping locations: {e}")
        # Return all reads with unknown status
        mapping_verification = {}
        for read_name in all_inserted_nova_read_names:
            mapping_verification[read_name] = 'unknown'
        mapping_report = {'unknown_mapping': len(all_inserted_nova_read_names),
                          'unknown_original_mapping': len(all_inserted_nova_read_names)}
        return mapping_verification, mapping_report

def create_insertion_lookup(insertions_data) -> Tuple[Dict[str, Dict], Set[str]]:
    """Create lookup dictionary for insertion metadata.
    
    Args:
        insertions_data: List of insertion records
        include_size: Whether to include insertion size info (default: True)
    
    Returns:
        Dict mapping modified read names to insertion metadata
    """
    read_to_insertion = {}
    all_inserted_nova_read_names = set()
    
    for insertion in insertions_data:
        read_name = insertion.get('modified_read_name', '')
        if read_name:
            metadata = {
                'base_read_name': insertion.get('base_read_name', ''),
                'modified_read_name': insertion.get('modified_read_name', ''),
                'original_chr': insertion.get('original_chr', ''),
                'original_pos': insertion.get('original_pos', 0),
                'insertion_id': insertion.get('insertion_id', ''),
                'insertion_type': insertion.get('insertion_type', 'unknown'),
                'insertion_length': insertion.get('insertion_length', 0),
                'insertion_pos': insertion.get('insertion_pos', 0)
            }
            
            read_to_insertion[read_name] = metadata
            all_inserted_nova_read_names.add(read_name)
    
    return read_to_insertion, all_inserted_nova_read_names

def calculate_observed_vars_by_ins_type(read_to_mapped_pos, read_to_insertion) -> Dict[str, List[Dict]]:
    """Obtain variant indices of observed insertion types (ex. random, simple repeat, predefined)."""
    try:
        observed_vars_by_ins_type = defaultdict(list)
        observed_reads = set(read_to_mapped_pos.keys())

        for read_name in observed_reads:
            if read_name in read_to_insertion:
                insertion_info = read_to_insertion[read_name]
                insertion_type = insertion_info['insertion_type']
                observed_vars_by_ins_type[insertion_type].append(read_name)
        return observed_vars_by_ins_type
        
    except Exception as e:
        print(f"ERROR: analyzing insertion types: {e}")
        return defaultdict(list)
    
def _calculate_observed_var_metrics(nova_variants) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """Obtain variant indices for SVTYPEs and precision flags."""
    try:
        observed_vars_by_svtype = defaultdict(list)
        observed_vars_by_precision = defaultdict(list)
        for variant in nova_variants:
            svtype = variant.get('SVTYPE', '')
            precision = variant.get('PRECISE', '')
            observed_vars_by_svtype[svtype].append(variant.get('ID', ''))
            observed_vars_by_precision[precision].append(variant.get('ID', ''))
        return observed_vars_by_svtype, observed_vars_by_precision
        
    except Exception as e:
        print(f"ERROR: Error analyzing observed variant metrics: {e}")
        return defaultdict(list), defaultdict(list)

def track_nova_read_usage(nova_variants) -> Dict[str, List[int]]:
    """Track which nova reads appear in multiple variant calls.
    
    Returns:
        Dict mapping nova read names to list of variant indices where they appear
    """
    read_to_variants = defaultdict(list)
    
    for variant in nova_variants:
        for read_name in variant['nova_read_names']:
            read_to_variants[read_name].append(variant['ID'])
    
    return read_to_variants

def analyze_false_negatives(all_inserted_nova_read_names, read_to_variants) -> Dict[str, Dict]:
    """Identify nova reads that were inserted but not found in VCF.
    
    Args:
        all_inserted_nova_read_names: Set of nova read names that were inserted
        read_to_variants: Dict mapping nova read names to variant indices
        
    Returns:
        Dict with false negative analysis
    """
    # Get reads that appear in variants
    vcf_reads = set(read_to_variants.keys())
    
    # Find missing reads
    missing_reads = all_inserted_nova_read_names - vcf_reads
    
    return {
        'false_negatives': len(missing_reads),
        'missing_read_names': sorted(list(missing_reads)),
        'false_negative_rate': (len(missing_reads) / len(all_inserted_nova_read_names) * 100) if len(all_inserted_nova_read_names) > 0 else 0,
        'filter_status': 'N/A'  # False negatives have no FILTER as they weren't called
    }

def calculate_expected_ins_type_counts(insertions_data) -> Dict[str, int]:
    """Calculate expected counts directly from insertions data."""
    type_counts = Counter()
    
    for insertion in insertions_data:
        insertion_type = insertion.get('insertion_type', 'unknown')
        type_counts[insertion_type] += 1
    
    return dict(type_counts)

def generate_detailed_summary_report(df, nova_variants, categories, expected_ins_type_counts, observed_vars_by_ins_type, 
                                     false_negatives, read_to_variants, mapping_report):
    """Generate a comprehensive summary report using detailed categorization."""
    print("\n" + "="*60)
    print("NOVA VARIANT DETECTION ANALYSIS")
    print("="*60)
    
    total_expected = sum(expected_ins_type_counts.values())

    print(f"\nTotal variants in VCF: {len(df)}")
    print(f"Expected variants with nova reads: {total_expected}")
    print(f"Total variants with nova reads: {len(nova_variants)}")
    
    # Calculate and display all VCF filter breakdown
    filter_col = 'FILTER' if 'FILTER' in df.columns else 'filter'
    if filter_col in df.columns:
        normalized_filters = df[filter_col].apply(normalize_filter_value)
        all_vcf_filter_counts = normalized_filters.value_counts().to_dict()
        
        print("\nALL VCF VARIANTS BY FILTER:")
        for filter_val, count in sorted(all_vcf_filter_counts.items()):
            print(f"   {filter_val}: {count}")
    else:
        print("\nALL VCF VARIANTS BY FILTER:")
        print(f"   UNKNOWN: {len(df)} (no FILTER column found)")
    
    print("\n1. TRUE POSITIVES (single nova with correct mapping):")
    tp_exclusive_vars = set(categories['true_positives']['exclusive'])
    tp_exclusive = len(tp_exclusive_vars)
    tp_shared_vars = set(categories['true_positives']['shared'])
    tp_shared = len(tp_shared_vars)
    tp_total = categories['true_positives']['total']
    
    print(f"   Exclusive (single nova read, unique to variant): {tp_exclusive}")
    print(f"   Shared (single nova read, present in other variants): {tp_shared}")
    print(f"   TOTAL TRUE POSITIVES: {tp_total}")
        
    print("\n2. FALSE POSITIVES:")
    fp_mapping_vars = set(categories['false_positives']['mapping_error'])
    fp_mapping = len(fp_mapping_vars)
    fp_multi_majority_nova_vars = set(categories['false_positives']['multi_read_majority_nova'])
    fp_multi_majority_nova = len(fp_multi_majority_nova_vars)
    fp_multi_minority_nova_vars = set(categories['false_positives']['multi_read_minority_nova'])
    fp_multi_minority_nova = len(fp_multi_minority_nova_vars)
    fp_total = categories['false_positives']['total']
    
    print(f"   Mapping errors (single nova, wrong position): {fp_mapping}")
    print(f"   Multi-read variants (majority nova): {fp_multi_majority_nova}")
    print(f"   Multi-read variants (minority nova): {fp_multi_minority_nova}")
    print(f"   TOTAL FALSE POSITIVES: {fp_total}")

    print("\n3. FALSE NEGATIVES:")
    fn_total = false_negatives['false_negatives']
    fn_rate = false_negatives['false_negative_rate']
    print(f"   Total false negatives: {fn_total}")
    print(f"   False negative rate: {fn_rate}")

    print("\n4A. PERFORMANCE METRICS (nova variants only):")
    print(f"   Recall: {tp_total}/{total_expected} ({tp_total/(tp_total+fn_total)*100:.2f}%)")
    print(f"   Precision: {tp_total}/{tp_total+fp_total} ({tp_total/(tp_total+fp_total)*100:.2f}%)")
    print(f"   F1 score: {2*tp_total/(2*tp_total+fp_total+fn_total):.2f}")
    print(f"   False positives: {fp_total}/{len(nova_variants)} ({fp_total/(len(nova_variants))*100:.2f}%)")

    print("\n4B. PERFORMANCE METRICS (all VCF variants):")
    print(f"   Precision: {tp_total}/{len(df)} ({tp_total/len(df)*100:.2f}%)")
    print(f"   False positives: {len(df)-tp_total}/{len(df)} ({(len(df)-tp_total)/len(df)*100:.2f}%)")

    print("\n5. NOVA READ REUSE ANALYSIS:")
    # Count reads that appear in multiple variants
    reused_reads = sum(1 for read, variants in read_to_variants.items() if len(variants) > 1)
    total_unique_nova_reads = len(read_to_variants)
    
    print(f"   Total unique nova reads: {total_unique_nova_reads}")
    print(f"   Reads appearing in multiple variants: {reused_reads} ({reused_reads/total_unique_nova_reads*100:.1f}%)")
    
    # Distribution of read usage
    usage_counts = Counter(len(variants) for variants in read_to_variants.values())
    print("   Read usage distribution:")
    for count, num_reads in sorted(usage_counts.items()):
        print(f"     Used in {count} variants: {num_reads} reads")
    
    print("\n6. PERFORMANCE BY INSERTION TYPE:")
    for ins_type in sorted(expected_ins_type_counts.keys()):
        expected_count = expected_ins_type_counts[ins_type]
        observed_variants = set(observed_vars_by_ins_type.get(ins_type, []))

        if len(observed_variants) > 0:
            observed_tp_exclusive_count = len(tp_exclusive_vars.intersection(observed_variants))
            observed_tp_shared_count = len(tp_shared_vars.intersection(observed_variants))
            observed_tp_total_count = observed_tp_exclusive_count + observed_tp_shared_count

            observed_fp_mapping_count = len(fp_mapping_vars.intersection(observed_variants))
            observed_fp_multi_majority_nova_count = len(fp_multi_majority_nova_vars.intersection(observed_variants))
            observed_fp_multi_minority_nova_count = len(fp_multi_minority_nova_vars.intersection(observed_variants))
            observed_fp_total_count = observed_fp_mapping_count + observed_fp_multi_majority_nova_count + observed_fp_multi_minority_nova_count

            print("Recall = observed reads / expected reads")
            print("Precision = observed TP reads / (observed TP reads + observed FP reads)")
            print (f"   {ins_type}:")
            print (f"      Read utilization (recall): {len(observed_variants)}/{expected_count} ({len(observed_variants)/expected_count*100:.1f}%)")
            print (f"      Read utilization (precision): {observed_tp_total_count}/{observed_tp_total_count+observed_fp_total_count} ({observed_tp_total_count/(observed_tp_total_count+observed_fp_total_count)*100:.1f}%)")
        else:
            continue

    print("\n7. VARIANT CHARACTERISTICS:")
    print("   By SVTYPE:")
    observed_svtypes, observed_precisions = _calculate_observed_var_metrics(nova_variants)
    for svtype, vars in sorted(observed_svtypes.items()):
        print(f"     {svtype}: {len(vars)}")
    
    print("   By PRECISION:")
    for precision, vars in sorted(observed_precisions.items()):
        print(f"     {precision}: {len(vars)}")

    print("\n8. FILTER BREAKDOWN:")
    filter_counts = Counter(categories['variant_filters'])
    tp_by_filter = dict(categories['true_positives']['by_filter'])
    fp_by_filter = dict(categories['false_positives']['by_filter'])
    
    print("   All nova variants by filter:")
    for filter_val, count in sorted(filter_counts.items()):
        print(f"     {filter_val}: {count}")
    
    print("   True positives by filter:")
    for filter_val, count in sorted(tp_by_filter.items()):
        print(f"     {filter_val}: {count}")
    
    print("   False positives by filter:")
    for filter_val, count in sorted(fp_by_filter.items()):
        print(f"     {filter_val}: {count}")
    
    print(f"   False negatives: {false_negatives['false_negatives']} (FILTER: N/A)")

    print("\n9. MAPPING VERIFICATION:")
    for read_name, status in mapping_report.items():
        print(f"   {read_name}: {status}")

def normalize_filter_value(filter_val):
    """Normalize a single FILTER value to handle lists, None, etc."""
    if isinstance(filter_val, list):
        return ';'.join(str(f) for f in filter_val) if filter_val else 'UNKNOWN'
    elif filter_val is None or filter_val == '':
        return 'UNKNOWN'
    else:
        return str(filter_val)

def generate_json_report(df, nova_variants, categories, expected_ins_type_counts, observed_vars_by_ins_type,
                        false_negatives, read_to_variants, mapping_report) -> dict:
    """Generate a JSON report matching the console output structure.
    
    Args:
        df: VCF dataframe
        nova_variants: List of variants containing nova reads
        categories: Categorized variants (true/false positives)
        expected_ins_type_counts: Expected counts by insertion type
        observed_vars_by_ins_type: Observed variants by insertion type
        false_negatives: False negative analysis results
        read_to_variants: Mapping of reads to variants
        mapping_report: Mapping verification results
        
    Returns:
        Dictionary containing all report data in structured format
    """
    total_expected = sum(expected_ins_type_counts.values())
    
    # Extract category data
    tp_exclusive_vars = set(categories['true_positives']['exclusive'])
    tp_shared_vars = set(categories['true_positives']['shared'])
    tp_total = categories['true_positives']['total']
    
    fp_mapping_vars = set(categories['false_positives']['mapping_error'])
    fp_multi_majority_nova_vars = set(categories['false_positives']['multi_read_majority_nova'])
    fp_multi_minority_nova_vars = set(categories['false_positives']['multi_read_minority_nova'])
    fp_total = categories['false_positives']['total']
    
    fn_total = false_negatives['false_negatives']
    fn_rate = false_negatives['false_negative_rate']
    
    # Calculate nova read reuse metrics
    reused_reads = sum(1 for read, variants in read_to_variants.items() if len(variants) > 1)
    total_unique_nova_reads = len(read_to_variants)
    usage_counts = Counter(len(variants) for variants in read_to_variants.values())
    
    # Calculate variant characteristics
    observed_svtypes, observed_precisions = _calculate_observed_var_metrics(nova_variants)
    
    # Calculate filter breakdowns
    filter_counts = Counter(categories['variant_filters'])
    tp_by_filter = dict(categories['true_positives']['by_filter'])
    fp_by_filter = dict(categories['false_positives']['by_filter'])
    
    # Build performance by insertion type
    performance_by_type = {}
    for ins_type in sorted(expected_ins_type_counts.keys()):
        expected_count = expected_ins_type_counts[ins_type]
        observed_variants = set(observed_vars_by_ins_type.get(ins_type, []))
        
        observed_tp_exclusive_count = len(tp_exclusive_vars.intersection(observed_variants))
        observed_tp_shared_count = len(tp_shared_vars.intersection(observed_variants))
        observed_tp_total_count = observed_tp_exclusive_count + observed_tp_shared_count
        
        observed_fp_mapping_count = len(fp_mapping_vars.intersection(observed_variants))
        observed_fp_multi_majority_nova_count = len(fp_multi_majority_nova_vars.intersection(observed_variants))
        observed_fp_multi_minority_nova_count = len(fp_multi_minority_nova_vars.intersection(observed_variants))
        observed_fp_total_count = observed_fp_mapping_count + observed_fp_multi_majority_nova_count + observed_fp_multi_minority_nova_count
        
        performance_by_type[ins_type] = {
            "expected": expected_count,
            "observed": len(observed_variants),
            "recall": round(len(observed_variants) / expected_count * 100, 2) if expected_count > 0 else 0,
            "precision": round(observed_tp_total_count / (observed_tp_total_count + observed_fp_total_count) * 100, 2)
                        if (observed_tp_total_count + observed_fp_total_count) > 0 else 0,
            "true_positives": {
                "exclusive": observed_tp_exclusive_count,
                "shared": observed_tp_shared_count,
                "total": observed_tp_total_count
            },
            "false_positives": {
                "mapping_errors": observed_fp_mapping_count,
                "multi_read_majority": observed_fp_multi_majority_nova_count,
                "multi_read_minority": observed_fp_multi_minority_nova_count,
                "total": observed_fp_total_count
            }
        }
    
    # Calculate filter breakdown for ALL VCF variants
    filter_col = 'FILTER' if 'FILTER' in df.columns else 'filter'
    if filter_col in df.columns:
        normalized_filters = df[filter_col].apply(normalize_filter_value)
        all_vcf_filter_counts = normalized_filters.value_counts().to_dict()
    else:
        all_vcf_filter_counts = {'UNKNOWN': len(df)}
    
    # Build complete report structure
    report = {
        "summary": {
            "total_variants_in_vcf": len(df),
            "expected_variants_with_nova": total_expected,
            "total_variants_with_nova": len(nova_variants),
            "all_vcf_filters": {
                "total_by_filter": all_vcf_filter_counts
            }
        },
        "true_positives": {
            "exclusive": {
                "count": len(tp_exclusive_vars)
            },
            "shared": {
                "count": len(tp_shared_vars)
            },
            "total": tp_total
        },
        "false_positives": {
            "mapping_errors": {
                "count": len(fp_mapping_vars)
            },
            "multi_read_majority_nova": {
                "count": len(fp_multi_majority_nova_vars)
            },
            "multi_read_minority_nova": {
                "count": len(fp_multi_minority_nova_vars)
            },
            "total": fp_total
        },
        "false_negatives": {
            "count": fn_total,
            "rate_percent": fn_rate
        },
        "performance_metrics": {
            "nova_variants_only": {
                "recall": round(tp_total / (tp_total + fn_total) * 100, 2) if (tp_total + fn_total) > 0 else 0,
                "recall_fraction": f"{tp_total}/{total_expected}",
                "precision": round(tp_total / (tp_total + fp_total) * 100, 2) if (tp_total + fp_total) > 0 else 0,
                "precision_fraction": f"{tp_total}/{tp_total + fp_total}",
                "f1_score": round(2 * tp_total / (2 * tp_total + fp_total + fn_total), 2) if (2 * tp_total + fp_total + fn_total) > 0 else 0,
                "false_positive_rate": round(fp_total / len(nova_variants) * 100, 2) if len(nova_variants) > 0 else 0
            },
            "all_vcf_variants": {
                "precision": round((tp_total / len(df)) * 100, 2) if len(df) > 0 else 0,
                "precision_fraction": f"{tp_total}/{len(df)}",
                "false_positive_rate": round((len(df) - tp_total) / len(df) * 100, 2) if len(df) > 0 else 0
            }
        },
        "nova_read_reuse": {
            "total_unique_reads": total_unique_nova_reads,
            "reads_in_multiple_variants": reused_reads,
            "reuse_rate_percent": round(reused_reads / total_unique_nova_reads * 100, 2) if total_unique_nova_reads > 0 else 0,
            "usage_distribution": {str(k): v for k, v in sorted(usage_counts.items())}
        },
        "performance_by_insertion_type": performance_by_type,
        "variant_characteristics": {
            "by_svtype": {svtype: len(vars) for svtype, vars in sorted(observed_svtypes.items())},
            "by_precision": {str(precision): len(vars) for precision, vars in sorted(observed_precisions.items())}
        },
        "filter_breakdown": {
            "total_by_filter": dict(filter_counts),
            "pass_count": filter_counts.get('PASS', 0),
            "aln_nm_count": filter_counts.get('ALN_NM', 0),
            "true_positives_by_filter": tp_by_filter,
            "false_positives_by_filter": fp_by_filter,
            "false_negatives": {
                "count": fn_total,
                "filter_status": false_negatives.get('filter_status', 'N/A')
            }
        },
        "mapping_verification": mapping_report
    }
    
    return report

def generate_csv_report(json_report: dict, output_file: Path) -> None:
    """Convert JSON report to a single pandas-readable CSV format.
    
    Creates a unified table with performance metrics by insertion type as the main structure,
    with summary metrics included as additional columns.
    
    Args:
        json_report: The JSON report dictionary
        output_file: Path to write the CSV file
    """
    # Extract summary metrics once
    summary = json_report['summary']
    tp = json_report['true_positives']
    fp = json_report['false_positives']
    fn = json_report['false_negatives']
    nova_metrics = json_report['performance_metrics']['nova_variants_only']
    all_vcf_metrics = json_report['performance_metrics']['all_vcf_variants']
    reuse_data = json_report['nova_read_reuse']
    mapping = json_report['mapping_verification']
    filter_data = json_report.get('filter_breakdown', {})
    all_vcf_filter_data = summary.get('all_vcf_filters', {})
    
    # Build main dataframe from performance by insertion type
    rows = []
    for ins_type, metrics in json_report['performance_by_insertion_type'].items():
        row = {
            # Insertion type metrics
            'insertion_type': ins_type,
            'expected': metrics['expected'],
            'observed': metrics['observed'],
            'recall': metrics['recall'],
            'precision': metrics['precision'],
            'tp_exclusive': metrics['true_positives']['exclusive'],
            'tp_shared': metrics['true_positives']['shared'],
            'tp_total': metrics['true_positives']['total'],
            'fp_mapping_errors': metrics['false_positives']['mapping_errors'],
            'fp_multi_read_majority': metrics['false_positives']['multi_read_majority'],
            'fp_multi_read_minority': metrics['false_positives']['multi_read_minority'],
            'fp_total': metrics['false_positives']['total'],
            
            # Add summary metrics (same for all rows)
            'total_variants_in_vcf': summary['total_variants_in_vcf'],
            'expected_variants_with_nova': summary['expected_variants_with_nova'],
            'total_variants_with_nova': summary['total_variants_with_nova'],
            'overall_tp_exclusive': tp['exclusive']['count'],
            'overall_tp_shared': tp['shared']['count'],
            'overall_tp_total': tp['total'],
            'overall_fp_mapping_errors': fp['mapping_errors']['count'],
            'overall_fp_multi_majority': fp['multi_read_majority_nova']['count'],
            'overall_fp_multi_minority': fp['multi_read_minority_nova']['count'],
            'overall_fp_total': fp['total'],
            'overall_fn': fn['count'],
            'overall_fn_rate_percent': fn['rate_percent'],
            'overall_recall': nova_metrics['recall'],
            'overall_precision': nova_metrics['precision'],
            'overall_f1_score': nova_metrics['f1_score'],
            'overall_precision_all_vcf': all_vcf_metrics['precision'],
            'total_unique_nova_reads': reuse_data['total_unique_reads'],
            'reads_in_multiple_variants': reuse_data['reads_in_multiple_variants'],
            'nova_read_reuse_rate_percent': reuse_data['reuse_rate_percent'],
            'mapping_unknown': mapping.get('unknown_mapping', 0),
            'mapping_unknown_original': mapping.get('unknown_original_mapping', 0),
            'mapping_correct': mapping.get('correct_mapping', 0),
            'mapping_incorrect': mapping.get('incorrect_mapping', 0),
            # Add filter columns
            'filter_pass_count': filter_data.get('pass_count', 0),
            'filter_aln_nm_count': filter_data.get('aln_nm_count', 0),
            'tp_pass': filter_data.get('true_positives_by_filter', {}).get('PASS', 0),
            'tp_aln_nm': filter_data.get('true_positives_by_filter', {}).get('ALN_NM', 0),
            'fp_pass': filter_data.get('false_positives_by_filter', {}).get('PASS', 0),
            'fp_aln_nm': filter_data.get('false_positives_by_filter', {}).get('ALN_NM', 0),
            # Add all-VCF filter columns
            'all_vcf_filter_pass_count': all_vcf_filter_data.get('total_by_filter', {}).get('PASS', 0),
            'all_vcf_filter_aln_nm_count': all_vcf_filter_data.get('total_by_filter', {}).get('ALN_NM', 0)
        }
        rows.append(row)
    
    # Create dataframe and save
    df = pd.DataFrame(rows)
    df.to_csv(output_file, index=False)
    
    # Also create a separate summary CSV that's more compact
    summary_file = output_file.parent / f"{output_file.stem}_compact.csv"
    summary_rows = [{
        'total_variants_in_vcf': summary['total_variants_in_vcf'],
        'expected_variants_with_nova': summary['expected_variants_with_nova'],
        'total_variants_with_nova': summary['total_variants_with_nova'],
        'tp_exclusive': tp['exclusive']['count'],
        'tp_shared': tp['shared']['count'],
        'tp_total': tp['total'],
        'fp_mapping_errors': fp['mapping_errors']['count'],
        'fp_multi_majority': fp['multi_read_majority_nova']['count'],
        'fp_multi_minority': fp['multi_read_minority_nova']['count'],
        'fp_total': fp['total'],
        'fn': fn['count'],
        'fn_rate_percent': fn['rate_percent'],
        'recall': nova_metrics['recall'],
        'precision': nova_metrics['precision'],
        'f1_score': nova_metrics['f1_score'],
        'precision_all_vcf': all_vcf_metrics['precision'],
        'total_unique_nova_reads': reuse_data['total_unique_reads'],
        'reads_in_multiple_variants': reuse_data['reads_in_multiple_variants'],
        'nova_read_reuse_rate_percent': reuse_data['reuse_rate_percent'],
        'mapping_unknown': mapping.get('unknown_mapping', 0),
        'mapping_unknown_original': mapping.get('unknown_original_mapping', 0),
        'mapping_correct': mapping.get('correct_mapping', 0),
        'mapping_incorrect': mapping.get('incorrect_mapping', 0),
        'variant_svtype_counts': ', '.join([f"{k}:{v}" for k, v in json_report['variant_characteristics']['by_svtype'].items()]),
        'variant_precision_counts': ', '.join([f"{k}:{v}" for k, v in json_report['variant_characteristics']['by_precision'].items()]),
        # Add filter summary
        'filter_pass_count': filter_data.get('pass_count', 0),
        'filter_aln_nm_count': filter_data.get('aln_nm_count', 0),
        'tp_pass': filter_data.get('true_positives_by_filter', {}).get('PASS', 0),
        'tp_aln_nm': filter_data.get('true_positives_by_filter', {}).get('ALN_NM', 0),
        'fp_pass': filter_data.get('false_positives_by_filter', {}).get('PASS', 0),
        'fp_aln_nm': filter_data.get('false_positives_by_filter', {}).get('ALN_NM', 0),
        'filter_all_counts': ', '.join([f"{k}:{v}" for k, v in filter_data.get('total_by_filter', {}).items()]),
        # Add all-VCF filter summary
        'all_vcf_filter_pass_count': all_vcf_filter_data.get('total_by_filter', {}).get('PASS', 0),
        'all_vcf_filter_aln_nm_count': all_vcf_filter_data.get('total_by_filter', {}).get('ALN_NM', 0),
        'all_vcf_filter_counts': ', '.join([f"{k}:{v}" for k, v in all_vcf_filter_data.get('total_by_filter', {}).items()])
    }]
    
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(summary_file, index=False)

def export_nova_variants_csv(df, nova_variants, categories, output_file):
    """Export all nova variants with identity classification to CSV.
    
    Args:
        df: Original VCF dataframe
        nova_variants: List of variants containing nova reads (with original_index)
        categories: Categorization results from categorize_variants()
        output_file: Path to write the nova variants CSV
    """
    # Extract all TP and FP nova read names from existing categorization
    tp_reads = set(categories['true_positives']['exclusive'] + 
                   categories['true_positives']['shared'])
    
    fp_reads = set(categories['false_positives']['mapping_error'] + 
                   categories['false_positives']['multi_read_majority_nova'] + 
                   categories['false_positives']['multi_read_minority_nova'])
    
    # Build list of nova variant indices and their identities
    nova_indices = []
    identities = []
    
    for variant in nova_variants:
        nova_reads_list = variant['nova_read_names']
        original_idx = variant['original_index']
        
        # Determine identity based on nova reads categorization
        # If any nova read is TP, mark variant as TP; otherwise FP
        has_tp_read = any(read in tp_reads for read in nova_reads_list)
        
        identity = "true_positive" if has_tp_read else "false_positive"
        
        nova_indices.append(original_idx)
        identities.append(identity)
    
    # Filter dataframe to nova variants and add identity column
    nova_df = df.loc[nova_indices].copy()
    nova_df['identity'] = identities
    
    # Save to CSV
    nova_df.to_csv(output_file, index=False)

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
    
    # Create insertion lookup
    print("Creating insertion lookup...")
    with open(insertions_json, 'r') as f:
        insertions_data = json.load(f)
    read_to_insertion, all_inserted_nova_read_names = create_insertion_lookup(insertions_data)

    print("Calculating expected insertion type counts...")
    expected_ins_type_counts = calculate_expected_ins_type_counts(insertions_data)
    
    print("Identifying variants with nova reads...")
    nova_variants = identify_nova_variants(df)

    if not nova_variants:
        print("No variants found with nova reads!")
        return
    
    print("Verifying mapping locations...")
    read_to_mapped_pos = find_mapping_locations(all_inserted_nova_read_names, modified_bam)
    
    mapping_verification, mapping_report = verify_mapping_locations(all_inserted_nova_read_names, 
                                                    read_to_insertion, 
                                                    read_to_mapped_pos,
                                                    position_tolerance=1000)
    
    print("Analyzing insertion types...")
    observed_vars_by_ins_type = calculate_observed_vars_by_ins_type(read_to_mapped_pos, read_to_insertion)
    
    # Track nova read usage
    print("Tracking nova read usage across variants...")
    read_to_variants = track_nova_read_usage(nova_variants)
    
    # Perform categorization
    print("Categorizing variants...")
    categories = categorize_variants(nova_variants, mapping_verification, read_to_variants)
    
    # Analyze false negatives
    print("Analyzing false negatives...")
    false_negatives = analyze_false_negatives(all_inserted_nova_read_names, read_to_variants)
    
    # Generate reports
    generate_detailed_summary_report(df, nova_variants, categories, expected_ins_type_counts, observed_vars_by_ins_type,
                                      false_negatives, read_to_variants, mapping_report)

    summary_json = output_dir / f"{args.output_prefix}_analysis_summary.json"
    json_report = generate_json_report(df, nova_variants, categories, expected_ins_type_counts, 
                                         observed_vars_by_ins_type, false_negatives, read_to_variants, mapping_report)
    
    with open(summary_json, 'w') as f:
        json.dump(json_report, f, indent=2)
        print(f"\nJSON report written to {summary_json}")
    
    # Generate CSV report
    summary_csv = output_dir / f"{args.output_prefix}_analysis_summary.csv"
    compact_csv = output_dir / f"{args.output_prefix}_analysis_summary_compact.csv"
    generate_csv_report(json_report, summary_csv)
    print(f"CSV report written to {summary_csv}")
    print(f"Compact CSV report written to {compact_csv}")
    
    # Export nova variants with identity
    nova_variants_csv = output_dir / f"{args.output_prefix}_analyzed_variants.csv"
    export_nova_variants_csv(df, nova_variants, categories, nova_variants_csv)
    print(f"Nova variants CSV written to {nova_variants_csv}")

if __name__ == "__main__":
    main()