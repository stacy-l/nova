#!/usr/bin/env python3
"""
VCF intersection analysis script.
Derived from truvari/consistency.py.

Identifies matching variants across multiple VCFs based on:
CHROM, POS, REF, ALT, SVTYPE, and SVLEN fields.
"""

import gzip
import json
import argparse
import re
import os
from collections import defaultdict, Counter


def parse_info_field(info_str):
    """
    Extract SVTYPE and SVLEN from INFO field.
    
    Args:
        info_str: INFO field string (e.g., "SVTYPE=INS;SVLEN=123;...")
    
    Returns:
        tuple: (svtype, svlen) - empty strings if not found
    """
    svtype_match = re.search(r'SVTYPE=([^;]+)', info_str)
    svlen_match = re.search(r'SVLEN=(-?\d+)', info_str)
    
    svtype = svtype_match.group(1) if svtype_match else ""
    svlen = svlen_match.group(1) if svlen_match else ""
    
    return svtype, svlen


def parse_vcf(vcf_path, ignore_alt):
    """
    Parse VCF file and yield variant keys and IDs.
    
    Args:
        vcf_path: Path to VCF file (can be gzipped)
        ignore_alt: Whether to ignore ALT field in key
    
    Yields:
        tuple: (variant_key, variant_id)
    """
    opener = gzip.open if vcf_path.endswith('.gz') else open
    
    with opener(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            
            chrom, pos, variant_id, ref, alt = fields[:5]
            info = fields[7]
            
            svtype, svlen = parse_info_field(info)
            if ignore_alt:
                key = f"{chrom}:{pos}:{ref}:{svtype}:{svlen}"
            else:
                key = f"{chrom}:{pos}:{ref}:{alt}:{svtype}:{svlen}"
            yield key, variant_id


def read_vcf_files(vcf_paths, ignore_alt=False):
    """
    Read all VCF files and track variant presence across files.
    
    Args:
        vcf_paths: List of VCF file paths
        ignore_alt: Whether to ignore ALT field in matching
    
    Returns:
        tuple: (all_presence, n_calls_per_vcf, variant_ids)
            all_presence: dict mapping variant keys to presence bit flags
            n_calls_per_vcf: list of variant counts per VCF
            variant_ids: dict mapping variant keys to sets of IDs per VCF
    """
    n_vcfs = len(vcf_paths)
    all_presence = defaultdict(lambda: 0)
    n_calls_per_vcf = [0] * n_vcfs
    variant_ids = defaultdict(lambda: [set() for _ in range(n_vcfs)])
    
    for i, vcf_path in enumerate(vcf_paths):
        flag = 1 << (n_vcfs - i - 1)
        seen_in_this_vcf = set()
        
        for key, variant_id in parse_vcf(vcf_path, ignore_alt):
            variant_ids[key][i].add(variant_id)
            if key not in seen_in_this_vcf:
                seen_in_this_vcf.add(key)
                n_calls_per_vcf[i] += 1
                all_presence[key] |= flag
    
    return all_presence, n_calls_per_vcf, variant_ids


def make_report(vcf_paths, all_presence, n_calls_per_vcf):
    """
    Generate intersection report data.
    
    Args:
        vcf_paths: List of VCF file paths
        all_presence: Dict mapping variant keys to presence bit flags
        n_calls_per_vcf: List of variant counts per VCF
    
    Returns:
        dict: Report data structure
    """
    output = {}
    output['vcfs'] = vcf_paths
    total_unique_calls = len(all_presence)
    output['total_calls'] = total_unique_calls
    n_vcfs = len(vcf_paths)
    output['num_vcfs'] = n_vcfs
    
    # VCF call counts
    output['vcf_counts'] = {}
    for i, vcf_path in enumerate(vcf_paths):
        output['vcf_counts'][vcf_path] = n_calls_per_vcf[i]
    
    # Initialize shared counts
    output['shared'] = []
    for i in reversed(range(n_vcfs)):
        output['shared'].append({'vcf_count': i + 1, 'num_calls': 0, 'call_pct': 0})
    
    # Process presence data
    output['detailed'] = []
    all_overlap = Counter(all_presence.values())
    
    for group_flag, ncalls in sorted(all_overlap.items(), key=lambda x: (-x[1], x[0])):
        group_bin = bin(group_flag)[2:].rjust(n_vcfs, '0')
        shared_calls_vcf_count = group_bin.count("1")
        output['shared'][-shared_calls_vcf_count]["num_calls"] += ncalls
        
        tot_pct = ncalls / total_unique_calls if total_unique_calls > 0 else 0
        
        # Calculate per-VCF percentages
        row = {
            'group': group_bin,
            'total': ncalls,
            'total_pct': tot_pct
        }
        
        for i in range(n_vcfs):
            if group_bin[i] == "1" and n_calls_per_vcf[i] > 0:
                row[vcf_paths[i]] = ncalls / n_calls_per_vcf[i]
            else:
                row[vcf_paths[i]] = 0
        
        output['detailed'].append(row)
    
    # Calculate shared percentages
    for i in range(n_vcfs):
        if total_unique_calls > 0:
            output['shared'][i]['call_pct'] = output['shared'][i]['num_calls'] / total_unique_calls
    
    return output


def write_console_report(output):
    """
    Write human-readable console report.
    
    Args:
        output: Report data structure
    """
    print(f"#\n# Total {output['total_calls']} calls across {output['num_vcfs']} VCFs\n#")
    
    print("#File\tNumCalls")
    for vcf_path in output['vcfs']:
        print(f"{vcf_path}\t{output['vcf_counts'][vcf_path]}")
    
    print("#\n# Summary of consistency\n#")
    print("#VCFs\tCalls\tPct")
    for row in output['shared']:
        print(f"{row['vcf_count']}\t{row['num_calls']}\t{row['call_pct'] * 100:.2f}%")
    
    print("#\n# Breakdown of VCFs' consistency\n#")
    print("#Group\tTotal\tTotalPct")
    for row in output['detailed']:
        group = row['group']
        ncalls = row['total']
        tot_pct = row['total_pct'] * 100
        
        print(f"{group}\t{ncalls}\t{tot_pct:.2f}%")


def write_private_variant_files(report, output_prefix, all_presence, variant_ids):
    """
    Write private variant ID files for each VCF.
    Files can be used with bcftools to filter private variants.

    Example:
    bcftools view -i 'ID=@01.txt' input.vcf.gz -o private_variants.vcf
    
    Args:
        report: Report data structure
        output_prefix: Output directory/prefix for files
        all_presence: Dict mapping variant keys to presence bit flags
        variant_ids: Dict mapping variant keys to IDs per VCF
    
    Returns:
        dict: Mapping of group patterns to output file paths
    """
    if output_prefix and os.path.dirname(output_prefix):
        os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
    
    n_vcfs = report['num_vcfs']
    private_files = {}
    
    variant_groups = defaultdict(list)
    
    for variant_key, presence_flag in all_presence.items():
        group_bin = bin(presence_flag)[2:].rjust(n_vcfs, '0')
        
        if group_bin.count('1') == 1:
            variant_groups[group_bin].append(variant_key)
    
    for group_pattern, variant_keys in variant_groups.items():
        vcf_index = group_pattern.index('1')
        
        variant_id_list = []
        for variant_key in variant_keys:
            ids_for_vcf = variant_ids[variant_key][vcf_index]
            variant_id_list.extend(ids_for_vcf)
        
        filename = f"{output_prefix}{group_pattern}.txt"
        with open(filename, 'w') as f:
            for variant_id in sorted(set(variant_id_list)):
                f.write(f"{variant_id}\n")
        
        private_files[group_pattern] = {
            'filename': filename,
            'count': len(set(variant_id_list))
        }
    
    return private_files


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        'vcf_files', 
        nargs='+',
        help='VCF files to analyze for intersections'
    )
    parser.add_argument(
        '--ignore-alt', 
        action='store_true',
        help='Do not match on ALT sequence (ex. compare mutated insertion sequences)'
    )    
    parser.add_argument(
        '--json', 
        action='store_true',
        help='Output results in JSON format'
    )
    parser.add_argument(
        '--report-pids',
        type=str,
        help='Creates private variant ID files and saves them with the given prefix/directory.'
    )
    
    return parser.parse_args()


def main():
    """Main entry point."""
    args = parse_args()
    
    if args.report_pids:
        all_presence, n_calls_per_vcf, variant_ids = read_vcf_files(args.vcf_files, args.ignore_alt)
        report = make_report(args.vcf_files, all_presence, n_calls_per_vcf)
        private_files = write_private_variant_files(report, args.report_pids, all_presence, variant_ids)
    else:
        all_presence, n_calls_per_vcf, _ = read_vcf_files(args.vcf_files, args.ignore_alt)
        report = make_report(args.vcf_files, all_presence, n_calls_per_vcf)
        private_files = {}
    
    if args.json:
        if private_files:
            report['private_variant_files'] = private_files
        print(json.dumps(report, indent=4))
    else:
        write_console_report(report)
        if private_files:
            print("\n# Private variant files created:")
            for group, info in private_files.items():
                print(f"#{group}\t{info['filename']}\t{info['count']} variants")


if __name__ == '__main__':
    main()