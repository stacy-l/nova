#!/usr/bin/env python3
"""
Large-scale memory usage measurement script for Nova simulation.
Uses 500 variants to test memory optimization at larger scales.
"""

import os
import sys
import json
import time
import resource
import psutil
from pathlib import Path
from typing import Dict, Any

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

from nova.variant_registry import VariantRegistry
from nova.variant_generator import VariantGenerator
from nova.read_selector import ReadSelector
from nova.read_inserter import ReadInserter


def get_memory_usage() -> Dict[str, float]:
    """Get current memory usage in MB."""
    process = psutil.Process()
    memory_info = process.memory_info()
    system_memory = psutil.virtual_memory()
    return {
        'rss_mb': memory_info.rss / 1024 / 1024,  # Resident Set Size
        'vms_mb': memory_info.vms / 1024 / 1024,  # Virtual Memory Size
        'system_used_gb': system_memory.used / 1024 / 1024 / 1024,  # System memory
        'system_percent': system_memory.percent
    }

def print_memory_stats(stage: str, memory_stats: Dict[str, float]):
    """Print formatted memory statistics."""
    print(f"{stage:.<40} RSS: {memory_stats['rss_mb']:>7.1f} MB | "
          f"VMS: {memory_stats['vms_mb']:>7.1f} MB | "
          f"System: {memory_stats['system_used_gb']:>6.1f} GB")


def run_memory_optimized_simulation(bam_file: str, config_file: str, output_dir: Path) -> Dict[str, Any]:
    """Run simulation using memory-optimized methods with memory tracking."""
    print("=== Large-Scale Memory-Optimized Simulation (500 variants) ===")
    
    # Track memory at each stage
    memory_log = []
    
    # Initial memory
    initial_memory = get_memory_usage()
    print_memory_stats("Initial memory", initial_memory)
    memory_log.append(("initial", initial_memory))
    
    # Setup registry and generator
    registry = VariantRegistry()
    generator = VariantGenerator(registry, random_seed=42)
    
    setup_memory = get_memory_usage()
    print_memory_stats("After setup", setup_memory)
    memory_log.append(("setup", setup_memory))
    
    # Load config and generate variants
    with open(config_file, 'r') as f:
        config = json.load(f)
    
    insertion_ids = generator.generate_from_config(config)
    print(f"Generated {len(insertion_ids)} variants")
    
    variant_memory = get_memory_usage()
    print_memory_stats("After variant generation", variant_memory)
    memory_log.append(("variants", variant_memory))
    
    # LAZY read selection
    selector = ReadSelector(bam_file, min_mapq=20, max_soft_clip_ratio=0.1)
    lazy_reads = selector.select_lazy_reads(len(insertion_ids))  # Use all variants
    print(f"Selected {len(lazy_reads)} lazy read references")
    
    selection_memory = get_memory_usage()
    print_memory_stats("After lazy read selection", selection_memory)
    memory_log.append(("lazy_selection", selection_memory))
    
    # STREAMING insertion with progress monitoring
    inserter = ReadInserter(registry, random_seed=42)
    num_to_insert = min(len(lazy_reads), len(insertion_ids))
    
    print(f"Starting streaming insertion for {num_to_insert} pairs...")
    start_time = time.time()
    
    # Use streaming method with direct file output
    output_prefix = output_dir / "nova_large"
    stats = inserter.save_streaming_results(
        lazy_reads[:num_to_insert], 
        insertion_ids[:num_to_insert], 
        str(output_prefix)
    )
    
    end_time = time.time()
    
    insertion_memory = get_memory_usage()
    print_memory_stats("After streaming insertion", insertion_memory)
    memory_log.append(("streaming_insertion", insertion_memory))
    
    # Save registry
    registry_file = output_dir / "nova_large_registry.json"
    registry.save_to_json(str(registry_file))
    
    final_memory = get_memory_usage()
    print_memory_stats("Final memory", final_memory)
    memory_log.append(("final", final_memory))
    
    # Calculate peak memory usage
    peak_rss = max(entry[1]['rss_mb'] for entry in memory_log)
    peak_vms = max(entry[1]['vms_mb'] for entry in memory_log)
    peak_system = max(entry[1]['system_used_gb'] for entry in memory_log)
    system_increase = final_memory['system_used_gb'] - initial_memory['system_used_gb']
    
    print(f"\nProcessing time: {end_time - start_time:.2f} seconds")
    print(f"Peak RSS memory: {peak_rss:.1f} MB")
    print(f"Peak VMS memory: {peak_vms:.1f} MB")
    print(f"Peak system memory: {peak_system:.1f} GB")
    print(f"System memory increase: {system_increase:.1f} GB")
    
    return {
        'method': 'memory_optimized',
        'memory_log': memory_log,
        'peak_rss_mb': peak_rss,
        'peak_vms_mb': peak_vms,
        'peak_system_gb': peak_system,
        'system_memory_increase_gb': system_increase,
        'processing_time': end_time - start_time,
        'stats': stats,
        'variants_generated': len(insertion_ids),
        'reads_selected': len(lazy_reads),
        'successful_insertions': stats['total_insertions']
    }

def run_traditional_simulation(bam_file: str, config_file: str, output_dir: Path) -> Dict[str, Any]:
    """Run simulation using traditional methods with memory tracking."""
    print("\n=== Large-Scale Traditional Simulation (500 variants) ===")
    
    # Track memory at each stage
    memory_log = []
    
    # Initial memory
    initial_memory = get_memory_usage()
    print_memory_stats("Initial memory", initial_memory)
    memory_log.append(("initial", initial_memory))
    
    # Setup registry and generator
    registry = VariantRegistry()
    generator = VariantGenerator(registry, random_seed=42)
    
    setup_memory = get_memory_usage()
    print_memory_stats("After setup", setup_memory)
    memory_log.append(("setup", setup_memory))
    
    # Load config and generate variants
    with open(config_file, 'r') as f:
        config = json.load(f)
    
    insertion_ids = generator.generate_from_config(config)
    print(f"Generated {len(insertion_ids)} variants")
    
    variant_memory = get_memory_usage()
    print_memory_stats("After variant generation", variant_memory)
    memory_log.append(("variants", variant_memory))
    
    # TRADITIONAL read selection (loads full pysam objects)
    selector = ReadSelector(bam_file, min_mapq=20, max_soft_clip_ratio=0.1)
    selected_reads = selector.select_reads(len(insertion_ids))  # Use all variants
    print(f"Selected {len(selected_reads)} full read objects")
    
    selection_memory = get_memory_usage()
    print_memory_stats("After traditional read selection", selection_memory)
    memory_log.append(("traditional_selection", selection_memory))
    
    # BATCH insertion (loads all results into memory)
    inserter = ReadInserter(registry, random_seed=42)
    num_to_insert = min(len(selected_reads), len(insertion_ids))
    
    print(f"Starting batch insertion for {num_to_insert} pairs...")
    start_time = time.time()
    
    insertion_records, modified_sequences, skip_stats = inserter.insert_random_mode(
        selected_reads[:num_to_insert], 
        insertion_ids[:num_to_insert]
    )
    
    batch_memory = get_memory_usage()
    print_memory_stats("After batch insertion", batch_memory)
    memory_log.append(("batch_insertion", batch_memory))
    
    # Save results using traditional methods
    registry_file = output_dir / "nova_large_traditional_registry.json"
    records_file = output_dir / "nova_large_traditional_insertions.json"
    sequences_file = output_dir / "nova_large_traditional_modified_reads.fasta"
    
    registry.save_to_json(str(registry_file))
    inserter.save_insertion_records(insertion_records, str(records_file))
    inserter.save_modified_sequences(modified_sequences, str(sequences_file))
    
    end_time = time.time()
    
    final_memory = get_memory_usage()
    print_memory_stats("Final memory", final_memory)
    memory_log.append(("final", final_memory))
    
    # Calculate peak memory usage
    peak_rss = max(entry[1]['rss_mb'] for entry in memory_log)
    peak_vms = max(entry[1]['vms_mb'] for entry in memory_log)
    peak_system = max(entry[1]['system_used_gb'] for entry in memory_log)
    system_increase = final_memory['system_used_gb'] - initial_memory['system_used_gb']
    
    print(f"\nProcessing time: {end_time - start_time:.2f} seconds")
    print(f"Peak RSS memory: {peak_rss:.1f} MB")
    print(f"Peak VMS memory: {peak_vms:.1f} MB")
    print(f"Peak system memory: {peak_system:.1f} GB")
    print(f"System memory increase: {system_increase:.1f} GB")
    
    return {
        'method': 'traditional',
        'memory_log': memory_log,
        'peak_rss_mb': peak_rss,
        'peak_vms_mb': peak_vms,
        'peak_system_gb': peak_system,
        'system_memory_increase_gb': system_increase,
        'processing_time': end_time - start_time,
        'variants_generated': len(insertion_ids),
        'reads_selected': len(selected_reads),
        'successful_insertions': len(insertion_records)
    }


def main():
    """Main function to run large-scale memory comparison."""
    print("Nova Large-Scale Memory Usage Measurement Script")
    print("=" * 60)
    
    # Configuration - using large config with 500 variants
    bam_file = 'tests/test_data/test_reads.bam'
    config_file = 'tests/test_data/test_config_large.json'
    output_dir = Path('test_output')
    
    # Verify files exist
    if not Path(bam_file).exists():
        print(f"ERROR: BAM file {bam_file} not found!")
        print("This script requires real BAM files for accurate memory measurement.")
        print("Please ensure test_reads.bam is available in tests/test_data/")
        return 1
    
    if not Path(config_file).exists():
        print(f"ERROR: Config file {config_file} not found!")
        return 1
    
    # Create output directory
    output_dir.mkdir(exist_ok=True)
    
    print(f"BAM file: {bam_file}")
    print(f"Config file: {config_file}")
    print(f"Output directory: {output_dir}")
    print(f"Scale: 500 variants (large test)")
    print()
    
    # Run both simulations for comparison
    memory_optimized_results = run_memory_optimized_simulation(bam_file, config_file, output_dir)
    traditional_results = run_traditional_simulation(bam_file, config_file, output_dir)
    
    # Compare results
    print("\n" + "=" * 80)
    print("LARGE-SCALE MEMORY USAGE COMPARISON (500 variants)")
    print("=" * 80)
    
    print(f"{'Method':<20} {'Peak RSS (MB)':<15} {'Peak VMS (MB)':<15} {'Peak Sys (GB)':<15} {'Time (s)':<10}")
    print("-" * 80)
    print(f"{'Memory-Optimized':<20} {memory_optimized_results['peak_rss_mb']:<15.1f} "
          f"{memory_optimized_results['peak_vms_mb']:<15.1f} "
          f"{memory_optimized_results['peak_system_gb']:<15.1f} "
          f"{memory_optimized_results['processing_time']:<10.2f}")
    print(f"{'Traditional':<20} {traditional_results['peak_rss_mb']:<15.1f} "
          f"{traditional_results['peak_vms_mb']:<15.1f} "
          f"{traditional_results['peak_system_gb']:<15.1f} "
          f"{traditional_results['processing_time']:<10.2f}")
    
    # Calculate memory savings
    rss_savings = ((traditional_results['peak_rss_mb'] - memory_optimized_results['peak_rss_mb']) 
                   / traditional_results['peak_rss_mb'] * 100)
    vms_savings = ((traditional_results['peak_vms_mb'] - memory_optimized_results['peak_vms_mb']) 
                   / traditional_results['peak_vms_mb'] * 100)
    system_savings = ((traditional_results['peak_system_gb'] - memory_optimized_results['peak_system_gb']) 
                      / traditional_results['peak_system_gb'] * 100)
    
    print(f"\nMemory Savings at 500-variant scale:")
    print(f"  RSS: {rss_savings:.1f}% {'reduction' if rss_savings > 0 else 'increase'}")
    print(f"  VMS: {vms_savings:.1f}% {'reduction' if vms_savings > 0 else 'increase'}")
    print(f"  System: {system_savings:.1f}% {'reduction' if system_savings > 0 else 'increase'}")
    
    # Performance comparison
    performance_ratio = memory_optimized_results['processing_time'] / traditional_results['processing_time']
    print(f"\nPerformance Trade-off:")
    print(f"  Memory-optimized is {performance_ratio:.1f}x slower than traditional")
    
    # Verify results consistency
    print(f"\nResults Verification:")
    print(f"  Memory-optimized: {memory_optimized_results['successful_insertions']} successful insertions")
    print(f"  Traditional: {traditional_results['successful_insertions']} successful insertions")
    
    # Save detailed results
    results_file = output_dir / "memory_large_comparison_results.json"
    with open(results_file, 'w') as f:
        json.dump({
            'scale': '500_variants',
            'memory_optimized': memory_optimized_results,
            'traditional': traditional_results,
            'comparison': {
                'rss_savings_percent': rss_savings,
                'vms_savings_percent': vms_savings,
                'system_savings_percent': system_savings,
                'performance_ratio': performance_ratio
            }
        }, f, indent=2)
    
    print(f"\nDetailed results saved to: {results_file}")
    print(f"Output files saved to: {output_dir}/nova_large_*")
    
    return 0


if __name__ == "__main__":
    exit(main())