"""
Command line interface for Nova variant insertion simulator.
"""

import click
import json
import logging
import random
import sys
from pathlib import Path
from typing import Optional, Dict, Any

from .variant_registry import VariantRegistry
from .variant_generator import VariantGenerator
from .read_selector import ReadSelector
from .read_inserter import ReadInserter
from .region_utils import RegionFilter, create_exclusion_filter_from_config


def setup_logging(verbose: bool = False) -> None:
    """
    Configure logging for the application.
    
    Args:
        verbose: Enable debug level logging if True
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

def _validate_config(config_file, logger):
    """
    Validate configuration file format and contents.
    
    CONFIG_FILE: Path to JSON configuration file to validate
    """
    try:
        logger.info(f"Loading configuration from {config_file}")
        with open(config_file, 'r') as f:
            config = json.load(f)
    except Exception as e:
        logger.error(f"Failed to load variant config file: {e}")
        sys.exit(1)
    try:
        registry = VariantRegistry()
        generator = VariantGenerator(registry)
        
        errors = generator.validate_config(config)
        
        if errors:
            logger.error("Configuration validation failed:")
            for error in errors:
                logger.error(f"  - {error}")
            sys.exit(1)
        else:
            logger.info("Configuration is valid")
            total_variants = 0
            
            if 'random' in config:
                random_configs = generator._normalize_to_list(config['random'])
                for random_config in random_configs:
                    total_variants += random_config['n']
            
            if 'simple' in config:
                simple_configs = generator._normalize_to_list(config['simple'])
                for simple_config in simple_configs:
                    total_variants += simple_config['n']
            
            if 'predefined' in config:
                predefined_configs = generator._normalize_to_list(config['predefined'])
                for pred_config in predefined_configs:
                    for type_config in pred_config.values():
                        total_variants += sum(type_config['spec'].values())
            
            logger.info(f"Configuration would generate {total_variants} total variants")        
    except Exception as e:
        logger.error(f"Failed to validate configuration: {e}")
        sys.exit(1)


@click.group()
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose logging')
@click.pass_context
def cli(ctx, verbose):
    """nova: de novo variant insertion simulator."""
    ctx.ensure_object(dict)
    ctx.obj['verbose'] = verbose
    setup_logging(verbose) 


@cli.command()
@click.argument('bam_file', type=click.Path(exists=True))
@click.argument('config_file', type=click.Path(exists=True))
@click.option('--output-dir', '-o', default='.', type=click.Path(),
              help='Output directory for results')
@click.option('--output-prefix', '-p', default='nova_sim',
              help='Prefix for output files')
@click.option('--min-mapq', default=20, type=int,
              help='Minimum mapping quality (default: 20)')
@click.option('--max-soft-clip-ratio', default=0.1, type=float,
              help='Maximum soft clipping ratio (default: 0.1)')
@click.option('--min-read-length', default=10000, type=int,
              help='Minimum read length (default: 10000)')
@click.option('--max-read-length', default=20000, type=int,
              help='Maximum read length (default: 20000)')
@click.option('--min-distance-from-ends', default=1000, type=int,
              help='Minimum distance from read ends for insertion (default: 1000)')
@click.option('--random-seed', type=int,
              help='Random seed for reproducibility')
@click.option('--max-reads-per-window', default=1, type=int,
              help='Maximum reads per genomic window (default: 1 for de novo simulation)')
@click.option('--disable-exclusion-regions', is_flag=True,
              help='Disable exclusion region filtering (overrides exclusion_regions from config)')
@click.pass_context
def simulate(ctx, bam_file, config_file, output_dir, output_prefix,
             min_mapq, max_soft_clip_ratio, min_read_length, max_read_length,
             min_distance_from_ends, random_seed, max_reads_per_window,
             disable_exclusion_regions):
    """
    Run de novo variant insertion simulation.
    Produces four output files:
    - *_insertions.json: Insertion records
    - *_modified_reads.fasta: Modified reads
    - *_registry.json: Sequence registry
    - *_statistics.json: Simulation statistics

    Positional arguments:
    bam_file: Path to input BAM file
    config_file: Path to JSON configuration file (determines variants to simulate)
    """
    logger = logging.getLogger(__name__)

    # Validate configuration before running - will exit if invalid
    _validate_config(config_file, logger)

    try:
        with open(config_file, 'r') as f:
            config = json.load(f)

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        registry = VariantRegistry()
        generator = VariantGenerator(registry, random_seed)
        
        # Setup exclusion regions from config
        exclusion_filter = None
        config_exclusions = config.get('exclusion_regions', [])
        if not disable_exclusion_regions and config_exclusions:
            logger.info(f"Loading exclusion regions from {len(config_exclusions)} BED files")
            exclusion_filter = create_exclusion_filter_from_config(config_exclusions)
        
        logger.info("Generating insertion sequences grouped by target regions")
        region_groups = generator.generate_from_config_with_regions(config)
        
        if not region_groups:
            logger.error("No insertion sequences generated")
            sys.exit(1)
        
        total_insertions = sum(len(ids) for ids in region_groups.values())
        logger.info(f"Generated {total_insertions} insertion sequences across {len(region_groups)} region groups")
        
        # Process each region group separately
        all_reads_with_metadata = []
        all_insertion_ids = []
        
        for region_bed_path, insertion_ids in region_groups.items():
            logger.info(f"Processing {len(insertion_ids)} insertions for region group: {region_bed_path}")
            
            # Setup target regions filter
            target_filter = None
            if region_bed_path != 'global':
                try:
                    target_filter = RegionFilter(bed_path=region_bed_path)
                    logger.info(f"Loaded target regions from {region_bed_path}")
                except Exception as e:
                    logger.error(f"Failed to load target regions from {region_bed_path}: {e}")
                    sys.exit(1)
            
            # Create ReadSelector for this region group
            read_selector = ReadSelector(
                bam_file, min_mapq, max_soft_clip_ratio,
                min_read_length, max_read_length, max_reads_per_window,
                target_regions=target_filter, exclusion_regions=exclusion_filter
            )
            
            # Select reads for this group
            logger.info(f"Selecting {len(insertion_ids)} reads for region group {region_bed_path}")
            reads_with_metadata = read_selector.select_reads(len(insertion_ids))
            
            if len(reads_with_metadata) < len(insertion_ids):
                logger.warning(f"Only found {len(reads_with_metadata)} suitable reads for {region_bed_path}, "
                             f"but need {len(insertion_ids)}. Truncating insertion list.")
                insertion_ids = insertion_ids[:len(reads_with_metadata)]
            
            all_reads_with_metadata.extend(reads_with_metadata)
            all_insertion_ids.extend(insertion_ids)
        
        if not all_insertion_ids:
            logger.error("No reads could be selected for any insertion sequences")
            sys.exit(1)
        
        logger.info(f"Total reads selected: {len(all_reads_with_metadata)} for {len(all_insertion_ids)} insertions")
        
        # Shuffle insertion IDs to randomize distribution of insertion types across chromosomes
        if random_seed is not None:
            random.seed(random_seed)
        random.shuffle(all_insertion_ids)
        
        records_file = output_path / f"{output_prefix}_insertions.json"
        sequences_file = output_path / f"{output_prefix}_modified_reads.fasta"
        registry_file = output_path / f"{output_prefix}_registry.json"
        stats_file = output_path / f"{output_prefix}_statistics.json"
        read_inserter = ReadInserter(registry, min_distance_from_ends, random_seed)
        logger.info("Inserting sequences into reads")
        insertion_stats = read_inserter.insert_streaming(
            all_reads_with_metadata, all_insertion_ids, str(records_file), str(sequences_file)
        ) # insertion records and sequences are written to files directly
        
        # read_inserter.save_insertion_records(insertion_records, str(records_file))
        # read_inserter.save_modified_sequences(modified_sequences, str(sequences_file))
        registry.save_to_json(str(registry_file))
        
        # insertion_stats = read_inserter.get_insertion_statistics(insertion_records)
        
        all_stats = {
            'input_parameters': {
                'bam_file': bam_file,
                'config_file': config_file,
                'total_insertions_requested': total_insertions,
                'reads_selected': len(all_reads_with_metadata),
                'insertions_completed': insertion_stats['total_insertions'],
                'exclusion_regions_used': config_exclusions if config_exclusions and not disable_exclusion_regions else [],
                'exclusion_regions_disabled': disable_exclusion_regions,
                'region_groups': {k: len(v) for k, v in region_groups.items()}
            },
            'insertion_statistics': insertion_stats
        }
        
        with open(stats_file, 'w') as f:
            json.dump(all_stats, f, indent=2)
        
        logger.info(f"Simulation completed successfully!")
        logger.info(f"Results saved to {output_dir}:")
        logger.info(f"  - Insertion records: {records_file}")
        logger.info(f"  - Modified reads: {sequences_file}")
        logger.info(f"  - Sequence registry: {registry_file}")
        logger.info(f"  - Statistics: {stats_file}")
        
    except Exception as e:
        logger.error(f"Simulation failed: {e}")
        if ctx.obj.get('verbose'):
            raise
        sys.exit(1)


@cli.command()
@click.argument('config_file', type=click.Path(exists=True))
def validate_config(config_file):
    """
    cli wrapper around _validate_config(), for convenience.
    """
    logger = logging.getLogger(__name__)
    _validate_config(config_file, logger)

def main():
    cli()

if __name__ == '__main__':
    main()