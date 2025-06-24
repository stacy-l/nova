"""
Command line interface for Nova variant insertion simulator.
"""

import click
import json
import logging
import sys
from pathlib import Path
from typing import Optional, Dict, Any

from .variant_registry import VariantRegistry
from .variant_generator import VariantGenerator
from .read_selector import ReadSelector
from .read_inserter import ReadInserter


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
@click.pass_context
def simulate(ctx, bam_file, config_file, output_dir, output_prefix,
             min_mapq, max_soft_clip_ratio, min_read_length, max_read_length,
             min_distance_from_ends, random_seed):
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
    
    try:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Loading configuration from {config_file}")
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        registry = VariantRegistry()
        generator = VariantGenerator(registry, random_seed)
        
        # Validate configuration
        validation_errors = generator.validate_config(config)
        if validation_errors:
            logger.error("Configuration validation failed:")
            for error in validation_errors:
                logger.error(f"  - {error}")
            sys.exit(1)
        
        logger.info("Generating insertion sequences")
        insertion_ids = generator.generate_from_config(config)
        
        if not insertion_ids:
            logger.error("No insertion sequences generated")
            sys.exit(1)
        
        total_insertions = len(insertion_ids)
        logger.info(f"Generated {total_insertions} insertion sequences")
        read_selector = ReadSelector(
            bam_file, min_mapq, max_soft_clip_ratio,
            min_read_length, max_read_length
        )
        
        logger.info(f"Selecting {total_insertions} reads from BAM file")
        reads_with_metadata = read_selector.select_reads(total_insertions)
        
        if len(reads_with_metadata) < total_insertions:
            logger.warning(f"Only found {len(reads_with_metadata)} suitable reads, "
                         f"but need {total_insertions}")
            insertion_ids = insertion_ids[:len(reads_with_metadata)]
        
        read_inserter = ReadInserter(registry, min_distance_from_ends, random_seed)
        logger.info("Inserting sequences into reads")
        insertion_records, modified_sequences, skip_stats = read_inserter.insert_random_mode(
            reads_with_metadata, insertion_ids
        )
        
        records_file = output_path / f"{output_prefix}_insertions.json"
        sequences_file = output_path / f"{output_prefix}_modified_reads.fasta"
        registry_file = output_path / f"{output_prefix}_registry.json"
        stats_file = output_path / f"{output_prefix}_statistics.json"
        read_inserter.save_insertion_records(insertion_records, str(records_file))
        read_inserter.save_modified_sequences(modified_sequences, str(sequences_file))
        registry.save_to_json(str(registry_file))
        
        insertion_stats = read_inserter.get_insertion_statistics(insertion_records)
        registry_stats = registry.get_statistics()
        
        all_stats = {
            'input_parameters': {
                'bam_file': bam_file,
                'config_file': config_file,
                'total_insertions_requested': total_insertions,
                'reads_selected': len(reads_with_metadata),
                'insertions_completed': len(insertion_records)
            },
            'insertion_statistics': insertion_stats,
            'registry_statistics': registry_stats,
            'insertion_success_summary': skip_stats
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
    Validate configuration file format and contents.
    
    CONFIG_FILE: Path to JSON configuration file to validate
    """
    logger = logging.getLogger(__name__)
    
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
        
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




def main():
    cli()


if __name__ == '__main__':
    main()