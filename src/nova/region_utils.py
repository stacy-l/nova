"""
Region filtering utilities for genomic coordinate operations using pyranges.

This module provides utilities for BED file operations and genomic overlap
checking to support region-targeted read selection in Nova.
"""

import logging
import random
from typing import List, Optional, Dict, Any, Tuple
from pathlib import Path
import pyranges as pr
import pandas as pd


class RegionFilter:
    """
    Manages genomic regions for filtering operations using pyranges.
    
    Supports both inclusion (target regions) and exclusion (avoidance regions)
    filtering modes for read selection.
    """
    
    def __init__(self, bed_path: Optional[str] = None, regions: Optional[pr.PyRanges] = None):
        """
        Initialize RegionFilter.
        
        Args:
            bed_path: Path to BED file to load regions from
            regions: Pre-loaded PyRanges object
        """
        self.logger = logging.getLogger(__name__)
        self.regions = None
        
        if bed_path is not None:
            self.load_from_bed(bed_path)
        elif regions is not None:
            self.regions = regions
        else:
            # Empty filter
            self.regions = pr.PyRanges()
    
    def load_from_bed(self, bed_path: str) -> None:
        """
        Load regions from a BED file.
        
        Args:
            bed_path: Path to BED file
            
        Raises:
            FileNotFoundError: If BED file doesn't exist
            ValueError: If BED file format is invalid
        """
        bed_file = Path(bed_path)
        
        if not bed_file.exists():
            raise FileNotFoundError(f"BED file not found: {bed_path}")
        
        try:
            self.regions = pr.read_bed(bed_path)
            self.logger.info(f"Loaded {len(self.regions)} regions from {bed_path}")
            
            # Log basic statistics
            if len(self.regions) > 0:
                total_length = (self.regions.End - self.regions.Start).sum()
                chromosomes = len(self.regions.Chromosome.unique())
                self.logger.debug(f"Regions span {chromosomes} chromosomes, total length: {total_length:,} bp")
            
        except Exception as e:
            self.logger.error(f"Failed to parse BED file {bed_path}: {e}")
            raise ValueError(f"Invalid BED file format: {e}")
    
    def overlaps_read(self, chromosome: str, start: int, end: int) -> bool:
        """
        Check if a read interval overlaps with any managed region.
        
        Args:
            chromosome: Chromosome name
            start: Read start position (0-based, inclusive)
            end: Read end position (0-based, exclusive)
            
        Returns:
            True if read overlaps with any region
        """
        if len(self.regions) == 0:
            return False
        
        # Create PyRanges object for the read
        read_interval = pr.from_dict({
            'Chromosome': [chromosome],
            'Start': [start],
            'End': [end]
        })
        
        # Check for overlap
        overlaps = self.regions.overlap(read_interval)
        return len(overlaps) > 0
    
    def contains_point(self, chromosome: str, position: int) -> bool:
        """
        Check if a genomic position is contained within any managed region.
        
        Args:
            chromosome: Chromosome name
            position: Genomic position (0-based)
            
        Returns:
            True if position is within any region
        """
        if len(self.regions) == 0:
            return False
        
        # Create PyRanges object for the position (single base interval)
        point_interval = pr.from_dict({
            'Chromosome': [chromosome],
            'Start': [position],
            'End': [position + 1]
        })
        
        # Check for overlap
        overlaps = self.regions.overlap(point_interval)
        return len(overlaps) > 0
    
    def filter_regions_by_chromosome(self, chromosome: str) -> pr.PyRanges:
        """
        Get all regions for a specific chromosome.
        
        Args:
            chromosome: Chromosome name
            
        Returns:
            PyRanges object containing regions for the chromosome
        """
        if len(self.regions) == 0:
            return pr.PyRanges()
        
        return self.regions[self.regions.Chromosome == chromosome]
    
    def get_chromosome_names(self) -> List[str]:
        """
        Get all chromosome names that have regions.
        
        Returns:
            List of chromosome names
        """
        if len(self.regions) == 0:
            return []
        
        return sorted(self.regions.Chromosome.unique().tolist())
    
    def get_statistics(self) -> Dict[str, Any]:
        """
        Get statistics about the managed regions.
        
        Returns:
            Dictionary with region statistics
        """
        if len(self.regions) == 0:
            return {
                'total_regions': 0,
                'total_length': 0,
                'chromosomes': 0,
                'regions_per_chromosome': {}
            }
        
        total_length = (self.regions.End - self.regions.Start).sum()
        chromosomes = self.regions.Chromosome.unique()
        
        # Count regions per chromosome
        chr_counts = self.regions.df.groupby('Chromosome', observed=True).size().to_dict()
        
        return {
            'total_regions': len(self.regions),
            'total_length': int(total_length),
            'chromosomes': len(chromosomes),
            'regions_per_chromosome': chr_counts
        }
    
    def merge_overlapping(self) -> None:
        """
        Merge overlapping regions within the same chromosome.
        
        This operation modifies the regions in-place.
        """
        if len(self.regions) > 0:
            self.regions = self.regions.merge()
            self.logger.debug(f"Merged overlapping regions, now have {len(self.regions)} regions")
    
    def __len__(self) -> int:
        """Return number of regions managed."""
        return len(self.regions) if self.regions is not None else 0
    
    def __bool__(self) -> bool:
        """Return True if filter contains any regions."""
        return len(self) > 0


def load_multiple_bed_files(bed_paths: List[str]) -> RegionFilter:
    """
    Load and combine regions from multiple BED files.
    
    Args:
        bed_paths: List of paths to BED files
        
    Returns:
        RegionFilter object containing combined regions
        
    Raises:
        FileNotFoundError: If any BED file doesn't exist
        ValueError: If any BED file format is invalid
    """
    logger = logging.getLogger(__name__)
    
    if not bed_paths:
        return RegionFilter()
    
    all_regions = []
    
    for bed_path in bed_paths:
        try:
            regions = pr.read_bed(bed_path)
            all_regions.append(regions)
            logger.debug(f"Loaded {len(regions)} regions from {bed_path}")
        except Exception as e:
            logger.error(f"Failed to load BED file {bed_path}: {e}")
            raise
    
    if all_regions:
        # Concatenate all regions
        combined_regions = pr.concat(all_regions)
        logger.info(f"Combined {len(combined_regions)} regions from {len(bed_paths)} BED files")
        return RegionFilter(regions=combined_regions)
    else:
        return RegionFilter()


def validate_bed_file(bed_path: str) -> tuple[bool, List[str]]:
    """
    Validate a BED file format and return validation results.
    
    Args:
        bed_path: Path to BED file to validate
        
    Returns:
        Tuple of (is_valid, list_of_error_messages)
    """
    errors = []
    
    try:
        bed_file = Path(bed_path)
        if not bed_file.exists():
            errors.append(f"BED file not found: {bed_path}")
            return (False, errors)
        
        # Try to read the BED file
        regions = pr.read_bed(bed_path)
        
        if len(regions) == 0:
            errors.append("BED file contains no valid regions")
        else:
            # Basic sanity checks for issues that would break region filtering
            df = regions.df
            
            # Check for invalid intervals (end <= start) - this would break overlap logic
            if (df['End'] <= df['Start']).any():
                errors.append("BED file contains invalid intervals where end <= start")
                
    except Exception as e:
        errors.append(f"Invalid BED file format - pyranges error: {e}")
    
    return (len(errors) == 0, errors)


def create_exclusion_filter_from_config(exclusion_bed_paths: List[str]) -> Optional[RegionFilter]:
    """
    Create a combined exclusion filter from multiple BED files specified in config.
    
    Args:
        exclusion_bed_paths: List of paths to exclusion BED files
        
    Returns:
        RegionFilter object for exclusion, or None if no valid files
    """
    logger = logging.getLogger(__name__)
    
    if not exclusion_bed_paths:
        return None
    
    valid_paths = []
    for bed_path in exclusion_bed_paths:
        is_valid, errors = validate_bed_file(bed_path)
        if is_valid:
            valid_paths.append(bed_path)
        else:
            logger.warning(f"Skipping invalid exclusion BED file {bed_path}: {errors}")
    
    if valid_paths:
        exclusion_filter = load_multiple_bed_files(valid_paths)
        exclusion_filter.merge_overlapping()  # Optimize for faster lookups
        logger.info(f"Created exclusion filter from {len(valid_paths)} BED files "
                   f"with {len(exclusion_filter)} total regions")
        return exclusion_filter
    else:
        logger.warning("No valid exclusion BED files found")
        return None


class RegionAwareWindowGenerator:
    """
    Generates sampling windows that respect inclusion/exclusion regions.
    
    This class pre-computes valid genomic regions based on BED files and generates
    random sampling windows only within those valid regions, eliminating the need
    for per-read region checking.
    """
    
    def __init__(self, bam_path: str):
        """
        Initialize the window generator.
        
        Args:
            bam_path: Path to BAM file for chromosome information
        """
        self.bam_path = bam_path
        self.logger = logging.getLogger(__name__)
        self._valid_regions_cache = {}  # Cache computed valid regions
        self._genome_regions = None     # Lazy-loaded genome-wide regions
    
    def _load_genome_regions(self) -> pr.PyRanges:
        """
        Load genome-wide regions from BAM file.
        
        Returns:
            PyRanges object covering entire genome
        """
        if self._genome_regions is not None:
            return self._genome_regions
        
        import pysam
        genome_data = {'Chromosome': [], 'Start': [], 'End': []}
        
        with pysam.AlignmentFile(self.bam_path, "rb") as bam:
            for chrom, length in zip(bam.references, bam.lengths):
                genome_data['Chromosome'].append(chrom)
                genome_data['Start'].append(0)
                genome_data['End'].append(length)
        
        self._genome_regions = pr.from_dict(genome_data)
        self.logger.debug(f"Loaded genome regions for {len(self._genome_regions)} chromosomes")
        return self._genome_regions
    
    def _make_cache_key(self, target_regions: Optional[RegionFilter], 
                       exclusion_regions: Optional[RegionFilter]) -> str:
        """
        Create a cache key for the region combination.
        
        Args:
            target_regions: Target/inclusion regions
            exclusion_regions: Exclusion regions
            
        Returns:
            String cache key
        """
        # Handle None and empty regions consistently
        def get_region_key(region_filter):
            if region_filter is None or len(region_filter) == 0:
                return "none"
            else:
                return str(id(region_filter))
        
        target_id = get_region_key(target_regions)
        exclusion_id = get_region_key(exclusion_regions)
        return f"{target_id}_{exclusion_id}"
    
    def get_valid_regions(self, target_regions: Optional[RegionFilter], 
                         exclusion_regions: Optional[RegionFilter]) -> pr.PyRanges:
        """
        Compute valid regions based on inclusion/exclusion filters.
        
        Args:
            target_regions: Optional inclusion regions
            exclusion_regions: Optional exclusion regions
            
        Returns:
            PyRanges object containing all valid regions
        """
        cache_key = self._make_cache_key(target_regions, exclusion_regions)
        
        if cache_key in self._valid_regions_cache:
            self.logger.debug("Using cached valid regions")
            return self._valid_regions_cache[cache_key]
        
        self.logger.info("Computing valid regions from filters")
        
        # Determine base regions
        if target_regions and len(target_regions) > 0:
            # Start with target regions
            valid_regions = target_regions.regions
            self.logger.debug(f"Starting with {len(valid_regions)} target regions")
        else:
            # Start with whole genome
            valid_regions = self._load_genome_regions()
            self.logger.debug(f"Starting with whole genome ({len(valid_regions)} chromosomes)")
        
        # Apply exclusions
        if exclusion_regions and len(exclusion_regions) > 0:
            original_length = (valid_regions.End - valid_regions.Start).sum()
            valid_regions = valid_regions.subtract(exclusion_regions.regions)
            final_length = (valid_regions.End - valid_regions.Start).sum() if len(valid_regions) > 0 else 0
            
            self.logger.info(f"After applying exclusions: {final_length:,}/{original_length:,} bp remain "
                           f"({100*final_length/original_length:.1f}%)")
        
        # Cache the result
        self._valid_regions_cache[cache_key] = valid_regions
        
        # Log statistics
        if len(valid_regions) > 0:
            total_length = (valid_regions.End - valid_regions.Start).sum()
            num_regions = len(valid_regions)
            chromosomes = len(valid_regions.Chromosome.unique())
            self.logger.info(f"Valid regions: {num_regions} regions across {chromosomes} chromosomes, "
                           f"total {total_length:,} bp")
        else:
            self.logger.warning("No valid regions after applying filters!")
        
        return valid_regions
    
    def calculate_required_windows(self, n_reads: int, reads_per_window: int) -> int:
        """
        Calculate the number of windows needed for the requested reads.
        
        Args:
            n_reads: Number of reads to select
            reads_per_window: Number of reads per window
            
        Returns:
            Number of windows required
        """
        import math
        return math.ceil(n_reads / reads_per_window)
    
    def generate_windows_for_regions(self, valid_regions: pr.PyRanges, n_windows: int, 
                                   window_size: int, seed: Optional[int] = None) -> List[Tuple[str, int, int]]:
        """
        Generate sampling windows within valid regions.
        
        Args:
            valid_regions: PyRanges object with valid genomic regions
            n_windows: Number of windows to generate
            window_size: Size of each window in base pairs
            seed: Random seed for reproducibility
            
        Returns:
            List of (chromosome, start, end) tuples
        """
        if len(valid_regions) == 0:
            self.logger.error("Cannot generate windows from empty regions")
            return []
        
        # Set random seed if provided
        if seed is not None:
            random.seed(seed)
        
        # Convert PyRanges to list of regions for easier manipulation
        regions_df = valid_regions.df
        regions_list = []
        for _, row in regions_df.iterrows():
            length = row['End'] - row['Start']
            regions_list.append({
                'chr': row['Chromosome'],
                'start': row['Start'],
                'end': row['End'],
                'length': length
            })
        
        # Calculate total valid length
        total_length = sum(r['length'] for r in regions_list)
        
        # Allocate windows proportionally to each region
        windows_per_region = []
        allocated = 0
        
        for region in regions_list:
            proportion = region['length'] / total_length
            region_windows = int(proportion * n_windows)
            windows_per_region.append(region_windows)
            allocated += region_windows
        
        # Distribute remaining windows to largest regions
        remaining = n_windows - allocated
        if remaining > 0:
            # Sort by region length descending
            sorted_indices = sorted(range(len(regions_list)), 
                                  key=lambda i: regions_list[i]['length'], 
                                  reverse=True)
            
            for i in range(min(remaining, len(sorted_indices))):
                windows_per_region[sorted_indices[i]] += 1
        
        # Generate windows within each region
        all_windows = []
        
        for region, num_windows in zip(regions_list, windows_per_region):
            if num_windows == 0:
                continue
            
            region_length = region['length']
            
            # Handle regions smaller than window size
            if region_length < window_size:
                # Use the entire region as a window
                for _ in range(num_windows):
                    all_windows.append((region['chr'], region['start'], region['end']))
                self.logger.debug(f"Region {region['chr']}:{region['start']}-{region['end']} "
                               f"smaller than window size, using entire region")
            else:
                # Generate random windows within the region
                for _ in range(num_windows):
                    max_start = region['end'] - window_size
                    start = random.randint(region['start'], max_start)
                    end = start + window_size
                    all_windows.append((region['chr'], start, end))
        
        # Shuffle windows to avoid chromosome ordering bias
        random.shuffle(all_windows)
        
        self.logger.info(f"Generated {len(all_windows)} windows from {len(regions_list)} valid regions")
        
        # Log chromosome distribution
        chr_counts = {}
        for chr, _, _ in all_windows:
            chr_counts[chr] = chr_counts.get(chr, 0) + 1
        
        self.logger.debug(f"Window distribution: {dict(sorted(chr_counts.items()))}")
        
        return all_windows