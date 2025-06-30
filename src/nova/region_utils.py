"""
Region filtering utilities for genomic coordinate operations using pyranges.

This module provides utilities for BED file operations and genomic overlap
checking to support region-targeted read selection in Nova.
"""

import logging
from typing import List, Optional, Dict, Any
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