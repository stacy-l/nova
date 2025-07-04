"""
Read selection module for filtering and stratifying reads from BAM files.
"""

import pysam
import random
from typing import List, Dict, Optional, Tuple
import logging
from dataclasses import dataclass

from .region_utils import RegionFilter, RegionAwareWindowGenerator


@dataclass
class ReadMetadata:
    """Metadata for a selected read."""
    read_name: str
    original_chr: str
    original_pos: int
    read_length: int
    mapq: int
    gc_content: float
    soft_clip_ratio: float


@dataclass
class LazyReadReference:
    """Lazy reference to a read that can fetch sequence on-demand."""
    bam_path: str
    read_name: str
    original_chr: str
    original_pos: int
    read_length: int
    mapq: int
    gc_content: float
    soft_clip_ratio: float
    
    def get_metadata(self) -> ReadMetadata:
        """Get ReadMetadata for this read."""
        return ReadMetadata(
            read_name=self.read_name,
            original_chr=self.original_chr,
            original_pos=self.original_pos,
            read_length=self.read_length,
            mapq=self.mapq,
            gc_content=self.gc_content,
            soft_clip_ratio=self.soft_clip_ratio
        )
    
    def get_sequence(self) -> Optional[str]:
        """Fetch the actual read sequence from BAM file on-demand."""
        try:
            with pysam.AlignmentFile(self.bam_path, "rb") as bam:
                # Search for the read by name and position
                for read in bam.fetch(self.original_chr, self.original_pos, self.original_pos + 1):
                    if read.query_name == self.read_name and read.reference_start == self.original_pos:
                        return read.query_sequence
        except Exception:
            pass
        return None


class ReadSelector:
    """
    Select appropriate reads for insertion based on filtering criteria.
    """
    
    def __init__(self, bam_path: str, min_mapq: int = 20, max_soft_clip_ratio: float = 0.1,
                 min_read_length: int = 10000, max_read_length: int = 20000, 
                 reads_per_window: int = 1, target_regions: Optional[RegionFilter] = None,
                 exclusion_regions: Optional[RegionFilter] = None):
        """
        Initialize ReadSelector.
        
        Args:
            bam_path: Path to BAM file
            min_mapq: Minimum mapping quality
            max_soft_clip_ratio: Maximum soft clipping ratio
            min_read_length: Minimum read length
            max_read_length: Maximum read length
            reads_per_window: Number of reads per genomic window (default: 1 for de novo simulation)
            target_regions: Optional RegionFilter for inclusion filtering (only select reads overlapping these regions)
            exclusion_regions: Optional RegionFilter for exclusion filtering (avoid reads overlapping these regions)
        """
        self.bam_path = bam_path
        self.min_mapq = min_mapq
        self.max_soft_clip_ratio = max_soft_clip_ratio
        self.min_read_length = min_read_length
        self.max_read_length = max_read_length
        self.reads_per_window = reads_per_window
        self.target_regions = target_regions
        self.exclusion_regions = exclusion_regions
        self.logger = logging.getLogger(__name__)
        
        # Log region filtering configuration
        if self.target_regions:
            stats = self.target_regions.get_statistics()
            self.logger.info(f"Targeting {stats['total_regions']} regions across {stats['chromosomes']} chromosomes")
        
        if self.exclusion_regions:
            stats = self.exclusion_regions.get_statistics()
            self.logger.info(f"Excluding {stats['total_regions']} regions across {stats['chromosomes']} chromosomes")
        
    def _calculate_gc_content(self, sequence: str) -> float:
        """
        Calculate GC content of a sequence.
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            GC content as a float between 0 and 1
        """
        if not sequence:
            return 0.0
        gc_count = sequence.count('G') + sequence.count('C')
        return gc_count / len(sequence)
    
    def _calculate_soft_clip_ratio(self, read: pysam.AlignedSegment) -> float:
        """
        Calculate soft clipping ratio from CIGAR.
        
        Args:
            read: Aligned read segment
            
        Returns:
            Ratio of soft-clipped bases to total read length
        """
        if not read.cigartuples:
            return 0.0
        
        total_soft_clipped = 0
        for op, length in read.cigartuples:
            if op == 4:  # BAM_CSOFT_CLIP
                total_soft_clipped += length
        
        return total_soft_clipped / read.query_length if read.query_length > 0 else 0.0
    
    def _passes_filters(self, read: pysam.AlignedSegment) -> bool:
        """
        Check if read passes basic filtering criteria (primary alignment, mapq, read length, soft clipping).
        
        Note: Region filtering is now handled by pre-computed windows, not per-read checks.
        """
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            return False
  
        if read.mapping_quality < self.min_mapq:
            return False

        if read.query_length is None or not (self.min_read_length <= read.query_length <= self.max_read_length):
            return False

        soft_clip_ratio = self._calculate_soft_clip_ratio(read)
        if soft_clip_ratio > self.max_soft_clip_ratio:
            return False
        
        return True
    
    def _lazy_select(self, bam: pysam.AlignmentFile, n_reads: int) -> List[LazyReadReference]:
        """
        Read selection method using region-aware window generation.
        
        Pre-generates windows that respect inclusion/exclusion regions, then samples
        reads from these valid windows, eliminating per-read region checks.
        
        Args:
            bam: Open BAM file handle
            n_reads: Number of reads to select
            
        Returns:
            List of LazyReadReference objects
        """
        # Initialize region-aware window generator
        window_generator = RegionAwareWindowGenerator(self.bam_path)
        
        # Get valid regions based on filters (cached operation)
        valid_regions = window_generator.get_valid_regions(
            self.target_regions,
            self.exclusion_regions
        )
        
        # Calculate required number of windows
        n_windows = window_generator.calculate_required_windows(n_reads, self.reads_per_window)
        
        # Generate windows from valid regions
        window_size = 1000000  # 1MB windows
        windows = window_generator.generate_windows_for_regions(
            valid_regions,
            n_windows,
            window_size,
            seed=random.getstate()[1][0] if hasattr(random.getstate()[1], '__getitem__') else None
        )
        
        if not windows:
            self.logger.error("No valid windows could be generated from the specified regions")
            return []
        
        if len(windows) < n_windows:
            self.logger.warning(f"Could only generate {len(windows)} windows, "
                              f"but {n_windows} were requested for {n_reads} reads")
        
        # Sample reads from windows
        selected_reads = []
        windows_used = 0
        
        # Process windows in order (already shuffled by generator)
        for chrom, start, end in windows:
            if len(selected_reads) >= n_reads:
                break
                
            try:
                window_reads = []
                for read in bam.fetch(chrom, start, end):
                    # Only check basic filters, no region checking needed
                    if not self._passes_filters(read):
                        continue
                    
                    lazy_ref = LazyReadReference(
                        bam_path=self.bam_path,
                        read_name=read.query_name,
                        original_chr=read.reference_name,
                        original_pos=read.reference_start,
                        read_length=read.query_length,
                        mapq=read.mapping_quality,
                        gc_content=self._calculate_gc_content(read.query_sequence),
                        soft_clip_ratio=self._calculate_soft_clip_ratio(read)
                    )
                    
                    window_reads.append(lazy_ref)
                    
                    # Stop if we have enough reads from this window
                    if len(window_reads) >= self.reads_per_window:
                        break
                
                # Add reads from this window
                if window_reads:
                    selected_reads.extend(window_reads[:self.reads_per_window])
                    windows_used += 1
                
            except Exception as e:
                self.logger.warning(f"Error fetching reads from {chrom}:{start}-{end}: {e}")
                continue
        
        self.logger.info(f"Region-aware sampling completed: {len(selected_reads)} reads selected "
                        f"from {windows_used} windows")
        
        return selected_reads[:n_reads]  # Ensure we don't return more than requested
    
    def select_reads(self, n_reads: int, use_legacy_methods: bool = False) -> List[LazyReadReference]:
        """
        Select reads from BAM file using unified sampling strategy, returning lazy references.
        
        Args:
            n_reads: Number of reads to select
            
        Returns:
            List of LazyReadReference objects
        """
        with pysam.AlignmentFile(self.bam_path, "rb") as bam:
            selected_reads = self._lazy_select(bam, n_reads)
            
            # Calculate chromosome distribution for logging
            chrom_counts = {}
            for lazy_ref in selected_reads:
                chrom = lazy_ref.original_chr
                chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
            
            self.logger.info(f"Selected {len(selected_reads)} reads from {self.bam_path} sampled across {len(chrom_counts)} chromosomes")
            
            if len(chrom_counts) > 0:
                self.logger.debug(f"Chromosome distribution: {dict(sorted(chrom_counts.items()))}")
            
            return selected_reads