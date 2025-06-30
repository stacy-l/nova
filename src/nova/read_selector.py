"""
Read selection module for filtering and stratifying reads from BAM files.
"""

import pysam
import random
from typing import List, Dict, Optional, Tuple
import logging
from dataclasses import dataclass

from .region_utils import RegionFilter


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
                 max_reads_per_window: int = 1, target_regions: Optional[RegionFilter] = None,
                 exclusion_regions: Optional[RegionFilter] = None):
        """
        Initialize ReadSelector.
        
        Args:
            bam_path: Path to BAM file
            min_mapq: Minimum mapping quality
            max_soft_clip_ratio: Maximum soft clipping ratio
            min_read_length: Minimum read length
            max_read_length: Maximum read length
            max_reads_per_window: Maximum reads per genomic window (default: 1 for de novo simulation)
            target_regions: Optional RegionFilter for inclusion filtering (only select reads overlapping these regions)
            exclusion_regions: Optional RegionFilter for exclusion filtering (avoid reads overlapping these regions)
        """
        self.bam_path = bam_path
        self.min_mapq = min_mapq
        self.max_soft_clip_ratio = max_soft_clip_ratio
        self.min_read_length = min_read_length
        self.max_read_length = max_read_length
        self.max_reads_per_window = max_reads_per_window
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
    
    def _passes_region_filters(self, read: pysam.AlignedSegment) -> bool:
        """
        Check if read passes region filtering criteria.
        
        Args:
            read: Aligned read segment
            
        Returns:
            True if read passes region filters (or no region filters are set)
        """
        # Get read coordinates
        read_chr = read.reference_name
        read_start = read.reference_start
        read_end = read.reference_end
        
        # Skip reads with invalid coordinates
        if read_chr is None or read_start is None or read_end is None:
            return False
        
        # Apply exclusion filters first (these take precedence)
        if self.exclusion_regions:
            if self.exclusion_regions.overlaps_read(read_chr, read_start, read_end):
                return False
        
        # Apply target region filters if specified
        if self.target_regions:
            if not self.target_regions.overlaps_read(read_chr, read_start, read_end):
                return False
        
        return True
    
    def _passes_filters(self, read: pysam.AlignedSegment) -> bool:
        """
        Check if read passes all filtering criteria (primary alignment, mapq, read length, soft clipping, regions)
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
        
        # Apply region filtering
        if not self._passes_region_filters(read):
            return False
        
        return True
    
    def _get_chromosome_sampling_regions(self, bam: pysam.AlignmentFile, 
                                       window_size: int = 1000000) -> List[Tuple[str, int, int]]:
        """
        Get random sampling regions across all chromosomes proportional to their lengths.
        
        Args:
            bam: Open BAM file handle
            window_size: Size of sampling windows
            
        Returns:
            List of tuples (chromosome, start, end) for random sampling
        """
        regions = []
        
        index_stats = bam.get_index_statistics()
        
        for stat in index_stats:
            chrom = stat.contig
            chrom_length = bam.get_reference_length(chrom)
            
            if chrom_length is None:
                continue
                
            num_windows = max(1, chrom_length // window_size)
            
            for i in range(num_windows):
                start = random.randint(0, max(0, chrom_length - window_size))
                end = min(start + window_size, chrom_length)
                regions.append((chrom, start, end))
        
        return regions

    def _get_chromosome_proportional_targets(self, bam: pysam.AlignmentFile, n_reads: int) -> Dict[str, int]:
        """
        Calculate target number of reads per chromosome proportional to chromosome length.
        
        Args:
            bam: Open BAM file handle
            n_reads: Total number of reads to select
            
        Returns:
            Dictionary mapping chromosome names to target read counts
        """
        index_stats = bam.get_index_statistics()
        
        # Calculate total genome length
        total_length = 0
        chrom_lengths = {}
        
        for stat in index_stats:
            chrom = stat.contig
            chrom_length = bam.get_reference_length(chrom)
            if chrom_length is not None:
                chrom_lengths[chrom] = chrom_length
                total_length += chrom_length
        
        # Allocate reads proportionally
        targets = {}
        allocated = 0
        
        for chrom, length in chrom_lengths.items():
            target = int((length / total_length) * n_reads)
            targets[chrom] = target
            allocated += target
        
        # Distribute remaining reads to largest chromosomes
        remaining = n_reads - allocated
        largest_chroms = sorted(chrom_lengths.keys(), key=lambda x: chrom_lengths[x], reverse=True)
        
        for i in range(remaining):
            if i < len(largest_chroms):
                targets[largest_chroms[i]] += 1
        
        return {k: v for k, v in targets.items() if v > 0}

    def _select_reads_unified(self, bam: pysam.AlignmentFile, n_reads: int) -> List[Tuple[pysam.AlignedSegment, ReadMetadata]]:
        """
        Unified read selection method combining proportional allocation with window limits.
        
        Uses chromosome-proportional allocation for fair distribution while strictly enforcing
        max_reads_per_window to prevent genomic clustering that leads to false positives.
        
        Args:
            bam: Open BAM file handle
            n_reads: Number of reads to select
            
        Returns:
            List of tuples containing (read, metadata)
        """
        selected_reads = []
        
        # Use proportional allocation across chromosomes
        chromosome_targets = self._get_chromosome_proportional_targets(bam, n_reads)
        
        self.logger.info(f"Using unified sampling for {n_reads} reads across {len(chromosome_targets)} chromosomes "
                        f"(max {self.max_reads_per_window} reads per window)")
        
        for chrom, target_reads in chromosome_targets.items():
            if target_reads == 0:
                continue
                
            # Get chromosome length for window creation
            chrom_length = bam.get_reference_length(chrom)
            if chrom_length is None:
                continue
            
            # Create sampling windows
            window_size = 1000000  # 1MB windows
            num_windows_needed = max(1, (target_reads + self.max_reads_per_window - 1) // self.max_reads_per_window)
            total_possible_windows = max(1, chrom_length // window_size)
            
            # Create more windows than we need for better randomization
            num_windows_to_create = min(total_possible_windows, num_windows_needed * 3)
            
            # Generate random windows for this chromosome
            chrom_windows = []
            for _ in range(num_windows_to_create):
                start = random.randint(0, max(0, chrom_length - window_size))
                end = min(start + window_size, chrom_length)
                chrom_windows.append((chrom, start, end))
            
            # Shuffle windows for random sampling
            random.shuffle(chrom_windows)
            
            # Collect reads from this chromosome
            chrom_reads = []
            windows_tried = 0
            max_windows_to_try = len(chrom_windows) * 2
            
            window_idx = 0
            while len(chrom_reads) < target_reads and windows_tried < max_windows_to_try:
                if window_idx >= len(chrom_windows):
                    # Reshuffle and restart if we've tried all windows
                    random.shuffle(chrom_windows)
                    window_idx = 0
                
                _, start, end = chrom_windows[window_idx]
                windows_tried += 1
                window_idx += 1
                
                try:
                    window_reads = 0
                    for read in bam.fetch(chrom, start, end):
                        # Stop if we've reached chromosome target or window limit
                        if len(chrom_reads) >= target_reads or window_reads >= self.max_reads_per_window:
                            break
                        
                        if not self._passes_filters(read):
                            continue
                        
                        metadata = ReadMetadata(
                            read_name=read.query_name,
                            original_chr=read.reference_name,
                            original_pos=read.reference_start,
                            read_length=read.query_length,
                            mapq=read.mapping_quality,
                            gc_content=self._calculate_gc_content(read.query_sequence),
                            soft_clip_ratio=self._calculate_soft_clip_ratio(read)
                        )
                        
                        chrom_reads.append((read, metadata))
                        window_reads += 1
                        
                except Exception as e:
                    self.logger.warning(f"Skipping region {chrom}:{start}-{end} due to error: {e}")
                    continue
            
            selected_reads.extend(chrom_reads)
            self.logger.debug(f"Selected {len(chrom_reads)} reads from {chrom} (target: {target_reads})")
        
        self.logger.info(f"Unified sampling completed: {len(selected_reads)} reads selected")
        return selected_reads

    def _select_lazy_reads_unified(self, bam: pysam.AlignmentFile, n_reads: int) -> List[LazyReadReference]:
        """
        Unified read selection method returning lazy references.
        
        Args:
            bam: Open BAM file handle
            n_reads: Number of reads to select
            
        Returns:
            List of LazyReadReference objects
        """
        selected_reads = []
        
        # Use proportional allocation across chromosomes
        chromosome_targets = self._get_chromosome_proportional_targets(bam, n_reads)
        
        self.logger.info(f"Using unified sampling for {n_reads} reads across {len(chromosome_targets)} chromosomes "
                        f"(max {self.max_reads_per_window} reads per window)")
        
        for chrom, target_reads in chromosome_targets.items():
            if target_reads == 0:
                continue
                
            # Get chromosome length for window creation
            chrom_length = bam.get_reference_length(chrom)
            if chrom_length is None:
                continue
            
            # Create sampling windows
            window_size = 1000000  # 1MB windows
            num_windows_needed = max(1, (target_reads + self.max_reads_per_window - 1) // self.max_reads_per_window)
            total_possible_windows = max(1, chrom_length // window_size)
            
            # Create more windows than we need for better randomization
            num_windows_to_create = min(total_possible_windows, num_windows_needed * 3)
            
            # Generate random windows for this chromosome
            chrom_windows = []
            for _ in range(num_windows_to_create):
                start = random.randint(0, max(0, chrom_length - window_size))
                end = min(start + window_size, chrom_length)
                chrom_windows.append((chrom, start, end))
            
            # Shuffle windows for random sampling
            random.shuffle(chrom_windows)
            
            # Collect reads from this chromosome
            chrom_reads = []
            windows_tried = 0
            max_windows_to_try = len(chrom_windows) * 2
            
            window_idx = 0
            while len(chrom_reads) < target_reads and windows_tried < max_windows_to_try:
                if window_idx >= len(chrom_windows):
                    # Reshuffle and restart if we've tried all windows
                    random.shuffle(chrom_windows)
                    window_idx = 0
                
                _, start, end = chrom_windows[window_idx]
                windows_tried += 1
                window_idx += 1
                
                try:
                    window_reads = 0
                    for read in bam.fetch(chrom, start, end):
                        # Stop if we've reached chromosome target or window limit
                        if len(chrom_reads) >= target_reads or window_reads >= self.max_reads_per_window:
                            break
                        
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
                        
                        chrom_reads.append(lazy_ref)
                        window_reads += 1
                        
                except Exception as e:
                    self.logger.warning(f"Skipping region {chrom}:{start}-{end} due to error: {e}")
                    continue
            
            selected_reads.extend(chrom_reads)
            self.logger.debug(f"Selected {len(chrom_reads)} reads from {chrom} (target: {target_reads})")
        
        self.logger.info(f"Unified sampling completed: {len(selected_reads)} reads selected")
        return selected_reads

    def _select_reads_with_proportional_sampling(self, bam: pysam.AlignmentFile, n_reads: int) -> List[Tuple[pysam.AlignedSegment, ReadMetadata]]:
        """
        Select reads using chromosome-proportional sampling for large simulations.
        
        Args:
            bam: Open BAM file handle
            n_reads: Number of reads to select
            
        Returns:
            List of tuples containing (read, metadata)
        """
        selected_reads = []
        chromosome_targets = self._get_chromosome_proportional_targets(bam, n_reads)
        
        self.logger.info(f"Using proportional sampling for {n_reads} reads across {len(chromosome_targets)} chromosomes")
        
        for chrom, target_reads in chromosome_targets.items():
            if target_reads == 0:
                continue
                
            # Get sampling regions for this chromosome
            chrom_length = bam.get_reference_length(chrom)
            if chrom_length is None:
                continue
                
            # Create windows for this chromosome
            window_size = 1000000
            num_windows = max(1, chrom_length // window_size)
            
            chrom_regions = []
            for i in range(num_windows):
                start = random.randint(0, max(0, chrom_length - window_size))
                end = min(start + window_size, chrom_length)
                chrom_regions.append((chrom, start, end))
            
            # Shuffle regions for this chromosome
            random.shuffle(chrom_regions)
            
            # Collect reads from this chromosome
            chrom_reads = []
            regions_tried = 0
            max_regions_to_try = len(chrom_regions) * 2
            
            region_idx = 0
            while len(chrom_reads) < target_reads and regions_tried < max_regions_to_try:
                if region_idx >= len(chrom_regions):
                    random.shuffle(chrom_regions)
                    region_idx = 0
                
                _, start, end = chrom_regions[region_idx]
                regions_tried += 1
                region_idx += 1
                
                try:
                    for read in bam.fetch(chrom, start, end):
                        if len(chrom_reads) >= target_reads:
                            break
                        
                        if not self._passes_filters(read):
                            continue
                        
                        metadata = ReadMetadata(
                            read_name=read.query_name,
                            original_chr=read.reference_name,
                            original_pos=read.reference_start,
                            read_length=read.query_length,
                            mapq=read.mapping_quality,
                            gc_content=self._calculate_gc_content(read.query_sequence),
                            soft_clip_ratio=self._calculate_soft_clip_ratio(read)
                        )
                        
                        chrom_reads.append((read, metadata))
                        
                except Exception as e:
                    self.logger.warning(f"Skipping region {chrom}:{start}-{end} due to error: {e}")
                    continue
            
            selected_reads.extend(chrom_reads)
            self.logger.debug(f"Selected {len(chrom_reads)} reads from {chrom} (target: {target_reads})")
        
        return selected_reads
    
    def _select_lazy_reads_with_proportional_sampling(self, bam: pysam.AlignmentFile, n_reads: int) -> List[LazyReadReference]:
        """
        Select reads using chromosome-proportional sampling, returning lazy references.
        
        Args:
            bam: Open BAM file handle
            n_reads: Number of reads to select
            
        Returns:
            List of LazyReadReference objects
        """
        selected_reads = []
        chromosome_targets = self._get_chromosome_proportional_targets(bam, n_reads)
        
        self.logger.info(f"Using proportional sampling for {n_reads} reads across {len(chromosome_targets)} chromosomes")
        
        for chrom, target_reads in chromosome_targets.items():
            if target_reads == 0:
                continue
                
            # Get sampling regions for this chromosome
            chrom_length = bam.get_reference_length(chrom)
            if chrom_length is None:
                continue
                
            # Create windows for this chromosome
            window_size = 1000000
            num_windows = max(1, chrom_length // window_size)
            
            chrom_regions = []
            for i in range(num_windows):
                start = random.randint(0, max(0, chrom_length - window_size))
                end = min(start + window_size, chrom_length)
                chrom_regions.append((chrom, start, end))
            
            # Shuffle regions for this chromosome
            random.shuffle(chrom_regions)
            
            # Collect reads from this chromosome
            chrom_reads = []
            regions_tried = 0
            max_regions_to_try = len(chrom_regions) * 2
            
            region_idx = 0
            while len(chrom_reads) < target_reads and regions_tried < max_regions_to_try:
                if region_idx >= len(chrom_regions):
                    random.shuffle(chrom_regions)
                    region_idx = 0
                
                _, start, end = chrom_regions[region_idx]
                regions_tried += 1
                region_idx += 1
                
                try:
                    for read in bam.fetch(chrom, start, end):
                        if len(chrom_reads) >= target_reads:
                            break
                        
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
                        
                        chrom_reads.append(lazy_ref)
                        
                except Exception as e:
                    self.logger.warning(f"Skipping region {chrom}:{start}-{end} due to error: {e}")
                    continue
            
            selected_reads.extend(chrom_reads)
            self.logger.debug(f"Selected {len(chrom_reads)} reads from {chrom} (target: {target_reads})")
        
        return selected_reads

    def _select_reads_with_window_limits(self, bam: pysam.AlignmentFile, n_reads: int) -> List[Tuple[pysam.AlignedSegment, ReadMetadata]]:
        """
        Select reads using window-limited sampling for small simulations.
        
        Args:
            bam: Open BAM file handle
            n_reads: Number of reads to select
            
        Returns:
            List of tuples containing (read, metadata)
        """
        selected_reads = []
        
        # Get all sampling regions and shuffle upfront
        sampling_regions = self._get_chromosome_sampling_regions(bam)
        random.shuffle(sampling_regions)
        
        # Limit reads per window to ensure better distribution
        max_reads_per_window = max(1, n_reads // 10)  # Max 10% of reads per window
        
        self.logger.info(f"Using window-limited sampling for {n_reads} reads (max {max_reads_per_window} per window)")
        
        regions_tried = 0
        max_regions_to_try = len(sampling_regions) * 3
        
        region_idx = 0
        while len(selected_reads) < n_reads and regions_tried < max_regions_to_try:
            if region_idx >= len(sampling_regions):
                # Reshuffle and restart if we've gone through all regions
                random.shuffle(sampling_regions)
                region_idx = 0
            
            chrom, start, end = sampling_regions[region_idx]
            regions_tried += 1
            region_idx += 1
            
            try:
                window_reads = 0
                for read in bam.fetch(chrom, start, end):
                    if len(selected_reads) >= n_reads or window_reads >= max_reads_per_window:
                        break
                    
                    if not self._passes_filters(read):
                        continue
                    
                    metadata = ReadMetadata(
                        read_name=read.query_name,
                        original_chr=read.reference_name,
                        original_pos=read.reference_start,
                        read_length=read.query_length,
                        mapq=read.mapping_quality,
                        gc_content=self._calculate_gc_content(read.query_sequence),
                        soft_clip_ratio=self._calculate_soft_clip_ratio(read)
                    )
                    
                    selected_reads.append((read, metadata))
                    window_reads += 1
                    
            except Exception as e:
                self.logger.warning(f"Skipping region {chrom}:{start}-{end} due to error: {e}")
                continue
        
        return selected_reads
    
    def _select_lazy_reads_with_window_limits(self, bam: pysam.AlignmentFile, n_reads: int) -> List[LazyReadReference]:
        """
        Select reads using window-limited sampling, returning lazy references.
        
        Args:
            bam: Open BAM file handle
            n_reads: Number of reads to select
            
        Returns:
            List of LazyReadReference objects
        """
        selected_reads = []
        
        # Get all sampling regions and shuffle upfront
        sampling_regions = self._get_chromosome_sampling_regions(bam)
        random.shuffle(sampling_regions)
        
        # Limit reads per window to ensure better distribution
        max_reads_per_window = max(1, n_reads // 10)  # Max 10% of reads per window
        
        self.logger.info(f"Using window-limited sampling for {n_reads} reads (max {max_reads_per_window} per window)")
        
        regions_tried = 0
        max_regions_to_try = len(sampling_regions) * 3
        
        region_idx = 0
        while len(selected_reads) < n_reads and regions_tried < max_regions_to_try:
            if region_idx >= len(sampling_regions):
                # Reshuffle and restart if we've gone through all regions
                random.shuffle(sampling_regions)
                region_idx = 0
            
            chrom, start, end = sampling_regions[region_idx]
            regions_tried += 1
            region_idx += 1
            
            try:
                window_reads = 0
                for read in bam.fetch(chrom, start, end):
                    if len(selected_reads) >= n_reads or window_reads >= max_reads_per_window:
                        break
                    
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
                    
                    selected_reads.append(lazy_ref)
                    window_reads += 1
                    
            except Exception as e:
                self.logger.warning(f"Skipping region {chrom}:{start}-{end} due to error: {e}")
                continue
        
        return selected_reads
    
    def select_reads(self, n_reads: int, use_legacy_methods: bool = False) -> List[Tuple[pysam.AlignedSegment, ReadMetadata]]:
        """
        Select reads from BAM file using unified sampling strategy.
        
        By default uses the new unified method that combines proportional allocation
        with strict window limits to prevent false positives. Legacy methods can
        be used for comparison testing.
        
        Args:
            n_reads: Number of reads to select
            use_legacy_methods: If True, use original hybrid strategy for comparison
            
        Returns:
            List of tuples containing (read, metadata)
        """
        with pysam.AlignmentFile(self.bam_path, "rb") as bam:
            if use_legacy_methods:
                # Legacy hybrid approach for comparison
                if n_reads < 500:
                    selected_reads = self._select_reads_with_window_limits(bam, n_reads)
                    strategy = "legacy-window-limited"
                else:
                    selected_reads = self._select_reads_with_proportional_sampling(bam, n_reads)
                    strategy = "legacy-proportional"
            else:
                # New unified approach (default)
                selected_reads = self._select_reads_unified(bam, n_reads)
                strategy = "unified"
            
            # Calculate chromosome distribution for logging
            chrom_counts = {}
            for _, metadata in selected_reads:
                chrom = metadata.original_chr
                chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
            
            self.logger.info(f"Selected {len(selected_reads)} reads from {self.bam_path} "
                           f"using {strategy} sampling across {len(chrom_counts)} chromosomes")
            
            if len(chrom_counts) > 0:
                self.logger.debug(f"Chromosome distribution: {dict(sorted(chrom_counts.items()))}")
            
            return selected_reads
    
    def select_lazy_reads(self, n_reads: int, use_legacy_methods: bool = False) -> List[LazyReadReference]:
        """
        Select reads from BAM file using unified sampling strategy, returning lazy references.
        
        By default uses the new unified method that combines proportional allocation
        with strict window limits to prevent false positives. Legacy methods can
        be used for comparison testing.
        
        Args:
            n_reads: Number of reads to select
            use_legacy_methods: If True, use original hybrid strategy for comparison
            
        Returns:
            List of LazyReadReference objects
        """
        with pysam.AlignmentFile(self.bam_path, "rb") as bam:
            if use_legacy_methods:
                # Legacy hybrid approach for comparison
                if n_reads < 500:
                    selected_reads = self._select_lazy_reads_with_window_limits(bam, n_reads)
                    strategy = "legacy-window-limited"
                else:
                    selected_reads = self._select_lazy_reads_with_proportional_sampling(bam, n_reads)
                    strategy = "legacy-proportional"
            else:
                # New unified approach (default)
                selected_reads = self._select_lazy_reads_unified(bam, n_reads)
                strategy = "unified"
            
            # Calculate chromosome distribution for logging
            chrom_counts = {}
            for lazy_ref in selected_reads:
                chrom = lazy_ref.original_chr
                chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
            
            self.logger.info(f"Selected {len(selected_reads)} reads from {self.bam_path} "
                           f"using {strategy} sampling across {len(chrom_counts)} chromosomes")
            
            if len(chrom_counts) > 0:
                self.logger.debug(f"Chromosome distribution: {dict(sorted(chrom_counts.items()))}")
            
            return selected_reads