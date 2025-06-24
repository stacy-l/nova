"""
Read selection module for filtering and stratifying reads from BAM files.
"""

import pysam
import random
from typing import List, Dict, Optional, Tuple
import logging
from dataclasses import dataclass


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


class ReadSelector:
    """
    Select appropriate reads for insertion based on filtering criteria.
    """
    
    def __init__(self, bam_path: str, min_mapq: int = 20, max_soft_clip_ratio: float = 0.1,
                 min_read_length: int = 10000, max_read_length: int = 20000):
        """
        Initialize ReadSelector.
        
        Args:
            bam_path: Path to BAM file
            min_mapq: Minimum mapping quality
            max_soft_clip_ratio: Maximum soft clipping ratio
            min_read_length: Minimum read length
            max_read_length: Maximum read length
        """
        self.bam_path = bam_path
        self.min_mapq = min_mapq
        self.max_soft_clip_ratio = max_soft_clip_ratio
        self.min_read_length = min_read_length
        self.max_read_length = max_read_length
        self.logger = logging.getLogger(__name__)
        
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
        Check if read passes all filtering criteria (primary alignment, mapq, read length, soft clipping)
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
    
    def select_reads(self, n_reads: int) -> List[Tuple[pysam.AlignedSegment, ReadMetadata]]:
        """
        Select reads from BAM file using random sampling across all chromosomes.
        
        Args:
            n_reads: Number of reads to select
            
        Returns:
            List of tuples containing (read, metadata)
        """
        selected_reads = []
        
        with pysam.AlignmentFile(self.bam_path, "rb") as bam:
            # Get random sampling regions across all chromosomes
            sampling_regions = self._get_chromosome_sampling_regions(bam)
            
            # Keep track of regions we've tried to avoid infinite loops
            regions_tried = 0
            max_regions_to_try = len(sampling_regions) * 3  # Allow multiple passes
            
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
                    for read in bam.fetch(chrom, start, end):
                        if len(selected_reads) >= n_reads:
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
                        
                except Exception as e:
                    # Skip regions that cause errors (e.g., invalid coordinates)
                    self.logger.warning(f"Skipping region {chrom}:{start}-{end} due to error: {e}")
                    continue
        
        self.logger.info(f"Selected {len(selected_reads)} reads from {self.bam_path} "
                        f"after sampling {regions_tried} regions")
        return selected_reads