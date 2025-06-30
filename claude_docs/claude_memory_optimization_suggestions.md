Memory Optimization Plan for Nova

    1. Lazy Read Processing (High Impact)

    - Modify ReadSelector to return minimal read metadata instead of full
    pysam.AlignedSegment objects
    - Store only: read_name, chr, pos, query_length, mapq in selection phase
    - Re-fetch sequences on-demand during insertion using pysam.fetch() by
    coordinates

    2. Streaming Insertion Pipeline (High Impact)

    - Convert ReadInserter.insert_random_mode() to use generator/iterator
    pattern
    - Process reads one-at-a-time instead of loading all into memory
    - Yield (InsertionRecord, SeqRecord) pairs as they're created
    - Write results to disk immediately instead of accumulating in lists

    3. Registry Memory Management (Medium Impact)

    - Add VariantRegistry.get_sequence_by_id() that can optionally clear after
     retrieval
    - Implement registry chunking for very large simulations
    - Consider disk-backed storage for large sequence sets

    4. Batch Processing Architecture (Medium Impact)

    - Add configurable batch size parameter (e.g., 50-100 reads per batch)
    - Process simulation in batches with memory cleanup between batches
    - Append results to output files progressively

    5. CLI Integration (Low Impact)

    - Update CLI to use new streaming architecture
    - Add --batch-size parameter for memory tuning
    - Add progress reporting for large simulations

    This should reduce memory usage from ~75MB+ to ~5-10MB for 500-read
    simulations while maintaining identical functionality.
