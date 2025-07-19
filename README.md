Introduction
============

Benchmark the speed and memory footprint of common genomic-interval operations across several libraries (Python's bioframe and pyranges; R's GenomicRanges; and bash's bedtools).

The workflow:

1. Grabs real data describing the size of sequences (human genome, human protein)
2. Creates synthetic BED files of varying sizes so that run-time can be swept over
many problem sizes.
3. Runs the same operation with every library and measures wall-clock time, and peak resident memory (RSS)
4. Collects all benchmarks into a single CSV.

The JSON and tabular files in this repo define the range of dataset parameters to be tested, as well as which operations are tested.

The code under the folder scripts/ defines how each operation is implemented for each library.


Usage
=====
- Install Snakemake and all sequence-interval libraries, then run  ```snakemake ```. Check the log for errors, and find tabular output in the ```results``` folder.


Notes on the different behavior of methods/operations in different libraries
============================================================================

Find nearest interval:
- Overlaps: GenomicRanges does not allow you to ignore the overlapping intervals. The other libraries ignore overlaps.
- Ties: Bioframe does not get all ties (nearest at same distance), and instead it only keeps one. 
- k-nearest: In the benchmark, Bioframe, Pyranges, and BEDTools are requested to get the two nearest intervals per interval, but GenomicRanges cannot; thus, GenomicRanges only finds the single nearest interval.
