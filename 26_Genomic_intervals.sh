18th October 2019

# Interval file types:
  Bed files = 0-based e.g peak calling output (regions of H3K27ac, promoter regions -5kb from TSS)
  GTF files = 1-based e.g annotations (protein coding genes, exons, UTRs)
  VCF files = 1-based e.g variant files (SNPs, INDELs, translocations)

Genome browser track
Bigwig or bam file

# Uses of intervals - use BEDTOOLS
  find which promoters overlap H3K27ac peaks and other marks of interest
  find promoters which don't overlap
  count number of different histone peaks in promoter
  identify genes overlapping SNPs/variants
  find nearest gene to peak/SNP
  
# Bedtools (https://bedtools.readthedocs.io/en/latest/):
  A = query file - processed line by line (biggest file)
  B = target/subject file - loaded into memory (smaller file)
  > 30 different functions - closest gene to regions, find overlaps, find regions specific to only one file
  Gzip all files
  Need to redirect file to output using > 
  -s flag tells it to take strandedness into account
  
# Example of usage of Bedtools:
Let’s imagine you have a BED file of ChiP-seq peaks from two different experiments. You want to identify peaks that were observed in 
both experiments (requiring 50% reciprocal overlap) and for those peaks, you want to find to find the closest, non-overlapping gene. 
Such an analysis could be conducted with two, relatively simple bedtools commands.

## intersect the peaks from both experiments.
## -f 0.50 combined with -r requires 50% reciprocal overlap between the peaks from each experiment.
$ bedtools intersect -a exp1.bed -b exp2.bed -f 0.50 -r > both.bed

## find the closest, non-overlapping gene for each interval where both experiments had a peak
## -io ignores overlapping intervals and returns only the closest, non-overlapping interval (in this case, genes)
$ bedtools closest -a both.bed -b genes.bed -io > both.nearest.genes.txt
  
Bedtools intersect - find/exclude overlaps # (http://quinlanlab.org/tutorials/bedtools/bedtools.html#bedtools-intersect)
Only overlapping intervals
OR annotate with coverage/overlaps from other BED/VCF/GTF OR complement - return all intervals in genome not covered

# Statistical tests - is proportion of overlapping statistically significant? 
# (https://bedtools.readthedocs.io/en/latest/content/tools/fisher.html)
Perform Fisher’s Exact Test on the number of overlaps/unique intervals between 2 files. Traditionally, in order to test whether 2 sets of 
intervals are related spatially, we resort to shuffling the genome and checking the simulated (shuffled) versus the observed. We can do 
the same analytically for many scenarios using Fisher’s Exact Test. This implementation can calculate the number of overlaps and the 
number of intervals unique to each file and it infers (or accepts) the number that are not present in each file.

Merge intervals that are too close - takes 1 file containing all intervals and merge all overlapping intervals together

Closest - reports nearest interval even if not overlapping eg finding closest gene to SNP or enhancer

Genome coverage - counts across genome how many intervals overlap across genome
Multicoverage - counts alignments from multiple position sorted and indexed BAM files that overlap BED file
Maskfasta, shift, flank, slop
