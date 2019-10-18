18th October 2019

# Interval file types:
  # Bed files = 0-based e.g peak calling output (regions of H3K27ac, promoter regions -5kb from TSS)
  # GTF files = 1-based e.g annotations (protein coding genes, exons, UTRs)
  # VCF files = 1-based e.g variant files (SNPs, INDELs, translocations)

# Genome browser track
# Bigwig or bam file

# Uses of intervals - use BEDTOOLS
  # find which promoters overlap H3K27ac peaks and other marks of interest
  # find promoters which don't overlap
  # count number of different histone peaks in promoter
  # identify genes overlapping SNPs/variants
  # find nearest gene to peak/SNP
  
# Bedtools (https://bedtools.readthedocs.io/en/latest/):
  # A = query file - processed line by line (biggest file)
  # B = target/subject file - loaded into memory (smaller file)
  # > 30 different functions - closest gene to regions, find overlaps, find regions specific to only one file
  # Gzip all files
  # Need to redirect file to output using > 
  # -s flag tells it to take strandedness into account
  
# Example of usage of Bedtools:
# Let’s imagine you have a BED file of ChiP-seq peaks from two different experiments. You want to identify peaks that were observed in 
# both experiments (requiring 50% reciprocal overlap) and for those peaks, you want to find to find the closest, non-overlapping gene. 
# Such an analysis could be conducted with two, relatively simple bedtools commands.

## intersect the peaks from both experiments.
## -f 0.50 combined with -r requires 50% reciprocal overlap between the peaks from each experiment.
bedtools intersect -a exp1.bed -b exp2.bed -f 0.50 -r > both.bed

## find the closest, non-overlapping gene for each interval where both experiments had a peak
## -io ignores overlapping intervals and returns only the closest, non-overlapping interval (in this case, genes)
bedtools closest -a both.bed -b genes.bed -io > both.nearest.genes.txt
  
# Bedtools intersect - find/exclude overlaps # (http://quinlanlab.org/tutorials/bedtools/bedtools.html#bedtools-intersect)
# Only overlapping intervals
# OR annotate with coverage/overlaps from other BED/VCF/GTF OR complement - return all intervals in genome not covered

# Statistical tests - is proportion of overlapping statistically significant? 
# (https://bedtools.readthedocs.io/en/latest/content/tools/fisher.html)
# Perform Fisher’s Exact Test on the number of overlaps/unique intervals between 2 files. Traditionally, in order to test whether 2 sets of 
# intervals are related spatially, we resort to shuffling the genome and checking the simulated (shuffled) versus the observed. We can do 
# the same analytically for many scenarios using Fisher’s Exact Test. This implementation can calculate the number of overlaps and the 
# number of intervals unique to each file and it infers (or accepts) the number that are not present in each file.

# Merge intervals that are too close - takes 1 file containing all intervals and merge all overlapping intervals together

# Closest - reports nearest interval even if not overlapping eg finding closest gene to SNP or enhancer

# Genome coverage - counts across genome how many intervals overlap across genome
# Multicoverage - counts alignments from multiple position sorted and indexed BAM files that overlap BED file

# Can also use Pybedtools - use bedtools in Python

# Exercise: Bedtools tutorial (http://quinlanlab.org/tutorials/bedtools/bedtools.html)
# Look at tutorial for full help, here be notes.

# looking at IGV browser in terminal will automatically create indexes for files that it needs
module load bio/IGV
igv.sh

# find intersects and show original records of where CpG islands and exons are, and number of base pairs that overlap
# -wo only overlaps with A are reported (number of bp)
# default is to report overlaps between features in A and B so long as at least one basepair of overlap exists
bedtools intersect -a cpg.bed -b exons.bed -wo | head -5
chr1	28735	  29810	  CpG:_116chr1	29320	  29370	  NR_024540_exon_10_0_chr1_29321_r	0	-	50
chr1	135124	135563	CpG:_30	chr1	134772	139696	NR_039983_exon_0_0_chr1_134773_r	0	-	439
chr1	327790	328229	CpG:_29	chr1	324438	328581	NR_028322_exon_2_0_chr1_324439_f	0	+	439
chr1	327790	328229	CpG:_29	chr1	324438	328581	NR_028325_exon_2_0_chr1_324439_f	0	+	439
chr1	327790	328229	CpG:_29	chr1	327035	328581	NR_028327_exon_3_0_chr1_327036_f	0	+	439

bedtools intersect -sorted # use if input files are sorted, will be quicker
bedtools intersect -names # if intersecting multiple B files, name them so can tell which is which in output
bedtools merge # for example alignments (narrowPeaks) from ChIP-Seq - would need to be (Unix) sorted first
bedtools merge -i exons.bed | head -n 20 # will not save over original file, only output to terminal
bedtools complement # For ChIP-seq peaks, may want to know which regions of the genome are not bound by the factor assayed

bedtools genomecov
# By default, bedtools genomecov will compute a histogram of coverage for the genome file provided. The default output format is as follows:
# 1. chromosome (or entire genome)
# 2. depth of coverage from features in input file (eg number of exons that overlap with genomic region)
# 3. number of bases on chromosome (or genome) with depth equal to column 2.
# 4. size of chromosome (or entire genome) in base pairs
# 5. fraction of bases on chromosome (or entire genome) with depth equal to column 2.
bedtools genomecov -i exons.bed -g genome.txt
chr1	0	241996316	249250621	0.970896 # 97% of chromosome 1 have no overlapping exons
chr1	1	4276763	249250621	0.0171585 # 1.7% of chromosome 1 has 1 overlapping exon
chr1	2	1475526	249250621	0.00591985 # 0.59% of chromosome 1 has 2 overlapping exons
chr1	3	710135	249250621	0.00284908 # ... and so on ...
