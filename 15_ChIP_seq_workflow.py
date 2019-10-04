"""
15_ChIP-seq workflow - 3rd October 2019

ChIP-seq analysis: identify TF binding sites (promoters/enhancers, target genes, sequence motifs, interacting partners, 
differential binding), chromatin state (active, poised, repressed)

Sources of bias:
- starting material - need many cells for robust results, consistent between replicates
- fragmentation of library (sonication - open chromatin shears more easily; enzymatic - sequence bias)
- antibody specificity. Use cocktail of antibodies - target different epitopes, avoid epitope masking. Use polyclonal Ab
- batch effects. 
ChIP-seq controls:
- mock (just Ab) - very little material - worth sequencing?
- no Ab (chromatin but IP with IgG) - very little material - worth sequencing?
- input control (chromatin not IP) - fragmentation bias
- spike in external control eg drosophila DNA
Sequencing considerations:
- single/paired end - paired end better for removal of duplicates
- read length long enough for unique mapping
- sequencing depth - greater for input (2 x IP), greater for broad peaks, Encode guidelines.
Mapping QC:
- > 90% mapping
- read filtering - mitochondrial reads, duplicates, multi-mapping reads, reads not properly paired
Encode metrics - PCR bottlenecking coefficients, non-redundant fraction
Peak calling:
- pool sample - better coverage
- caveats - samples should be of similar depth, otherwise one sample might influence peak calling more than others (down-sample?)
Peak calling tools - MACS2 (most widely used), SICER, LanceOtron, Homer, SPP, PeakRanger
Viewing peaks (in IGV) - different file formats:
- Bam (check raw reads, memory intensive) 
- Bedcoverage (breaks genome into regions and gives coverage score, make with Bedtools, plain text, less memory intensive)
- Bigwig (like bedcoverage, but compressed and indexed; faster, less memory intensive)
Peak calling QC:
- reproducibility (overlap replicates Bedtools, keep peaks 2+ samples)
- irreproducible discovery rate (Encode project), rank peaks on p-value, identify reproducible peaks on correlation of ranks
- black-listed regions (sticky regions that always appear in ChIP-seq)
Peak annotation:
- genomic context (promoter, enhancer etc) - GREAT online tool
- target genes (nearest neighbour) - GREAT
Differential binding (2 different conditions)
- featureCounts (count reads under peaks in each condition)
- DESeq2/edgeR (treat like RNA-seq data, must have >= 3 replicates)
Motif analysis: MEME suit, Homer motif analysis
"""

#download files in command line -q quiet & send to background
#wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/001/SRR1592591/SRR1592591_1.fastq.gz &
#wget -q ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/001/SRR1592591/SRR1592591_2.fastq.gz &
#remember to check md5sum!

#create soft links in same folder - used Filezilla to transfer over and edit as shell script in nano
#have to use # at end as Windows text editor adds $ r at end of line
#change name to something meaningful and easy to use regex later!
#ln -s /ifs/obds-training/exercises/chipseq_pipeline/SRR1592591_1.fastq.gz mm_rep2_dex_IP_1.fastq.gz #
#ln -s /ifs/obds-training/exercises/chipseq_pipeline/SRR1592591_2.fastq.gz mm_rep2_dex_IP _2.fastq.gz #

#run shell script to create links with new meaningful names
#sh <shell script>

#copy links to new own directory
#cp -d <link name> <directory>

"""
chipseq_pipeline.yml

queue: all.q
threads: 12
memory: 8G
bowtie2:
    options:
    ref: /ifs/mirror/genomes/bowtie/mm10
picard:
    ref: /ifs/mirror/genomes/plain/mm10.fasta
"""
import sys
import gzip
from cgatcore import pipeline as P
from ruffus import *

P.get_parameters('chipseq_pipeline.yml')

@follows(mkdir('fastqc')) #same as before, can run independently of other processes
@transform('*.fastq.gz', regex(r'(.*).fastq.gz'),r'fastqc/\1_fastqc.html')
def qc_reads(infile, outfile):
    statement = 'fastqc -q -t %(threads)s --nogroup %(infile)s --outdir fastqc'
    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_memory = P.PARAMS['memory'])

@follows(mkdir('sam'))
@collate('*.fastq.gz', regex(r'(.*)_[1-2].fastq.gz'), r'sam/\1.sam')
def align_reads(infiles, outfile):
    read1, read2 = infiles
    options = ''
    if P.PARAMS['bowtie2_options']: #if there are bowtie2 options specified in .yml
        options = P.PARAMS['bowtie2_options']
    statement = '''bowtie2 -p %(threads)s %(options)s -x %(bowtie2_ref)s #have to use bowtie2 instead of hisat2
    -1 %(read1)s -2 %(read2)s -S %(outfile)s''' 
    P.run(statement,
          job_queue   = P.PARAMS['queue'],
          job_threads = P.PARAMS['threads'],
          job_memory  = P.PARAMS['memory'])

@follows(mkdir('bam'))
@transform(align_reads, regex(r'sam/(.*).sam'), r'bam/\1.bam') #can skip this step and go straight from .fastq.gz to bam - see later
def sort_reads(infile, outfile):
    """Sort aligned reads."""
    statement = 'samtools sort %(infile)s -o %(outfile)s --output-fmt BAM -@ %(threads)s'
    P.run(statement,
          job_queue   = P.PARAMS['queue'],
          job_threads = P.PARAMS['threads'],
          job_memory  = P.PARAMS['memory'])

@transform(sort_reads, regex(r'bam/(.*).bam'), r'bam/\1.bam.bai')
def create_bam_file_index(infile, outfile):
    """Bam files are compressed. The index allows fast access to different
    slices of the file."""
    statement = 'samtools index %(infile)s %(outfile)s'
    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_memory = P.PARAMS['memory'])

@follows(mkdir('picard'), create_bam_file_index)
@transform(sort_reads, regex(r'bam/(.*).bam'), r'picard/\1.txt')
def get_alignment_statistics(infile, outfile):
    """Quality control alignments and write a report."""
    statement = 'picard CollectAlignmentSummaryMetrics R=%(picard_ref)s I=%(infile)s O=%(outfile)s'
    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_memory = P.PARAMS['memory'])
    
@follows(create_bam_file_index, mkdir('deduplicated')) #remove duplicates    
@transform(align_reads, regex(r'bam/(.*.bam)'), r'deduplicated/\1')
def remove_duplicates(infile, outfile):
    stats = outfile.replace('.bam', '_marked_duplicates.txt')
    statement = ' '.join(['picard MarkDuplicates',
                          'I=%(infile)s',
                          'O=%(outfile)s',
                          'M=%(stats)s',
                          'REMOVE_DUPLICATES=true',
                          'REMOVE_SEQUENCING_DUPLICATES=true',
                          'CREATE_INDEX=true'])
    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_memory = P.PARAMS['memory'])

@follows(mkdir('peaks')) #call peaks - generates .narrowPeak file which is like a .bed file
@collate(remove_duplicates, regex(r'deduplicated/(.*)_(gr|ip).bam'), r'peaks/\1_peaks.narrowPeak')
def call_peaks (infiles, outfile):
    treatment, control = infiles
    file_base = outfile.replace('_peaks.narrowPeak','')
    statement = '''macs2 callpeak -t %(treatment)s -c %(control)s -n %(file_base)s'''
    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_memory = P.PARAMS['memory'],
          job_condaenv = 'macs2-env') #separate conda environment

"""
Next steps:
Use bedtools merge to combine narrowPeak files
Use awk to change narrowPeak files into .saf format
Use featureCounts to count reads under peaks - .saf file will take place of .gft file in command, will give information on peaks,
and use .bam file as normal to give information on counts.
"""
    
#use this to combine aligning reads and sorting, so go straight from .fastq.gz to .bam
@follows(mkdir('bam'))
@collate('*.fastq.gz', regex(r'(.*)_[1-2].fastq.gz'), r'bam/\1.bam')
def align_reads(infiles, outfile):
    
    ''' Aligns fq files using bowtie2 before conversion to bam file using
        Samtools view. Bam file is then sorted and the unsorted bam file is replaced'''
    
    fq1, fq2 = infiles   
    sorted_bam = outfile.replace('.bam', '_sorted.bam')
    options = ''
        
    if P.PARAMS['bowtie2_options']:
        options = P.PARAMS['bowtie2_options']
    
    cmd = '''bowtie2 -x %(bowtie2_index)s -1 %(fq1)s -2 %(fq2)s -p %(threads)s %(options)s |
             samtools view -b - > %(outfile)s &&
             samtools sort -@ %(threads)s -m 5G -o %(sorted_bam)s %(outfile)s &&
             mv %(sorted_bam)s %(outfile)s'''

    P.run(cmd, 
          job_queue=P.PARAMS['queue'], 
          job_threads=P.PARAMS['threads'])
          
     
if __name__ == "__main__":
    sys.exit( P.main(sys.argv) )
