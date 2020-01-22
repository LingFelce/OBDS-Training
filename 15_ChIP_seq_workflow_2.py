"""
Created on Thu May 30 13:08:52 2019
@author: rrigby
From other OBDS training cohort

"""
from ruffus import *
from cgatcore import pipeline as P
import sys

PARAMS = P.get_parameters("chipseq_pipeline.yml")

@follows(mkdir('fastqc'))
@transform('*.fastq.gz', regex(r'(.*).fastq.gz'), r'fastqc/\1_fastqc.zip')
def fastqc(input_file, output_file):
    statement = 'fastqc -t 8 %(input_file)s -o fastqc'
    P.run(statement, job_queue=PARAMS['q'], job_threads=8)

@transform('*_1.fastq.gz', suffix('_1.fastq.gz'), '.sam') 
def bowtie2(input_file, output_file):
    read_2 = input_file.replace('1.fastq.gz', '2.fastq.gz')
    statement = '''bowtie2 -x %(bowtie2_ref)s 
    -1 %(input_file)s -2 %(read_2)s -p %(bowtie2_threads)s -S %(output_file)s 
    2> %(output_file)s.log'''
    P.run(statement, job_queue=PARAMS['q'], job_threads=PARAMS['bowtie2_threads'], job_memory = '8G')
    
@transform(bowtie2, suffix('.sam'), '.bam')     
def sam_to_bam(input_file, output_file):
    statement = '''samtools sort -@ 12 -o %(output_file)s  %(input_file)s 
    2> %(output_file)s.log'''
    P.run(statement, job_queue=PARAMS['q'], job_threads=12, job_memory = '8G')
    
@transform(sam_to_bam, suffix('.bam'), '.bam.bai')     
def bam_index(input_file, output_file):
    statement = 'samtools index %(input_file)s -@ 12 2> %(output_file)s.log'
    P.run(statement, job_queue=PARAMS['q'], job_threads=12, job_memory = '8G')  
 
@follows(bam_index)    
@transform(sam_to_bam, suffix('.bam'), '.idxstat')     
def samtools_idxstat(input_file, output_file):
    statement = '''samtools idxstats %(input_file)s -@ 12 > %(output_file)s
    2> %(output_file)s.log'''
    P.run(statement, job_queue=PARAMS['q'], job_threads=12, job_memory = '8G') 

@follows(bam_index)    
@transform(sam_to_bam, suffix('.bam'), '_dna_metrics.txt')     
def alignment_summary_metrics(input_file, output_file):
    statement = '''picard CollectAlignmentSummaryMetrics R=%(picard_ref)s I=%(input_file)s O=%(output_file)s
    2> %(output_file)s.log'''
    P.run(statement, job_queue=PARAMS['q'], job_memory = '32G') 

@follows(bam_index)
@transform(sam_to_bam, suffix('.bam'), '.flagstat')     
def samtools_flagstat(input_file, output_file):
    statement = '''samtools flagstat %(input_file)s -@ 12 > %(output_file)s
    2> %(output_file)s.log'''
    P.run(statement, job_queue=PARAMS['q'], job_threads=12, job_memory = '8G')   
    
@follows(samtools_flagstat, samtools_idxstat, alignment_summary_metrics, fastqc)
@merge(samtools_flagstat, 'mapping_qc.html')    
def multiqc(input_file, output_file):
    statement = '''multiqc -o %(output_file)s .'''
    P.run(statement, job_queue=PARAMS['q'], job_threads=12, job_memory = '8G') 
    
@follows(bam_index)    
@transform(sam_to_bam, suffix('.bam'), '.rmdup.bam')     
def remove_duplicates(input_file, output_file):
    statement = '''picard MarkDuplicates I=%(input_file)s O=%(output_file)s M=%(output_file)s.markdups.metrics.txt REMOVE_DUPLICATES=TRUE
    2> %(output_file)s.log'''
    P.run(statement, job_queue=PARAMS['q'], job_memory = '32G')
       
@transform(remove_duplicates, suffix('.bam'), '.bam.bai')     
def index_rmdup(input_file, output_file):
    statement = '''samtools index %(input_file)s'''
    P.run(statement, job_queue=PARAMS['q'], job_memory = '8G')

@collate(remove_duplicates, regex(r'(.*)_.*.rmdup.bam'), r'\1_peaks.narrowPeak')
def macs2(input_files, output_file):
    output_name = output_file.replace("_peaks.narrowPeak", "")  
    chip_file, input_file = input_files
    statement = '''macs2 callpeak -t %(chip_file)s -c %(input_file)s -f BAMPE -g mm -n 
    %(output_name)s -B -q 0.01'''
    P.run(statement, job_condaenv = "macs2-env", job_queue=PARAMS['q'], job_threads=12, job_memory = '8G')

@merge(macs2, 'merged_peaks.bed')
def merge_bed(input_files, output_file):
    input_file_list = ' '.join(input_files)
    statement = '''cat %(input_file_list)s | sort -k1,1 -k2,2n | bedtools merge > %(output_file)s'''   
    P.run(statement, job_condaenv = "obds_env", job_queue=PARAMS['q'], job_memory = '8G')

@follows(index_rmdup) 
@merge(remove_duplicates, 'count_peaks_coverage.tsv')
def count_peaks(input_files, output_file):
    input_file_list = ' '.join(input_files)
    output_header= "chr\tstart\tstop\t" + "\t".join(input_files) +"\n"
    output_header=output_header.replace(".rmdup.bam", "")
    statement = '''bedtools multicov -bams %(input_file_list)s -bed merged_peaks.bed > %(output_file)s 
    && sed -i '1i %(output_header)s' %(output_file)s '''   
    P.run(statement, job_condaenv = "obds_env", job_queue=PARAMS['q'], job_memory = '8G')



   
if __name__=="__main__":
    sys.exit(P.main(sys.argv))
