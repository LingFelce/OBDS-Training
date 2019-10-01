#13_RNA-Seq Workflow 1st October 2019

"""
rna_seq_pipeline.yaml
queue: all.q
threads: 12
hisat:
	options: --fr #forward strand
	ref: /ifs/mirror/genomes/hisat2/mm10 #genome reference
picard:
	ref: /ifs/mirror/genomes/plain/mm10.fasta

featurecounts:
	gtf: /ifs/obds-training/exercises/rnaseq/genes.gtf.gz
	options: -s 1 #strand option (see featureCounts --help)


Regex expression regex101
(.*)_[1-2].fastq.gz
gm
1st Capturing Group (.*)
.* matches any character (except for line terminators)
* Quantifier — Matches between zero and unlimited times, as many times as possible, giving back as needed (greedy)
_ matches the character _ literally (case sensitive)
Match a single character present in the list below [1-2]
1-2 a single character in the range between 1 (index 49) and 2 (index 50) (case sensitive)
. matches any character (except for line terminators)
fastq matches the characters fastq literally (case sensitive)
. matches any character (except for line terminators)
gz matches the characters gz literally (case sensitive)
Global pattern flags
g modifier: global. All matches (don't return after first match)
m modifier: multi line. Causes ^ and $ to match the begin/end of each line (not only begin/end of string)

"""

import gzip
from ruffus import *
from cgatcore import pipeline as P
import sys
import statistics

P.get_parameters('rnaseq_pipeline.yaml') #.yml for cbrg, .yaml for cgat

@follows(mkdir('fastqc')) #make fastqc folder before running code below
@transform('*.fastq.gz', regex(r'(.*).fastq.gz'),r'fastqc/\1_fastqc.html') #find all fastq.gz files, save to fastqc folder and use name_fastqc.html
def fastqc (infile, outfile): #next time call functions run_fastqc etc so don't have confused!
    statement = '''fastqc -q -t %(threads)s --nogroup %(infile)s --outdir fastqc''' #need to direct output
    P.run(statement, job_queue=P.PARAMS['queue'], job_threads=P.PARAMS['threads']) #only for setting run parameters eg memory, threads, queue

@follows(mkdir('sam')) #make sam folder before running code below
@collate('*.fastq.gz', regex(r'(.*)_[1-2].fastq.gz'), r'sam/\1.sam') #look for fastq.gz files with same name ending in _1 or _2, output as sam file
def hisat2 (infiles, outfile):
    read1, read2 = infiles #2 infiles
    statement = '''hisat2 -p %(threads)s %(hisat_option)s -x %(hisat_ref)s #hisat options in yaml file
    -1 %(read1)s -2 %(read2)s -S %(outfile)s'''
    P.run(statement, job_queue=P.PARAMS['queue'], job_threads=P.PARAMS['threads'], job_memory ='8G')

@follows(mkdir('bam')) #make folder bam before running code below
@transform(hisat2, regex(r'sam/(.*).sam'),r'bam/\1.bam') #dependent on hisat2 running properly.Take all .sam files, convert to .bam 
def samtools_sort (infile, outfile):
    statement = '''samtools sort %(infile)s -o %(outfile)s --output-fmt BAM -@ %(threads)s''' #probably don't need --output bit as already specified above
    P.run(statement, job_queue=P.PARAMS['queue'], job_threads=P.PARAMS['threads'], job_memory='8G')

@transform(samtools_sort, regex(r'bam/(.*).bam'), r'bam/\1.bam.bai') #dependent on samtools_sort above. Take all bam and output as .bam.bai
def samtools_index  (infile, outfile):
    statement = '''samtools index %(infile)s %(outfile)s'''
    P.run(statement, job_queue=P.PARAMS['queue'], job_memory='8G')

@follows(mkdir('picard'), samtools_index) #create picard folder, run after samtools_index
@transform(samtools_sort, regex(r'bam/(.*).bam'), r'picard/\1.txt') #use output from samtools_sort (.bam files), save in picard folder as .txt
def picard (infile, outfile):
	statement = ' '.join(['picard CollectAlignmentSummaryMetrics', 'R=%(picard_ref)s', 'I=%(infile)s', 'O=%(outfile)s',])
	P.run(statement, job_queue=P.PARAMS['queue'], job_memory='16G')

@merge(picard, 'multiqc_report.html') #run after picard, output is .html file
def multiqc (infile, outfile):
	statement = '''multiqc .'''
	P.run(statement)

@follows(mkdir('featureCounts'),samtools_index)
@merge(samtools_sort, 'featureCounts/total.counts')
def featureCounts (infiles, outfile):
	infile_string = ' '.join(infiles) #converts list bam files into string so can be read properly by featureCounts
	statement = '''featureCounts %(featurecounts_options)s -T %(threads)s -a %(featurecounts_gtf)s -o %(outfile)s %(infile_string)s'''
	P.run(statement, job_queue=P.PARAMS['queue'], job_threads=P.PARAMS['threads'], job_memory='4G'
	

if __name__ == "__main__":
    sys.exit( P.main(sys.argv) )
