14_Alignment free RNA-Seq 2nd October 2019

"""
Kallisto pseudo alignment, skips a lot of steps - mapping, alignment, counts etc.
For each read, finds compatible transcripts by building De-Bruijn graph from reference transcriptome.
Bootstrapping - sample reads randomly with replacement, estimate transcript abundances, repeat 100x
- get confidence interval for abundance of each transcript, identify transcripts with good quality quantification estimates
Use Ensembl geneset (manually curated); Refseq better for small RNAs, longer first exons (TSS)
Output from Kallisto = TPM (transcripts per kilobase million) counts table, use R txiimport package to convert to counts for DESeq2/edgeR
Can use Sleuth for DE directly from pseudoalignments
Transcript level differential expression is hard, recommend quantifying at gene level

kallisto_pipeline_cgat.yml
queue: all.q
threads:12
index:
	fasta: /ifs/mirror/genomes/plain/mm10.fasta
	gtf: /ifs/obds-training/exercises/rnaseq/genes.gtf.gz
quant:
	option: --fr-stranded #first read forward

"""
#To avoid conflicts with github, in terminal create new branch "ling"
git branch ling
git checkout ling
#Then make edits to files on branch locally
git add 
git commit -m "comment"
git push origin ling
#On github website, select branch and go to new pull request. Merge branch with master. 
#Back in terminal, switch back to master, pull from github and delete branch.
git branch master
git pull
git branch -d ling

from ruffus import *
from cgatcore import pipeline as P
import sys

P.get_parameters('kallisto_pipeline_cgat.yml') 

@follows(mkdir('fastqc')) #same as yesterday, fastq quality control
@transform('*.fastq.gz', regex(r'(.*).fastq.gz'),r'fastqc/\1_fastqc.html')
def fastqc (infile, outfile):
    statement = '''fastqc -q -t %(threads)s --nogroup %(infile)s --outdir fastqc'''
    P.run(statement, job_queue=P.PARAMS['queue'], job_threads=P.PARAMS['threads'])

@merge([P.PARAMS['index_fasta'],P.PARAMS['index_gtf']], 'transcripts.fasta') #input files from .yml put as a list, output name
def create_fasta (infiles,outfile):
    fasta,gtf = infiles #unpack list
    statement = '''zcat %(gtf)s | gffread -w %(outfile)s -g %(fasta)s - ''' #unzip gtf.gz first, then use - as standard input after -g
    P.run(statement, job_queue=P.PARAMS['queue'], job_memory ='8G')
   
@transform(create_fasta, suffix('.fasta'),'.index') #follow from create_fasta, use output as input
def create_index (infile, outfile):
    statement = '''kallisto index -i %(outfile)s %(infile)s'''
    P.run(statement, job_queue=P.PARAMS['queue'], job_memory ='8G')

@follows(mkdir('quant'), create_index) #create quant folder, follow from create_index to use output as input
@collate('*.fastq.gz', regex(r'(.*)_[1-2].fastq.gz'), r'quant/\1/abundances.h5') #find all fastq.gz files, output into quant folder/name of sample/abundances.h5 (one of 3 output files generated by default)
def run_quant(infiles, outfile):
    read1, read2 = infiles #set variables
    outdir = outfile.replace('abundances.h5','') #set directory for output, remove file name
    statement = ''' kallisto quant -t %(threads)s %(quant_option)s -i transcripts.index -o %(outdir)s %(read1)s %(read2)s''' #-o is only for output directory
    P.run(statement, job_queue=P.PARAMS['queue'], job_memory ='8G', job_threads=P.PARAMS['threads'])

@merge(run_quant, 'multiqc_report.html')
def multiqc (infile, outfile):
    statement = '''multiqc .'''
    P.run(statement)