"""
16_Capture-C 4th October 2019

Chromosome conformation capture techniques (often abbreviated to 3C technologies or 3C-based methods) are a set of molecular biology 
methods used to analyze the spatial organization of chromatin in a cell. 
These methods quantify the number of interactions between genomic loci that are nearby in 3-D space, but may be separated by many 
nucleotides in the linear genome. Such interactions may result from biological functions, such as promoter-enhancer interactions, 
or from random polymer looping, where undirected physical motion of chromatin causes loci to collide.
Interaction frequencies may be analyzed directly, or they may be converted to distances and used to reconstruct 3-D structures.

First, the cell genomes are cross-linked with formaldehyde, which introduces bonds that "freeze" interactions between genomic loci. 
Treatment of cells with 1-3% formaldehyde, for 10-30min at room temperature is most common, however, standardization for preventing high 
protein-DNA cross linking is necessary, as this may negatively affect the efficiency of restriction digestion in the subsequent step. 
The genome is then cut into fragments with a restriction endonuclease. The size of restriction fragments determines the resolution of 
interaction mapping. Restriction enzymes (REs) that make cuts on 6bp recognition sequences, such as EcoR1 or HindIII, are used for this 
purpose, as they cut the genome once every 4000bp, giving ~ 1 million fragments in the human genome. For more precise interaction 
mapping, a 4bp recognizing RE may also be used. The next step is, random ligation. This takes place at low DNA concentrations in the 
presence of T4 DNA ligase, such that ligation between cross-linked interacting fragments is favored over ligation between fragments 
that are not cross-linked. 
Subsequently, interacting loci are quantified by amplifying ligated junctions by PCR methods

Sequence capture-based methods
A number of methods use oligonucleotide capture to enrich 3C and Hi-C libraries for specific loci of interest.
These methods include Capture-C, NG Capture-C, Capture-3C, and Capture Hi-C.
These methods are able to produce higher resolution and sensitivity than 4C based methods.

Capture-C similar to 4C:
Chromosome conformation capture-on-chip (4C) captures interactions between one locus and all other genomic loci. 
It involves a second ligation step, to create self-circularized DNA fragments, which are used to perform inverse PCR. 
Inverse PCR allows the known sequence to be used to amplify the unknown sequence ligated to it.
In contrast to 3C and 5C, the 4C technique does not require the prior knowledge of both interacting chromosomal regions. 
Results obtained using 4C are highly reproducible with most of the interactions that are detected between regions proximal to one another. 
On a single microarray, approximately a million interactions can be analyzed.

Use biotinylated probes to enrich for captured sequences, don't sequence whole library, just sites of interest.

See Jelena's notes:
http://userweb.molbiol.ox.ac.uk/public/telenius/CCseqBasicManual/ppMan/CCseqBasic/2_workflow/index.html
http://userweb.molbiol.ox.ac.uk/public/telenius/CCseqBasicManual/ppMan/CCseqBasic/DOCS/readsFragments_multipage.pdf
http://userweb.molbiol.ox.ac.uk/public/telenius/captureManual/UserManualforCaptureCanalysis.pdf

capturec_pipeline.yml
queue: all.q
threads: 12
memory: 8G
bowtie2:
    options: --reorder
    ref: /ifs/mirror/genomes/bowtie/mm10
picard:
    ref: /ifs/mirror/genomes/plain/mm10.fasta
ccanalyser:
    ref: /ifs/obds-training/exercises/capturec/mm10.txt
    genome: mm10
    oligo: /ifs/obds-training/exercises/capturec/fragmentmm10.txt
    pu: "http://www.cgat.org/downloads/"
    pf: /ifs/obds-training/lingf/week3/capture_c
    memory: 128G

"""
#Capture-C pipeline exercise - files required: read1.fastq and read2.fastq, fragment.txt Hba-1 mouse globin locus

#fastqc code as before

#trimming using trim-galore - need to use collate decorator as need to put in reads as pair for trimming
#FLASH (Fast Length Adjustment of SHort reads) is a very fast and accurate software tool to merge paired-end reads from 
#next-generation sequencing experiments. FLASH is designed to merge pairs of reads when the original DNA fragments are shorter than 
#twice the length of reads. The resulting longer reads can significantly improve genome assemblies. 
#@collate is for collating pairs or groups; @merge is to merge all samples, irrespective of grouping

import sys
from cgatcore import pipeline as P
from ruffus import *

P.get_parameters('capturec_pipeline.yml')

@follows(mkdir('fastqc'))
@transform('*.fastq.gz', regex(r'(.*).fastq.gz'),r'fastqc/\1_fastqc.html')
def qc_reads(infile, outfile):
    statement = 'fastqc -q -t %(threads)s --nogroup %(infile)s --outdir fastqc'
    P.run(statement,
          job_queue  = P.PARAMS['queue'],
          job_memory = P.PARAMS['memory']
          job_threads = P.PARAMS['threads'])

@follows (mkdir('trim'))
@collate('*.fastq.gz', regex(r'(.*)_[1-2].fastq.gz'), r'trim/\1_1_val_1.fq.gz')
def trim(infiles, outfile):
    ''' Trim fastq files'''
    fq1, fq2 = infiles
    cmd = '''trim_galore --paired %(fq1)s  %(fq2)s -o trim --cores %(threads)s'''
    P.run(cmd, 
          job_queue=P.PARAMS['queue'], 
          job_threads=P.PARAMS['threads'],
          job_memory=P.PARAMS['memory'])
    
@follows(mkdir('flash'), trim)
@collate ('trim/*.fq.gz', regex(r'trim/(.*)_[1-2]_.*'), r'flash/\1_extendedFrags.fastq.gz')
def combine_reads (infiles, outfile):
    read1, read2 = infiles
    out_prefix = outfile.replace('_extendedFrags.fastq.gz','') #otherwise will call file out.extendedFrags.fastq.gz, want read.extendedFrags.fastq.gz
    statement = '''flash %(read1)s %(read2)s -t %(threads)s -o %(out_prefix)s -z'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_threads=P.PARAMS['threads'],
          job_memory=P.PARAMS['memory'])
    
@follows(mkdir('digest'))          
@transform(flash_reads, regex(r'flash/(.*)_extendedFrags.fastq.gz'), r'digest/\1_digest.fastq.gz')
def digest(infile, outfile):
    '''in silico restriction enzyme digest'''
    statement = '''python digest_fastq.py -i %(infile)s -o %(outfile)s ''' #Alistair's script (see bottom)
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_memory=P.PARAMS['memory'])
   
@follows(mkdir('sam'))
@collate(digest, regex(r'digest/(.*).fastq.gz'), r'sam/\1.sam')
def align_reads(infiles, outfile):
    ''' Aligns digested fq files using bowtie2, output is .sam'''
    options = ''
    if P.PARAMS['bowtie2_options']:
        options = P.PARAMS['bowtie2_options'] #setting --reorder as an option
    cmd = '''bowtie2 -x %(bowtie2_ref)s -U %(infile)s -p %(threads)s %(options)s -S %(outfile)s'''
    P.run(cmd, 
          job_queue=P.PARAMS['queue'], 
          job_threads=P.PARAMS['threads'])
   
@follows(mkdir('ccanalyser'))
@transform(align_reads, regex(r'sam/(.*).sam'), r'ccanalyser/\1.sam') 
def ccanalyser (infile, outfile): #Jelena's script
    statement = '''perl CCanalyser2.pl -f %(infile) -r %(ccanalyser_ref)s #ref is DpnII cut sites throughout genome
    --genome %(ccanalyser_genome)s -o %(ccanalyser_oligo)s -s miseq  #genome used, oligo coordinates (pull down from at least one end of fragment) and sample name
    --pu %(ccanalyser_pu)s --pf %(ccanalyser_pf)s''' 
    P.run(statement,
          job_queue=P.PARAMS['queue'], 
          job_memory=P.PARAMS['ccanalyser_memory'])
    
#separate python script for in silico digest written by Alistair
import argparse
import os
import sys
import pysam
import gzip
import re

p = argparse.ArgumentParser()
p.add_argument('-i', '--input_fn', help='fastq file to parse')
p.add_argument('-o', '--output', help='output file prefix', default='output.fq.gz')
p.add_argument('-s', '--cutsite', help='Sequence of restriction site', 
               default='GATC')
args = p.parse_args()


def main():
    
    with gzip.open(f'{args.output}', 'wb') as w:
        
        cut_site = re.compile(f'{args.cutsite}')
        
        for record in pysam.FastqFile(args.input_fn):
            fragments = cut_site.split(record.sequence)
            lengths = [len(frag) for frag in fragments]
    
            last_slice = 0
            for ii, length in enumerate(lengths):
                current_slice = last_slice + length
                
                w.write(f'@{record.name}:PE:1:{ii}\n'.encode())
                w.write(f'{record.sequence[last_slice:current_slice]}\n'.encode())
                w.write('+\n'.encode())
                w.write(f'{record.quality}\n'.encode())
                
                last_slice = current_slice

if __name__ == '__main__':
    main()
