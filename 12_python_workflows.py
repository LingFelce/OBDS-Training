#12_Reproducible computational biology workflows 30th September 2019

"""
Perform steps in parallel for multiple input files
Pipeline scripts (bash, perl, python), makefiles, GUI workbenches (free and commercial), programmatic workflow systems
github.com/pditommaso/awesome-pipeline
Ruffus - open source workflow library for Python
CGATcore - wrapper for ruffus, can add parameters eg in a .yml configuration file so don't need to retype code, sent jobs to cluster in parallel
Building & running CGAT pipelines:
	1. write pipline script (series of python functions/tasks, import ruffus - task dependencies, and CGATcore - )
	2. write configuration file ((pipeline.yml)
	3. run from command line
	4. examine output
Defining tasks/functions - split, transform, merge
Python decorators - @symbol - describes input and output files e.g:
	@split(<input files>'all.gtf',<output files>'chr*.gtf') 
	def split_chrom(input, output):
	@transform(split_chrom, suffix(<input files>), <output files>) 
	def count_genes(input, output):
	@merge(count_genes, suffix(<input files>), <output files>)
General pipeline structure: documentation, import section, read parameters (.yml), define tasks (functions), _main_ (CGATcore command - run from command line)
Can have different conda environment for each section of pipeline
Look at slide for example of inside a task - very detailed!
CGATcore - uses actions to interact with ruffus pipeline 
	$python <script>mapping.py <action>plot <target function name>buildBigWig 
	$python pipeline_test.py make split_chrom -v 5
- other actions: show, plot (schematic of pipeline), make (run the pipeline), touch, config. Verbosity level (how much info for log file - recommend 5)
Grouping tasks - pipelines can have many independent branches, helpful to add dummy ruffus functions after each section of pipeline and "full"at end to run whole pipeline
	@follows(mapping,qc,views,duplication)
	def full():
	$python mapping.py make full -v 5
Pipeline status - presence of files, time stamps
Decorators can only be split, transform or merge (for now)
"""

#Exercise
#Write split-transform-merge pipeline
#Split file by chromosome
#Count number of transcripts for each chromosome
#Read all count files, calculate average and write to file

#$zcat genes.gtf.gz | cut -f1 | uniq #look at first column (chromosomes), see how many unique occurrences on command line

import gzip
from ruffus import *
from cgatcore import pipeline as P
import sys
import statistics

def split_chrom(infile): #define function
    with gzip.open(infile, 'rt') as inf: #read zipped file as text
        previous_chr = ''
        for i, line in enumerate(inf): #index and line - 2 variables
            col_list = line.split('\t') #set columns
            current_chr = col_list[0] #chromosome number/name is in first column
            if current_chr == previous_chr:
                outfile.write(line)
            elif i == 0: #if first line open file for first time
                outfile = open (f'{current_chr}.gtf', 'w')
                outfile.write(line)
	    else: #covers first line and subsequent chromosome number changes
                outfile.close() #not first line - close previous file
		outfile =  open(f'{current_chr}.gtf','w') #write as gtf file named after chr 
                outfile.write(line)     
            previous_chr = current_chr #resets chromosome number
                
split_chrom('test.gtf.gz') #run function

@transform(split_chrom, suffix('.gtf'),'.counts') #decorator/wrapper
def count_genes(infile, outfile):
    statement = '''wc -l %(infile)s > %(outfile)s''' #counts no. of lines, s specifies string. Each line = transcript
    P.run(statement) #if on cbrg have to specify job queue. Will run and save to single file per chromosome
	#or can do len(chr1.gtf) - count number of lines in file


@merge(count_genes,'all.average')
def average (infiles, outfile): #each count file has 2 items e.g 100 (no. of transcripts) chr1.gtf (original file name)
    total_counts = {} #create dictionary
    for infile in infiles:
        with open (infile) as inf:
            count, chrom = inf.read().strip().split(' ') #.read = reading the line, .strip = taking away white space, .split = setting 2 variables
            total_counts[chrom] = int(count) #key = chrom, count (make integer) = value
    median = statistics.median(total_counts.values()) #calculate median from integers only
    with open(outfile,'w') as count:
        for key, value in total_counts.items():
            count.write(f'{key}\t{value}\n') #write dictionary
        count.write(f'Median\t{median}\n') #write median

if __name__ == "__main__":
    sys.exit( P.main(sys.argv) ) #if have this at bottom can run in command line eg $python workflow.py make for full pipeline

#If want to set input file using .yml then make .yml file with (save .yml file in same directory):
gtf = /ifs/obds-training/lingf/obds/devel/OBDS_Training_Sep_2019/genes.gtf.gz

#Then under import statements
P.get_parameters('workflow.yml')

#And have to change decorators e.g
@split(P.PARAMS['gtf'], 'chr*.gtf') #end up with dictionary, gtf is key and value is file path
