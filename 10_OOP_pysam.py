#10_Object oriented programming - Python for Computational Genomics 26th September 2019

"""
Procedural v object oriented (procedural - write functions for each sub-problem, task oriented; object oriented focus on data not logic, use objects to model real-world concepts)
Pipelines use a mixture of both types of programming
classes represent real world objects, contain data (attributes) and have functions (methods)
Classes communicate with each other by calling methods (methods are functions are associated with class)
Inheritance - new classes inherit code, single, multi-level, hierarchal structure of classes allows for efficient code reuse (so don't have to make edits to same bit of code repeated)
Encapsulation improves security (can't tamper with code), keeps private
Abstraction - main function or call or method will stay the same (like write etc) but background code can change due to updates etc - won't see that, hidden
Polymorphism - one interface, many forms?
Use Pandas to make dataframes in Python (like R)
Computational genomic workflow: read standard file format (FASTA/FASTQ-sequence data; SAM/BAM/CRAM-genomic alignments; BED/GTF/GFF-genomic annotations; VCF-variant calls), process data
(quality control stats, trimming reads, genomic alignment, remove duplicates, count reads, overlapping intervals, extending intervals, filtering, annotation), write standard file format
(libraries - Pysam, HT Seq, PyBedtools)
Pysam - lightweight Python wrapper of C library HTSlib, used by Samtools. Use with genomic file formats. Google Pysam read the docs for help.

"""
#pandas example
import pandas as pd
df = pd.DataFrame()
print(df) #prints empty dataframe

#FASTQ in Pysam
with pysam.FastxFile(filename) as fh: #fh object of class FastxFile
	for entry in fh: #object of class Fastqproxy #loop, iterate through each read entry
		print(entry.name) #public attributes of class FastqProxy #object within object
		print(entry.sequence)
		print(entry.comment)
		print(entry.quality)

with pysam.FastxFile(filename) as fin, #maybe use to subset a FASTQ file?
open (out_filename, mode='w') as fout:
	for entry in fin:
		fout.write(str(entry)) #str method of class FastqProxy) #writing every entry out to a file

#SAM in Pysam
import pysam
samfile = pysam.AlignmentFile("ex1.bam", "rb") #AlignmentFile class representing SAM file, r=read SAM, rb=read BAM, rc=read CRAM
iter = samfile.fetch("chr1", 1000, 2000) #fetch is method of class AlignmentFile, returns as iterator. Fetch all reads from chrom 1 overlapping position 1000-2000

for aln in iter: #class representing alignment
	if aln.is_paired: #iterate over file print paired alignments #is_paired is an attribute
		print(aln)

#AlignmentFile Class - represents SAM/BAM/CRAM file. Methods include:
count() #count reads in region
count_coverage() #coverage of region
find_introns() #find introns based on NNN in CIGAR
get_index_statistics() #number of mapped and unmapped reads
head() #first n alignments

#AlignedSegment Class - single line in SAM file - an alignment
#attributes include query_name, flag, reference_name etc (column headings in SAM file)
#methods include:
infer_read_length() #from SIGAR
get_tags() #optional tags
get_overlap(start,stop)

#Pileup command allows iteration over each base in specific region (for SNP, variant calling etc)
iter = samfile.pileup('chr1', 100, 200) #method of class AlignmentFile, returns an iterator
for base in iter:
	print(str(x))

#writing SAM files - like a dictionary
header = { 'HD': {'VN': '1.0'}, 'SQ': [{'LN': 1575, 'SN': 'chr1'}, {'LN': 1584, 'SN': 'chr2'}] }
with pysam.AlignmentFile(“out.bam”, "wb", header=header) as outf:
	a = pysam.AlignedSegment()
	a.query_name= "read_28833_29006_6945"
	a.query_sequence="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
	a.flag= 99
	a.reference_id= 0
	a.reference_start= 32
	a.mapping_quality= 20
	a.cigar= ((0,10), (2,1), (0,25))
	a.next_reference_id= 0
	a.next_reference_start=199
	a.template_length=167
	a.query_qualities= pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
	a.tags= (("NM", 1), ("RG", "L1"))
	outf.write(a)

#using samtools in command line (csamtools)
pysam.sort("-o", "output.bam", "ex1.bam") #in python - allows for more investigation
samtools sort -o output.bam ex1.bam #command line

#Pysam works with Tabix indexed files (BED-user defined/GFF/GTF)
#Tabix is generic indexer for TAB-delimited genome position files, gives random access to compressed files

tbx = pysam.TabixFile("example.bed.gz") #tbx object of class TabixFile

for row in tbx.fetch("chr11", 100, 200, parser = pysam.asBed()): #pysam.asGTF() or pysam.asTuple()*
	print("name is", row.name) #if BED or GTF parser used, fields are accessible by name

#Using Pysam to convert SAM to BED
import pysam
import argparse 

parser = argparse.ArgumentParser(description='converting a sam/bam to bed with pysam')
parser.add_argument('--input', help='input sam or bam file required argument')
parser.add_argument('--output', default='output.bed', 
                    help='output file required argument')
parser.add_argument('--min_mapq', help='Minimum mapping quality',
                     default=0, type=int)

args = parser.parse_args()

print(args.input)
print(args.output)
print(args.min_mapq)

with open (args.output, 'w') as bed:    
    with pysam.AlignmentFile(args.input, "r") as sambam:
        for read in sambam.fetch():
          if not read.reference_name == None: #remove unmapped reads
              if read.mapping_quality >= args.min_mapq: #remove low quality reads (user defined)
                  if read.tlen > 0: #configure strandness
                      strand = '+'
                  elif read.tlen < 0:
                      strand = '-'
                  
                  bed.write(f"{read.reference_name}\t{read.reference_start}\t{read.reference_end}\t{read.query_name}\t{read.mapping_quality}\t{strand}\n")

#updated with additional options
import pysam
import argparse 

parser = argparse.ArgumentParser(description='converting a sam/bam to bed with pysam')
parser.add_argument('--input', help='input sam or bam file required argument')
parser.add_argument('--output', default='output.bed', 
                    help='output file required argument')
parser.add_argument('--min_mapq', help='Minimum mapping quality - input integer',
                     default=0, type=int)
parser.add_argument('--truncate', action="store_true", help='Truncate read coordinates to first base - leave out if don't want to truncate')
parser.add_argument('--padding', help='padding length to add to aligment/fragment - input integer',
                     default=0, type=int)

args = parser.parse_args()

print(args.input)
print(args.output)
print(args.min_mapq)
print(args.truncate)


#to be able to input BAM/SAM/CRAM file
if args.input.split('.')[-1] == 'bam':
    filetype = 'rb'
elif args.input.split('.')[-1] == 'sam':
    filetype = 'r'
elif args.input.split('.')[-1] == 'cram':
    filetype = 'rc'
else:
    print("Error: file needs to be .bam, .sam, .cram to work!")

with open (args.output, 'w') as bed:    
    with pysam.AlignmentFile(args.input, filetype) as sambam:
        for read in sambam.fetch():
            if not read.reference_name == None: #check if mapped
                if read.is_proper_pair: #check if proper pair
                    if read.mapping_quality >= args.min_mapq: #filter out low quality mapped reads (user defined)
                        if read.tlen > 0: #set strandness
                            strand = '+'
                        elif read.tlen < 0:
                            strand = '-'
                        
                        if args.truncate:
                            if strand == '+':
                                new_start = read.reference_start
                                new_end = read.reference_start + 1
                            else:
                                new_start = read.reference_end - 1
                                new_end = read.reference_end
                            bed.write(f"{read.reference_name}\t{new_start}\t{new_end}\t{read.query_name}\t{read.mapping_quality}\t{strand}\n")
     
                        else: #padding and truncation are separate
                            bed.write(f"{read.reference_name}\t{read.reference_start - args.padding}\t{read.reference_end + args.padding}\t{read.query_name}\t{read.mapping_quality}\t{strand}\n")
                    
#in command line can type:
python pysam_sambam_to_bed.py --help                    
python pysam_sambam_to_bed.py --input inputfile.sam --output outputfile.bed --min_mapq 1 --truncate <and/or> --padding 100
