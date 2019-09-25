#8_Biological algorithms - genomic file formats 24th September 2019

#SAM file (sequence alignment) -> BED file (chromosomal information)

#don't have to git push all my files to OBDS github!
#format string - pull out terms that you want and put them along in a single string. Use f-string
#once convert to .bed file, will have unmapped reads at the top with * symbol - remove * lines from SAM file first
#failed multi-mapped reads (?) have * in CIGAR column - need to remove too.
#output fragments instead of reads - action store_true (default is false)
#can also add padding (don't know why but you can)
#only know read sequence, not fragment sequence

"""
#information we need for BED file and what column it's in SAM file (remember column 0)
chrom   - RNAME         - 3rd column
Start   - POS           - 4th column*
End     - POS + SEQ     - 4th column + 10th column*
name    - QNAME         - 1st column
score   - MAPQ          - 5th column
Strand  - +
Fragment- TLEN		- 9th column

*need to convert to integer and do something else!

"""
import argparse #to enable running Python code in terminal
#within Python -> run configuration per file -> command line options: specify destination and output folder above
#--sam /ifs/obds-training/lingf/obds/devel/OBDS_Training_Sep_2019/ERR177084.test1000.sam --bed (need to add new arguments!)
#can set parser arguments with default integer values, then specify value in command line option eg --min_mapq 1 (if default set to 0 in code below)
#so in terminal can type:
python sam_to_bed.py --sam /ifs/obds-training/lingf/obds/devel/OBDS_Training_Sep_2019/ERR177084.test1000.sam --bed padding.bed --fragments --min_mapq 1 --padding 100

parser = argparse.ArgumentParser(description='sam to bed')
parser.add_argument('--sam', dest='filename',
                    help='input sam file required argument') 
parser.add_argument('--bed', dest='output', default='output.bed',
                    help='output bed file')
parser.add_argument('--min_mapq', help='Minimum mapping quality',
                     default=0, type=int)
parser.add_argument('--fragments', help='output fragments rather than reads', 
		    action='store_true')
parser.add_argument('--padding',help='padding length to add to alignment/fragment',
		    default=0, type=int) #padding length user defined

args = parser.parse_args()
#print(args.filename)
#print(args.output)
#print(args.min_mapq)
#print(args.fragments)

#filename = '/t1-data/user/wchen/obds/devel/OBDS_Training_Sep_2019/ERR177084.test1000.sam' - not needed as already specified
#output = '/home/cgat/wchen/obds/week2/ERR177084.test1000.bed' - not needed as already specified

with open (args.output, 'w') as bed:    
    with open(args.filename, 'r') as sam:
        
	processed_lines = 0
        processed_headers = 0
        processed_reads = 0
        failed_mapping = 0
        failed_mapqthresh = 0
	fragment_count = 0
	        
	for line in sam:
            processed_lines += 1
            if line.startswith ("@"): #find line that starts with @, and do nothing
                processed_headers += 1
            else: #find the content
                processed_reads += 1
                col_list=line.split('\t')  #set a list
                if not col_list[2] == '*': #or if col_list[2] =! '*':
                    if int(col_list[4]) >= args.min_mapq:
                        start = int(col_list[3]) -1 - args.padding
                        if args.fragments:
				if int(col_list[8]) > 0:	
					end = int(col_list[8])+ int(col_list[3]) + args.padding #adding template length to chromosome start to get different chromosome end (template should be longer than read)
					bed_line = f"{col_list[2]}\t{start}\t{end}\t{col_list[0]}\t{col_list[4]}\t+\n"
                        		bed.write(bed_line) 
					fragment_count += 1 #count fragments
			else:						#if fail > 0 then will not go through to else.
				end = len(col_list[9])+ int(col_list[3]) + args.padding #if not output the read length instead
                        	bed_line = f"{col_list[2]}\t{start}\t{end}\t{col_list[0]}\t{col_list[4]}\t+\n"
                        	bed.write(bed_line) 
                    else:
                        failed_mapqthresh += 1
                else:
                    failed_mapping += 1
                    
print(f'Lines processed: {processed_lines}',
      f'Header lines processed: {processed_headers}',
      f'Alignments processed: {processed_reads}',
      f'Reads not mapped: {failed_mapping}', 
      f'Reads failing mapping quailty filter {failed_mapqthresh}',
      f'Fragment counts {fragment_count},
      sep='\n')





#repeated but with my useful comments (and Wentao's)
with open (args.output, 'w') as bed: #when looking at sam file, will save output into bed
    with open(args.filename, 'r') as sam:
        for line in sam:
            if line.startswith ("@"): #find line that starts with @, and do nothing
               pass
            else: #look at everything else
                col_list=line.split('\t')  #creates a list
                if not col_list[2] == '*': #or if col_list[2] =! '*': #unmapped reads in sam file have *; if line has no * then continue
                    start = int(col_list[3]) -1 #sam to bed offset -1 for start, have to convert to integer for arithmetic
                    end = len(col_list[9])+int(col_list[3]) #sequence length + chrom start (convert to integer) = chrom end
		    #form string made up of required columns with tab space \t and line break \n at end
                    bed_line = f"{col_list[2]}\t{start}\t{end}\t{col_list[0]}\t{col_list[4]}\t+\n"
                    bed.write(bed_line) #write string to bed file 
