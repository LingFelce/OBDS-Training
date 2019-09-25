#9_Biological algorithms - overlaps (reproducing BEDtools flow)

"""
BED tools bed and sam input, output will be how many times a read/alignment in sam file overlaps with feature (e.g gene) in bed file - use bed as reference
comparing sets of coordinates

Looking for - 100% overlap, partial overlap, no overlap

1. load sam and bed files
2. specify files from command line
3. iterate over files (for line etc)
4. compare for given line in bed line (outer loop), for inner loop extract reads at location
5. assess overlap, store counts and write out results once per bed feature


#bed	=====================
#sam
1.---------------
		2.--------
			3.-------------
4.-------------------------------------

1. single overlap
2. double internal
3. single overlap
4. double overlap

"""

parser = argparse.ArgumentParser(description='comparing sam and bed')
parser.add_argument('--sam', help='input sam file required argument') 
parser.add_argument('--bed', help='input bed file required argument')
parser.add_argument('--output', help='output file required argument')
parser.add_argument('--min_overlap', help='minimal length of bases overlapped', default=0, type=int') #find overlaps with user defined minimum eg 10 base pairs

args = parser.parse_args()

print(args.sam)
print(args.bed)
print(args.output)
print(args.min_overlap) #will show user defined minimum number of base pairs for overlap

with open (args.output, 'w') as output:

	with open (args.bed, 'r') as bed:
				
		for bed_index, bed_line in enumerate (bed): #gives line number and counts at same time
			chrom, start, end, name, score, strand = bed_line.split('\t') #naming the 6 columns in bed file
			start = int(start) #convert to integer for comparison
			end = int(end) #convert to integer for comparison
			#print(chrom, start, end, name, score, strand)
			count_overslaps = 0 #will still list negative counts
			if end - start >= args.min_overlap: #only check of minimum overlap if bed fragment is actually larger than specified overlap
				
				with open (args.sam, 'r') as sam: #open sam file within loop so that will continue iterating over file multiple times (rather than once and stop)
					

					for sam_index, sam_line in enumerate (sam):
						if not sam_line.startswith ("@"): #pass over headers
							col_list = sam_line.split('\t') #set a list
							sam_start = int(col_list[3]) -1 #sam start position is 1 based, bed start position is 0 based
							sam_end = len(col_list[9])+ int(col_list[3]) #calculate chromosome end
							sam_chrom = col_list[2] #chromosome name
							sam_length = sam_end - sam_start

							if chrom == sam_chrom and sam_length >= args.min_overlap: #compare chromsome (strings) first
								if sam_start >= start and sam_start <= end and sam_start <= end - args.min_overlap #covers 2. and 3. - see diagram above. Also checks if reads are x bp in bed file (min overlap)
									count_overlaps +=1 #counter
								elif sam_start <= start + args.min_overlap and sam_end >= start: #covers 1. and 4. - see diagram above. Also checks if reads are x bp in bed file (min overlap)
									count_overlaps +=1 #counter
				
				
			output.write(f'{bed_line.strip()}\t{count_overlaps}\n') #save to an output file (user specifies file name and path in command line directory). Strips anyway white space	
			#print(bed_line, count_overlaps)
