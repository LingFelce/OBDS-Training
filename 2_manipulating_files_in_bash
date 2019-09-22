#redirecting input, output and error
command1 < file1 #input file1 to command1
command1 < (command2) #output of command 2 as file input to command1
command1 > file1 #standard output of command1 to file1
command1 >> file1 #append standard output of command1 to file1
command1 2> file2 #error output of command1 to file2
command1 1>&2 #standard output to same place as standard error
command1 2>&1 #standard error to same place as standard output
command1 > /dev/null #discard standard output of command1

#combining commands
command1; command2 #run cmd1 then cmd2
command1 && command2 #run cmd2 if cmd1 is successful
command1 || command2 #run cmd2 if cmd1 is not successful
command1 | command2 #Pipe stderr from cmd1 to stdin cmd2
command1 |& command2 #Pipe stdout from cmd1 to stdin cmd2

#loops
for i in {1..5}; do COMMAND; done
for (i=1;i<=10;i+=2); do COMMAND; done
for i in *.txt; do COMMAND; done
for x in *.bed; do wc –l $x; done #* Glob for files, use filename as variable

#regular expressions www.regexe.com and regex101.com

#search within files using grep (find things that match search term)
grep [options] regex [file...]
grep –i “exception” pipeline.log #ignore case
grep –c “chr1” p300.bed #print number of matches instead of lines (count)
grep –v “#” pedigree.vcf #invert match - prints lines that don't match

#manipulating text files in Linux
paste #displays the corresponding lines of multiple files side-by-side
paste file1.txt file2.tsv 

join #joining lines of two files on a common field
join file1.txt file2.tsv #merge files line by line

cut #Remove sections of each line of a file
cut -f3 file1.tsv #remove sections from each line of files. Extract only the third field from a tab delimited file
#(the default field delimiter is tab - \t)

sort #sorts files using one or more keys
sort file1.txt file2.txt file3.txt > sorted.txt
sort --key=1,1 --key=2n filename #start at field one, end at field one. Field 2 is 2nd sort key. n=numeric sort

uniq #reports or filters out repeated lines in a file
uniq filename #remove duplicate lines
uniq -c filename #prefix lines by number of occurances

tr #translate or delete characters

grep #search for lines which match a specified pattern

sed #a stream editor - text substitution. http://sed.sourceforge.net/sed1line.txt
sed s/regexp/replacement/[flags]
sed s/chr1/1/g
#d delete pattern space; g apply replacement to all matches; 2 replace 2nd occurance of pattern in line

#Print a specific line (e.g. line 42) from a file:
sed -n 42p <file>
#Extract every 4th line starting at the second line (extract the sequence from FASTQ file):
sed -n '2~4p' file.fastq
#Replace all occurrences of foo with bar in file.txt: (s = substitute)
sed 's/foo/bar/g' file.txt
#Trim leading whitespaces and tabs in file.txt:
$ sed 's/^[ \t]*//' file.txt
#Trim trailing whitespaces and tabs in file.txt:
sed 's/[ \t]*$//' file.txt
#Trim leading and trailing whitespaces and tabulations in file.txt:
sed 's/^[ \t]*//;s/[ \t]*$//' file.txt
#Delete blank lines in file.txt:
sed '/^$/d' file.txt

AWK #programming language designed for tabular text processing, arithmetic operations
#https://www.tutorialspoint.com/awk/index.htm
#https://www.shortcutfoo.com/app/dojos/awk/cheatsheet
awk ‘Pattern { action }’ <file>
awk ‘{ print $0 }’ <file> #$ denotes position $0 entire line, $1 first variable, $2 second variable
awk ‘($0 ~ /chr1/) {print $0}’ peaks.bed
awk ‘/chr1/ {print}’ peaks.bed
awk ‘{print $3-$2}’ peaks.bed
awk ‘{print $1 ”:” $2 ”-” $3}’ peaks.bed
awk ‘BEGIN{SUM=0}{SUM+=$3-$2}END{print SUM}’ peaks.bed

rev filename #reverses order of characters in every line

#transliteration
echo "lowercase letters" | tr a-z A-Z #convert text to uppercase
echo AATGATACGGCGA | rev | tr ATGC TACG #reverse complement
