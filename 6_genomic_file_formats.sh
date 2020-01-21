##Genomics File Formats 20th September 2019

#Fastq md5 and file name (from ENA)

0d6d6d9075aa03c75b54d38845e99f02 ERR1755084_1.fastq.gz
0aa75e821d1a58a40422970b7182e754 ERR1755084_2.fastq.gz
f5944eece094f7bb2ff1d5c9d3954aee ERR1755087_1.fastq.gz
3ac3829e93c4e236df6125938451a9f3 ERR1755087_2.fastq.gz

md5sum -c hash.md5

#FASTQC - output is .html file that you can download using filezilla

fastqc --nogroup -t 12 -q ERR1755084_1.fastq.gz &

#HiSat2 for mapping and alignment

#--fr is first strand (strand specific)

hisat2 --fr -p 8 --quiet -x /ifs/mirror/genomes/hisat2/mm10 -1 ERR1755084_1.fastq.gz -2 ERR1755084_2.fastq.gz -S ERR1755084.sam &
hisat2 --fr -p 12 --quiet -x /ifs/mirror/genomes/hisat2/mm10 -1 ERR1755087_1.fastq.gz -2 ERR1755087_2.fastq.gz -S ERR1755087.sam &

#QC for alignment

samtools sort ERR1755084.sam -o ERR1755084.bam --output-fmt BAM
samtools sort ERR1755087.sam -o ERR1755087.bam --output-fmt BAM

samtools index ERR1755084.bam ERR1755084.bam.bai &
samtools index ERR1755087.bam ERR1755087.bam.bai &

picard CollectAlignmentSummaryMetrics R=/ifs/mirror/genomes/plain/mm10.fasta I=ERR1755084.bam O=ERR1755084.txt & 
picard CollectAlignmentSummaryMetrics R=/ifs/mirror/genomes/plain/mm10.fasta I=ERR1755087.bam O=ERR1755087.txt &

#MultiQC

multiqc . #(make sure in correct directory)

#to go in .bashrc to make sure multiqc runs

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

#Counting

featureCounts -s 1 -T 12 -a *.gtf.gz -o featureCounts.summary *.bam

#to open IGV genome browser in terminal in X2Go

module load bio/IGV

igv.sh


