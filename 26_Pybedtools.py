""""""""""""""""""""""
18th October 2019
Pybedtools wraps and extends bedtools and offers feature level manipulations from within Python.
https://daler.github.io/pybedtools/

Done in Jupyter Lab, see exported HTML

Three brief examples: https://daler.github.io/pybedtools/3-brief-examples.html

"""""""""""""""""""""
!date
!pwd
!conda env list | grep '*'
!python --version

import pybedtools
import pandas as pd

# create some bedtool objects using .bed files from bedtools tutorial
a = pybedtools.BedTool('cpg.bed')
b = pybedtools.BedTool('exons.bed')

# intersection with BedTool.intersect() method
a_and_b = a.intersect(b) 

#a_and_b is a new BedTool instance, pointing to temporary file

a.head()

b.head()

# shows all CpG islands that overlap exons
a_and_b.head()

# intersection using -u switch (True/False to indicate we want to see features in a that overlapped with b)
a_with_b = a.intersect(b, u=True)
a_with_b.head()
# only lists unique overlapped entries (collapsed form of above)

# save files with trackline - customisable, can add extra info here
c = a_with_b.saveas('intersection_of_cpg_and_exons.bed', trackline='track name="cpg and exons"')
# or d = a_with_b.moveto('another_location.bed') if don't want to save trackline

print(c.fn) # .fn is temporary file
# printing will show track name

# open file to show track line
c.head()

# convert intersection file into dataframe using pandas
df_c = c.to_dataframe()
df_c

import pybedtools
exons = pybedtools.BedTool('exons.bed')
cpg = pybedtools.BedTool('cpg.bed')
exons_with_cpg = exons.intersect(cpg, u=True)
exons_with_cpg.head() # lists which exons have CpG islands in

x1 = exons.intersect(cpg)
x1.head()

x2 = exons.intersect(a=exons.fn, b=cpg.fn)
x2.head()

x3 = exons.intersect(b=cpg.fn)
x3.head()

x4 = exons.intersect(cpg, a=exons.fn)
x4.head()

# four ways of getting same output
x1 == x2 == x3 == x4

# chain methods together (pipe)
x1 = a.intersect(b, u=True)
x2 = x1.merge() # merge overlapping intervals

# combine into one line
x3 = a.intersect(b, u=True).merge()

# even longer command on separate lines
x4 = a\
    .intersect(b, u=True)\
    .saveas('a-with-bed.bed')\
    .merge()\
    .saveas('a-with-b-merged.bed')

# quicker way to do intersection
x5 = a.intersect(b, u=True)
x6 = a + b
x5 == x6

# find no overlaps
x7 = a.intersect(b, v=True)
x8 = a - b
x7 == x8

# if want to do something after addition/subtraction
x9 = (a + b).merge()
x2 == x3 == x9

# create BedTool object (again)
a = pybedtools.BedTool('exons.bed')

# access Intervals by indexing BedTool object
feature = a[0]
feature

# slicing
features = a[1:3]
list(features)

# printing converts to original line from file
print(feature)

feature.chrom
feature.start

# BED is 0-based, GFF, GTF and VCF are 1-based

# create GFF Interval from scratch
gff = ["chr1",
      "fake",
      "mRNA",
      "51", # start is greater than start for BED below
      "300",
      ".",
      "+",
      ".",
      "ID=mRNA1;Parent=gene1;"]
gff = pybedtools.create_interval_from_list(gff)

# make sure keys are sorted - only for testing
gff.attrs.sort_keys = True
print(gff)

# create corresponding BED Interval; will need to subtract 1 from start
bed = ["chr1",
      "50",
      "300",
      "mRNA1",
      ".",
      "+"]
bed = pybedtools.create_interval_from_list(bed)
print(bed)

# format auto-detected based on position of chrom/start/stop provided
gff.file_type

# format auto-detected based on position of chrom/start/stop provided
bed.file_type

bed.start == gff.start == 50

# add some new attributes
gff.attrs['Awesomeness'] = "99"
gff.attrs['ID'] = 'transcript1'
gff.attrs.sort_keys = True
assert gff.attrs.sort_keys # what is assert?
print(gff)
x2 == x3
