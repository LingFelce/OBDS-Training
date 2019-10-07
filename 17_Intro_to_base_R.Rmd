#Introduction to base R - 7th October 2019

R Studio - interactive development environment

Tools > Global options > Git/SVN > Type in Git details > Enable version control interface for R Studio projects
Should then get Git tab on top right for commiting etc.
New files will appear in environment, tick to add then press commit. Can type commit message. Can push separately once 
set up with online repository.
Can also set up project by cloning from Git repository (eg template files)

Packages for R: CRAN (package repository, general use), Bioconductor (bioinformatics, more documentation/vignettes)
Github also has some packages but not as rigorously tested as CRAN or Bioconductor, so be wary.

Install via Conda on command line if want to keep track of package versions etc for example:
$conda install r-tidyverse
$conda install bioconductor-deseq2

Better to use <- notation not = for assigning variables. (alt + - as a shortcut)

Attributes associated with each value eg names or other metadata - similar to python dictionaries
>attributes(vector_name)

R is nested functions, can use tidyverse to make code look nicer e.g tibbles are tidyverse dataframes (simple and better functionality)

```{r examples}
df[1] #1st column only
df[1,] #1st row only
df[5:10,1:2] #rows 5:10, columns 1:2

df[grep("TSS", df&=$column1),] #get rows containing TSS in column 1

df[df$column1 %in% c('promoter', 'enhancer'),] #rows that are enhancer or promoter values in column 1
df[df$column1 == 'promoter' | df$column1 == 'enhancer',] #rows that are enhancer or promoter values in column 1
```
