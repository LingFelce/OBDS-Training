# Conda environments
# general purpose packaging system (Python, R, etc)
# Cross platform (Windows, Linux)
# Create separate software environments
# Some packages have different requirements so cannot be used in same environment
# e.g macs2 (peak caller software) requires python 2.7, current version is python 3.6
# Can easily share environment versions with others

# Conda docs: https://conda.io

# Conda workshop: https://github.com/OBDS-Training/Conda_Workshops/blob/master/1_Conda_intro.md

# log into cluster using ssh
ssh cgath1

# move onto computer node on cluster
qrsh

# go to working directory
cd /ifs/obds-training/lingf
mkdir conda
cd conda

#install copy of conda install script
curl -o Miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

# run install script to install conda
bash Miniconda.sh -b -p obds_conda_install

# add location to path variable so that terminal knows where to find conda software so we can use it
# activate conda installation
source /full/file/path/to/where/you/have/installed/obds_conda_install/etc/profile.d/conda.sh

# activate base environment to move into default conda environment
conda activate base

# test if source command has worked
conda --help
which conda
conda info

# can use conda to search for software packages, but conda needs channels to get software
# add channels in correct order!
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

# check channels installed
conda info

# list currently installed packages
conda list

#check packages up to date
conda update --all

# install a new package (if from default, conda-forge or bioconda channel - check conda website)
conda install <packagename>

# check if installed correctly
conda list <packagename>
which <packagename>
<packagename> --help

# remove package
conda remove <packagename>

# useful to have separate environments to share/export
conda env -h
conda env list 
conda env list -h

# create new environment with desired package and activate
conda create -n <environmentname> <packagename>
conda activate <environmentname>

# check current activate conda environment
conda env list

# check if software installed properly and what verison of Python have installed
<environmentname> -h
python --version

# can export conda environment and check .yml file
conda env export -n <environmentname> > env.yml
cat env.yml

# can recreate environment in another conda installation
conda env create -n <environmentname> -f env.yml

# Base is default environment every time you install conda, contains latest version of Python 
# and a few basic Python packages

# track what you've added (full history of packages that you can revert to)
conda list --revisions

# currently using obds_full_env.yml file, can use this to create new environment, containing:
# Python & associated libraries
# Python - version 3.6 (you should specify version number in the yaml file)
    # numpy  - (a python package for doing fast mathematical calculations and manipulations)
    # pandas - (a python package for making/using dataframes)
    # matplotlib - (a python package for plotting)
    # seaborn -  (a much prettier python package for plotting)
    # scipy - (a collection of python packages for data analysis, includes ipython, pandas etc)
    # spyder - an interactive development environment (IDE) for python similar to Rstudio 
    # ruffus - a python pipelining program that we will use to write pipelines
    # cgatcore - a library from cgat to make pipelines usable with a computer cluster
    # pysam - a python package for working with bam/sam alignment files 
    # pybedtools - a python wrapper for bedtools meaning you can use bedtools functionality in python scripts
    # drmaa - for the management of submitting jobs to the cluster
    # ggplot  - python version of ggplot
    # jupyter - interactive notebooks for python - install jupyterlab instead
# R  - note that r is called `r-base` in conda - search it and check you can find it. Might have to Google to get correct names
    # rstudio - IDE for R 
    # tidyverse library - family of r packages that have usefull data processing functionality
    # deseq2 - r bioconductor statistical package for differential expression
    # edger - r bioconductor statistical package for differential expression
    # seurat - r Cran package for single cell analysis 
    # goseq - r bioconductor package for gene ontology analysis
    # gsea - r package for gene set enrichment analysis
# Bioinformatics software  
    # fastqc - QC of fastq raw sequence files
    # multiqc - collects summary statistics from other bioinformatic programs
    # hisat2 - super quick read aligner (mapper) for spliced sequencing reads
    # bowtie2 - slower read aligner for unspliced sequencing reads
    # samtools - manipulate bam/sam alignment files 
    # picard - QC of alignment files 
    # subread - counting of reads in features
 
 # modify .bashrc so that can use alias (keyboard shortcut) to activate conda environments
# load conda env
DEST="/ifs/obds/${USER}/conda"
source /ifs/obds-training/lingf/conda/obds_conda_install/etc/profile.d/conda.sh
alias cc='conda activate base && conda activate obds_full_env'

