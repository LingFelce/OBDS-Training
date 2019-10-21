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




# Base is default environment

# track what you've added (full history of packages that you can revert to)
conda list --revisions
