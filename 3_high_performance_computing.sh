#job profiling
top
qstat -j <job_id> | grep usage

qhost #computational resources available
qconf -shgrpl #show all host groups
qconf -shgrp <@name> #show host group conf.
qconf -sql #show list of queues
qconf -sq <queuename> #show queue config
qconf -spl #parallel environments
qhost -j #cluster
qstat -f -u "*" #cluster
qsub [sge-options] <batch-prog>[prog-args] #eg pipline.py
qrsh [sge-options][commands] #use qrsh to switch to computer node, so not doing too much processing on head node

#submit jobs to cluster
qsub –cwd # use current directory as working one
qsub –i <input.file> # redirect standard input
qsub –o <output.file> # redirect standard output
qsub –e <error.file> # redirect standard error
qsub –j y # redirect standard error to output
qsub –p <priority> # lowest: -1024, highest: 1023
qsub –N <job-name> # change job name
qsub –M <email> -m e # email me when done
qsub –q <queue-name> # select specific queue
qsub –v name=value # export env variables
qsub –hold_jid <job-id> # start job after <job-id>
qsub –l h_rt=1:30 #run time limit
qsub –l mem_free=2G #memory usage per slot
qsub –pe dedicated <n-slots> #shared memory multi-processing
qsub –pe mpi <n-slots> #message passing interface
qsub –t n[-m[:s]] <script.sh> #array jobs
qsub sge-job.sh #batch jobs - use submission script

#sge-job.sh script
#!/bin/bash
# The SGE batch system uses the current directory as working directory.
# Both files (output.dat and error.dat) will be placed in the current
# directory. The batch system assumes to find the executable in this directory.
#$ -cwd
# Redirect output stream to this file.
#$ -o output.dat
# Redirect error stream to this file.
#$ -e error.dat
# Send status information to this email address.
#$ -M sebastian.lunavalero@imm.ox.ac.uk
# Send an e-mail when the job is done.
#$ -m e
# For example an additional script file to be executed in the current
# working directory. In such a case assure that script.sh has
# execution rights (chmod +x script.sh).
./script.sh

qsub sge-serial.sh; qstat; cat serial.output #simple serial job
qsub -l h=cgat022 sge-serial.sh; cat serial output #select specific execution host
qsub -N myjob sge-parallel.sh #run parallel job with customised name

qdel <job-id>
