#$ -cwd

#$   -cwd

#$  -N Sim_multiseed

#$   -l   mem_free=2.0G

#$   -l   h_vmem=2G

#$  -t 1-200

#  -m  e  -M  ddey1@jhu.edu

#$ -j y

# -o $HOME/dagar/nimble_outputs/o.$JOB_ID.$TASK_ID

#$ -o output/o.$JOB_ID.$TASK_ID

# -pe local 1



date

echo "I am running on compute node:"

hostname


# Run the "R" program to read in the "plot1.r" script.


R CMD BATCH --no-save Sim_notourmod.R log/sim${SGE_TASK_ID}.log

echo "Finished running. Good workout!"