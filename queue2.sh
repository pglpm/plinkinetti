#!/bin/bash
#PBS -l walltime=150:00:00
#PBS -l nodes=1:ppn=12

#! General jobscript. Do not modify below
#!echo $SLURM_ARRAY_TASK_ID
#!CHANGE="./${TEXT}${SLURM_ARRAY_TASK_ID}"
#!cd $CHANGE
#!echo $CHANGE
#! NODES='wc -l <$PBS_NODEFILE'
## module load anaconda
## source activate r
Rscript job.R
#!cp $DATASET $SCRATCH
#!cd $SCRATCH
#!chkfile $OUTFILE
#!YourProgram $DATASET > $OUTFILE
