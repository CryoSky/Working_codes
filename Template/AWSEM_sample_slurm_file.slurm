#!/bin/bash

#For introduction of NOTS and DAVinCI, please check https://docs.rice.edu/confluence/display/CD/Getting+Started+on+NOTS and https://docs.rice.edu/confluence/display/CD/Getting+Started+on+DAVinCI

#SBATCH --account=ctbp-common
#SBATCH --partition=ctbp-common #You need to specify the name of the condo account to use a condo on the cluster.
#SBATCH --job-name=T0815-exHOER #Assigns a job name.
#SBATCH -o trial.out 
#SBATCH -e trial.err
#SBATCH --ntasks=1 #Required: The maximum number of tasks per job. Usually used for MPI jobs.
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=1024M
#SBATCH --mail-user=1964569056@qq.com
#SBATCH --mail-type=FAIL
#SBATCH --export=ALL #Optional:  Exports all environment variables to the job.
#SBATCH --nodes=1

#############
#date
#echo My job ran on: $SLURM_NODELIST
#echo My jobid was $SLURM_JOBID
#echo My taskid was $SLURM_ARRAY_TASK_ID (when multiple jobs)https://docs.rice.edu/confluence/display/CD/Using+Job+Arrays+to+submit+large+number+of+jobs+via+SLURM
#echo My job array jobid was $SLURM_ARRAY_JOB_ID
#echo My submission directory was $SLURM_SUBMIT_DIR
#module purge
#module load XXX
#The command used in Linux system.
#############
