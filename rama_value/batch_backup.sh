#!/bin/bash
# Written by Shikai Jin on 2018-01-30, this file is used to backup required files for data analysis 
# for one protein (totally 100 runs) and reappear the task
# Example in Linux: ./batch_backup.sh T0784-3u6g-HO+ER
dir=$1
backupdir='/scratch/sj52/backup'
pathway=$(cd `dirname $0`;pwd)
ramalist="0.25 0.5 1.0 1.5 2.0"
task=`echo "${dir}" | cut -d \- -f 1`

mkdir -vp ${backupdir}/${task} #-v means print a message, -p means make parent directories as needed
for rama in $ramalist
do
  for i in $(seq 1 20)
    do
      cd "${pathway}"/${dir}_${rama}_${i}
      mkdir -vp ${backupdir}/${task}/${dir}_${rama}_${i}
      cp energy* ${backupdir}/${task}/${dir}_${rama}_${i}/
      cp dump* ${backupdir}/${task}/${dir}_${rama}_${i}/
      cp wham* ${backupdir}/${task}/${dir}_${rama}_${i}/
      cp *.in ${backupdir}/${task}/${dir}_${rama}_${i}/
      cp *.slurm ${backupdir}/${task}/${dir}_${rama}_${i}/
    done
done
