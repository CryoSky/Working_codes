#!/bin/bash
# Written by Shikai Jin on 2018-01-18, submit 20 jobs simultaneously for specified RAMA value with random velocity seed number.
# Example in Linux: ./batch_create_directory.sh T0784-3u6g-HO+ER 0.25
dir=$1 # Notice: there is no space before and after =
rama=$2
task=`echo "${dir}" | cut -d \- -f 1`

for i in $(seq 1 20)
do
  mkdir ${dir}_${rama}_${i} 
  cp -r ./${dir}/. ./${dir}_${rama}_${i}/

  sed -i -e "/[Rama]/{n;s/2.0/${rama}/g}" ./${dir}_${rama}_${i}/fix_backbone_coeff.data
  sed -i -e "/velocity/{s/2349852/$RANDOM/g}" ./${dir}_${rama}_${i}/${task}-HO.in
  sed -i -e "`cat ./${dir}_${rama}_${i}/${task}-HO.in | grep -n "langevin" | sed -n 2p | cut -d : -f1`s/500/800/g" ./${dir}_${rama}_${i}/${task}-HO.in
  
  cd `dirname $0`/${dir}_${rama}_${i}
  echo `pwd`
  sbatch ${task}-HO.slurm
  cd ..
done
