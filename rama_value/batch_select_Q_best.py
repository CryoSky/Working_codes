#! usr/bin/env python3
# Written by Shikai Jin, 2018-01-26, modified from Xun Chen's unpublished code
# This file is used for extracting Q max value and steps from 20 cycles, save on a new txt file.
# Example in Linux: python batch_select_Q_best.py /scratch/sj52/RAMA_Tempfree/T0784 T0784-3u6g-HO+ER_0.25

import os
import sys


def selectmaxQ(wham):  # extract the max Qw value from wham.dat file
    with open (wham, 'r+') as fp:
        max = 0.0
        steps = 0.0
        for line in fp.readlines():
            if line.split()[1] != 'timestep':
                curr = float(line.split()[1])
                if max < curr:
                    max = curr
                    steps = float(line.split()[0])
    return (max, steps)

def savemaxQ(Qlist, dir, prefix): # save max Qw value to "prefix"_maxQ.txt
    with open ('%s/data_analysis/%s_maxQ.txt' %(dir, prefix), 'w+') as fp:
      for i in range(20):
        fp.write(str(Qlist[i]) + '\n')
      fp.close()
    
def savesteps(steplist, dir, prefix): # save max Qw value to "prefix"_maxQ.txt
    with open ('%s/data_analysis/%s_steps.txt' %(dir, prefix), 'w+') as fp:
      for i in range(20):
        fp.write(str(steplist[i]) + '\n')
      fp.close()

    
def main():
    awsemdir = '/scratch/sj52/Rotation3/awsemmd-master/results_analysis_tools'
    dir = sys.argv[1]
    prefix = sys.argv[2]
    list_Q = [0 for j in range(20)] # initialize list for saving Q value
    list_step = [0 for j in range(20)] # initialize list for saving steps

# create "current directory"/data_analysis/"protein prefix"_Best_structure as the saving directory
    if os.path.exists('%s/data_analysis' %(dir)) == False:
        os.mkdir('%s/data_analysis' %(dir))
    os.chdir('%s/data_analysis' %(dir))
    os.system('mkdir %s_Best_structure' %(prefix)) # two ways to create directory
    datadir = dir + '/' + 'data_analysis' + '/' + prefix + '_Best_structure'
    print(datadir)

    for i in range(1, 21):
        branch = dir + '/' + prefix + '_' + str(i)
        wham = branch + '/wham.dat'
        pdbfile = prefix + '_' + str(i) + '_best'
        (list_Q[i-1], list_step[i-1]) = selectmaxQ(wham) # file start from 1 while the list start from 0
        snapshot = round(list_step[i-1])/1000 # snapshot means frame, not the timestep, here I set 1000 timestep per frame, also noticed the snapshot is float and should be convert to int to use
        print(snapshot)
        # Noticed this maybe not effective when you input the restart file
        os.chdir('%s' %(datadir))
        print(branch)
        os.system('python2 %s/BuildAllAtomsFromLammps.py %s/dump.lammpstrj %s.pdb %s' % (awsemdir, branch, pdbfile, int(snapshot))) # Build file coded by Python2
        
    savemaxQ(list_Q, dir, prefix)
    savesteps(list_step, dir, prefix)

main()
