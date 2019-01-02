#! usr/bin/env python3
# Written by Shikai Jin, 2018-01-26, modified from Xun Chen's unpublished code
# This file is used for extracting Q max value and steps from 20 cycles, save on a new txt file.
# Example in Linux: python batch_select_Q_best.py /scratch/sj52/RAMA_Tempfree/T0784 T0784-3u6g-HO+ER_0.25

import os
import sys

'''
def selectmaxQ(wham):  # extract the max Qw value from wham.dat file
    with open (wham, 'r+') as fp:
        max = 0.0
        steps = 0.0
        for line in fp.readlines():
            if line.split()[1] != 'Fix':
                curr = float(line.split()[1])
                if max < curr:
                    max = curr
                    steps = float(line.split()[0])
    return (max, steps)
'''

def selectmaxQw(qw_list):  # extract the max Qw value from wham.dat file
    with open (qw_list, 'r+') as fp:
        max = 0.0
        steps = 0
        maxstep = 0
        for line in fp.readlines():
            steps = steps + 1
            if line.split()[0] != 'Fix' and steps > 3:
                curr = float(line.split()[0])
                if max < curr:
                    max = curr
                    maxstep = steps
    maxstep = maxstep - 3
    print("The max step is %s" %maxstep)
    return (max, maxstep)


def selectenergy(energy,steps):  # extract the max Qw value from wham.dat file
    qgo = 0
    er = 0
    #print (steps)
    with open (energy, 'r+') as fp:
        for line in fp.readlines():
            if line.split()[0] != 'Step':
                filestep = float(line.split()[0])
                #print(filestep)
                if filestep == steps:
                    print ("everything is ok")
                    qgo = float(line.split()[-2])
                    er = float(line.split()[11])
    print (er)
    print (qgo)
    return (qgo, er, steps)

def savemaxQ(Qlist, dir, prefix): # save max Qw value to "prefix"_maxQ.txt
    with open ('%s/data_analysis/%s_maxQ.txt' %(dir, prefix), 'w+') as fp:
      for i in range(20):
        fp.write(str(Qlist[i]) + '\n')
      fp.close()

def savemaxQgo(Qlist, dir, prefix): # save max Qw value to "prefix"_maxQ.txt
    with open ('%s/data_analysis/%s_maxQgo.txt' %(dir, prefix), 'w+') as fp:
      for i in range(20):
        fp.write(str(Qlist[i]) + '\n')
      fp.close()

def savemaxer(Qlist, dir, prefix): # save max Qw value to "prefix"_maxQ.txt
    with open ('%s/data_analysis/%s_maxer.txt' %(dir, prefix), 'w+') as fp:
      for i in range(20):
        fp.write(str(Qlist[i]) + '\n')
      fp.close()
    
def savesteps(steplist, dir, prefix): # save max Qw value to "prefix"_maxQ.txt
    with open ('%s/data_analysis/%s_steps.txt' %(dir, prefix), 'w+') as fp:
      for i in range(20):
        fp.write(str(steplist[i]) + '\n')
      fp.close()

    
def main():
    awsemdir = '/home/sj52/awsemmd-amylometer/results_analysis_tools'
    dir = sys.argv[1]
    prefix = sys.argv[2]
    list_Q = [0 for j in range(20)]
    list_Qgo = [0 for j in range(20)] # initialize list for saving Q value
    list_step = [0 for j in range(20)] # initialize list for saving steps
    list_er = [0 for j in range(20)]



# create "current directory"/data_analysis/"protein prefix"_Best_structure as the saving directory
    if os.path.exists('%s/data_analysis' %(dir)) == False:
        os.mkdir('%s/data_analysis' %(dir))
    os.chdir('%s/data_analysis' %(dir))
    os.system('mkdir %s_Best_structure' %(prefix)) # two ways to create directory
    datadir = dir + '/' + 'data_analysis' + '/' + prefix + '_Best_structure'
    print(datadir)

    for i in range(1, 21):
        branch = dir + '/' + prefix + str(i)
        qw_list = branch + '/qw_list'
        #wham = branch + '/wham.dat'
        energyfile = branch + '/energy.log'
        pdbfile = prefix + '_' + str(i) + '_best'

        #dir = os.getcwd()

        os.chdir(branch)
        os.system("pwd")
        os.system(
            "python /home/sj52/awsemmd-amylometer/create_project_tools/CalcQValue_structure.py ../5fjl_A_selected.pdb dump.lammpstrj qw_list 0")
        os.chdir(dir)


        (list_Q[i-1], list_step[i-1]) = selectmaxQw(qw_list) # file start from 1 while the list start from 0
        #snapshot = round(list_step[i-1])/1000 # snapshot means frame, not the timestep, here I set 1000 timestep per frame, also noticed the snapshot is float and should be convert to int to use
        snapshot = list_step[i-1]
        print(snapshot)
        # Noticed this maybe not effective when you input the restart file
        (list_Qgo[i - 1],list_er[i - 1], list_step[i - 1]) = selectenergy(energyfile, snapshot*1000)
        os.chdir('%s' %(datadir))
        print(branch)
        os.system('python2 %s/BuildAllAtomsFromLammps.py %s/dump.lammpstrj %s.pdb %s' % (awsemdir, branch, pdbfile, int(snapshot))) # Build file coded by Python2

    savemaxQ(list_Q, dir, prefix)
    savemaxQgo(list_Qgo, dir, prefix)
    savemaxer(list_er, dir, prefix)
    savesteps(list_step, dir, prefix)

main()
