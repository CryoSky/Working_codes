#! usr/bin/env python3
# Written by Shikai Jin, 2018-02-08, latest modified 2018-02-15
# This file is used for extracting last 20 frames of energy.log, calculate the average and standard deviation
# Title is Chain, Chi, RAMA, DSSP, P_AP, Water, Burial, Helix, AMH_GO, Frag_Mem, Q_GO, V_Total
# Example in Linux: python batch_extract_latest_energy.py /scratch/sj52/Rotation3/RAMA_Tempfree/1mba 1mba-3pt8-HO+ER_1.5

import os
import string
import re
import numpy as np
import sys


def extractenergy(energy):  # extract all energy value from energy.log file
    last_energy = np.zeros(shape=(20,12))
    with open (energy, 'rb+') as fp: # you MUST add b in open because all other modes only allow pointer starts from the beginning
        offset = -5000    # set offset value for file pointer, for details see this page: http://blog.csdn.net/binchasing/article/details/51067923
        while True:
            fp.seek(offset, 2)  # seek(off, 2)means file pointer, from the end (2) go forward 5000 characters(offset)
            lines = fp.readlines()  # load all lines within the file pointer range
            if len(lines) >= 20:  # judge whether it has 20 lines
                j = 0 # pointer for the matrix line
                for i in range(-20,0):  # use range(-20,0) to avoid losing the last line
                    temp1 = lines[i].decode().split()  # split current lines by space ' ',remember to add () after split
                    if temp1[0] != 'Step' and temp1[0] != '1000':
                        temp2 = list(map(float, temp1)) # use map(float, ) function to convert all elements in temp1 to float
                        #print(temp2) 
                        temp3 = []
                        for k in (1,3,4,6,7,8,9,10,11,12,17,18): # k means useful energy term in each line
                            temp3.append(temp2[k])
                        print (temp3)
                        last_energy[j] = temp3
                        j += 1
                break
            else:
                offset *= 2
        #print(last_energy)
        return last_energy

def calculate(total_list):
    aver = np.mean(total_list, axis=0) # axis = 0 means for each column
    aver = aver.astype('U32') # I don't know why but you must change the dtype to U32 or you will get error like "TypeError: ufunc 'add' did not contain a loop with signature matching types dtype('<U32') dtype('<U32') dtype('<U32')"
    deviation = np.std(total_list, axis=0)
    deviation = deviation.astype('U32')
    return (aver, deviation)

def export(aver,deviation,dir,prefix):
    with open('%s/data_analysis/%s_energyaverage.txt' % (dir, prefix), 'w+') as fp:
        for i in aver:
            fp.write(i + '\n')
        fp.close()
    with open('%s/data_analysis/%s_energystd.txt' % (dir, prefix), 'w+') as fp:
        for j in deviation:
            fp.write(j + '\n')
        fp.close()


def main():
    dir = sys.argv[1]
    prefix = sys.argv[2]
    total_list = np.array([])
    # for rama in [0.25, 0.5, 1.0, 1.5, 2.0]
    for i in range(1,21):
        branch = dir + '/' + prefix + '_' + str(i)
        logfile = branch + '/energy.log'
        print (i)
        one_list = extractenergy(logfile)
        if i == 1:
           total_list = one_list
        else:
           np.concatenate((total_list, one_list)) # https://docs.scipy.org/doc/numpy-dev/reference/generated/numpy.concatenate.html you should add two parentheses
    (aver, deviation) = calculate(total_list)
    export(aver, deviation, dir, prefix)

main()
