#! usr/bin/env python3
# Written by Shikai Jin on 2018-03-10, latest modified on 2018-03-12
# This file is used for -ing
# Example in Linux: Python ContactMap.py 1mba.pdb
# Only works on AWSEM genrated files due to 'T' on each line

import string
import numpy as np
import sys
import math

pdbfile = sys.argv[1]

with open(pdbfile, 'r') as fp:
    lines = fp.readlines()
    i = -1
    aanum = 0
    while True:  # find total amino acids number by residue number of last 'ATOM' line
        temp1 = lines[i].split()
        if temp1[0] == 'TER':
            aanum = temp1[3]
            break
        elif temp1[0] == 'ATOM':
            aanum = temp1[4]
            if aanum == 'T':
                aanum = temp1[5]
            break
        else:
            i -= 1
    aanum = int(aanum)
    dismatrix = np.zeros(shape=(aanum, aanum))

    x, y, z = 0, 0, 0
    d = []
    residue = 0
    for line in lines:
        line = line.strip()
        temp2 = line.split()

        if temp2[0] == 'ATOM' and temp2[2] == 'CA':
            if temp2[4] == 'T':
                x = float(temp2[6])
                y = float(temp2[7])
                z = float(temp2[8])
            else:
                x = float(temp2[5])
                y = float(temp2[6])
                z = float(temp2[7])
            residue += 1

            if residue == 1:
                var = []
                var.append(x) # append method, change the original list, return None and dont generate new list
                var.append(y)
                var.append(z)
                d = np.array(var, dtype=complex)
            else:
                var = []
                var.append(x)
                var.append(y)
                var.append(z)
                newline = np.array(var, dtype=complex)
                temp = np.row_stack((d, newline))
                d = temp
    print(d)

for i in range(1,aanum+1): # generate number from 1 to aanum
    for j in range(1,aanum+1):
        distance = math.sqrt((d[j-1,0]-d[i-1,0])**2+(d[j-1,1]-d[i-1,1])**2+(d[j-1,2]-d[i-1,2])**2)
        print(distance)
        dismatrix[j-1,i-1]=distance
np.savetxt('contact_matrix.txt',dismatrix)
