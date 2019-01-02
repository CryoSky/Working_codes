#!/usr/bin/python3
# -*- coding: utf-8 -*-
# Written by Xun Chen, latest modified 2019-01-01
# This file is used for generating extended structure by given a sequence file using Pymol

import os
import string
import re
import sys

squence_file = sys.argv[1] + '.seq'
pml_file = sys.argv[2] + '.pml'
pdb_file = sys.argv[3] + '.pdb'


# title of fasta >1QYS:A|PDBID|CHAIN|SEQUENCE

with open(squence_file,'r') as fopen:
    lines = fopen.readlines()
    for line1 in lines:
       line = line1.rstrip('\n')
       break 
with open (pml_file,'w') as fwrite:
    data = ''
    data += 'sequence = \"' + line + '\" \n'
    data += 'for aa in sequence: cmd._alt(string.lower(aa))\n'
    data += 'alter (all),resi=str(int(resi)-1)\n' 
    data += 'save ' + pdb_file + '\n'
    data += 'quit \n'
    fwrite.writelines(data)
os.system("pymol %s"%(pml_file))
     

