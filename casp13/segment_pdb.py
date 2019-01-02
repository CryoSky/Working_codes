#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2018-May-24, latest modified on 2018-Jun-13

import prody
import sys

def get_chain(pdbfile, chain_id):
    prefix = pdbfile.split('.')[0]
    
    structure = prody.parsePDB(pdbfile)
    hv = structure.getHierView()
    temp = hv[chain_id]
    protein = temp.select('protein') # Select protein part of pdb
    prody.writePDB('%s_%s.pdb' % (prefix, chain_id), protein)

def main():
    pdbfile = sys.argv[1]
    chain_id = sys.argv[2]
    if pdbfile[-4:].lower()!=".pdb":
        pdbfile = pdbfile + ".pdb"
    if len(sys.argv)!= 3: # 2 arguments and #0
        print ("\n" + sys.argv[0] + " PDB_file Chain_ID\n")
	exit()
    get_chain(pdbfile, chain_id)


if __name__== '__main__':
    main()
