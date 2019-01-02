#!/usr/bin/python3
# -*- coding: utf-8 -*-
# Written by Shikai Jin, 2019-01-01, latest modified 2019-01-01, from Mingchen's code
# This file is used for checking the gro files indexed by frag.mem, are they contains all residues
# Example in Linux: python *this_file*.py frags.mem
import os
import sys

mem = sys.argv[1]

with open("tmp.mem", "w") as out:
    with open(mem, "r") as f:
        for i in range(4):  # Pass the first 4 lines
            line = next(f)
            out.write(line)
        for line in f:
            gro = line.split()[0]
            i = line.split()[2]  # Starting residue of fragment
            n = line.split()[3]  # How many residues loaded in the gro file
            delete = False

            with open(gro, "r") as one:
                next(one)
                next(one)  # Pass the first 2 lines
                all_residues = set()  # Initialize a set
                for atom in one:
                    residue = atom.split()[0]
                    # print(residue)
                    all_residues.add(int(residue))  # Save all residue indexes in gro file
                for test in range(int(i), int(i) + int(n)):  # Check whether are there some missing residues in gro file
                    if test not in all_residues:
                        print("ATTENTION", gro, i, n, "missing:", test)
                        delete = True
            if not delete:  # If delete=true then discard this gro file in output
                out.write(line)

os.system("mv %s %s_back" % (mem, mem))
os.system("mv tmp.mem %s" % (mem))
