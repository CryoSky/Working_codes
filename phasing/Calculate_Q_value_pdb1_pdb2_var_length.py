#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2018-Nov-17, latest modified on 2018-Nov-17
# A big extension of calculate Q script
# Example in Linux: python /mnt/e/linux/script/Calculate_Q_value_pdb1_pdb2_var_length.py T1003/3D-JIGSAW_SL1_TS1.pdb T1003_reference_6hrh_A.pdb 0 --start1 30 --end1 465


import argparse
import os
import sys
from Bio.PDB.PDBParser import PDBParser
from math import *


def vector(p1, p2):
    return [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]]


def vabs(a):
    return sqrt(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2))


cutoff = 9.5


# For the equation of calculate the Q value read dx.doi.org/10.1021/jp212541y

def compute_Q_pdb1_pdb2(ca_atoms_pdb1, ca_atoms_pdb2, q_type, sigma_sq):
    Q_value = 0
    count_a = 0
    ca_length = len(ca_atoms_pdb2)
    if q_type == 1:
        minimum_separation = 4
    elif q_type == 0:
        minimum_separation = 3
    else:
        print("The Q value type %s should be either 1 or 0" % q_type)
        exit()

    if len(ca_atoms_pdb1) != len(ca_atoms_pdb2):
        print("The number of C alpha atoms in pdb1 file %s doesn't match the number in pdb2 file %s" % (
            len(ca_atoms_pdb1), len(ca_atoms_pdb2)))
        exit()

    for ia in range(ca_length):
        for ja in range(ia + minimum_separation, ca_length):
            if (ia + 1) in ca_atoms_pdb1 and (ja + 1) in ca_atoms_pdb1:
                rij_N = vabs(vector(ca_atoms_pdb1[ia + 1], ca_atoms_pdb1[ja + 1]))
                if q_type == 1 and rij_N >= cutoff:  # What's the function of this line?
                    continue
                rij = vabs(vector(ca_atoms_pdb2[ia + 1], ca_atoms_pdb2[ja + 1]))
                dr = rij - rij_N
                Q_value = Q_value + exp(-dr * dr / (2 * sigma_sq[ja - ia]))
                count_a = count_a + 1
    Q_value = Q_value / count_a
    return Q_value


def calc_sigma_sq(ca_atoms_pdb):
    sigma = []
    sigma_sq = []
    sigma_exp = 0.15

    for i in range(0, len(ca_atoms_pdb) + 1):
        sigma.append((1 + i) ** sigma_exp)
        sigma_sq.append(sigma[-1] * sigma[-1])

    return sigma_sq


def pdb_load(pdb_file, start, end):
    p = PDBParser(PERMISSIVE=1)  # Useful when find errors in PDB
    s = p.get_structure('pdb', pdb_file)
    chains = s[0].get_list()  # Compartible for multiple chains, but only for first model
    chain_id = 0
    ca_atoms_pdb = {}  # A new dictionary to record the C_alpha atom coordinate in pdb file
    pdb_chain_id = []

    for chain in chains:
        chain_id = chain_id + 1
        first_residue = chain.get_unpacked_list()[0]
        last_residue = chain.get_unpacked_list()[-1]
        sequence_id_flag = 0
        if end is None:  # Default end point is last residue
            end = int(last_residue.get_id()[1])
            # print("End value is %s" %end)
            # print(chain[353]['CA'].get_coord())
        if start is None:  # Some fxxking pdbs start from -7
            start = int(first_residue.get_id()[1])
        for res in chain:
            is_regular_res = res.has_id('CA') and res.has_id('O') or (
            res.get_id()[1] == last_residue or res.has_id('CA'))
            hetero_flag = res.get_id()[
                0]  # Get a list for each residue, include hetero flag, sequence identifier and insertion code
            # https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
            if (hetero_flag == ' ' or hetero_flag == 'H_MSE' or hetero_flag == 'H_M3L' or hetero_flag == 'H_CAS') \
                    and is_regular_res:
                # The residue_id indicates that there is not a hetero_residue or water ('W')
                sequence_id = res.get_id()[1]
                # print(sequence_id)
                if sequence_id_flag == 0:
                    last_sequence_id = sequence_id
                else:
                    if last_sequence_id != int(sequence_id) - 1:
                        print ("WARNING: the residues between %s and %s are lost" % (last_sequence_id, sequence_id))
                        # Some fxxking pdbs lost residues halfway
                    last_sequence_id = sequence_id
                    # print(res.get_id())
            if sequence_id >= start and sequence_id <= end:
                ca_atoms_pdb[sequence_id] = res['CA'].get_coord()
                # print(ca_atoms_pdb.keys())
                # ca_atoms_pdb.append(res['CA'].get_coord())
                pdb_chain_id.append(chain_id)
            sequence_id_flag = sequence_id_flag + 1
    #print(ca_atoms_pdb.keys())
    return ca_atoms_pdb


# def check_residue(ca_atoms_pdb1, ca_atoms_pdb2):




def main():
    parser = argparse.ArgumentParser(
        description="This script calculates Q value for all chains in two pdbs.")
    parser.add_argument("PDB_filename1", help="The file name of input pdb1", type=str)
    parser.add_argument("--start1", help="The start residue of protein1", type=int)
    parser.add_argument("--end1", help="The end residue of protein1", type=int)
    parser.add_argument("PDB_filename2", help="The file name of input pdb2", type=str)
    parser.add_argument("--start2", help="The start residue of protein2", type=int)
    parser.add_argument("--end2", help="The end residue of protein2", type=int)
    # parser.add_argument("output_filename", help="The file name of output file", type=str)
    parser.add_argument("q_type", help="The Q value type, 0 for Q_wolynes, 1 for Q_onuchic", type=int)
    args = parser.parse_args()
    pdb1_file = args.PDB_filename1
    pdb2_file = args.PDB_filename2
    q_type = args.q_type
    start1 = args.start1
    end1 = args.end1
    start2 = args.start2
    end2 = args.end2
    print(end2)

    if pdb1_file[-4:].lower() != ".pdb" or pdb2_file[-4:].lower() != ".pdb":
        print("It must be a pdb file.")
        exit()

    ca_atoms_pdb1 = pdb_load(pdb1_file, start1, end1)
    ca_atoms_pdb2 = pdb_load(pdb2_file, start2, end2)
    #print(ca_atoms_pdb2)

    sorted_keys1 = sorted(ca_atoms_pdb1.keys())
    sorted_keys2 = sorted(ca_atoms_pdb2.keys())

    range_pdb1 = int(list(sorted_keys1)[-1]) - int(list(sorted_keys1)[0])
    range_pdb2 = int(list(sorted_keys2)[-1]) - int(list(sorted_keys2)[0])

    print("first range %s" %range_pdb1)
    print("second range %s"%range_pdb2)

    if len(ca_atoms_pdb1) != len(ca_atoms_pdb2):
        if range_pdb1 == range_pdb2:
            print("The two pdb have same residue range but different length.\
             So there are some residues lost in the halfway. Now checking...")
            difference = sorted_keys1[0] - sorted_keys2[0]
            for i in sorted_keys1:
                if i - difference in sorted_keys2:
                    continue
                else:
                    print(
                    "The residue %s in protein1 has no corresponding residue in protein2, the value of it %s has been removed" % (
                    i, ca_atoms_pdb1[i]))
                    ca_atoms_pdb1.pop(i)
            for j in sorted_keys2:
                if j + difference in sorted_keys1:
                    continue
                else:
                    print(
                    "The residue %s in protein2 has no corresponding residue in protein1, the value of it %s has been removed" % (
                    j, ca_atoms_pdb2[j]))
                    ca_atoms_pdb2.pop(j)
            if len(ca_atoms_pdb1) != len(ca_atoms_pdb2):
                print("No idea.")
                exit()
        else:
            print ("Error: Two PDB structures have different lengths and ranges!")
            print (ca_atoms_pdb1.keys())
            print (len(ca_atoms_pdb1))
            print (ca_atoms_pdb2.keys())
            print (len(ca_atoms_pdb2))
            exit()

    new_ca_atoms_pdb1 = {i + 1: ca_atoms_pdb1[k] for i, k in enumerate(sorted(ca_atoms_pdb1.keys()))}
    #print(new_ca_atoms_pdb1.keys())  # Check the number whether start from 1 and match the length
    new_ca_atoms_pdb2 = {i + 1: ca_atoms_pdb2[k] for i, k in enumerate(sorted(ca_atoms_pdb2.keys()))}
    # Change the keys of ca_atoms_pdb dict into natural increment\
    # Since we trim the model so the real residue number may not work, if we have same length just reorder from 0
    # https://stackoverflow.com/questions/39126272/reset-new-keys-to-a-dictionary

    sigma_sq = calc_sigma_sq(new_ca_atoms_pdb1)

    if len(ca_atoms_pdb1) > 0:
        q = compute_Q_pdb1_pdb2(new_ca_atoms_pdb1, new_ca_atoms_pdb2, q_type, sigma_sq)  # For last frame
        # with open(output_file, 'a+') as out:
        #    out.write(str(round(q, 3)))
        #    out.write(' ')
        #    out.write('\n')
        print(str(round(q, 3)), end="")  # Cancel line break

if __name__ == '__main__':
    main()
