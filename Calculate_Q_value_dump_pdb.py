#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2018-Nov-7, latest modified on 2018-Nov-7
# Modified from Mingchen's python2 version CalcQValue_structure.py
# Example in Linux: Python Calculate_Q_value_dump_pdb.py test.pdb dump.lammpstrj q_output 0

import argparse
from math import *
import os
import sys
from Bio.PDB.PDBParser import PDBParser


def vector(p1, p2):
    return [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]]


def vabs(a):
    return sqrt(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2))


atom_type = {'1': 'C', '2': 'N', '3': 'O', '4': 'C', '5': 'H', '6': 'C'}
atom_desc = {'1': 'C-Alpha', '2': 'N', '3': 'O',
             '4': 'C-Beta', '5': 'H-Beta', '6': 'C-Prime'}
PDB_type = {'1': 'CA', '2': 'N', '3': 'O', '4': 'CB', '5': 'HB', '6': 'C'}
cutoff = 9.5


# For the equation of calculate the Q value read dx.doi.org/10.1021/jp212541y

def compute_Q_dump_pdb(ca_atoms_pdb, ca_atoms_dump, q_type, sigma_sq):
    Q_value = 0
    count_a = 0
    ca_length = len(ca_atoms_dump)
    if q_type == 1:
        minimum_separation = 4
    elif q_type == 0:
        minimum_separation = 3
    else:
        print("The Q value type %s should be either 1 or 0" % q_type)
        exit()

    if len(ca_atoms_dump) != len(ca_atoms_pdb):
        print("The number of C alpha atoms in dump file %s doesn't match the number in pdb file %s" %(len(ca_atoms_dump), len(ca_atoms_pdb)))
        exit()

    for ia in range(ca_length):
        for ja in range(ia + minimum_separation, ca_length):
            if (ia + 1) in ca_atoms_pdb and (ja + 1) in ca_atoms_pdb:
                rij_N = vabs(vector(ca_atoms_pdb[ia + 1], ca_atoms_pdb[ja + 1]))
                if q_type == 1 and rij_N >= cutoff: # What's the function of this line?
                    continue
                rij = vabs(vector(ca_atoms_dump[ia], ca_atoms_dump[ja]))
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

def pdb_load(pdb_file):
    p = PDBParser(PERMISSIVE=1)  # Useful when find errors in PDB
    s = p.get_structure('pdb', pdb_file)
    chains = s[0].get_list()  # Compartible for multiple chains, but only for first model
    chain_id = 0
    ca_atoms_pdb = {}  # A new dictionary to record the C_alpha atom coordinate in pdb file
    pdb_chain_id = []


    for chain in chains:
        chain_id = chain_id + 1
        for res in chain:
            is_regular_res = res.has_id('CA') and res.has_id('O')
            hetero_flag = res.get_id()[0]  # Get a list for each residue, include hetero flag, sequence identifier and insertion code
            # https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
            if (hetero_flag == ' ' or hetero_flag == 'H_MSE' or hetero_flag == 'H_M3L' or hetero_flag == 'H_CAS') \
                    and is_regular_res:
                # The residue_id indicates that there is not a hetero_residue or water ('W')
                sequence_id = res.get_id()[1]
                ca_atoms_pdb[sequence_id] = res['CA'].get_coord()
                # ca_atoms_pdb.append(res['CA'].get_coord())
                pdb_chain_id.append(chain_id)
    return ca_atoms_pdb


def lammps_load_and_calc(lammpsdump_file, ca_atoms_pdb, q_type, sigma_sq, output_file):
    ca_atoms_dump = []
    box_type = []
    box_boundary = []
    with open(lammpsdump_file, 'r') as dump_file:
        for dump_line in dump_file:
            dump_line = dump_line.strip()
            if dump_line[:5] == "ITEM:":
                item = dump_line[6:]  # Then check next line what happened
            else:
                if item == "TIMESTEP":
                    if len(ca_atoms_dump) > 0:
                        q = compute_Q_dump_pdb(ca_atoms_pdb, ca_atoms_dump, q_type, sigma_sq) # Calculate Q online, one frame by one frame
                        with open(output_file, 'a+') as out:
                            out.write(str(round(q, 3)))
                            out.write(' ')
                            out.write('\n')
                    # After calculation for each frame, refresh the three lists
                        ca_atoms_dump = []
                        box_type = []
                        box_boundary = []

                elif item == "NUMBER OF ATOMS":
                    continue
                elif item[:10] == "BOX BOUNDS":
                    box_type.append(dump_line)
                    dump_line = dump_line.split()
                    box_boundary.append([float(dump_line[0]), float(dump_line[1])])
                elif item[:5] == "ATOMS":
                    dump_line = dump_line.split()
                    i_atom = dump_line[0]
                    x = float(dump_line[2])
                    y = float(dump_line[3])
                    z = float(dump_line[4])
                    x_position = (box_boundary[0][1] - box_boundary[0][0]) * x + box_boundary[0][0]
                    y_position = (box_boundary[1][1] - box_boundary[1][0]) * y + box_boundary[1][0]
                    z_position = (box_boundary[2][1] - box_boundary[2][0]) * z + box_boundary[2][0]
                    desc = atom_desc[dump_line[1]]
                    if desc == 'C-Alpha':
                        residue_id = int((int(i_atom) + 2) / 3)
                        atom = [x_position, y_position, z_position, residue_id] # record C_alpha atom coordinate
                        ca_atoms_dump.append(atom)

    q = compute_Q_dump_pdb(ca_atoms_pdb, ca_atoms_dump, q_type, sigma_sq) # For last frame
    with open(output_file, 'a+') as out:
        out.write(str(round(q, 3)))
        out.write(' ')
        out.write('\n')

def main():

    parser = argparse.ArgumentParser(
        description="This script calculates Q value for each chain in a pdb with every frame of a dump.lammpstrj from AWSEM")
    parser.add_argument("PDB_filename", help="The file name of input pdb", type=str)
    parser.add_argument("dump_filename", help="The file name of dump file, usually dump.lammpstrj", type=str)
    parser.add_argument("output_filename", help="The file name of output file", type=str)
    parser.add_argument("q_type", help="The Q value type, 0 for Q_wolynes, 1 for Q_onuchic", type=int, default=1)
    args = parser.parse_args()
    pdb_file = args.PDB_filename
    dump_file = args.dump_filename
    output_file = args.output_filename
    q_type = args.q_type

    if pdb_file[-4:].lower() != ".pdb":
        print("It must be a pdb file.")
        exit()

    ca_atoms_pdb = pdb_load(pdb_file)
    sigma_sq = calc_sigma_sq(ca_atoms_pdb)
    lammps_load_and_calc(dump_file, ca_atoms_pdb, q_type, sigma_sq, output_file)

if __name__ == '__main__':
    main()
