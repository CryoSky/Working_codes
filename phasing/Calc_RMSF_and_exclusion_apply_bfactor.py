#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2018-Oct-22, latest modified on 2018-Nov-11
# Due to limitation of MDAnalysis, this python code must be run on Linux system and python3

import MDAnalysis as mda
import matplotlib
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSF

matplotlib.use('Agg')  # Avoid _tkinter.TclError: no display name and no $DISPLAY environment variable
import matplotlib.pyplot as plt
import numpy
import Bio.PDB as bpdb
from Bio.PDB.PDBParser import PDBParser
import os
import argparse
import math

'''
def exclusion(exclusion):
    if id in exclusion:
        return id
    return -1    # Avoid returning NoneType
'''


class Backboneselect(bpdb.Select):
    def __init__(self):
        super().__init__()  # Inherit attributes from parents

    def accept_atom(self, atom):
        if atom.get_name() == 'CA' or atom.get_name() == 'CB' or atom.get_name() == 'O' or atom.get_name() == 'HB':
            return 1
        else:
            return 0


class RMSFexclusion(bpdb.Select):
    def __init__(self, exclusion):
        super().__init__()  # Inherit attributes from parents
        self.exclusion = exclusion

    def accept_residue(self, residue):
        flag = 1
        for exclusion_id in self.exclusion:
            if residue.get_id()[1] == exclusion_id:
                flag = 0
        return flag


def plot(resnums, rmsf):
    plt.plot(resnums, rmsf)
    plt.savefig("rmsf_per_residue.png")


def topology_backbone(topology_prefix):
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure('topology', '%s.pdb' % topology_prefix)

    output = bpdb.PDBIO()
    output.set_structure(structure)
    output_filename = str(topology_prefix) + '_backbone.pdb'
    output.save('%s' % output_filename, Backboneselect())
    print("The backbone of selected topology file has been saved as %s" % output_filename)
    return output_filename


def RMSFcalculation(topology, trajectory, threshold):
    # https://www.mdanalysis.org/docs/documentation_pages/analysis/rms.html#MDAnalysis.analysis.rms.RMSF
    u = mda.Universe("%s" % topology, "%s" % trajectory, in_memory=True)  # Initialization in MDAnalysis
    ref = mda.Universe("%s" % topology)  # The backbone file has 3 atoms per residue
    prealigner = align.AlignTraj(u, ref, select="protein and name CA",
                                 in_memory=True).run()  # Fit to the best frame to get a better average structure
    protein = u.select_atoms("protein")
    reference_coordinates = u.trajectory.timeseries(asel=protein).mean(axis=1)
    print ("The atomic coordinates are\n%s" % reference_coordinates)
    reference = mda.Merge(protein).load_new(reference_coordinates[:, None, :], order="afc")
    aligner = align.AlignTraj(u, reference, select="protein and name CA", in_memory=True).run()
    calphas = protein.select_atoms("name CA")
    rmsfer = RMSF(calphas, verbose=True).run()
    # rmsf_data = numpy.zeros(rmsfer.rmsf.size)
    resnums = calphas.resnums
    rmsf_data = rmsfer.rmsf
    plot(resnums, rmsf_data)
    selection_residue = []
    index = 1
    for rmsf_pointer in rmsf_data:
        if rmsf_pointer >= threshold:
            selection_residue.append(index)
        index += 1

    return selection_residue, rmsf_data


def pdbBfactor(topology_prefix, rmsf_data):
    '''Learned from https://github.com/harmslab/pdbtools/blob/master/pdbtools/bfactor.py
    Goes throungh pdb line by line. The b-factor of pdb file is replaced by RMSF converted
    value'''
    pdbfile = topology_prefix + '.pdb'
    output = []
    with open(pdbfile, 'r') as f1:
        for line in f1:
            if line[0:6] == "ATOM  ":
                resnum = line[23:26].strip()
                index = int(resnum) - 1
                bfactor = (8 * pow(math.pi, 2) / 3) * pow(rmsf_data[index], 2)  # B=[(8*PI**2)/3] * (RMSF)**2
                output.append("%s%6.2F%s" % (line[:60], bfactor, line[66:]))
            elif line[0:6] == "HETATM":
                output.append("%s%6s%s" % (line[:60], "NA", line[66:]))
            else:
                output.append(line)
    with open('%s_rmsf_contained.pdb' % topology_prefix, 'w+') as f2:
        for line in output:
            f2.write(line)
    print("The pdb file with rmsf value has been saved as %s_rmsf_contained.pdb" % topology_prefix)
    rmsf_topology_prefix = '%s_rmsf_contained' % topology_prefix
    return rmsf_topology_prefix
    # return output


def topology_rmsfexclusion(topology_prefix, selection_residue):
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure('topology', '%s.pdb' % topology_prefix)
    # residues = structure.get_residues()

    output = bpdb.PDBIO()
    output.set_structure(structure)
    output.save('%s_rmsf_exclusion.pdb' % topology_prefix, RMSFexclusion(selection_residue))


def main():
    parser = argparse.ArgumentParser(
        description="This script calculates the RMSF value of each residue based on C_alpha atom and exclude\
        the atoms based on specified threshold. The script must run under Linux system and Python3.")
    parser.add_argument("trajectory", help="The file name of trajectory, in AWSEM usually dump.lammpstrj", type=str)
    parser.add_argument("topology", help="The file name of topology for align, any frame of run.pdb is ok", type=str)
    parser.add_argument("threshold", help="The RMSF exclusion threshold", type=float)

    args = parser.parse_args()
    trajectory = args.trajectory
    topology = args.topology
    threshold = args.threshold

    if topology[-4:].lower() != ".pdb":
        topology = topology + ".pdb"
    topology_prefix = topology[:-4]

    topology_filename = topology_backbone(topology_prefix)
    (selection_residue, rmsf_data) = RMSFcalculation(topology_filename, trajectory, threshold)
    os.rename('rmsf_per_residue.png', "%s_rmsf_per_residue.png" % topology_prefix)
    print("The figure has been saved as %s_rmsf_per_residue.png" % topology_prefix)
    rmsf_topology_prefix = pdbBfactor(topology_prefix, rmsf_data)
    topology_rmsfexclusion(rmsf_topology_prefix, selection_residue)


if __name__ == '__main__':
    main()
