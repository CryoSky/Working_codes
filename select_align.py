#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2018-Jun-23, latest modified on 2018-Jun-28
# This is a revised Python3 version of old script (Weihua's version)
# This script is used for getting the aligned regions from the sequence alignment algorithm.


import Bio.PDB as bpdb
import argparse


def read_alignment(alignment):
    with open(alignment, 'r') as f:
        strs = f.readline().split()
        target = strs[1]
        targ_start = int(strs[0])
        strs = f.readline().split()
        template = strs[1]
        temp_start = int(strs[0])

        print('The target length is: ' + str(len(target)))
        print('The template length is: ' + str(len(template)))

        if len(target) != len(template):
            print('ERROR: The target length is not equal to template length!')
            exit()

        targ_ind = targ_start - 1
        temp_ind = temp_start - 1
        targ_seq = []
        temp_seq = []

        for i in range(len(target)):
            if template[i] != '-':
                temp_ind = temp_ind + 1
            if target[i] != '-':
                targ_ind = targ_ind + 1
            if target[i] != '-' and template[i] != '-':
                targ_seq.append(targ_ind)
                temp_seq.append(temp_ind)

        print('The target sequence is ' + str(targ_seq))
        print('The template sequence is ' + str(temp_seq))

        return targ_seq, temp_seq

def find(resid, seqind):
    for i in range(len(seqind)):
        if resid == seqind[i]:
            return (i + 1)
    return -1    # Avoid returning NoneType


class ResSelect(bpdb.Select):
    def __init__(self, temp_seq, chain_id):
        super().__init__() # Inherit attributes from parents
        self.temp_seq = temp_seq
        self.chain_id = chain_id

    def accept_residue(self, res):
        ind = find(res.id[1], self.temp_seq)
        if ind > 0 and res.parent.id == self.chain_id:
            return True
        else:
            return False

def main():
    parser = argparse.ArgumentParser(
        description="This script gets the aligned regions from the sequence alignment algorithm")
    parser.add_argument("input", help="The file name of input pdb", type=str)
    parser.add_argument("output", help="The file name of output pdb", type=str)
    parser.add_argument("align", help="The alignment file name", type=str)
    parser.add_argument("chain", help="The chain ID", type=str)
    args = parser.parse_args()
    input = args.input
    output = args.output
    alignment = args.align
    chain_id = args.chain

    if input[-4:].lower() != ".pdb":
        input = input + ".pdb"
    if output[-4:].lower() != ".pdb":
        output = output + ".pdb"

    (targ_seq, temp_seq) = read_alignment(alignment)

    s = bpdb.PDBParser().get_structure('temp', input)
    io1 = bpdb.PDBIO()
    io1.set_structure(s)
    io1.save('temp1.pdb', ResSelect(temp_seq, chain_id))

    s = bpdb.PDBParser().get_structure('temp', 'temp1.pdb')
    model = s[0]
    chain = model[chain_id]

    for residue in chain:
        resind = targ_seq[find(residue.id[1], temp_seq) - 1]
        residue.id = (' ', resind, ' ')

    io2 = bpdb.PDBIO()
    io2.set_structure(s)
    io2.save(output)

if __name__ == '__main__':
    main()
