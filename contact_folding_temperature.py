### Copyright By Mingchen Chen, Oct/31st/2017, revised by Shikai Jin
### Wolynes Research Group

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas
import prody
from matplotlib.pyplot import cm as cm


def trajectory_contact(trajectory, n_frame, skip):
    distances = []
    for snapshot in trajectory.getCoordsets()[np.arange(0, n_frame, skip)]:
        snap = snapshot[trajectory.getNames() == 'CA']
        d = []
        name = []
        i = 0
        for x1, y1, z1 in snap:
            i = i + 1
            j = 0
            for x2, y2, z2 in snap:
                j = j + 1
                if j > i + 1:
                    if np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2) <= 10:
                        d += [1]
                    else:
                        d += [0]
                    name += [(i, j)]
        distances += [pandas.Series(d, index=name)]
    Distances = pandas.concat(distances, axis=1)  # contact information for each selected frame
    return Distances


def calcTF(array, start, end):
    alen = len(array)  # import the one specific contact for all frames
    accuracy = np.zeros((1, alen))
    for i in range(alen):
        for j in range(alen):
            if j < i and array[j] == 0:
                accuracy[0, i] = accuracy[0, i] + 1 # Add 1 if the frame doesn't form the contact before the index
            if j > i and array[j] == 1:
                accuracy[0, i] = accuracy[0, i] + 1 # Add 1 if the frame forms the contact after the index
    pos = np.argmax(accuracy)  # returns the indices of the maximum values along an axis.
    return -pos * 1.0 / alen * (start - end) + start


def calc_temperature(trjnative, matrix, contactmp, Distances, Tstart, Tend):
    for snapshot in trjnative.getCoordsets()[np.arange(0, 1, 1)]:
        snap = snapshot[trjnative.getNames() == 'CA']
        i = 0
        count = 0
        for x1, y1, z1 in snap:
            i = i + 1
            j = 0
            for x2, y2, z2 in snap:
                j = j + 1
                if j > i + 1:
                    count = count + 1
                    if np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2) <= 10:
                        matrix[j - 1, i - 1] = calcTF(Distances.values[count - 1], Tstart, Tend)
                        # print Distances.values[count-1]
                        matrix[i - 1, j - 1] = matrix[j - 1, i - 1]
                        # print matrix[i-1,j-1], i, j
                        contactmp[j - 1, i - 1] = 1
    return matrix


def main():
    parser = argparse.ArgumentParser(
        description="This script calculates the folding temperature of contact pairs and visualize on contact map."
                    "Example Input: Python ContactFoldingTemperature.py n_frame seqlen input_trj.pdb native.pdb Tstart Tend SkipRate")
    parser.add_argument("n_frame", help="The number of frames in trajectory", type=int)
    parser.add_argument("seqlen", help="The sequence length", type=int)
    parser.add_argument("trajectory_file", help="The file name of input pdb", type=str)
    parser.add_argument("native_file", help="The file name of input pdb", type=str)
    parser.add_argument("Tstart", help="The start temperature of trajectory", type=float)
    parser.add_argument("Tend", help="The ending temperature of trajectory", type=float)
    parser.add_argument("skip", help="The skip steps", type=int)
    parser.add_argument("title", help="The file name of output", type=str)
    args = parser.parse_args()
    n_frame = args.n_frame
    seqlen = args.seqlen
    trajectory_file = args.trajectory_file
    native_file = args.native_file
    Tstart = args.Tstart
    Tend = args.Tend
    skip = args.skip
    title = args.title

    trajectory = prody.parsePDB(trajectory_file)
    trjnative = prody.parsePDB(native_file)

    Distances = trajectory_contact(trajectory, n_frame, skip)

    matrix = np.zeros((seqlen, seqlen))
    contactmp = np.zeros((seqlen, seqlen))

    new_matrix = calc_temperature(trjnative, matrix, contactmp, Distances, Tstart, Tend)

    plt.figure(figsize=(seqlen / 10, seqlen / 10))
    cmap = cm.get_cmap('jet', 30)
    cmap.set_under(color='white')
    cax = plt.imshow(new_matrix, interpolation="nearest", cmap=cmap, vmin=Tend, vmax=Tstart)
    plt.colorbar()
    plt.xticks(fontsize=20, rotation=90)
    plt.yticks(fontsize=20)
    # plt.show(cax)
    plt.savefig('%s_contact_TF.png' % title)


if __name__ == '__main__':
    main()
