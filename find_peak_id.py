#!/usr/bin/python
import argparse
# import metrics
import mdtraj as md
import numpy as np
import os,h5py
# import matplotlib.pyplot as plt
from progressbar import ProgressBar
import progressbar
from mdtraj import io

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i',
        '--input',
        help = "input assignment directory from lfdp"
        )

    parser.add_argument('-o',
        '--output',
        default = 'Data/leader.txt',
        help ="output file for leader.txt default is Data/leader.txt"
        )

    parser.add_argument('-p',
        '--percent',
        default = 0.0,
        type = float,
        help = "Percent of neighbours(default = 0)"
        )

    args = parser.parse_args()


    input_file = os.path.join(os.getcwd(),args.input,'lf_leader.h5')
    lf_leader = io.loadh(input_file)['arr_0']
    input_file = os.path.join(os.getcwd(),args.input,'lf_assign.h5')
    lf_assign = io.loadh(input_file)['arr_0']

    input_file = os.path.join(os.getcwd(),args.input,'dp_peaks.h5')
    dp_peaks = io.loadh(input_file)['arr_0']
    input_file = os.path.join(os.getcwd(),args.input,'dp_assign.h5')
    dp_assign = io.loadh(input_file)['arr_0']

    input_file = os.path.join(os.getcwd(),args.input,'lfdp_peaks.h5')
    lfdp_leader = io.loadh(input_file)['arr_0']
    input_file = os.path.join(os.getcwd(),args.input,'lfdp_assign.h5')
    lfdp_assign = io.loadh(input_file)['arr_0']


    Ass = lfdp_assign
    Leader = lfdp_leader
    size = Ass.size

    NumStates = max(Ass) + 1
    state = 0

    for i in range(NumStates):
        cluster = np.where(Ass == i)[0]
        if len(cluster) > size * args.percent :
            Ass[cluster] = state
            Leader[state] = Leader[i]
            state += 1
        else:
            Ass[cluster] = -1
            Leader[i] = -1


    NumStates = max(Ass) + 1

    if not os.path.isdir(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))

    with open(args.output,'w') as f:
        for i in range(NumStates):
            temp = np.where(lf_assign == lf_assign[Leader[i]])[0].tolist()
            f.write('\t'.join(map(str,temp)))
            f.write('\n')

if __name__ == '__main__':
    main()
