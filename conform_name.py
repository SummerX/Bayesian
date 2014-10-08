#!/usr/bin/python
import argparse
import mdtraj as md
import numpy as np
import os,h5py
from progressbar import ProgressBar
import progressbar
from mdtraj import io
from conformation_detect import *
import scipy.stats.mstats

def load_file(jcoup_file):
    fileName, fileExtension = os.path.splitext(jcoup_file)
    if fileExtension == '.h5':
        return io.loadh(jcoup_file)['arr_0']
    else:
        return np.loadtxt(jcoup_file, dtype = float)

def filter_assignment(assignment,percent):
    numstates = max(assignment) + 1   
    state = 0
    size = len(assignment)

    for i in range(numstates):
        cluster = np.where(assignment == i)[0]
        if len(cluster) > size * percent:
            assignment[cluster] = state
            state += 1
        else:
            assignment[cluster] = -1
    return assignment


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('assignment',
        help = 'Path to an assignment file.')

    parser.add_argument('-o',
        '--output',
        default = 'Data/conform.dat',
        help = "output file for conformation file default is Data/conform.dat")
    
    parser.add_argument('-t',
        '--trajectory',
        default='trajectory',
        help ='Path to the trajectory directory or project file (default=trajectory)'
        )

    parser.add_argument('-p',
        '--percent',
        help = "conformation above percent to keep for Bayesian",
        default = 0.00,
        type = float
        )

    parser.add_argument('-d',
        '--angle',
        default = 'conform_angle.dat',
        help = "angle used to calculate Jcoupling default is conform_angle.dat (phi first, psi second)")

    args = parser.parse_args()

    traj = md.load(os.path.join(os.getcwd(),args.trajectory,'traj.h5'))
    
    phi = md.compute_phi(traj)[1] * 180.0 / np.pi
    psi = md.compute_psi(traj)[1] * 180.0 / np.pi

    if os.path.isfile(args.angle):
        angle = np.loadtxt(args.angle,dtype = int,unpack=True).reshape((2,-1))
        phi = np.hstack((phi[:,angle[0,:]],phi[:,[0]]))
        psi = np.hstack((psi[:,[0]],psi[:,angle[1,:]]))


    assignment = load_file(args.assignment)

    assignment = filter_assignment(assignment,args.percent)

    numstates = max(assignment) + 1

    with open(args.output, 'w') as output:
        for i in range(numstates):
            cluster = np.where(assignment == i)[0]
            conform = ''
            for i in range(len(phi[0,:])-1):
                conform = conform +'-'+detect_conf_print(scipy.stats.mstats.mode(map(int,phi[cluster,i]))[0],scipy.stats.mstats.mode(map(int,psi[cluster,i+1]))[0])
            output.write(conform[1:]+'\n')



if __name__ == '__main__':
    main()
