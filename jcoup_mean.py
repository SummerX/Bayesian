#!/usr/bin/python
import argparse
import mdtraj as md
import numpy as np
import os,h5py
from progressbar import ProgressBar
import progressbar
from mdtraj import io


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

    parser.add_argument('-j',
    	'--jcoup',
    	help = 'all jcoupling file',
    	default = 'Data/J_coup_kalplus.dat'
    )

    parser.add_argument('-o',
        '--output',
        default = 'Data/jcoup_mean.txt',
        help ="output file for leader.txt default is Data/jcoup_mean.txt"
        )

    parser.add_argument('-omd',
        '--mdoutput',
        default = 'Data/md_mean.dat',
        help = "output file for md mean default is Data/md_mean.dat")

    parser.add_argument('-omdsd',
        '--mdsdoutput',
        default = 'Data/md_sd.dat',
        help = "output file for md sd default is Data/md_sd.dat")

    subparsers = parser.add_subparsers(help = 'types of input assignment')

    leader_parser = subparsers.add_parser('leader')
    leader_parser.set_defaults(which = 'leader')

    leader_parser.add_argument('-i',
        '--input',
        help = "input file for list of structure from each cluster",
        default = 'Data/leader.txt'
        )


    all_parser = subparsers.add_parser('all')
    all_parser.set_defaults(which = 'all')

    all_parser.add_argument('-i',
        '--input',
        help = "input file for list of structure from each cluster",
        default = 'Data/lfdp_assign.h5'
        )

    all_parser.add_argument('-p',
        '--percent',
        help = "conformation above percent to keep for Bayesian",
        default = 0.00,
        type = float
        )


    args = parser.parse_args()

    Jcoup = load_file(args.jcoup)

    if args.which == 'leader':
        
        output = open(args.output, 'w')

        with open(args.input,'r') as f:
        	for line in f:
        		cluster = np.array(map(int,line.split()))
        		J_mean = np.mean(Jcoup[cluster,:], axis = 0)
        		print len(cluster)
        		output.write('\t'.join(map(str,J_mean)))
        		output.write('\n')
        output.close()

    elif args.which == 'all':

        assignment = load_file(args.input)

        assignment = filter_assignment(assignment,args.percent)

        numstates = max(assignment) + 1

        md_output = open(args.mdoutput,'w')

        mdsd_output = open(args.mdsdoutput,'w')

        with open(args.output, 'w') as output:
            for i in range(numstates):
                cluster = np.where(assignment == i)[0]
                J_mean = np.mean(Jcoup[cluster,:], axis = 0)
                print len(cluster)
                output.write('\t'.join(map(str,J_mean)))
                output.write('\n')
                md_output.write(str(float(len(cluster))/len(assignment))+'\n')
                mdsd_output.write('0.20'+'\n')


if __name__ == '__main__':
    main()
