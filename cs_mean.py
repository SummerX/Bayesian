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

def filter_cs(CS,atomindices):
    if not os.path.isfile(atomindices):
        return CS
    else:
        CS_new = np.empty(shape = (len(CS),0))
        with open(atomindices,'r') as f:
            for line in f:
                line = map(int,line.split())
                # print np.reshape(np.mean(CS[:,line],axis = 1),(-1,1)).shape
                CS_new = np.hstack((CS_new,np.reshape(np.mean(CS[:,line],axis = 1),(-1,1))))
    return CS_new

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-cs',
    	'--cs',
    	help = 'all jcoupling file',
    	default = 'Data/cs.dat'
    )

    parser.add_argument('-o',
        '--output',
        default = 'Data/cs_mean.txt',
        help ="output file for leader.txt default is Data/jcoup_mean.txt"
        )

    parser.add_argument('-a',
        '--atomindices',
        default = 'cs_atomindices.dat',
        help = 'atomindices file for calculate chemical shifts default is cs_atomindices.dat'
        )

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
        default = 'Data/lfdp_assignment.h5'
        )

    all_parser.add_argument('-p',
        '--percent',
        help = "conformation above percent to keep for Bayesian",
        default = 0.00,
        type = float
        )


    args = parser.parse_args()

    CS = load_file(args.cs)

    CS = filter_cs(CS,args.atomindices)




    if args.which == 'leader':
        
        output = open(args.output, 'w')

        with open(args.input,'r') as f:
        	for line in f:
        		cluster = np.array(map(int,line.split()))
        		CS_mean = np.mean(CS[cluster,:], axis = 0)
        		print len(cluster)
        		output.write('\t'.join(map(str,CS_mean)))
        		output.write('\n')
        output.close()

    elif args.which == 'all':

        assignment = load_file(args.input)

        assignment = filter_assignment(assignment,args.percent)

        numstates = max(assignment) + 1

        with open(args.output, 'w') as output:
            for i in range(numstates):
                cluster = np.where(assignment == i)[0]
                CS_mean = np.mean(CS[cluster,:], axis = 0)
                print len(cluster)
                output.write('\t'.join(map(str,CS_mean)))
                output.write('\n')













if __name__ == '__main__':
    main()
