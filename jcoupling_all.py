#!/usr/bin/python
from karplus import *
import mdtraj as md
import argparse
import os
import numpy as np

def main():
	parser = argparse.ArgumentParser()

	parser.add_argument('-o',
	    '--output',
	    default = 'Data/J_coup_kalplus.dat',
	    help = "output file name for trajectory (default = Data/J_coup_kalplus.dat)"
	    )

	parser.add_argument('-t',
	    '--trajectory',
	    default='trajectory',
	    help ='Path to the trajectory directory or project file (default=trajectory)'
	    )

	parser.add_argument('-f',
		'--function',
		default = 'karp_function.dat',
		help = "Function used to calculate Jcoupling default is karp_function.dat")

	parser.add_argument('-d',
		'--angle',
		default = 'karp_angle.dat',
		help = "angle used to calculate Jcoupling default is karp_angle.dat (phi first, psi second)")

	args = parser.parse_args()
	traj = md.load(os.path.join(os.getcwd(),args.trajectory,'traj.h5'))

	function = np.loadtxt(args.function,dtype = str)
	dispatcher = {'J3HNHa':J3HNHa,'J3HNC':J3HNC,'J3HaC':J3HaC,'J3CC':J3CC,'J3HNCb':J3HNCb,'J1NCa':J1NCa,'J2NCa':J2NCa,'J3HNCa':J3HNCa}

	angles = np.loadtxt(args.angle,dtype = int)

	phi = md.compute_phi(traj)[1] * 180.0 / np.pi
	psi = md.compute_psi(traj)[1] * 180.0 / np.pi

	snap=traj.n_frames
	nj=len(function)

	print 'There are '+str(nj)+' jcoupling and totally '+str(snap)+' snapshot'

	J = np.zeros((snap, nj))

	for i in range(nj):
		print i
		J[:,i] = map(dispatcher[function[i]],phi[:,angles[i,0]],psi[:,angles[i,1]])

	if not os.path.isdir(os.path.dirname(args.output)):
		os.mkdir(os.path.dirname(args.output))
		
	np.savetxt(args.output,J,newline='\n',fmt='%f')


if __name__ == '__main__':
	main()