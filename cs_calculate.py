#!/usr/bin/python
import mdtraj as md
import os
from distutils.spawn import find_executable
import numpy as np
from multiprocessing import Pool
import argparse
import multiprocessing
import shutil


def convert_traj_to_pdb(traj,out_dir):
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	for i in range(traj.n_frames):
		traj[i].save_pdb(os.path.join(out_dir,'traj%d.pdb'%i))

def calculate_shiftx2(out_dir,pH,temperature):
	files = os.path.join(out_dir, 'traj*.pdb')
	shiftx2 = find_executable('shiftx2.py')
	if not shiftx2:
		print "Can not find shiftx2!"
	
	file_ref = os.path.join(out_dir, 'traj0.pdb')
	tmpfolder = os.path.join(out_dir,'tmp')
	if not os.path.exists(tmpfolder):
		os.makedirs(tmpfolder)
	cmd = "%s -b '%s' -p %.1f -t %.2f -i %s -f csv -z %s" %(shiftx2, files, pH, temperature,file_ref, tmpfolder)
	print cmd
	os.system(cmd)


def gathering_cs(traj,out_dir):
	atom_indices = np.loadtxt('nmr_atomindices.dat',dtype = np.int).flatten()
	cs = np.empty([0,len(atom_indices)],dtype = float)
	for i in range(traj.n_frames):
		filename = os.path.join(out_dir,'traj%d.pdb.cs'%i)
		data = np.genfromtxt(filename,delimiter=",",skip_header = 1, usecols = 3)[atom_indices]
		cs = np.vstack((cs,data))	
	filename = os.path.basename(out_dir)+'_cs.txt'
	print filename
	np.savetxt(filename, cs,fmt = '%.3f')
	shutil.rmtree(out_dir)


def cal_cs_helper(args):
	'''traj,out_dir,pH,temperature'''
	traj, out_dir, pH, temperature = args
	convert_traj_to_pdb(traj,out_dir)
	calculate_shiftx2(out_dir,pH,temperature)
	return gathering_cs(traj,out_dir)

def combine_data(output,ind_max):
	data = np.loadtxt('PDB_0_cs.txt')
	for i in range(1,ind_max):
		data = np.vstack((data, np.loadtxt('PDB_%d_cs.txt'%i)))
	np.savetxt(output,data,fmt = '%.3f')



def main():
	parser = argparse.ArgumentParser()

	parser.add_argument('-ph',
		'--pH',
		default = 5.0,
		type = float,
		help = " PH used to calculate cs, default = 5.0"
		)

	parser.add_argument('-T',
		'--temperature',
		default = 298.0,
		type = float,
		help = " temperature, default = 298.0"
		)

	parser.add_argument('-t',
		'--trajectory',
		default = 'trajectory/traj.h5',
		help = "trajectory file used to calculate cs, default is trajectory/traj.h5"
		)

	parser.add_argument('-o',
		'--output',
		default = 'cs.dat',
		help = "output file to save output chemical shifts, default is cs.dat"
		)

	parser.add_argument('-m',
		'--mode',
		choices = ['serial','parallel'],
		default = 'parallel',
		help = 'mode used to calculate shiftx2, single or parallel'
		)

	parser.add_argument('-b',
		'--begin',
		default = 0,
		type = int,
		help = "start point of trajectory to calculate shiftx2"
		)

	args = parser.parse_args()

	

	if parser.mode == 'parallel':
		traj = md.load(os.path.join(os.getcwd(),args.trajectory))

		chuck_size = 300
		
		ind_max = traj.n_frames/chuck_size

		variable = [(traj[i*chuck_size:(i+1)*chuck_size],'PDB_%d'%i,args.pH,args.temperature) for i in range(ind_max)]

		if ind_max*chuck_size < traj.n_frames:
			variable.append((traj[ind_max*chuck_size:],'PDB_%d'%ind_max,args.pH,args.temperature))
			ind_max += 1

		nprocess = multiprocessing.cpu_count()
		pool = Pool(processes=nprocess)   
		print 'You have {0:1d} CPUs to run MCMC'.format(nprocess)
		pool.map_async(cal_cs_helper, variable).get(9999999)
		combine_data(args.output,ind_max)

	else:
		traj = md.load(os.path.join(os.getcwd(),args.trajectory))[args.begin:]
		cal_cs_helper((traj,'PDB',args.pH, args.temperature))



if __name__ == "__main__":
	main()
