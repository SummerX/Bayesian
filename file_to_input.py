#!/usr/bin/python

import numpy as np
import os
import json
import argparse

def experiment_input(in_dir,out_dir):
	print "Start generate experimental file for bayesian"
	experimental = {}
	experimental['exp_mean'] = np.loadtxt(os.path.join(in_dir,"exp_mean.dat")).tolist()
	experimental['exp_sd'] = np.loadtxt(os.path.join(in_dir,"exp_sd.dat"))[0,:].tolist()
	with open(os.path.join(out_dir,'experimental.json'),'w') as output:
		json.dump(experimental,output,indent = 4, sort_keys = True)
	print "Finish generate experimental.json!"  

def theoretical_input(in_dir,out_dir):
	print "Start generate theoretical file for bayesian"
	theoretical = {}
	theoretical['md_mean'] = np.loadtxt(os.path.join(in_dir,'md_mean.dat')).tolist()
	theoretical['md_sd'] = np.loadtxt(os.path.join(in_dir,'md_sd.dat')).tolist()
	theoretical['theo_mean'] = np.loadtxt(os.path.join(in_dir,'theo_mean.dat')).tolist()
	theoretical['theo_sd'] = np.loadtxt(os.path.join(in_dir,"exp_sd.dat"))[1,:].tolist()
	with open(os.path.join(out_dir,'theoretical.json'),'w') as output:
		json.dump(theoretical,output,indent = 4, sort_keys = True)
	print "Finished generate theoretical.json!" 

def bayesian_input(in_dir,out_dir,prior,nc,nd,ns,ss):
	print "Start generate bayesian input file"
	bayesian = {}
	bayesian['prior'] = prior
	bayesian['number_states'] = nc
	bayesian['number_experimental_data'] = nd
	bayesian['number_steps'] = ns
	bayesian['step_width'] = ss
	bayesian['conformation_name'] = np.genfromtxt(os.path.join(in_dir,'conform.dat'),dtype = str).flatten().tolist()
	with open(os.path.join(out_dir,'bayesian.json'),'w') as output:
		json.dump(bayesian, output, indent = 4, sort_keys = True)
	print "Finished generate bayesian.json!"



def main():
	parser = argparse.ArgumentParser()

	parser.add_argument('-i',
		'--in_dir',
		default = '.',
		help = "input directory, (default is local directory))"
		)

	parser.add_argument('-o',
		'--out_dir',
		default = '.',
		help = "output directory, (default is local directory))"
		)

	parser.add_argument('-p',
		'--prior',
		choices = ['RC','MD'],
		default = 'RC',
		help = "Prior used for Bayesian, default is RC prior"
		)

	parser.add_argument('-nc',
		'--num_states',
		default = 0,
		type = int,
		help = "number of states"
		)

	parser.add_argument('-nd',
		'--num_data',
		type = int,
		default = 0,
		help = "number of experimental data"
		)

	parser.add_argument('-ns',
		'--num_steps',
		type = int,
		default = 100000,
		help = "number of steps for MCMC, default = 100000"
		)

	parser.add_argument('-ss',
		'--step_size',
		type = float,
		default = 0.02,
		help = "step width for MCMC, default = 0.02"
		)

	args = parser.parse_args()

	experiment_input(args.in_dir,args.out_dir)
	theoretical_input(args.in_dir,args.out_dir)
	bayesian_input(args.in_dir,args.out_dir,args.prior,args.num_states,args.num_data,args.num_steps,args.step_size)


if __name__ == "__main__":
	main()