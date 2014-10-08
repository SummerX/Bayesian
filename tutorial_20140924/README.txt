Bayesian: 

prerequrest:
	progressbar
	json
	numpy
	matplotlib
	multiprocessing
	scipy

input file:
	bayesian.json:
		controling file for bayesian in json format
		"conformation_name" : list of conformation names according to theoretical data
		"number_experimental_data" : number of experimental data included
		"number_states" : number of sub-states
		"number_steps" : number of steps for MCMC
		"prior" : prior type (RC or MD)
		"step_width" : step size for MCMC

	experimental.json:
		experimental results in json format
		"exp_mean" : list of experimental data
		"exp_sd" : standard deviation for experimental data

	theoretical.json:
		theoretical results in json format
		"md_mean" : list of MD prior for each sub-states
		"md_sd" : list of MD stadard deviation for each sub-states
		"theo_mean" : matrixs of theoretical calculated experimental variables for each conformation
		"theo_sd" : standard deviation for each theoretical calculated experimental data


usage: Bayesian.py [-h] [-p {serial,parallel}] [-i IN_DIR] [-o OUT_DIR]

optional arguments:
  -h, --help            show this help message and exit
  -p {serial,parallel}, --mp {serial,parallel}
                        Using single processor or parallel, default is single
  -i IN_DIR, --in_dir IN_DIR
                        input directory, (default is local directory))
  -o OUT_DIR, --out_dir OUT_DIR
                        output directory, (default is local directory))


output:
	posterior_MCMC.txt.npy:
		MCMC output in numpy format
	Bayesian_posterior_hist.png:
		figure of posterior distribution

	print out on screen:
		summary of input data, progress of process, summary of posterior, chi^2





