#!/usr/bin/python

from progressbar import ProgressBar
import progressbar
import argparse
import numpy as np
import json
import os,sys
from scipy.stats import norm
from matplotlib import pyplot as plt
from matplotlib import cm
import multiprocessing
from multiprocessing import Pool
import warnings

widgets = ['Progress: ', progressbar.Percentage(), ' ', progressbar.Bar(marker=progressbar.RotatingMarker()),' ', progressbar.ETA(), ' ',progressbar.SimpleProgress()]

class Bayesian_control:
    ''' the control variables of MCMC'''
    def __init__(self, bayesian_input):
        try:
            self.Prior_type = bayesian_input['prior']
        except KeyError:
            print "'prior:' does not exist in bayesian.json"
        try:
            self.Number_states = bayesian_input['number_states']
        except KeyError:
            print "'number_states:' does not exist in bayesian.json"
        try:
            self.Step_width = bayesian_input['step_width']
        except KeyError:
            print "'step_width:' does not exist in bayesian.json"
        try:
            self.Number_steps = bayesian_input['number_steps']
        except KeyError:
            print "'number_steps:' does not exist in bayesian.json"
        try:
            self.Number_experimental_data = bayesian_input['number_experimental_data']
        except KeyError:
            print "'number_experimental_data:' does not exist in bayesian.json"
        try:
            self.Confor_name = bayesian_input['conformation_name']
        except KeyError:
            print "'conformation_name:' does not exist in bayesian.json"
    def print_control(self):
        ''' print out the control information'''
        print "*"*20+"Control Parameters"+"*"*20
        print "\nPrior type of Bayesian: " + self.Prior_type
        print "\nNumber of states: "+ str(self.Number_states)
        print "\nConformation Names: " + '\t'.join(name for name in self.Confor_name)
        print "\nNumber of experimental data: " + str(self.Number_experimental_data)
        print "\nNumber of steps: " + str(self.Number_steps)
        print "\nStep width of Monte Carlo: " + str(self.Step_width)
        print "-"*60    

class Prior_parameters:
    ''' the parameters for prior '''
    def __init__(self,control,theoretical_input):
        if control.Prior_type == 'RC':
            self.Ini_weight = np.array([1.0/control.Number_states]*control.Number_states)
            self.Ini_weight_sd = np.array([0.2]*control.Number_states)
        else:
            try:
                self.Ini_weight = np.array(theoretical_input['md_mean'])
            except KeyError:
                print "'md_mean:' does not exist in theoretical.json"
            try:
                self.Ini_weight_sd = np.array(theoretical_input['md_sd'])
            except KeyError:
                print "'Initial weight uncertainty:' does not exist in theoretical.json"

    def print_prior(self,names):
        ''' print out prior information'''
        print "*"*20+"\tPrior Parameters\t"+"*"*20
        print '\t\t'+'\t'.join(name for name in names)
        print 'w_ini\t\t'+'\t'.join("%.3f" %weight for weight in self.Ini_weight)
        print 'w_ini sd\t'+'\t'.join("%.2f" %weight_sd for weight_sd in self.Ini_weight_sd)
        print "-"*60

class Likelihood_parameters:
    ''' the parameters for likelihood'''
    def __init__(self,theoretical_input, experimental_input):
        ''' obtain likelihood parameters and check'''
        try:
            self.Comp_data = np.array(theoretical_input['theo_mean'])
        except KeyError:
            print "'theo_mean:' does not exist in theoretical.json"
        try:
            self.Comp_sd = np.array(theoretical_input['theo_sd'])
        except KeyError:
            print "'theo_sd:' does not exist in theoretical.json"
        try:
            self.Exp_mean = np.array(experimental_input['exp_mean'])
        except KeyError:
            print "'exp_mean:' does not exist in experimental.json'"
        try:
            self.Exp_sd = np.array(experimental_input['exp_sd'])
        except KeyError:
            print "'exp_sd:' does not exist in experimental.json'"

    def print_likelihood(self,names):
        '''print out likelihood information'''
        print "*"*20+"\tLikelihood Parameters\t"+"*"*20
        print "\nComputational data:"
        for i in range(len(self.Comp_data)):
            row = self.Comp_data[i]
            print names[i]+'\t'+'\t'.join("%.3f" %item for item in row)
        
        print "\nExperimental data:"
        print '\t'+'\t'.join("%.3f" %exp_data for exp_data in self.Exp_mean)
        
        print "\nComputational data SD:"
        print '\t'+'\t'.join("%.3f" %comp_data_sd for comp_data_sd in self.Comp_sd)

        print "\nExperimental data SD:"
        print '\t'+'\t'.join("%.3f" %exp_data_sd for exp_data_sd in self.Exp_sd)                
        print "-"*60

class MCMC_data:
    ''' the class to save MCMC results '''
    def __init__(self,*args):
        if type(args[0]) is not str : # load matrix
            self.Weight = args[0] # matrix save mcmc sampled weight
            self.get_mean_sd() # calculate mean and sd for weight
            if len(args) > 1:
                self.Accept_rate = args[1]
        else: # load weight from .npy file
            self.load(args[0])
    
    def get_mean_sd(self):
        self.Mean = self.Weight.mean(axis = 0)
        self.SD = self.Weight.std(axis = 0)
    
    def save(self, out_dir, file_name = 'posterior_MCMC.txt'):
        '''save mcmc results to file_name.npy file'''
        print "saved MCMC results to file 'posterior_MCMC.txt' "
        np.save(os.path.join(out_dir,file_name), self.Weight)
        
    def load(self, in_dir, file_name):
        '''load mcmc results from file_name'''
        self.Weight = np.load(file_name)
        self.get_mean_sd()
        
    def print_mcmc(self,names):
        '''print a summary of MCMC results'''
        print "*"*20+"\tPosterior:\t"+"*"*20
        print '\t\t'+'\t'.join(name for name in names)
        print 'w_posterior\t'+'\t'.join("%.3f" %weight for weight in self.Mean)
        print 'w_posterior SD\t'+'\t'.join("%.3f" %weight_sd for weight_sd in self.SD)
        print "-"*60

    def plot_mcmc(self,names,out_dir,plot_file = "Bayesian_posterior_hist.png"):
        '''plot the distribution of MCMC'''
        colors = cm.rainbow(np.linspace(0, 1, len(names)))
        fig = plt.figure(figsize=(12,8))
        fontsize = 15

        for i in range(len(names)):            
            binEdges= plt.hist(self.Weight[:,i],bins=20,facecolor=colors[i],alpha =0.15,histtype = 'bar', edgecolor = colors[i], rwidth = 0.9, normed = 1);
            bincenters = 0.5*(binEdges[1][1:]+binEdges[1][:-1])
            plt.plot(bincenters,binEdges[0],lw = 2,color = colors[i],label = names[i]+'\t('+"%.2f"%(self.Mean[i]*100)+r'$\pm$'+"%.2f"%(self.SD[i]*100)+')%')
        plt.xlabel('Conformation Ratio',fontsize = fontsize)
        plt.ylabel('Probability',fontsize = fontsize)
        plt.ylim(ymin = 0.0 , ymax = 20)
        plt.xlim(xmin = 0.0, xmax = 1.0)
        plt.title("Posterior",position = (0.3,0.9),fontsize = fontsize+20)
        plt.legend(loc= 1, ncol = 1, prop={'size':fontsize},frameon = False)
        fig.savefig(os.path.join(out_dir,plot_file))


def bayesian_serial(in_dir,out_dir):
    ''' serial bayesian '''
    
    # read and parse control
    with open(os.path.join(in_dir,'bayesian.json'),'r') as input_file:
        bayesian_input = json.load(input_file)
    control = Bayesian_control(bayesian_input)
    control.print_control()

    # read experimental file
    with open(os.path.join(in_dir, 'experimental.json'), 'r') as input_file:
        experimental_input = json.load(input_file)

    # read theoretical file
    with open(os.path.join(in_dir, 'theoretical.json'), 'r') as input_file:
        theoretical_input = json.load(input_file)

    # generate prior variables
    prior_variables = Prior_parameters(control,theoretical_input)
    prior_variables.print_prior(control.Confor_name)

    # generate likelihood variables
    likelihood_variables = Likelihood_parameters(theoretical_input, experimental_input)
    likelihood_variables.print_likelihood(control.Confor_name)

    #running MCMC for Bayesian
    MCMC_result = MCMC(control,prior_variables,likelihood_variables)
    MCMC_result = MCMC_data(*MCMC_result)

    # save and plot MCMC results
    MCMC_result.save(out_dir)
    MCMC_result.plot_mcmc(control.Confor_name,out_dir)

    MCMC_result.print_mcmc(control.Confor_name)
    cal_chi_square(prior_variables,likelihood_variables,MCMC_result)

def bayesian_paralle(in_dir,out_dir):

    # read and parse control
    with open(os.path.join(in_dir,'bayesian.json'),'r') as input_file:
        bayesian_input = json.load(input_file)
    control = Bayesian_control(bayesian_input)
    control.print_control()

    # read experimental file
    with open(os.path.join(in_dir, 'experimental.json'), 'r') as input_file:
        experimental_input = json.load(input_file)

    # read theoretical file
    with open(os.path.join(in_dir, 'theoretical.json'), 'r') as input_file:
        theoretical_input = json.load(input_file)

    # generate prior variables
    prior_variables = Prior_parameters(control,theoretical_input)
    prior_variables.print_prior(control.Confor_name)

    # generate likelihood variables
    likelihood_variables = Likelihood_parameters(theoretical_input, experimental_input)
    likelihood_variables.print_likelihood(control.Confor_name)

     #running MCMC for Bayesian parallel
    nprocess = multiprocessing.cpu_count()
    pool = Pool(processes=nprocess)   

    control.Number_steps /= nprocess

    print 'You have {0:1d} CPUs to run MCMC'.format(nprocess)
    print "each threat run ",control.Number_steps," steps"

    MCMC_result = pool.map_async(MCMC_helper, [(control,prior_variables,likelihood_variables)]*nprocess).get(9999999)

    MCMC_result = Combine_MCMC_result(MCMC_result,nprocess)

    MCMC_result = MCMC_data(*MCMC_result)

    # save and plot MCMC results
    MCMC_result.save(out_dir)
    MCMC_result.plot_mcmc(control.Confor_name,out_dir)

    MCMC_result.print_mcmc(control.Confor_name)
    cal_chi_square(prior_variables,likelihood_variables,MCMC_result)



def Cal_Comp_Data(weight, Comp_data):
    '''calculate computaional mean from theoretical data and weight'''
    return np.dot(weight, Comp_data)

def Prior(weight, prior_variables):
    ''' prior distribution '''
    return np.prod(norm.pdf(weight,prior_variables.Ini_weight,prior_variables.Ini_weight_sd))

def Likelihood(weight, likelihood_variables):
    ''' likelihood distribution '''
    Comp_mean = Cal_Comp_Data(weight, likelihood_variables.Comp_data)
    Comp_Exp_sd = np.sqrt(likelihood_variables.Comp_sd**2+likelihood_variables.Exp_sd**2)
    return np.prod(norm.pdf(likelihood_variables.Exp_mean, Comp_mean, Comp_Exp_sd))

def Posterior(weight, prior_variables, likelihood_variables):
    ''' Posterior distribution '''
    v_Prior = Prior(weight,prior_variables)
    v_Likelihood = Likelihood(weight,likelihood_variables)
    if min(weight) < 0:
        return -np.inf
    else:
        return v_Prior * v_Likelihood

##### MCMC Part #####
def MCMC(control, prior_variables, likelihood_variables):
    ''' the main part of MCMC processes'''
    np.random.seed()
    accept_criteria = np.random.uniform(0,1,control.Number_steps)
    random_walk_steps = np.reshape(np.random.uniform(-control.Step_width,control.Step_width,control.Number_steps*(control.Number_states-1)),(control.Number_steps,control.Number_states-1))
    weight_mcmc = np.zeros((control.Number_steps,control.Number_states))
    weight_mcmc[0,:] = prior_variables.Ini_weight
    temp_steps = np.zeros((1,control.Number_states))
    accept_steps = 0
    print "\nRunning MCMC....."
    pbar = ProgressBar(widgets=widgets)
    for i in pbar(range(1, control.Number_steps)):
        temp_steps[0,:control.Number_states-1] = weight_mcmc[i-1, : control.Number_states-1] + random_walk_steps[i,:]
        temp_steps[0,control.Number_states-1] = 1 - sum( temp_steps[0,:control.Number_states-1])
        if Posterior(weight_mcmc[i-1,:],prior_variables,likelihood_variables) * accept_criteria[i] <= Posterior(temp_steps[0,:],prior_variables,likelihood_variables):
            weight_mcmc[i,] = temp_steps[0,]
            accept_steps += 1
        else:
            weight_mcmc[i,] = weight_mcmc[i-1,]

    print "\nFinished MCMC"
    accept_rate = float(accept_steps)/control.Number_steps * 100
    print "\naccept rates: \t", accept_rate,"%"
    if accept_rate > 70:
        warnings.warn('The accept rate is too high, please increase the step size.')
    elif accept_rate < 30:
        warnings.warn('The accept rate is too low, please decrease the step size.')
    return (weight_mcmc, accept_rate)    

def MCMC_helper(args):
    '''help function for parallel MCMC'''
    return MCMC(*args)

def Combine_MCMC_result(MCMC_result,nprocess):
    ''' combining results from parallel MCMC'''
    weight_mcmc = MCMC_result[0][0]
    accept_rate = MCMC_result[0][1]
    for i in range(1,nprocess):
        weight_mcmc = np.vstack((weight_mcmc,MCMC_result[i][0]))
        accept_rate += MCMC_result[i][1]
    accept_rate /=nprocess
    return (weight_mcmc,accept_rate)


def cal_chi_square(prior_variables,likelihood_variables,MCMC_result):
    ''' calculate chi square and output'''
    Comp_Exp_var = likelihood_variables.Comp_sd**2+likelihood_variables.Exp_sd**2
    chi_square_prior = np.mean((Cal_Comp_Data(prior_variables.Ini_weight, likelihood_variables.Comp_data)-likelihood_variables.Exp_mean)**2/Comp_Exp_var)
    chi_square_posterior = np.mean((Cal_Comp_Data(MCMC_result.Mean, likelihood_variables.Comp_data)-likelihood_variables.Exp_mean)**2/Comp_Exp_var)    
    print "*"*20+"\tX^2\t"+"*"*20
    print "X^2"
    print "Prior\tPosterior"
    print "%.2f\t%.2f"%(chi_square_prior,chi_square_posterior)
    print "-"*60


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-p',
        '--mp',
        choices = ['serial', 'parallel'],
        default = 'serial',
        help = "Using single processor or parallel, default is single"
        )

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

    args = parser.parse_args()

    if args.mp == "serial":
        bayesian_serial(args.in_dir, args.out_dir)
    else:
        bayesian_paralle(args.in_dir,args.out_dir)





if __name__  == '__main__':
    main()
