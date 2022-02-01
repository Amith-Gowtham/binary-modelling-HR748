 # -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 10:26:52 2020

@author: Amith
"""

import phoebe
import numpy as np
import matplotlib.pyplot as plt
import emcee
import sys
import corner
import time
from phoebe import u,c

logger = phoebe.logger(clevel='WARNING')

b = phoebe.default_binary()
b.add_constraint('semidetached', 'secondary')
b['period@orbit'] = 7.6409
b['ecc@orbit'] = 0.53193
b['q@orbit'] = 0.9707
b['vgamma@system'] = -15.9991
b['t0@system'] = 39255.796
b['asini@constraint']=22.45689
b['sma@orbit'] = 22.63
b['teff@primary'] = 6673.45
b['teff@secondary'] = 6379.71
b['incl@orbit'] = 85.33

# Loading data into PHOEBE
logfile = '/mcmc.txt'
lcv = np.loadtxt('CygVMag.txt') # .txt file of the observation set)
b.add_dataset('lc', times=lcv[:,0],
              fluxes=lcv[:,1], 
              sigmas=0.05*np.ones(len(lcv)),
              passband='Johnson:V')
phoebe.interactive_checks_off()
phoebe.interactive_constraints_off() 
b.set_value_all('irrad_method', 'horvat')
b.flip_constraint('compute_phases', 'compute_times')
b['compute_phases@lc@dataset'] = np.linspace(-0.5,0.5,21)
print(b)

# Defining essential functions and routines
def lnprob(x, adjpars, priors):
    # Check to see that all values are within the allowed limits:
    # if not np.all([priors[i][0] < x[i] < priors[i][1] for i in range(len(priors))]):
    #   return -np.inf
    for i in range(len(adjpars)):
        b[adjpars[i]] = x[i]   
    # Let's assume that our priors are uniform on the range of the physical parameter combinations.
    # This is already handled in Phoebe, which will throw an error if the system is not physical,
    # therefore it's easy to implement the lnprior as =0 when system checks pass and =-inf if they don't.
    # Here we'll 'package' this in a simple try/except statement:    
    try:
        b.run_compute(irrad_method='none')
        # sum of squares of the residuals
        fluxes_model = b['fluxes@model'].interp_value(times=lcv[:,0])
        lnp = -0.5*np.sum((fluxes_model-b['value@fluxes@dataset'])**2 / b['value@sigmas@dataset']**2) 
    except:
        lnp = -np.inf
    sys.stderr.write("lnp = %e\n" % (lnp))
    return lnp

def run(adjpars, priors, nwalkers, niter):
    ndim = len(adjpars)
    with open(logfile, "w") as f:
        f.write('# Number of parameters being adjusted: %d\n' % (len(adjpars)))
        f.write('# \n')
        f.write('# %15s %14s %14s\n' % ('Parameter:', 'Lower limit:', 'Upper limit:'))
        for i in range(len(adjpars)):
            f.write('# %15s %14.5f %14.5f\n' % (adjpars[i], priors[i][0], priors[i][1]))
        f.write('# \n')
    p0 = np.array([[p[0] + (p[1]-p[0])*np.random.rand() for p in priors] for i in range(nwalkers)])
#     pool = MPIPool()
#     if not pool.is_master():
#         pool.wait()
#         sys.exit(0)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[adjpars, priors])
    for result in sampler.sample(p0, iterations=niter, storechain=False):
        position = result[0]
        with open(logfile, 'a') as f:
            for k in range(position.shape[0]):
                f.write("%d %s %f\n" % (k, " ".join(['%.12f' % i for i in position[k]]), result[1][k]))
#     pool.close()

chain = np.loadtxt(logfile) # Text file to record process logs
def plot_posteriors(chain, paramnames, truths, skipiters=10):
    samples = np.delete(chain, 0, 1)
    samples = np.delete(samples, -1, 1)
    fig = corner.corner(samples[nwalkers*skipiters:], labels=paramnames, quantiles=[0.16, 0.5, 0.84],
                       show_titles=True, title_kwargs={"fontsize": 12}, truths=truths )
    plt.show()
    
# Code to choose parameters and their priors whose values are to be estimated 
adjpars = ['sma@orbit', 'teff@primary', 'teff@secondary', 'incl@orbit']
priors = [(22.4, 22.9), (6000, 7000), (6000, 7000), (80, 90)]
nwalkers = 32
niters = 10   
state = None

# Code to run the PHOEBE-MCMC routine
time1 = time.time()
run(adjpars, priors, nwalkers, niters)
time2 = time.time()

# Code to print results 
plot_posteriors(chain, ['a', 'T\u209A', 'T\u209B', 'i'], truths = [22.44, 6550, 6496, 87.035], skipiters=0)
print('Time taken to run:', (time2-time1)/60)   

