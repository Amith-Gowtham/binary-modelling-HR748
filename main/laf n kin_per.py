# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 20:14:04 2020

@author: Amith
"""

import numpy as np
import csv
from matplotlib import pyplot as plt
import time

# Setting up Arrays for calculations
ti, mag, ph = [], [], []
with open('CygVMag.txt','r') as rawfile:            # Dataset .txt file            
    plots = csv.reader(rawfile, delimiter='\t')
    for row in plots:
        ti.append(float(row[0]))
        mag.append(float(row[1]))
        
# Defining Constants
LIMITS = [7.0001, 8.000]
STEP_SIZE = 0.0001
T_ST = min(ti)        

# Defining functions
def laf_n_kin(l, ss, t_st, mag, ti):
    """
    This routine prints the periodogram for LK statistic as defined in Lafler 
    and Kinman (1965) and finds the most probable period for the given dataset
    within the range.

    Parameters
    ----------
    l : array or ndarray
        List whose elements are the lower and upper limit of periodogram.
    ss : float
        Step-size between any two test periods within the periodogram.
    t_st : float
        T_0 from which the origin of phase is determined.
    mag, ti : array or ndarray
        Lists containing the time-varying values of magnitudes and thier 
        corresponding time value in either M.J.D or H.J.D 
        
    Returns
    -------
    None

    """
    l = np.arange(l[0], l[1], ss)      
    theta = []
    for j in range(0, len(l)):
        def phase(x):
            return (x - t_st)/per - np.floor((x - t_st)/per)
        per = l[j]
        ph = list(map(phase, ti))
        c = []
        for i in range(0, len(mag)):
            d = [mag[i], ph[i]]
            c.append(d)
        c.sort(key = lambda x: x[1])
        xs = [x[0] for x in c]
        m = np.mean(xs) 
        a1, b1 = [], []
        for i in range(0, len(xs)): 
            a = (xs[i-1] - xs[i])**2
            b = (xs[i] - m)**2
            a1.append(a)
            b1.append(b)
        a2, b2 = sum(a1), sum(b1)
        theta.append(a2/b2)
        e1 = []
        for i in range(0, len(theta)): 
            e = [l[i], theta[i]]
            e1.append(e)
        e1.sort(key = lambda x: x[1])
    print("Minimised theta value:", e1[0][1])
    print("Corresponding Period value:", e1[0][0], '+/-', ss)  
    plt.plot(l,theta, '.r', label = 'Trial Periods')
    plt.xlabel(r'Test Period')
    plt.ylabel(r'$\Theta$')
    plt.plot(e1[0][0], e1[0][1], '*b', label = r'Period with the least $\Theta$')
    plt.legend()
    plt.show()  
    return None      
    
# Code to print Graphical Output
st = time.time()
laf_n_kin(LIMITS, STEP_SIZE, T_ST, mag, ti)
et = time.time()
print("Program Run-time: ", et-st, ' seconds')    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




    
