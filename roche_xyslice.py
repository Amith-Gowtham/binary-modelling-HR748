# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 16:45:07 2020

@author: Amith
"""

import matplotlib.pyplot as plt
import mpmath as mp
from cmath import *
import numpy as np
from sympy import *
import time

# Initialising Constants and Symbols
x, y, z = symbols('x y z')
e = 0.53193
nu = 0
nud = 0
q = 1/0.9707

# Specify the upper and lower limits of x and y-axis and the resolution for the roche graph
xlims = [-3.0, 13.0]
ylims = [-10.0, 10.0]
res = [5000, 5000]

# Necessary Functions
def f(e):
    """
    Function to return the value of 'f' as defined in (Sepinsky, 2007, Pg. 4)
    Redefined in terms in 'e'. It shows whether the system is Synchronous/
    Asynchronous.
    (Assuming Synchronicity params F_1 = F_2 = 1). 
    
    Parameters
    ----------
    e : float
        Eccentricity of the binary system in question. Unitless.

    Returns
    -------
    f : float
        Rotational angular velocity of star 1 normalised to the value of 
        w_p. Unitless.

    """
    return ((1-e)**(3/2))/((1+e)**(1/2))

def A(e, nu):
    """
    To define the factor which depends on f, e, v(true anomaly) as defined in 
    (Sepinsky, 2007)

    Parameters
    ----------
    nu : float
        True Anomaly of Star 2 wrt Star 1. Radians.

    Returns
    -------
    float
        Equation 21 in (Sepinsky, 2007)

    """
    a = (f(e)**2)*(((1+e)**4)/(1+e*np.cos(nu))**3)
    return a


def effpot(x, y, z, q, e, nu):
    """
    To return the value of Effective Potential at a point in space(x,y,z)

    Parameters
    ----------
    x, y, z : float
        Cartesian Coordinates at which Eff.Pot is to be calculated.
    q : float
        Mass ratio of the Binary System (M_1/M_2). Unitless.
    e : float
        Eccentricity of the binary system in question. Unitless.    
    nu : float
        True Anomaly of Star 2 wrt Star 1. Radians.    

    Returns
    -------
    V : float
        Effective Potential as defined by Equation 20 in (Sepinsky, 2007).

    """
    r = (x**2 + y**2 + z**2)**(1/2)
    s = (r**2 -2*x + 1)**(1/2)
    F = x - (q/r) - (1/s)
    G = -(1/2)*(x**2 + y**2)*(1+q)*A(e,nu)
    V = (F+G)
    return V

def roche_plot(xlims, ylims, res):
    
    """
    Function that plots the Equipotential curves that passes through Lagrange
    Points (L_[1-5])on the x-axis.

    Parameters
    ----------
    xlims, ylims : array or ndarray
        Upper limit and Lower limit of caryesian x and y axis
    res : array or ndarray
        An array whose first element is no. of points between upper and lower
        limit of x-axis and likewise for the second element in y-axis.

    Returns
    -------
    None.

    """
    if len(xlims) == 2 and len(ylims) == 2:
        vx = diff(effpot(x, 0, 0, q, e, nu), x)                                                                                                                                                                                                                                                             
        sol = solve(simplify(vx), x)
        real_sol = [i for i in sol if phase(i) == 0 or phase(i) == 3.141592653589793]                                                
        for j in real_sol:
            if 0<j<1: L1 = j
            if j>1: L2 = j
            if j<0: L3 = j
        alp = (1+q)*A(e, nu)   
        xf = mp.power(mp.power(q/(alp-1), 2), 1/3)/2
        yf = mp.power(xf*(2 - xf), 1/2)
        levels = []
        for i in real_sol:
            d = effpot(i, 0, 0, q, e, nu)                                              
            levels.append(d)
        levels.sort()   
        print(levels)                                                           
        X = np.linspace(xlims[0], xlims[1], res[0])
        Y = np.linspace(ylims[0], ylims[1], res[1])
        X, Y = np.meshgrid(X, Y)
        Z = effpot(X, Y, 0, q, e, nu)
        plt.contour(X, Y, Z, levels=levels, linestyles='solid', colors='k')
        plt.annotate(r'$\nu$ = %r degrees'% nud, xy=(0.05, 0.95), xycoords='axes fraction')
        plt.annotate(r'$\mathcal{A}$ = %e'% A(e, nu), xy=(0.05, 0.9), xycoords='axes fraction')
        plt.plot(L3, 0, '.r', label='L3')                                  
        plt.plot(L1, 0, '.g', label='L1')                                      
        plt.plot(L2, 0, '.b', label='L2') 
        plt.plot(xf, yf, '.y', label='L4') 
        plt.plot(xf, -yf, '.m', label='L5')          
        plt.plot(0, 0, 'ok', label=r'$M_p$')                            
        plt.plot(1, 0, 'xk', label=r'$M_s$')
        plt.legend()
        plt.show()
    else:
        print("Invalid entry of axes' limits")                
    return None

# Code to plot Roche Potential Graph on z=0 xy plane 
st = time.time()
roche_plot(xlims, ylims, res)
et = time.time()
print('Programme Runtime: ',et-st)
