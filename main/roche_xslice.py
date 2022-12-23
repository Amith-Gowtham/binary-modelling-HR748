# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 12:12:47 2020

@author: Amith
"""
import numpy as np
from numpy import linspace, meshgrid, array
import matplotlib.pyplot as plt
from sympy import *
from cmath import *

x = Symbol('x')

# Constants
y, z = 0, 0 
e , q = 0.53193, 1/0.9707  
nu = 0
nud = 0
xlims, res = [-3.0, 10], 1000

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
    return (f(e)**2)*((1 + e)**4)/(1 + e*np.cos(nu))**3

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
    s = (r**2 - 2*x + 1)**(1/2)
    F = x - (q/r) - (1/s)
    G = -(1/2)*(x**2 + y**2)*(1+q)*A(e,nu)
    V = (F+G)
    return V

def x_slice(x, y, z, xlims, res):
    """
    Function to plot the x-slice of the equipotential curve is to be plotted.
    
    Parameters
    ----------
    y, z : float
        Specific y and z-coordinates 
    xlims : array or ndarray
        Upper limit and Lower limit of cartesian x.
    res : int
        No. of points between upper and lower limit of x-axis.
    
    Returns
    -------
    None.

    """
    vx = diff(effpot(x, 0, 0, q, e, nu), x)                                                                                                                                                                                                                                                             
    sol = solve(simplify(vx), x)
    real_sol = [i for i in sol if phase(i) == 0 or phase(i) == 3.141592653589793]                                                
    for j in real_sol:
        if 0<j<1: L1 = j
        if j>1: L2 = j
        if j<0: L3 = j
    x = linspace(xlims[0], xlims[1], res)
    v = effpot(x, y, z, q, e, nu)
    plt.annotate(r'$\nu$ = %r degrees'% nud, xy=(0.05, 0.95), xycoords='axes fraction')
    plt.annotate(r'$\mathcal{A}$ = %e'% A(e, nu), xy=(0.05, 0.9), xycoords='axes fraction')
    plt.plot(x, v, 'k')
    plt.plot(L1, effpot(L1, 0, 0, q, e, nu), 'og', label='L1')
    plt.plot(L2, effpot(L2, 0, 0, q, e, nu), 'ob', label='L2')
    plt.plot(L3, effpot(L3, 0, 0, q, e, nu), 'or', label='L3')
    plt.legend()
    plt.show()
    return None

# Command to plot the graph in question.
x_slice(x, y, z, xlims, res)



