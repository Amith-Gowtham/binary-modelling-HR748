# Calling necessary modules
from matplotlib.pyplot import plot
import matplotlib.pylab as plt
from matplotlib.axes import Axes
import numpy as np
import csv
from math import sqrt

# Constants (Period and T_0s)
per = 7.6409
tsv = 9385.7 + (per/2)
tsb = 9385.7  + (per/2)
tsu = 9278.7 + (per/2) + 0.0275


# Defining necessary functions

def phase_v(x):        
    return (x - tsv)/per - np.floor((x - tsv)/per)
def phase_b(x):        
    return (x - tsb)/per - np.floor((x - tsb)/per)
def phase_u(x):        
    return (x - tsu)/per - np.floor((x - tsu)/per)

# Importing Data 
t1, t2, t3 = [], [], []
vmag, bmag, umag = [], [], []
fig, ax = plt.subplots()
with open('CygVMag.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        t1.append(float(row[0]))
        vmag.append(float(row[1]))
with open('CygBMag.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        t2.append(float(row[0]))
        bmag.append(float(row[1]))
with open('CygUMag.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        t3.append(float(row[0]))
        umag.append(float(row[1]))

ph_1 = list(map(phase_v, t1))           
ph_2 = list(map(phase_b, t2))
ph_3 = list(map(phase_u, t3))

ax.plot(ph_1, vmag, '.k', label = 'V Magnitudes')
ax.plot(ph_2, bmag, '.b', label = 'B Magnitudes')
ax.plot(ph_3, umag, '.r', label = 'U Magnitudes')

ax.set_xlabel('Phase')
ax.set_ylabel('Mag')
ax.invert_yaxis()
ax.legend()
ax.grid(True)
plt.show()
