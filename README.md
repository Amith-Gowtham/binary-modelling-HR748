# binary-modelling-HR7484
A repository for codes written in Python to analyze, study and model the Eclipsing Variable V1143 Cygni (HR 7484) for my Bachelor of Science thesis in Physics. 
This repository contains five Python programmes in all and they are listed as they would appear chronologically in the aforesaid thesis and are briefly described.


1. laf n kin_per.py 

	Periodogram code for Lafler and Kinman Method of finding periodicity of a periodic/semi-periodic variable star. This program basically accepts observable 
datasets (tested only for small-medium sized sampling datasets in Photometry like magnintudes) and gives out the most likely value for periodicity and the graphical output of the 
periodogram for the given range.


2. ph_plot.py

	This programme is more of a companion to the previous code. This programme accepts the value of periodicity of a variable star (here, V1143 Cygni) and plots its 
Phased Light Curve (provided that one has supplied the data points).


3. fitting.py

	This is more of a modification to the source code for estimating values of certain photometric parameters of the variable star in question, obtained from a workshop repository 
in PHOEBE project's GitHub page. The only changes which are made are the application of constraints, selecting parameters and setting their priors for the MCMC run etc. among some 
notes and docstrings.


4. roche_xslice.py and roche_xyslice.py

	These two codes are used for the sole purpose of plotting the Equilibrium points and their corresponding zero-velocity curves on the specified z-level and y-level(for the former) 
and specified z-level (for the latter). The first programme describes the change in potential of the system along any horizontal parallel to the arbitrary x-axis and the second 
programme does the same but on a plane parallel to the arbitrary xy-plane.

The reader is welcome to download these programmes and can suggest any improvements (if any) for which I will be grateful.



