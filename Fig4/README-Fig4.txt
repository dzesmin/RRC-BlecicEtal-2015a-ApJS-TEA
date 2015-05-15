
          ===== Figure 4 =====

This directory contains instructions, data, and the code to reproduce 
Fig. 4 from Blecic et al. (2015a), TEA: A code for calculating
thermochemical equilibrium abundances.

The directory contains the temperature and pressure data of WASP12b 
and WASP43b, and the code to plot the figures:
 - PTprofiles.py - code to plot the figures
 - tp_orp.dat - P-T data for WASP12b, O-rich case, third and fourth columns
 - tp_crp.dat - P-T data for WASP12b, C-rich case, third and fourth columns
 - WASP43b_PHASE7_NH3_H2S_TP.txt - WASP43b 

The T-P profiles are provided by Kevin Stevenson and Michael Line and are 
used in the papers by Stevenson et al. (2014a):  
http://iopscience.iop.org/0004-637X/791/1/36/pdf/apj_791_1_36.pdf
DOI:10.1088/0004-637X/791/1/36
and Stevenson et al. (2014b):
http://www.sciencemag.org/content/346/6211/838.full.pdf
DOI: 10.1126/science.1256758

To reproduce Fig. 4 run PTprofiles.py in terminal as:
PTprofiles.py

