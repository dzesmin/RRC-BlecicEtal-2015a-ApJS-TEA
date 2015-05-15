
          ===== Figure 5 =====

This directory contains instructions, data, and the code to reproduce 
Fig. 5 from Blecic et al. (2015a), TEA: A code for calculating
thermochemical equilibrium abundances.

The directory consists the following:
 - /WASP12b-O-rich  - TEA run with input and output files for
                      WASP12b O-rich case, left panel of Figure 5
 - /WASP12b-C-rich  - TEA run with input and output files for
                      WASP12b C-rich case, right panel of Figure 5
 - tp_orp.dat - P-T data for WASP12b, O-rich case, third and fourth columns
 - tp_crp.dat - P-T data for WASP12b, C-rich case, third and fourth columns
 - plotTEA.py - code to plot the figures

The T-P profiles are provided by Kevin Stevenson and are used in the
paper by Stevenson et al. (2014): 
http://iopscience.iop.org/0004-637X/791/1/36/pdf/apj_791_1_36.pdf
DOI:10.1088/0004-637X/791/1/36

To reproduce each TEA run, open atm_inputs/ subdirectory from the 
respective directories, make run directory, copy TEA.cfg to that directory, 
edit TEA.cfg with paths to the files provided in the respective atm_inputs/
directory and run TEA for each case separately as: 

<path_to_TEA>/tea/runatm.py <path_to_preatm.atm> <output_dir>

To plot Fig. 5 use plotTEA.py function provided.
Edit the code between lines 171-197 and 242-280 use a set of code lines for 
the desired C/O case that you want to plot. Keep the rest of the lines
commented out. Currently, the WASP12b O-rich set of code lines is active in 
the code, the rest is commented out.
 
To run plotTEA.py, provide the path to the each result file and type
in terminal:
plotTEA.py <../results/NAME.tea> H,CO,CO2,CH4,H2O,HCN,C2H2,C2H4,N2,NH3,HS,H2S

