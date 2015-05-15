
          ===== Figure 6 =====

This directory contains instructions, data, and the code to reproduce 
Fig. 6 from Blecic et al. (2015a), TEA: A code for calculating
thermochemical equilibrium abundances. 

The directory consists the following:
 - WASP-43b-solar\    - TEA run with input and output files for WASP43b
                        solar metallicity case
 - WASP-43b-10xsolar\ - TEA run with input and output files for WASP43b
                        10 times solar metallicity case
 - WASP-43b-50xsolar\ - TEA run with input and output files for WASP43b
                        50 times solar metallicity case
 - WASP43b_PHASE7_NH3_H2S_TP.txt - P-T data for WASP43b,first and second
                        column
 - plotTEA.py - code to plot the figures

The T-P profile is provided by Kevin Stevenson and Mike Line, and is used in 
the paper by Stevenson et al. (2014): 
http://www.sciencemag.org/content/346/6211/838.full.pdf
DOI: 10.1126/science.1256758

To reproduce each TEA run, open atm_inputs/ subdirectory from the 
respective directories, make run directory, copy TEA.cfg to that directory, 
edit TEA.cfg with paths to the files provided in the respective atm_inputs/
directory and run TEA for each case separately as: 

<path_to_TEA>/tea/runatm.py <path_to_preatm.atm> desired_name


To plot Fig. 6 use plotTEA.py function provided.
Edit the code between lines 171 and 212 and use a set of code lines for the 
desired metallicity case that you want to plot. Keep the rest of the lines
commented out. Currently, the WASP43b-solar set of code lines is active in 
the code, the rest is commented out.
 
To run plotTEA.py, provide the path to the each result file and type 
in terminal:
plotTEA.py <../results/NAME.tea> H,CO,CO2,CH4,H2O,HCN,C2H2,C2H4,N2,NH3,HS,H2S


