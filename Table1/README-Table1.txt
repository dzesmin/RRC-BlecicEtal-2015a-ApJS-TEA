
          ===== Table 1 =====

This directory contains instructions, data, and the code to reproduce 
Table 1 from Blecic et al. (2015a), TEA: A code for calculating
thermochemical equilibrium abundances. 

The directory consists the following: 

 - TEA.cfg - config file to run TEA 
 - input_white.txt - input file for TEA
 - header_white.txt - prepared header file for TEA
 - format_white.py - source code adjusted to give moles instead of mole 
         fractions as the final output
 - runsingle.py - source code adjusted to use header_white.txt with 
         White et al. (1958) free energy values, instead of JANAF 
          thermo data
 - white/ - directory containing full TEA run 

To reproduce this run, and get TEA results using White's free energy values do
the following:
 - go to the TEA/tea/ source directory and rename format.py to format_original.py
 - rename format_white.py to format.py and place in the TEA/tea/ source directory 
 - do the same procedure with runsingle.py
 - edit runsingle.py line 138 and give the path to the 'header_white.txt' file
 - make run directory
mkdir run

 - cd to run
cd run

 - copy TEA.cfg to run/ directory
 - copy input_white.txt to run/ directory
 - edit TEA.cfg with correct paths to the inputs files
 - run TEA as:
<path_to_TEA>/tea/runsingle.py input_white.txt white

 - 'white' subdirctory will be created and in the results/ subdirectory
   'results-visual.txt' will list final mole numbers in the column 'Final Abun'
 - those are the values in the 4th column of the Table 1.
 - first, second and third column are from White et al. (1958) paper, and the
   last column is difference between White's and TEA results


