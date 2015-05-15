
          ===== Figure 7 =====

This directory contains instructions, data, and the code to reproduce 
Fig. 7 from Blecic et al. (2015a), TEA: A code for calculating
thermochemical equilibrium abundances. 

The directory consists the following:
 - Left_panel/  - directory containing data to reproduce left panel of Fig 7
 - Right_panel/ - directory containing data to reproduce right panel of Fig 7


 ####### Left_panel #######

This directory contains the CEA-TEA-TEBS/ subdirectory that carries all 
information needed to reproduce left panel of Fig. 7. It contains the 
following subdirectories and files:

 - runCEA-TEA-TEBS-comparison.py - code to produce Fig. 7 left panel
 - CEA/ - directory that contains the CEA source files used to produce
          this figure
 - CEA-inputs/ - directory that contains CEA input files produced after
          running the first section of runCEA-TEA-TEBS-comparison.py
 - CEA-outputs/ - directory that contains CEA output files produced after
          running the first section of the runCEA-TEA-TEBS-comparison.py
 - TEA/ - directory that contains the TEA run with input and output files
 - run/ - directory that contains the CEA result file produced after running
          the second section of the runCEA-TEA-TEBS-comparison.py code

To reproduce all above directories and their data do the following:

 === run TEA ===
Make run directory, open TEA/atm_inputs/ subdirectory, copy TEA.cfg to the
run/ directory, edit TEA.cfg with paths to the files provided in atm_inputs/
directory and run TEA as: 

<path_to_TEA>/tea/runatm.py <path_to_preatm.atm> desired_name

 === make plots ===
 - open runCEA-TEA-TEBS-comparison.py in text editor 
 - open terminal and type:
ipython

 - edit line 135 with the correct location of the FCEA2-cmd file
 - run section == run CEA ==, copy and paste lines 1-140 into the 
   ipython session
 - multiple .inp and .out files will be made in the CEA/ directory
 - exit the ipython session and make two directories outside of the CEA/:
   CEA-inputs/ and CEA-outputs/
 - move all *.inp files into ./CEA-inputs directory
mv *.inp ../CEA-inputs/

 - move all *.out files into ./CEA-outputs directory
 mv *.out ../CEA-outputs/

 - edit line 148 with the correct path to the ./CEA-outputs directory
 - make directory run
mkdir run

 - cd to run
cd run

 - open terminal and type
ipython

 - run all other sections from the runCEA-TEA-TEBS-comparison.py 
 - copy and paste lines 142 - 510 into the ipython session
 - the CEATEATEBS.dat file is made inside the run/ directory
 - the plots CEA-TEA-TEBS.png and CEA-TEA-TEBS.ps will be made


 ####### Right_panel ######

This directory contains the CEA-TEA subdirectory that carries all information
needed to reproduce right panel of Fig. 7. It contains the following 
subdirectories and files:

 - runCEA-TEA-comparison.py - code to produce Fig. 7 right panel
 - CEA/ - directory that contains the CEA source files used to produce
          this figure
 - CEA-inputs/ - directory that contains CEA input files produced after
          running the first section of runCEA-TEA-comparison.py
 - CEA-outputs/ - directory that contains CEA output files produced after
          running the first section of the runCEA-TEA-comparison.py
 - TEA/ - directory that contains the TEA run with input and output files
 - run/ - directory that contains the CEA result file produced after running
          the second section of the runCEA-TEA-comparison.py code

To reproduce all above directories and their data do the following:

 === run TEA ===
Make run directory, open TEA/atm_inputs/ subdirectory, copy TEA.cfg to the
run/ directory, edit TEA.cfg with paths to the files provided in atm_inputs/
directory and run TEA as: 

<path_to_TEA>/tea/runatm.py <path_to_preatm.atm> desired_name

 === make plots ===
 - open in text editor runCEA-TEA-comparison.py 
 - open terminal and type:
ipython

 - edit line 184 with the correct location of the FCEA2-cmd file
 - run section == run CEA ==, copy and paste lines 1-188 into the
   ipython session
 - exit the ipython session and make two directories outside of the CEA/:
   CEA-inputs/ and CEA-outputs/
 - move all *.inp files into ./CEA-inputs directory
mv *.inp ../CEA-inputs/

 - move all *.out files into ./CEA-outputs directory
 mv *.out ../CEA-outputs/

 - edit line 197 with the correct path to the ./CEA-outputs directory
 - make directory run
mkdir run

 - cd to run
cd run

 - open terminal and type
ipython

 - run all other sections from the runCEA-TEA-comparison.py 
 - copy and paste lines 191 - 571 into the ipython session
 - the CEA-TEA.dat file is made inside the run directory
 - the plots CEA-TEA-WASP-43b.png and CEA-TEA-WASP-43b.ps will be made



