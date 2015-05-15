
          ===== Figure 8 =====

This directory contains instructions, data, and the code to reproduce 
Fig. 8 from Blecic et al. (2015a), TEA: A code for calculating
thermochemical equilibrium abundances. 

The directory consists the following:
 - Left_panel/  - directory containing data to reproduce left panel of Fig 8
 - Right_panel/ - directory containing data to reproduce right panel of Fig 8


 ####### Left_panel ######

This directory contains the CEA-TEA-freeEnergiesCEA/ subdirectory that 
carries all information needed to reproduce left panel of Fig. 8. 
It contains the following subdirectories and files:

 - runCEA-TEA-comparison.py - code to produce Fig. 8 left panel
 - CEA/ - directory that contains the CEA source files used to produce
          this figure
 - CEA-inputs/ - directory that contains CEA input files produced after
          running the first section of runCEA-TEA-comparison.py
 - CEA-outputs/ - directory that contains CEA output files produced after
          running the first section of the runCEA-TEA-comparison.py
 - TEArun-freeEnergiesCEA/ - directory that contains the TEA run with input
          and output files that uses CEA free energies
 - run/ - directory that contains the CEA result file produced after running
          the second section of the runCEA-TEA-comparison.py code
 - makeheader-freeEnergiesCEA.py - makeheader.py source file with CEA free 
          energies calculated using thermo.inp and the instructions given
          in the CEA theory document, Section 4, and in the TEA theory
          document, Section 7.1
 - TEA-CEA-thermo.py - code to calculate CEA free energies based on the 
          thermo.inp and the instructions given in the CEA theory document,
          Section 4, and in the TEA theory document, Section 7.1  
        

To reproduce all above directories and their data do the following:

 ========= make new makeheader.py with CEA thermo data ===========
 - open terminal and type
ipython
 
 - go to TEA/tea/ source directory and rename makeheader.py to 
   makeheader_original.py
 - open TEA-CEA-thermo.py and copy and paste all lines into the ipython session
 - copy final output (CEA free energies) and the temperature range  in the  
   the makeheader.py line 548 
 - TEA is run with this new makeheader.py to produce final abundances based
   on the CEA thermochemical data given
 - makeheader-freeEnergiesCEA.py is the file with the CEA free energies already
   in
 - so to reproduce, take makeheader-freeEnergiesCEA.py rename it to 
   makeheader.py
 - copy newly made makeheader.py to the TEA/tea/ source

 ========= run TEA  ========= 
 - make run/ directory
 - open TEArun-freeEnergiesCEA/inputs/ directory and copy TEA.cfg to the
   run/ directory 
 - edit TEA.cfg with paths to the files provided in inputs/ directory
 - run TEA as: 

<path_to_TEA>/tea/runatm.py <path_to_preatm.atm> desired_name

 ========== make plots ============
 - open runCEA-TEA-comparison.py in text editor 
 - open terminal and type:
ipython

 - edit line 135 with the correct location of the FCEA2-cmd file
 - run section == run CEA ==, copy and paste lines 1-140 into the
   ipython session
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

 - run all other sections from the runCEA-TEA-comparison.py 
 - copy and paste lines 142 - 437
 - the CEA-TEA.dat file is made inside the run directory
 - the plots CEA-TEA-freeEnergiesCEA.png and CEA-TEA-freeEnergiesCEA.ps 
   will be made



 ####### Right_panel ######

This directory contains the CEA-TEA-freeEnergiesCEA subdirectory that 
carries all information needed to reproduce right panel of Fig. 8. 
It contains the following subdirectories and files:

 - runCEA-TEA-comparison.py - code to produce Fig. 8 right panel
 - CEA/ - directory that contains the CEA source files used to produce
          this figure
 - CEA-inputs/ - directory that contains CEA input files produced after
          running the first section of runCEA-TEA-comparison.py
 - CEA-outputs/ - directory that contains CEA output files produced after
          running the first section of the runCEA-TEA-comparison.py
 - TEArun-freeEnergiesCEA/ - directory that contains the TEA run with input
          and output files that uses CEA free energies
 - run/ - directory that contains the CEA result file produced after running
          the second section of the runCEA-TEA-comparison.py code
 - makeheader-freeEnergiesCEA.py - makeheader.py source file with CEA free 
          energies calculated using thermo.inp and the instructions given
          in the CEA theory document, Section 4, and in the TEA theory
          document, Section 7.1
 - TEA-CEA-thermo.py - code to calculate CEA free energies based on the 
          thermo.inp and the instructions given in the CEA theory document,
          Section 4, and in the TEA theory document, Section 7.1  

To reproduce all above directories and their data do the following:

 ========= make new makeheader.py with CEA thermo data ===========
 - open terminal and type
ipython
 
 - go to TEA/tea/ source directory and rename makeheader.py to 
   makeheader_original.py
 - open TEA-CEA-thermo.py and copy and paste all lines into the ipython session
 - copy final output (CEA free energies) and the temperature range  in the  
   the makeheader.py line 548 
 - TEA is run with this new makeheader.py to produce final abundances based
   on the CEA thermochemical data given
 - makeheader-freeEnergiesCEA.py is the file with the CEA free energies already
   in
 - so to reproduce, take makeheader-freeEnergiesCEA.py rename it to 
   makeheader.py
 - copy newly made makeheader.py to the TEA/tea/ source

 ========= run TEA  ========= 
 - make run/ directory
 - open TEArun-freeEnergiesCEA/inputs/ directory and copy TEA.cfg to the
   run/ directory 
 - edit TEA.cfg with paths to the files provided in inputs/ directory
 - run TEA as: 

<path_to_TEA>/tea/runatm.py <path_to_preatm.atm> desired_name


 ========== make plots =====================
 - open in text editor runCEA-TEA-comparison.py 
 - open terminal and type:
ipython

 - edit line 184 with the correct location of the FCEA2-cmd file
 - run section == run CEA ==, copy and paste lines 1-188 into 
   the ipython session
 - exit the ipython session and make two directories outside of the CEA/:
   CEA-inputs and CEA-outputs
 - move all *.inp files into ./inputs directory
mv *.inp ../CEA-inputs/

 - move all *.out files into ./outputs directory
 mv *.out ../CEA-outputs/

 - edit line 197 with the correct path to the ./CEA-outputs directory
 - make directory run
mkdir run

 - cd to run
cd run

 - open terminal and type
ipython

 - run all other sections from the runCEA-TEA-comparison.py 
 - copy and paste lines 191 - 571
 - the CEA-TEA.dat file is made inside the run directory
 - the plots CEA-TEA-WASP43b-freeEnergiesCEA.png and
   CEA-TEA-WASP43b-freeEnergiesCEA.ps will be made


