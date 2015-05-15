
          ===== Table 2 =====

This directory contains instructions how to reproduce Table 2
from Blecic et al. (2015a), TEA: A code for calculating
thermochemical equilibrium abundances. 

Data listed in the Table 2 are taken from the headers/ and results/ 
directories of the respected TEA directories (runs), and the CEA 
result file (run) placed  in the run/ directory, done to produce
Fig. 7 left panel and Fig. 8 left panel.

CEA free energies - Fig8/Left_panel/CEA-TEA-freeEnergiesCEA/
                    /TEArun-freeEnergiesCEA/headers/
                    take headers with temperatures 2500, 2700, 2900 K
                    and read the values of the species free energies 
                    in the column 'Chemical potential'
TEA (JANAF) free energies - Fig7/Left_panel/CEA-TEA-TEBS/TEA/headers/
                    take headers with temperatures 2500, 2700, 2900 K 
                    and read the values of the species free energies 
                    in the column 'Chemical potential'
CEA final abundances- Fig8/Left_panel/CEA-TEA-freeEnergiesCEA/run/
                     open CEA-TEA.dat find 2500, 2700, 2900 K rows
                     and read the values of the species final abundances
TEA final abundances using CEA free energies - Fig8/Left_panel/
                    /CEA-TEA-freeEnergiesCEA/TEArun-freeEnergiesCEA/results/
                    open CEA-TEA-comp.tea and read rows with 2500, 2700, 2900 K
TEA final abundances using JANAF free energies - Fig7/Left_panel/CEA-TEA-TEBS/
                    /TEA/results/
                    open CEA-TEA.tea and read rows with 2500, 2700, 2900 K
