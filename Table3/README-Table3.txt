
          ===== Table 3 =====

This directory contains instructions how to reproduce Table 3
from Blecic et al. (2015a), TEA: A code for calculating
thermochemical equilibrium abundances. 

Data listed in the Table 3 are taken from the headers/ and results/ 
directories of the respected TEA directories (runs), and the CEA 
result file (run) placed  in the run/ directory, done to produce
Fig. 7 right panel and Fig. 8 right panel.

CEA free energies - Fig8/Right_panel/CEA-TEA-freeEnergiesCEA/
                    /TEArun-freeEnergiesCEA/headers/
                    take headers with temperatures 1720, 1805, 1810 K
                    and respective pressures of 3.80e-01, 1.66e+00, 
                    2.19e+01 bar, and read the values of the species free
                    energies from the column 'Chemical potential'
TEA (JANAF) free energies - Fig7/Right_panel/CEA-TEA/TEA/headers/
                    take headers with temperatures 1720, 1805, 1810 K
                    and respective pressures of 3.80e-01, 1.66e+00, 
                    2.19e+01 bar, and read the values of the species free 
                    energies from the column 'Chemical potential'
CEA final abundances - Fig8/Right_panel/CEA-TEA-freeEnergiesCEA/run/
                     open CEA-TEA.dat find temperatures 1720, 1805, 1810 K
                     and respective pressures of 3.80e-01, 1.66e+00, 
                     2.19e+01 bar and read the values of the species final 
                     abundances
TEA final abundances using CEA free energies - Fig8/Right_panel/
                    /CEA-TEA-freeEnergiesCEA/TEArun-freeEnergiesCEA/results/
                    open TEArun-FreeEnergiesCEA.tea and read rows with 
                    temperatures 1720, 1805, 1810 K and respective pressures 
                    of 3.80e-01, 1.66e+00, 2.19e+01 bar
TEA final abundances using JANAF free energies - Fig7/Right_panel/CEA-TEA/
                    /TEA/results/
                    open CEA-TEA-freeEnergiesTEA.tea and read rows with 
                    temperatures 1720, 1805, 1810 K and respective pressures 
                    of 3.80e-01, 1.66e+00, 2.19e+01 bar
