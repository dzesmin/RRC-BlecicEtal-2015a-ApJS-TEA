#! /usr/bin/env python

# ******************************* START LICENSE *******************************
# Thermal Equilibrium Abundances (TEA), a code to calculate gaseous molecular
# abundances in planetary atmospheres under thermochemical equilibrium
# conditions.
#
# This project was completed with the support of the NASA Earth and Space 
# Science Fellowship Program, grant NNX12AL83H, held by Jasmina Blecic, 
# PI Joseph Harrington. Project developers included graduate student 
# Jasmina Blecic and undergraduate M. Oliver Bowman. 
# 
# Copyright (C) 2014 University of Central Florida.  All rights reserved.
# 
# This is a test version only, and may not be redistributed to any third
# party.  Please refer such requests to us.  This program is distributed
# in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.
# 
# We welcome your feedback, but do not guarantee support.  Please send
# feedback or inquiries to both:
# 
# Jasmina Blecic <jasmina@physics.ucf.edu>
# Joseph Harrington <jh@physics.ucf.edu>
# 
# or alternatively,
# 
# Jasmina Blecic and Joseph Harrington
# UCF PSB 441
# 4000 Central Florida Blvd
# Orlando, FL 32816-2385
# USA
# 
# Thank you for testing TEA!
# ******************************* END LICENSE *******************************

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from PIL import Image


def plotTEA():
    '''
    This code plots a figure of temperature vs. abundances for any multiTP 
    final output produced by TEA, RESULT_ATM_FILE. It accepts 2 arguments 
    on the command line: the path to the multiTP TEA output, and the names
    of the species user wants to plot. See notes for the full description
    of arguments. To run the code do:
    ../TEA/tea/plotTEA.py <RESULT_ATM_FILE> <SPECIES_NAMES>
    See 'Notes' section for an example.

    Parameters
    ----------
    None

    Returns
    -------
    plot_out - string
               Name and location of output plot.

    Notes
    -----
    To plot the results, go to the results/ directory and open final
    .tea file. Read the species names and column numbers you want to plot.
    Arguments given should be in the following format:
    ../TEA/tea/plotTEA.py <RESULT_ATM_FILE> <SPECIES_NAMES>

    RESULT_ATM_FILE - string
               Path to the atm_file.           
    SPECIES_NAMES  - list of strings
               List of species that user wants to plot. Species names should 
               be given with their symbols (without their states) and no
               breaks between species names (e.g., CH4,CO,H2O).

    Example: ../TEA/tea/plotTEA.py ../TEA/doc/examples/multiTP/results/multiTP_Example.tea CO,CH4,H2O,NH3
    The plot is opened once execution is completed and saved in the ./plots/
    subdirectory.

    The lower range on the y axis (mixing fraction) should be set at most -14 
    (the maximum precision of the TEA code). The lower range on the x axis
    (temperature) is currently set to 200 K, although the best precision is
    reached above 1000 K.
    '''
    
    # Get plots directory, create if non-existent
    plots_dir = "plots/"
    if not os.path.exists(plots_dir): os.makedirs(plots_dir)
    
    # Counts number of arguments given
    noArguments = len(sys.argv)

    # Prints usage if number of arguments different from 4 or adiabatic profile
    if noArguments != 3:
        print("\nUsage  : ../TEA/tea/plotTEA.py atmfile species(divided by comma, no breaks)")
        print("Example: ../TEA/tea/plotTEA.py ../TEA/doc/examples/multiTP/results/multiTP_Example.tea CO,CH4,H2O,NH3\n")
     
    # Sets the first argument given as the atmospheric file
    filename = sys.argv[1]

    # Sets the second argument given as the species names
    species = sys.argv[2]

    # Open the atmospheric file and read
    f = open(filename, 'r')
    lines = np.asarray(f.readlines())
    f.close()

    # Get molecules names
    imol = np.where(lines == "#SPECIES\n")[0][0] + 1
    molecules = lines[imol].split()
    nmol = len(molecules)
    for m in np.arange(nmol):
        molecules[m] = molecules[m].partition('_')[0]

    # Take user input for species and split species strings into separate strings 
    #      convert the list to tuple
    species = tuple(species.split(','))
    nspec   = len(species)

    # Populate column numbers for requested species and 
    #          update list of species if order is not appropriate
    columns = []
    spec    = []
    for i in np.arange(nmol):
        for j in np.arange(nspec):
            if molecules[i] == species[j]:
                columns.append(i+2)
                spec.append(species[j])

    # Convert spec to tuple
    spec = tuple(spec)

    # Concatenate spec with pressure for data and columns
    data    = tuple(np.concatenate((['p'], spec)))
    usecols = tuple(np.concatenate(([0], columns)))

    # Load all data for all interested species
    data = np.loadtxt(filename, dtype=float, comments='#', delimiter=None,    \
                    converters=None, skiprows=8, usecols=usecols, unpack=True)

    # Open a figure
    plt.figure(1)
    plt.clf()
 
    # WASP43b all plots (please DO NOT USE THE SAME COLORS FOR YOUR PAPER!)
    colors = 'b', '#FF1CAE','#FF0000' , '#FFAA00', '#00FFFF', '#00FF00', '#91219E', '#BCEE68' , 'g', '#ffc3a0', 'c','m'
    color_index = 0

    # Plot all specs with different colours and labels
    for i in np.arange(nspec):
        # addition for WASP-43b species that does not change with metallicity
        if i==0 or i==3 or i==5 or i==6 or i==7 or i==9: 
            plt.loglog(data[i+1], data[0], '--', color=colors[color_index], \
                                        linewidth=2, label=str(spec[i]))
        else:
            plt.loglog(data[i+1], data[0], '-', color=colors[color_index], \
                                        linewidth=3, label=str(spec[i]))
        color_index += 1

    # Label the plot
    plt.xlabel('Mixing Fraction'    , fontsize=14)
    plt.ylabel('Pressure [bar]'    , fontsize=14)

    # Font-size of x and y ticks
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    # WASP43b solar
    plt.text(2e-8, 1e-4, 'H', color='b', fontsize=14)
    plt.text(8e-4, 1, 'CO', color='#FF1CAE', fontsize=14)
    plt.text(5e-7,1e-4, 'CO$_{2}$', color='#FF0000', fontsize=14)
    plt.text(2e-12,3e-5, 'CH$_{4}$', color='#FFAA00', fontsize=14)
    plt.text(8e-4,1e-3, 'H$_{2}$O', color='#00FFFF', fontsize=14)
    plt.text(7e-14,2e-4, 'HCN', color='#00FF00', fontsize=14)
    plt.text(2e-19,2e-5, 'C$_{2}$H$_{2}$', color='#91219E', fontsize=14)
    plt.text(4e-19,6e-4, 'C$_{2}$H$_{4}$', color='#BCEE68', fontsize=14)
    plt.text(7e-5,1e-4, 'N$_{2}$', color='g', fontsize=14)
    plt.text(1.5e-10,2e-5, 'NH$_{3}$', color='#ffc3a0', fontsize=14)
    plt.text(2.5e-9, 6e-6, 'HS', color='c', fontsize=14)
    plt.text(1e-6,1e-5, 'H$_{2}$S', color='m', fontsize=14)

    # WASP43b 10xsolar
    #plt.text(2e-6, 2e-2, 'H', color='b', fontsize=14)
    #plt.text(8e-3, 1, 'CO', color='#FF1CAE', fontsize=14)
    #plt.text(1.5e-6,1e-4, 'CO$_{2}$', color='#FF0000', fontsize=14)
    #plt.text(9e-12,8e-5, 'CH$_{4}$', color='#FFAA00', fontsize=14)
    #plt.text(8e-3,1e-3, 'H$_{2}$O', color='#00FFFF', fontsize=14)
    #plt.text(3.5e-13,4e-4, 'HCN', color='#00FF00', fontsize=14)
    #plt.text(2e-19,2e-5, 'C$_{2}$H$_{2}$', color='#91219E', fontsize=14)
    #plt.text(4e-21,1.5e-4, 'C$_{2}$H$_{4}$', color='#BCEE68', fontsize=14)
    #plt.text(8e-4,1e-4, 'N$_{2}$', color='g', fontsize=14)
    #plt.text(4e-10,7e-4, 'NH$_{3}$', color='#ffc3a0', fontsize=14)
    #plt.text(3e-7, 3e-5, 'HS', color='c', fontsize=14)
    #plt.text(1e-5,3, 'H$_{2}$S', color='m', fontsize=14)

    # WASP43b 50xsolar
    #plt.text(4e-8, 2e-5, 'H', color='b', fontsize=14)
    #plt.text(4e-2, 1, 'CO', color='#FF1CAE', fontsize=14)
    #plt.text(3e-5,1e-4, 'CO$_{2}$', color='#FF0000', fontsize=14)
    #plt.text(4e-12,6e-6, 'CH$_{4}$', color='#FFAA00', fontsize=14)
    #plt.text(4e-2,1e-3, 'H$_{2}$O', color='#00FFFF', fontsize=14)
    #plt.text(5e-13,3e-4, 'HCN', color='#00FF00', fontsize=14)
    #plt.text(2e-19,2e-5, 'C$_{2}$H$_{2}$', color='#91219E', fontsize=14)
    #plt.text(4e-21,1.5e-4, 'C$_{2}$H$_{4}$', color='#BCEE68', fontsize=14)
    #plt.text(3.5e-3,1e-4, 'N$_{2}$', color='g', fontsize=14)
    #plt.text(8e-7,3e-1, 'NH$_{3}$', color='#ffc3a0', fontsize=14)
    #plt.text(8.5e-7, 5e-5, 'HS', color='c', fontsize=14)
    #plt.text(3e-5,1e-5, 'H$_{2}$S', color='m', fontsize=14)

    # ======================= Kevin's PT profile ========================== #

    # Kevin's PT profile from WASP43b_PHASE7_NH3_H2S_TP.txt
    # WASP43b, Stevenson et al 2014
    pres = np.array([  3.16227770e+01,   2.63026800e+01,   2.18776160e+01,
         1.81970090e+01,   1.51356120e+01,   1.25892540e+01,
         1.04712850e+01,   8.70963590e+00,   7.24435960e+00,
         6.02559590e+00,   5.01187230e+00,   4.16869380e+00,
         3.46736850e+00,   2.88403150e+00,   2.39883290e+00,
         1.99526230e+00,   1.65958690e+00,   1.38038430e+00,
         1.14815360e+00,   9.54992590e-01,   7.94328230e-01,
         6.60693450e-01,   5.49540870e-01,   4.57088190e-01,
         3.80189400e-01,   3.16227770e-01,   2.63026800e-01,
         2.18776160e-01,   1.81970090e-01,   1.51356120e-01,
         1.25892540e-01,   1.04712850e-01,   8.70963590e-02,
         7.24435960e-02,   6.02559590e-02,   5.01187230e-02,
         4.16869380e-02,   3.46736850e-02,   2.88403150e-02,
         2.39883290e-02,   1.99526230e-02,   1.65958690e-02,
         1.38038430e-02,   1.14815360e-02,   9.54992590e-03,
         7.94328230e-03,   6.60693450e-03,   5.49540870e-03,
         4.57088190e-03,   3.80189400e-03,   3.16227770e-03,
         2.63026800e-03,   2.18776160e-03,   1.81970090e-03,
         1.51356120e-03,   1.25892540e-03,   1.04712850e-03,
         8.70963590e-04,   7.24435960e-04,   6.02559590e-04,
         5.01187230e-04,   4.16869380e-04,   3.46736850e-04,
         2.88403150e-04,   2.39883290e-04,   1.99526230e-04,
         1.65958690e-04,   1.38038430e-04,   1.14815360e-04,
         9.54992590e-05,   7.94328230e-05,   6.60693450e-05,
         5.49540870e-05,   4.57088190e-05,   3.80189400e-05,
         3.16227770e-05,   2.63026800e-05,   2.18776160e-05,
         1.81970090e-05,   1.51356120e-05,   1.25892540e-05,
         1.04712850e-05,   8.70963590e-06,   7.24435960e-06,
         6.02559590e-06,   5.01187230e-06,   4.16869380e-06,
         3.46736850e-06,   2.88403150e-06,   2.39883290e-06])

    # Kevin's PT profile from WASP43b_PHASE7_NH3_H2S_TP.txt
    # WASP43b, Stevenson et al 2014
    temp = np.array([ 1811.8938 ,  1810.9444 ,  1810.1535 ,  1809.4948 ,  1808.9463 ,
        1808.4898 ,  1808.1098 ,  1807.7936 ,  1807.5304 ,  1807.3114 ,
        1807.1291 ,  1806.9766 ,  1806.8464 ,  1806.7212 ,  1806.53   ,
        1806.1269 ,  1805.2849 ,  1803.7403 ,  1800.5841 ,  1794.9518 ,
        1786.6255 ,  1775.3705 ,  1761.0973 ,  1742.3631 ,  1719.6396 ,
        1694.0976 ,  1666.1517 ,  1636.2055 ,  1603.0265 ,  1567.6227 ,
        1531.3326 ,  1494.5529 ,  1457.6432 ,  1419.5923 ,  1381.4921 ,
        1344.4864 ,  1308.8483 ,  1274.7949 ,  1241.6997 ,  1210.3302 ,
        1181.2851 ,  1154.5844 ,  1130.2066 ,  1107.7735 ,  1087.5396 ,
        1069.5885 ,  1053.7529 ,  1039.8606 ,  1027.6493 ,  1017.057  ,
        1007.9573 ,  1000.1681 ,   993.52384,   987.85759,   983.05871,
         979.01311,   975.60886,   972.74937,   970.3489 ,   968.33983,
         966.66146,   965.26054,   964.09216,   963.11826,   962.30738,
         961.63264,   961.0714 ,   960.60473,   960.2169 ,   959.89468,
         959.627  ,   959.40467,   959.22003,   959.06677,   958.93955,
         958.83394,   958.74627,   958.6735 ,   958.61313,   958.56303,
         958.52146,   958.48696,   958.45833,   958.43458,   958.41487,
         958.39852,   958.38496,   958.3737 ,   958.36437,   958.35662])

    # ======================= Kevin's PT profile ========================== #

    # WASP-43b Kevin
    plt.xlim(1e-21, 1)

    # Pressure limits
    plt.ylim(max(pres), min(pres))

    # ================ inset plot with PT profile ========================== #

    # WASP-43 Kevin all plots
    b = plt.axes([.21, .25, .10, .19], axisbg='w')

    # WASP-43b all metallicities
    plt.semilogy(temp, pres, color='r', linestyle='-', linewidth=1)
    plt.xlim(700, 2000)    
    plt.xticks(np.arange(1000, 2001, 500))   

    plt.xlabel('T [K]'  , fontsize=8)
    plt.ylabel('P [bar]', fontsize=8)
    plt.yticks(fontsize=8) 
    plt.xticks(fontsize=8)
    plt.ylim(max(pres), min(pres))

    # ================ inset plot with PT profile ========================== #

    # Save plots
    plot_out1 = plots_dir + filename.split("/")[-1][:-4] +'.ps'
    plot_out = plots_dir  + filename.split("/")[-1][:-4] +'.png'
    plt.savefig(plot_out1, dpi=300) 
    plt.savefig(plot_out , dpi=300)
    plt.close()

    return plot_out

# Call the function to execute
if __name__ == '__main__':
    # Make plot and retrieve plot's name
    plot_out = plotTEA()
    
    # Open plot
    plot = Image.open(plot_out)
    plot.show()
