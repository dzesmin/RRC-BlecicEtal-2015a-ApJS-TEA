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
 
    # WASP-12b all plots (please DO NOT USE THE SAME COLORS FOR YOUR PAPER!)
    colors = 'b', '#FF1CAE','#FF0000' , '#FFAA00', '#00FFFF', '#00FF00', '#91219E', '#BCEE68' , 'g', '#ffc3a0', 'c','m'
    color_index = 0

    # Plot all specs with different colours and labels
    for i in np.arange(nspec):
        # addition for WASP-43b species that does not change with metalicity
        if i==0 or i==1 or i==8 or i==9 or i==10 or i==11: 
            plt.loglog(data[i+1], data[0], '--', color=colors[color_index], \
                                        linewidth=1, label=str(spec[i]))
        else:
            plt.loglog(data[i+1], data[0], '-', color=colors[color_index], \
                                        linewidth=2, label=str(spec[i]))
        color_index += 1

    # Label the plot
    plt.xlabel('Mixing Fraction'    , fontsize=14)
    plt.ylabel('Pressure [bar]'    , fontsize=14)

    # Font-size of x and y ticks
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    # WASP-12b solar
    plt.text(8e-4, 1e-3, 'CO', color='#FF1CAE', fontsize=14)
    plt.text(6e-10,3e-2, 'CO$_{2}$', color='#FF0000', fontsize=14)
    plt.text(1e-16,1e-3, 'CH$_{4}$', color='#FFAA00', fontsize=14)
    plt.text(8e-4,0.2, 'H$_{2}$O', color='#00FFFF', fontsize=14)
    plt.text(2e-13,5e-3, 'HCN', color='#00FF00', fontsize=14)
    plt.text(8e-23,1.5e-4, 'C$_{2}$H$_{2}$', color='#91219E', fontsize=14)
    plt.text(2e-23,7e-3, 'C$_{2}$H$_{4}$', color='#BCEE68', fontsize=14)
    plt.text(1e-1, 1e-2, 'H', color='b', fontsize=14)
    plt.text(4e-6,1e-4, 'N$_{2}$', color='g', fontsize=14)
    plt.text(7e-13,3e-5, 'NH$_{3}$', color='#ffc3a0', fontsize=14)
    plt.text(1e-6, 1e-3, 'HS', color='c', fontsize=14)
    plt.text(4e-8,7e-1, 'H$_{2}$S', color='m', fontsize=14)

    # WASP-12b C/O=1.2
    #plt.text(8e-4, 10, 'CO', color='#FF1CAE', fontsize=14)
    #plt.text(6e-17,1e-4, 'CO$_{2}$', color='#FF0000', fontsize=14)
    #plt.text(2e-9,6e-5, 'CH$_{4}$', color='#FFAA00', fontsize=14)
    #plt.text(1e-9,6e-3, 'H$_{2}$O', color='#00FFFF', fontsize=14)
    #plt.text(5e-7,3e-2, 'HCN', color='#00FF00', fontsize=13)
    #plt.text(5e-5,7e-2, 'C$_{2}$H$_{2}$', color='#91219E', fontsize=14)
    #plt.text(3e-11, 1, 'C$_{2}$H$_{4}$', color='#BCEE68', fontsize=14)
    #plt.text(2e-2, 1e-4, 'H', color='b', fontsize=14)
    #plt.text(5e-5,8e-3, 'N$_{2}$', color='g', fontsize=14)
    #plt.text(5e-12,2e-3, 'NH$_{3}$', color='#ffc3a0', fontsize=14)
    #plt.text(4.5e-7, 9e-2, 'HS', color='c', fontsize=13)
    #plt.text(2e-7,8e-5, 'H$_{2}$S', color='m', fontsize=13)

    # ======================= Kevin's PT profiles ========================== #

    # Kevin's PT profile from tp_crp.dat and tp_orp.dat
    # WASP12b, Stevenson et al 2014

    pres = np.array([  1.00000000e+02,   8.49753000e+01,   7.22081000e+01,
         6.13591000e+01,   5.21401000e+01,   4.43062000e+01,
         3.76494000e+01,   3.19927000e+01,   2.71859000e+01,
         2.31013000e+01,   1.96304000e+01,   1.66810000e+01,
         1.41747000e+01,   1.20450000e+01,   1.02353000e+01,
         8.69749000e+00,   7.39072000e+00,   6.28029000e+00,
         5.33670000e+00,   4.53488000e+00,   3.85353000e+00,
         3.27455000e+00,   2.78256000e+00,   2.36449000e+00,
         2.00923000e+00,   1.70735000e+00,   1.45083000e+00,
         1.23285000e+00,   1.04762000e+00,   8.90215000e-01,
         7.56463000e-01,   6.42807000e-01,   5.46228000e-01,
         4.64159000e-01,   3.94421000e-01,   3.35160000e-01,
         2.84804000e-01,   2.42013000e-01,   2.05651000e-01,
         1.74753000e-01,   1.48497000e-01,   1.26186000e-01,
         1.07227000e-01,   9.11163000e-02,   7.74264000e-02,
         6.57933000e-02,   5.59081000e-02,   4.75081000e-02,
         4.03702000e-02,   3.43047000e-02,   2.91505000e-02,
         2.47708000e-02,   2.10490000e-02,   1.78865000e-02,
         1.51991000e-02,   1.29155000e-02,   1.09750000e-02,
         9.32603000e-03,   7.92483000e-03,   6.73415000e-03,
         5.72237000e-03,   4.86260000e-03,   4.13201000e-03,
         3.51119000e-03,   2.98365000e-03,   2.53536000e-03,
         2.15443000e-03,   1.83074000e-03,   1.55568000e-03,
         1.32194000e-03,   1.12332000e-03,   9.54548000e-04,
         8.11131000e-04,   6.89261000e-04,   5.85702000e-04,
         4.97702000e-04,   4.22924000e-04,   3.59381000e-04,
         3.05386000e-04,   2.59502000e-04,   2.20513000e-04,
         1.87382000e-04,   1.59228000e-04,   1.35305000e-04,
         1.14976000e-04,   9.77010000e-05,   8.30218000e-05,
         7.05480000e-05,   5.99484000e-05,   5.09414000e-05,
         4.32876000e-05,   3.67838000e-05,   3.12572000e-05,
         2.65609000e-05,   2.25702000e-05,   1.91791000e-05,
         1.62975000e-05,   1.38489000e-05,   1.17681000e-05,
         1.00000000e-05])


    # Kevin tp_orp.dat solar
    # WASP12b, Stevenson et al 2014
    temp = np.array([ 2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,
        2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,
        2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,
        2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,
        2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,
        2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,
        2948.3 ,  2948.3 ,  2948.3 ,  2948.3 ,  2948.26,  2947.68,
        2946.6 ,  2945.04,  2943.03,  2940.59,  2937.75,  2934.53,
        2930.96,  2927.08,  2922.89,  2918.48,  2914.35,  2910.51,
        2906.96,  2903.69,  2900.71,  2898.01,  2895.59,  2891.16,
        2883.82,  2873.64,  2860.69,  2845.05,  2826.77,  2805.92,
        2782.58,  2756.81,  2728.67,  2698.24,  2667.89,  2638.57,
        2610.28,  2583.01,  2556.78,  2531.57,  2507.39,  2484.25,
        2462.13,  2441.04,  2420.97,  2401.94,  2383.94,  2366.96,
        2351.02,  2336.1 ,  2322.21,  2309.35,  2297.52,  2286.71,
        2276.94,  2268.2 ,  2260.48,  2253.79,  2248.13,  2238.36,
        2234.76,  2232.19,  2230.64,  2230.13])

    '''
    # Kevin tp_crp.dat C/O=1.2
    # WASP12b, Stevenson et al 2014
    temp = np.array([ 2964.76,  2964.76,  2964.76,  2964.76,  2964.76,  2964.76,
        2964.76,  2964.76,  2964.76,  2964.76,  2964.76,  2964.76,
        2964.76,  2964.76,  2964.76,  2964.76,  2964.76,  2964.76,
        2964.76,  2964.76,  2964.76,  2964.76,  2964.76,  2964.76,
        2964.76,  2964.76,  2964.76,  2964.76,  2964.76,  2964.76,
        2964.76,  2964.76,  2964.76,  2964.76,  2964.76,  2964.76,
        2964.76,  2964.51,  2963.85,  2962.85,  2961.53,  2959.96,
        2957.23,  2949.72,  2937.53,  2920.75,  2899.48,  2873.8 ,
        2844.07,  2810.51,  2773.19,  2732.13,  2687.4 ,  2639.98,
        2593.58,  2548.2 ,  2503.84,  2460.5 ,  2418.18,  2376.87,
        2336.59,  2297.33,  2259.09,  2221.86,  2185.66,  2150.48,
        2116.31,  2083.17,  2051.05,  2019.94,  1989.86,  1960.79,
        1932.75,  1905.72,  1879.72,  1854.73,  1830.77,  1807.82,
        1785.9 ,  1764.99,  1745.1 ,  1726.24,  1708.39,  1691.56,
        1675.76,  1660.97,  1647.2 ,  1634.46,  1622.73,  1612.02,
        1602.33,  1593.66,  1586.01,  1579.39,  1573.78,  1564.09,
        1560.52,  1557.97,  1556.44,  1555.93])
    '''

    # ======================= Kevin's PT profiles ========================== #

    # WASP-12b
    plt.xlim(1e-24, 1)

    # Pressure limits
    plt.ylim(max(pres), min(pres))

    # ================ inset plot with PT profile ========================== #

    # WASP12b all plots
    b = plt.axes([.21, .22, .10, .19], axisbg='w')

    # WASP12b solar
    plt.semilogy(temp, pres, color='k', linestyle='-', linewidth=2)
    plt.xlim(2000, 3100)   # WASP-12b
    plt.xticks(np.arange(2000, 3100, 500))

    # WASP12b C/O=1.2
    #plt.semilogy(temp, pres, color='orange', linestyle='--', linewidth=3)
    #plt.xlim(1400, 3100)   # WASP-12b
    #plt.xticks(np.arange(1500, 3100, 1500))  # WASP-12b

    plt.xlabel('T [K]'     , fontsize=8)
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
