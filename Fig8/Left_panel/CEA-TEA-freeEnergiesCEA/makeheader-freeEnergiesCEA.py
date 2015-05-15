
# ******************************* START LICENSE *******************************
# Thermal Equilibrium Abundances (TEA), a code to calculate gaseous molecular
# abundances under thermochemical equilibrium conditions.
# 
# This project was completed with the support of the NASA Earth and
# Space Science Fellowship Program, grant NNX12AL83H, held by Jasmina
# Blecic, Principal Investigator Joseph Harrington. Project developers
# included graduate student Jasmina Blecic and undergraduate M. Oliver
# Bowman.
# 
# Copyright (C) 2014 University of Central Florida.  All rights reserved.
# 
# This is a test version only, and may not be redistributed to any third
# party.  Please refer such requests to us.  This program is distributed
# in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.
# 
# Our intent is to release this software under an open-source,
# reproducible-research license, once the code is mature and the first
# research paper describing the code has been accepted for publication
# in a peer-reviewed journal.  We are committed to development in the
# open, and have posted this code on github.com so that others can test
# it and give us feedback.  However, until its first publication and
# first stable release, we do not permit others to redistribute the code
# in either original or modified form, nor to publish work based in
# whole or in part on the output of this code.  By downloading, running,
# or modifying this code, you agree to these conditions.  We do
# encourage sharing any modifications with us and discussing them
# openly.
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
# 4111 Libra Drive
# Orlando, FL 32816-2385
# USA
# 
# Thank you for testing TEA!
# ******************************* END LICENSE *******************************

from readconf import *

import numpy as np
import re
import os

from scipy.interpolate import UnivariateSpline

# =============================================================================
# This module contains functions to write headers for a single T-P and multiple
# T-P runs. It consists of two main functions, make_singleheader() and 
# make_atmheader() called by the runsingle.py and runatm.py modules
# respectively. The header_setup(), atm_headarr(), single_headarr(), and
# write_header() are the supporting functions for the main functions. 
# Imported by runatm.py and runsingle.py to create header files.
# =============================================================================

# Correct directory names
if location_TEA[-1] != '/':
    location_TEA += '/'

if location_out[-1] != '/':
    location_out += '/'

def header_setup(temp, pressure, spec_list,                      \
                 thermo_dir, stoich_file = 'lib/stoich.txt'):
    '''
    This function is a common setup for both single T-P and multiple T-P runs. 
    Given the thermochemical data and stoichiometric table, this function 
    returns stoichiometric values and an array of chemical potentials for the 
    species of interest at the current temperature and pressure.  It also 
    returns an array of booleans that marks which information should be read 
    from stoich_file for the current species. It is executed by the 
    make_atmheader() and make_singleheader() functions.

    Parameters
    ----------
    temp: float
         Current temperature value.
    pressure: float
         Current pressure value.
    spec_list: string array
         Array containing names of molecular species.
    thermo_dir = 'lib/gdata':  string
         Name of directory containing thermodynamic data.
    stoich_file = 'lib/stoich.txt': string
         Name of file containing stoichiometric data.

    Returns
    -------
    stoich_data: array
         Full stoichiometric information from stoich_file for species used.
    spec_stoich: float array
         Array containing values from stoich_file that correspond to the 
         species used.
    g_RT: float array
         Array containing species' chemical potentials.
    is_used: boolean array
         Array containing booleans to trim stoichiometric data to only the 
         species of interest.
    '''

    # Ensure that inputs are floats
    temp     = np.float(temp)
    pressure = np.float(pressure)

    # Obtain thermo_dir files, and count both files and species
    gdata_files = os.listdir(thermo_dir)
    n_gdata     = np.size(gdata_files)
    n_spec      = np.size(spec_list)

    # Create index of where species listed in the input file are in thermo_dir
    spec_ind = np.zeros(n_spec)
    for i in np.arange(n_spec):
        spec_file = spec_list[i] + '.txt'
        if spec_file in gdata_files:
            spec_ind[i] = gdata_files.index(spec_file)

    # Create array of G/RT
    g_RT = np.zeros(n_spec)
    R = 8.3144621 # J/K/mol
    for i in np.arange(n_spec):
        spec_file = thermo_dir + '/' + gdata_files[np.int(spec_ind[i])]
        f = open(spec_file, 'r')
        data = []
        headerline = True
        for line in f.readlines():
            if headerline:
                headerline = False
            else:
                l = [np.float(value) for value in line.split()]
                data.append(l)
        
        f.close()
       
        # Convert data to an array
        data = np.asarray(data)

        # Equation for g_RT term equation (11) in TEA Document
        #  G        G-H(298)      delta-f H(298)
        # ---  =    -------  +    -------------
        # R*T         R*T              R*T

        # First term is divided by T in JANAF   
        # Second term needs to be converted to Joules (JANAF gives kJ)
        #  G     G-H(298)             1                                    1000
        # ---  = ------- [J/K/mol] * ---         + delta-f H(298) [kJ/mol] ------ 
        # R*T      T                  R [J/K/mol]                          R*T [J/K/mol][K] 

        # Spline interpolation of JANAF term1 values
        spline_term1 = UnivariateSpline(data[:, 0], data[:,1], s=1)
        gdata_term1  = spline_term1(temp)
       
        # Take term2 from gdata/ tables
        for j in np.arange(len(data[:,0])):
            if data[:,0][j] == 298.15:
                gdata_term2  = data[:,2][j]   

        # Calculate the above equation   
        g_RT[i] = - (gdata_term1 / R) + (gdata_term2 * 1000 / (temp * R))
        
    # Get number of elements that occur in species of interest
    nostate = np.copy(spec_list)
    for i in np.arange(n_spec):
        nostate[i] = re.search('(.*?)_', spec_list[i]).group(1)

    # Location of the stoich_file
    stoich_file = location_TEA + stoich_file

    # Get stoichiometric information for species of interest
    f = open(stoich_file, 'r')
    stoich_data = []
    bline = True
    for line in f.readlines():
        l = [value for value in line.split()]
        stoich_data.append(l)

    # Store information in stoich_data array
    stoich_data = np.asarray(stoich_data)
    f.close()
    
    # Count total number of elements in stoich_file
    n_ele = np.size(stoich_data[0, 1:])
    
    # Allocate array to store current element and species stoichiometric values
    spec_stoich = np.empty((n_spec+1, n_ele), dtype=np.float)
    spec_stoich[0] = stoich_data[0,1:]

    # Place species' stoichiometric data into array
    for i in np.arange(n_spec):
        idx = np.where(stoich_data[:, 0] == nostate[i])[0]
        if np.size(idx) != 1:
            idx = idx[0]
        spec_stoich[i+1] = stoich_data[idx,1:]
    
    # Determine which elements are used to trim down final stoichiometric table
    is_used = np.empty(n_ele, dtype=np.bool)
    for j in np.arange(n_ele):
        if np.sum(spec_stoich[1:, j]) != 0:
            is_used[j] = True
        else:
            is_used[j] = False

    return stoich_data, spec_stoich, g_RT, is_used


def single_headarr(spec_list, stoich_data, spec_stoich, is_used):
    '''
    This function gathers data needed for TEA to run in a single T-P case. 
    These are: elemental abundances, species names, and their stoichiometric
    values. For the list of elements and species used, it takes the abundances
    and stoichiometric values and puts them in the final array. This function
    is run by make_singleheader() and is dependent on results from 
    header_setup().

    Parameters
    ----------
    spec_list: string array
         Array containing names of molecular species.
    stoich_data: array
         Full stoichiometric information from header_setup() for species of 
         interest.
    spec_stoich: float array
         Array containing values from header_setup() for species of interest.
    is_used: boolean array
         Array containing booleans to trim stoichiometric data to only the 
         species of interest.

    Returns
    -------
    stoich_arr: array
         Array containing elemental abundances, species names, and their 
         stoichiometric values.
    '''
    
    # Get number of species used
    n_spec = np.size(spec_list)
    
    # Allocate final header array
    stoich_arr = np.empty((n_spec + 2, np.sum(is_used) + 1), dtype=np.object)
    
    # First row is 'b' values (elemental abundances)
    stoich_arr[0,0] = 'b'
    stoich_arr[0,1:] = stoich_data[0, np.where(is_used)[0] + 1]
    
    # Second row is list of species
    stoich_arr[1,0] = 'Species'
    stoich_arr[1,1:] = stoich_data[1, np.where(is_used)[0] + 1]
    
    # Each row after contains the weights of each element referred to as
    #      stoichiometric coefficients
    for i in np.arange(n_spec):
        stoich_arr[i+2, 0] = spec_list[i]
        stoich_arr[i+2, 1:] = map(int,spec_stoich[i+1, np.where(is_used)[0]])
    
    # Convert dex abundances into number densities
    finalstoich_conv = np.empty((n_spec + 2, np.sum(is_used) + 1), \
                                                 dtype=np.object)
    finalstoich_conv[0,0] = 'b'
    finalstoich_conv[0,1:] = map(float, stoich_arr[0,1:])

    # Number densities for elements in the system are equal to 10**(dex)
    finalstoich_conv[0,1:] = 10**(finalstoich_conv[0,1:])
    
    # Elemental abundance is equal to elemental number density divided by
    #           total sum of all elemental number densities in the system
    finalstoich_conv[0,1:] /= sum(finalstoich_conv[0,1:])
    
    # Place converted values back into final header array
    stoich_arr[0] = finalstoich_conv[0]
    
    return stoich_arr


def atm_headarr(spec_list, stoich_data, spec_stoich, atom_arr, q, is_used):
    '''
    This function gathers data needed for TEA to run in a multiple T-P case. 
    These are: elemental abundances, species names, and their stoichiometric
    values. For the list of elements and species used, it takes the abundances
    and stoichiometric values and puts them in the final array. This function
    is run by make_atmheader() and is dependent on results from header_setup().

    Parameters
    ----------
    spec_list: string array
         Array containing names of molecular species.
    stoich_data: array
         Full stoichiometric information from header_setup() for species of 
         interest.
    spec_stoich: float array
         Array containing values from header_setup() for species of interest.
    atom_arr: string array
         Array containing elemental symbols and abundances.
    q: integer
         Current line number in pre-atm file.
    is_used: boolean array
         Array containing booleans to trim stoichiometric data to only the 
         species of interest.

    Returns
    -------
    stoich_arr: array
         Array containing elemental abundances, species names, and their 
         stoichiometric values.
    '''

    # Get number of species and elements used
    n_spec = np.size(spec_list)
    n_atom = np.size(atom_arr[0])
    
    # Allocate final abundance array
    stoich_arr = np.empty((n_spec + 2, np.sum(is_used) + 1), dtype=np.object)
    stoich_arr[0,0] = 'b'
    
    # Get only the abundances used in the species list
    for n in np.arange(np.sum(is_used)):
        # Get list of used elements from species list
        cur_ele = stoich_data[1, np.where(is_used)[0][n] + 1]
        
        # Place used elemental abundances into final abundance array
        for m in np.arange(n_atom):
            if atom_arr[0][m] == cur_ele:
                cur_abn = atom_arr[q][m]
        stoich_arr[0,1+n] = cur_abn

    # Fill in final abundance array
    stoich_arr[1,0] = 'Species'
    stoich_arr[1,1:] = stoich_data[1, np.where(is_used)[0] + 1]
    for i in np.arange(n_spec):
        stoich_arr[i+2, 0] = spec_list[i]
        stoich_arr[i+2, 1:] = map(int,spec_stoich[i+1, np.where(is_used)[0]])
    
    return stoich_arr


def write_header(desc, pressure, temp, stoich_arr, n_spec, g_RT):
    '''
    This function writes a header file that contains all necessary data
    for TEA to run. It is run by make_atmheader() and make_singleheader().

    Parameters
    ----------
    desc: string
         Name of output directory given by user. 
    pressure: float
         Current pressure value.
    temp: float
         Current temperature value.
    stoich_arr: array
         Array containing elemental abundances, species names, and their 
         stoichiometric values.
    n_spec: integer
         Number of species. 
    g_RT: float array
         Array containing chemical potentials for each species at the 
         current T-P.

    Returns
    -------
    None
    '''

    # Comment at top of header file
    header_comment = ("# This is a header file for one T-P. It contains the following data:\n"
    "# pressure (bar), temperature (K), elemental abundances (b, unitless),\n"
    "# species names, stoichiometric values of each element in the species (a),\n"
    "# and chemical potentials.\n\n")

    # Make header directory if it does not exist to store header files
    if not os.path.exists(location_out + desc + '/headers/'): os.makedirs(location_out + desc + '/headers/')
    #if not os.path.exists(location_out + 'headers/' + desc + '/'): os.makedirs(location_out + '/headers/' + desc + '/')

    # Create header file to be used by the main pipeline
    outfile = location_out + desc + '/headers/' + 'header_' + desc + ".txt"

    # Open file to write
    f = open(outfile, 'w+')

    # Write comments and data
    f.write(header_comment)
    f.write(np.str(pressure) + '\n')
    f.write(np.str(temp)     + '\n')

    b_line = stoich_arr[0][0]
    for i in np.arange(np.shape(stoich_arr)[1] - 1):
        b_line += ' ' + np.str(stoich_arr[0][i + 1])
    f.write(b_line + '\n')

    # Retrieve the width of the stoichiometric array
    width = np.shape(stoich_arr)[1]


    print g_RT
    for i in np.arange(len(g_RT)):
        g_RT[i]=str(g_RT[i])
    # Add title for list of chemical potentials
    g_RT = np.append(["Chemical potential"], g_RT)

    # Add title for species names
    stoich_arr[1][0] = "# Species"

    # Loop over all species and titles
    for i in np.arange(n_spec + 1):
        # Loop over species and number of elements
        for j in np.arange(width):
            # Write species names
            if j == 0:
                f.write(np.str(stoich_arr[i+1][j]).rjust(9) + "  ")
            # Write elemental number density
            else:
                f.write(np.str(stoich_arr[i+1][j]).ljust(3))
        # Write chemical potentials, adjust spacing for minus sign
        if g_RT[i][0] == '-':
            f.write(np.str(g_RT[i]).rjust(13) + '\n')
        else:
            f.write(np.str(g_RT[i]).rjust(14) + '\n')

    f.close()


def make_singleheader(infile, desc, thermo_dir):
    '''
    This is the main function that creates single-run TEA header. It reads
    the input T-P file and retrieves necessary data. It calls the
    header_setup(), single_headarr(), and write_header() functions to create a
    header for the single T-P point. This function is called by the 
    runsingle.py module. 

    Parameters
    ----------
    infile: string
         Name of the input file (single T-P file).         
    desc: string
         Name of output directory given by user.    
    thermo_dir = 'lib/gdata':  string
         Name of directory containing thermodynamic data.

    Returns
    -------
    desc: string
         Name of output directory given by user. 
    pressure: float
         Current pressure value.
    temp: float
         Current temperature value.
    stoich_arr: array
         Array containing elemental abundances, species names, and their 
         stoichiometric values.
    n_spec: integer
         Number of species. 
    g_RT: float array
         Array containing chemical potentials for each species at the 
         current T-P.
    '''

    # Read single input file
    f = open(infile, 'r')

    # Allocate species list
    spec_list = []
    
    # Begin by reading first line of file (l = 0)
    l = 0
    
    # Loop over all lines in input file to retrieve appropriate data
    for line in f.readlines():
        # Retrieve temperature
        if (l == 0):
            temp = np.float([value for value in line.split()][0])
        # Retrieve pressure
        if (l == 1):
            pressure = np.float([value for value in line.split()][0])
        # Retrieve list of species
        if (l > 1):
            val = [value for value in line.split()][0]
            spec_list = np.append(spec_list, val)
        # Update line number
        l += 1
    
    f.close()
    
    # Retrieve number of species used
    n_spec = np.size(spec_list)
    
    # Execute header setup
    stoich_data, spec_stoich, g_RT, is_used =                              \
                        header_setup(temp, pressure, spec_list, thermo_dir)
    
    # Execute single-specific header setup
    stoich_arr = single_headarr(spec_list, stoich_data, spec_stoich, is_used)
    
    # Write header array to file
    write_header(desc, pressure, temp, stoich_arr, n_spec, g_RT)


def make_atmheader(q, spec_list, pressure, temp, atom_arr, desc, \
                   thermo_dir = 'lib/gdata'):
    '''
    This is the main function that creates a header for one T-P of a 
    pre-atm file. It retrieves number of elements and species used for 
    only the q-th T-P point in the pre-atm file. It then calls the 
    header_setup(), atm_headarr(), and write_header() functions to create a 
    header for this point. This function is called by the runatm.py module. 

    Parameters
    ----------
    q: integer
         Current line number in pre-atm file.
    spec_list: string array
         Array containing names of molecular species.
    pressure: float
         Current pressure value.
    temp float
         Current temperature value.
    atom_arr: string array
         Array containing elemental symbols and abundances.
    desc: string
         Name of output directory given by user.    
    thermo_dir = 'lib/gdata':  string
         Name of directory containing thermodynamic data.

    Returns
    -------
    None
    '''

    # Retrieve number of elements and species used
    n_spec   = np.size(spec_list)
    n_atom   = np.size(atom_arr[0])
    
    # Execute header setup
    stoich_data, spec_stoich, g_RT, is_used =                             \
                       header_setup(temp, pressure, spec_list, thermo_dir)
    
    # Execute multiple-specific header setup
    stoich_arr = atm_headarr(spec_list, stoich_data, spec_stoich, atom_arr,\
                                                               q, is_used)

    ########################################## JB added for CEA-TEA comparison, calculated CEA free energies based on thermo.inp and CEA theory doc

    # C H O N H2 CO CH4 H2O N2 NH3 used both CEA temperature ranges 100 - 1000 and 1000 - 6000
    # values produced with TEA-CEA-thermo.py
    g_RT_CEA = \
    np.array([[ 246.17083749, 840.71265488, 547.8395475, 277.90072018, -18.67806007, -159.83085247, -115.71539877, -317.11219575, -26.16486495, -82.07110085], \
         [ 117.06942949,  411.73699577,  265.58349055,  130.22870673,  -16.03051564, -90.56455389,  -67.65797087, -168.5031827,   -23.36551338,  -51.19461686], \
         [  73.5995867,   268.30588659,  171.06209811,   80.52458487,  -15.71729001, -68.08711861,  -52.32279721, -119.66025955,  -23.04528681,  -41.60248807], \
         [  51.65229151,  196.37737345,  123.5890281,    55.44806008,  -15.85408413, -57.14625053,  -45.02087002,  -95.58221113,  -23.1827982,   -37.17045471], \
         [  38.35757863,  153.09374437,   94.97885032,   40.27140132,  -16.113398,   -50.76011218,  -40.88570144,  -81.34357282,  -23.44312835,  -34.7456328 ], \
         [  29.41053873,  124.15401684,   75.82150049,   30.0677758,   -16.40434374, -46.62293047,  -38.31615827,  -71.99322719,  -23.73608877,  -33.29809417], \
         [  22.95999346,  103.42294291,   62.0778767,    22.71867616,  -16.69651887, -43.75533653,  -36.63181566,  -65.41883956,  -24.03195154,  -32.39396037], \
         [  18.07728855,   87.82982105,   51.72536289,   17.16149948,  -16.97910706, -41.67177719,  -35.49436123,  -60.56880984,  -24.32009708,  -31.81969987], \
         [  14.24481596,   75.66701203,   43.6385945,    12.80410298,  -17.24852198, -40.10462491,  -34.71664127,  -56.86138868,  -24.59682952,  -31.4585344 ], \
         [ 11.15100266,  65.9089218,   37.14134454,   9.29013704, -17.50405995, -38.89449725, -34.18678245, -53.94898416, -24.86116642, -31.24152113], \
         [  8.59693558,  57.90225699,  31.80264907,   6.39216742, -17.74620822, -37.94071192, -33.83383152, -51.61130364, -25.11323503, -31.12549977], \
         [  6.44958077,  51.21105962,  27.33477597,   3.95813229, -17.97588694, -37.17663579, -33.61065069, -49.70206249, -25.3536106, -31.08216729], \
         [  4.616544,    45.53322305,  23.53822804,   1.8824519,  -18.19414182, -36.5564762,  -33.48474839, -48.12037511, -25.58302726, -31.09227234], \
         [  3.03161979,  40.65274791,  20.27029589,   0.0894962,  -18.40199883, -36.04777144, -33.43303356, -46.79443971, -25.8022518,  -31.1423038 ], \
         [  1.6461039,   36.41107921,  17.42617253,  -1.4763533,  -18.60040403, -35.62690248, -33.4386786,  -45.67178209, -26.01202766, -31.22253014], \
         [  0.4233531,   32.68918307,  14.92713832,  -2.85692712, -18.79020655, -35.27629861, -33.489162,   -44.71317823, -26.21305155, -31.32579042], \
         [ -0.66474184,  29.39594313,  12.71290865,  -4.08430225, -18.97215951, -34.9826362, -33.57500118, -43.88873285, -26.40596533, -31.44672   ], \
         [ -1.64011199,  26.46042065,  10.73652763,  -5.18349662, -19.14692792, -34.73564246, -33.68890665, -43.1752718,  -26.59135523, -31.58123873], \
         [ -2.52012531,  23.82655565,   8.96087117,  -6.17431587, -19.31509864, -34.52728011, -33.82520226, -42.5545618,  -26.76975448, -31.72620325], \
         [ -3.31871928,  21.44945576,   7.35619689,  -7.07264834, -19.47719042, -34.35117886, -33.97941856, -42.01206695, -26.94164725, -31.87916521], \
         [ -4.04721155,  19.29274364,   5.89839257,  -7.89139154, -19.63366302, -34.20223103, -34.148002,   -41.53606214, -27.10747318, -32.03819962], \
         [ -4.71489055,  17.32662682,   4.56770118,  -8.64112716, -19.78492532, -34.07629865, -34.32810326, -41.11698949, -27.26763173, -32.2017809 ], \
         [ -5.32945292,  15.52647012,   3.34777769,  -9.33062071, -19.93134226, -33.96999782, -34.51742065, -40.74698351, -27.42248635, -32.3686917 ], \
         [ -5.89733218,  13.87172439,   2.22498108,  -9.96719694, -20.07324082, -33.88053783, -34.71408283, -40.41951556, -27.57236828, -32.53795517], \
         [ -6.42394896,  12.34511174,   1.18783584, -10.55702561, -20.2109150,  -33.8055993,  -34.9165597,  -40.12912409, -27.71757996, -32.70878374], \
         [ -6.91390395,  10.93199857,   0.22661734, -11.10534164, -20.34463035, -33.74324099, -35.12359403, -39.87120748, -27.85839804, -32.88054008], \
         [ -7.3711281,    9.61990742,  -0.66697085, -11.61661669, -20.47462715, -33.69182758, -35.33414843, -39.6418632,  -27.9950761,  -33.05270703], \
         [ -7.79900102,   8.39813342,  -1.5000508 , -12.0946941,  -20.60112399, -33.64997328, -35.54736383, -39.43776158, -28.12784702, -33.22486413], \
         [ -8.20044493,   7.25743967,  -2.27877163, -12.54289616, -20.72432014, -33.61649725, -35.76252669, -39.25604592, -28.2569251,  -33.39666924], \
         [ -8.57800005,   6.18981358,  -3.0084708,  -12.96411005, -20.84439788, -33.59038819, -35.97904295, -39.09425268, -28.38250791, -33.56784401]])


    # Divide temperature range
    T = np.linspace(100.0, 3000.0, 30)

    for i in np.arange(len(T)):
        if float(temp) == T[i]:
            g_RT = g_RT_CEA[i]

    # Write header array to file
    write_header(desc, pressure, temp, stoich_arr, n_spec, g_RT)

