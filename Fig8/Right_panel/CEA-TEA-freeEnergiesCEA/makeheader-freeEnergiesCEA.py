
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
    np.array([[ 12.36414503, -16.36954367,  69.72808408,  39.68537711,  10.66759989, \
   13.27864516, -17.39932449, -39.3624422 , -77.5876507 , -34.38252737, \
  -55.08134714,  -9.91622307, -23.19291862, -24.7526095 , -31.31756266,\
  -29.46973513],\
 [ 12.36414503, -16.36954367,  69.72808408,  39.68537711,  10.66759989,\
   13.27864516, -17.39932449, -39.3624422 , -77.5876507 , -34.38252737,\
  -55.08134714,  -9.91622307, -23.19291862, -24.7526095 , -31.31756266,\
  -29.46973513],\
 [ 12.36384159, -16.36956164,  69.72712761,  39.68474017,  10.66725528,\
   13.27826285, -17.39934993, -39.36232405, -77.58717619, -34.38247637,\
  -55.08106228,  -9.91643528, -23.19304865, -24.75263582, -31.31754211,\
  -29.46974243],\
 [ 12.36353816, -16.36957961,  69.72617116,  39.68410325,  10.66691068,\
   13.27788055, -17.39937537, -39.36220591, -77.5867017 , -34.38242538,\
  -55.08077744,  -9.91664748, -23.19317868, -24.75266215, -31.31752156,\
  -29.46974973],\
 [ 12.36293131, -16.36961555,  69.72425832,  39.68282943,  10.66622149,\
   13.27711598, -17.39942624, -39.36196963, -77.58575273, -34.38232339,\
  -55.08020776,  -9.91707188, -23.19343873, -24.75271481, -31.31748047,\
  -29.46976434],\
 [ 12.3626279,  -16.36963352,  69.72330193,  39.68219255,  10.66587691,\
   13.2767337,  -17.39945167, -39.36185149, -77.58527827, -34.3822724,\
  -55.07992293,  -9.91728407, -23.19356876, -24.75274114, -31.31745993,\
  -29.46977164],\
 [ 12.36202109, -16.36966946,  69.72138921,  39.68091881,  10.66518776,\
   13.27596917, -17.39950254, -39.36161523, -77.58432937, -34.38217043,\
  -55.0793533 ,  -9.91770845, -23.19382881, -24.75279379, -31.31741885,\
  -29.46978625],\
 [ 12.36111093, -16.36972337,  69.71852028,  39.6790083 ,  10.6641541,\
   13.27482244, -17.39957885, -39.36126086, -77.58290611, -34.38201749,\
  -55.07849889,  -9.91834499, -23.19421886, -24.75287278, -31.31735723,\
  -29.46980816],\
 [ 12.36020081, -16.36977728,  69.71565153,  39.67709792,  10.6631205,\
   13.27367577, -17.39965516, -39.36090651, -77.58148294, -34.38186456,\
  -55.07764455,  -9.91898149, -23.19460891, -24.75295176, -31.31729562,\
  -29.46983008],\
 [ 12.35929076, -16.36983118,  69.71278295,  39.67518765,  10.66208696,\
   13.27252918, -17.39973146, -39.36055219, -77.58005986, -34.38171165,\
  -55.07679026,  -9.91961795, -23.19499894, -24.75303074, -31.31723403,\
  -29.469852  ],\
 [ 12.35807743, -16.36990306,  69.70895845,  39.6726408 ,  10.660709,\
   13.27100049, -17.39983319, -39.3600798 , -77.57816257, -34.38150781,\
  -55.0756513 ,  -9.92046652, -23.19551896, -24.75313605, -31.31715191,\
  -29.46988124],\
 [ 12.35656091, -16.3699929 ,  69.70417828,  39.66945753,  10.65898671,\
   13.2690898 , -17.39996036, -39.35948938, -77.57579121, -34.38125304,\
  -55.07422774,  -9.92152715, -23.19616895, -24.75326767, -31.3170493,\
  -29.46991779],\
 [ 12.35474128, -16.37010071,  69.69844271,  39.66563805,  10.65692019,\
   13.26679723, -17.40011296, -39.35878098, -77.57294591, -34.38094739,\
  -55.07251969,  -9.92279978, -23.19694889, -24.75342563, -31.3169262,\
  -29.46996167],\
 [ 12.35231546, -16.37024445,  69.6907964 ,  39.66054612,  10.65416522,\
   13.26374091, -17.40031641, -39.3578366 , -77.56915279, -34.38053996,\
  -55.07024264,  -9.92449641, -23.19798873, -24.75363622, -31.31676213,\
  -29.47002021],\
 [ 12.34989002, -16.37038818,  69.68315134,  39.65545503,  10.65141069,\
   13.26068507, -17.40051985, -39.3568924 , -77.56536033, -34.38013267,\
  -55.06796602,  -9.92619279, -23.19902847, -24.75384681, -31.31659814,\
  -29.47007877],\
 [ 12.34655567, -16.37058579,  69.67264144,  39.64845614,  10.64762393,\
   13.2564841 , -17.40079957, -39.35559444, -77.56014682, -34.37957284,\
  -55.06483633,  -9.92852493, -23.20045796, -24.75413635, -31.31637277,\
  -29.47015933],\
 [ 12.34261602, -16.37081933,  69.66022372,  39.64018675,  10.64314975,\
   13.25152051, -17.40113013, -39.35406093, -77.55398705, -34.37891154,\
  -55.06113862,  -9.9312805 , -23.20214712, -24.75447852, -31.3161066,\
  -29.47025462],\
 [ 12.33807154, -16.37108877,  69.6458997 ,  39.63064788,  10.63798868,\
   13.2457949 , -17.40151152, -39.3522921 , -77.54688184, -34.37814893,\
  -55.05687338,  -9.93445921, -23.20409583, -24.7548733 , -31.31579972,\
  -29.47036465],\
 [ 12.33261996, -16.37141206,  69.62871671,  39.61920508,  10.63179745,\
   13.23892646, -17.40196914, -39.35017037, -77.53835873, -34.37723438,\
  -55.05175702,  -9.93827255, -23.20643383, -24.75534701, -31.31543181,\
  -29.47049682],\
 [ 12.3256569,  -16.37182511,  69.60676992,  39.60458986,  10.62388968,\
   13.23015373, -17.40255381, -39.34746063, -77.52747307, -34.37606674,\
  -55.04522249,  -9.94314338, -23.20942058, -24.75595224, -31.31496224,\
  -29.47066592],\
 [ 12.31778946, -16.37229198,  69.581973  ,  39.58807655,  10.61495485,\
   13.22024162, -17.40321466, -39.34439928, -77.51517428, -34.37474806,\
  -55.03783974,  -9.94864713, -23.21279594, -24.75663633, -31.31443216,\
  -29.47085735],\
 [ 12.30780973, -16.37288444,  69.55051906,  39.56712998,  10.6036212,\
   13.20766831, -17.40405329, -39.34051653, -77.49957454, -34.3730763,\
  -55.0284756 ,  -9.95562899, -23.21707859, -24.75750448, -31.31376047,\
  -29.47110075],\
 [ 12.29632605, -16.37356652,  69.51432577,  39.54302714,  10.5905796,\
   13.19320028, -17.40501881, -39.33604937, -77.48162544, -34.37115393,\
  -55.01770134,  -9.96366365, -23.22200811, -24.758504  , -31.31298856,\
  -29.47138161],\
 [ 12.28213458, -16.37440995,  69.46959944,  39.51324156,  10.57446292,\
   13.17532089, -17.40621272, -39.33052993, -77.45944627, -34.36878023,\
  -55.00438813,  -9.97359376, -23.22810214, -24.75974   , -31.31203612,\
  -29.47172987],\
 [ 12.26524292, -16.37541459,  69.4163648 ,  39.47778959,  10.55527986,\
   13.15403978, -17.40763487, -39.32396187, -77.43305034, -34.36595773,\
  -54.9885441 ,  -9.98541457, -23.23535881, -24.76121233, -31.31090457,\
  -29.47214608],\
 [ 12.24475691, -16.37663407,  69.35180488,  39.43479511,  10.53201501,\
   13.12823056, -17.40936118, -39.31599845, -77.4010424 , -34.36253876,\
  -54.96933195,  -9.99975264, -23.24416417, -24.76299962, -31.30953536,\
  -29.47265331],\
 [ 12.22039058, -16.37808607,  69.27501998,  39.38365869,  10.50434369,\
   13.09753309, -17.41141669, -39.30652986, -77.36297849, -34.35847811,\
  -54.94648549, -10.01680928, -23.25464391, -24.76512786, -31.30791124,\
  -29.47326014],\
 [ 12.1912627,  -16.379824,    69.18323513,  39.32253199,  10.47126529,\
   13.06083739, -17.41387704, -39.29521555, -77.31748597, -34.35363239,\
  -54.91918122, -10.03720299, -23.26718086, -24.76767541, -31.30597609,\
  -29.47399057],\
 [ 12.1562036,  -16.38191898,  69.07276765,  39.24896195,  10.43145166,\
   13.01667022, -17.41684298, -39.28160399, -77.26274391, -34.34781216,\
  -54.88632683, -10.06175509, -23.28228417, -24.77074667, -31.30365608,\
  -29.47487701],\
 [ 12.11435925, -16.38442395,  68.94093095,  39.16115864,  10.38393323,\
   12.96395612, -17.42038949, -39.2653676 , -77.19742724, -34.34088303,\
  -54.84712787, -10.09106711, -23.30032988, -24.77441945, -31.30090022,\
  -29.47594542],\
 [ 12.06429884, -16.38742727,  68.7832232 ,  39.0561229 ,  10.32708553,\
   12.90089326, -17.42464178, -39.24595686, -77.11931445, -34.33261852,\
  -54.80025231, -10.1261462 , -23.32194658, -24.77882357, -31.29762214,\
  -29.47723856],\
 [ 12.00462302, -16.39101672,  68.5952448 ,  38.93092315,  10.25932003,\
   12.82571974, -17.42972427, -39.22283728, -77.02623906, -34.32280252,\
  -54.74440174, -10.16797987, -23.34775498, -24.78408818, -31.29374139,\
  -29.47880142],\
 [ 11.93367734, -16.39529721,  68.37179652,  38.78209463,  10.17875871,\
   12.73635273, -17.43578568, -39.19537921, -76.91564447, -34.31118386,\
  -54.67804439, -10.21773766, -23.37849373, -24.79036771, -31.28916612,\
  -29.48068978],\
 [ 11.84956974, -16.40039045,  68.10693646,  38.60567704,  10.08325413,\
   12.63041071, -17.44299862, -39.16286622, -76.78461458, -34.297482,\
  -54.59943391, -10.27676021, -23.41501496, -24.79784148, -31.2837963,\
  -29.48297138],\
 [ 11.75019319, -16.40643448,  67.79405334,  38.39726302,   9.97041517,\
   12.50524213, -17.45155897, -39.12450583, -76.62991318, -34.28139451,\
  -54.50663335, -10.34654508, -23.45827861, -24.80671319, -31.27752799,\
  -29.48572757],\
 [ 11.63267608, -16.4136186 ,  67.4241396 ,  38.15084736,   9.83698324,\
   12.35723364, -17.46173535, -39.07922043, -76.44713492, -34.26251383,\
  -54.39700656, -10.42913536, -23.50959802, -24.81726227, -31.27022299,\
  -29.4890719 ],\
 [ 11.49458965, -16.4221115 ,  66.98959629,  37.86136089,   9.68020343,\
   12.18333112, -17.47376747, -39.02611668, -76.23259183, -34.24052903,\
  -54.26835052, -10.52627473, -23.57012122, -24.82973861, -31.26178964,\
  -29.49312034],\
 [ 11.33340164, -16.43209599,  66.48251657,  37.52352656,   9.49720439,\
   11.98035201, -17.48791535, -38.96427783, -75.98246964, -34.2151435,  -54.11839,\
  -10.63979388, -23.64107556, -24.84441372, -31.25215296, -29.49801002],\
[  1.11462479e+01,  -1.64437854e+01,   6.58939724e+01,   3.71313834e+01,\
    9.28473941e+00,   1.17446983e+01,  -1.75044827e+01,  -3.88926810e+01,\
   -7.56924837e+01,  -3.41860484e+01,  -5.39445722e+01,  -1.07717757e+01,\
   -2.37238784e+01,  -2.48616052e+01,  -3.12412479e+01,  -2.95039116e+01],\
 [  1.09300120e+01,  -1.64574221e+01,   6.52142724e+01,   3.66784572e+01,\
    9.03927710e+00,   1.14724573e+01,  -1.75238148e+01,  -3.88102346e+01,\
   -7.53580167e+01,  -3.41529481e+01,  -5.37441512e+01,  -1.09245049e+01,\
   -2.38201167e+01,  -2.48816747e+01,  -3.12290337e+01,  -2.95110344e+01],\
 [  1.06814205e+01,  -1.64732747e+01,   6.44332695e+01,   3.61579641e+01,\
    8.75711013e+00,   1.11595228e+01,  -1.75462954e+01,  -3.87158224e+01,\
   -7.49742831e+01,  -3.41155910e+01,  -5.35142875e+01,  -1.11004086e+01,\
   -2.39315216e+01,  -2.49050248e+01,  -3.12155112e+01,  -2.95196325e+01],\
 [  1.03974193e+01,  -1.64916187e+01,   6.35415531e+01,   3.55636039e+01,\
    8.43478307e+00,   1.08020697e+01,  -1.75723184e+01,  -3.86084557e+01,\
   -7.45369285e+01,  -3.40738427e+01,  -5.32524101e+01,  -1.13017969e+01,\
   -2.40598189e+01,  -2.49320709e+01,  -3.12007554e+01,  -2.95300015e+01],\
 [  1.00768271e+01,  -1.65126306e+01,   6.25356439e+01,   3.48930218e+01,\
    8.07096929e+00,   1.03986357e+01,  -1.76021389e+01,  -3.84879009e+01,\
   -7.40445817e+01,  -3.40279361e+01,  -5.29577446e+01,  -1.15296936e+01,\
   -2.42059889e+01,  -2.49630851e+01,  -3.11850063e+01,  -2.95424213e+01],\
 [  9.72074395e+00,  -1.65363549e+01,   6.14192591e+01,   3.41486544e+01,\
    7.66693300e+00,   9.95063329e+00,  -1.76358256e+01,  -3.83548198e+01,\
   -7.34994540e+01,  -3.39785062e+01,  -5.26316695e+01,  -1.17835348e+01,\
   -2.43700542e+01,  -2.49981474e+01,  -3.11686711e+01,  -2.95571266e+01],\
 [  9.32543579e+00,  -1.65631798e+01,   6.01810113e+01,   3.33228591e+01,\
    7.21845683e+00,   9.45339870e+00,  -1.76739369e+01,  -3.82081146e+01,\
   -7.28964554e+01,  -3.39256120e+01,  -5.22712074e+01,  -1.20662449e+01,\
   -2.45543685e+01,  -2.50378495e+01,  -3.11520043e+01,  -2.95746055e+01],\
 [  8.89073677e+00,  -1.65932849e+01,   5.88207656e+01,   3.24154862e+01,\
    6.72537602e+00,   8.90676463e+00,  -1.77167374e+01,  -3.80480843e+01,\
   -7.22360911e+01,  -3.38699361e+01,  -5.18767428e+01,  -1.23782621e+01,\
   -2.47597809e+01,  -2.50824790e+01,  -3.11355173e+01,  -2.95952710e+01],\
 [  8.41690682e+00,  -1.66268454e+01,   5.73397776e+01,   3.14273042e+01,\
    6.18801235e+00,   8.31110567e+00,  -1.77644879e+01,  -3.78752396e+01,\
   -7.15196291e+01,  -3.38123323e+01,  -5.14491304e+01,  -1.27197696e+01,\
   -2.49870658e+01,  -2.51323223e+01,  -3.11198213e+01,  -2.96195802e+01],\
 [  7.91020180e+00,  -1.66636221e+01,   5.57580663e+01,   3.03715978e+01,\
    5.61348899e+00,   7.67433681e+00,  -1.78168615e+01,  -3.76923025e+01,\
   -7.07574535e+01,  -3.37544414e+01,  -5.09946736e+01,  -1.30866535e+01,\
   -2.52341831e+01,  -2.51870531e+01,  -3.11057661e+01,  -2.96477127e+01],\
 [  7.38113448e+00,  -1.67030349e+01,   5.41088663e+01,   2.92704840e+01,\
    5.01375139e+00,   7.00971519e+00,  -1.78730464e+01,  -3.75034653e+01,\
   -6.99662119e+01,  -3.36982742e+01,  -5.05233976e+01,  -1.34716621e+01,\
   -2.54968865e+01,  -2.52458362e+01,  -3.10942311e+01,  -2.96795407e+01],\
 [  6.82509060e+00,  -1.67456138e+01,   5.23782203e+01,   2.81145772e+01,\
    4.38359514e+00,   6.31149148e+00,  -1.79338146e+01,  -3.73074879e+01,\
   -6.91398584e+01,  -3.36441844e+01,  -5.00318069e+01,  -1.38785257e+01,\
   -2.57783826e+01,  -2.53094927e+01,  -3.10857250e+01,  -2.97158126e+01],\
 [  6.24697321e+00,  -1.67911892e+01,   5.05818596e+01,   2.69143104e+01,\
    3.72860575e+00,   5.58587483e+00,  -1.79989438e+01,  -3.71065495e+01,\
   -6.82866273e+01,  -3.35935956e+01,  -4.95249204e+01,  -1.43040726e+01,\
   -2.60772137e+01,  -2.53778038e+01,  -3.10810044e+01,  -2.97567337e+01],\
 [  5.65127839e+00,  -1.68396016e+01,   4.87342000e+01,   2.56792457e+01,\
    3.05390459e+00,   4.83855846e+00,  -1.80682275e+01,  -3.69026424e+01,\
   -6.74140564e+01,  -3.35478224e+01,  -4.90073366e+01,  -1.47453943e+01,\
   -2.63920488e+01,  -2.54505653e+01,  -3.10807655e+01,  -2.98024913e+01],\
 [  5.05301110e+00,  -1.68897700e+01,   4.68820992e+01,   2.44406551e+01,\
    2.37650642e+00,   4.08840294e+00,  -1.81401392e+01,  -3.67012133e+01,\
   -6.65447762e+01,  -3.35087127e+01,  -4.84925744e+01,  -1.51916715e+01,\
   -2.67157007e+01,  -2.55261822e+01,  -3.10855061e+01,  -2.98523058e+01],\
 [  4.46708661e+00,  -1.69404729e+01,   4.50717984e+01,   2.32294516e+01,\
    1.71330363e+00,   3.35411925e+00,  -1.82129421e+01,  -3.65073579e+01,\
   -6.57006127e+01,  -3.34774624e+01,  -4.79935957e+01,  -1.56318630e+01,\
   -2.70403261e+01,  -2.56028286e+01,  -3.10952542e+01,  -2.99050393e+01],\
 [  3.88400189e+00,  -1.69925414e+01,   4.32739567e+01,   2.20260011e+01,\
    1.05354132e+00,   2.62380201e+00,  -1.82878429e+01,  -3.63179623e+01,\
   -6.48679360e+01,  -3.34536989e+01,  -4.75023715e+01,  -1.60731537e+01,\
   -2.73713197e+01,  -2.56817737e+01,  -3.11102510e+01,  -2.99616012e+01],\
 [  3.30868770e+00,  -1.70455630e+01,   4.15038405e+01,   2.08405126e+01,\
    4.02802826e-01,   1.90363559e+00,  -1.83642642e+01,  -3.61346994e+01,\
   -6.40539358e+01,  -3.34378460e+01,  -4.70231836e+01,  -1.65118959e+01,\
   -2.77061068e+01,  -2.57624075e+01,  -3.11305162e+01,  -3.00216160e+01],\
 [  2.74578654e+00,  -1.70990960e+01,   3.97757022e+01,   1.96825363e+01,\
   -2.33662386e-01,   1.19942930e+00,  -1.84415828e+01,  -3.59590290e+01,\
   -6.32651388e+01,  -3.34300663e+01,  -4.65598905e+01,  -1.69445504e+01,\
   -2.80420057e+01,  -2.58440685e+01,  -3.11558974e+01,  -3.00845972e+01],\
 [  2.21124202e+00,  -1.71515163e+01,   3.81382427e+01,   1.85847444e+01,\
   -8.37842033e-01,   5.31102707e-01,  -1.85174568e+01,  -3.57956999e+01,\
   -6.25234087e+01,  -3.34301690e+01,  -4.61252883e+01,  -1.73586697e+01,\
   -2.83690413e+01,  -2.59242742e+01,  -3.11853674e+01,  -3.01485141e+01],\
 [  1.71713176e+00,  -1.72014038e+01,   3.66279201e+01,   1.75716655e+01,\
   -1.39611811e+00,  -8.63032378e-02,  -1.85898198e+01,  -3.56478917e+01,\
   -6.18444303e+01,  -3.34371175e+01,  -4.57284317e+01,  -1.77444353e+01,\
   -2.86786963e+01,  -2.60008264e+01,  -3.12175085e+01,  -3.02113380e+01],\
 [  1.24697491e+00,  -1.72502038e+01,   3.51938580e+01,   1.66092532e+01,\
   -1.92714233e+00,  -6.73437263e-01,  -1.86607572e+01,  -3.55101986e+01,\
   -6.12045572e+01,  -3.34501695e+01,  -4.53553689e+01,  -1.81142817e+01,\
   -2.89802386e+01,  -2.60759186e+01,  -3.12526865e+01,  -3.02746171e+01],\
 [  8.03323669e-01,  -1.72974905e+01,   3.38434723e+01,   1.57025435e+01,\
   -2.42805467e+00,  -1.22715315e+00,  -1.87296423e+01,  -3.53830178e+01,\
   -6.06065264e+01,  -3.34685378e+01,  -4.50075991e+01,  -1.84658829e+01,\
   -2.92712504e+01,  -2.61488778e+01,  -3.12901912e+01,  -3.03376048e+01],\
 [  3.88529666e-01,  -1.73428300e+01,   3.25835022e+01,   1.48561290e+01,\
   -2.89622567e+00,  -1.74456265e+00,  -1.87958315e+01,  -3.52666225e+01,\
   -6.00526657e+01,  -3.34912851e+01,  -4.46863581e+01,  -1.87970056e+01,\
   -2.95492793e+01,  -2.62190129e+01,  -3.13292184e+01,  -3.03995010e+01],\
 [  1.50659196e-02,  -1.73846184e+01,   3.14512853e+01,   1.40951804e+01,\
   -3.31761115e+00,  -2.21016883e+00,  -1.88569616e+01,  -3.51639817e+01,\
   -5.95585150e+01,  -3.35165836e+01,  -4.44004909e+01,  -1.90971957e+01,\
   -2.98047345e+01,  -2.62838104e+01,  -3.13677758e+01,  -3.04578181e+01],\
 [ -3.10029892e-01,  -1.74217648e+01,   3.04674636e+01,   1.34336831e+01,\
   -3.68431323e+00,  -2.61527566e+00,  -1.89114043e+01,  -3.50763558e+01,\
   -5.91319726e+01,  -3.35424764e+01,  -4.41543413e+01,  -1.93601587e+01,\
   -3.00312188e+01,  -2.63415345e+01,  -3.14040824e+01,  -3.05106550e+01],\
 [ -6.03864015e-01,  -1.74559741e+01,   2.95796984e+01,   1.28365379e+01,\
   -4.01566273e+00,  -2.98126375e+00,  -1.89616298e+01,  -3.49985781e+01,\
   -5.87494281e+01,  -3.35690947e+01,  -4.39340935e+01,  -1.95992027e+01,\
   -3.02393350e+01,  -2.63947978e+01,  -3.14391730e+01,  -3.05601278e+01],\
 [ -8.64596921e-01,  -1.74868479e+01,   2.87931263e+01,   1.23072680e+01,\
   -4.30961134e+00,  -3.30588930e+00,  -1.90070315e+01,  -3.49307244e+01,\
   -5.84124145e+01,  -3.35953560e+01,  -4.37404858e+01,  -1.98124391e+01,\
   -3.04268044e+01,  -2.64429521e+01,  -3.14721776e+01,  -3.06054338e+01],\
 [ -1.09057429e+00,  -1.75140097e+01,   2.81123253e+01,   1.18490207e+01,\
   -4.56431955e+00,  -3.58713898e+00,  -1.90470328e+01,  -3.48728217e+01,\
   -5.81222238e+01,  -3.36201857e+01,  -4.35741149e+01,  -1.99981278e+01,\
   -3.05914727e+01,  -2.64853830e+01,  -3.15022436e+01,  -3.06457997e+01],\
 [ -1.27289455e+00,  -1.75362027e+01,   2.75636865e+01,   1.14796289e+01,\
   -4.76978052e+00,  -3.81398197e+00,  -1.90797577e+01,  -3.48267319e+01,\
   -5.78894088e+01,  -3.36416547e+01,  -4.34408756e+01,  -2.01485502e+01,\
   -3.07258455e+01,  -2.65200975e+01,  -3.15275142e+01,  -3.06791286e+01],\
 [ -1.40933361e+00,  -1.75529764e+01,   2.71534915e+01,   1.12033884e+01,\
   -4.92351338e+00,  -3.98369716e+00,  -1.91045162e+01,  -3.47926132e+01,\
   -5.77159627e+01,  -3.36585777e+01,  -4.33417557e+01,  -2.02614804e+01,\
   -3.08273073e+01,  -2.65463623e+01,  -3.15470293e+01,  -3.07045236e+01],\
 [ -1.51556825e+00,  -1.75661361e+01,   2.68343311e+01,   1.09884173e+01,\
   -5.04319927e+00,  -4.11581594e+00,  -1.91239556e+01,  -3.47662715e+01,\
   -5.75813826e+01,  -3.36722704e+01,  -4.32649334e+01,  -2.03496282e+01,\
   -3.09068520e+01,  -2.65669846e+01,  -3.15625878e+01,  -3.07245696e+01],\
 [ -1.59327215e+00,  -1.75758174e+01,   2.66010126e+01,   1.08312445e+01,\
   -5.13073401e+00,  -4.21243823e+00,  -1.91382649e+01,  -3.47471296e+01,\
   -5.74832086e+01,  -3.36825755e+01,  -4.32089417e+01,  -2.04142248e+01,\
   -3.09653392e+01,  -2.65821649e+01,  -3.15741719e+01,  -3.07393850e+01],\
 [ -1.64547606e+00,  -1.75823481e+01,   2.64443227e+01,   1.07256820e+01,\
   -5.18953880e+00,  -4.27734531e+00,  -1.91479218e+01,  -3.47343294e+01,\
   -5.74173778e+01,  -3.36896375e+01,  -4.31714199e+01,  -2.04576812e+01,\
   -3.10047788e+01,  -2.65924096e+01,  -3.15820521e+01,  -3.07494116e+01],\
 [ -1.67464491e+00,  -1.75860065e+01,   2.63567938e+01,   1.06667101e+01,\
   -5.22239455e+00,  -4.31360965e+00,  -1.91533328e+01,  -3.47271985e+01,\
   -5.73806392e+01,  -3.36936322e+01,  -4.31504882e+01,  -2.04819830e+01,\
   -3.10268671e+01,  -2.65981500e+01,  -3.15864895e+01,  -3.07550397e+01],\
 [ -1.68882589e+00,  -1.75877875e+01,   2.63142456e+01,   1.06380426e+01,\
   -5.23836764e+00,  -4.33123960e+00,  -1.91559675e+01,  -3.47237371e+01,\
   -5.73627896e+01,  -3.36955870e+01,  -4.31403205e+01,  -2.04938031e+01,\
   -3.10376190e+01,  -2.66009451e+01,  -3.15886558e+01,  -3.07577826e+01],\
 [ -1.69664347e+00,  -1.75887700e+01,   2.62907914e+01,   1.06222398e+01,\
   -5.24717306e+00,  -4.34095833e+00,  -1.91574210e+01,  -3.47218305e+01,\
   -5.73529527e+01,  -3.36966681e+01,  -4.31347178e+01,  -2.05003207e+01,\
   -3.10435501e+01,  -2.66024871e+01,  -3.15898525e+01,  -3.07592966e+01],\
 [ -1.70031998e+00,  -1.75892323e+01,   2.62797616e+01,   1.06148082e+01,\
   -5.25131412e+00,  -4.34552888e+00,  -1.91581049e+01,  -3.47209342e+01,\
   -5.73483274e+01,  -3.36971775e+01,  -4.31320835e+01,  -2.05033862e+01,\
   -3.10463403e+01,  -2.66032125e+01,  -3.15904159e+01,  -3.07600090e+01],\
 [ -1.70206579e+00,  -1.75894518e+01,   2.62745242e+01,   1.06112792e+01,\
   -5.25328052e+00,  -4.34769923e+00,  -1.91584296e+01,  -3.47205087e+01,\
   -5.73461312e+01,  -3.36974195e+01,  -4.31308328e+01,  -2.05048420e+01,\
   -3.10476654e+01,  -2.66035571e+01,  -3.15906836e+01,  -3.07603475e+01],\
 [ -1.70326010e+00,  -1.75896020e+01,   2.62709412e+01,   1.06088651e+01,\
   -5.25462573e+00,  -4.34918396e+00,  -1.91586519e+01,  -3.47202176e+01,\
   -5.73446289e+01,  -3.36975852e+01,  -4.31299772e+01,  -2.05058379e+01,\
   -3.10485720e+01,  -2.66037928e+01,  -3.15908667e+01,  -3.07605790e+01],\
 [ -1.70445425e+00,  -1.75897522e+01,   2.62673588e+01,   1.06064513e+01,\
   -5.25597076e+00,  -4.35066849e+00,  -1.91588741e+01,  -3.47199266e+01,\
   -5.73431268e+01,  -3.36977509e+01,  -4.31291217e+01,  -2.05068337e+01,\
   -3.10494786e+01,  -2.66040286e+01,  -3.15910499e+01,  -3.07608106e+01],\
 [ -1.70583192e+00,  -1.75899255e+01,   2.62632259e+01,   1.06036666e+01,\
   -5.25752249e+00,  -4.35238115e+00,  -1.91591304e+01,  -3.47195909e+01,\
   -5.73413939e+01,  -3.36979422e+01,  -4.31281349e+01,  -2.05079826e+01,\
   -3.10505245e+01,  -2.66043006e+01,  -3.15912613e+01,  -3.07610777e+01],\
 [ -1.70748484e+00,  -1.75901334e+01,   2.62582673e+01,   1.06003256e+01,\
   -5.25938425e+00,  -4.35443600e+00,  -1.91594381e+01,  -3.47191882e+01,\
   -5.73393149e+01,  -3.36981717e+01,  -4.31269509e+01,  -2.05093610e+01,\
   -3.10517796e+01,  -2.66046269e+01,  -3.15915150e+01,  -3.07613983e+01],\
 [ -1.70950466e+00,  -1.75903875e+01,   2.62522080e+01,   1.05962429e+01,\
   -5.26165927e+00,  -4.35694696e+00,  -1.91598140e+01,  -3.47186962e+01,\
   -5.73367746e+01,  -3.36984524e+01,  -4.31255042e+01,  -2.05110455e+01,\
   -3.10533133e+01,  -2.66050257e+01,  -3.15918251e+01,  -3.07617902e+01],\
 [ -1.71189114e+00,  -1.75906878e+01,   2.62450489e+01,   1.05914192e+01,\
   -5.26434726e+00,  -4.35991372e+00,  -1.91602583e+01,  -3.47181149e+01,\
   -5.73337733e+01,  -3.36987842e+01,  -4.31237951e+01,  -2.05130359e+01,\
   -3.10551258e+01,  -2.66054971e+01,  -3.15921917e+01,  -3.07622533e+01],\
 [ -1.71482748e+00,  -1.75910573e+01,   2.62362405e+01,   1.05854841e+01,\
   -5.26765457e+00,  -4.36356402e+00,  -1.91608050e+01,  -3.47173998e+01,\
   -5.73300808e+01,  -3.36991928e+01,  -4.31216924e+01,  -2.05154850e+01,\
   -3.10573562e+01,  -2.66060771e+01,  -3.15926429e+01,  -3.07628233e+01],\
 [ -1.71831313e+00,  -1.75914960e+01,   2.62257844e+01,   1.05784388e+01,\
   -5.27158058e+00,  -4.36789717e+00,  -1.91614542e+01,  -3.47165512e+01,\
   -5.73256979e+01,  -3.36996783e+01,  -4.31191967e+01,  -2.05183925e+01,\
   -3.10600043e+01,  -2.66067658e+01,  -3.15931790e+01,  -3.07635001e+01],\
 [ -1.72253078e+00,  -1.75920270e+01,   2.62131328e+01,   1.05699142e+01,\
   -5.27633106e+00,  -4.37314029e+00,  -1.91622399e+01,  -3.47155247e+01,\
   -5.73203952e+01,  -3.37002665e+01,  -4.31161773e+01,  -2.05219108e+01,\
   -3.10632092e+01,  -2.66075993e+01,  -3.15938280e+01,  -3.07643195e+01],\
 [ -1.72747942e+00,  -1.75926502e+01,   2.61982890e+01,   1.05599123e+01,\
   -5.28190483e+00,  -4.37929205e+00,  -1.91631621e+01,  -3.47143206e+01,\
   -5.73141744e+01,  -3.37009575e+01,  -4.31126353e+01,  -2.05260393e+01,\
   -3.10669706e+01,  -2.66085776e+01,  -3.15945902e+01,  -3.07652813e+01],\
 [ -1.73352404e+00,  -1.75934117e+01,   2.61801582e+01,   1.05476955e+01,\
   -5.28871300e+00,  -4.38680620e+00,  -1.91642889e+01,  -3.47128505e+01,\
   -5.73065771e+01,  -3.37018030e+01,  -4.31083097e+01,  -2.05310828e+01,\
   -3.10715665e+01,  -2.66097730e+01,  -3.15955222e+01,  -3.07664570e+01],\
 [ -1.74075392e+00,  -1.75943229e+01,   2.61584731e+01,   1.05330837e+01,\
   -5.29685610e+00,  -4.39579365e+00,  -1.91656373e+01,  -3.47110929e+01,\
   -5.72974918e+01,  -3.37028162e+01,  -4.31031373e+01,  -2.05371160e+01,\
   -3.10770656e+01,  -2.66112035e+01,  -3.15966383e+01,  -3.07678641e+01],\
 [ -1.74944035e+00,  -1.75954182e+01,   2.61324205e+01,   1.05155288e+01,\
   -5.30663965e+00,  -4.40659161e+00,  -1.91672583e+01,  -3.47089826e+01,\
   -5.72865788e+01,  -3.37040365e+01,  -4.30969249e+01,  -2.05443659e+01,\
   -3.10836757e+01,  -2.66129232e+01,  -3.15979813e+01,  -3.07695563e+01]])

    T = np.array([  958.36,   958.36,   958.37,   958.38,   958.4 ,   958.41, \
         958.43,   958.46,   958.49,   958.52,   958.56,   958.61, \
         958.67,   958.75,   958.83,   958.94,   959.07,   959.22, \
         959.4 ,   959.63,   959.89,   960.22,   960.6 ,   961.07, \
         961.63,   962.31,   963.12,   964.09,   965.26,   966.66, \
         968.34,   970.35,   972.75,   975.61,   979.01,   983.06, \
         987.86,   993.52, 1000.17,  1007.96,  1017.06,  1027.65, \
        1039.86,  1053.75,  1069.59,  1087.54,  1107.77,  1130.21, \
        1154.58,  1181.29,  1210.33,  1241.7 ,  1274.79,  1308.85, \
        1344.49,  1381.49,  1419.59,  1457.64,  1494.55,  1531.33, \
        1567.62,  1603.03,  1636.21,  1666.15,  1694.1 ,  1719.64, \
        1742.36,  1761.1 ,  1775.37,  1786.63,  1794.95,  1800.58, \
        1803.74,  1805.28,  1806.13,  1806.53,  1806.72,  1806.85, \
        1806.98,  1807.13,  1807.31,  1807.53,  1807.79,  1808.11, \
        1808.49,  1808.95,  1809.49,  1810.15,  1810.94,  1811.89])

    for i in np.arange(len(T)):
        if float(temp) == T[i]:
            g_RT = g_RT_CEA[i]

    # Write header array to file
    write_header(desc, pressure, temp, stoich_arr, n_spec, g_RT)

