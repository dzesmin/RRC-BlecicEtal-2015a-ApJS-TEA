

# ================================ Run CEA ===============================#

import os
import sys
import numpy as np

# Number of layers   
steps = 30
pres = np.linspace( 1,    1, steps)
temp = np.linspace( 100, 3000, steps)

# Input elements
in_elem  = 'H C N O'
# Output species
out_spec = 'H C N O H2 CO CH4 H2O N2 NH3'

# The abundance profile file
abun     = 'abundances.txt'

def autoCEA(abun, pres, temp, in_elem, out_spec):
    ''' 
    Function that runs CEA with multiple (T-P) points.
    It creates the multiple .inp files for a CEA run and calls CEA
    and produces output files.
        
    Modification history: 06-01-2012  Oliver Bowman    intial version
                          11-01-2014  Jasmina Blecic   
                                      made the function work with multiple
                                      T-P points, introduce new arguments
                                      rewrote data handling
    '''
                                                                                                           
    if abun == 'abundances.txt':
        file = str("TEA")

        for data in np.arange(len(temp)): 
                file_string = (np.str(temp[data]) + "K_" + np.str(pres[data]) + \
                              "bar_" + file)
                desc = str(file_string)

                # Import abundance profile
                f = open(abun, 'r')
                abundata = []
                for line in f.readlines():
                    l = [value for value in line.split()]
                    abundata.append(l)               
                abundata = np.asarray(abundata)
                f.close()
            
                # Set up abundance conversions
                n_ele = abundata.shape[0]
                sum = 0
                for i in np.arange(n_ele):
                    sum += 10**np.float(abundata[i][2])

                # Create class for data storage
                class Element():
                    def __init__(self):
                        self.foo = 'DUMMY STRING'
                        self.id  = 'Element Symbol'
                        self.name = 'Element Name'
                        self.atnum  = []
                        self.weight = []
                        self.lgdex  = []
                        self.abun   = []
                        self.relw   = []
                    def __call__(self):
                        print('Element Name:        ' +      self.id     )
                        print('Element Symbol:      ' +      self.name   )
                        print('Atomic Number:       ' + str( self.atnum ))
                        print('Atomic Weight:       ' + str( self.weight))
                        print('Log abundance (dex): ' + str( self.lgdex ))
                        print('Solar Abundance:     ' + str( self.abun  ))
                        print('Relative Weight:     ' + str( self.relw  ))
                    
                # Fill classes for each element
                for i in np.arange(n_ele):
                    atnum  = np.float(abundata[i][0])
                    s      =          abundata[i][1]
                    lgdex  = np.float(abundata[i][2])
                    name   =          abundata[i][3]
                    weight = np.float(abundata[i][4])
                    globals()[s] = Element()
                    globals()[s].id     = s
                    globals()[s].name   = name
                    globals()[s].atnum  = atnum
                    globals()[s].weight = weight
                    globals()[s].lgdex  = lgdex
                    globals()[s].abun   = 10**lgdex / 10**12
                    globals()[s].relw   = 10**lgdex / sum * weight
    
                # == Begin writing .inp file == #
            
                # Create desired   ------ pressure range string
                pres_range = np.array(pres[data])
                pres_list = pres_range.tolist() 
                pres_list = '%.8f' % pres_list
                pres_str = np.str(pres_list) + ','
                #                  ------ temperature range string
                temp_range = np.array(temp[data])
                temp_list = temp_range.tolist()       
                if temp[data] != temp_list:
                    for j in np.arange(np.size(temp_list)):
                        temp_list[j] = '%.8f' % temp_list[j]
                        temp_str = ','.join(temp_list) + ','
                else:
                    temp_str = np.str(temp[data]) + ','           
                #                  ------ abundances
                if in_elem != '':
                    speclist = in_elem.split(' ')
                    n_spec   = np.size(speclist)
                    abunlist = np.zeros(n_spec)
                    for i in np.arange(n_spec):
                        abunlist[i] = globals()[speclist[i]].abun

                # Write to .inp file
                f = open(desc + '.inp', 'wb')
                f.write('problem    case=' + desc + '\n')
                f.write('     tp    t,k=' + temp_str + '  p,bar=' + pres_str + '\n') 
                f.write('react\n')
                for i in np.arange(n_spec):
                    f.write('  name=' + speclist[i] + ' moles=' + \
                            '%.10f'%abunlist[i] + '\n')
                f.write('only\n')
                f.write('  ' + out_spec + '\n')
                f.write('output    trace= 1.0E-100\n')
                f.write('end')
                f.close()
            
                # Run CEA
                import subprocess
                print("  PYTHON:  Wrapper is passing \"" + desc + ".inp\" to Fortran...")
                subprocess.call(['/home/jasmina/Work/TEA/TEA-APJSupp-paper/TEA-paper-runs/test-compendium/CEA-TEA-TEBS/CEA/FCEA2-cmd', desc.encode()])

# Call the function
autoCEA(abun, pres, temp, in_elem, out_spec)


# ===================== Make results file from CEA out files ======================== #

import numpy as np
import glob
import string  

# Location of the CEA output files
out  = '/home/jasmina/Work/TEA/TEA-APJSupp-paper/TEA-paper-runs/test-compendium/CEA-TEA/CEA-outputs'
path = out + '/*.out'

# choose species of interest with the exact names as they appear in CEA out files
# NO COMMAS in between species
specs = '*C CH4 *CO *H *H2 H2O *N NH3 *N2 *O'      

# Desired name of the results file
res_file = 'CEA-TEA.dat'

def CEAresults(path, specs, res_file):
    '''
    This routine makes a result file from CEA output files. 
    It takes the directory name, desired species, and the desired name of the
    result file, opens the results file, writes header and species abundances
    as they appaear in the CEA output files. For species not listed in the CEA
    output files, it sets zeroes
    '''

    # Splits the specs string to separate species strings
    specs = specs.split()

    # Open file to write CEA results
    results = open(res_file, 'w')

    # Make header
    header_file = '#Source directory:   ' + str(path.split("*")[0])
    results.write(header_file)
    results.write('\n')
    header = ['#Pressure'.ljust(12)] + [' Temp'.ljust(12)]
    for i in np.arange(len(specs)):
        header = header + [specs[i].ljust(14)]
    
    # Write header into the results file
    for i in np.arange(len(header)):
        results.write(header[i])
    results.write('\n')
    
    # Take all files from the CEA out directory and sort them 
    files=glob.glob(path)
    files.sort()
 
    for file in files:     
        f=open(file, 'r')  
    
        # Initiate abundances array with zeros fos species not occuring in CEA out
        abundances = ['0.0000e+0'  for i in np.arange(len(specs))]

        # allocates arrays of species and abundances
        species = []
        abun    = []
    
        # Take temp and pressure from CEA out file
        for line in f:
            if string.find(line, '      tp    t,k=') >-1:
                info      = line.split()
                pres_info = info[2]
                pres      = (pres_info.partition('=')[-1]).partition(',')[0]
                temp_info = info[1]
                temp      = (temp_info.partition('=')[-1]).partition(',')[0]
                temp      = '%.2f'%float(temp)

            # Read data in between two markers
            if string.find(line, ' MOLE FRACTIONS') > -1:
                # reads the atmospheric file from the row below the header
                for line in f:
                    if string.find(line, "  * THERMODYNAMIC PROPERTIES FITTED") > -1:
                        for line in f:
                            pass
                    else:
                        data = line.split()
                        if data != []:
                            # Replace - with e-
                            data[1] = data[1].replace('-', 'e-')
    
                            # append the data to species and abundances
                            species = (np.append(species, data[0])).tolist()
                            abun    = (np.append(abun, data[1])).tolist()

        # Write pressure, temperature,
        results.write(pres.ljust(12) + ' ')
        results.write(temp.ljust(10) + ' ')
    
        # Check whether the desired species appear in CEA out and write in the results file
        for i in np.arange(len(specs)):
            for j in np.arange(len(species)):
                if specs[i] == species[j]:
                    abundances[i] = abun[j]
            results.write(abundances[i].ljust(13) + ' ')
    
        results.write('\n')
         
    # Close results file
    results.close()
    print 'The result file is created.'
    
# Call the function
CEAresults(path, specs, res_file)


# ================================ Plot CEA ============================== #

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from PIL import Image

plt.ion()

# Name of the results file
res_file = 'CEA-TEA.dat'

# Species names. MUST BE IN THE CORRECT ORDER
species = 'CH4 CO H2O NH3 N2'  

# Open a figure
plt.figure(1)
plt.clf()


def plotCEA(res_file, species):
    '''
    This function reads the CEA results file and plots the results.
    '''


    # Open the atmospheric file and read
    f = open(res_file, 'r')
    lines = np.asarray(f.readlines())
    f.close()

    # Get molecules names without *
    header = lines[1].split()
    molecules = header[2:]
    nmol = len(molecules)
    for m in np.arange(nmol):
        if not molecules[m].find('*'):
            molecules[m] = molecules[m].replace('*', '')

    # Take user input for species and split species strings into separate strings 
    #      convert the list to tuple
    species = tuple(species.split())
    nspec = len(species)
    molecules = tuple(molecules)

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
    data    = tuple(np.concatenate((['T'], species)))
    usecols = tuple(np.concatenate(([1], columns)))

    # Load all data for all interested species
    data = np.loadtxt(res_file, dtype=float, comments='#', delimiter=None,    \
                    converters=None, skiprows=2, usecols=usecols, unpack=True)
 
    # Set different colours of lines
    colors = 'krcbg'
    color_index = 0

    # Sort pressure array and abundances array
    temp = sorted(data[0])
    abun = np.empty((nspec, len(data[0])), dtype=np.object)
    for i in np.arange(nspec):
        abun[i] = [x for (y,x) in sorted(zip(data[0],data[i+1]))]

    abun = abun.tolist()

    # Plot all specs with different colours and labels
    for i in np.arange(nspec):
        if i==0:
            plt.semilogy(temp, abun[i], ':', color=colors[color_index], \
                                                       linewidth = 2, label='CEA')
        else:
            plt.semilogy(temp, abun[i], ':', color=colors[color_index], \
                                        linewidth=2)
        color_index += 1

    # Label the plot
    plt.xlabel('Temperature [K]', fontsize=14)
    plt.ylabel('Log10 Mixing Fraction' , fontsize=14)
    plt.legend(loc='best', prop={'size':10})
 
    return

# Call the function
plotCEA(res_file, species)



# ================================ Plot TEA ============================== #

filename = '/home/jasmina/Work/TEA/TEA-APJSupp-paper/TEA-paper-runs/test-compendium/CEA-TEA/TEArun-freeEnergiesCEA/results/CEA-TEA-comp.tea'
species = 'CH4,CO,H2O,NH3,N2'

def plotTEA(filename, species):
    '''
    This function plots TEA results.
    '''
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
    nspec = len(species)

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

    # Concatenate spec with temperature for data and columns
    data    = tuple(np.concatenate((['T'], spec)))
    usecols = tuple(np.concatenate(([1], columns)))

    # Load all data for all interested specs
    data = np.loadtxt(filename, dtype=float, comments='#', delimiter=None,    \
                    converters=None, skiprows=8, usecols=usecols, unpack=True)

    # Set different colours of lines
    colors = 'rkcgb'
    color_index = 0

    # Plot all specs with different colours and labels
    for i in np.arange(nspec):
        if i==1:
            plt.semilogy(data[0], data[i+1], '-', color=colors[color_index], \
                                                       linewidth=1, label='TEA')

        plt.semilogy(data[0], data[i+1], '-', color=colors[color_index], linewidth=1)
        color_index += 1

# Call the function
plotTEA(filename, species)



# ================================ Annotation============================== #

# Attotation
plt.annotate('H$_{2}$O' , (850, 1.3e-3), fontsize=14, color='cyan')
plt.annotate('CH$_{4}$', (1450, 1e-5), fontsize=14, color='k')
plt.annotate('CO', (2500, 8e-4), fontsize=14, color='r')
plt.annotate('N$_{2}$', (2100, 1e-4), fontsize=14, color='g')
plt.annotate('NH$_{3}$', (1140, 2.5e-6), fontsize=14, color='b')

# Label
plt.legend(loc='lower center', prop={'size':12})

plt.xlabel('Temperature [K]' , fontsize=14)
plt.ylabel('Mixing Fraction', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlim(600, 3000)
plt.ylim(3e-10, 3e-3)

# Save figures
plt.savefig("CEA-TEA-freeEnergiesCEA.png", dpi=300)
plt.savefig("CEA-TEA-freeEnergiesCEA.ps", dpi=300)



