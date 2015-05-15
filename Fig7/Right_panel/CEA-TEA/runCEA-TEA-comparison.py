

# ================================ Run CEA ===============================#


import os
import sys
import numpy as np

         
# Kevin's PT profile from WASP43b_PHASE7_NH3_H2S_TP.txt, cut to 1e-5 bars
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
         1.04712850e-05, 1e-5])

# Kevin's PT profile from WASP43b_PHASE7_NH3_H2S_TP.txt, cut to 1e-5 bars
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
         958.52146,   958.48696, 958.48])


# Input elements
in_elem  = 'H He C O N S'
# Output species
out_spec = 'H He C N O S H2 CO CO2 CH4 H2O HCN C2H4 N2 NH3 H2S'
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
                subprocess.call(['/home/jasmina/Work/TEA/TEA-APJSupp-paper/TEA-paper-runs/test-compendium/WASP-43b-comp/CEA/FCEA2-cmd', desc.encode()])

# Call the function
autoCEA(abun, pres, temp, in_elem, out_spec)


# ===================== Make results file from CEA out files ======================== #

import numpy as np
import glob
import string  

# Location of the CEA output files
out  = '/home/jasmina/Work/TEA/TEA-APJSupp-paper/TEA-paper-runs/test-compendium/WASP-43b-comp/CEA-outputs'
path = out + '/*.out'

# choose species of interest with the exact names as they appear in CEA out files
# NO COMMAS in between species
specs = '*H *CO *CO2 CH4 H2O HCN C2H4 *N2 NH3 H2S'      

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
species = 'H CO CO2 CH4 H2O HCN C2H4 N2 NH3 H2S' 


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
    data    = tuple(np.concatenate((['p'], species)))
    usecols = tuple(np.concatenate(([0], columns)))

    # Load all data for all interested species
    data = np.loadtxt(res_file, dtype=float, comments='#', delimiter=None,    \
                    converters=None, skiprows=1, usecols=usecols, unpack=True)
 
    # Sort pressure array and abundances array
    pres = sorted(data[0])
    abun = np.empty((nspec, len(data[0])), dtype=np.object)
    for i in np.arange(nspec):
        abun[i] = [x for (y,x) in sorted(zip(data[0],data[i+1]))]
    
    # Return name of plot created 
    return abun, pres

# Run plotCEA
abun, pres = plotCEA(res_file, species)

# Initiate figure
plt.ion()
plt.figure(1)

# Set different colours of lines (PLEASE DO NOT USE THE SAME COLORS FOR YOUR
#     PAPER!)
colors = 'b', '#FF1CAE','#FF0000' , '#FFAA00', '#00FFFF', '#00FF00','#ffc3a0', '#91219E', 'g', '#BCEE68', 'm' ,'c'
color_index = 0

# Plot all specs with different colours and labels
species = species.split()
for i in np.arange(len(species)):
    plt.loglog(abun[i], pres, '-', color=colors[color_index], \
                                    linewidth=2, label=str(species[i]))
    color_index += 1

# Label the plot
plt.xlabel('Mixing Fraction', fontsize=14)
plt.ylabel('Pressure [bar]' , fontsize=14)

# Temperature range (plt.xlim) and pressure range (plt.ylim)
plt.ylim(max(pres), min(pres))     
plt.xlim(1e-21, 2e-2)


# ================================ Plot TEA ============================== #

filename = '/home/jasmina/Work/TEA/TEA-APJSupp-paper/TEA-paper-runs/test-compendium/WASP-43b-comp/TEA/results/CEA-TEA-freeEnergiesTEA.tea'

species = 'H,CO,CO2,CH4,H2O,HCN,C2H4,N2,NH3,H2S'

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


# 'H,CO,CO2,CH4,H2O,HCN,C2H4,N2,NH3,H2S'
colors = 'k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k' 

color_index = 0

# Plot all specs with different colours and labels
for i in np.arange(nspec):
    if i==0 or i==3 or i==5 or i==6 or i==7 or i==9: 
        plt.loglog(data[i+1], data[0], ':', color=colors[color_index], \
                                        linewidth=2, label=str(spec[i]))
    else:
        plt.loglog(data[i+1], data[0], ':', color=colors[color_index], \
                                        linewidth=2, label=str(spec[i]))
    color_index += 1


# Label the plot
plt.xlabel('Mixing Fraction'    , fontsize=14)
plt.ylabel('Pressure [bar]'    , fontsize=14)

# WASP43b solar
plt.text(1.5e-8, 5e-5, 'H', color='b', fontsize=14)
plt.text(8e-4, 1, 'CO', color='#FF1CAE', fontsize=14)
plt.text(6e-7,5e-4, 'CO$_2$', color='#FF0000', fontsize=14)
plt.text(4e-12,4e-5, 'CH$_4$', color='#FFAA00', fontsize=14)
plt.text(8e-4,1e-3, 'H$_2$O', color='#00FFFF', fontsize=14)
plt.text(2.5e-13,8e-4, 'HCN', color='#00FF00', fontsize=14)
plt.text(7e-19,8e-4, 'C$_2$H$_4$', color='#ffc3a0', fontsize=14)
plt.text(7e-5,1e-2, 'N$_2$', color='#91219E', fontsize=14)
plt.text(1.5e-10,1e-3, 'NH$_3$', color='g', fontsize=14)
plt.text(1e-6,3e-5, 'H$_2$S', color='#BCEE68', fontsize=14)
plt.text(2e-15, 5e-5, 'CEA', color='#91219E',  fontsize=14)
plt.text(2e-15, 1.05e-4, 'TEA', color='#91219E', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

# Annotation for CEA and TEA
x1 = [5e-17, 1e-15] 
y1 = [5e-5, 5e-5]
y2 = [8e-5, 8e-5]

plt.plot(x1, y1,  '-', color = 'k', linewidth=2)
plt.plot(x1, y2,  ':', color = 'k', linewidth=2)


# ====================== Inset plot ================================

plt.xlim(1e-23, 2e-2)
b = plt.axes([.195, .20, .24, .30], axisbg='w')

plt.ylabel('Pressure [bar]', fontsize=10)
plt.xlabel('Mixing Fraction', fontsize=10)
plt.yticks(fontsize=9) 
plt.xticks(fontsize=9)

plt.ylim(max(pres), 1)
plt.xlim(1e-8, 1e-5) 

############## CEA #########################
species = 'H CO CO2 CH4 H2O HCN C2H4 N2 NH3 H2S'

# Set different colours of lines
colors = 'b', '#FF1CAE','#FF0000' , '#FFAA00', '#00FFFF', '#00FF00','#ffc3a0', '#91219E', 'g', '#BCEE68', 'm' ,'c'
color_index = 0

# Plot all specs with different colours and labels
species = species.split()
for i in np.arange(len(species)):
    plt.loglog(abun[i], pres, '-', color=colors[color_index], \
                                    linewidth=1, label=str(species[i]))
    color_index += 1

############## TEA #########################
species = 'H,CO,CO2,CH4,H2O,HCN,C2H4,N2,NH3,H2S'


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


# 'H,CO,CO2,CH4,H2O,HCN,C2H4,N2,NH3,H2S'
colors = 'k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k' 

color_index = 0

# Plot all specs with different colours and labels
for i in np.arange(nspec):
    if i==0 or i==3 or i==5 or i==6 or i==7 or i==9: 
        plt.loglog(data[i+1], data[0], ':', color=colors[color_index], \
                                        linewidth=1, label=str(spec[i]))
    else:
        plt.loglog(data[i+1], data[0], ':', color=colors[color_index], \
                                        linewidth=1, label=str(spec[i]))
    color_index += 1

plt.title('Detail', fontsize=10)
plt.text(2e-8,20, 'CO$_2$', color='#FF0000', fontsize=8)
plt.text(7e-7,2, 'CH$_4$', color='#FFAA00', fontsize=8)
plt.text(1.5e-7,10, 'HCN', color='#00FF00', fontsize=8)
plt.text(3e-7,4, 'NH$_3$', color='g', fontsize=8)


# ====================== Inset plot ================================

# free energies TEA
plt.savefig('CEA-TEA-WASP-43b.png', dpi=300)
plt.savefig('CEA-TEA-WASP-43b.ps', dpi=300)

