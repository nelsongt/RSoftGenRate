from __future__ import division
import os
import csv
import glob
import math
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# ---- PARAMETERS TO ADJUST ----
prefixname = 'Pyra2NoOx'
detectornumber = '2' # Component number for rsoft monitor
outputfile = '3x3.plx'


## Get all files in the output directory
filename_to_glob = '3x3Rand_m2_f1_absorption.vtk' #FW output
filelist = glob.glob(filename_to_glob)



## Collect information on the simulation, like dimensions, mesh grid, etc.
meshpoints = pd.read_csv(filelist[0],delim_whitespace=True,header=None,skiprows=5,nrows=1) #RSoft writes info in header
meshorigin = pd.read_csv(filelist[0],delim_whitespace=True,header=None,skiprows=6,nrows=1) 
meshspacing = pd.read_csv(filelist[0],delim_whitespace=True,header=None,skiprows=7,nrows=1) 
xmin = meshorigin.iloc[0,1]
xmax = xmin + (meshpoints.iloc[0,1]-1)*meshspacing.iloc[0,1] #don't include xmin in # of points
xs = np.linspace(xmin,xmax,meshpoints.iloc[0,1])
ymin = meshorigin.iloc[0,2]
ymax = ymin + (meshpoints.iloc[0,2]-1)*meshspacing.iloc[0,2]
ys = np.linspace(ymin,ymax,meshpoints.iloc[0,2])
zmin = meshorigin.iloc[0,3]
zmax = zmin + (meshpoints.iloc[0,3]-1)*meshspacing.iloc[0,3]
zs = np.linspace(zmin,zmax,meshpoints.iloc[0,3])

print xmin
print xmax
print ymin
print ymax
print zmin
print zmax

integrateddata = np.zeros(meshpoints.iloc[0,3])
slice = meshpoints.iloc[0,1] * meshpoints.iloc[0,2] #number of points on an xy plane


## Collect all data
data = pd.read_csv(filelist[0],delim_whitespace=True,header=None,skiprows=12) #RSoft writes 12 header lines
wavelength = 900

integrateddata[:] = np.mean(data.values.reshape(-1, slice), axis=1) #Get averaged absorption across xy plane   

print integrateddata

xax = zs-zs[0]



yax = integrateddata

fdtd=np.column_stack((xax,yax))


## Plot generation profile
plt.figure()
plt.semilogy(xax, yax)
plt.show()


## Write the result to Sentaurus PLX format
outdata = pd.DataFrame(data=fdtd) # can't include header as it contains the delimiter (space)
outheader = '\n'.join(
    [unicode(line, 'utf8') for line in 
        ['# from Sentaurus', 'Theta = 0 [deg] Intensity = 1.0 [W*cm^-2]\n']
    ]
)
with open(outputfile, 'w') as ict:
    for line in outheader:
        ict.write(line)
    outdata.to_csv(ict,sep=' ',float_format='%6.4e',index=False,header=None) # turn off data indexing and header, header written above directly
