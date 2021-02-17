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
projectname = 'ProjectName-ShouldNotHaveUnderscores'
detectornumber = '2' # Component number for rsoft monitor
filename = 'ProjectGenerationRate.plx'


## Get all files in the work directory
filenames_to_glob = projectname + '_work/raw/' + projectname + '_*_m' + detectornumber + '_f*_absorption.vtk' #FW output
filelist = glob.glob(filenames_to_glob)
numfiles = len(filelist)
sortedfiles = [0]*numfiles


## Make a sorted list to work on, rsoft unfortunately uses 0, 1, 2,... 10,
## 11, .. which messes up the order. This fixes it.
for i in xrange(numfiles):
    nameparts = filelist[i].split('_')
    fileidx = nameparts[2]   # 1=projectname, 2='work/raw/projectname', 3=fileidx, 4=mdetectoridx, 5=freqidx, 6='absorption.dat'
    idxnum = int(fileidx)
    sortedfiles[idxnum] = filelist[i]
    

## Collect information on the simulation, like dimensions, mesh grid, etc.
meshpoints = pd.read_csv(sortedfiles[0],delim_whitespace=True,header=None,skiprows=5,nrows=1) #RSoft writes info in header
meshorigin = pd.read_csv(sortedfiles[0],delim_whitespace=True,header=None,skiprows=6,nrows=1) 
meshspacing = pd.read_csv(sortedfiles[0],delim_whitespace=True,header=None,skiprows=7,nrows=1) 
xmin = meshorigin.iloc[0,1]
xmax = xmin + (meshpoints.iloc[0,1]-1)*meshspacing.iloc[0,1] #don't include xmin in # of points
xs = np.linspace(xmin,xmax,meshpoints.iloc[0,1])
ymin = meshorigin.iloc[0,2]
ymax = ymin + (meshpoints.iloc[0,2]-1)*meshspacing.iloc[0,2]
ys = np.linspace(ymin,ymax,meshpoints.iloc[0,2])
zmin = meshorigin.iloc[0,3]
zmax = zmin + (meshpoints.iloc[0,3]-1)*meshspacing.iloc[0,3]
zs = np.linspace(zmin,zmax,meshpoints.iloc[0,3])


#alldata = np.zeros((simsize[0], simsize[1], numfiles))
integrateddata = np.zeros((numfiles,meshpoints.iloc[0,3]))
wavelengths = np.zeros(numfiles)
slice = meshpoints.iloc[0,1] * meshpoints.iloc[0,2] #number of points on an xy plane


## Collect all data
for i in xrange(numfiles):
    data = pd.read_csv(sortedfiles[i],delim_whitespace=True,header=None,skiprows=12) #RSoft writes 12 header lines
    wavelengthstr = pd.read_csv(projectname + '_work/raw/' + projectname + '_' + str(i) + '.syms',sep='=',header=None,nrows=1)
    wavelength = wavelengthstr.iloc[0,1]  # "excitlambda=..."
    wavelengths[i] = wavelength
    #wavelengths(i) = 0.3 + (i-1)*.85/99  # for RCWA output only
    
    #alldata[:,:,i] = data
    integrateddata[i,:] = np.mean(data.values.reshape(-1, slice), axis=1) #Get averaged absorption across xy plane   
    print i+1


xax = zs-zs[0]


## Setup the spectrum
#am15gspectrum = pd.read_csv('ASTMG173.csv',sep=',',header=None,skiprows=2)
#am15gfunc = interp1d(am0spectrum.loc[:,0]/1000, am0spectrum.loc[:,2])
#am15g = amfunc(wavelengths)

am0spectrum = pd.read_csv('ASTMG173.csv',sep=',',header=1)
am0func = interp1d(am0spectrum.iloc[:,0]/1000, am0spectrum.iloc[:,1])
am0 = am0func(wavelengths)


## Integrate data against spectrum
h=6.626e-34 # Js Planck's constant
c=2.998e8 #m/s speed of light
deltaWL = np.mean(np.diff(wavelengths))  #average wavelength step, must have at least 2 data files or line88 div0
gax = integrateddata*(np.tile((wavelengths*am0),(xax.size,1)).transpose())/(h*c)
yax = np.sum(gax,axis=0)*deltaWL/1000

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
with open(filename, 'w') as ict:
    for line in outheader:
        ict.write(line)
    outdata.to_csv(ict,sep=' ',float_format='%6.4e',index=False,header=None) # turn off data indexing and header, header written above directly
