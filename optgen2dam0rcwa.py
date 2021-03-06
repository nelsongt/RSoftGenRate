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
projectname = '3JGeCell'
detectornumber = '1'
filename = '3JGeCell.plx'


## Get all files in the output directory
#filenames_to_glob = projectname + '_work/raw/' + projectname + '_*_m' + detectornumber + '_f*_absorption.dat' #FW output - FDTD
filenames_to_glob = projectname + '_work/raw/' + projectname + '_*_m' + detectornumber + '_absorption.dat' #FW output - RCWA
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
idata = pd.read_csv(sortedfiles[0],delim_whitespace=True,header=None,skiprows=4) # RSoft writes 4 header lines
simsize = idata.shape
xinfo = pd.read_csv(sortedfiles[0],delim_whitespace=True,header=None,skiprows=2,nrows=1)
xmin = xinfo.iloc[0,1]
xmax = xinfo.iloc[0,2]
xs = np.linspace(xmin,xmax,simsize[0])
yinfo = pd.read_csv(sortedfiles[0],delim_whitespace=True,header=None,skiprows=3,nrows=1)
ymin = yinfo.iloc[0,1]
ymax = yinfo.iloc[0,2]
ys = np.linspace(ymin,ymax,simsize[1])


xax = ys-ys[0]


#alldata = np.zeros((simsize[0], simsize[1], numfiles))
integrateddata = np.zeros((numfiles,simsize[1]))
wavelengths = np.zeros(numfiles)


## Collect all data
for i in xrange(numfiles):
    data = pd.read_csv(sortedfiles[i],delim_whitespace=True,header=None,skiprows=4) # RSoft writes 4 header lines
    #textdata = pd.read_csv(sortedfiles[i],delim_whitespace=True,header=None,skiprows=2,nrows=1)  #FDTD
    textdata = pd.read_csv(projectname + '_work/raw/' + projectname + '_' + str(i) + '.syms',delim_whitespace=True,header=None,nrows=1) #RCWA
    #wavelengthstr = textdata.loc[0,5].split('=')  # "wavelength = ..." #FDTD
    wavelengthstr = textdata.loc[0,0].split('=')  # "wavelength = ..." #RCWA
    wavelength = wavelengthstr[1]
    wavelengths[i] = wavelength

    #alldata[:,:,i] = data
    integrateddata[i,:] = np.sum(data,0)/simsize[0]   #Get averaged absorption across X direction
    print i+1

	

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
deltaWL = np.mean(np.diff(wavelengths))   # average wavelength step
gax = integrateddata*(np.tile((wavelengths*am0),(xax.size,1)).transpose())/(h*c)
yax = np.sum(gax,axis=0)*deltaWL/1000

rcwa=np.column_stack((xax,yax))


## Plot generation profile
plt.figure()
plt.semilogy(xax, yax)
plt.show()


## Write the result to Sentaurus PLX format
outdata = pd.DataFrame(data=rcwa) # can't include header as it contains the delimiter (space)
outheader = '\n'.join(
    [unicode(line, 'utf8') for line in 
        ['# from Sentaurus', 'Theta = 0 [deg] Intensity = 1.0 [W*cm^-2]\n']
    ]
)
with open(filename, 'w') as ict:
    for line in outheader:
        ict.write(line)
    outdata.to_csv(ict,sep=' ',float_format='%6.4e',index=False,header=None) # turn off data indexing and header, header written above directly
