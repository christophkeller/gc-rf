import argparse
import sys
import numpy as np
import random
import datetime
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from netCDF4 import Dataset
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, RandomForestRegressor
from sklearn.externals import joblib
from sklearn.tree import _tree

# time step in seconds
dt_chem = 15.0*60.0

def read_files(filenames,ispec='O3',fraction=1.0,nbins=10,outtype='tend'):
   '''
   Prepares the data arrays for the random forest (input and output features).
   
   Parameters
   ----------
   filenames : str
       list of pickle file names to be read
   ispec : str 
       Species name of target variable
   tend : bool
       Target variable is tendency or absolute concentration
   fraction : float
       Fraction of data to be used
   nbins : int
       Number of bins to evenly draw samples from 
   outtype : str
       Output type
   '''
   first = True
   nfiles = len(filenames)
   for filename in filenames:
       print 'Reading '+filename
       [fvars, pvars, farr, parr] = joblib.load(filename)
       # reduce output to species of interest
       if ispec=='NOx' or ispec=='NOratio':
           # get NO and NO2 in kg N
           pidx1 = pvars.index('KPP_AFTER_CHEM_NO')
           pidx2 = pvars.index('KPP_AFTER_CHEM_NO2')
           ptmp1 = np.copy(parr[:,pidx1]) * (14./30.)
           ptmp2 = np.copy(parr[:,pidx2]) * (14./46.)
           # get NOx
           if ispec=='NOx':
               ptmp  = ptmp1 + ptmp2
           if ispec=='NOratio':
               ptmp  = ptmp1 / ( ptmp1 + ptmp2 )
       else:
           pidx = pvars.index('KPP_AFTER_CHEM_'+ispec)
           ptmp = np.copy(parr[:,pidx])
       if outtype=='tend':
           if ispec=='NOx':
               fidx1 = fvars.index('KPP_BEFORE_CHEM_NO')
               fidx2 = fvars.index('KPP_BEFORE_CHEM_NO2')
               ftmp1 = np.copy(farr[:,fidx1]) * (14./30.)
               ftmp2 = np.copy(farr[:,fidx2]) * (14./46.)
               ftmp  = ftmp1 + ftmp2
           else:
               fidx = fvars.index('KPP_BEFORE_CHEM_'+ispec)
               ftmp = farr[:,fidx]
           ptmp = ( ptmp - ftmp ) / dt_chem
       if first:
           lens = farr.shape[0]
           nrows_per_file = int(lens*fraction)
           nrows_total    = nrows_per_file*nfiles
           nrows_per_bin  = nrows_per_file/nbins
           megaf = np.zeros((nrows_total,farr.shape[1]),dtype='float32')
           megap = np.zeros((nrows_total,1),dtype='float32')
           first=False
           i1 = 0
       # get percentiles
       ranges = [np.percentile(ptmp,i) for i in range(0,101,100/nbins)]
       for ibin in range(0,nbins):
           # get all indeces that are within this percentile
           idxs = np.where( (ptmp >= ranges[ibin]) & (ptmp <= ranges[ibin+1]) )[0]
           # randomly pick values
           idx = np.random.choice(idxs,nrows_per_bin,replace=False)
           # pass to master array
           i0 = i1
           i1 = i0 + nrows_per_bin
           megaf[i0:i1,:] = farr[idx,:]
           megap[i0:i1,0] = ptmp[idx]
   # remove entries where NUMDEN is zero
   idx = fvars.index('KPP_AIRDEN')
   msk = np.where(megaf[:,idx]>0.0)[0]
   megaf = megaf[msk,:]
   megap = megap[msk,:]
   # scale inputs & output. This improves the quality of the fit
   idx = fvars.index('KPP_AIRDEN')
   if ispec=='NOratio':
       megap[:,0] = megap[:,0] * 1.0e3
   else:
       megap[:,0] = megap[:,0]/megaf[:,idx]*1.0e15
   idx2 = [i for i in fvars if 'KPP_BEFORE_CHEM' in i]
   for ivar in idx2:
       ii = fvars.index(ivar)
       megaf[:,ii] = megaf[:,ii]/megaf[:,idx]*1.0e15
   if outtype == 'logconc': 
       msk = np.where(megap[:,0]<=0.0)[0]
       if len(msk) > 0:
           megap[msk,0] = 1.0e-20*megaf[msk,idx]
       megap[:,0] = np.log10(megap[:,0])
   if outtype=='ratio':
       idx = fvars.index('KPP_BEFORE_CHEM_'+ispec)
       megap[:,0] = ( megap[:,0] / megaf[:,idx] ) * 1000.0
   # eventually drop target variable from input features
   if outtype == 'drop': 
       idx = fvars.index('KPP_BEFORE_CHEM_'+ispec)
       megaf[:,idx] = 0.0
   return megaf,megap,fvars,pvars

def GCspecies(template):
    '''
    Returns the list of all species in the template file, plus NOx and NOratio.
    '''
    fid = Dataset(template)
    specs = [i.split('_')[3] for i in fid.variables if 'KPP_BEFORE_CHEM' in i]
    specs.append('NOx')
    specs.append('NOratio')
    return(specs)

def getMeta(args):
    '''
    Wrapper routine to get meta data based of the species number (species name and file prefix).
    '''
    specs = GCspecies(args.template)
    ispec = specs[args.specnumber-1]
    pfix  = 'spc'+'{:03d}'.format(args.specnumber)+'_'+ispec+'_'+'{:02d}'.format(args.treenumber)
    return ispec,pfix
