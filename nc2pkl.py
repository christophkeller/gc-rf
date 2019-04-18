import numpy as np
import glob
import sys
import os
import myrf
import argparse
from netCDF4 import Dataset
from sklearn.externals import joblib

def nc2pkl(args): 
    '''
    Extracts all features from GEOS-Chem netCDF files and packs them into pickled files.

    Parameters
    ----------
    Parser object with the following arguments:
    template : str
        filename of template netCDF file (needed for variable selection)
    files : str
        path with filename(s) to GEOS-Chem netCDF files
    outdir : str
        location for writing the pickled files into

    Returns
    -------
    None
    '''
    gc_specs = myrf.GCspecies(args.template)
    files = glob.glob(args.files)
    for file in files:
        filename = file.split('/')[-1]
        print 'Reading ',file
        sys.stdout.flush()
        fid = Dataset(file,'r')
        # Select all features: J-values, species concentrations before
        # chemistry, and some meteorological variables 
        vars1 = [var for var in fid.variables if 'KPP_JVAL' in var]
        jvals = []
        for ivar in vars1:
            if np.sum(fid.variables[ivar]) > 0.0:
                jvals.append(ivar)
        vars2 = []
        for spec in gc_specs:
            if 'KPP_BEFORE_CHEM_'+spec in fid.variables:
                vars2.append('KPP_BEFORE_CHEM_'+spec)
        vars = jvals + vars2
        vars.append('KPP_SUNCOS')
        vars.append('KPP_AIRDEN')
        vars.append('KPP_TEMP')
        vars.append('KPP_PRESS')
        vars.append('KPP_RH')
        vars.append('KPP_QLIQ')
        vars.append('KPP_QICE')
        nrow = fid.variables['KPP_TEMP'].shape[1] * fid.variables['KPP_TEMP'].shape[2] * fid.variables['KPP_TEMP'].shape[3]
        in_arr = np.zeros((nrow,len(vars)),dtype='float32')
        ivar = 0
        for var in vars:
            iarr = fid.variables[var][0,:,:,:]
            in_arr[:,ivar] = np.reshape(iarr,(nrow))
            ivar = ivar+1
        # Select all target variables: these are the species concentrations after chemistry 
        pvars = []
        for spec in gc_specs:
            if 'KPP_AFTER_CHEM_'+spec in fid.variables:
                pvars.append('KPP_AFTER_CHEM_'+spec)
        pred_arr = np.zeros((nrow,len(pvars)),dtype='float32')
        ivar = 0
        for var in pvars:
            iarr = fid.variables[var][0,:,:,:]
            pred_arr[:,ivar] = np.reshape(iarr,(nrow))
            ivar = ivar+1
        # Write out everything 
        if not os.path.exists(args.outdir):
            os.mkdir(args.outdir)
        ofile = args.outdir+'/'+filename+'.pkl'
        # print statements for debugging 
        #print 'Shape of in_arr: ',in_arr.shape[0],in_arr.shape[1],in_arr.dtype
        #print 'Shape of pred_arr: ',pred_arr.shape[0],pred_arr.shape[1],pred_arr.dtype
        sys.stdout.flush()
        joblib.dump([vars,pvars,in_arr,pred_arr], ofile)
        print 'Written to '+ofile
        # all done
        fid.close()

def parse_args():
    '''
    Parser for nc2pkl
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-t','--template',type=str,default="data/nc4/c48-ml.tavg1_3d_kpp_Nv.20130702_0000z.nc4")
    parser.add_argument('-f','--files',type=str,default="data/nc4/c48-ml.tavg1_3d_kpp_Nv.*.nc4")
    parser.add_argument('-o','--outdir',type=str,default="data/pkl")
    return parser.parse_args()

if __name__ == '__main__':
    nc2pkl(parse_args())
