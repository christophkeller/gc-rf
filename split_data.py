import numpy as np
import glob
import sys
import os
import shutil
import argparse
import ntpath

def main(args):
    '''
    Splits files into training and validation data.
    
    Parameters
    ----------
    source : str
        file names (including path)
    training : str
        target location for training data
    validation : str
        target location for validation data
    fraction : float
        fraction of files used for validation
    '''
    if not os.path.exists(args.training):
        os.mkdir(args.training)
    if not os.path.exists(args.validation):
        os.mkdir(args.validation)
    filenames = glob.glob(args.source)
    train_files = np.random.choice(filenames,int(float(len(filenames))*args.fraction),replace=False)
    for ii in range(0,len(train_files)):
        shutil.move(train_files[ii],args.validation+'/'+ntpath.basename(train_files[ii]))
    print('Moved '+str(len(train_files))+ ' files to '+args.validation)
    filenames = glob.glob(args.source)
    for ii in range(0,len(filenames)):
        shutil.move(filenames[ii],args.training+'/'+ntpath.basename(filenames[ii]))
    print('Moved '+str(len(filenames))+ ' files to '+args.training)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--source',type=str,default="data/pkl/*.pkl")
    parser.add_argument('-t','--training',type=str,help='location for training data',default='data/pkl/training')
    parser.add_argument('-v','--validation', type=str,help='location for validation data',default='data/pkl/validation')
    parser.add_argument('-f','--fraction',type=float,help='fraction of validation data',default=0.25)
    return parser.parse_args() 

if __name__ == '__main__':
    main(parse_args())

# eof
