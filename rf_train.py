import argparse
import numpy as np
import glob
import sys
import random
import glob
import myrf
import os
from scipy import stats 
import sklearn
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, RandomForestRegressor
from sklearn.ensemble import ExtraTreesRegressor 
from sklearn.externals import joblib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from sklearn.metrics import mean_squared_error, r2_score
from math import sqrt
from scipy import stats

def rf_train(args):
    '''
    Train random forest for GEOS-Chem chemistry mechanism

    This routine creates a random forest for a given GEOS-Chem species, 
    using as inputs and outputs the data from the pickled GEOS-Chem data
    files. The various parameters for training (tendencies vs. absolute
    concentrations, data thinning, etc.) are passed via the argument object.

    Parameters
    ----------
    template : str
        filename of template netCDF file (needed for variable selection)
    specnumber : int
        species number of target species.
    treenumber : int
        tree number to be given to forest. Ignore if not training trees separately
    nbins : int
        Number of bins to evenly draw samples from 
    fraction : float
        Fraction of data to be used
    source-training : str
        File names (incl. path) used for training
    source-validation: str
        File names (incl. path) used for validation 
    outtype : str 
        Target variable type: 'tend', 'conc', 'logconc', 'ratio', or 'drop'
    nleaves : int
        number of leaves per tree
    ntrees : int
        number of trees per forest
    treepath : str
        path to save tree file

    Returns
    -------
    None
    '''
    # get metadata
    ispec,pfix = myrf.getMeta(args)
    # check if file already exists - skip if so
    ofile = args.treepath+'/tree_'+pfix+'_'+args.outtype+'.sav'
    if os.path.isfile(ofile):
        print('Tree already exists - leaving!')
        sys.stdout.flush()
        sys.exit()     
    # set flags based on outtype
    if args.outtype=="tend":
        action  = 'T'
    elif args.outtype=="logconc":
        action  = 'L'
    elif args.outtype=="conc":
        action  = 'C'
    elif args.outtype=="ratio":
        action  = 'R'
    elif args.outtype=="drop":
        action  = 'D'
    # now read full array of files
    filenames = glob.glob(args.source_training)
    farr,parr,fvars,pvars = myrf.read_files(filenames,ispec,args.fraction,args.nbins,args.outtype)
    print 'Shape of farr: ',farr.shape
    print 'Shape of parr: ',parr.shape
    # train forest
    print('training tree...')
    rf = RandomForestRegressor(n_estimators=args.ntrees,n_jobs=1, max_leaf_nodes=args.nleaves, criterion='mse', bootstrap=True, random_state=2, verbose=5)
    rf.fit(farr,parr[:,0])
    ofile = 'figures/'+pfix+'_scatter_training.png'
    plot_true_vs_pred(rf,farr,parr,fvars,ispec,args.outtype,ofile,False)
    if args.outtype == 'tend':
        ofile = 'figures/'+pfix+'_scatter_training_tend+conc.png'
        plot_true_vs_pred(rf,farr,parr,fvars,ispec,args.outtype,ofile,True)
    # validation
    filenames = glob.glob(args.source_validation)
    farr,parr,fvars,pvars = myrf.read_files(filenames,ispec,args.fraction,args.nbins,args.outtype)
    ofile = 'figures/'+pfix+'_scatter_validation.png'
    plot_true_vs_pred(rf,farr,parr,fvars,ispec,args.outtype,ofile,False)
    if args.outtype == 'tend':
        ofile = 'figures/'+pfix+'_scatter_validation_tend+conc.png'
        plot_true_vs_pred(rf,farr,parr,fvars,ispec,args.outtype,ofile,True)
    # plot features
    ofile = 'figures/'+pfix+'_feature_importances.png'
    plot_features(rf,ispec,fvars,ofile)
    # write forest to ascii file 
    opath = args.treepath+'/'+args.outtype+'/'+ispec
    if not os.path.exists(opath):
        os.mkdir(opath)
    # correction factor
    corr = 1.0
    ofile = opath + '/'+ispec
    forest_write(rf,fvars,ispec,action,ofile,corr)
    # to save python file
    #joblib.dump(rf,ofile)
    #print('Written to '+ofile)
    return

def forest_write(forest,vars,spec,action,ofile,corr=1.0):
    '''
    Helper routine to write forest to ascii file. This ascii
    file can then be read by the Fortran routine.
    '''
    ntrees    = len(forest.estimators_)
    maxnodes  = 0
    for tree in forest.estimators_:
        if tree.tree_.node_count > maxnodes:
            maxnodes = tree.tree_.node_count
    hdr       = 'Trees: '+str(ntrees)+'\nNodes: '+str(maxnodes)
    left      = np.zeros((maxnodes,ntrees))
    right     = np.zeros((maxnodes,ntrees))
    feature   = np.zeros((maxnodes,ntrees))
    threshold = np.zeros((maxnodes,ntrees))
    value     = np.zeros((maxnodes,ntrees))
    itree     = 0
    # testing only
    for tree in forest.estimators_:
        nn                    = tree.tree_.node_count 
        left[0:nn,itree]      = tree.tree_.children_left[0:nn]
        right[0:nn,itree]     = tree.tree_.children_right[0:nn]
        feature[0:nn,itree]   = tree.tree_.feature[0:nn]
        threshold[0:nn,itree] = tree.tree_.threshold[0:nn]
        value[0:nn,itree]     = tree.tree_.value[0:nn,0,0]
        itree                 = itree + 1
    fo = open(ofile+"_header.txt",'w')
    fo.write('Features: '+str(len(vars))+'\n')
    fo.write('Predictor: '+spec+'\n')
    fo.write('Action: '+action+'\n')
    fo.write(hdr+'\n')
    fo.write('corr_factor: '+'{0:.6f}'.format(corr)+'\n')
    fo.write('dt_chem: '+'{0:.2f}'.format(myrf.dt_chem)+'\n')
    fo.close()
    print('Header written to '+ofile+'_header.txt')
    fo = open(ofile+"_feature_names.txt",'w')
    fo.write('Features: '+str(len(vars))+'\n')
    for var in vars:
        fo.write(var+'\n')
    fo.close()
    print('Feature names written to '+ofile+'_feature_names.txt')
    np.savetxt(ofile+"_lefts.txt",left,fmt='%06d',header=hdr)
    print('Lefts written to '+ofile+'_lefts.txt')
    np.savetxt(ofile+"_rights.txt",right,fmt='%06d',header=hdr)
    print('Rights written to '+ofile+'_rights.txt')
    np.savetxt(ofile+"_features.txt",feature,fmt='%06d',header=hdr)
    print('Features written to '+ofile+'_features.txt')
    np.savetxt(ofile+"_thresholds.txt",threshold,fmt='%.4e',header=hdr)
    print('Thresholds written to '+ofile+'_thresholds.txt')
    np.savetxt(ofile+"_values.txt",value,fmt='%.4e',header=hdr)
    print('Values written to '+ofile+'_values.txt')
    return

def plot_true_vs_pred(forest,farr,parr,fvars,ispec,outtype,ofile,plot_conc):
    '''
    Helper routine to plot prediction vs. true values for a given forest and features and prediction arrays.
    '''
    # compute predictions
    pred = forest.predict(farr)
    # convert units
    true = parr * 28.96/48.0/1.0e3
    pred = pred * 28.96/48.0/1.0e3
    # eventually add tendency to concentration
    if outtype == 'tend' and plot_conc:
        idx  = fvars.index('KPP_AIRDEN')
        idx2 = fvars.index('KPP_BEFORE_CHEM_'+ispec)
        conc = farr[:,idx2]/farr[:,idx] * 28.96 / 48.0 * 1.0e-6
        true[:,0] = conc + true[:,0]
        pred[:] = conc + pred[:]
    # statistics
    R2    = r2_score(true[:,0],pred[:])
    nrmse = sqrt(mean_squared_error(true[:,0],pred[:]))/np.std(true[:,0])
    nmb   = np.sum(pred[:]-true[:,0])/np.sum(true[:,0])
    slope, intercept, r_value, p_value, std_err = stats.linregress(true[:,0],pred[:])
    # scatter plot
    if outtype == 'tend':
        if plot_conc:
            title = ispec+'; concentration + tendency'
            xlabel='true concentration + tendency [ppbv]'
            ylabel='predicted concentration + tendency [ppbv]'
        else: 
            title = ispec+'; tendency'
            xlabel='true tendency [ppbv]'
            ylabel='predicted tendency [ppbv]'
    if outtype=='conc' or outtype=='drop':
        title = ispec+'; concentration'
        xlabel='true concentration [ppbv]'
        ylabel='predicted concentration [ppbv]'
    if outtype=='logconc':
        title = ispec+'; log concentration'
        xlabel='log of true concentration [log10( conc / molec cm-3 )]'
        ylabel='log of predicted concentration [log10( conc / molec cm-3 )]'
    if outtype=='ratio':
        title = ispec+'; ratio'
        xlabel='true concentration ratio [-]'
        ylabel='predicted concentration ratio [-]'
    fig, ax = plt.subplots()
    ax.hexbin(true[:,0],pred[:],cmap=plt.cm.gist_earth_r,bins='log',gridsize=300)
    # 1:1 line
    minval = np.min((np.min(true),np.min(pred)))
    maxval = np.max((np.max(true),np.max(pred)))
    ax.set_xlim(minval,maxval)
    ax.set_ylim(minval,maxval)
    ax.plot((0.95*minval,1.05*maxval),(0.95*minval,1.05*maxval),color='grey',linestyle='dashed')
    # regression line
    ax.plot((0.95*minval,1.05*maxval),(intercept+(0.95*minval*slope),intercept+(1.05*maxval*slope)),color='blue',linestyle='dashed')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    str_ncell = '{:,}'.format(pred.shape[0])
    istr = ' '.join(('N =',str_ncell))
    ax.text(0.05,0.95,istr,transform=ax.transAxes)
    str_R2    = '{0:.2f}'.format(R2)
    istr = ' '.join(('R$^{2}$=',str_R2))
    ax.text(0.05,0.90,istr,transform=ax.transAxes)
    str_rmse  = '{0:.2f}'.format(nrmse*100)
    istr = ' '.join(('NRMSE [%] =',str_rmse))
    ax.text(0.05,0.85,istr,transform=ax.transAxes)
    str_nmb  = '{0:.2f}'.format(nmb*100)
    istr = ' '.join(('NMB [%] =',str_nmb))
    ax.text(0.05,0.80,istr,transform=ax.transAxes)
    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(ofile)
    plt.close()
    print('figure written to '+ofile)
    return

def plot_features(forest,ispec,fvars,ofile,nbars=20):
    '''
    Helper routine to plot importance of input features for a given forest 
    '''
    importances = forest.feature_importances_
    std = np.std([tree.feature_importances_ for tree in forest.estimators_],axis=0)
    indices = np.argsort(importances)[::-1]
    plt.figure(1)
    plt.title(ispec+" - feature importances")
    idx = [indices[i] for i in range(nbars)]
    fvs = [fvars[i] for i in idx]
    for ivar in range(len(fvs)):
        fvs[ivar] = fvs[ivar].replace('KPP_','')
        fvs[ivar] = fvs[ivar].replace('BEFORE_CHEM_','')
    pos = np.arange(nbars)
    idx = idx[::-1]
    plt.barh(pos,importances[idx]*100.0,xerr=std[idx]*100.0,align="center",color=plt.cm.jet(1.*pos/float(nbars)))
    plt.yticks(pos,fvs[::-1])
    plt.ylim([-1,nbars])
    plt.xlabel('Feature importance [%]')
    plt.tight_layout()
    plt.savefig(ofile)
    plt.close()
    print('figure written to '+ofile)
    sys.stdout.flush()
    return
 
def parse_args():
    '''
    Parser for rf_train 
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-t','--template',type=str,default="data/nc4/c48-ml.tavg1_3d_kpp_Nv.20130702_0000z.nc4")
    parser.add_argument('-n','--specnumber',type=int,help='species number',default=41)
    parser.add_argument('-o','--outtype', type=str,help='output type',default='tend')
    parser.add_argument('-tr','--treenumber',type=int,help='tree number',default=0)
    parser.add_argument('-b','--nbins',type=int,help='number of bins',default=10)
    parser.add_argument('-f','--fraction',  type=float,help='fraction of data to use',default=0.1)
    parser.add_argument('-st','--source-training', type=str,help='data source for training',default="data/pkl/training/*.pkl")
    parser.add_argument('-sv','--source-validation', type=str,help='data source for validation',default="data/pkl/validation/*.pkl")
    parser.add_argument('-nl','--nleaves', type=int,help='number of leaves',default=10000)
    parser.add_argument('-nt','--ntrees', type=int,help='number of trees',default=10)
    parser.add_argument('-tp','--treepath', type=str,help='path to save tree file',default='trees')
    return parser.parse_args() 

if __name__ == '__main__':
    rf_train(parse_args())

# eof
