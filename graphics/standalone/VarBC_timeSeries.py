#################################################################################################################################
#
## Script to plot the time series of VarBC predictor coefficients and coefficients errors
#
# Usage:
#
# python VarBC_timeSeries.py -exp EXP -expName shortEXPNAME -dateIni DATEINI -dateEnd DATEEND -mhsType 'ncdiag' -removeEmiss True
#################################################################################################################################

import numpy as np
import pandas as pd
import time, os, argparse
import h5py as h5
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from datetime import datetime, timedelta
from pathlib import Path
from VarBC_dict import VarBCDict
import math
from concurrent.futures import ProcessPoolExecutor

def plot(date_list, sensor, data, expname, typ, title, channel, nobs):

    time = pd.to_datetime(date_list, format='%Y%m%d%H')

    fig = plt.figure()
    gs = GridSpec(ncols=1, nrows=1, height_ratios=[1], wspace=0, hspace=0)
    ax1 = fig.add_subplot(gs[0])

    predictors = list(data.keys())

    ax1.set_title(title+'\n', loc='center', fontweight='bold')
    ax1.set_title('Sensor: '+sensor, loc='left')
    ax1.set_title('nObsUsed @last cycle='+str(nobs), loc='right')
    ax1.xaxis.grid(linestyle=":", alpha=0.2, color='grey')
    ax1.yaxis.grid(linestyle=":", alpha=0.2, color='grey')

    ax1.set_xticks(time[::8])  # every 8 analysis (2 days)
    ax1.set_xticklabels([date.strftime('%d') for date in time[::8]])

    ax1.set_ylabel('Standard Deviation' if typ == 'predcov' else 'Beta', fontsize=12)

    for pred in predictors:
        if len(data[pred]) > 0:
          ax1.plot(time, data[pred], ls='-', label=pred)

    ax1.legend(loc="best", framealpha=0, ncol=1)
    ax1.text(x=0.94, y=-0.12, s=time[-1].strftime('%Y-%b'), transform=ax1.transAxes)

    folder_name = 'predictors'
    Path(folder_name).mkdir(parents=True, exist_ok=True)
    fig.savefig(f'{folder_name}/{expname}_{typ}_{sensor}_{channel}.png', dpi=300, format='png', bbox_inches='tight')
    plt.close(fig)

def main(args_tuple):
    h0 = time.time()

    main_path, date_list, exp, expname, sensorSat, prefix, mhsType, removeEmiss = args_tuple

    nchans = sensorSat['channels']
    chans = sensorSat['analyzed channels']
    n_dates = len(date_list)

    for c in nchans:
        if c not in chans:
            continue
        print(sensorSat['name'], 'channel: ', c)

        ind = chans.index(c)
        ss = sensorSat['name']
        data_coeff, data_cov, predlist = {}, {}, None

        for date_idx, datestr in enumerate(date_list):
            datadir = main_path+exp+'/CyclingDA/'+datestr+'/dbOut/'
            ss = sensorSat['name']
            satbias_file = 'satbias_'+ss+'.h5'
            satbias_cov_file = 'satbias_cov_'+ss+'.h5'
            obsout = 'obsout_da_'+prefix+'.h5'

            if not all(os.path.exists(os.path.join(datadir, f)) for f in [obsout, satbias_file, satbias_cov_file]):
                continue

            with h5.File(os.path.join(datadir, satbias_file), "r") as f_coeff, \
                 h5.File(os.path.join(datadir, satbias_cov_file), "r") as f_cov:
                coeff = f_coeff['BiasCoefficients']
                cov = f_cov['BiasCoefficientErrors']
                nobs_ = f_cov['numberObservationsUsed']

                if predlist is None:
                    predlist = [k for k in coeff.keys() if not removeEmiss or k != 'emissivityJacobian']
                    data_coeff = {pred: np.full(n_dates, np.nan) for pred in predlist}
                    data_cov = {pred: np.full(n_dates, np.nan) for pred in predlist}

                for pred in predlist:
                    # this is because we specify channels differently for ABI/AHI/MHS-ncdiag (a subset)
                    if (prefix.startswith('abi') or prefix.startswith('ahi') or
                        (prefix.startswith('mhs') and mhsType == 'ncdiag')):
                        data_coeff[pred][date_idx] = coeff[pred][0][ind]
                        data_cov[pred][date_idx] = math.sqrt(cov[pred][0][ind])
                        nobs = nobs_[0][ind]
                    else:
                        # here c-1 to extract the correct predictor value for the specific channel
                        # as indices in python start at zero
                        data_coeff[pred][date_idx] = coeff[pred][0][c-1]
                        data_cov[pred][date_idx] = math.sqrt(cov[pred][0][c-1])
                        nobs = nobs_[0][c-1]

        if predlist and any(~np.isnan(data_coeff[pred]).all() for pred in predlist):
           plot(date_list, ss, data_coeff, expname, 'predcoeff', 'Bias coefficients @ch'+str(c), c, nobs)
           plot(date_list, ss, data_cov, expname, 'predcov', 'Bias coefficients errors @ch'+str(c), c, nobs)

    print(f'[DONE] {prefix} in {time.time() - h0:.2f} seconds')

if __name__ == '__main__':
    h = time.time()

    parser = argparse.ArgumentParser(description='Plot time series of VarBC predictor coefficients and coefficients errors', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-exp', '--exp', type=str, help='Experiment name', required=True)
    parser.add_argument('-expName', '--expName', type=str, help='Experiment short name for convenience', required=True)
    parser.add_argument('-dateIni', '--dateIni', type=str, help='First analysis time', required=True)
    parser.add_argument('-dateEnd', '--dateEnd', type=str, help='Last analysis time', required=True)
    parser.add_argument('-mhsType',  '--mhsType', type=str, help='MHS data type: raw or ncdiag', default='raw')
    parser.add_argument('-removeEmiss', '--removeEmiss', type=bool, help='Remove emissivity Jacobian from predictors list', default=False)
    args = parser.parse_args()

    user = os.environ['USER']
    main_path = '/glade/derecho/scratch/'+user+'/pandac/'

    delta = 6
    start = datetime.strptime(args.dateIni, "%Y%m%d%H")
    end = datetime.strptime(args.dateEnd, "%Y%m%d%H")
    date_list = [(start + timedelta(hours=i)).strftime("%Y%m%d%H")
                 for i in range(0, int((end - start).total_seconds() // 3600) + 1, delta)]

    tasks = []
    for prefix in VarBCDict:
        tasks.append((main_path, date_list, args.exp, args.expName, VarBCDict[prefix], prefix, args.mhsType, args.removeEmiss))

    with ProcessPoolExecutor() as executor:
        executor.map(main, tasks)

    print(f'Total time elapsed: {time.time() - h:.2f} seconds')
