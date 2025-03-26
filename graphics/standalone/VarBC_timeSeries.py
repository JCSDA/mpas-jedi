#############################################################################################
#
## Script to plot the time series of VarBC predictor coefficients and coefficients errors
#
# Usage:
#
# python VarBC_timeSeries.py -exp EXP -expName shortEXPNAME -dateIni DATEINI -dateEnd DATEEND
#############################################################################################

import numpy as np
import pandas as pd
import time, os, argparse
import h5py as h5
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from datetime import datetime, timedelta
from pathlib import Path
from VarBC_dict import VarBCDict
import math

def plot(days, sensor, data, expname, typ, title, channel, nobs):

    time = pd.to_datetime(days, format='%Y-%m-%d %H:%M:%S')

    fig = plt.figure()
    gs = GridSpec(ncols=1, nrows=1, height_ratios=[1], wspace=0, hspace=0)
    ax1 = fig.add_subplot(gs[0])

    predictors = list(data.keys())

    ax1.set_title(title+'\n', loc='center', fontweight='bold')
    ax1.set_title('Sensor: '+sensor+'  nObsUsed @last cycle='+str(nobs) , loc='left')
    ax1.xaxis.grid(linestyle=":", alpha=0.2, color='grey')
    ax1.yaxis.grid(linestyle=":", alpha=0.2, color='grey')

    ax1.set_xticks(time[::8])  # every 8 analysis (2 days)
    ax1.set_xticklabels([date.strftime('%d') for date in time[::8]])

    if (typ == 'predcov'):
      ax1.set_ylabel('Standard Deviation',fontsize=12)
    if (typ == 'predcoeff'):
      ax1.set_ylabel('Beta',fontsize=12)

    for pred in predictors:
      if len(data[pred]) > 0:
        ax1.plot(time, data[pred],  ls='-', label=pred)
        ax1.legend(loc="best",framealpha=0, ncol=1)
        ax1.text(x=0.94, y =-0.12, s=time[-1].strftime('%Y')+'-'+time[-1].strftime('%b'), transform=ax1.transAxes)

    folder_name = 'predictors'
    folder_path = Path(folder_name)
    folder_path.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder_name+'/'+expname+'_'+typ+'_'+sensor+'_'+str(channel)+'.png', dpi=300, format='png', bbox_inches='tight')
    plt.close(fig)

def main(main_path, dateIni, dateEnd, exp, expname, sensorSat, prefix):
    h0 = time.time()

    nchans = sensorSat['channels']
    chans = sensorSat['analyzed channels']

    for c in nchans:
      if c in chans:
        print(sensorSat['name'], 'channel: ', c)
        data_coeff = {v: [] for v in predlist}
        data_cov   = {v: [] for v in predlist}
        delta='6'
        datei = datetime.strptime(str(dateIni), "%Y%m%d%H")
        datef = datetime.strptime(str(dateEnd), "%Y%m%d%H")
        date  = datei
        date_list  = []

        while (date <= datef):
          datestr  = date.strftime("%Y%m%d%H")
          day  = date.strftime("%d")
          date_list.append(date)

          datadir      = main_path+exp+'/CyclingDA/'+datestr+'/dbOut/'
          ss = sensorSat['name']
          satbias_file = 'satbias_'+ss+'.h5'
          satbias_cov_file = 'satbias_cov_'+ss+'.h5'
          obsout = 'obsout_da_'+prefix+'.h5'
          if os.path.exists(datadir+obsout):
            if os.path.exists(datadir+satbias_file) and os.path.exists(datadir+satbias_cov_file):
              coeff = h5.File(datadir+satbias_file, "r")['BiasCoefficients']
              cov   = h5.File(datadir+satbias_cov_file, "r")['BiasCoefficientErrors']
              nobs  = h5.File(datadir+satbias_cov_file, "r")['numberObservationsUsed']

              for fp in range(len(predlist)):
                # this is because we specify channels differently for ABI/AHI/MHS (a subset)
                if ('abi' in prefix) or ('ahi' in prefix) or ('mhs' in prefix):
                  ind = chans.index(c)
                  data_coeff[predlist[fp]] = np.append(data_coeff[predlist[fp]], coeff[predlist[fp]][0][ind])
                  data_cov[predlist[fp]]   = np.append(data_cov[predlist[fp]],   cov[predlist[fp]][0][ind])
                else:
                  data_coeff[predlist[fp]] = np.append(data_coeff[predlist[fp]], coeff[predlist[fp]][0][c-1])
                  data_cov[predlist[fp]]   = np.append(data_cov[predlist[fp]],   math.sqrt(cov[predlist[fp]][0][c-1]))

          date = date + timedelta(hours=int(delta))
        non_empty = all(len(data_coeff.get(fp, [])) > 0 and len(data_cov.get(fp, [])) > 0 for fp in predlist)

        if non_empty:
          plot(date_list, ss, data_coeff, expname, 'predcoeff', 'Bias coefficients @ch'+str(c), c, nobs[0][c-1])
          plot(date_list, ss, data_cov, expname, 'predcov', 'Bias coefficients errors @ch'+str(c), c, nobs[0][c-1])

    hf = time.time()
    print('Time elapsed: ',hf  - h0)


if __name__ == '__main__':

   parser = argparse.ArgumentParser(description='Plot time series of VarBC predictor coefficients and coefficients errors', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   required = parser.add_argument_group(title='required arguments')
   required.add_argument('-exp', '--exp',type=str, help='Experiment name', required=True)
   required.add_argument('-expName', '--expName',type=str, help='Experiment short name for convenience', required=True)
   required.add_argument('-dateIni', '--dateIni',type=str, help='First analysis time', required=True)
   required.add_argument('-dateEnd', '--dateEnd',type=str, help='Last analysis time', required=True)
   args = parser.parse_args()

   exp = args.exp
   shortName = args.expName
   dateIni   = args.dateIni
   dateEnd   = args.dateEnd

   user = os.environ['USER']
   main_path = '/glade/derecho/scratch/'+user+'/pandac/'

   predlist = ['constant', 'emissivityJacobian', 'lapseRate', 'lapseRate_order_2', 'sensorScanAngle', 'sensorScanAngle_order_2', 'sensorScanAngle_order_3', 'sensorScanAngle_order_4']

   obsPrefix = VarBCDict.keys()

   for prefix in obsPrefix:
     main(main_path, dateIni, dateEnd, exp, shortName, VarBCDict[prefix], prefix)
