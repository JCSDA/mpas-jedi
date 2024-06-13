import sys
import argparse
from datetime import datetime, timedelta
import multiprocessing
from multiprocessing import Pool
import numpy as np
import os, time, itertools
from itertools import repeat
import logging
import plot
import utils as ut
ut.logg()

class ConditionalVerification:
  '''
  Conditional brigthness temperature forecasts verification by cloud fraction and cloud layering against ABI products
  Optional:
  - conditional hydrometeors verification
  - domain specification
    default: fullDisk
    options: finer (refinement area), coaser (outside refinement area)
  '''
  def __init__(self,args):
      self.name = 'ConditionalVerification'
      logging.info(self.name)

      self.main_dir = args.rootDir
      self.ctp_path = args.ctpDir
      self.static_path = args.staticDir
      self.experiments = args.experiments
      self.labelL = args.expnames
      dateIni = args.dateIni
      dateEnd = args.dateEnd
      delta = args.delta
      self.hr = args.forecastHours
      self.hydro = args.hydro
      self.domain = args.domain

      self.Nhours = range(self.hr+1)
      self.dateList, self.ddList = ut.dateList(dateIni, dateEnd, delta)

      self.layerDefinition = 'ERA5'

      self.workers = ut.available_processors()
      logging.info('Number of cpus available: {}'.format(str(self.workers)))


  def driver(self,binningMethod,binningList):

      self.binningMethod = binningMethod
      self.binningList   = binningList

      # Initialize lists/dictionaries
      varDict = ut.varDict3L(self.experiments,self.Nhours,self.binningList)
      countT  = ut.varDict3L(self.experiments,self.Nhours,self.binningList,0)
      self.binList = list(itertools.product(self.binningList,ut.varsList))

      # read / process files in parallel
      pool_input  = list(itertools.product(self.dateList,self.experiments))
      mworkers = len(pool_input)

      # create a thread pool
      for hour in self.Nhours:
        with Pool(processes=self.workers) as pool:
          args = [(hour,datestr,exp) for datestr,exp in pool_input]
          items_dataDict = pool.map(self.read_input, args)

          logging.info('Aggregating binned values - %s',self.binningMethod)
          for hr,datestr,exp,data in items_dataDict:
            for c,var in self.binList:
              varDict[var][exp][hr][c] = np.append( varDict[var][exp][hr][c], data[var][c] )

        pool.close()
        pool.join()

      # get stats
      logging.info('Computing stats - %s',self.binningMethod)
      continous_stats, countT = self.stats_continuous(varDict)
      categorical_stats = self.stats_categorical(countT)

      # make plots
      logging.info('Making plots - %s',self.binningMethod)
      self.make_plots(countT, continous_stats, categorical_stats)

    ###################### END OF DRIVE #########################


  def read_input(self,args):
      hr = args[0]; datestr = args[1]; exp = args[2]

      logging.info('Working on: %s %s %d',exp,datestr,hr)
      datestrptime = datetime.strptime(str(datestr), "%Y%m%d%H")
      nso_date = datestrptime + timedelta(hours=hr)
      suffix = 'abi_g16'
      app = 'hofx'
      obsoutfile  = 'obsout_'+app+'_'+suffix+'.h5'
      geovalfile  = 'geoval_'+app+'_'+suffix+'.nc4'
      nCellsVMesh = 835586
      static_file = 'x20.'+str(nCellsVMesh)+'.static.nc'

      fcst_path = self.main_dir+exp+'/ExtendedFC/'
      veri_path = self.main_dir+exp+'/Verification/fc_abiOnly_model/mean/'

      obsoutfile  = veri_path+datestr+'/'+str(hr)+'hr/dbOut/'+ obsoutfile
      geovalsfile = veri_path+datestr+'/'+str(hr)+'hr/dbOut/'+ geovalfile
      staticfile  = self.static_path+static_file

      so_date = nso_date.strftime("%Y-%m-%d_%H")
      so_mpas_10kmL2file = 'saca_obs.'+so_date+'.00.00.nc'
      sofile = self.ctp_path+so_mpas_10kmL2file

      binnedData = self.get_binned_data(datestr, obsoutfile, geovalsfile, sofile, staticfile)

      return hr, datestr, exp, binnedData


  def get_binned_data(self, datestr, obsoutfile, geovalsfile, sofile, staticfile):

      # Initialize variable/dictionaries
      binnedData = ut.varDict_binned(self.binningList)
      varsDict = {}
      ch = 6
      threshold = 78.0

      if os.path.exists(obsoutfile) and os.path.exists(geovalsfile) and os.path.exists(sofile):
        try:
          obsvalue = ut.open_dataset(obsoutfile, enginetype='h5', group='ObsValue')
          hofx     = ut.open_dataset(obsoutfile, enginetype='h5', group='hofx')
          metadata = ut.open_dataset(obsoutfile, enginetype='h5', group='MetaData')
          geovalnc = ut.open_dataset(geovalsfile, enginetype='nc')
          obsctpnc = ut.open_dataset(sofile, enginetype='nc')
          staticnc = ut.open_dataset(staticfile, enginetype='nc')

          matchIds, obsIds  = ut.getDomainIndexes(threshold,self.domain,metadata,staticnc)

          varsDict['obs'] = ut.get_obsvalue(obsvalue,ch,obsIds)
          varsDict['bkg'] = ut.get_hofx(hofx,ch,obsIds)
          ocldfrac = ut.get_cldfrac(metadata,obsIds)
          fcldfrac, pressure, psfc = ut.get_geoval(geovalnc,obsIds)
          octp = ut.get_ctp(obsctpnc,matchIds)

          if self.hydro:
            varsDict['ihumr'], varsDict['imclw'], varsDict['imcli'], varsDict['imclr'], varsDict['imcls'], varsDict['imclg'] = ut.get_hydrometeors(geovalnc,obsIds)

          if np.count_nonzero(varsDict['obs']) != 0 or (len(varsDict['obs']) != 0) or (len(varsDict['bkg']) != 0) or ut.count_non_nan(ocldfrac) != 0 or ut.count_non_nan(octp) != 0: #if obs exist and ctp have non nan data

            if self.layerDefinition == 'ERA5':
              PTOP_LOW,PTOP_MID,PTOP_HIGH = ut.get_layerDefinitions_era5(psfc)
            elif self.layerDefinition == 'UPP':
              PTOP_LOW,PTOP_MID,PTOP_HIGH = ut.get_layerDefinitions_upp(psfc)

            fcstLow  = ut.get_FcstCloudFrac_low(fcldfrac, pressure, psfc, PTOP_LOW)
            fcstMid  = ut.get_FcstCloudFrac_mid(fcldfrac, pressure, psfc, PTOP_LOW, PTOP_MID)
            fcstHigh = ut.get_FcstCloudFrac_high(fcldfrac, pressure, psfc, PTOP_MID, PTOP_HIGH)

            tmp = np.vstack( (fcstLow,fcstMid,fcstHigh) ) # stack the rows into one 2d array
            fcfrmax = np.max(tmp,axis=0) # get maximum value across low/mid/high for each pixel (minimum overlap assumption)

            bins_ids = ut.get_binsIndexes(self.binningMethod, ocldfrac, fcfrmax, fcstLow, fcstMid, fcstHigh, octp, PTOP_LOW, PTOP_MID, PTOP_HIGH)
            for c,var in self.binList:
              if len(bins_ids[c]) > 0:
                binnedData[var].update({c: ut.varbinning(bins_ids, c, varsDict[var])})
              else:
                continue

            # close and clean
            obsvalue.close()
            hofx.close()
            metadata.close()
            geovalnc.close()
            obsctpnc.close()
            staticnc.close()
            del varsDict, psfc, ocldfrac, fcfrmax, fcstLow, fcstMid, fcstHigh, octp, PTOP_LOW, PTOP_MID, PTOP_HIGH

        except OSError as error:
          print(error)

      return binnedData


  def stats_continuous(self,varDict):
      ##########################################################
      # Calculate continuous stats (mean, count, bias, rmse, std)
      ##########################################################

      # Initialize lists/dictionaries
      countT     = ut.createNestedDictT(self.experiments,self.Nhours,self.binningList)
      obsDict    = ut.varContStatDict(self.experiments,self.binningList)
      bkgDict    = ut.varContStatDict(self.experiments,self.binningList)
      ombDict    = ut.varContStatDict(self.experiments,self.binningList)
      hydrosDict = ut.varContHydroStatDict(self.experiments,self.binningList)
      ilists     = list(itertools.product(ut.cont_stats,self.experiments,self.Nhours,self.binningList))

      for stat,exp,hr,c in ilists:
        obs_ = varDict['obs'][exp][hr][c]
        bkg_ = varDict['bkg'][exp][hr][c]
        omb_ = ut.omb(obs_,bkg_)
        countT[exp][hr][c] = len(omb_)

        ombDict[stat][exp][c] = np.append( ombDict[stat][exp][c], ut.calc_stat(stat,omb_) )

        if stat in ut.cont_stats[:2]:
          obsDict[stat][exp][c] = np.append( obsDict[stat][exp][c], ut.calc_stat(stat,obs_,omb_) )
          bkgDict[stat][exp][c] = np.append( bkgDict[stat][exp][c], ut.calc_stat(stat,bkg_,omb_) )

          if self.hydro:
            for hy in ut.hydrosList:
              hydrosDict[hy][stat][exp][c] = np.append( hydrosDict[hy][stat][exp][c], ut.calc_stat(stat,varDict[hy][exp][hr][c],omb_))

        # plot PDF
        if c in ut.criteria_only and hr <= 6:
          if ut.count_non_nan(obs_) > 0 and ut.count_non_nan(bkg_) > 0:
            plot.pdf(bkg_,obs_,title=str(hr)+'hr '+c,savename='PDF_'+str(hr)+'hr_'+c+'_'+exp)

        del obs_, bkg_, omb_

      continous_stats = {'obs': obsDict, 'bkg': bkgDict, 'omb': ombDict}
      if self.hydro:
        continous_stats.update(hydrosDict)

      return continous_stats, countT


  def stats_categorical(self,countT):
      #########################################################
      # Calculate categorical stats (FBIAS, FAR, POD, GSS, etc)
      #########################################################

      # Initialize lists/dictionaries
      categorical_stats = ut.createNestedDictT(self.experiments,ut.cat_stats,ut.layering_list)
      clists = list(itertools.product(self.experiments,self.Nhours))

      for exp,hr in clists:
        if self.binningMethod == 'CloudFraction':
          for st in range(len(ut.cat_stats)):
            stat = ut.cat_stats[st]
            stats_2by2 = ut.contigency_table_2by2(countT[exp][hr])
            categorical_stats[exp][stat]['all'] = np.append( categorical_stats[exp][stat]['all'], stats_2by2[st] )
        else:
          for st in range(len(ut.cat_stats)):
            stat = ut.cat_stats[st]
            if stat in ut.cat_stats_reduced:
              stat_4by4 = ut.contigency_table_4by4(countT[exp][hr], stat)

              for layer in ut.layering_list:
                categorical_stats[exp][stat][layer] = np.append( categorical_stats[exp][stat][layer], stat_4by4[layer] )

      return categorical_stats


  def make_plots(self,countT,continous_stats,categorical_stats):

      plotList, statLayer  = ut.get_list4ploting(self.binningList)

      ############ LINE PLOTS ####################################
      # Continuous
      for statContinuous in ut.cont_stats:
        ## omb
        plot.line(statContinuous, plotList, self.hr, self.experiments, data1=continous_stats['omb'][statContinuous], labelsL=self.labelL, savename=statContinuous+'_omb_full-table')
        plot.line(statContinuous, ut.criteria_3p, self.hr, self.experiments, data1=continous_stats['omb'][statContinuous], labelsL=self.labelL, savename=statContinuous+'_omb_cloudy-clear-all')

        if statContinuous in ut.cont_stats[:2]:
          ## obs, bkg
          plot.line(statContinuous, plotList, self.hr, self.experiments, continous_stats['bkg'][statContinuous], continous_stats['obs'][statContinuous], self.labelL, statContinuous+'_obs-bkg_full-table')
          plot.line(statContinuous, ut.criteria_3p, self.hr, self.experiments, continous_stats['bkg'][statContinuous], continous_stats['obs'][statContinuous], self.labelL, statContinuous+'_obs-bkg_cloudy-clear-all')

          ## hydros
          if self.hydro:
            for hy in ut.hydrosList:
              plot.line(statContinuous, plotList, self.hr, self.experiments, data1=continous_stats[hy][statContinuous], labelsL=self.labelL, savename=statContinuous+'_'+hy+'_full-table')
              plot.line(statContinuous, ut.criteria_3p, self.hr, self.experiments, data1=continous_stats[hy][statContinuous], labelsL=self.labelL, savename=statContinuous+'_'+hy+'_cloudy-clear-all')

      # Categorical
      for stat,layer in statLayer:
        plot.line(stat, [layer], self.hr, self.experiments, data1=categorical_stats, labelsL=self.labelL, savename=stat+'_'+layer)

#=========================================================================
# main program
def main(args):

    h0 = time.time()
    logging.info('Starting '+__name__)

    verification = ConditionalVerification(args)

    binningList = ut.binningDict[args.binningMethod]
    binningMethod = args.binningMethod
    verification.driver(binningMethod,binningList)

    hf = time.time()
    ftime = hf - h0
    logging.info('Finished '+__name__+' successfully')
    logging.info('Time elapsed: {}'.format(str(ftime)))


if __name__ == "__main__":

    # Get command line arguments
    parser = argparse.ArgumentParser(
        description=(
            'Reads input file(s),'
            ' aggregates and bins BT values by binning criteria. '
            ' Compute statistics and makes plots.')
    )
    required = parser.add_argument_group(title='required arguments')
    required.add_argument(
        '-dir', '--rootDir',
        help="experiments root directory",
        type=str, required=True)
    required.add_argument(
        '-ctp', '--ctpDir',
        help="directory for superobbed ABI CTP observations",
        type=str, required=True)
    required.add_argument(
        '-static', '--staticDir',
        help="static file directory",
        type=str, required=True)
    required.add_argument(
        '-exp', '--experiments',
        help="experiments list",
        type=ut.list_of_strings, required=True)
    required.add_argument(
        '-n', '--expnames',
        help="experiments names list",
        type=ut.list_of_strings, required=True)
    required.add_argument(
        '-di', '--dateIni',
        metavar="YYYYMMDDHH",
        help="initial analysis time",
        type=str, required=True)
    required.add_argument(
        '-de', '--dateEnd',
        metavar="YYYYMMDDHH",
        help="final analysis time",
        type=str, required=True)
    required.add_argument(
        '-dt', '--delta',
        metavar="HH",
        help="interval between analysis times",
        type=int, required=True)
    required.add_argument(
        '-fhr', '--forecastHours',
        metavar="HH",
        help="length of the forecast",
        type=int, required=True)
    required.add_argument(
        '-b', '--binningMethod',
        help="binning method",
        type=str, required=True)
    optional = parser.add_argument_group(title='optional arguments')
    optional.add_argument(
        '-hy', '--hydro',
        help="conditional binning for hydrometeors",
        type=bool, default=True)
    optional.add_argument(
        '-r', '--domain',
        help="area of interest (fullDisk, finer, coarser)",
        type=str, default="fullDisk")

    args = parser.parse_args()

    main(args)
