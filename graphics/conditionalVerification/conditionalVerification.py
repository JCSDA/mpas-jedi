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
  - cloud fraction threshold
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
      step = args.fcstFrequency
      self.hr = args.forecastHours
      self.doHydro = args.doHydro
      self.domain = args.domain
      self.plotScatter = args.scatter
      self.cldfraThresh = args.cldfraThresh
      self.meshRes = args.meshRes

      self.Nhours = range(0,self.hr+1,step)
      self.dateList, self.ddList = ut.dateList(dateIni, dateEnd, delta)

      self.layerDefinition = 'ERA5'

      self.workers = ut.available_processors()
      logging.info('Number of cpus available: {}'.format(str(self.workers)))


  def driver(self,binningMethod,binningList):

      self.binningMethod = binningMethod
      self.binningList   = binningList

      # Initialize lists/dictionaries
      self.varsList = ut.varsList_
      if self.doHydro:
        self.varsList = self.varsList + ut.hydrosList
      self.binList = list(itertools.product(self.binningList,self.varsList))
      varDict = ut.varDict3L(self.experiments,self.Nhours,self.binningList,self.varsList)
      countT  = ut.varDict3L(self.experiments,self.Nhours,self.binningList,self.varsList,0)

      # read / process files in parallel
      pool_input  = list(itertools.product(self.dateList,self.experiments))
      if self.doHydro:
        mpworkers = len(pool_input)
      else:
        mpworkers = self.workers
      logging.info('Number of cpus to be used: {}'.format(str(mpworkers)))
      all_results = []

      # create pool
      with Pool(processes=mpworkers, maxtasksperchild=1) as pool:
        for hour in self.Nhours:
          args = [(hour,datestr,exp) for datestr,exp in pool_input]
          all_results.append(pool.map(self.read_input, args, chunksize=1))

      logging.info('Aggregating binned values - %s',self.binningMethod)
      for items in all_results:
        for hr,datestr,exp,data in items:
          for c,var in self.binList:
            varDict[var][exp][hr][c] = np.append( varDict[var][exp][hr][c], data[var][c] )

      # get stats
      logging.info('Computing stats - %s',self.binningMethod)
      continous_stats, countT = self.stats_continuous(varDict)
      categorical_stats = self.stats_categorical(countT)

      # make plots
      logging.info('Making plots - %s',self.binningMethod)
      self.make_linePlots(countT, continous_stats, categorical_stats)
      self.make_pdfPlots(varDict)
      if self.plotScatter:
        self.make_scatterPlots(varDict)


    ###################### END OF DRIVE #########################


  def read_input(self,args):
      hr = args[0]; datestr = args[1]; exp = args[2]

      logging.info('Working on: %s %s %d',exp,datestr,hr)
      datestrptime = datetime.strptime(str(datestr), "%Y%m%d%H")
      nso_date = datestrptime + timedelta(hours=hr)
      suffix = 'abi_g16'
      app = 'hofx'
      obsoutfile_ = 'obsout_'+app+'_'+suffix+'.h5'
      geovalfile_ = 'geoval_'+app+'_'+suffix+'.nc4'

      nCellsVMesh, meshRatio = ut.get_meshSpec(self.meshRes)
      static_file = 'x'+str(meshRatio)+'.'+str(nCellsVMesh)+'.static.nc'


      fcst_path = self.main_dir+exp+'/ExtendedFC/'
      veri_path = self.main_dir+exp+'/Verification/fc/mean/'

      obsoutfile  = veri_path+datestr+'/'+str(hr)+'hr/dbOut/'+ obsoutfile_
      geovalfile = veri_path+datestr+'/'+str(hr)+'hr/dbOut/'+ geovalfile_
      staticfile  = self.static_path+static_file

      so_date = nso_date.strftime("%Y-%m-%d_%H")
      so_mpas_10kmL2file = 'saca_obs.'+so_date+'.00.00.nc'
      sofile = self.ctp_path+so_mpas_10kmL2file

      binnedData = self.get_binned_data(datestr, obsoutfile, geovalfile, sofile, staticfile)

      return hr, datestr, exp, binnedData


  def get_binned_data(self, datestr, obsoutfile, geovalsfile, sofile, staticfile):

      # Initialize variable/dictionaries
      binnedData = ut.varDict_binned(self.binningList,self.varsList)
      varsDict = dict.fromkeys(self.varsList)
      varsDictClean = {}
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

          obs = ut.get_obsvalue(obsvalue,ch,obsIds)
          bkg = ut.get_hofx(hofx,ch,obsIds)
          ocldfrac = ut.get_cldfrac(metadata,obsIds)
          fcldfrac3d, pressure, psfc = ut.get_geoval(geovalnc,obsIds)
          octp = ut.get_ctp(obsctpnc,matchIds)
          lat, lon = ut.get_latlon(metadata,obsIds)

          varsDict = {'obs':obs, 'bkg':bkg, 'ocldfrac':ocldfrac, 'fcldfrac3d':fcldfrac3d,
                     'pressure':pressure, 'psfc':psfc, 'octp':octp, 'lat':lat, 'lon':lon}
          if self.doHydro:
            for hy in ut.hydrosList:
              varsDict[hy] = ut.get_hydrometeors(geovalnc,obsIds,hy)

          # clean data
          varsDictClean = ut.cleanData(self.binningMethod, varsDict)
          lenghtList = [len(varsDictClean[k]) for  k in varsDictClean.keys()]
          identical = ut.check(lenghtList)

          if identical: #if all variables have the same lenght

            if self.layerDefinition == 'ERA5':
              PTOP_LOW,PTOP_MID,PTOP_HIGH = ut.get_layerDefinitions_era5(varsDictClean['psfc'])
            elif self.layerDefinition == 'UPP':
              PTOP_LOW,PTOP_MID,PTOP_HIGH = ut.get_layerDefinitions_upp(varsDictClean['psfc'])

            fcstLow  = ut.get_FcstCloudFrac_low(varsDictClean['fcldfrac3d'], varsDictClean['pressure'], varsDictClean['psfc'], PTOP_LOW)
            fcstMid  = ut.get_FcstCloudFrac_mid(varsDictClean['fcldfrac3d'], varsDictClean['pressure'], varsDictClean['psfc'], PTOP_LOW, PTOP_MID)
            fcstHigh = ut.get_FcstCloudFrac_high(varsDictClean['fcldfrac3d'],varsDictClean['pressure'], varsDictClean['psfc'], PTOP_MID,PTOP_HIGH)

            tmp = np.vstack( (fcstLow,fcstMid,fcstHigh) ) # stack the rows into one 2d array from python_stuff.py
            varsDictClean['fcldfrac'] = np.max(tmp,axis=0) # get maximum value across low/mid/high for each pixel (minimum overlap assumption) from python_stuff.py

            bins_ids = ut.get_binsIndexes(self.binningMethod, self.cldfraThresh, varsDictClean['ocldfrac'], varsDictClean['fcldfrac'], fcstLow, fcstMid, fcstHigh, varsDictClean['octp'], PTOP_LOW, PTOP_MID, PTOP_HIGH)

            for c,var in self.binList:
              if len(bins_ids[c]) > 0:
                binnedData[var].update({c: ut.varbinning(bins_ids, c, varsDictClean[var])})

            # close and clean
            obsvalue.close()
            hofx.close()
            metadata.close()
            geovalnc.close()
            obsctpnc.close()
            staticnc.close()

        except OSError as error:
          print(error)

      del varsDict, varsDictClean

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
        countT[exp][hr][c] = ut.count_non_nan(obs_)

        ombDict[stat][exp][c] = np.append( ombDict[stat][exp][c], ut.calc_stat(stat,omb_,obs_) )

        if stat in ut.cont_stats[:2]:
          obsDict[stat][exp][c] = np.append( obsDict[stat][exp][c], ut.calc_stat(stat,obs_,obs_) )
          bkgDict[stat][exp][c] = np.append( bkgDict[stat][exp][c], ut.calc_stat(stat,bkg_,obs_) )

          if self.doHydro:
            for hy in ut.hydrosList:
              hydrosDict[hy][stat][exp][c] = np.append( hydrosDict[hy][stat][exp][c], ut.calc_stat(stat,varDict[hy][exp][hr][c],obs_))

      continous_stats = {'obs': obsDict, 'bkg': bkgDict, 'omb': ombDict}
      if self.doHydro:
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
        for st in range(len(ut.cat_stats)):
          stat = ut.cat_stats[st]
          if self.binningMethod == 'CloudFraction':
            categorical_stats[exp][stat]['all'] = np.append( categorical_stats[exp][stat]['all'], ut.contigency_table_2by2(countT[exp][hr],stat) )
          else:
            if stat in ut.cat_stats_reduced:
              stat_4by4 = ut.contigency_table_4by4(countT[exp][hr], stat)
              for layer in ut.layering_list:
                categorical_stats[exp][stat][layer] = np.append( categorical_stats[exp][stat][layer], stat_4by4[layer] )
      return categorical_stats


  def make_linePlots(self,countT,continous_stats,categorical_stats):

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
          if self.doHydro:
            for hy in ut.hydrosList:
              plot.line(statContinuous, plotList, self.hr, self.experiments, data1=continous_stats[hy][statContinuous], labelsL=self.labelL, savename=statContinuous+'_'+hy+'_full-table')
              plot.line(statContinuous, ut.criteria_3p, self.hr, self.experiments, data1=continous_stats[hy][statContinuous], labelsL=self.labelL, savename=statContinuous+'_'+hy+'_cloudy-clear-all')

      # Categorical
      for stat,layer in statLayer:
        plot.line(stat, [layer], self.hr, self.experiments, data1=categorical_stats, labelsL=self.labelL, savename=stat+'_'+layer)


  def make_scatterPlots(self,varDict):
      ilists     = list(itertools.product(self.experiments,self.Nhours,self.binningList))
      for exp,hr,c in ilists:
        if hr <=3:
          plot.scatter(varDict['lon'][exp][hr][c], varDict['lat'][exp][hr][c], varDict['octp'][exp][hr][c], 'Obs CTP '+str(hr)+'hr '+c,'jet',savename='scatter_octp_'+str(hr)+'hr_'+c+'_'+exp)
          plot.scatter(varDict['lon'][exp][hr][c], varDict['lat'][exp][hr][c], varDict['ocldfrac'][exp][hr][c], 'Obs cldfrac '+str(hr)+'hr '+c,'jet',savename='scatter_ocldfrac_'+str(hr)+'hr_'+c+'_'+exp)


  def make_pdfPlots(self,varDict):
      ilists     = list(itertools.product(self.experiments,self.Nhours,self.binningList))
      for exp,hr,c in ilists:
        # plot PDF
        if c in ut.criteria_only and hr <= 3:
          bkg_ = varDict['bkg'][exp][hr][c]
          obs_ = varDict['obs'][exp][hr][c]
          if ut.count_non_nan(obs_) > 0 and ut.count_non_nan(bkg_) > 0:
            plot.pdf(bkg_,obs_,title=str(hr)+'hr '+c,savename='PDF_'+str(hr)+'hr_'+c+'_'+exp)

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
        '-freq', '--fcstFrequency',
        metavar="HH",
        help="interval between forecast times (output frequency)",
        type=int, required=True)
    required.add_argument(
        '-b', '--binningMethod',
        help="binning method",
        type=str, required=True)
    required.add_argument(
        '-mres', '--meshRes',
        help="mesh resolution",
        type=str, required=True)
    optional = parser.add_argument_group(title='optional arguments')
    optional.add_argument(
        '-doHy', '--doHydro',
        help="whether to do conditional binning for hydrometeors",
        type=ut.str2bool, default=False)
    optional.add_argument(
        '-r', '--domain',
        help="area of interest (fullDisk, finer, coarser)",
        type=str, default="fullDisk")
    optional.add_argument(
        '-scatter', '--scatter',
        help="optionally make scatter plots",
        type=ut.str2bool, default=False)
    optional.add_argument(
        '-thres', '--cldfraThresh',
        help="cloud fraction threshold for clear/cloudy conditions",
        type=float, default=20.0)

    args = parser.parse_args()

    main(args)
