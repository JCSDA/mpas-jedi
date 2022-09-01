import PlotModelSpreadArgs
import basic_plot_functions as bpf
import binning_utils as bu
import config as conf
from copy import deepcopy
from datetime import datetime
import datetime as dt
import logging
import logsetup
import multiprocessing as mp
import numpy as np
import modelsp_utils as mu
import plot_utils as pu
import predefined_configs as pconf
import var_utils as vu

_logger = logging.getLogger(__name__)

class PlotModelSpread():
  '''
  Diagnose model-space statistics
  Driven by
    - static selections in conf
    - command-line arguments in PlotModelSpreadArgs
  '''
  def __init__(self):
    self.name = 'PlotModelSpread'
    self.args = PlotModelSpreadArgs.args
    self.logger = logging.getLogger(self.name)
    #self.nprocs = min(mp.cpu_count(), self.args.nprocs)

    # construct mean DB into 0th member slot
    date = str(self.args.date)
    initDate = datetime.strptime(date,'%Y%m%d%H')
    self.fileDate= initDate.strftime('%Y-%m-%d_%H.%M.%S')

    # mean/deterministic states
    #TODO: allow for entire ReferenceStateFile and MeanStateFile as args instead of date
    self.ReferenceStateFile = str(self.args.referenceState)+'.'+self.fileDate+'.nc'
    self.MeanStateFile = str(self.args.meanState)+'.'+self.fileDate+'.nc'

    self.bgEnsemblePath = (str(self.args.backgroundEnsemblePath)+'.'+self.fileDate+'.nc').replace('YYYYMMDDHH', self.args.date)
    self.anEnsemblePath = (str(self.args.analysisEnsemblePath)+'.'+self.fileDate+'.nc').replace('YYYYMMDDHH', self.args.date)
    self.inflatedEnsemblePath = (str(self.args.inflatedEnsemblePath)+'.'+self.fileDate+'.nc').replace('YYYYMMDDHH', self.args.date)

    self.logger.info('ReferenceState: '+self.ReferenceStateFile)

  def diagnose(self):
    def RMS(x):
      return np.sqrt(np.nanmean(np.square(x)))

    self.logger.info('Populating binning metadata')

    # initialize netcdf4 Dataset's
    ReferenceState = mu.getNCData(self.ReferenceStateFile)

    #TODO: for when MeanState is relevant
    #self.logger.info('MeanState: '+self.MeanStateFile)
    #MeanState = mu.getNCData(self.MeanStateFile)

    MeanStateFile = 'bg/mem001/bg'+'.'+self.fileDate+'.nc'
    self.logger.info('MeanState: '+MeanStateFile)
    MeanState = mu.getNCData(MeanStateFile)

    # bg ens member states
    bgEnsemble = []
    for member in list(range(1, self.args.nMembers+1)):
      ensembleFile = str(self.bgEnsemblePath).format(member)
      self.logger.info('adding background member file: '+ensembleFile)
      bgEnsemble.append(mu.getNCData(ensembleFile))

    # an ens member states
    anEnsemble = []
    for member in list(range(1, self.args.nMembers+1)):
      ensembleFile = str(self.anEnsemblePath).format(member)
      self.logger.info('adding analysis member file: '+ensembleFile)
      anEnsemble.append(mu.getNCData(ensembleFile))

    # an ens member states
    inflatedEnsemble = []
    for member in list(range(1, self.args.nMembers+1)):
      ensembleFile = str(self.inflatedEnsemblePath).format(member)
      self.logger.info('adding inflated member file: '+ensembleFile)
      inflatedEnsemble.append(mu.getNCData(ensembleFile))

    # read "metadata" fields that are used for binning and
    #  create reshaped copies to fit 1D and 2D arrays
    pressure = mu.varRead('pressure', ReferenceState)
    nCells = pressure.shape[0]
    nVertLevelsP1 = pressure.shape[1]+1
    nVertLevels = nVertLevelsP1-1

    dbVals = {}
    for nn in [1, nVertLevels, nVertLevelsP1]:
      dbVals[(nn,)] = {}

    # modVarLev (unique for variables with different numbers of levels)
    for nn in [1, nVertLevels, nVertLevelsP1]:
      dbVals[(nn,)][vu.modVarLev] = np.arange(1, nn+1)

    # modVarLat, modVarLon (same for all level-distributions)
    grid = mu.readGrid(gridFile=self.ReferenceStateFile)
    for nn in [1, nVertLevels, nVertLevelsP1]:
      dbVals[(nn,)][vu.modVarLat] = grid['latitude']
      dbVals[(nn,)][vu.modVarLon] = grid['longitude']

    # add pressure as a binning variable
    # TODO: enable binning by altitude
    # TODO: enable binning by pressure
    # TODO: enable pressure on W levels, if desired
    #dbVals[(nVertLevels,)][vu.modVarPrs] = np.array(pressure)
    #dbVals[(1,)][vu.modVarPrs] = np.array(pressure[:,0]).flatten()


    ###############################################
    ## Extract constructor info about the DiagSpace
    ###############################################

    DiagSpaceName = 'mpas'
    DiagSpaceInfo = conf.DiagSpaceConfig[DiagSpaceName]
    DiagSpaceGrp = DiagSpaceInfo['DiagSpaceGrp']
    binVarConfigs = DiagSpaceInfo.get('binVarConfigs',{})
    selectDiagNames = DiagSpaceInfo.get('diagNames',{})

    modelVars = vu.modVarNames2d+vu.modVarNamesBase3d

    ######################################
    ## Collect statistics for all modelVars
    ######################################

    nsubplots = len(modelVars)
    nxplots = np.int(np.ceil(np.sqrt(nsubplots)))
    nyplots = np.int(np.ceil(np.true_divide(nsubplots, nxplots)))
    subplotWidth = 1.9
    subplotAspect = 0.75

    srFig = pu.setup_fig(nxplots, nyplots, subplotWidth, subplotAspect)
    spreadFig = pu.setup_fig(nxplots, nyplots, subplotWidth, subplotAspect)

    fieldsDB = {
      'ReferenceState': ReferenceState,
      'MeanState': MeanState,
      'bgEnsemble': bgEnsemble,
      'anEnsemble': anEnsemble,
      'infEnsemble': inflatedEnsemble,
    }

    subStats = []
    iplot = 0

    self.logger.info('Calculating then plotting ensemble spread')

    for varName in modelVars: 
      sigmax = {}
      for ee, ensembleState in enumerate(['b', 'a', 'inf']):
        diagName = 'sigmax'+ensembleState

        modelVarName = mu.aggModelVariable(varName)

        diagFunction = mu.diagnosticFunctions[diagName]
        diagnostic = diagFunction.evaluate(varName, fieldsDB)

        dShape = diagnostic.shape
        nDims = len(dShape)
        assert (nDims==1 or nDims==2), 'Number of dimensions must be 1 or 2'

        nn = 1
        if nDims==2:
          nn = dShape[1]

        sigmax[ensembleState] = np.full(nn, np.NaN)

        for lev in np.arange(0, mu.aggMaxLevel(varName, nn)):
          if nDims==2:
            sigmax[ensembleState][lev] = RMS(diagnostic[:,lev])
          else:
            sigmax[ensembleState][lev] = RMS(diagnostic)

      SpreadRatios = []
      self.logger.info('field = '+varName)

      variational = np.divide(sigmax['a'], sigmax['b'])
      SpreadRatios.append(variational)

      RTPP = np.divide(sigmax['inf'], sigmax['a'])
      SpreadRatios.append(RTPP)

      Total = np.divide(sigmax['inf'], sigmax['b'])
      SpreadRatios.append(Total)

      levels = np.arange(1, nn+1)

      marker = None
      if nn == 1: marker = 'o'
      bpf.plotProfile(
        spreadFig,
        [sigmax['b'], sigmax['a'], sigmax['inf']], levels,
        ['bg', 'an', 'relaxed'],
        varName, 'model level', 'Ensemble Spread ($\sigma_x$)',
        bpf.defaultIndepConfig,
        False, False, None,
        nyplots, nxplots, nsubplots, iplot,
        marker=marker)

      bpf.plotProfile(
        srFig,
        SpreadRatios, levels,
        ['EDA', 'RTPP', 'Net'],
        varName, 'model level', 'Spread Ratio',
        bpf.defaultIndepConfig,
        False, True, 1.0,
        nyplots, nxplots, nsubplots, iplot,
        dmin=0.5, dmax=2.0,
        marker=marker)


      iplot += 1

    filename = 'sigmax_'+self.fileDate
    pu.finalize_fig(spreadFig, str(filename), 'pdf', True, 0.6)

    filename = 'SpreadRatio_'+self.fileDate
    pu.finalize_fig(srFig, str(filename), 'pdf', True, 0.6)

    ReferenceState.close()
    MeanState.close()
    for f in bgEnsemble+anEnsemble+inflatedEnsemble:
      f.close()

#=========================================================================
# main program
def main():

  statistics = PlotModelSpread()
  statistics.diagnose()

  _logger.info('Finished '+__name__+' successfully')

if __name__ == '__main__': main()
