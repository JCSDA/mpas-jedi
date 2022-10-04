#!/usr/bin/env python3

import GenerateABEIFactorsArgs

from basic_plot_functions import plotDistri
import binning_utils as bu
from distutils.util import strtobool
import predefined_configs as pconf
from collections import defaultdict
from collections.abc import Iterable
import config as conf
from copy import deepcopy
import diag_utils as du
import fnmatch
import glob
from JediDB import JediDB
from JediDBArgs import obsFKey
import Interpolate
import logging
import logsetup
import modelsp_utils as modelUtils
import multiprocessing as mp
from netCDF4 import Dataset
import numpy as np
import os
#from scipy.spatial import cKDTree
import stat_utils as su
import var_utils as vu

_logger = logging.getLogger(__name__)

class GenerateABEIFactors:
  '''
  Diagnose observation-space adaptive background error inflation
  Driven by
    - static selections in conf
    - command-line arguments in GenerateABEIFactorsArgs
  '''
  def __init__(self):
    global argsProcessor
    self.name = 'GenerateABEIFactors'
    self.args = GenerateABEIFactorsArgs.args
    self.logger = logging.getLogger(self.name)
    self.nprocs = min(mp.cpu_count(), self.args.nprocs)

    # observation batch size when updating local model cell inflation values
    # trade-off is between memory overhead and wall-time when self.nprocs>1
    # 5000 works well for MPAS 30km mesh and 15x15 super-ob subgrids
    self.batchSize = 5000

    # localization radius for inflation projection
    self.localizationRadius = self.args.localizationRadius * 1.0e3 #convert from km to meters
    self.cutoffRadius = 3.5 #normalized distance where localization is set to zero

    # construct mean DB into 0th member slot
    self.logger.info('database path: '+self.args.dbPath)
    self.db = JediDB(self.args.dbPath)
    self.channels = self.args.channels.split(',')
    instNames = self.args.IRInstruments.split(',')
    dbOSNames = list(self.db.ObsSpaceName.values())
    self.osKeys = []
    for inst in instNames:
      if inst == '': continue
      message = self.__class__.__name__+'.__init__: '+inst+' not included in JediDB'
      assert inst in dbOSNames, message
      for osKey, dsName in self.db.ObsSpaceName.items():
        if inst == dsName:
          self.osKeys.append(osKey)
          break

    self.plotLambda = bool(strtobool(self.args.plotLambda))

#  def getWorkers(self):
#    if self.nprocs > 1:
#      workers = mp.Pool(processes = self.nprocs)
#    else:
#      workers = None
#    return workers

#  def finalizeWorkers(self, workers):
#    if workers is not None:
#      workers.close()
#      workers.join()

  def diagnose(self):
    '''
    conducts diagnoseObsSpace across multiple ObsSpaces in parallel
    '''
    # Loop over all experiment+observation combinations (keys) alphabetically
    for kk, osKey in enumerate(self.osKeys):
      self.logger.info(osKey)
      #Note: processing must be done serially due to the nature of the inflation
      #      update procedure... could process in parallel across channels/obsVars
      self.diagnoseObsSpace(self.db, osKey, {'resetLambda': kk==0})

  def getObsVars(self, db, osKey):
    # create observed variable list by selecting those variables in the
    # obs feedback files (obsFKey) with the proper suffix
#    if self.args.jediAppName == 'variational':
#        markerSuffix = vu.ombgGroup
#    elif self.args.jediAppName == 'hofx':
#        markerSuffix = vu.hofxGroup
#    else:
#        logger.error('JEDI Application is not supported:: '+self.args.jediAppName)
#
#    return db.varList(osKey, obsFKey, markerSuffix)
    return ['brightness_temperature_'+str(c) for c in self.channels]

  @staticmethod
  def getBinVarConfigs():
    binVarConfigs = defaultdict(list)
    binVarConfigs[vu.obsVarQC] += [bu.passQCMethod]
    return binVarConfigs

  @staticmethod
  def getBinMethods(binVarConfigs, ObsSpaceName, fileFormat):
    binMethods = {}
    for binVarKey, binMethodKeys in binVarConfigs.items():
        binVarConfig = pconf.binVarConfigs.get(binVarKey,pconf.nullBinVarConfig)
        for binMethodKey in binMethodKeys:
            config = binVarConfig.get(binMethodKey,pconf.nullBinMethod).copy()

            if (len(config['values']) < 1 or
                len(config['filters']) < 1): continue

            config['dsName'] = ObsSpaceName
            config['fileFormat'] = fileFormat
            binMethods[(binVarKey,binMethodKey)] = bu.BinMethod(config)

    return binMethods

  @staticmethod
  def getDiagNames():
    return ['ABEILambda']

  @staticmethod
  def getDBVars(obsVars, binMethods, diagnosticConfigs, fileFormat):
    dbVars = []
    for varName in obsVars:
        for diagName, diagnosticConfig in diagnosticConfigs.items():
            if 'BinFunction' not in diagnosticConfig: continue

            # variables for diagnostics
            for grpVar in diagnosticConfig['BinFunction'].dbVars(
                varName, diagnosticConfig['outerIter']):
                if diagnosticConfig[vu.mean]:
                    dbVars.append(grpVar)

            # variables for binning
            for (binVarKey,binMethodKey), binMethod in binMethods.items():
                for grpVar in binMethod.dbVars(
                    varName, diagnosticConfig['outerIter']):
                    dbVars.append(grpVar)

    # manually add latitude and longitude to dbVars
    for var in [vu.latMeta, vu.lonMeta]:
      dbVars.append(vu.base2dbVar(var, var, fileFormat))

    return dbVars


  def diagnoseObsSpace(self, db, osKey, localConfig):
    #  osKey - key of db dict to process
    logger = logging.getLogger(self.name+'.diagnoseObsSpace('+osKey+')')

    # initialize mean db file handles
    db.initHandles(osKey)

    ###############################################
    ## Extract constructor info about the ObsSpace
    ###############################################

    ObsSpaceName  = db.ObsSpaceName[osKey]
    ObsSpaceInfo  = conf.DiagSpaceConfig[ObsSpaceName]
    ObsSpaceGrp   = ObsSpaceInfo['DiagSpaceGrp']

    obsVars = self.getObsVars(db, osKey)
    fileFormat = db.fileFormat(osKey, obsFKey)

    ########################################################
    ## Construct dictionary of binMethods for this ObsSpace
    ########################################################

    binVarConfigs = self.getBinVarConfigs()
    binMethods = self.getBinMethods(binVarConfigs, ObsSpaceName, fileFormat)

    ######################################
    ## Construct diagnostic configurations
    ######################################

    diagnosticConfigs = du.diagnosticConfigs(
        self.getDiagNames(), ObsSpaceName,
        includeEnsembleDiagnostics = False,
        fileFormat = fileFormat)

    ######################################
    ## Generate list of required variables
    ######################################

    dbVars = self.getDBVars(obsVars, binMethods, diagnosticConfigs, fileFormat)


    ##################################
    ## Read required variables from db
    ##################################

    # read mean database variable values into memory
    dbVals = db.readVars(osKey, dbVars)

    # destroy mean file handles
    db.destroyHandles(osKey)


    #################################################
    ## Calculate and plot diagnostics for all obsVars
    ## Calculate and write out inflation factors
    #################################################

    # read model variables used for inflation calculation/templating
    modelGridFile = self.args.modelGridFile
    modelGrid = modelUtils.getNCData(modelGridFile)

    grid = modelUtils.readGrid(gridFile=modelGridFile)
    modelLatsDeg = grid['latitude']
    modelLonsDeg = grid['longitude']
    modelAreas = grid['area']
    EarthSphereR = grid['R']

    self.modelLats = np.multiply(modelLatsDeg, np.pi / 180.0)
    self.modelLons = np.multiply(modelLonsDeg, np.pi / 180.0)
    self.EarthSphereR = EarthSphereR

    modelTemplateVars = {}
    for templateType, templateVar in modelUtils.templateVariables.items():
      modelTemplateVars[templateType] = {
        'values': modelUtils.varRead(templateVar, modelGrid),
        'dims': modelUtils.varDims(templateVar, modelGrid),
      }

    # setup observation lat/lon
    latMeta = vu.base2dbVar(vu.latMeta, vu.latMeta, fileFormat)
    obsLatsDeg = dbVals[latMeta]
    obsLats = np.multiply(obsLatsDeg, np.pi / 180.0)

    lonMeta = vu.base2dbVar(vu.lonMeta, vu.lonMeta, fileFormat)
    obsLonsDeg = dbVals[lonMeta]
    obsLons = np.multiply(obsLonsDeg, np.pi / 180.0)

    obsnLocs = len(obsLons)

    # setup interpolation objects
    logger.info('Setting up interpolators')

    # interpolator from model to observation locations
    model2obs = Interpolate.InterpolateLonLat(
      self.modelLons, self.modelLats,
      #weightMethod = 'unstinterp',
      weightMethod = 'barycentric',
      Radius = EarthSphereR)
    model2obs.initWeights(obsLons, obsLats)

    # interpolator for calculating local distances between model and observation locations
    localRadius = self.cutoffRadius*self.localizationRadius
    localArea = np.pi * localRadius**2

    meanCellArea = modelAreas.mean()
    nLocalCells = np.floor(localArea / meanCellArea * 1.2).astype(int)

    #TODO: improve area strategy for regional or variable resolution meshes
    #minCellArea = modelAreas.min()
    #nLocalCells = np.floor(localArea / minCellArea * 1.2).astype(int)

    logger.info('Maximum number of cells per obs localization radius: '+str(nLocalCells))

    modelLocal2obs = Interpolate.InterpolateLonLat(
      self.modelLons, self.modelLats,
      #weightMethod = 'unstinterp',
      weightMethod = 'barycentric',
      nInterpPoints = nLocalCells, #TODO: use reasonable number of local model points per observation based on mesh area and 3.5x localization radius
      Radius = EarthSphereR)

    obsInds = np.arange(0, len(obsLats))

    # plotting attributes
    # TODO: move these attributes to diag_utils
    diagColors = {}
    diagColors['ABEILambda'] = 'YlOrBr'

    minmaxValue = {}
#    minmaxValue['ABEILambda'] = [
#      bu.ABEILambda().minLambda,
#      bu.ABEILambda().maxLambda,
#    ]

    for diagName, diagnosticConfig in diagnosticConfigs.items(): # only 1 diagnosticConfig
        if 'BinFunction' not in diagnosticConfig: continue

        logger.info('Calculating/writing diagnostic stats for:')
        logger.info('Diagnostic => '+diagName)
        Diagnostic = diagnosticConfig['BinFunction']
        outerIter = diagnosticConfig['outerIter']

        for varName in obsVars:
            logger.info('Variable => '+varName)

            varShort, varUnits = vu.varAttributes(varName)

            Diagnostic.evaluate(dbVals, varName, outerIter)
            diagValues = Diagnostic.result
            nLocs = len(diagValues)-np.isnan(diagValues).sum()

            if nLocs == 0:
                logger.warning('All missing values for diagnostic: '+diagName)

            for (binVarKey,binMethodKey), binMethod in binMethods.items(): #only 1 binMethod
                if binMethod.excludeDiag(diagName): continue

                binVarName, binGrpName = vu.splitObsVarGrp(binVarKey)
                binVarShort, binVarUnits = vu.varAttributes(binVarName)

                # initialize binMethod filter function result
                # NOTE: binning is performed using mean values
                #       and not ensemble member values
                binMethod.evaluate(dbVals, varName, outerIter)

                for binVal in binMethod.getvalues(): #only 1 binVal
                    # apply binMethod filters for binVal
                    binnedDiagnostic = binMethod.apply(diagValues,diagName,binVal)

                    # Plot horizontal distribution of binnedDiagnostic
                    if self.plotLambda and False:
                      logger.info('plotting obs-space diagnostic: '+diagName)
                      dotsize = 9.0
                      color = diagColors[diagName]
                      minmax = minmaxValue.get(diagName, [None, None])
                      nLocs = len(binnedDiagnostic)-np.isnan(binnedDiagnostic).sum()
                      plotDistri(
                        obsLatsDeg, obsLonsDeg, binnedDiagnostic,
                        ObsSpaceName, varShort, varUnits,
                        ObsSpaceName+'_'+binVarName+'='+binVal,
                        nLocs, diagName,
                        minmax[0], minmax[1], dotsize, color)

                    #TODO: create a base class to be used by both this and DiagnoseObsStatistics
                    #TODO: encapsulate below behavior in a member function
                    if diagName == 'ABEILambda':
                      lambdaVarName = 'modelLambda_'+varShort
                      resetLambda = localConfig.get('resetLambda', True)
                      inflationOutFile = varShort+'_'+self.args.inflationOutFile

                      if not os.path.exists(inflationOutFile) or resetLambda:
                        # create new NC file with same header info as grid file
                        logger.info('creating fresh inflation file')
                        modelUtils.headerOnlyFileFromTemplate(
                          modelGrid,
                          inflationOutFile,
                          date = self.args.datetime,
                        )

                      # create the inflation NC DataSet
                      inflationOut = modelUtils.getNCData(inflationOutFile, 'a')

                      if modelUtils.hasVar(lambdaVarName, inflationOut):
                        # read previous values for lambdaVarName (i.e., from previous sensor)
                        modelLambda = modelUtils.varRead(lambdaVarName, inflationOut)
                      else:
                        modelLambda = np.full_like(modelTemplateVars['1D-c']['values'],
                          bu.ABEILambda().minLambda)

                      # work in the range lambda >= 0.0
                      modelLambda -= bu.ABEILambda().minLambda

                      # alternative method using weight sum, not valid where NLocalObs is small
                      #modelWeightedSumNum = np.full_like(modelTemplateVars['1D-c']['values'],
                        #0.0)
                      #modelWeightedSumDenom = np.full_like(modelTemplateVars['1D-c']['values'],
                        #0.0)
                      #modelNLocalObs = np.full_like(modelTemplateVars['1D-c']['values'], 0, dtype=int)

                      ## Project observation-space inflation factors to
                      #  model-space using localization function
                      #  in sequential loop over the observation locations
                      #  Minamide and Zhang (2018), eq. 10
                      logger.info('projecting inflation factors to model space')

                      valid = np.isfinite(binnedDiagnostic)
                      s=0
                      nValid = valid.sum()
                      while s < nValid:
                        e=np.min([s+self.batchSize, nValid])

                        logger.info('  observation batch index limits => ['+str(s)+','+str(e)+')')
                        obsInds_ = obsInds[valid][s:e]
                        obsLats_ = obsLats[valid][s:e]
                        obsLons_ = obsLons[valid][s:e]

                        # work in the range lambda >= 0.0
                        obsLambdas_ = binnedDiagnostic[valid][s:e] - bu.ABEILambda().minLambda

                        #logger.info('initializing local distances')
                        modelLocal2obs.initNeighbors(obsLons_, obsLats_, self.nprocs)
                        normalizedDistance = modelLocal2obs.nnInterpDistances / self.localizationRadius

                        #See Greybush et al. (2011) for discussion
                        # https://doi.org/10.1175/2010MWR3328.1
                        # Greybush attributes Gaspari and Cohn (1999) with prescribing
                        # the Gaussian localization function in model-space and
                        # relates it to Hunt et al. (2007) application to obs-space
                        # ... could revisit rho function for more optimal results

                        # cut off when normalized distance is less than self.cutoffRadius
                        isLocal = bu.lessBound(normalizedDistance, self.cutoffRadius, False)

                        # localization, rho
                        rho = np.exp(-(normalizedDistance ** 2) / 2.0)

                        for ii, (obsInd, obsLambda) in enumerate(list(zip(
                             obsInds_, obsLambdas_))):

                          localModelInds = modelLocal2obs.nnInterpInds[ii, isLocal[ii,:]]

                          # interpolate current lambda field to this obs location
                          interpLambda = model2obs.applyAtIndex(modelLambda, obsInd)

                          # update lambda field
                          modelLambda[localModelInds] = \
                            modelLambda[localModelInds] + \
                            rho[ii, isLocal[ii,:]] * (obsLambda - interpLambda)

                          #modelWeightedSumNum[localModelInds] = \
                            #modelWeightedSumNum[localModelInds] + \
                            #rho[ii, isLocal[ii,:]] * (obsLambda - bu.ABEILambda().minLambda)

                          #modelWeightedSumDenom[localModelInds] = \
                            #modelWeightedSumDenom[localModelInds] + \
                            #rho[ii, isLocal[ii,:]]

                          #modelNLocalObs[localModelInds] = \
                            #modelNLocalObs[localModelInds] + 1

                        s=e

                      #maxLocalObs = float(modelNLocalObs.max())
                      #impactedCells = bu.greatBound(modelWeightedSumDenom, 0.0, False)

                      #modelLambda[impactedCells] = \
                        #modelLambda[impactedCells] + \
                        #modelWeightedSumNum[impactedCells] / modelWeightedSumDenom[impactedCells] * \
                        #modelNLocalObs[impactedCells].astype(float) / maxLocalObs

                      # truncate negative values caused by one of
                      # + non-conservative interpolation
                      # + obs-space discretization
                      # + numerical precision loss due to addition/subtraction
                      modelLambda[modelLambda < 0.0] = 0.0

                      # convert back to lambda >= 1.0
                      modelLambda += bu.ABEILambda().minLambda

                      if self.plotLambda:
                        # Plot horizontal distribution of modelLambda
                        logger.info('plotting model-space inflation factors')
                        dotsize = 9.0
                        color = diagColors[diagName]
                        minmax = minmaxValue.get(diagName, [None, None])
                        plotDistri(
                          modelLatsDeg, modelLonsDeg, modelLambda,
                          'model-'+ObsSpaceName.upper()+' '+diagName, varShort, '',
                          'model-'+ObsSpaceName+'_'+binVarName+'='+binVal,
                          0, diagName,
                          minmax[0], minmax[1], dotsize, color)

                      logger.info('writing inflation factors to NC file')

                      # write universal inflation factor used for all state variables below
                      attrs = {
                        'units': 'unitless',
                        'long_name': 'Column-wise Inflation factor derived from '+varShort,

                      }
                      modelUtils.varWrite(lambdaVarName, modelLambda, inflationOut,
                        attrs, modelTemplateVars['1D-c']['dims'], 'double')

                      # write state-variable-specific inflation factors
                      for modelInflateVarName in modelUtils.inflationVariables:
                        templateInfo = modelUtils.variableTraits[modelInflateVarName]
                        if modelUtils.hasVar(modelInflateVarName, modelGrid):
                          templateVarName = modelInflateVarName
                          templateVar = {
                            'attrs': modelUtils.varAttrs(modelInflateVarName, modelGrid),
                          }
                        else:
                          templateVarName = modelUtils.templateVariables[templateInfo['templateVar']]
                          templateVar = {
                            'attrs': {
                              'units': templateInfo['units'],
                              'long_name': templateInfo['long_name'],
                            },
                          }
                        templateVar['values'] = modelUtils.varRead(templateVarName, modelGrid)
                        templateVar['dims'] = deepcopy(modelUtils.varDims(templateVarName, modelGrid))
                        templateVar['datatype'] = modelUtils.varDatatype(templateVarName, modelGrid)

                        shape = templateVar['values'].shape
                        nDims = len(shape)
                        if nDims == 1:
                          modelInflateVar = modelLambda
                        elif nDims == 2:
                          modelInflateVar = np.empty_like(templateVar['values'])
                          for level in np.arange(0, shape[1]):
                            modelInflateVar[:,level] = modelLambda

                        modelUtils.varWrite(modelInflateVarName, modelInflateVar, inflationOut,
                          templateVar['attrs'], templateVar['dims'], templateVar['datatype'],
                        )

                      # close the inflation NC DataSet
                      inflationOut.close()

                #END binMethod.values LOOP
            #END binMethods tuple LOOP
        #END obsVars LOOP
    #END diagnosticConfigs LOOP

    modelGrid.close()


#=========================================================================
# main program
def main():
  _logger.info('Starting '+__name__)

  abei = GenerateABEIFactors()

  # diagnose abei
  abei.diagnose()

  _logger.info('Finished '+__name__+' successfully')

if __name__ == '__main__': main()
