#!/usr/bin/env python3

import GenerateABEIFactorsArgs

from basic_plot_functions import plotDistri
import binning_utils as bu
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

    # localization length scale for inflation projection
    self.localizationLength = 2.0e5 #meters ... same as all-sky gaussian thinning distance
    #self.localizationLength = 7.2e5 #meters ... 6x grid spacing
    #self.localizationLength = 2.e6 #meters same as background ensemble localization length

    # construct mean DB into 0th member slot
    self.logger.info('database path: '+self.args.dbPath)
    self.db = JediDB(self.args.dbPath, 'nc4')
    self.osNames = self.args.IRInstruments.split(',')
    dbOSNames = list(self.db.ObsSpaceName.values())
    self.osKeys = []
    for inst in self.osNames:
      assert inst in dbOSNames, 'GenerateABEIFactors::__init__ '+inst+' not included in JediDB'
      for osKey, osName in self.db.ObsSpaceName.items():
        if inst == osName:
          self.osKeys.append(osKey)
          break

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
#        markerSuffix = vu.depbgGroup
#    elif self.args.jediAppName == 'hofx':
#        markerSuffix = vu.hofxGroup
#    else:
#        logger.error('JEDI Application is not supported:: '+self.args.jediAppName)
#
#    return db.varList(osKey, obsFKey, markerSuffix)
    return [
      'brightness_temperature_8',
      'brightness_temperature_9',
      'brightness_temperature_10',
    ]
  @staticmethod
  def getBinVarConfigs():
    binVarConfigs = defaultdict(list)
    binVarConfigs[vu.obsVarQC] += [bu.goodQCMethod]
    return binVarConfigs

  @staticmethod
  def getBinMethods(binVarConfigs, ObsSpaceName):
    binMethods = {}
    for binVarKey, binMethodKeys in binVarConfigs.items():
        binVarConfig = pconf.binVarConfigs.get(binVarKey,pconf.nullBinVarConfig)
        for binMethodKey in binMethodKeys:
            config = binVarConfig.get(binMethodKey,pconf.nullBinMethod).copy()

            if (len(config['values']) < 1 or
                len(config['filters']) < 1): continue

            config['osName'] = ObsSpaceName
            binMethods[(binVarKey,binMethodKey)] = bu.BinMethod(config)

    return binMethods

  @staticmethod
  def getDiagNames():
    return ['ABEILambda']

  @staticmethod
  def getDBVars(obsVars, binMethods, diagnosticConfigs):
    dbVars = []
    for varName in obsVars:
        for diagName, diagnosticConfig in diagnosticConfigs.items():
            if 'ObsFunction' not in diagnosticConfig: continue

            # variables for diagnostics
            for varGrp in diagnosticConfig['ObsFunction'].dbVars(
                varName, diagnosticConfig['outerIter']):
                if diagnosticConfig[vu.mean]:
                    dbVars.append(varGrp)

            # variables for binning
            for (binVarKey,binMethodKey), binMethod in binMethods.items():
                for varGrp in binMethod.dbVars(
                    varName, diagnosticConfig['outerIter']):
                    dbVars.append(varGrp)

    dbVars.append(vu.latMeta)
    dbVars.append(vu.lonMeta)

    return dbVars


  def diagnoseObsSpace(self, db, osKey, localConfig):
    #  osKey - key of db dict to process
    logger = logging.getLogger(self.name+'.diagnoseObsSpace('+osKey+')')
    osName = self.osNames[self.osKeys.index(osKey)]

    # initialize mean db file handles
    db.initHandles(osKey)

    ###############################################
    ## Extract constructor info about the ObsSpace
    ###############################################

    ObsSpaceName  = db.ObsSpaceName[osKey]
    ObsSpaceInfo  = conf.DiagSpaceConfig[ObsSpaceName]
    ObsSpaceGrp   = ObsSpaceInfo['DiagSpaceGrp']

    obsVars = self.getObsVars(db, osKey)

    ########################################################
    ## Construct dictionary of binMethods for this ObsSpace
    ########################################################

    binVarConfigs = self.getBinVarConfigs()
    binMethods = self.getBinMethods(binVarConfigs, ObsSpaceName)

    ######################################
    ## Construct diagnostic configurations
    ######################################

    diagnosticConfigs = du.diagnosticConfigs(
        self.getDiagNames(), ObsSpaceName)
# TODO: change above line after bugfix/hofxDiagsCleanup is merged into develop
#        self.getDiagNames(), ObsSpaceName, False)

    ######################################
    ## Generate list of required variables
    ######################################

    dbVars = self.getDBVars(obsVars, binMethods, diagnosticConfigs)


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
    (modelLatsDeg, modelLonsDeg, EarthSphereR) = modelUtils.readGrid(gridFile=self.args.modelGridFile, returnR=True)
    modelLats = np.multiply(modelLatsDeg, np.pi / 180.0)
    modelLons = np.multiply(modelLonsDeg, np.pi / 180.0)

    modelTemplateVars = {}
    for templateType, templateVar in modelUtils.templateVariables.items():
      modelTemplateVars[templateType] = {
        'values': modelUtils.varRead(templateVar, self.args.modelGridFile),
        'dims': modelUtils.varDims(templateVar, self.args.modelGridFile),
      }

    # setup observation lat/lon
    obsLatsDeg = dbVals[vu.latMeta]
    obsLonsDeg = dbVals[vu.lonMeta]
    obsLats = np.multiply(obsLatsDeg, np.pi / 180.0)
    obsLons = np.multiply(obsLonsDeg, np.pi / 180.0)
    obsnLocs = len(obsLons)

    # setup interpolation object
    model2obs = Interpolate.InterpolateLonLat(
      modelLons, modelLats,
      weightMethod = 'barycentric',
      Radius = EarthSphereR)
    model2obs.initWeights(obsLons, obsLats)

    #TODO: move these attributes to diag_utils
    diagColors = {}
    minmaxValue = {}
    diagColors['ABEILambda'] = 'BuPu'
    minmaxValue['ABEILambda'] = [
      bu.ABEILambda().minLambda,
      bu.ABEILambda().maxLambda,
    ]

    for diagName, diagnosticConfig in diagnosticConfigs.items():
        if 'ObsFunction' not in diagnosticConfig: continue

        logger.info('Calculating/writing diagnostic stats for:')
        logger.info('Diagnostic => '+diagName)
        Diagnostic = diagnosticConfig['ObsFunction']
        outerIter = diagnosticConfig['outerIter']

        for varName in obsVars:
            logger.info('Variable => '+varName)

            varShort, varUnits = vu.varAttributes(varName)

            Diagnostic.evaluate(dbVals, varName, outerIter)
            diagValues = Diagnostic.result
            nLocs = len(diagValues)-np.isnan(diagValues).sum()

            if nLocs == 0:
                logger.warning('All missing values for diagnostic: '+diagName)

            for (binVarKey,binMethodKey), binMethod in binMethods.items():
                if diagName in binMethod.excludeDiags: continue

                binVarName, binGrpName = vu.splitObsVarGrp(binVarKey)
                binVarShort, binVarUnits = vu.varAttributes(binVarName)

                # initialize binMethod filter function result
                # NOTE: binning is performed using mean values
                #       and not ensemble member values
                binMethod.evaluate(dbVals, varName, outerIter)

                for binVal in binMethod.values:
                    # apply binMethod filters for binVal
                    binnedDiagnostic = binMethod.apply(diagValues,diagName,binVal)

                    # Plot horizontal distribution of binnedDiagnostic
                    if self.args.plotLambda:
                      logger.info('plotting obs-space diagnostic')
                      dotsize = 9.0
                      color = diagColors[diagName]
                      minmax = minmaxValue.get(diagName, [None, None])
                      nLocs = len(binnedDiagnostic)-np.isnan(binnedDiagnostic).sum()
                      plotDistri(
                        obsLatsDeg, obsLonsDeg, binnedDiagnostic,
                        osName, varShort, varUnits, osName+'_'+binVarName+'='+binVal,
                        nLocs, diagName,
                        minmax[0], minmax[1], dotsize, color)

                    #TODO: create a base class to be used by both this and DiagnoseObsStatistics
                    #TODO: encapsulate below behavior in a member function
                    if diagName == 'ABEILambda':
                      lambdaVarName = 'modelLambda_'+varShort
                      resetLambda = localConfig.get('resetLambda', True)
                      inflationFile = varShort+'_'+self.args.inflationFile

                      if not os.path.exists(inflationFile) or resetLambda:
                        # create new NC file with same header info as grid file
                        logger.info('creating fresh inflation file')
                        modelUtils.createHeaderOnlyFile(
                          self.args.modelGridFile,
                          inflationFile,
                          date = self.args.datetime,
                        )

                      if modelUtils.hasVar(lambdaVarName, inflationFile):
                        # read previous values for lambdaVarName (i.e., from previous sensor)
                        modelLambda = modelUtils.varRead(lambdaVarName, inflationFile)
                      else:
                        modelLambda = np.full_like(modelTemplateVars['1D-c']['values'],
                          bu.ABEILambda().minLambda)

                      ## Project observation-space inflation factors to
                      #  model-space using localization function
                      #  in sequential loop over the observation locations
                      #  Minamide and Zhang (2018), eq. 10
                      logger.info('projecting inflation factors to model space')
                      for obsInd, (obsLambda, obsLat, obsLon) in enumerate(list(zip(
                           binnedDiagnostic,  obsLats, obsLons))):
                        if np.isfinite(obsLambda):
                          gcDistances = Interpolate.HaversineDistance(
                                          obsLon, obsLat,
                                          modelLons, modelLats,
                                          EarthSphereR)
                          scales = gcDistances / self.localizationLength
                          rho = np.zeros(scales.shape)

                          #See Greybush et al. (2011) for discussion
                          # https://doi.org/10.1175/2010MWR3328.1
                          # Greybush attributes Gaspari and Cohn (1999) with prescribing
                          # the Gaussian localization function in model-space and
                          # relates it to Hunt et al. (2007) application to obs-space
                          # ... could revisit rho function for more optimal results
                          crit = scales < 3.5
                          rho[crit] = np.exp(- (scales[crit] ** 2) / 2.0)
                          crit = (np.abs(rho) > 0.0)
                          modelLocInds = np.arange(0, len(rho))[crit]

                          # interpolate current lambda field to this obs location
                          interpLambda = model2obs.applyAtIndex(modelLambda, obsInd)

                          # update lambda field
                          modelLambda[modelLocInds] = modelLambda[modelLocInds] + \
                            rho[crit] * (obsLambda - interpLambda)

                      # truncate negative valuess caused by obs-space discretization
                      modelLambda[modelLambda < bu.ABEILambda().minLambda] = \
                        bu.ABEILambda().minLambda

                      if self.args.plotLambda:
                        # Plot horizontal distribution of modelLambda
                        logger.info('plotting model-space inflation factors')
                        dotsize = 9.0
                        color = diagColors[diagName]
                        minmax = minmaxValue.get(diagName, [None, None])
                        plotDistri(
                          modelLatsDeg, modelLonsDeg, modelLambda,
                          'model-'+osName.upper()+' '+diagName, varShort, '',
                          'model-'+osName+'_'+binVarName+'='+binVal,
                          0, diagName,
                          minmax[0], minmax[1], dotsize, color)

                      logger.info('writing inflation factors to NC file')

                      # write universal inflation factor used for all state variables below
                      attrs = {
                        'units': 'unitless',
                        'long_name': 'Column-wise Inflation factor derived from '+varShort,

                      }
                      modelUtils.varWrite(lambdaVarName, modelLambda, inflationFile,
                        attrs, modelTemplateVars['1D-c']['dims'], 'double')

                      # write state-variable-specific inflation factors
                      for modelInflateVarName in modelUtils.inflationVariables:
                        templateInfo = modelUtils.variableTraits[modelInflateVarName]
                        if modelUtils.hasVar(modelInflateVarName, self.args.modelGridFile):
                          templateVarName = modelInflateVarName
                          templateVar = {
                            'attrs': modelUtils.varAttrs(modelInflateVarName, self.args.modelGridFile),
                          }
                        else:
                          templateVarName = modelUtils.templateVariables[templateInfo['templateVar']]
                          templateVar = {
                            'attrs': {
                              'units': templateInfo['units'],
                              'long_name': templateInfo['long_name'],
                            },
                          }
                        templateVar['values'] = modelUtils.varRead(templateVarName, self.args.modelGridFile)
                        templateVar['dims'] = modelUtils.varDims(templateVarName, self.args.modelGridFile)
                        templateVar['datatype'] = modelUtils.varDatatype(templateVarName, self.args.modelGridFile)

                        shape = templateVar['values'].shape
                        nDims = len(shape)
                        if nDims == 1:
                          modelInflateVar = modelLambda
                        elif nDims == 2:
                          modelInflateVar = np.empty_like(templateVar['values'])
                          for level in np.arange(0, shape[1]):
                            modelInflateVar[:,level] = modelLambda

                        modelUtils.varWrite(modelInflateVarName, modelInflateVar, inflationFile,
                          templateVar['attrs'], templateVar['dims'], templateVar['datatype'],
                        )

                #END binMethod.values LOOP
            #END binMethods tuple LOOP
        #END obsVars LOOP
    #END diagnosticConfigs LOOP


#=========================================================================
# main program
def main():
  _logger.info('Starting '+__name__)

  abei = GenerateABEIFactors()

  # diagnose abei
  abei.diagnose()

  _logger.info('Finished '+__name__+' successfully')

if __name__ == '__main__': main()
