#!/usr/bin/env python3

import basic_plot_functions as bpf
import binning_utils as bu
import predefined_configs as pconf
from collections.abc import Iterable
from collections import defaultdict, OrderedDict
from config import DiagSpaceConfig
from copy import deepcopy
import collections
import datetime as dt
import diag_utils as du
from fit2D import fit2D, poly2DEquation
import inspect
import itertools as itt
import logging
import multiprocessing as mp
import numpy as np
from pathlib import Path
import plot_utils as pu
import re
import os
import stat_utils as su
import StatisticsDatabase as sdb
from textwrap import indent
import var_utils as vu
import yaml

## plot settings
figureFileType = 'pdf' #['pdf','png']

interiorLabels = True

dataYAMLMissingFloat = 999. # fill value when ~np.isfinite
dataYAMLPrecision = 4 # number of significant digits
def dataYAMLFmtFloat(f):
  if np.isfinite(f):
    return float(('{:.'+str(dataYAMLPrecision-1)+'e}').format(f))
  else:
    return dataYAMLMissingFloat

def anWorkingDir(DiagSpace, analysisType):
  return DiagSpace+'_analyses'+'/'+analysisType

###################################
## Base class for all analysisTypes
###################################
class AnalysisBase():
    def __init__(self, db, analysisType, diagnosticGroupings = {}):
        self.analysisType = analysisType
        self.DiagSpaceName = db.DiagSpaceName
        self.diagnosticGroupings = diagnosticGroupings

        self.logger = logging.getLogger(__name__+'.'+self.DiagSpaceName+'.'+self.analysisType)

        '''
        use relativeErrorType to choose between three relativeErrorFunction options:
          disable - use absolute diagnostics for all experiments

          one hundred centered - use oneHundredCenteredPercentDifference of
            positive-semidefinite diagnostics for non-control experiments

          zero centered (default) - use zeroCenteredPercentDifference of
            positive-semidefinite diagnostics for non-control experiments
        '''

        self.relativeErrorType = 'zero centered'

        assert self.relativeErrorType in ['zero centered', 'one hundred centered', 'disable'], (
          self.__class__.__name__+': invalid relativeErrorType: '+str(self.relativeErrorType))

        relativeErrorConfigurations = {
          'disable': {
            'center': None,
            'function': None,
            'labeler': None,
            'limiter': None,
          },
          'one hundred centered': {
            'center': 100.,
            'function': self.oneHundredCenteredPercentDifference,
            'labeler': self.oneHundredCenteredLabeler,
            'limiter': self.oneHundredCenteredLimiter,
          },
          'zero centered': {
            'center': 0.,
            'function': self.zeroCenteredPercentDifference,
            'labeler': self.zeroCenteredLabeler,
            'limiter': self.zeroCenteredLimiter,
          },
        }
        rconf = relativeErrorConfigurations[self.relativeErrorType]
        self.relativeErrorCenter = rconf['center']
        self.relativeErrorFunction = rconf['function']
        self.relativeErrorLabeler = rconf['labeler']
        self.relativeErrorLimiter = rconf['limiter']

        ## Extract useful variables from the database
        self.db = db
        self.diagnosticConfigs = db.diagnosticConfigs
        self.availableDiagnostics = list(self.diagnosticConfigs.keys())

        self.expNames = self.db.expNames
        self.nExp = len(self.expNames)
        self.cntrlExpName = self.db.cntrlExpName
        self.noncntrlExpNames = self.db.noncntrlExpNames

        self.fcTDeltas = self.db.fcTDeltas
        self.fcTDeltas_totmin = self.db.fcTDeltas_totmin
        self.fcMap = list(zip(self.fcTDeltas, self.fcTDeltas_totmin))
        self.nFC = len(self.fcTDeltas)

        self.cyDTimes = self.db.cyDTimes
        self.nCY = len(self.cyDTimes)

        varNames = np.asarray(self.db.varNames)
        chlist = np.asarray(self.db.chlist)
        varUnitss = np.asarray(self.db.varUnitss)

        # use DiagSpaceConfig settings in config.py to select plotting order for variables/channels
        # only those variables/channels that are in self.db and in DiagSpaceConfig will be
        # plotted
        # TODO: alternatively could do this in StatisticsDatabase.py to reduce memory usage
        # allow for two optionl settings for each DiagSpaceConfig element
        #  + 'diagnosed channels' or 'diagnosed variables' (when statistics are generated, TODO)
        #  + 'analyzed channels' or 'analyzed variables' (here, when plots are created)
        if any(chlist < 0):
          dbVarsList = list(varNames)
          orderConf = 'analyzed variables'
        else:
          dbVarsList = list(chlist)
          orderConf = 'analyzed channels'

        selectOrder = DiagSpaceConfig[self.DiagSpaceName].get(orderConf, dbVarsList)
        allVarsIndex = np.full_like(selectOrder, -1, dtype=int)

        for ii, value in enumerate(selectOrder):
          if value in dbVarsList:
            allVarsIndex[ii] = dbVarsList.index(value)

        allVarsIndex = allVarsIndex[allVarsIndex>=0]
        varNames = varNames[allVarsIndex]
        varUnitss = varUnitss[allVarsIndex]
        chlist = chlist[allVarsIndex]

        varLabels = []
        # kludge to override some labels
        # TODO: extend varDictObs to include labels similar to du.availableDiagnostic
        self.labelReplacements = {
            'Ref (%)': 'refractivity (N-units)',
            'Bnd (%)': 'bending angle (rad)',
            'ModLev': 'model level',
            'alt (m)': 'altitude (m)',
            'lat ('+vu.degree+')': 'latitude ('+vu.degree+')',
            'impact_height (m)': 'impact height (m)',
        }
        for (varName, varUnits) in zip(varNames, varUnitss):
            label = varName
            if varUnits != vu.miss_s:
                label = label+' ('+varUnits+')'
            for orig, sub in self.labelReplacements.items():
              label = label.replace(orig, sub)
            varLabels.append(label)

        self.varNames = list(varNames)
        self.chlist = list(chlist)
        self.varMap = list(zip(self.varNames, varLabels))
        self.nVars = len(self.varNames)

        self.allBinStrVals = self.db.allBinStrVals
        self.allBinNumVals = self.db.allBinNumVals
        self.allBinNumVals2DasStr = self.db.allBinNumVals2DasStr

        ## Establish default configuration
        self.blocking = False
        self.parallelism = False

        # TODO(JJG): decide if nproc is needed
        # nproc could be used to initialize workers within a single
        # analysisType to be distributed however they would be most useful.
        # That would allow for more granularity of computational effort
        # but requires that analysisType to be "blocking" due to the nature
        # of multiprocessing Pool's.
        # self.nproc = nproc

        self.requiredStatistics = []
        self.requestAggDFW = False
        self.blankBinMethodFile = bu.identityBinMethod

        self.subplotWidth = 2.5
        self.subplotAspect = 1.0

        self.MAX_FC_SUBFIGS = 6
        self.MAX_FC_LINES = 6

        ## Setup paths
        CWD = os.getcwd()
        workDir = anWorkingDir(self.DiagSpaceName, self.analysisType)
        nWDPieces = len(workDir.split('/'))
        localWD = '/'.join(CWD.split('/')[-nWDPieces:])
        self.WorkPath = Path('./')
        if localWD != workDir:
            self.WorkPath = self.WorkPath/workDir
        self.myFigPath = self.WorkPath
        self.myFigPath.mkdir(parents=True, exist_ok=True)

    def oneHundredCenteredPercentDifference(self,
      experiment,
      reference,
      dmin0 = np.NaN,
      dmax0 = np.NaN,
    ):
      exp = experiment.astype(float)
      ref = reference.astype(float)

      out = np.full_like(exp, np.NaN)

      validDenom = bu.greatBound(np.abs(ref), 0., False)

      out[validDenom] = 100. * (
        np.divide(
          np.subtract(
            exp[validDenom],
            ref[validDenom]
          ),
          ref[validDenom]
        ) + 1.0
      )

      dmin, dmax = self.oneHundredCenteredLimiter(dmin0, dmax0, out)

      centralValue = 100.0

      label = self.oneHundredCenteredLabeler()

      return out, dmin, dmax, centralValue, label

    def oneHundredCenteredLabeler(self):
      return '% of '+self.cntrlExpName
      #return '100 x [1+\n(EXP-'+self.cntrlExpName+') / \n'+self.cntrlExpName+']'
      #return '100 x [1+\n(EXP-CONTROL)/\nCONTROL]'

    @staticmethod
    def initLimits(dmin0=np.NaN, dmax0=np.NaN, d=None):
      dmin = dmin0
      dmax = dmax0
      if d is not None and not (np.isfinite(dmin) or np.isfinite(dmax)):
        dmin = np.nanmin(d)
        dmax = np.nanmax(d)
      return dmin, dmax

    def oneHundredCenteredLimiter(self, dmin0=np.NaN, dmax0=np.NaN, d=None):
      dmin, dmax = self.initLimits(dmin0, dmax0, d)

      # update dmax, conditional on dmin
      dmax = np.nanmax([dmax, 100.0/(dmin/100.0)])

      # set absolute min/max in case there are large outliers
      #if np.isfinite(dmin):
      dmin = np.nanmax([dmin, 66.7])
      dmin = np.nanmin([dmin, 98.0])

      #if np.isfinite(dmax):
      dmax = np.nanmin([dmax, 150.0])
      dmax = np.nanmax([dmax, 102.0])

      return dmin, dmax

    def zeroCenteredPercentDifference(self,
      experiment,
      reference,
      dmin0 = np.NaN,
      dmax0 = np.NaN,
    ):
      exp = experiment.astype(float)
      ref = reference.astype(float)

      out = np.full_like(exp, np.NaN)

      validDenom = bu.greatBound(np.abs(ref), 0., False)

      out[validDenom] = 100. * (
        np.divide(
          np.subtract(
            exp[validDenom],
            ref[validDenom]
          ),
        ref[validDenom]
        )
      )

      dmin, dmax = self.zeroCenteredLimiter(dmin0, dmax0, out)

      centralValue = 0.0

      label = self.zeroCenteredLabeler()

      return out, dmin, dmax, centralValue, label

    def zeroCenteredLabeler(self):
      return '% diff. from '+self.cntrlExpName
      #return '100 x (EXP-'+self.cntrlExpName+') / \n'+self.cntrlExpName+''
      #return '100 x (EXP-CONTROL)/\nCONTROL'

    def zeroCenteredLimiter(self, dmin0=np.NaN, dmax0=np.NaN, d=None):
      dmin, dmax = self.initLimits(dmin0, dmax0, d)

      # update dmax, conditional on dmin
      #dmax = np.nanmax([dmax, 100.0/((100.0+dmin)/100.0) - 100.0])
      dmax = np.nanmax([dmax, 1.0/dmin])

      # set absolute min/max in case there are large outliers
      #if np.isfinite(dmin):
      dmin = np.nanmax([dmin, -33.3])
      dmin = np.nanmin([dmin, -2.0])

      #if np.isfinite(dmax):
      dmax = np.nanmin([dmax, 50.0])
      dmax = np.nanmax([dmax, 2.0])

      return dmin, dmax


    def binMethodFile(self, binMethod, before = True):
        '''
        Format the binMethod for file naming
        '''
        binMethodFile = ''
        if binMethod != self.blankBinMethodFile:
            if before:
                binMethodFile = '_'+binMethod
            else:
                binMethodFile = binMethod+'_'

        return binMethodFile

    def fcName(self, diagnosticGroup):
        '''
        Format the diagnosticGroup for forecast analysisType's
        '''
        fcDiagName = diagnosticGroup
        if self.fcTDeltas[-1] > dt.timedelta(0):
            fcDiagName = fcDiagName.replace('omm','omf')
            fcDiagName = fcDiagName.replace('omb','omf')
            fcDiagName = fcDiagName.replace('bmo','fmo')
            fcDiagName = fcDiagName.replace('mmo','fmo')
            fcDiagName = fcDiagName.replace('hmo','fmo')
            fcDiagName = fcDiagName.replace('bak','fc')
        return fcDiagName

    def expDiagnosticLabel(self, expName, diagnosticName, allDiagnosticNames):
        if len(allDiagnosticNames) > 1:
            diagLabel = self.diagnosticConfigs[diagnosticName]['label']
            return r''+str(expName)+' – '+str(diagLabel)
        else:
            return str(expName)

    def statPlotAttributes(self, diagnosticGroup, statName,
                           allDiagnosticNames=None, isDifferencePlot=False, isPercentDiffPlot=False):
        '''
        Define plotting attributes for the combination of diagnosticGroup and statName
        '''
        ommDiagnostics = ['omb', 'oma', 'omm', 'omf', 'omb_nobc', 'oma_nobc', 'omm_nobc',
                          'rltv_omb', 'rltv_oma', 'rltv_omm', 'rltv_omf']
        mmoDiagnostics = ['bmo', 'amo', 'mmo', 'fmo', 'mmgfsan', 'rltv_mmgfsan', 'log_mogfsan']
        truncateDiagnostics = ommDiagnostics+mmoDiagnostics
        diagnosticGroup_ = diagnosticGroup
        #for diag in truncateDiagnostics:
        #    if pu.prepends(diag, diagnosticGroup_) or pu.postpends(diag, diagnosticGroup_):
        #        diagnosticGroup_ = diag

        allDiagnosticNames_ = deepcopy(allDiagnosticNames)
        if allDiagnosticNames_ is None:
            cntrlDiagnosticName = diagnosticGroup_
            allDiagnosticNames_ = [diagnosticGroup_]
        else:
            cntrlDiagnosticName = allDiagnosticNames_[0]

        # retrieve label
        diagLabel = self.diagnosticConfigs.get(diagnosticGroup_, {}).get('label', diagnosticGroup_)

        fcDiagName = self.fcName(diagnosticGroup_)
        fcDiagLabel = self.diagnosticConfigs.get(fcDiagName, {}).get('label', diagnosticGroup_)

        if statName in du.diagnosticIndependentStatistics:
            statDiagLabel = statName
            fcstatDiagLabel = statName
        elif statName == vu.miss_s:
            statDiagLabel = diagLabel
            fcstatDiagLabel = fcDiagLabel
        elif set(allDiagnosticNames_).issubset(set(du.statisticDependentDiagnostics)):
            statDiagLabel = diagLabel+'('+statName+')'
            fcstatDiagLabel = fcDiagLabel+'('+statName+')'
        else:
            statDiagLabel = statName+'('+diagLabel+')'
            fcstatDiagLabel = statName+'('+fcDiagLabel+')'

        #The following attributes apply to unbounded and/or symmetric plotted quantities
        # e.g., omb, oma, ana/bak for velocity, differences

        # default axis formatting
        sciTicks = False
        logScale = False
        centralValue = None

        #for diag in truncateDiagnostics:
        #    if pu.prepends(diag, cntrlDiagnosticName) or pu.postpends(diag, cntrlDiagnosticName):
        #        cntrlDiagnosticName = diag
        #    for idiag, adiag in enumerate(allDiagnosticNames_):
        #        if pu.prepends(diag, adiag) or pu.postpends(diag, adiag):
        #            allDiagnosticNames_[idiag] = diag

        oneCenteredRatioDiagPrefixes = ['CRx', 'CRy', 'InnovationRatio', 'SRx', 'SRy', 'OENI']
        allOneCenteredDiagnostics = statName not in ['Count', 'Mean', 'Skew', 'ExcessKurtosis']
        for diag in allDiagnosticNames_:
            oneCentered = False
            for diagPrefix in oneCenteredRatioDiagPrefixes:
                oneCentered = oneCentered or pu.prepends(diagPrefix, diag)
            allOneCenteredDiagnostics = allOneCenteredDiagnostics and oneCentered
        if allOneCenteredDiagnostics:
            centralValue = 1.
            logScale = True

        # TODO: add log-scaling for some diagnostics, need to ensure log range is reasonable,
        #       including ignoring outliers that over-smooth the rest of the colors
        # logScaledDiagPrefixes = ['sigmao', 'ideal-sigmao', 'sigma']
        # allLogScaledDiagnostics = statName not in ['Count', 'Mean']
        # for diag in allDiagnosticNames_:
        #     logScaled = False
        #     for diagPrefix in logScaledDiagPrefixes:
        #         logScaled = logScaled or pu.prepends(diagPrefix, diag)
        #     allLogScaledDiagnostics = allLogScaledDiagnostics and logScaled
        # if allLogScaledDiagnostics:
        #     logScale = True and not isDifferencePlot

        if statName=='Count':
            sciTicks = True
            logScale = True

        cntrlExpDiagnosticLabel = self.expDiagnosticLabel(self.cntrlExpName, cntrlDiagnosticName, allDiagnosticNames_)
        statDiagDiffLabel = statDiagLabel+': [EXP-'+cntrlExpDiagnosticLabel+']'
        fcstatDiagDiffLabel = fcstatDiagLabel+': [EXP-'+cntrlExpDiagnosticLabel+']'
        #statDiagPercentDiffLabel = '100 x [EXP-'+cntrlExpDiagnosticLabel+'] / \n'+cntrlExpDiagnosticLabel
        statDiagPercentDiffLabel = '% diff. from '+cntrlExpDiagnosticLabel

        if statName in ['Mean', 'Skew', 'ExcessKurtosis']:
            if diagnosticGroup_ in ommDiagnostics:
                centralValue = 0.
                statDiagDiffLabel = statDiagLabel+': ['+cntrlExpDiagnosticLabel+'-EXP]'
                fcstatDiagDiffLabel = fcstatDiagLabel+': ['+cntrlExpDiagnosticLabel+'-EXP]'
                #statDiagPercentDiffLabel = '100 x ['+cntrlExpDiagnosticLabel+'-EXP] / \n'+cntrlExpDiagnosticLabel
                statDiagPercentDiffLabel = '% diff. to '+cntrlExpDiagnosticLabel
            if diagnosticGroup_ in mmoDiagnostics:
                centralValue = 0.

            if 'OENI' in diagnosticGroup: centralValue = 0.
            if 'ACI' in diagnosticGroup: centralValue = 0.

        if 'rltv_mmgfsan' in diagnosticGroup_ and statName in ['STD', 'RMS']: logScale = True
        if 'log_mogfsan' in diagnosticGroup_ and statName in ['Mean']: centralValue = 0.

        #if 'rltv_mmgfsan' in diagnosticGroup_: logScale = True

        if isDifferencePlot:
          statDiagLabel = statDiagDiffLabel
          fcstatDiagLabel = fcstatDiagDiffLabel

        if isPercentDiffPlot:
          statDiagLabel = statDiagPercentDiffLabel
          fcstatDiagLabel = statDiagPercentDiffLabel


        pubConventions = {
          'Mean': 'mean',
          'RMS': 'rms',
          'rms(mmgfsan)': 'RMSE',
          'rms(omf)': 'RMSE',
          'rms(omm)': 'RMSE',
          'rms(dx)': 'rms',
          'rms(dy)': 'rms',
          'ddT': '$dd^T$',
        }
        for orig, sub in pubConventions.items():
          statDiagLabel = statDiagLabel.replace(orig, sub)
          fcstatDiagLabel = fcstatDiagLabel.replace(orig, sub)

        # change to raw string to allow for TeX parsing
        statDiagLabel = r''+statDiagLabel
        fcstatDiagLabel = r''+fcstatDiagLabel

        return statDiagLabel, fcstatDiagLabel, sciTicks, logScale, centralValue

    def UNIONcntrlANDexpCYDTimes(self, dfw, myLoc = {}):
        '''
        Determine the union of cyDTimes available between
        the control and each other experiment at each fcTDelta
        '''
        cntrlLoc = deepcopy(myLoc)
        cntrlLoc['expName'] = self.cntrlExpName

        expsCYDTimes = {}
        for fcTDelta in self.fcTDeltas:
            cntrlLoc['fcTDelta'] = fcTDelta
            cntrlCYDTimes = set(dfw.levels('cyDTime', cntrlLoc))

            expLoc = deepcopy(cntrlLoc)
            for expName in self.expNames:
                expLoc['expName'] = expName
                expCYDTimes = set(dfw.levels('cyDTime', expLoc))
                expsCYDTimes[(expName, fcTDelta)] = list(cntrlCYDTimes & expCYDTimes)
                if len(cntrlCYDTimes) != len(expCYDTimes):
                    self.logger.warning(
                      self.cntrlExpName+' and '+expName+
                      ' have different number of CYDTimes at forecast length '+
                      str(fcTDelta.total_seconds()/3600.)+
                      'hr. Only using common CYDTimes for CI calculation.')
        return expsCYDTimes

    def analyze(self, workers = None):
        self.logger.info('analyze()')
        if self.blocking or workers is None:
            # analyses with internal blocking
            self.analyze_()
        elif self.parallelism:
            # each analysis determines how to split external workers
            self.analyze_(workers)
        else:
            # divide workers acros analyses without internal parallelism/blocking
            workers.apply_async(self.analyze_)

    def analyze_(self, workers = None):
        '''
        virtual method
        '''
        raise NotImplementedError()




def categoryBinValsAttributes(dfw, fullBinVar, binMethod, options):
    '''
    Utility function for providing an ordered list of
    pairs of binVals and associated labels for
    category binMethods in the context of a DFWrapper
    '''

    binVar = vu.varDictAll.get(fullBinVar, [None, fullBinVar])[1]

    dbSelect1DBinVals = dfw.levels('binVal')
    binUnitss = dfw.uniquevals('binUnits')
    #if (len(binUnitss) == 0 or
    #    len(dbSelect1DBinVals) == 1): return None, None
    assert (len(binUnitss) != 0 and len(dbSelect1DBinVals) > 0), 'ERROR: categoryBinValsAttributes received invalid binVar/binMethod'

    binUnits = binUnitss[0]

    # reorder select1DBinVals to match binMethod definition
    # TODO(JJG): clean up for readability
    select1DBinVals = []
    binMethods = pconf.binVarConfigs.get(fullBinVar, {})
    for method, value in binMethods.items():
      if binMethod is None or method == binMethod:
        tmp = value.get('values', dbSelect1DBinVals)
        if (not isinstance(tmp, Iterable) or
          isinstance(tmp, str)):
          select1DBinVals += [tmp]
        else:
          select1DBinVals += list(tmp)

    for b in dbSelect1DBinVals:
        select1DBinVals.append(b)
    for b in list(select1DBinVals):
        if b not in dbSelect1DBinVals:
            select1DBinVals.remove(b)

    select1DBinVals = pu.uniqueMembers(select1DBinVals)

    binTitles = []
    for binVal in select1DBinVals:
        t = ''
        if pu.isfloat(binVal) or pu.isint(binVal):
            t += ' @ '+binVar+'='+binVal
            if binUnits != vu.miss_s:
                t += ' '+binUnits
        elif binVal not in pconf.goodFlagNames+[bu.allQCMethod]:
            t += ' @ '+binVal
        binTitles.append(t)

    binValsMap = list(zip(select1DBinVals, binTitles))

    return binValsMap


#=============
# 1-D figures
#=============
class CategoryBinMethodBase(AnalysisBase):
    '''
    Base class used to analyze statistics across binMethods with zero-dimensioned or
      category binValues, e.g., QC flag, named latitude band, cloudiness regime, surface type
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)
        self.parallelism = True
        self.maxBinVarTier = 2

        # default binVar/binMethod combinations
        self.binVarDict = {
            (vu.obsVarQC, bu.goodQCMethod): {},
            (vu.obsVarLat, bu.latbandsMethod): {},
            (vu.obsVarLat, bu.troplatbandsMethod): {},
            (vu.obsVarLat, bu.polarlatbandsMethod): {},
            (vu.obsVarCldFracX, bu.cloudbandsMethod): {},
            (vu.obsVarCldFracY, bu.cloudbandsMethod): {},
            # vu.modVarLat is redundant with vu.obsVarLat (both have varShort=="lat")
            #(vu.modVarLat, bu.latbandsMethod): {},
            (vu.noBinVar, bu.noBinMethod): {},
            (vu.obsRegionBinVar, bu.geoirlatlonboxMethod): {'binVarTier': 2},
            (vu.modelRegionBinVar, bu.geoirlatlonboxMethod): {'binVarTier': 2},
            (vu.obsVarPrs, bu.PjetMethod): {'binVarTier': 3},
            (vu.obsVarAlt, bu.altjetMethod): {'binVarTier': 3},
            (vu.obsVarImpact, bu.altjetMethod): {'binVarTier': 3},
            (vu.obsVarLandFrac, bu.surfbandsMethod): {'binVarTier': 3},
        }
        self.maxDiagnosticsPerAnalysis = 10 // self.nExp

    def subplotArrangement(self, binValsMap):
        # subplot configuration
        if len(binValsMap) > 1:
            nxplots = len(binValsMap)
            nyplots = self.nVars
            nsubplots = nxplots * nyplots
        else:
            nsubplots = self.nVars
            nxplots = np.int(np.ceil(np.sqrt(nsubplots)))
            while nsubplots%nxplots > 0 and nsubplots%nxplots / nxplots <= 0.5: nxplots += 1
            nyplots = np.int(np.ceil(np.true_divide(nsubplots, nxplots)))

        return nxplots, nyplots, nsubplots

    def analyze_(self, workers = None):
        useWorkers = (not self.blocking and self.parallelism and workers is not None)

        # TODO(JJG): construct member Diagnostic objects (create new class) from
        #            diagnosticConfigs instead of referencing dictionary
        #            entries below.
        # TODO(JJG): use same color, vary line style within diagnosticGroupings
        diagnosticGrouped = {}
        for diag in self.availableDiagnostics:
            diagnosticGrouped[diag] = False

        diagnosticGroupings = deepcopy(self.diagnosticGroupings)
        for group in list(diagnosticGroupings.keys()):
            diags = diagnosticGroupings[group]
            if (len(diags) > self.maxDiagnosticsPerAnalysis or
               not set(diags).issubset(set(list(self.availableDiagnostics)))):
                del diagnosticGroupings[group]
                continue
            for diag in diags: diagnosticGrouped[diag] = True

        for diag in self.availableDiagnostics:
            if not diagnosticGrouped[diag]:
                diagnosticGroupings[diag] = [diag]

        for diagnosticGroup, diagnosticNames in diagnosticGroupings.items():
            if len(diagnosticNames) > self.maxDiagnosticsPerAnalysis: continue
            if len(set(diagnosticNames) & set(self.availableDiagnostics)) == 0: continue
            diagnosticConfigs = {}
            selectedStatistics = set([])
            for diagnosticName in diagnosticNames:
                diagnosticConfigs[diagnosticName] = deepcopy(self.diagnosticConfigs[diagnosticName])
                selectedStatistics = set(list(selectedStatistics) +
                                         diagnosticConfigs[diagnosticName]['selectedStatistics'])
            availableStatistics = set([])
            for diagnosticName in diagnosticNames:
                diagnosticConfigs[diagnosticName] = deepcopy(self.diagnosticConfigs[diagnosticName])
                availableStatistics = set(list(availableStatistics) +
                                         diagnosticConfigs[diagnosticName]['availableStatistics'])
            if not set(self.requiredStatistics).issubset(availableStatistics): continue

            diagLoc = {'diagName': diagnosticNames}
            diagBinVars = self.db.dfw.levels('binVar', diagLoc)
            diagBinMethods = self.db.dfw.levels('binMethod', diagLoc)
            for (fullBinVar, binMethod), options in self.binVarDict.items():
                if options.get('binVarTier', 1) > self.maxBinVarTier: continue
                binVar = vu.varDictAll.get(fullBinVar, [None, fullBinVar])[1]
                if (binVar not in diagBinVars or
                    (binMethod is not None and binMethod not in diagBinMethods)): continue

                self.logger.info(diagnosticGroup+', '+binVar+', '+str(binMethod))

                if useWorkers:
                    workers.apply_async(self.innerloopsWrapper,
                        args = (diagnosticGroup, diagnosticConfigs, fullBinVar, binMethod, selectedStatistics, options))
                else:
                    self.innerloopsWrapper(
                        diagnosticGroup, diagnosticConfigs, fullBinVar, binMethod, selectedStatistics, options)

    def innerloopsWrapper(self,
        diagnosticGroup, diagnosticConfigs, fullBinVar, binMethod, selectedStatistics, options):

        binVar = vu.varDictAll.get(fullBinVar, [None, fullBinVar])[1]

        # narrow mydfwDict by binVar and binMethod to reduce run-time and memory
        myLoc = {}
        myLoc['binVar'] = binVar
        myLoc['binMethod'] = binMethod

        mydfwDict = {'dfw': self.db.loc(myLoc)}

        # aggregate statistics when requested
        if self.requestAggDFW:
            mydfwDict['agg'] = sdb.DFWrapper.fromAggStats(mydfwDict['dfw'], ['cyDTime'])
            sdb.createORreplaceDerivedDiagnostics(mydfwDict['agg'], diagnosticConfigs)

        # further narrow mydfwDict by diagName
        # NOTE: derived diagnostics may require multiple diagName values;
        # can only narrow by diagName after aggregation
        myLoc['diagName'] = list(diagnosticConfigs.keys())
        for key in mydfwDict.keys():
            mydfwDict[key] = sdb.DFWrapper.fromLoc(mydfwDict[key], myLoc)

        binValsMap = categoryBinValsAttributes(
            mydfwDict['dfw'], fullBinVar, binMethod, options)

        nxplots, nyplots, nsubplots = self.subplotArrangement(binValsMap)

        for statName in selectedStatistics:
            if statName not in options.get('onlyStatNames', selectedStatistics): continue

            self.innerloops(
                mydfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
                nsubplots, nxplots, nyplots)

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):
        '''
        virtual method
        '''
        raise NotImplementedError()


class CYAxisExpLines(CategoryBinMethodBase):
    '''
    Creates a timeseries figure between firstCycleDTime and lastCycleDTime
      for each forecast length between fcTDeltaFirst and fcTDeltaLast
      -  x-axis: cycle initial time
      -    line: per experiment
      - subplot: combination of DiagSpace variable and binVal
      -    file: combination of binVar, statistic, and FC lead time (if applicable)
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.subplotWidth = 2.4
        self.subplotAspect = 0.75

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nCY < 2: return
        assert myLoc['binMethod'] is not None, 'CYAxisExpLines.innerloops: binMethod cannot be None'

        bgstatDiagLabel, fcstatDiagLabel, sciTicks, logScale, centralValue = \
            self.statPlotAttributes(diagnosticGroup, statName, myLoc['diagName'])

        figPath = self.myFigPath/diagnosticGroup
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        lineLoc = {}
        axisLimitsLoc = {}

        #file loop 1
        for (fcTDelta, fcTDelta_totmin) in self.fcMap:
            # calculate valid time for x-axis
            xVals = []
            for cyDTime in self.cyDTimes:
                xVals.append(cyDTime+fcTDelta)

            lineLoc['fcTDelta'] = fcTDelta
            axisLimitsLoc['fcTDelta'] = fcTDelta

            # establish a new figure
            fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, interiorLabels)
            iplot = 0

            #subplot loop 1
            for (varName, varLabel) in self.varMap:
                lineLoc['varName'] = varName
                axisLimitsLoc['varName'] = varName

                #subplot loop 2
                for binVal, binTitle in binValsMap:
                    lineLoc['binVal'] = binVal
                    axisLimitsLoc['binVal'] = binVal

                    # use common y-axis limits across axisLimitsLoc database locations
                    if statName == 'Count':
                        dmin = 0.
                    else:
                        dmin = dfwDict['dfw'].min(axisLimitsLoc, statName)
                    dmax = dfwDict['dfw'].max(axisLimitsLoc, statName)

                    # collect statName for all lines on this subplot
                    linesVals = []
                    linesLabel = []
                    linesGroup = []
                    for expName in self.expNames:
                        lineLoc['expName'] = expName
                        for diagnosticName in myLoc['diagName']:
                            linesGroup.append(expName)
                            linesLabel.append(self.expDiagnosticLabel(
                                expName, diagnosticName, myLoc['diagName']))

                            lineLoc['diagName'] = diagnosticName

                            lineCYDTimes = dfwDict['dfw'].levels('cyDTime', lineLoc)

                            lineVals = np.full(self.nCY, np.NaN)
                            cyLoc = deepcopy(lineLoc)
                            for cyDTime in lineCYDTimes:
                                icy = self.cyDTimes.index(cyDTime)
                                cyLoc['cyDTime'] = cyDTime
                                lineVals[icy] = dfwDict['dfw'].loc1(cyLoc, statName)
                            linesVals.append(lineVals)

                    # define subplot title
                    title = varLabel+binTitle

                    # perform subplot agnostic plotting (all expNames)
                    bpf.plotTimeSeries(
                        fig,
                        xVals, linesVals, linesLabel,
                        title, None, bgstatDiagLabel,
                        sciTicks, logScale, centralValue,
                        nyplots, nxplots, nsubplots, iplot,
                        dmin = dmin, dmax = dmax,
                        interiorLabels = interiorLabels)

                    iplot = iplot + 1

                # end binVal loop

            # end varMap loop

            # save each figure
            filename = ('%s%s_TSeries_%smin_%s_%s_%s'%(
                       myLoc['binVar'], self.binMethodFile(myLoc['binMethod']), fcTDelta_totmin,
                       self.DiagSpaceName, diagnosticGroup, statName))

            pu.finalize_fig(fig, str(figPath/filename), figureFileType, interiorLabels, 0.35, 0.55)

        # end fcMap loop


class FCAxisExpLines(CategoryBinMethodBase):
    '''
    Creates a timeseries figure between fcTDeltaFirst and fcTDeltaLast containing
      aggregated statistics for the period between firstCycleDTime and lastCycleDTime
      -  x-axis: forecast duration
      -    line: per experiment
      - subplot: combination of DiagSpace variable and binVal
      -    file: combination of binVar and statistic
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.requestAggDFW = True

        self.subplotWidth = 2.4
        self.subplotAspect = 0.75

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nFC < 2: return
        assert myLoc['binMethod'] is not None, 'FCAxisExpLines.innerloops: binMethod cannot be None'

        bgstatDiagLabel, fcstatDiagLabel, sciTicks, logScale, centralValue = \
            self.statPlotAttributes(diagnosticGroup, statName, myLoc['diagName'])

        fcDiagName = self.fcName(diagnosticGroup)
        figPath = self.myFigPath/fcDiagName
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        lineLoc = {}
        axisLimitsLoc = {}

        # establish a new figure
        fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, interiorLabels)
        iplot = 0

        #subplot loop 1
        for (varName, varLabel) in self.varMap:
            lineLoc['varName'] = varName
            axisLimitsLoc['varName'] = varName

            #subplot loop 2
            for binVal, binTitle in binValsMap:
                lineLoc['binVal'] = binVal
                axisLimitsLoc['binVal'] = binVal

                # use common y-axis limits across axisLimitsLoc database locations
                if statName == 'Count':
                    dmin = 0.
                else:
                    dmin = dfwDict['agg'].min(axisLimitsLoc, statName)
                dmax = dfwDict['agg'].max(axisLimitsLoc, statName)

                #collect aggregated statNames, varying across fcTDelta
                linesVals = []
                linesLabel = []
                linesGroup = []
                for expName in self.expNames:
                    lineLoc['expName'] = expName
                    for diagnosticName in myLoc['diagName']:
                        linesGroup.append(expName)
                        linesLabel.append(self.expDiagnosticLabel(
                            expName, diagnosticName, myLoc['diagName']))
                        lineLoc['diagName'] = diagnosticName

                        lineFCTDeltas = dfwDict['agg'].levels('fcTDelta', lineLoc)

                        lineVals = np.full(self.nFC, np.NaN)
                        fcLoc = deepcopy(lineLoc)
                        for fcTDelta in lineFCTDeltas:
                            ifc = self.fcTDeltas.index(fcTDelta)
                            fcLoc['fcTDelta'] = fcTDelta
                            lineVals[ifc] = dfwDict['agg'].loc1(fcLoc, statName)
                        linesVals.append(lineVals)

                # define subplot title
                title = varLabel+binTitle

                # perform subplot agnostic plotting (all expNames)
                bpf.plotTimeSeries(
                    fig,
                    self.fcTDeltas, linesVals, linesLabel,
                    title, 'Lead Time', fcstatDiagLabel,
                    sciTicks, logScale, centralValue,
                    nyplots, nxplots, nsubplots, iplot,
                    dmin = dmin, dmax = dmax,
                    interiorLabels = interiorLabels)

                iplot = iplot + 1

            # end statMap loop

        # end varMap loop

        # save figure
        filename = ('%s%s_TSeries_%s-%smin_%s_%s_%s'%(
                   myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                   self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                   self.DiagSpaceName, fcDiagName, statName))

        pu.finalize_fig(fig, str(figPath/filename), figureFileType, interiorLabels, 0.35, 0.55)


class FCAxisExpLinesDiffCI(CategoryBinMethodBase):
    '''
    Similar to FCAxisExpLines, except
      - shows difference between experiment(s) and control
      - control is selected using cntrlExpIndex
      - statistics are narrowed down by su.bootStrapStats
      - confidence intervals (CI) are shown at each lead time
      -    line+shaded region: per experiment
      - subplot: combination of DiagSpace variable and binVal
      -    file: combination of binVar and statistic
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        # OPTIONAL: implement fine-grained parallelism for bootStrapping
        #self.blocking = True

        # used for percent ratio plots
        self.requestAggDFW = True

        self.subplotWidth = 2.4
        self.subplotAspect = 0.75

        for key in self.binVarDict:
            if 'onlyStatNames' in self.binVarDict[key]:
                self.binVarDict[key]['onlyStatNames'] += su.bootStrapStats
            else:
                self.binVarDict[key]['onlyStatNames'] = su.bootStrapStats

        self.requiredStatistics = ['Count']

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nFC < 2: return
        if self.nExp * len(myLoc['diagName']) < 2: return
        if self.cntrlExpName not in dfwDict['dfw'].levels('expName'): return
        assert myLoc['binMethod'] is not None, 'FCAxisExpLinesDiffCI.innerloops: binMethod cannot be None'

        bgstatDiagLabel_abs, fcstatDiagLabel_abs, sciTicks_abs, logScale_abs, centralValue_abs = \
            self.statPlotAttributes(diagnosticGroup, statName, myLoc['diagName'], isDifferencePlot=True)

        fcDiagName = self.fcName(diagnosticGroup)
        figPath = self.myFigPath/fcDiagName
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        # Only bootstrap over the union of cyDTimes available
        # from both experiments at each fcTDelta
        myExpsCYDTimes = self.UNIONcntrlANDexpCYDTimes(dfwDict['dfw'])

        # establish a new figure
        fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, interiorLabels)
        iplot = 0

        useRelativeDifference = (
          statName in su.posSemiDefiniteStats and
          self.relativeErrorType != 'disable' and
          not set(myLoc['diagName']).issubset(du.absoluteOnlyDiagnostics))

        binValLoc = {}
        #subplot loop 1
        for (varName, varLabel) in self.varMap:
            binValLoc['varName'] = varName

            #subplot loop 2
            for binVal, binTitle in binValsMap:
                binValLoc['binVal'] = binVal

                # intermediate tempdfw reduces extraction time in inner loops
                tempdfw = sdb.DFWrapper.fromLoc(dfwDict['dfw'], binValLoc)
                normdfw = sdb.DFWrapper.fromLoc(dfwDict['agg'], binValLoc)

                cntrlLoc = deepcopy(binValLoc)
                cntrlLoc['expName'] = self.cntrlExpName
                cntrlLoc['diagName'] = myLoc['diagName'][0]

                normLoc = deepcopy(cntrlLoc)

                # define subplot title
                if useRelativeDifference:
                    title = varName
                else:
                    title = varLabel
                title += binTitle

                linesVals = defaultdict(list)
                linesLabel = []
                linesGroup = []
                for expName in self.expNames:
                    for diagnosticName in myLoc['diagName']:
                        if (expName == cntrlLoc['expName'] and
                            diagnosticName == cntrlLoc['diagName']): continue
                        linesGroup.append(expName)
                        linesLabel.append(self.expDiagnosticLabel(
                            expName, diagnosticName, myLoc['diagName']))

                        lineVals = defaultdict(list)
                        for fcTDelta in self.fcTDeltas:
                            cntrlLoc['cyDTime'] = myExpsCYDTimes[(expName, fcTDelta)]
                            cntrlLoc['fcTDelta'] = fcTDelta

                            expLoc = deepcopy(cntrlLoc)
                            expLoc['diagName'] = diagnosticName
                            expLoc['expName'] = expName

                            X = tempdfw.loc(expLoc)
                            Y = tempdfw.loc(cntrlLoc)

                            ciVals = su.bootStrapClusterFunc(
                                         X, Y,
                                         n_samples = 10000,
                                         statNames = [statName])

                            # normalizing value for ratio
                            normLoc['fcTDelta'] = fcTDelta
                            normalizingStat = normdfw.loc1(normLoc, statName)

                            for trait in su.ciTraits:
                                t = float(ciVals[statName][trait])
                                # automatically generate relative difference plots for positive-semi-definite statistics
                                if useRelativeDifference:
                                  # divide by cntrlLoc aggregated statName
                                  t /= normalizingStat
                                  if self.relativeErrorType == 'one hundred centered':
                                    t += 1.0
                                  t *= 100.0
                                lineVals[trait] += [t]

                        for trait in su.ciTraits:
                            linesVals[trait].append(lineVals[trait])

                # use specific y-axis limits for each varName
                #dmin = np.nanmin(linesVals[su.cimin])
                #dmax = np.nanmax(linesVals[su.cimax])
                dmin = np.nanmin(linesVals[su.cimean])
                dmax = np.nanmax(linesVals[su.cimean])

                if not (np.isfinite(dmin) or np.isfinite(dmax)):
                  iplot = iplot + 1
                  continue

                fcstatDiagLabel = fcstatDiagLabel_abs
                sciTicks = sciTicks_abs
                logScale = logScale_abs
                centralValue = 0.

                if useRelativeDifference:

                  fcstatDiagLabel = self.relativeErrorLabeler()
                  fcstatDiagLabel = statName.replace('RMS','rms').replace('Mean','mean')+': '+fcstatDiagLabel

                  #dmin, dmax = self.relativeErrorLimiter(dmin, dmax)

                  centralValue = self.relativeErrorCenter
                  sciTicks = False
                  logScale = False

                # perform subplot agnostic plotting (all expNames)
                bpf.plotTimeSeries(
                    fig,
                    self.fcTDeltas, linesVals[su.cimean],
                    linesLabel,
                    title, 'Lead Time', fcstatDiagLabel,
                    sciTicks, logScale, centralValue,
                    nyplots, nxplots, nsubplots, iplot,
                    linesValsMinCI = linesVals[su.cimin],
                    linesValsMaxCI = linesVals[su.cimax],
                    dmin = dmin, dmax = dmax,
                    lineAttribOffset = 1,
                    interiorLabels = interiorLabels)

                iplot = iplot + 1

            # end binValsMap loop

        # end varName loop

        # save figure
        filename = ('%s%s_TSeries_%s-%smin_%s_%s_%s'%(
                   myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                   self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                   self.DiagSpaceName, fcDiagName, statName))

        pu.finalize_fig(fig, str(figPath/filename), figureFileType, interiorLabels, 0.35, 0.55)


class FCScoreCard(CategoryBinMethodBase):
    '''
    Similar to FCAxisExpLinesDiffCI, except
      - instead of individual time-series plots, data is displayed in a 2D score card
        with a separate row in 2D subplot placing each FCAxisExpLinesDiffCI subplot axis
      - subplot: column by experiment, row by DiagSpace variable
      -    file: combination of binVar, binMethod, and statistic
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        # OPTIONAL: implement fine-grained parallelism for bootStrapping
        #self.blocking = True

        # default binVar/binMethod combinations
        # TODO: extend quads to dictionary with each key as subfigure title
        # TODO: extend quads to include varName key as alterantive to using separate list
        self.binVarDict = {
            ## vertical-latitudinal bins
            # model-space any
            (vu.modVarDiagPrs, None): {
              'binMethods': [
                ('NXTro', 'NX'),
                ('Tro', 'TR'),
                ('SXTro', 'SX'),
              ],
              'binVals': [
                ('50', '50'),
                ('100', '100'),
                ('200', '200'),
                #('250', '250'),
                ('500', '500'),
                #('700', '700'),
                ('850', '850'),
                #('925', '925'),
              ],
              'varNames': [
                'T', 'RH', 'U', 'V', 'height', #'Tdp',
              ],
            },
            # equivalent example using 'quads'
#            (vu.modVarDiagPrs, None): {
#              'quads': [
#                ('NXTro', 'NX', '50', '50'),
#                ('NXTro', 'NX', '100', '100'),
#                ('NXTro', 'NX', '200', '200'),
#                #('NXTro', 'NX', '250', '250'),
#                ('NXTro', 'NX', '500', '500'),
#                #('NXTro', 'NX', '700', '700'),
#                ('NXTro', 'NX', '850', '850'),
#                #('NXTro', 'NX', '925', '925'),
#                ('Tro', 'TR', '50', '50'),
#                ('Tro', 'TR', '100', '100'),
#                ('Tro', 'TR', '200', '200'),
#                #('Tro', 'TR', '250', '250'),
#                ('Tro', 'TR', '500', '500'),
#                #('Tro', 'TR', '700', '700'),
#                ('Tro', 'TR', '850', '850'),
#                #('Tro', 'TR', '925', '925'),
#                ('SXTro', 'SX', '50', '50'),
#                ('SXTro', 'SX', '100', '100'),
#                ('SXTro', 'SX', '200', '200'),
#                #('SXTro', 'SX', '250', '250'),
#                ('SXTro', 'SX', '500', '500'),
#                #('SXTro', 'SX', '700', '700'),
#                ('SXTro', 'SX', '850', '850'),
#                #('SXTro', 'SX', '925', '925'),
#              ],
#            },
            # observation-space pressure
            (vu.obsVarPrs, None): {
              'binMethods': [
                ('NXTro', 'NX'),
                ('Tro', 'TR'),
                ('SXTro', 'SX'),
              ],
              'binVals': [
                #('5.0', '5'),
                #('10.0', '10'),
                #('25.0', '25'),
                ('50.0', '50'),
                ('100.0', '100'),
                ('200.0', '200'),
                #('250.0', '250'),
                ('500.0', '500'),
                #('700.0', '700'),
                ('850.0', '850'),
                #('925.0', '925'),
                #('1000.0', '1000'),
              ],
              #'varNames': {},
            },
            # observation-space altitude
            (vu.obsVarAlt, None): {
              'binMethods': [
                ('NPol', 'NP'),
                ('NMid', 'NM'),
                ('Tro', 'TR'),
                ('SMid', 'SM'),
                ('SPol', 'SP'),
              ],
              'binVals': [
                ('28000', '28km'),
                #('25000', '25km'),
                ('20000', '20km'),
                ('15000', '15km'),
                ('10000', '10km'),
                ('5000', '5km'),
                ('2000', '2km'),
                ('1000', '1km'),
              ],
              #'varNames': {},
            },
            # observation-space impact height
            (vu.obsVarImpact, None): {
              'binMethods': [
                ('NPol', 'NP'),
                ('NMid', 'NM'),
                ('Tro', 'TR'),
                ('SMid', 'SM'),
                ('SPol', 'SP'),
              ],
              'binVals': [
                ('28000', '28km'),
                #('25000', '25km'),
                ('20000', '20km'),
                ('15000', '15km'),
                ('10000', '10km'),
                ('5000', '5km'),
                ('2000', '2km'),
                ('1000', '1km'),
              ],
              #'varNames': {},
            },

            ## vertically-aggregated latitude bins
            (vu.obsVarLat, bu.latbandsMethod): {
              'subplot index': 'varName',
              'binVals': [
                ('NXTro', 'NX'),
                ('Tro', 'TR'),
                ('SXTro', 'SX'),
              ],
            },
            (vu.obsVarLat, bu.troplatbandsMethod): {
              'subplot index': 'varName',
              'binVals': [
                ## 5-bin
                ('NXTro', 'NX'),
                ('NTro', 'NT'),
                ('ITCZ', 'ITCZ'),
                ('STro', 'ST'),
                ('SXTro', 'SX'),
                ## 3-bin
                #('NXTro', 'NX'),
                #('Tro', 'TR'),
                #('SXTro', 'SX'),
              ],
            },
            (vu.obsVarLat, bu.polarlatbandsMethod): {
              'subplot index': 'varName',
              'binVals': [
                ('NPol', 'NP'),
                ('NMid', 'NM'),
                ('Tro', 'TR'),
                ('SMid', 'SM'),
                ('SPol', 'SP'),
              ],
            },

        }

        # used for percent ratio plots
        self.requestAggDFW = True

        for key in self.binVarDict:
            if 'onlyStatNames' in self.binVarDict[key]:
                self.binVarDict[key]['onlyStatNames'] += su.bootStrapStats
            else:
                self.binVarDict[key]['onlyStatNames'] = su.bootStrapStats

        self.requiredStatistics = ['Count']

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
        a, b, c):

        availableVarNames = dfwDict['dfw'].levels('varName')
        varNames = options.get('varNames', [v for (v, l) in self.varMap])
        usableVarNames = list(set(varNames) & set(availableVarNames))
        varNames = [v for v in varNames if v in usableVarNames]
        if len(varNames) > 1:
          variables = pu.uniqueMembers(list(zip(varNames, varNames)))
        else:
          variables = pu.uniqueMembers(list(zip(varNames, [''])))

        # TODO: make units look right OR remove lines below
        #varUnitss = []
        #for varName in varNames:
        #    units = dfwDict['dfw'].uniquevals('varUnits', {'varName': varName})
        #    varUnitss.append(units[0])
        #varLabels = []
        #for (varName, varUnits) in zip(varNames, varUnitss):
        #    label = varName
        #    if varUnits != vu.miss_s:
        #        label = label+' ('+varUnits+')'
        #    for orig, sub in self.labelReplacements.items():
        #      label = label.replace(orig, sub)
        #    varLabels.append(label)
        #variables = pu.uniqueMembers(list(zip(varNames, varLabels)))

        sorters = {}
        sorters['varName'] = variables

        if 'quads' in options and myLoc['binMethod'] is None:
          quads = options['quads']

          indexGenerator = itt.product(variables, quads)

          varMap = []
          binMethods = []
          binVals = []
          for (var, quad) in indexGenerator:
            varMap.append(var)
            binMethods.append(tuple(quad[0:1]))
            binVals.append(tuple(quad[2:]))

          binMethodFile = self.binMethodFile('-'.join(pu.uniqueMembers(binMethods)))
        else:
          if myLoc['binMethod'] is None:
            methods = options['binMethods']
          else:
            methods = [(None, '')]

          vals = options['binVals']

          indexGenerator = itt.product(variables, methods, vals)
          varMap = []
          binMethods = []
          binVals = []
          for (var, method, val) in indexGenerator:
            varMap.append(var)
            binMethods.append(method)
            binVals.append(val)

        sorters['binMethod'] = pu.uniqueMembers(binMethods)
        sorters['binVal'] = pu.uniqueMembers(binVals)

        if myLoc['binMethod'] is None:
          binMethodFile = self.binMethodFile('-'.join([m[0] for m in sorters['binMethod']]))
        else:
          binMethodFile = self.binMethodFile(myLoc['binMethod'])

        ## subplot index: select the index that varies across subplots
        # OPTIONS: 'varName', 'binMethod'
        # TODO: allow for subplotIndex==None (one subplot for all combinations)
        subplotIndex = options.get('subplot index', 'binMethod')
        assert subplotIndex in ['varName', 'binMethod'], 'FCAxisExpLinesDiffCI.innerloops: invalid subplot index'

        indices = {
          'varName': varMap,
          'binMethod': binMethods,
          'binVal': binVals,
        }
        # sort on binVal first (maintains detail-level ordering)
        sortOrder = ['binVal']
        nRow = len(sorters['binVal'])

        # sort on non-subplotIndex indices next
        for key, sorter in sorters.items():
          if key != subplotIndex and key not in sortOrder:
            sortOrder.append(key)
            nRow *= len(sorter)

        # subplot rowIndex is identical to sortOrder without subplotIndex
        rowIndex = deepcopy(sortOrder)

        # sort on subplotIndex last
        sortOrder.append(subplotIndex)

        def customSortIndex(ll, sorter):
          order = {key: i for i, key in enumerate(sorter)}
          ii = list(range(len(ll)))
          ii.sort(key=lambda i: order[ll[i]])
          return ii


        # sort row values across all indices according to sortOrder
        allRowIndices = []
        for r in sortOrder:
          allRowIndices.append(indices[r])

        allRowIndices = list(zip(*allRowIndices))

        for sort in sortOrder:
          # reverse order to place top of list at top of each subfigure,
          # keeping subplot order the same
          sorter = deepcopy(sorters[sort])
          if sort != subplotIndex: sorter.reverse()

          # get index order for indices[sort] when sorted against sorter
          order = customSortIndex(indices[sort], sorter)

          # sort allRowIndices according to order
          allRowIndices = list(map(allRowIndices.__getitem__, order))

          # unzip
          allRowIndices = list(zip(*allRowIndices))

          # restore allRowIndices back to indices
          for jj, r in enumerate(sortOrder):
            indices[r] = allRowIndices[jj]

          # re-unzip
          allRowIndices = list(zip(*allRowIndices))

        # get list of all row names
        allRowNames = []
        for row in allRowIndices:
            name = ''
            for ii in range(len(rowIndex)):
                sub = row[ii][1]
                if ii > 0 and sub is not None and sub != '':
                  name +=','
                name += sub
            allRowNames.append(name)

        # axes settings (override a, b, c args above)
        nSubplotsPerExp = len(sorters[subplotIndex])
        if self.nExp > 2:
            nxplots = self.nExp-1
            nyplots = nSubplotsPerExp
            nsubplots = nxplots * nyplots
        else:
            nsubplots = nSubplotsPerExp
            nxplots = np.int(np.ceil(np.sqrt(nsubplots)))
            while nsubplots%nxplots > 0 and nsubplots%nxplots / nxplots <= 0.5: nxplots += 1
            nyplots = np.int(np.ceil(np.true_divide(nsubplots, nxplots)))

        xVals = self.fcTDeltas

        # assume tick labels are self explanatory and no x-label or y-label are needed
        #xLabel = 'Lead Time'
        #yLabel = myLoc['binVar']
        #for orig, sub in self.labelReplacements.items():
        #    yLabel = yLabel.replace(orig, sub)
        xLabel = None
        yLabel = None

        if self.nFC < 2: return
        if self.cntrlExpName not in dfwDict['dfw'].levels('expName'): return

        bgstatDiagLabel_abs, fcstatDiagLabel_abs, sciTicks_abs, logScale_abs, centralValue_abs = \
            self.statPlotAttributes(diagnosticGroup, statName, myLoc['diagName'], isDifferencePlot=True)

        fcDiagName = self.fcName(diagnosticGroup)
        figPath = self.myFigPath/fcDiagName
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        # Only bootstrap over the union of cyDTimes available
        # from both experiments at each fcTDelta
        myExpsCYDTimes = self.UNIONcntrlANDexpCYDTimes(dfwDict['dfw'])

        useRelativeDifference = (
          statName in su.posSemiDefiniteStats and
          self.relativeErrorType != 'disable' and
          not set(myLoc['diagName']).issubset(du.absoluteOnlyDiagnostics))

        # approximate componenent spacing (all in inches)
        boxEdge = 0.18

        # horizontal dimensions
        nCol = self.nFC
        colWidth = boxEdge
        colTotal = nCol*colWidth
        maxStrLength = max([len(str(s)) for s in allRowNames])
        yLabelWidth = 0.06*maxStrLength #approx. remainder
        outerWidth = yLabelWidth + 0.05*colTotal + 0.65
        slopWidth = 0.00
        subplotWidth = outerWidth + colTotal + slopWidth
        outerWidthFraction = outerWidth/subplotWidth

        # vertical dimensions
        rowHeight = boxEdge
        rowTotal = nRow*rowHeight
        outerHeight = 0.70
        slopHeight = 0.00
        subplotHeight = outerHeight + rowTotal + slopHeight
        outerHeightFraction = outerHeight/subplotHeight

        # establish a new figure
        subplotAspect = subplotHeight / subplotWidth
        fig = pu.setup_fig(nxplots, nyplots, subplotWidth, subplotAspect, True)
        iplot = 0

        planeLoc = {}
        # file loop
        for diagnosticName in myLoc['diagName']:
            planeLoc['diagName'] = diagnosticName

            # subplot loop 2
            for subplotName, subplotLabel in pu.uniqueMembers(indices[subplotIndex]):
                planeLoc[subplotIndex] = subplotName

                rowNames = []
                rowIndices = []
                for ii, row in enumerate(allRowIndices):
                    if row[-1][0] != subplotName: continue
                    rowNames.append(allRowNames[ii])
                    rowIndices.append(allRowIndices[ii])

                # intermediate tempdfw reduces extraction time in inner loops
                tempdfw = sdb.DFWrapper.fromLoc(dfwDict['dfw'], planeLoc)
                normdfw = sdb.DFWrapper.fromLoc(dfwDict['agg'], planeLoc)

                # subplot loop 1
                for expName in self.expNames:
                    if expName == self.cntrlExpName: continue

                    cntrlLoc = {}
                    cntrlLoc['expName'] = self.cntrlExpName

                    normLoc = {}
                    normLoc['expName'] = self.cntrlExpName

                    # define subplot title
                    # which is correct for non-varName subplotIndex?
                    #title = subplotName
                    title = subplotLabel
                    if self.nExp > 2:
                        title = expName+'\n'+title

                    planeVals = {}
                    for trait in su.ciTraits:
                        planeVals[trait] = np.full((nRow, nCol), np.NaN)

                    # row loop
                    for iy, row in enumerate(rowIndices):
                        for ii, index in enumerate(rowIndex):
                            if index != 'binMethod' or myLoc['binMethod'] is None:
                                cntrlLoc[index] = row[ii][0]
                                normLoc[index] = row[ii][0]

                        # column loop
                        for ix, fcTDelta in enumerate(self.fcTDeltas):
                            cntrlLoc['cyDTime'] = myExpsCYDTimes[(expName, fcTDelta)]
                            cntrlLoc['fcTDelta'] = fcTDelta

                            expLoc = deepcopy(cntrlLoc)
                            expLoc['expName'] = expName

                            X = tempdfw.loc(expLoc)
                            Y = tempdfw.loc(cntrlLoc)

                            ciVals = su.bootStrapClusterFunc(
                                         X, Y,
                                         n_samples = 10000,
                                         statNames = [statName])

                            # normalizing value for ratio
                            normLoc['fcTDelta'] = fcTDelta
                            normalizingStat = normdfw.loc1(normLoc, statName)

                            for trait in su.ciTraits:
                                t = float(ciVals[statName][trait])
                                # automatically generate relative difference plots for positive-semi-definite statistics
                                if useRelativeDifference:
                                  # divide by cntrlLoc aggregated statName
                                  t /= normalizingStat
                                  if self.relativeErrorType == 'one hundred centered':
                                    t += 1.0
                                  t *= 100.0
                                planeVals[trait][iy, ix] = t

                    # use specific y-axis limits for each subplot
                    dmin = np.nanmin(planeVals[su.cimean])
                    dmax = np.nanmax(planeVals[su.cimean])

                    if not (np.isfinite(dmin) or np.isfinite(dmax)):
                      iplot = iplot + 1
                      continue

                    fcstatDiagLabel = fcstatDiagLabel_abs
                    sciTicks = sciTicks_abs
                    logScale = logScale_abs
                    centralValue = 0.

                    if useRelativeDifference:

                      fcstatDiagLabel = self.relativeErrorLabeler()
                      fcstatDiagLabel = statName.replace('RMS','rms').replace('Mean','mean')+': '+fcstatDiagLabel

                      #dmin, dmax = self.relativeErrorLimiter(dmin, dmax)

                      centralValue = self.relativeErrorCenter
                      sciTicks = False
                      logScale = False

                    # perform subplot agnostic plotting
                    bpf.scoreCard(
                        fig,
                        xVals, rowNames, planeVals[su.cimean],
                        title, xLabel, yLabel, fcstatDiagLabel,
                        bpf.defaultIndepConfig,
                        sciTicks, logScale, centralValue,
                        nyplots, nxplots, nsubplots, iplot,
                        contourValsMinCI = planeVals[su.cimin],
                        contourValsMaxCI = planeVals[su.cimax],
                        dmin = dmin, dmax = dmax)

                    iplot = iplot + 1

                # end subplotIndex loop

            # end expName loop

            # save figure
            filename = ('%s%s_TSeries_%s-%smin_%s_%s_%s'%(
                       myLoc['binVar'], binMethodFile,
                       self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                       self.DiagSpaceName, diagnosticName, statName))

        # end diagnosticName loop

        pu.finalize_fig(fig, str(figPath/filename), figureFileType, True, outerWidthFraction, outerHeightFraction)


class CYAxisFCLines(CategoryBinMethodBase):
    '''
    Similar to CYAxisExpLines, except
      each line is for a different forecast lead time and
      each experiment is in a different file
      -  x-axis: valid time of forecast
      -    line: per FC lead time
      - subplot: combination of DiagSpace variable and binVal
      -    file: combination of binVar, statistic, and experiment
      - self.MAX_FC_LINES determines number of FC lead time lines to include
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.subplotWidth = 2.4
        self.subplotAspect = 0.75

        self.maxDiagnosticsPerAnalysis = 1

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nFC < 2 or self.nCY < 2: return
        assert myLoc['binMethod'] is not None, 'CYAxisFCLines.innerloops: binMethod cannot be None'

        bgstatDiagLabel, fcstatDiagLabel, sciTicks, logScale, centralValue = \
            self.statPlotAttributes(diagnosticGroup, statName)

        fcDiagName = self.fcName(diagnosticGroup)
        figPath = self.myFigPath/fcDiagName
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        lineLoc = {}
        axisLimitsLoc = {}

        #file loop 1
        for expName in self.expNames:
            lineLoc['expName'] = expName

            # establish a new figure
            fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, interiorLabels)
            iplot = 0

            #subplot loop 1
            for (varName, varLabel) in self.varMap:
                lineLoc['varName'] = varName
                axisLimitsLoc['varName'] = varName

                #subplot loop 2
                for binVal, binTitle in binValsMap:
                    lineLoc['binVal'] = binVal
                    axisLimitsLoc['binVal'] = binVal

                    # use common y-axis limits across axisLimitsLoc database locations
                    if statName == 'Count':
                        dmin = 0.
                    else:
                        dmin = dfwDict['dfw'].min(axisLimitsLoc, statName)
                    dmax = dfwDict['dfw'].max(axisLimitsLoc, statName)

                    # collect statName for all lines on this subplot, letting cyDTime vary
                    xsVals = []
                    linesVals = []
                    self.fcTDeltas_labels = []
                    for fcTDelta in self.fcTDeltas:
                        lineLoc['fcTDelta'] = fcTDelta

                        # calculate valid time for x-axis
                        xVals = []
                        for cyDTime in self.cyDTimes:
                            xVals.append(cyDTime+fcTDelta)
                        xsVals.append(xVals)

                        #Setting to avoid over-crowding
                        if self.fcTDeltas.index(fcTDelta) > (self.MAX_FC_LINES-1): continue

                        self.fcTDeltas_labels.append(
                            pu.timeDeltaTicks(fcTDelta.total_seconds(),0))

                        lineCYDTimes = dfwDict['dfw'].levels('cyDTime', lineLoc)

                        lineVals = np.full(self.nCY, np.NaN)
                        cyLoc = deepcopy(lineLoc)
                        for cyDTime in lineCYDTimes:
                            icy = self.cyDTimes.index(cyDTime)
                            cyLoc['cyDTime'] = cyDTime
                            lineVals[icy] = dfwDict['dfw'].loc1(cyLoc, statName)

                        linesVals.append(lineVals)

                    # define subplot title
                    title = varLabel+binTitle

                    # perform subplot agnostic plotting (all expNames)
                    bpf.plotTimeSeries(
                        fig,
                        xsVals, linesVals, self.fcTDeltas_labels,
                        title, None, bgstatDiagLabel,
                        sciTicks, logScale, centralValue,
                        nyplots, nxplots, nsubplots, iplot,
                        dmin = dmin, dmax = dmax,
                        interiorLabels = interiorLabels)

                    iplot = iplot + 1

                # end binValsMap loop

            # end varMap loop

            expFileName = re.sub('\.', '', re.sub('\s+', '-', expName))
            filename = ('%s%s_TSeries_%s_%s_%s_%s'%(
                       myLoc['binVar'], self.binMethodFile(myLoc['binMethod']), expFileName,
                       self.DiagSpaceName, fcDiagName, statName))

            pu.finalize_fig(fig, str(figPath/filename), figureFileType, interiorLabels, 0.35, 0.55)

        # end expName loop


###########################################
## Figures with individual lines per binVal
###########################################

class BinValLinesAnalysisType(CategoryBinMethodBase):
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.binVarDict = {
            (vu.obsVarQC, bu.badQCMethod): {
                'onlyStatNames': ['Count'],
            },
            (vu.obsVarQC, bu.allQCMethod): {
                'onlyStatNames': ['Count'],
            },
            (vu.obsVarLat, bu.latbandsMethod): {},
            (vu.obsVarLat, bu.troplatbandsMethod): {},
            # vu.modVarLat is redundant with vu.obsVarLat (both have varShort=="lat")
            #(vu.modVarLat, bu.latbandsMethod): {'binVarTier': 1},
            (vu.obsVarCldFracX, bu.cloudbandsMethod): {},
            (vu.obsVarCldFracY, bu.cloudbandsMethod): {},
            (vu.obsVarLandFrac, bu.surfbandsMethod): {'binVarTier': 3},
        }

    def subplotArrangement(self, dummy):
        # subplot configuration
        return self.nExp, self.nVars, self.nExp * self.nVars


class CYAxisBinValLines(BinValLinesAnalysisType):
    '''
    Similar to CYAxisExpLines, except
      each line is for a different binVal (e.g., latitude band, cloudiness, etc.)
      -    line: binVals for named bins (e.g., NXTro, Tro, SXTro for latitude)
      - subplot: column by experiment, row by DiagSpace variable
      -    file: combination of statistic and forecast length
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.subplotWidth = 2.4
        self.subplotAspect = 0.75

        self.maxDiagnosticsPerAnalysis = 1

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nCY < 2: return

        bgstatDiagLabel, fcstatDiagLabel, sciTicks, logScale, centralValue = \
            self.statPlotAttributes(diagnosticGroup, statName)


        figPath = self.myFigPath/diagnosticGroup
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        lineLoc = {}
        binStrVals = [binVal for binVal, binTitle in binValsMap]
        lineLoc['binVal'] = binStrVals

        axisLimitsLoc = deepcopy(lineLoc)

        #file loop 1
        for (fcTDelta, fcTDelta_totmin) in self.fcMap:
            # calculate valid time for x-axis
            xVals = []
            for cyDTime in self.cyDTimes:
                xVals.append(cyDTime+fcTDelta)

            lineLoc['fcTDelta'] = fcTDelta
            axisLimitsLoc['fcTDelta'] = fcTDelta

            # establish a new figure
            fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, interiorLabels)

            iplot = 0

            #subplot loop 1
            for (varName, varLabel) in self.varMap:
                lineLoc['varName'] = varName
                axisLimitsLoc['varName'] = varName

                # use common y-axis limits across axisLimitsLoc database locations
                if statName == 'Count':
                    dmin = 0.
                else:
                    dmin = dfwDict['dfw'].min(axisLimitsLoc, statName)
                dmax = dfwDict['dfw'].max(axisLimitsLoc, statName)

                #subplot loop 2
                for expName in self.expNames:
                    lineLoc['expName'] = expName

                    # collect statName for all lines on this subplot, letting cyDTime vary
                    linesVals = []
                    for binVal in binStrVals:
                        lineLoc['binVal'] = binVal
                        lineCYDTimes = dfwDict['dfw'].levels('cyDTime', lineLoc)

                        lineVals = np.full(self.nCY, np.NaN)
                        cyLoc = deepcopy(lineLoc)
                        for cyDTime in lineCYDTimes:
                            icy = self.cyDTimes.index(cyDTime)
                            cyLoc['cyDTime'] = cyDTime
                            lineVals[icy] = dfwDict['dfw'].loc1(cyLoc, statName)

                        linesVals.append(lineVals)

                    # end binVal loop

                    # define subplot title
                    title = expName+'\n'+varLabel

                    # perform subplot agnostic plotting (all expNames)
                    bpf.plotTimeSeries(
                        fig,
                        xVals, linesVals, binStrVals,
                        title, None, bgstatDiagLabel,
                        sciTicks, logScale, centralValue,
                        nyplots, nxplots, nsubplots, iplot,
                        dmin = dmin, dmax = dmax,
                        interiorLabels = interiorLabels)

                    iplot = iplot + 1

                # end expName Loop

            # end varMap Loop

            filename = ('%s%s_TSeries_%smin_%s_%s_%s'%(
                       myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                       fcTDelta_totmin, self.DiagSpaceName,
                       diagnosticGroup, statName))

            pu.finalize_fig(fig, str(figPath/filename), figureFileType, interiorLabels, 0.35, 0.55)

        # end fcMap loop


# TODO(JJG): implement FCAxisBinValLines similar to FCAxisExpLines


#########################################################
## Figures with binVal on one axis, i.e., 2D and profiles
#########################################################
class MultiDimBinMethodBase(AnalysisBase):
    '''
    Base class used to analyze statistics across binMethods with one-dimensional binValues
      that are assigned numerical values, e.g., altitude, pressure, latitude, cloud fraction
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)
        self.parallelism = True
        self.maxBinVarTier = 2

        # default 1D binVars
        self.binVarDict = {
            vu.obsVarAlt: {
              'profilefunc': bpf.plotProfile,
              'binFilter': {
                # maximum altitude to show on all figures
                'maxvalue': 30000.,
              },
            },
            vu.obsVarImpact: {
              'profilefunc': bpf.plotProfile,
              'binFilter': {
                # maximum altitude to show on all figures
                'maxvalue': 30000.,
              },
            },
            vu.obsVarACI: {'profilefunc': bpf.plotSeries, 'binVarTier': 3},
            vu.obsVarCldFracX: {'profilefunc': bpf.plotSeries, 'binVarTier': 2},
            vu.obsVarCldFracY: {'profilefunc': bpf.plotSeries, 'binVarTier': 1},
            vu.obsVarLat: {'profilefunc': bpf.plotProfile},
            vu.obsVarPrs: {'profilefunc': bpf.plotProfile},
            vu.obsVarCI: {'profilefunc': bpf.plotSeries, 'binVarTier': 2},
            vu.obsVarLogCI: {'profilefunc': bpf.plotSeries, 'binVarTier': 3},
            vu.modVarDiagPrs: {'profilefunc': bpf.plotProfile},
            # vu.modVarLat is redundant with vu.obsVarLat (both have varShort=="lat")
            #vu.modVarLat: {'profilefunc': bpf.plotProfile, 'binVarTier': 1},
            vu.modVarLev: {
              'profilefunc': bpf.plotProfile,
              'binFilter': {
                # maximum model level to show on all figures
                #'maxvalue': 40,
              },
            },
            vu.obsVarGlint: {'profilefunc': bpf.plotSeries, 'binVarTier': 3},
            vu.obsVarLandFrac: {'profilefunc': bpf.plotSeries, 'binVarTier': 3},
            vu.obsVarLT: {'profilefunc': bpf.plotSeries, 'binVarTier': 3},
            vu.obsVarSenZen: {'profilefunc': bpf.plotSeries, 'binVarTier': 3},
        }
        self.maxDiagnosticsPerAnalysis = 10 // self.nExp

    def analyze_(self, workers = None):
        useWorkers = (not self.blocking and self.parallelism and workers is not None)
        diagnosticGrouped = {}
        for diag in self.availableDiagnostics:
            diagnosticGrouped[diag] = False

        diagnosticGroupings = deepcopy(self.diagnosticGroupings)
        for group in list(diagnosticGroupings.keys()):
            diags = diagnosticGroupings[group]
            if (len(diags) > self.maxDiagnosticsPerAnalysis or
               not set(diags).issubset(set(list(self.availableDiagnostics)))):
                del diagnosticGroupings[group]
                continue
            for diag in diags: diagnosticGrouped[diag] = True

        for diag in self.availableDiagnostics:
            if not diagnosticGrouped[diag]:
                diagnosticGroupings[diag] = [diag]

        for diagnosticGroup, diagnosticNames in diagnosticGroupings.items():
            if len(diagnosticNames) > self.maxDiagnosticsPerAnalysis: continue
            if len(set(diagnosticNames) & set(self.availableDiagnostics)) == 0: continue
            diagnosticConfigs = {}
            selectedStatistics = set([])
            for diagnosticName in diagnosticNames:
                diagnosticConfigs[diagnosticName] = deepcopy(self.diagnosticConfigs[diagnosticName])
                selectedStatistics = set(list(selectedStatistics) +
                                         diagnosticConfigs[diagnosticName]['selectedStatistics'])
            availableStatistics = set([])
            for diagnosticName in diagnosticNames:
                diagnosticConfigs[diagnosticName] = deepcopy(self.diagnosticConfigs[diagnosticName])
                availableStatistics = set(list(availableStatistics) +
                                         diagnosticConfigs[diagnosticName]['availableStatistics'])
            if not set(self.requiredStatistics).issubset(availableStatistics): continue

            diagBinVars = self.db.dfw.levels('binVar', {'diagName': diagnosticNames})

            for fullBinVar, options in self.binVarDict.items():
                if options.get('binVarTier', 1) > self.maxBinVarTier: continue
                binVar = vu.varDictAll.get(fullBinVar, [None, fullBinVar])[1]
                if (binVar not in diagBinVars): continue
                binVarLoc = {}
                binVarLoc['diagName'] = diagnosticNames
                binVarLoc['binVar'] = binVar
                binVarLoc['binVal'] = self.allBinNumVals2DasStr

                #Make figures for all binMethods
                binMethods = self.db.dfw.levels('binMethod', binVarLoc)
                for binMethod in binMethods:

                    #TODO: REMOVE, for testing only
                    #if binMethod != bu.identityBinMethod: continue

                    self.logger.info(diagnosticGroup+', '+binVar+', '+binMethod)

                    if useWorkers:
                        workers.apply_async(self.innerloopsWrapper,
                            args = (diagnosticGroup, diagnosticConfigs, binVar, binMethod, selectedStatistics, options))
                    else:
                        self.innerloopsWrapper(
                            diagnosticGroup, diagnosticConfigs, binVar, binMethod, selectedStatistics, options)

    def innerloopsWrapper(self,
        diagnosticGroup, diagnosticConfigs, binVar, binMethod, selectedStatistics, options):

        myLoc = {}
        myLoc['binVar'] = binVar
        myLoc['binVal'] = self.allBinNumVals2DasStr
        myLoc['binMethod'] = binMethod

        # narrow mydfwDict by binVar, binVal, and binMethod to reduce run-time and memory
        mydfwDict = {'dfw': self.db.loc(myLoc)}

        # aggregate statistics when requested
        if self.requestAggDFW:
            mydfwDict['agg'] = sdb.DFWrapper.fromAggStats(mydfwDict['dfw'], ['cyDTime'])
            sdb.createORreplaceDerivedDiagnostics(mydfwDict['agg'], diagnosticConfigs)

        # further narrow mydfwDict by diagName
        # NOTE: derived diagnostics may require multiple diagName values;
        # can only narrow by diagName after aggregation
        myLoc['diagName'] = list(diagnosticConfigs.keys())
        for key in mydfwDict.keys():
            mydfwDict[key] = sdb.DFWrapper.fromLoc(mydfwDict[key], myLoc)

        ## Get all float/int binVals associated with binVar
        binStrVals = mydfwDict['dfw'].levels('binVal')
        binUnits = mydfwDict['dfw'].uniquevals('binUnits')[0]

        # assume all bins represent same variable/units
        binLabel = binVar
        if binUnits != vu.miss_s:
            binLabel += ' ('+binUnits+')'
        for orig, sub in self.labelReplacements.items():
            binLabel = binLabel.replace(orig, sub)

        # bin info
        binNumVals = []
        for binVal in binStrVals:
            ibin = self.allBinStrVals.index(binVal)
            binNumVals.append(self.allBinNumVals[ibin])

        # filter out binVals less than (greater than) minvalue (maxvalue)
        binFilter = options.get('binFilter', None)
        # Example to include only values from 3 to 40:
        # 'binFilter': {
        #   'minvalue': 3,
        #   'maxvalue': 40,
        # }
        if binFilter is not None:
          remove = np.full_like(binNumVals, False, bool)
          minvalue = binFilter.get('minvalue', None)
          if minvalue is not None:
            less = bu.lessBound(np.asarray(binNumVals), minvalue)
            remove[less] = True
          maxvalue = binFilter.get('maxvalue', None)
          if maxvalue is not None:
            great = bu.greatBound(np.asarray(binNumVals), maxvalue)
            remove[great] = True

          binStrVals = list(np.asarray(binStrVals)[~remove])
          binNumVals = list(np.asarray(binNumVals)[~remove])

        # special independent variable axis configs
        binVarIs = {}
        specialBinVars = [
          vu.obsVarPrs,
          vu.obsVarMCI,
          vu.obsVarOCI,
          vu.obsVarLogCI,
          vu.obsVarCldFracX,
          vu.obsVarCldFracY,
          vu.modVarDiagPrs,
        ]
        for var in specialBinVars:
            var_dict = vu.varDictAll.get(var,['',''])
            binVarIs[var] = (var_dict[1] == binVar)

        binConfig = deepcopy(bpf.defaultIndepConfig)
        pCoord = binVarIs[vu.obsVarPrs] or binVarIs[vu.modVarDiagPrs]
        binConfig['invert'] = pCoord
        if pCoord: binConfig['transform'] = 'Pressure'
        if binVarIs[vu.obsVarMCI] or binVarIs[vu.obsVarOCI] or binVarIs[vu.obsVarLogCI]:
            binConfig['transform'] = 'CloudImpact'
        if binVarIs[vu.obsVarCldFracX] or binVarIs[vu.obsVarCldFracY]:
            binConfig['transform'] = 'logit'

        # sort bins by numeric value
        indices = list(range(len(binNumVals)))
        indices.sort(key=binNumVals.__getitem__)
        binNumVals = list(map(binNumVals.__getitem__, indices))
        binStrVals = list(map(binStrVals.__getitem__, indices))

        myBinConfigs = {
            'str': binStrVals,
            'num': binNumVals,
            'binLabel': binLabel,
            'binConfig': binConfig,
        }
        if len(binStrVals) < 2: return

        # only analyze variables that have non-zero Count when sliced by myLoc
        nVarsLoc = 0
        varMapLoc = []
        for (varName, varLabel) in self.varMap:
          if 'Count' in selectedStatistics:
              countDF = mydfwDict['dfw'].loc({'varName': varName}, 'Count')
              if countDF.shape[0] > 0:
                if np.nansum(countDF.to_numpy()) > 0:
                  nVarsLoc += 1
                  varMapLoc.append((varName, varLabel))
          else:
              statDF = mydfwDict['dfw'].loc({'varName': varName}, list(selectedStatistics)[0])
              if statDF.shape[0] > 0:
                if np.isfinite(statDF.to_numpy()).sum() > 0:
                  nVarsLoc += 1
                  varMapLoc.append((varName, varLabel))

        for statName in selectedStatistics:
            if statName not in options.get('onlyStatNames', selectedStatistics): continue

            self.innerloops(
                mydfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options)

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options):
        '''
        virtual method
        '''
        raise NotImplementedError()


class CYandBinValAxes2D(MultiDimBinMethodBase):
    '''
    Creates raster maps with binVar binVals on y-axis
      - only applicable to binned diagnostics (e.g., vertical dimension, latitude, zenith angle)
      - subplot: column by experiment, row by DiagSpace variable
      -    file: combination of binVar, binMethod, statistic, and FC lead time
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.subplotWidth = 2.5
        self.subplotAspect = 0.70

        self.maxDiagnosticsPerAnalysis = 1

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options):

        if self.nCY < 2: return

        bgstatDiagLabel_abs, fcstatDiagLabel_abs, sciTicks_abs, logScale_abs, centralValue_abs = \
            self.statPlotAttributes(diagnosticGroup, statName)

        figPath = self.myFigPath/diagnosticGroup
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        # axes settings
        nxplots = self.nExp
        nyplots = nVarsLoc
        nsubplots = nxplots * nyplots

        xLabel = None

        binNumVals = myBinConfigs['num']
        binStrVals = myBinConfigs['str']
        nBinVals = len(binStrVals)
        binLabel = myBinConfigs['binLabel']

        planeLoc = {}
        axisLimitsLoc = {}

        useRelativeDifference = (
          statName in su.posSemiDefiniteStats and
          self.relativeErrorType != 'disable' and
          not set(myLoc['diagName']).issubset(du.absoluteOnlyDiagnostics))

        #file loop 1
        for (fcTDelta, fcTDelta_totmin) in self.fcMap:
            # calculate valid time for x-axis
            xVals = []
            for cyDTime in self.cyDTimes:
                xVals.append(cyDTime+fcTDelta)

            planeLoc['fcTDelta'] = fcTDelta

            # establish a new figure
            fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, interiorLabels)

            iplot = 0

            #subplot loop 1
            for (varName, varLabel) in varMapLoc:
                planeLoc['varName'] = varName
                axisLimitsLoc['varName'] = varName

                # use common c-axis limits across axisLimitsLoc database locations
                if statName == 'Count':
                    dmin_abs = 0.
                else:
                    dmin_abs = dfwDict['dfw'].min(axisLimitsLoc, statName)
                dmax_abs = dfwDict['dfw'].max(axisLimitsLoc, statName)

                # letting cyDTime and binVal vary
                # extract control experiment
                cntrlLoc = deepcopy(planeLoc)
                cntrlLoc['expName'] = self.cntrlExpName
                cntrlPlaneCYDTimes = dfwDict['dfw'].levels('cyDTime', cntrlLoc)
                cntrlPlaneVals = np.full((nBinVals, self.nCY), np.NaN)

                for ibin, binVal in enumerate(binStrVals):
                    cntrlLoc['binVal'] = binVal
                    tmp = dfwDict['dfw'].loc(cntrlLoc, statName).to_numpy()
                    for jcy, cyDTime in enumerate(cntrlPlaneCYDTimes):
                        if jcy > len(tmp)-1: continue
                        icy = self.cyDTimes.index(cyDTime)
                        cntrlPlaneVals[ibin, icy] = tmp[jcy]

                #subplot loop 2
                dmin_relative = np.NaN
                dmax_relative = np.NaN
                for expName in self.expNames:
                    expLoc = deepcopy(planeLoc)
                    expLoc['expName'] = expName

                    # define subplot title
                    if useRelativeDifference:
                        title = varName
                    else:
                        title = varLabel
                    title = expName+'\n'+title

                    bgstatDiagLabel = bgstatDiagLabel_abs
                    sciTicks = sciTicks_abs
                    logScale = logScale_abs
                    centralValue = centralValue_abs
                    dmin = dmin_abs
                    dmax = dmax_abs
                    if (expName == cntrlLoc['expName']):
                        planeVals = deepcopy(cntrlPlaneVals)
                    else:
                        # letting cyDTime and binVal vary
                        # extract this experiment
                        expPlaneCYDTimes = dfwDict['dfw'].levels('cyDTime', expLoc)
                        expPlaneVals = np.full_like(cntrlPlaneVals, np.NaN)

                        for ibin, binVal in enumerate(binStrVals):
                            expLoc['binVal'] = binVal
                            tmp = dfwDict['dfw'].loc(expLoc, statName).to_numpy()
                            for jcy, cyDTime in enumerate(expPlaneCYDTimes):
                                if jcy > len(tmp)-1: continue
                                icy = self.cyDTimes.index(cyDTime)
                                expPlaneVals[ibin, icy] = tmp[jcy]

                        # automatically generate relative difference plots for positive-semi-definite statistics
                        if useRelativeDifference:
                            sciTicks = False
                            logScale = False

                            planeVals, dmin_relative, dmax_relative, centralValue, label = self.relativeErrorFunction(
                              expPlaneVals,
                              cntrlPlaneVals,
                              dmin_relative,
                              dmax_relative,
                            )
                            dmin = dmin_relative
                            dmax = dmax_relative
                            bgstatDiagLabel = statName.replace('RMS','rms').replace('Mean','mean')+': '+label
                        else:
                            planeVals = deepcopy(expPlaneVals)

                    cLabel = bgstatDiagLabel

                    # perform subplot agnostic plotting (all expNames)
                    bpf.plot2D(
                        fig,
                        xVals, binNumVals, planeVals,
                        title, xLabel, binLabel, cLabel,
                        bpf.defaultIndepConfig,
                        myBinConfigs['binConfig'],
                        sciTicks, logScale, centralValue,
                        nyplots, nxplots, nsubplots, iplot,
                        dmin = dmin, dmax = dmax,
                        interiorLabels = interiorLabels)

                    iplot = iplot + 1

            filename = ('%s%s_BinValAxisTSeries_%smin_%s_%s_%s'%(
                       myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                       fcTDelta_totmin, self.DiagSpaceName,
                       diagnosticGroup, statName))

            pu.finalize_fig(fig, str(figPath/filename), figureFileType, interiorLabels, 0.35, 0.55)

        # end fcTDelta loop


class FCandBinValAxes2D(MultiDimBinMethodBase):
    '''
    Creates raster maps with binVar binVals on y-axis
      - only applicable to binned diagnostics (e.g., vertical dimension, latitude, zenith angle)
      - subplot: column by experiment, row by DiagSpace variable
      -    file: combination of binVar, binMethod, and statistic
    '''

    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.requestAggDFW = True

        self.subplotWidth = 2.6
        self.subplotAspect = 0.75

        self.maxDiagnosticsPerAnalysis = 1

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options):

        if self.nFC < 2: return

        bgstatDiagLabel_abs, fcstatDiagLabel_abs, sciTicks_abs, logScale_abs, centralValue_abs = \
            self.statPlotAttributes(diagnosticGroup, statName)

        fcDiagName = self.fcName(diagnosticGroup)
        figPath = self.myFigPath/fcDiagName
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        # axes settings
        nxplots = self.nExp
        nyplots = nVarsLoc
        nsubplots = nxplots * nyplots

        xVals = self.fcTDeltas
        xLabel = 'Lead Time'

        binNumVals = myBinConfigs['num']
        binStrVals = myBinConfigs['str']
        nBinVals = len(binStrVals)
        binLabel = myBinConfigs['binLabel']


        planeLoc = {}
        axisLimitsLoc = {}

        useRelativeDifference = (
          statName in su.posSemiDefiniteStats and
          self.relativeErrorType != 'disable' and
          not set(myLoc['diagName']).issubset(du.absoluteOnlyDiagnostics))

        # Only bootstrap over the union of cyDTimes available
        # from both experiments at each fcTDelta
        myExpsCYDTimes = self.UNIONcntrlANDexpCYDTimes(dfwDict['dfw'])

        # establish a new figure
        fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, interiorLabels)

        iplot = 0
        #subplot loop 1
        for (varName, varLabel) in varMapLoc:
            planeLoc['varName'] = varName
            axisLimitsLoc['varName'] = varName

            # use common c-axis limits across axisLimitsLoc database locations
            if statName == 'Count':
                dmin_abs = 0.
            else:
                dmin_abs = dfwDict['agg'].min(axisLimitsLoc, statName)
            dmax_abs = dfwDict['agg'].max(axisLimitsLoc, statName)

            # letting fcTDelta and binVal vary
            # extract control experiment
            cntrlAggLoc = deepcopy(planeLoc)
            cntrlAggLoc['expName'] = self.cntrlExpName
            cntrlPlaneFCTDeltas = dfwDict['agg'].levels('fcTDelta', cntrlAggLoc)
            cntrlAggPlaneVals = np.full((nBinVals, self.nFC), np.NaN)

            for ibin, binVal in enumerate(binStrVals):
                cntrlAggLoc['binVal'] = binVal
                aggFCStats = dfwDict['agg'].loc(cntrlAggLoc, statName).to_numpy()
                for jfc, fcTDelta in enumerate(cntrlPlaneFCTDeltas):
                    if jfc > len(aggFCStats)-1: continue
                    ifc = self.fcTDeltas.index(fcTDelta)
                    cntrlAggPlaneVals[ibin, ifc] = aggFCStats[jfc]

            # intermediate tempdfw reduces extraction time in inner loops
            tempdfw = sdb.DFWrapper.fromLoc(dfwDict['dfw'], planeLoc)

            #subplot loop 2
            dmin_relative = np.NaN
            dmax_relative = np.NaN
            for expName in self.expNames:
                # define subplot title
                if useRelativeDifference:
                    title = varName
                else:
                    title = varLabel
                title = expName+'\n'+title

                fcstatDiagLabel = fcstatDiagLabel_abs
                sciTicks = sciTicks_abs
                logScale = logScale_abs
                centralValue = centralValue_abs
                dmin = dmin_abs
                dmax = dmax_abs

                planeVals = {}
                for trait in su.ciTraits:
                    planeVals[trait] = np.full_like(cntrlAggPlaneVals, np.NaN)

                if (expName == cntrlAggLoc['expName']):
                    planeVals[su.cimean] = deepcopy(cntrlAggPlaneVals)
                else:
                    expAggLoc = deepcopy(cntrlAggLoc)
                    expAggLoc['expName'] = expName

                    cntrlLoc = deepcopy(cntrlAggLoc)
                    expLoc = deepcopy(cntrlLoc)
                    expLoc['expName'] = expName

                    # letting fcTDelta and binVal vary
                    # extract this experiment
                    expPlaneFCTDeltas = dfwDict['agg'].levels('fcTDelta', expAggLoc)
                    expAggPlaneVals = np.full_like(cntrlAggPlaneVals, np.NaN)
                    for ibin, binVal in enumerate(binStrVals):
                        iy = ibin

                        cntrlLoc['binVal'] = binVal
                        expLoc['binVal'] = binVal

                        expAggLoc['binVal'] = binVal
                        aggFCStats = dfwDict['agg'].loc(expAggLoc, statName).to_numpy()

                        for jfc, fcTDelta in enumerate(expPlaneFCTDeltas):
                            if jfc > len(aggFCStats)-1: continue
                            ifc = self.fcTDeltas.index(fcTDelta)
                            ix = ifc
                            expAggPlaneVals[iy, ix] = aggFCStats[jfc]

                            # normalizing value for ratio
                            normalizingStat = cntrlAggPlaneVals[iy, ix]

                            cntrlLoc['cyDTime'] = myExpsCYDTimes[(expName, fcTDelta)]
                            expLoc['cyDTime'] = myExpsCYDTimes[(expName, fcTDelta)]

                            if useRelativeDifference:
                                if statName in su.bootStrapStats:
                                    expLoc['fcTDelta'] = fcTDelta
                                    cntrlLoc['fcTDelta'] = fcTDelta

                                    X = tempdfw.loc(expLoc)
                                    Y = tempdfw.loc(cntrlLoc)

                                    # confidence intervals
                                    ciVals = su.bootStrapClusterFunc(
                                                 X, Y,
                                                 n_samples = 10000,
                                                 statNames = [statName])

                                else:
                                    ciVals = {statName: {
                                        su.cimean: expAggPlaneVals[iy, ix] - normalizingStat,
                                        su.cimin: np.NaN,
                                        su.cimax: np.NaN,
                                    }}
                            else:
                                ciVals = {statName: {
                                    su.cimean: expAggPlaneVals[iy, ix],
                                    su.cimin: np.NaN,
                                    su.cimax: np.NaN,
                                }}

                            for trait in su.ciTraits:
                                t = float(ciVals[statName][trait])
                                # automatically generate relative difference plots for positive-semi-definite statistics
                                if useRelativeDifference:
                                  # divide by cntrlLoc aggregated statName
                                  t /= normalizingStat
                                  if self.relativeErrorType == 'one hundred centered':
                                    t += 1.0
                                  t *= 100.0
                                planeVals[trait][iy, ix] = t


                    # automatically generate relative difference plots for positive-semi-definite statistics
                    if useRelativeDifference:
                        sciTicks = False
                        logScale = False

                        notused, dmin_relative, dmax_relative, centralValue, label = self.relativeErrorFunction(
                          expAggPlaneVals,
                          cntrlAggPlaneVals,
                          dmin_relative,
                          dmax_relative,
                        )
                        dmin = dmin_relative
                        dmax = dmax_relative
                        fcstatDiagLabel = statName.replace('RMS','rms').replace('Mean','mean')+': '+label

                cLabel = fcstatDiagLabel

                # perform subplot agnostic plotting (all expNames)
                bpf.plot2D(
                    fig,
                    xVals, binNumVals, planeVals[su.cimean],
                    title, xLabel, binLabel, cLabel,
                    bpf.defaultIndepConfig,
                    myBinConfigs['binConfig'],
                    sciTicks, logScale, centralValue,
                    nyplots, nxplots, nsubplots, iplot,
                    contourValsMinCI = planeVals[su.cimin],
                    contourValsMaxCI = planeVals[su.cimax],
                    dmin = dmin, dmax = dmax,
                    interiorLabels = interiorLabels)

                iplot = iplot + 1

        # save figure
        filename = ('%s%s_BinValAxisTSeries_%s-%smin_%s_%s_%s'%(
                   myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                   self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                   self.DiagSpaceName, fcDiagName, statName))

        pu.finalize_fig(fig, str(figPath/filename), figureFileType, interiorLabels, 0.35, 0.55)


class BinValAxisProfile(MultiDimBinMethodBase):
    '''
    Similar to FCandBinValAxes2D, except
      - each vertical column of raster points is plotted as a profile on
        a separate set of axes instead of in 2-D color
      - therefore this is a valid plot even for a single forecast length (omb)
      -    line: per experiment
      - subplot: column by lead time, row by DiagSpace variable
      -    file: combination of binVar, binMethod, and statistic
      - self.MAX_FC_SUBFIGS determines number of FC lead times to include
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.requestAggDFW = True

        self.subplotWidth = 1.6
        self.subplotAspect = 1.0

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options):

        bgstatDiagLabel, fcstatDiagLabel, sciTicks, logScale, centralValue = \
            self.statPlotAttributes(diagnosticGroup, statName, myLoc['diagName'])

        fcDiagName = self.fcName(diagnosticGroup)
        figPath = self.myFigPath/fcDiagName
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        if self.nFC > 1:
            nxplots = min([self.nFC, self.MAX_FC_SUBFIGS])
            nyplots = nVarsLoc
            nsubplots = nxplots * nyplots
        else:
            nsubplots = np.nanmax([nVarsLoc, 1])
            nxplots = np.nanmax([np.int(np.ceil(np.sqrt(nsubplots))), 1])
            while nsubplots%nxplots > 0 and nsubplots%nxplots / nxplots <= 0.5: nxplots += 1
            nyplots = np.int(np.ceil(np.true_divide(nsubplots, nxplots)))

        ptLoc = {}
        axisLimitsLoc = {}

        # establish a new figure
        fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, interiorLabels)
        iplot = 0

        binNumVals = myBinConfigs['num']
        binStrVals = myBinConfigs['str']
        binLabel = myBinConfigs['binLabel']

        figureData = {}
        figureData['binNumVals'] = [float(f) for f in binNumVals]
        figureData['binLabel'] = binLabel
        figureData['dataLabel'] = fcstatDiagLabel
        figureData['subplots'] = []

        #subplot loop 1
        for (varName, varLabel) in varMapLoc:
            subplotData = {}
            subplotData['varName'] = str(varName)

            ptLoc['varName'] = varName
            axisLimitsLoc['varName'] = varName

            #subplot loop 2
            for fcTDelta in self.fcTDeltas:
                tDeltaDays = float(fcTDelta.total_seconds()) / 3600.0 / 24.0
                subplotData['tDeltaDays'] = tDeltaDays

                ptLoc['fcTDelta'] = fcTDelta
                axisLimitsLoc['fcTDelta'] = fcTDelta

                # define subplot title
                title = varLabel
                if len(self.fcTDeltas) > 1:
                  title += ' @ '+str(tDeltaDays)+' days'

                # use common x-axis limits across axisLimitsLoc database locations
                if statName == 'Count':
                    dmin = 0.
                else:
                    dmin = dfwDict['agg'].min(axisLimitsLoc, statName)
                dmax = dfwDict['agg'].max(axisLimitsLoc, statName)

                subplotData['title'] = title
                #subplotData['iplot'] = iplot
                subplotData['dmin'] = dataYAMLFmtFloat(dmin)
                subplotData['dmax'] = dataYAMLFmtFloat(dmax)

                #Setting to avoid over-crowding
                if self.fcTDeltas.index(fcTDelta) > (self.MAX_FC_SUBFIGS-1): continue

                #collect aggregated statNames, varying across fcTDelta
                linesVals = []
                linesLabel = []
                linesGroup = []
                linesData = {}
                for expName in self.expNames:
                    linesData[str(expName)] = {}
                    ptLoc['expName'] = expName
                    for diagnosticName in myLoc['diagName']:
                        linesGroup.append(expName)
                        lineLabel = self.expDiagnosticLabel(
                            expName, diagnosticName, myLoc['diagName'])
                        linesLabel.append(lineLabel)

                        ptLoc['diagName'] = diagnosticName

                        lineVals = []
                        for binVal in binStrVals:
                            ptLoc['binVal'] = binVal
                            pt = dfwDict['agg'].loc1(ptLoc, statName)
                            lineVals.append(pt)
                            #if len(pt) == 1:
                              #lineVals.append(pt[0])
                            #else:
                              #lineVals.append(np.NaN)
                        linesVals.append(lineVals)
                        key = str(diagnosticName)+' '+self.DiagSpaceName+str(myLoc['binMethod'])+'_'+str(varName)
                        linesData[str(expName)][key] = []
                        for d in lineVals: linesData[str(expName)][key].append(dataYAMLFmtFloat(d))

                subplotData['linesData'] = linesData

                figureData['subplots'].append(deepcopy(subplotData))

                # perform subplot agnostic plotting (all expNames)
                options['profilefunc'](
                    fig,
                    linesVals, binNumVals,
                    linesLabel,
                    title, binLabel, fcstatDiagLabel,
                    myBinConfigs['binConfig'],
                    sciTicks, logScale, centralValue,
                    nyplots, nxplots, nsubplots, iplot,
                    dmin = dmin, dmax = dmax,
                    interiorLabels = interiorLabels)

                iplot = iplot + 1

        # save figure
        filename = '%s%s_BinValAxis_%s-%smin_%s_%s_%s'%(
                   myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                   self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                   self.DiagSpaceName, fcDiagName, statName)

        pu.finalize_fig(fig, str(figPath/filename), figureFileType, interiorLabels, 0.50, 0.55)

        # save figure data as yaml
        figureYAML = yaml.safe_dump(
          figureData,
          indent=2,
          width=2147483647,
          allow_unicode=False,
          default_flow_style=None,
        )
        with open(str(dataPath/filename)+'.yaml', 'w') as file:
           file.write(figureYAML)


class BinValAxisProfileDiffCI(MultiDimBinMethodBase):
    '''
    Similar to BinValAxisProfile, except
      shows difference between experiment(s) and control
      - control is selected using cntrlExpIndex
      - statistics are narrowed down by su.bootStrapStats
      - confidence intervals (CI) are shown at each lead time and binVal
      -    line+shaded region: per experiment
      - subplot: column by lead time, row by DiagSpace variable
      -    file: combination of binVar, binMethod, and statistic
      - self.MAX_FC_SUBFIGS determines number of FC lead times to include
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        # OPTIONAL: implement fine-grained parallelism for bootStrapping
        #self.blocking = True

        # used for percent ratio plots
        self.requestAggDFW = True

        for key in self.binVarDict:
            if 'onlyStatNames' in self.binVarDict[key]:
                self.binVarDict[key]['onlyStatNames'] += su.bootStrapStats
            else:
                self.binVarDict[key]['onlyStatNames'] = su.bootStrapStats

        self.subplotWidth = 1.6
        self.subplotAspect = 1.0

        self.requiredStatistics = ['Count']

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options):

        if self.nExp * len(myLoc['diagName']) < 2: return
        if self.cntrlExpName not in dfwDict['dfw'].levels('expName'): return

        bgstatDiagLabel_abs, fcstatDiagLabel_abs, sciTicks_abs, logScale_abs, centralValue_abs = \
            self.statPlotAttributes(diagnosticGroup, statName, myLoc['diagName'], isDifferencePlot=True)

        fcDiagName = self.fcName(diagnosticGroup)
        figPath = self.myFigPath/fcDiagName
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        if self.nFC > 1:
            nxplots = min([self.nFC, self.MAX_FC_SUBFIGS])
            nyplots = nVarsLoc
            nsubplots = nxplots * nyplots
        else:
            nsubplots = np.nanmax([nVarsLoc, 1])
            nxplots = np.nanmax([np.int(np.ceil(np.sqrt(nsubplots))), 1])
            while nsubplots%nxplots > 0 and nsubplots%nxplots / nxplots <= 0.5: nxplots += 1
            nyplots = np.int(np.ceil(np.true_divide(nsubplots, nxplots)))

        # Only bootstrap over the union of cyDTimes available
        # from both experiments at each fcTDelta
        myExpsCYDTimes = self.UNIONcntrlANDexpCYDTimes(dfwDict['dfw'])

        # establish a new figure
        fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, interiorLabels)
        iplot = 0

        binNumVals = myBinConfigs['num']
        binStrVals = myBinConfigs['str']
        binLabel = myBinConfigs['binLabel']

        useRelativeDifference = (
          statName in su.posSemiDefiniteStats and
          self.relativeErrorType != 'disable' and
          not set(myLoc['diagName']).issubset(du.absoluteOnlyDiagnostics))

        fcLoc = {}
        #subplot loop 1
        for (varName, varLabel) in varMapLoc:
            fcLoc['varName'] = varName

            #subplot loop 2
            for fcTDelta in self.fcTDeltas:
                fcLoc['fcTDelta'] = fcTDelta

                # intermediate tempdfw reduces extraction time in inner loops
                tempdfw = sdb.DFWrapper.fromLoc(dfwDict['dfw'], fcLoc)
                normdfw = sdb.DFWrapper.fromLoc(dfwDict['agg'], fcLoc)

                cntrlLoc = deepcopy(fcLoc)
                #cntrlLoc = {}
                cntrlLoc['expName'] = self.cntrlExpName
                cntrlLoc['diagName'] = myLoc['diagName'][0]
                normLoc = deepcopy(cntrlLoc)

                #Setting to avoid over-crowding
                if self.fcTDeltas.index(fcTDelta) > (self.MAX_FC_SUBFIGS-1): continue

                linesVals = defaultdict(list)
                linesLabel = []
                linesGroup = []
                for expName in self.expNames:
                    for diagnosticName in myLoc['diagName']:
                        if (expName == cntrlLoc['expName'] and
                            diagnosticName == cntrlLoc['diagName']): continue
                        cntrlLoc['cyDTime'] = myExpsCYDTimes[(expName, fcTDelta)]
                        linesGroup.append(expName)
                        linesLabel.append(self.expDiagnosticLabel(
                            expName, diagnosticName, myLoc['diagName']))

                        lineVals = defaultdict(list)

                        for binVal in binStrVals:
                            cntrlLoc['binVal'] = binVal
                            expLoc = deepcopy(cntrlLoc)
                            expLoc['diagName'] = diagnosticName
                            expLoc['expName'] = expName

                            X = tempdfw.loc(expLoc)
                            Y = tempdfw.loc(cntrlLoc)

                            # confidence intervals
                            ciVals = su.bootStrapClusterFunc(
                                         X, Y,
                                         n_samples = 10000,
                                         statNames = [statName])

                            # normalizing value for ratio
                            normLoc['binVal'] = binVal
                            normalizingStat = normdfw.loc1(normLoc, statName)

                            for trait in su.ciTraits:
                                t = float(ciVals[statName][trait])
                                # automatically generate relative difference plots for positive-semi-definite statistics
                                if useRelativeDifference:
                                  # divide by cntrlLoc aggregated statName
                                  t /= normalizingStat
                                  if self.relativeErrorType == 'one hundred centered':
                                    t += 1.0
                                  t *= 100.0
                                lineVals[trait] += [t]

                        for trait in su.ciTraits:
                            linesVals[trait].append(lineVals[trait])

                # use specific y-axis limits for each varName
                #dmin = np.nanmin(linesVals[su.cimin])
                #dmax = np.nanmax(linesVals[su.cimax])
                dmin = np.nanmin(linesVals[su.cimean])
                dmax = np.nanmax(linesVals[su.cimean])

                if not (np.isfinite(dmin) or np.isfinite(dmax)):
                  iplot = iplot + 1
                  continue

                fcstatDiagLabel = fcstatDiagLabel_abs
                sciTicks = sciTicks_abs
                logScale = logScale_abs
                centralValue = 0.

                if useRelativeDifference:
                  fcstatDiagLabel = self.relativeErrorLabeler()
                  fcstatDiagLabel = statName.replace('RMS','rms').replace('Mean','mean')+': '+fcstatDiagLabel

                  #dmin, dmax = self.relativeErrorLimiter(dmin, dmax)

                  centralValue = self.relativeErrorCenter
                  sciTicks = False
                  logScale = False

                # define subplot title
                if useRelativeDifference:
                    title = varName
                else:
                    title = varLabel

                if len(self.fcTDeltas) > 1:
                  title += ' @ '+str(float(fcTDelta.total_seconds()) / 3600.0 / 24.0)+' days'


                # perform subplot agnostic plotting (all expNames)
                options['profilefunc'](
                    fig,
                    linesVals[su.cimean], binNumVals,
                    linesLabel,
                    title, binLabel, fcstatDiagLabel,
                    myBinConfigs['binConfig'],
                    sciTicks, logScale, centralValue,
                    nyplots, nxplots, nsubplots, iplot,
                    linesValsMinCI = linesVals[su.cimin],
                    linesValsMaxCI = linesVals[su.cimax],
                    dmin = dmin, dmax = dmax,
                    lineAttribOffset = 1,
                    interiorLabels = interiorLabels)

                iplot = iplot + 1

        # save figure
        filename = ('%s%s_BinValAxis_%s-%smin_%s_%s_%s'%(
                   myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                   self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                   self.DiagSpaceName, fcDiagName, statName))

        pu.finalize_fig(fig, str(figPath/filename), figureFileType, interiorLabels, 0.50, 0.55)


class BinValAxes2D(MultiDimBinMethodBase):
    '''
    Creates raster maps with binVar binVals on both x- and y-axis
      - applicable to pairs of binned diagnostics (e.g., latitude+vertical dimension, longitude+latitude)
      - this is a valid plot even for a single cycle and/or forecast length (omb)
      - subplot: column by experiment, row by DiagSpace variable
      -    file: combination of binVar2D, (binMethod), statistic, and FC lead time (if applicable)
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.maxBinVarTier = 2
        self.requestAggDFW = True

        # default binVars
        self.binVarDict = {
            pconf.LonLat2D: {
                'plotfunc': bpf.map2D,
                'binVarTier': 1,
                'subplotWidth': 2.0,
                'subplotAspect': {
                    'default': 0.55,
                    'abi_g16': 0.9,
                    'ahi_himawari8': 0.9,
                },
                'ybuffer': 0.45,
            },
            pconf.ModelLonLat2D: {
                'plotfunc': bpf.map2D,
                'binVarTier': 1,
                'subplotWidth': 2.0,
                'subplotAspect': 0.55,
                'ybuffer': 0.45,
            },
            pconf.LatAlt2D: {
                'plotfunc': bpf.plot2D,
                'binVarTier': 1,
                'subplotWidth': 2.0,
                'subplotAspect': 0.95,
                'xbuffer': 0.6,
                'ybuffer': 0.55,
            },
            pconf.LatImpact2D: {
                'plotfunc': bpf.plot2D,
                'binVarTier': 1,
                'subplotWidth': 2.0,
                'subplotAspect': 0.95,
                'xbuffer': 0.6,
            },
            pconf.LatPrs2D: {
                'plotfunc': bpf.plot2D,
                'binVarTier': 1,
                'subplotWidth': 2.0,
                'subplotAspect': 0.8,
                'xbuffer': 0.55,
                'ybuffer': 0.70,
            },
            pconf.ModelLatLev2D: {
                'plotfunc': bpf.plot2D,
                'binVarTier': 1,
                'subplotWidth': 2.0,
                'subplotAspect': 0.8,
                'xbuffer': 0.55,
                'ybuffer': 0.70,
            },
            pconf.ModelLatPrs2D: {
                'plotfunc': bpf.plot2D,
                'binVarTier': 1,
                'subplotWidth': 2.0,
                'subplotAspect': 0.8,
                'xbuffer': 0.55,
                'ybuffer': 0.70,
            },
            pconf.ObsModel2D: {
                'plotfunc': bpf.plot2D,
                'binVarTier': 3,
                'subplotWidth': 2.0,
                'xbuffer': 0.5,
                'ybuffer': 0.55,
            },
            pconf.CldFrac2D: {
                'plotfunc': bpf.plot2D,
                'binVarTier': 1,
                'subplotWidth': 2.0,
                # uncomment to enable STD fits (i.e., Poly2DLat)
                #'twoDFittingStatistics': ['RMS'],
                'xbuffer': 0.5,
                'ybuffer': 0.55,
            },
            #pconf.CloudImpact2D: {
            #    'plotfunc': bpf.plot2D,
            #    'binVarTier': 1,
            #    'subplotWidth': 3.0,
            #},
            #pconf.OkamotoCloudImpact2D: {
            #    'plotfunc': bpf.plot2D,
            #    'binVarTier': 1,
            #    'subplotWidth': 3.0,
            #},
        }
        self.maxDiagnosticsPerAnalysis = 1

        self.requiredStatistics = ['Count']

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options):

        # TODO: REMOVE, for testing only
        #if myLoc['binMethod'] != bu.identityBinMethod: return

        subplotWidth = options.get('subplotWidth', 3.5)
        subplotAspect = deepcopy(options.get('subplotAspect', 1.0))
        xbuffer = deepcopy(options.get('xbuffer', 0.35))
        ybuffer = deepcopy(options.get('ybuffer', 0.55))

        twoDFittingStatistics = options.get('twoDFittingStatistics', [])
        if isinstance(subplotAspect, dict):
            if self.DiagSpaceName in subplotAspect:
                subplotAspect = subplotAspect[self.DiagSpaceName]
            else:
                subplotAspect = subplotAspect.get('default', 1.0)

        bgstatDiagLabel_abs, fcstatDiagLabel_abs, sciTicks_abs, logScale_abs, centralValue_abs = \
            self.statPlotAttributes(diagnosticGroup, statName)

        figPath = self.myFigPath/diagnosticGroup
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        nxplots = self.nExp
        nyplots = nVarsLoc
        nsubplots = nxplots * nyplots

        planeLoc = {}
        axisLimitsLoc = {}

        # retrieve a list of coordinates for all X/Y locations, formatted as strings
        binCoordsLevels = dfwDict['agg'].levels('binVal')

        # parse comma-separated coordinates into a list of tuples
        binCoords = list([tuple(c.split(',')) for c in binCoordsLevels])

        # unzip binCoords into independent X tuple and Y tuple
        xCoords, yCoords = zip(*binCoords)

        # determine the coordinates of the structued X/Y grid points
        xUnique = np.array(pu.uniqueMembers(xCoords))
        xVals = np.asarray(xUnique, dtype=np.float)
        xSort = np.argsort(xVals)
        xVals = xVals[xSort]
        nXVals = len(xVals)
        xValsStr = list(xUnique[xSort])

        yUnique = np.array(pu.uniqueMembers(yCoords))
        yVals = np.asarray(yUnique, dtype=np.float)
        ySort = np.argsort(yVals)
        yVals = yVals[ySort]
        nYVals = len(yVals)
        yValsStr = list(yUnique[ySort])

        # determine the X and Y indices of each location on the structured grid
        xIndex = []
        yIndex = []
        for c in binCoords:
          xIndex.append(xValsStr.index(c[0]))
          yIndex.append(yValsStr.index(c[1]))

        xUnits, xLabel = tuple(vu.varDictAll[pconf.binVars2D[myLoc['binVar']][0]])
        xVar = xLabel
        if xUnits != vu.miss_s:
            xLabel += ' ('+xUnits+')'
        for orig, sub in self.labelReplacements.items():
            xLabel = xLabel.replace(orig, sub)

        yUnits, yLabel = tuple(vu.varDictAll[pconf.binVars2D[myLoc['binVar']][1]])
        yVar = yLabel
        if yUnits != vu.miss_s:
            yLabel += ' ('+yUnits+')'
        for orig, sub in self.labelReplacements.items():
            yLabel = yLabel.replace(orig, sub)

        # special independent variable axes configs
        xVarIs = {}
        yVarIs = {}
        specialBinVars = [
          vu.obsVarPrs,
          vu.obsVarMCI,
          vu.obsVarOCI,
          vu.obsVarLogCI,
          vu.obsVarCldFracX,
          vu.obsVarCldFracY,
          vu.modVarDiagPrs,
        ]
        for var in specialBinVars:
            var_dict = vu.varDictAll.get(var,['',''])
            xVarIs[var] = (var_dict[1] == xVar)
            yVarIs[var] = (var_dict[1] == yVar)

        # TODO: make these configuration choices in one place, possibly var_utils or bpf
        xConfig = deepcopy(bpf.defaultIndepConfig)
        pCoord = xVarIs[vu.obsVarPrs] or xVarIs[vu.modVarDiagPrs]
        xConfig['invert'] = pCoord
        if pCoord: xConfig['transform'] = 'Pressure'
        if xVarIs[vu.obsVarMCI] or xVarIs[vu.obsVarOCI] or xVarIs[vu.obsVarLogCI]:
            xConfig['transform'] = 'CloudImpact'
#        if xVarIs[vu.obsVarCldFracX] or xVarIs[vu.obsVarCldFracY]:
#            xConfig['transform'] = 'logit'

        yConfig = deepcopy(bpf.defaultIndepConfig)
        pCoord = yVarIs[vu.obsVarPrs] or yVarIs[vu.modVarDiagPrs]
        yConfig['invert'] = pCoord
        if pCoord: yConfig['transform'] = 'Pressure'
        if yVarIs[vu.obsVarMCI] or yVarIs[vu.obsVarOCI] or yVarIs[vu.obsVarLogCI]:
            yConfig['transform'] = 'CloudImpact'
#        if yVarIs[vu.obsVarCldFracX] or yVarIs[vu.obsVarCldFracY]:
#            yConfig['transform'] = 'logit'

        useRelativeDifference = (
          statName in su.posSemiDefiniteStats and
          self.relativeErrorType != 'disable' and
          not set(myLoc['diagName']).issubset(du.absoluteOnlyDiagnostics))

        # Only bootstrap over the union of cyDTimes available
        # from both experiments at each fcTDelta
        myExpsCYDTimes = self.UNIONcntrlANDexpCYDTimes(dfwDict['dfw'])

        #file loop 1
        for (fcTDelta, fcTDelta_totmin) in self.fcMap:
            planeLoc['fcTDelta'] = fcTDelta

            # establish a new figure
            fig = pu.setup_fig(nxplots, nyplots, subplotWidth, subplotAspect, interiorLabels)

            iplot = 0

            #subplot loop 1
            polynomialDegrees = np.asarray([4, 6, 8, 10, 12, 14, 16])

            fitEquationConfigs = {}
            for expName in self.expNames:
              fitEquationConfigs[expName] = {}
              for degree in polynomialDegrees:
                degStr = str(degree)
                fitEquationConfigs[expName][degStr] = OrderedDict()

            for (varName, varLabel) in varMapLoc:
                planeLoc['varName'] = varName
                axisLimitsLoc['varName'] = varName

                # use common c-axis limits across axisLimitsLoc database locations
                if statName == 'Count':
                    dmin_abs = 0.
                else:
                    dmin_abs = dfwDict['agg'].min(axisLimitsLoc, statName)
                dmax_abs = dfwDict['agg'].max(axisLimitsLoc, statName)

                # letting binVal vary
                # extract control experiment
                cntrlAggLoc = deepcopy(planeLoc)
                cntrlAggLoc['expName'] = self.cntrlExpName
                cntrlAggPlaneVals = np.full((nYVals, nXVals), np.NaN)
                for ibin, binVal in enumerate(binCoordsLevels):
                    cntrlAggLoc['binVal'] = binVal
                    cntrlAggPlaneVals[yIndex[ibin], xIndex[ibin]] = \
                        dfwDict['agg'].loc1(cntrlAggLoc, statName)

                # intermediate tempdfw reduces extraction time in inner loops
                tempdfw = sdb.DFWrapper.fromLoc(dfwDict['dfw'], planeLoc)

                #subplot loop 2
                dmin_relative = np.NaN
                dmax_relative = np.NaN
                for expName in self.expNames:
                    if useRelativeDifference:
                        title = varName
                    else:
                        title = varLabel
                    title = expName+'\n'+title
                    expFileName = re.sub('\.', '', re.sub('\s+', '-', expName))

                    bgstatDiagLabel = bgstatDiagLabel_abs
                    sciTicks = sciTicks_abs
                    logScale = logScale_abs
                    centralValue = centralValue_abs
                    dmin = dmin_abs
                    dmax = dmax_abs

                    planeVals = {}
                    for trait in su.ciTraits:
                        planeVals[trait] = np.full_like(cntrlAggPlaneVals, np.NaN)

                    expAggLoc = deepcopy(cntrlAggLoc)
                    expAggPlaneVals = deepcopy(cntrlAggPlaneVals)
                    if (expName == cntrlAggLoc['expName']):
                        planeVals[su.cimean] = deepcopy(cntrlAggPlaneVals)
                    else:
                        expAggLoc['expName'] = expName

                        cntrlLoc = deepcopy(cntrlAggLoc)
                        cntrlLoc['cyDTime'] = myExpsCYDTimes[(expName, fcTDelta)]
                        expLoc = deepcopy(cntrlLoc)
                        expLoc['expName'] = expName

                        # letting binVal vary
                        # extract this experiment
                        expAggPlaneVals.fill(np.NaN)
                        for ibin, binVal in enumerate(binCoordsLevels):
                            iy = yIndex[ibin]
                            ix = xIndex[ibin]

                            expAggLoc['binVal'] = binVal
                            expAggPlaneVals[iy, ix] = \
                                dfwDict['agg'].loc1(expAggLoc, statName)

                            cntrlLoc['binVal'] = binVal
                            expLoc['binVal'] = binVal

                            # normalizing value for ratio
                            normalizingStat = cntrlAggPlaneVals[iy, ix]

                            if useRelativeDifference:
                                if statName in su.bootStrapStats:
                                    X = tempdfw.loc(expLoc)
                                    Y = tempdfw.loc(cntrlLoc)

                                    # confidence interval on difference
                                    ciVals = su.bootStrapClusterFunc(
                                                 X, Y,
                                                 n_samples = 10000,
                                                 statNames = [statName])

                                else:
                                    ciVals = {statName: {
                                        su.cimean: expAggPlaneVals[iy, ix] - normalizingStat,
                                        su.cimin: np.NaN,
                                        su.cimax: np.NaN,
                                    }}
                            else:
                                ciVals = {statName: {
                                    su.cimean: expAggPlaneVals[iy, ix],
                                    su.cimin: np.NaN,
                                    su.cimax: np.NaN,
                                }}

                            for trait in su.ciTraits:
                                t = float(ciVals[statName][trait])
                                # automatically generate relative difference plots for positive-semi-definite statistics
                                if useRelativeDifference:
                                  # divide by cntrlLoc aggregated statName
                                  t /= normalizingStat
                                  if self.relativeErrorType == 'one hundred centered':
                                    t += 1.0
                                  t *= 100.0
                                planeVals[trait][iy, ix] = t

                        # automatically generate relative difference plots for positive-semi-definite statistics
                        if useRelativeDifference:
                            # TODO: replace with same mechanisms as other DiffCI classes
                            sciTicks = False
                            logScale = False

                            notused, dmin_relative, dmax_relative, centralValue, label = self.relativeErrorFunction(
                              expAggPlaneVals,
                              cntrlAggPlaneVals,
                              dmin_relative,
                              dmax_relative,
                            )
                            if statName in su.bootStrapStats:
                              if ~any(np.isfinite(np.array([dmin_relative, dmax_relative]))):
                                dmin = np.nanmin(planeVals[su.cimean])
                                dmax = np.nanmax(planeVals[su.cimean])
                            else:
                              dmin = dmin_relative
                              dmax = dmax_relative

                            bgstatDiagLabel = statName.replace('RMS','rms').replace('Mean','mean')+': '+label

                    cLabel = bgstatDiagLabel

                    if statName in twoDFittingStatistics:
                        countVals = np.full_like(expAggPlaneVals, 0, dtype=int)
                        Y = np.empty_like(expAggPlaneVals)
                        X = np.empty_like(expAggPlaneVals)
                        for ibin, binVal in enumerate(binCoordsLevels):
                            iy = yIndex[ibin]
                            ix = xIndex[ibin]
                            X[iy, ix] = xVals[ix]
                            Y[iy, ix] = yVals[iy]
                            expAggLoc['binVal'] = binVal
                            c = dfwDict['agg'].loc1(expAggLoc, 'Count')
                            if np.isfinite(c):
                              countVals[iy, ix] = c

                        counts = countVals.flatten()
                        xf = X.flatten()
                        yf = Y.flatten()
                        zf = expAggPlaneVals.flatten()

                        ns = counts.sum()
                        xs = np.empty(ns, dtype=float)
                        ys = np.empty(ns, dtype=float)
                        zs = np.empty(ns, dtype=float)

                        NormWeight = countVals.astype(float) / ns

                        ss = 0
                        for c, x, y, z in zip(counts, xf, yf, zf):
                            ee = ss+c
                            xs[ss:ee] = x
                            ys[ss:ee] = y
                            zs[ss:ee] = z
                            ss = ee

                        self.logger.info('varName,expName=>'+varName+','+expName)
                        L2Norms = {}
                        L2Norms['values'] = []
                        weightedL2Norms = {}
                        weightedL2Norms['values'] = []

                        nxFit = len(polynomialDegrees)+1
                        nyFit = 3
                        fitFig = pu.setup_fig(nxFit, nyFit, subplotWidth, subplotAspect, interiorLabels)
                        fplot = 0

                        # plot statName
                        bpf.plot2D(
                            fitFig,
                            xVals, yVals, expAggPlaneVals,
                            'Data, '+title, xLabel, yLabel, cLabel,
                            xConfig, yConfig,
                            sciTicks_abs, logScale_abs, centralValue_abs,
                            nyFit, nxFit, nyFit*nxFit, fplot,
                            dmin = dmin_abs, dmax = dmax_abs,
                            interiorLabels = interiorLabels)

                        # plot count
                        bpf.plot2D(
                            fitFig,
                            xVals, yVals, countVals,
                            'Data, '+title, xLabel, yLabel, 'Count',
                            xConfig, yConfig,
                            True, True, None,
                            nyFit, nxFit, nyFit*nxFit, fplot+nxFit,
                            dmin = np.NaN, dmax = np.NaN,
                            interiorLabels = interiorLabels)

                        delta = np.abs(np.nanmax([(dmax - dmin) / 5., dmin, dmax]))

                        for degree in polynomialDegrees:
                            equation = poly2DEquation(degree)
                            degStr = str(degree)

                            try:
                              # handle when fit2D raises exception
                              fit = fit2D(xs, ys, zs, equation=equation)
                            except:
                              if degStr in fitEquationConfigs[expName]: del fitEquationConfigs[expName][degStr]
                              continue

                            exponents, coeffStr = fit.terms(precision=6)
                            coeffs = fit.coeffs()

                            #self.logger.info('\n'+str(exponents))
                            self.logger.info('\n '+degStr+','+varName+' coeffs: '+str(coeffStr))
                            self.logger.info('\n '+degStr+','+varName+' coeffs[0]: '+str(coeffs[0]))
                            self.logger.info('\n '+degStr+','+varName+' sum(coeffs): '+str(coeffs.sum()))

                            try:
                              fitEquationConfigs[expName][degStr][varName] = fit.terms(precision=6, returnType='dict')['terms']
                            except:
                              continue

                            Z = fit.predict(X, Y)

                            #Plot Z
                            fplot = fplot+1
                            bpf.plot2D(
                                fitFig,
                                xVals, yVals, Z,
                                'Degree '+degStr+' fit, '+title, xLabel, yLabel, cLabel,
                                xConfig, yConfig,
                                sciTicks_abs, logScale_abs, centralValue_abs,
                                nyFit, nxFit, nyFit*nxFit, fplot,
                                dmin = dmin_abs, dmax = dmax_abs,
                                interiorLabels = interiorLabels)

                            Q = Z-expAggPlaneVals

                            #Plot Q
                            L2 = np.sqrt(np.nansum(Q**2))
                            L2Norms[degStr] = '{:0.8f}'.format(L2)
                            L2Norms['values'].append(L2)
                            fitInfo = 'L2 = '+L2Norms[degStr]

                            bpf.plot2D(
                                fitFig,
                                xVals, yVals, Q,
                                fitInfo, xLabel, yLabel, '('+statName+'_fit - '+statName+')',
                                xConfig, yConfig,
                                False, False, 0,
                                nyFit, nxFit, nyFit*nxFit, fplot+nxFit,
                                dmin = -delta, dmax = delta,
                                interiorLabels = interiorLabels)

                            #Weighted Q^2
                            weightedL2 = np.sqrt(np.nansum(Q**2 * NormWeight))
                            weightedL2Norms[degStr] = '{:0.8f}'.format(weightedL2)
                            weightedL2Norms['values'].append(weightedL2)

                            fitInfo = 'Count-weighted L2 = '+weightedL2Norms[degStr]
                            bpf.plot2D(
                                fitFig,
                                xVals, yVals, Q**2 * NormWeight,
                                fitInfo, xLabel, yLabel, '('+statName+'_fit - '+statName+')^2*W_n',
                                xConfig, yConfig,
                                False, True, None,
                                nyFit, nxFit, nyFit*nxFit, fplot+2*nxFit,
                                dmin = np.NaN, dmax = delta/5.,
                                interiorLabels = interiorLabels)

                        self.logger.info('\nfit2D L2 norms: '+str(L2Norms))
                        self.logger.info('\nfit2D count-weighted L2 norms: '+str(weightedL2Norms))

                        filename = ('%s_%s_fit2D_%s_%s_%s%s_%smin'%(
                                   varName, expFileName,
                                   myLoc['binVar'], diagnosticGroup, statName,
                                   self.binMethodFile(myLoc['binMethod']),
                                   fcTDelta_totmin))

                        pu.finalize_fig(fitFig, str(figPath/filename), figureFileType, interiorLabels, 0.6)

                        # generate plots of l-curves for L2 norms vs. degree to find optimal degree
                        LFig = pu.setup_fig(2, 1, inch_width=2.0, aspect=1.0, ybuffer=True)

                        bpf.plotSeries(
                            LFig,
                            [np.asarray(L2Norms['values'])], polynomialDegrees,
                            [varName],
                            'l-curve', 'polynomial degree', 'L2 norm',
                            bpf.defaultIndepConfig,
                            False, True, None,
                            1, 2, 2, 0,
                            interiorLabels=True)

                        bpf.plotSeries(
                            LFig,
                            [np.asarray(weightedL2Norms['values'])], polynomialDegrees,
                            [varName],
                            'l-curve', 'polynomial degree', 'count-weighted L2 norm',
                            bpf.defaultIndepConfig,
                            False, True, None,
                            1, 2, 2, 1,
                            interiorLabels=True)

                        filename = ('%s_%s_L-Curves_%s_%s_%s%s_%smin'%(
                                   varName, expFileName,
                                   myLoc['binVar'], diagnosticGroup, statName,
                                   self.binMethodFile(myLoc['binMethod']),
                                   fcTDelta_totmin))

                        pu.finalize_fig(LFig, str(figPath/filename), figureFileType, True, 0.6)

                    # perform subplot agnostic plotting (all expNames)
                    if options['plotfunc'] is bpf.map2D:
                        options['plotfunc'](
                            fig,
                            xVals, yVals, planeVals[su.cimean],
                            title, cLabel,
                            sciTicks, logScale, centralValue,
                            nyplots, nxplots, nsubplots, iplot,
                            contourValsMinCI = planeVals[su.cimin],
                            contourValsMaxCI = planeVals[su.cimax],
                            dmin = dmin, dmax = dmax,
                            interiorLabels = interiorLabels)

                    else:
                        options['plotfunc'](
                            fig,
                            xVals, yVals, planeVals[su.cimean],
                            title, xLabel, yLabel, cLabel,
                            xConfig, yConfig,
                            sciTicks, logScale, centralValue,
                            nyplots, nxplots, nsubplots, iplot,
                            contourValsMinCI = planeVals[su.cimin],
                            contourValsMaxCI = planeVals[su.cimax],
                            dmin = dmin, dmax = dmax,
                            interiorLabels = interiorLabels)

                    iplot = iplot + 1

            # save each figure
            filename = ('%s%s_BinValAxes2D_%smin_%s_%s_%s'%(
                       myLoc['binVar'],
                       self.binMethodFile(myLoc['binMethod']),
                       fcTDelta_totmin, self.DiagSpaceName,
                       diagnosticGroup, statName))

            pu.finalize_fig(fig, str(figPath/filename), figureFileType, interiorLabels, xbuffer, ybuffer)

            if statName in twoDFittingStatistics:

                for expName, e1 in fitEquationConfigs.items():
                    expFileName = re.sub('\.', '', re.sub('\s+', '-', expName))
                    for degree, e2 in e1.items():
                        degStr = str(degree)
                        self.logger.info('\n '+expName+', degree: '+degStr)
                        filename = ('%s_degree=%s_fit2D_%s_%s_%s%s_%smin_%s.yaml'%(
                                   expFileName, degStr,
                                   myLoc['binVar'], diagnosticGroup, statName,
                                   self.binMethodFile(myLoc['binMethod']),
                                   fcTDelta_totmin, self.DiagSpaceName))

                        fn = str(figPath/filename)
                        anchorKeyParts = [
                          '_options for',
                          myLoc['binVar'],
                          'Terms',
                          self.DiagSpaceName,
                          'degree'+degStr+':',
                          '&'+self.DiagSpaceName+'_fit2D_'+myLoc['binVar']+self.binMethodFile(myLoc['binMethod'])+'_degree'+degStr,
                        ]
                        with open(fn, 'w') as file:
                          file.write(' '.join(anchorKeyParts)+'\n')

                        # create compactAllVarDict entry for all variables
                        firstVarName = list(e2.keys())[0]
                        compactAllVarDict = {'polynomial terms': []}
                        for term in e2[firstVarName]:
                          t = deepcopy(term)
                          del t['value']
                          compactAllVarDict['polynomial terms'].append(t)

                        # append other varName terms to 'values'
                        compactAllVarDict['polynomial coefficients'] = []
                        for iVar, (varName, termList) in enumerate(e2.items()):
                            vv = varName.replace('BT', vu.obsVarBT+'_')
                            vv, ch = vu.splitIntSuffix(vv)
                            compactAllVarDict['polynomial coefficients'].append(
                              {'name': vv, 'values': []}
                            )
                            if ch != '':
                                compactAllVarDict['polynomial coefficients'][iVar]['channel'] = int(ch)
                            for term in termList:
                                compactAllVarDict['polynomial coefficients'][iVar]['values'].append(float(term['value']))
                        compactAllVarYAML = yaml.safe_dump(
                          compactAllVarDict,
                          indent=2,
                          width=80,
                          allow_unicode=False,
                          default_flow_style=None,
                        )
                        compactAllVarYAML = indent(compactAllVarYAML, '  ')
                        with open(fn, 'a') as file:
                           file.write(compactAllVarYAML)


class BinValAxisPDF(AnalysisBase):
    '''
    Similar to BinValAxisProfile, except
      uses Count statistic to analyze a PDF across binVals
      -  x-axis: binVal
      -    line: per binMethod
      - subplot: combination of FC lead time and DiagSpace variable
      -    file: per experiment (if applicable)
    '''
    statsToPlot = ['Count']
    def __init__(self, db, analysisType, diagnosticGroupings):
      super().__init__(db, analysisType, diagnosticGroupings)
      self.binVarDict = OrderedDict()

      ## vu.obsVarNormDep
      self.binVarDict[vu.obsVarNormDep] = {}

      # all methods together
      #self.binVarDict[vu.obsVarNormDep]['all'] = {
      #  'standard gaussian': True,
      #}

      self.binVarDict[vu.obsVarNormDep][bu.identityBinMethod] = {
        'standard gaussian': True,
        'binMethodEqualsAny': [bu.identityBinMethod],
      }

      # Okamoto
      self.binVarDict[vu.obsVarNormDep][bu.OkamotoMethod] = {
        'standard gaussian': True,
        'binMethodEqualsAny': [bu.identityBinMethod, bu.OkamotoMethod],
      }

      # ClearCloudMode
      self.binVarDict[vu.obsVarNormDep][bu.ClearCloudModeMethod] = {
        'standard gaussian': True,
        'binMethodPrefix': bu.ClearCloudModeMethod+'=',
      }

      # ClearCloudMode
      self.binVarDict[vu.obsVarNormDep][bu.OkamotoMethod+','+bu.ClearCloudModeMethod] = {
        'standard gaussian': True,
        'binMethodPrefix': bu.OkamotoMethod+','+bu.ClearCloudModeMethod+'=',
      }

      # Quadrature
      #self.binVarDict[vu.obsVarNormDep][bu.QuadratureMethod] = {
        #'standard gaussian': True,
        #'binMethodEqualsAny': [bu.identityBinMethod, bu.QuadratureMethod],
      #}

# useful for diagnosing departures across ClearCloudMode, for different DA cycling scenarios
#      ## vu.obsVarDep, vu.obsVarLogDepRatio, vu.obsVarClearSkyDep
#      for departureVar in [vu.obsVarDep]:#, vu.obsVarLogDepRatio, vu.obsVarClearSkyDep]:
#        self.binVarDict[departureVar] = OrderedDict()
#
#        # ClearCloudMode
#        self.binVarDict[departureVar][bu.ClearCloudModeMethod] = {
#          'binMethodPrefix': bu.ClearCloudModeMethod+'=',
#        }

#        # ClearCloud sub-groups
#        for (v, subgroupLabel, mode, k1, k2) in pconf.ClearCloudSubgroupCases:
#          self.binVarDict[departureVar][bu.ClearCloudModeMethod+'-'+mode+subgroupLabel] = {
#            'binMethodPrefix': bu.ClearCloudModeMethod+','+subgroupLabel+'='+mode+',',
#          }

# useful for diagnosing departures across CI predictor values (e.g., Okamoto CI)
#      for departureVar in [vu.obsVarDep]:
#        # Okamoto
#        self.binVarDict[departureVar][bu.OkamotoMethod] = {
#          'binMethodPrefix': bu.OkamotoMethod+',SCI=',
#        }

      self.requestAggDFW = True

      self.subplotWidth = 1.2
      self.subplotAspect = 1.3

      self.requiredStatistics = ['Count', 'Mean', 'STD']

    def analyze_(self, workers = None):

      myLoc = {}
      myLoc['binVal'] = self.allBinNumVals2DasStr

      for diagnosticName, diagnosticConfig in self.diagnosticConfigs.items():
        if diagnosticName not in self.db.dfw.levels('diagName'): continue
        selectedStatistics = diagnosticConfig['selectedStatistics']
        availableStatistics = diagnosticConfig['availableStatistics']
        if not set(self.requiredStatistics).issubset(availableStatistics): continue

        diagBinVars = self.db.dfw.levels('binVar', {'diagName': diagnosticName})

        myLoc['diagName'] = diagnosticName

        for fullBinVar, binMethodCases in self.binVarDict.items():
          binVar = vu.varDictAll.get(fullBinVar, [None, fullBinVar])[1]
          if binVar not in diagBinVars: continue

          myLoc['binVar'] = binVar

          # reducing to dfwDict speeds extractions in innerloops
          binVarDFW = self.db.loc(myLoc)

          ## Get applicable binMethods
          allBinMethods = binVarDFW.levels('binMethod')

          for binMethodCase, options in binMethodCases.items():

            pdffunc = options.get('pdffunc', bpf.plotPDF)
            binMethodPrefix = options.get('binMethodPrefix', None)
            binMethodEqualsAny = options.get('binMethodEqualsAny', [])
            binMethodContainsAny = options.get('binMethodContainsAny', [])
            standardGaussian = options.get('standard gaussian', False)
            createAllNumericMethod = options.get('create method from all numeric', True)
            allNumericMethod = 'all'

            # collect binMethod's for this case
            caseBinMethods = []
            for binMethod in allBinMethods:
              # include all by default
              skip = False

              # exclude if binMethodPrefix not in binMethod
              if binMethodPrefix is not None:
                if not pu.prepends(binMethodPrefix, binMethod): skip = True
                #if binMethodPrefix not in binMethod: skip = True

              # only include if binMethod equals one of binMethodEqualsAny
              if len(binMethodEqualsAny) > 0: skip = True
              for isThis in binMethodEqualsAny:
                if binMethod == isThis: skip = False

              # only include if binMethod contains one of binMethodContainsAny
              if len(binMethodContainsAny) > 0: skip = True
              for containsThis in binMethodContainsAny:
                if containsThis in binMethod: skip = False

              if skip: continue

              caseBinMethods.append(binMethod)

            if len(caseBinMethods) == 0: continue

            caseLoc = deepcopy(myLoc)
            caseLoc['binMethod'] = caseBinMethods

            mydfwDict = {'dfw': sdb.DFWrapper.fromLoc(binVarDFW, caseLoc)}

            # include aggregated statistics when requested
            if self.requestAggDFW:
              mydfwDict['agg'] = sdb.DFWrapper.fromAggStats(mydfwDict['dfw'], ['cyDTime'])

            binMethodLabels0 = OrderedDict()
            for binMethod in caseBinMethods:
              if binMethod == bu.identityBinMethod:
                if len(caseBinMethods) > 1:
                  binMethodLabels0[binMethod] = r'$\sigma_o$'
                else:
                  binMethodLabels0[binMethod] = ''
              else:
                label = binMethod
                if binMethodPrefix is not None:
                  label = label.replace(binMethodPrefix,'')
                binMethodLabels0[binMethod] = label

            # default attribute offset for line color and style
            lineAttribOffset=1

            # sort the binMethod lines if all labels are floats
            if all(pu.isfloat(l) for l in binMethodLabels0.values()):
              binMethodFloats = np.asarray(list(binMethodLabels0.values()), np.float64)
              binMethodIndices = list(np.argsort(binMethodFloats))

            else:
              createAllNumericMethod = False
              binMethodIndices = list(np.argsort(list(binMethodLabels0.values())))

            #binMethodLabels0 = list(map(binMethodLabels0.__getitem__, binMethodIndices))
            caseBinMethods = list(map(caseBinMethods.__getitem__, binMethodIndices))

            # create correctly ordered binMethodLabels
            binMethodLabels = OrderedDict()

            # create a pseudo-binMethod that combines counts from all numeric real caseBinMethods
            if createAllNumericMethod:
              binMethodLabels[allNumericMethod] = allNumericMethod
              lineAttribOffset=0

            for binMethod in caseBinMethods:
              binMethodLabels[binMethod] = binMethodLabels0[binMethod]
            del binMethodLabels0

            ## Get all float/int binVals associated with binVar
            binStrVals = mydfwDict['dfw'].levels('binVal')
            binUnits = mydfwDict['dfw'].uniquevals('binUnits')[0]

            # assume all bins represent same variable/units
            binLabel = binVar
            if binUnits != vu.miss_s:
              binLabel += ' ('+binUnits+')'
            for orig, sub in self.labelReplacements.items():
              binLabel = binLabel.replace(orig, sub)

            # bin info
            binNumVals = []
            for binVal in binStrVals:
              ibin = self.allBinStrVals.index(binVal)
              binNumVals.append(self.allBinNumVals[ibin])

            # sort bins by numeric value
            indices = list(range(len(binNumVals)))
            indices.sort(key=binNumVals.__getitem__)
            binNumVals = np.array(list(map(binNumVals.__getitem__, indices)))
            binStrVals = list(map(binStrVals.__getitem__, indices))

            if len(binStrVals) < 2: continue

            self.logger.info('binVar,binMethodCase=>'+binVar+','+binMethodCase)

            fcDiagName = self.fcName(diagnosticName)
            figPath = self.myFigPath/fcDiagName
            figPath.mkdir(parents=True, exist_ok=True)
            dataPath = figPath/'data'
            dataPath.mkdir(parents=True, exist_ok=True)

            if self.nFC > 1:
              nxplots = min([self.nFC, self.MAX_FC_SUBFIGS])
              nyplots = self.nVars
              nsubplots = nxplots * nyplots
            else:
              nsubplots = self.nVars
              nxplots = np.int(np.ceil(np.sqrt(nsubplots)))
              while nsubplots%nxplots > 0 and nsubplots%nxplots / nxplots <= 0.5: nxplots += 1
              nyplots = np.int(np.ceil(np.true_divide(nsubplots, nxplots)))

            subplotLoc = deepcopy(caseLoc)

            #file loop 1
            for expName in self.expNames:
              subplotLoc['expName'] = expName

              # establish a new figure
              fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, interiorLabels)
              fig_normalized = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, interiorLabels)

              iplot = 0

              #subplot loop 1
              for (varName, varLabel) in self.varMap:
                subplotLoc['varName'] = varName

                #subplot loop 2
                for fcTDelta in self.fcTDeltas:
                  subplotLoc['fcTDelta'] = fcTDelta

                  subplotdfw = sdb.DFWrapper.fromLoc(mydfwDict['agg'], subplotLoc)

                  # aggregate across binVal for each binMethod
                  eachBinMethod = sdb.DFWrapper.fromAggStats(subplotdfw, ['binVal'])

                  #Setting to avoid over-crowding
                  if self.fcTDeltas.index(fcTDelta) > (self.MAX_FC_SUBFIGS-1): continue

                  #collect aggregated statNames, varying across fcTDelta
                  statsVals = {}
                  for statName in self.statsToPlot:
                    statsVals[statName] = {}
                    #statsVals[statName] = OrderedDict()
                    for binMethod in caseBinMethods:
                      statsVals[statName][binMethod] = np.empty(len(binStrVals), dtype=su.statDtypes[statName])

                  ptLoc = deepcopy(subplotLoc)
                  binMethodAggStats = {}
                  for binMethod in caseBinMethods:
                    ptLoc['binMethod'] = binMethod

                    for ii, binVal in enumerate(binStrVals):
                      ptLoc['binVal'] = binVal
                      for statName in statsVals.keys():
                        statsVals[statName][binMethod][ii] = subplotdfw.loc1(ptLoc, statName)

                    # extract gross statistics for each binMethod
                    binMethodAggStats[binMethod] = {}
                    for statName in su.aggregatableFileStats:
                      binMethodAggStats[binMethod][statName] = eachBinMethod.loc1(
                        {'binMethod': binMethod}, statName)

                  if createAllNumericMethod:
                    # aggregate across binMethod for each binVal
                    eachBinVal = sdb.DFWrapper.fromAggStats(subplotdfw, ['binMethod'])

                    # extract self.statsToPlot at each binVal for allNumericMethod
                    for statName in self.statsToPlot:
                      statsVals[statName][allNumericMethod] = np.empty(len(binStrVals), dtype=su.statDtypes[statName])
                      for ii, binVal in enumerate(binStrVals):
                        statsVals[statName][allNumericMethod][ii] = eachBinVal.loc1(
                          {'binVal': binVal}, statName)

                    # extract gross statistics for allNumericMethod
                    allNumeric = sdb.DFWrapper.fromAggStats(subplotdfw, ['binMethod', 'binVal'])

                    binMethodAggStats[allNumericMethod] = {}
                    for statName in su.aggregatableFileStats:
                      s = allNumeric.var(statName).to_numpy()
                      if isinstance(s, Iterable):
                        if len(s) == 1:
                          binMethodAggStats[allNumericMethod][statName] = s[0]
                      else:
                        binMethodAggStats[allNumericMethod][statName] = s

                  # extract counts for easy plotting
                  countsVals = []
                  for binMethod in binMethodLabels.keys():
                    countsVals.append(statsVals['Count'][binMethod])

                  # append Mean, STD, Skew, ExcessKurtosis to labels
                  binMethodLabelsValues = list(binMethodLabels.values())

                  binMethodLabelsWithMetrics = []
                  label = []
                  for ii, (binMethod, binMethodLabel) in enumerate(list(binMethodLabels.items())):
                    if binMethodLabel != '':
                      label += [binMethodLabel]
                    label += [r'$\mu=${:.1f}'.format(binMethodAggStats[binMethod]['Mean'])]
                    label += [r'$\sigma=${:.2f}'.format(binMethodAggStats[binMethod]['STD'])]
                    #label += [r'$\nu=${:.2f}'.format(binMethodAggStats[binMethod]['Skew'])]
                    #label += [r'$\kappa=${:.2f}'.format(binMethodAggStats[binMethod]['ExcessKurtosis'])]

                    binMethodLabelsWithMetrics.append(','.join(label))

                  # define subplot title
                  title = varLabel
                  if len(self.fcTDeltas) > 1:
                    title += ' @ '+str(float(fcTDelta.total_seconds()) / 3600.0 / 24.0)+' days'

                  # perform subplot agnostic plotting (all expNames)
                  # raw counts
                  pdffunc(
                    fig,
                    countsVals, binNumVals,
                    binMethodLabelsValues,
                    title,
                    binLabel,
                    nyplots, nxplots, nsubplots, iplot,
                    lineAttribOffset = lineAttribOffset,
                    interiorLabels = interiorLabels,
                    normalized = False,
                    standardGaussian = False,
                  )

                  # counts normalized by bin-width and total count
                  pdffunc(
                    fig_normalized,
                    countsVals, binNumVals,
                    binMethodLabelsWithMetrics,
                    title,
                    binLabel,
                    nyplots, nxplots, nsubplots, iplot,
                    lineAttribOffset = lineAttribOffset,
                    interiorLabels = interiorLabels,
                    normalized = True,
                    standardGaussian = standardGaussian,
                  )

                  iplot = iplot + 1

                # end fcTDelta loop

              # end varMap loop

              # save each figure
              filename = ('%s_%s_BinValAxis_%s-%smin_%s_%s_%s'%(
                         binVar, binMethodCase, self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                         self.DiagSpaceName, fcDiagName, expName))

              pu.finalize_fig(fig, str(figPath/filename), figureFileType, interiorLabels, 0.6)

              filename = ('%s_%s-normalized_BinValAxis_%s-%smin_%s_%s_%s'%(
                         binVar, binMethodCase, self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                         self.DiagSpaceName, fcDiagName, expName))

              pu.finalize_fig(fig_normalized, str(figPath/filename), figureFileType, interiorLabels, 0.6)

            # end expName loop


class BinValAxisPDFMultiExp(AnalysisBase):
    '''
    Similar to BinValAxisProfile, except
      uses Count statistic to analyze a PDF across binVals
      -  x-axis: binVal
      -    line: per binMethod and experiment
      - subplot: combination of FC lead time and DiagSpace variable
    '''
    statsToPlot = ['Count']
    def __init__(self, db, analysisType, diagnosticGroupings):
      super().__init__(db, analysisType, diagnosticGroupings)
      self.binVarDict = OrderedDict()

      ## vu.obsVarNormDep
      self.binVarDict[vu.obsVarNormDep] = {}

      self.binVarDict[vu.obsVarNormDep][bu.identityBinMethod] = {
        'standard gaussian': True,
        'binMethodEqualsAny': [bu.identityBinMethod],
      }

      self.requestAggDFW = True

      self.subplotWidth = 1.2
      self.subplotAspect = 1.3

      self.requiredStatistics = ['Count', 'Mean', 'STD']

    def analyze_(self, workers = None):

      myLoc = {}
      myLoc['binVal'] = self.allBinNumVals2DasStr

      for diagnosticName, diagnosticConfig in self.diagnosticConfigs.items():
        if diagnosticName not in self.db.dfw.levels('diagName'): continue
        selectedStatistics = diagnosticConfig['selectedStatistics']
        availableStatistics = diagnosticConfig['availableStatistics']
        if not set(self.requiredStatistics).issubset(availableStatistics): continue

        diagBinVars = self.db.dfw.levels('binVar', {'diagName': diagnosticName})

        myLoc['diagName'] = diagnosticName

        for fullBinVar, binMethodCases in self.binVarDict.items():
          binVar = vu.varDictAll.get(fullBinVar, [None, fullBinVar])[1]
          if binVar not in diagBinVars: continue

          myLoc['binVar'] = binVar

          # reducing to dfwDict speeds extractions in innerloops
          binVarDFW = self.db.loc(myLoc)

          ## Get applicable binMethods
          allBinMethods = binVarDFW.levels('binMethod')

          for binMethodCase, options in binMethodCases.items():

            pdffunc = options.get('pdffunc', bpf.plotPDF)
            binMethodPrefix = options.get('binMethodPrefix', None)
            binMethodEqualsAny = options.get('binMethodEqualsAny', [])
            binMethodContainsAny = options.get('binMethodContainsAny', [])
            standardGaussian = options.get('standard gaussian', False)
            createAllNumericMethod = options.get('create method from all numeric', True)
            allNumericMethod = 'all'

            # collect binMethod's for this case
            caseBinMethods = []
            for binMethod in allBinMethods:
              # include all by default
              skip = False

              # exclude if binMethodPrefix not in binMethod
              if binMethodPrefix is not None:
                if not pu.prepends(binMethodPrefix, binMethod): skip = True
                #if binMethodPrefix not in binMethod: skip = True

              # only include if binMethod equals one of binMethodEqualsAny
              if len(binMethodEqualsAny) > 0: skip = True
              for isThis in binMethodEqualsAny:
                if binMethod == isThis: skip = False

              # only include if binMethod contains one of binMethodContainsAny
              if len(binMethodContainsAny) > 0: skip = True
              for containsThis in binMethodContainsAny:
                if containsThis in binMethod: skip = False

              if skip: continue

              caseBinMethods.append(binMethod)

            if len(caseBinMethods) == 0: continue

            caseLoc = deepcopy(myLoc)
            caseLoc['binMethod'] = caseBinMethods

            mydfwDict = {'dfw': sdb.DFWrapper.fromLoc(binVarDFW, caseLoc)}

            # include aggregated statistics when requested
            if self.requestAggDFW:
              mydfwDict['agg'] = sdb.DFWrapper.fromAggStats(mydfwDict['dfw'], ['cyDTime'])

            binMethodLabels0 = OrderedDict()
            for binMethod in caseBinMethods:
              if binMethod == bu.identityBinMethod:
                if len(caseBinMethods) > 1:
                  binMethodLabels0[binMethod] = r'$\sigma_o$'
                else:
                  binMethodLabels0[binMethod] = ''
              else:
                label = binMethod
                if binMethodPrefix is not None:
                  label = label.replace(binMethodPrefix,'')
                binMethodLabels0[binMethod] = label

            # default attribute offset for line color and style
            lineAttribOffset=1

            # sort the binMethod lines if all labels are floats
            if all(pu.isfloat(l) for l in binMethodLabels0.values()):
              binMethodFloats = np.asarray(list(binMethodLabels0.values()), np.float64)
              binMethodIndices = list(np.argsort(binMethodFloats))

            else:
              createAllNumericMethod = False
              binMethodIndices = list(np.argsort(list(binMethodLabels0.values())))

            #binMethodLabels0 = list(map(binMethodLabels0.__getitem__, binMethodIndices))
            caseBinMethods = list(map(caseBinMethods.__getitem__, binMethodIndices))

            # create correctly ordered binMethodLabels
            binMethodLabels = OrderedDict()

            # create a pseudo-binMethod that combines counts from all numeric real caseBinMethods
            if createAllNumericMethod:
              binMethodLabels[allNumericMethod] = allNumericMethod
              lineAttribOffset=0

            for binMethod in caseBinMethods:
              binMethodLabels[binMethod] = binMethodLabels0[binMethod]
            del binMethodLabels0

            ## Get all float/int binVals associated with binVar
            binStrVals = mydfwDict['dfw'].levels('binVal')
            binUnits = mydfwDict['dfw'].uniquevals('binUnits')[0]

            # assume all bins represent same variable/units
            binLabel = binVar
            if binUnits != vu.miss_s:
              binLabel += ' ('+binUnits+')'
            for orig, sub in self.labelReplacements.items():
              binLabel = binLabel.replace(orig, sub)

            # bin info
            binNumVals = []
            for binVal in binStrVals:
              ibin = self.allBinStrVals.index(binVal)
              binNumVals.append(self.allBinNumVals[ibin])

            # sort bins by numeric value
            indices = list(range(len(binNumVals)))
            indices.sort(key=binNumVals.__getitem__)
            binNumVals = np.array(list(map(binNumVals.__getitem__, indices)))
            binStrVals = list(map(binStrVals.__getitem__, indices))

            if len(binStrVals) < 2: continue

            self.logger.info('binVar,binMethodCase=>'+binVar+','+binMethodCase)

            fcDiagName = self.fcName(diagnosticName)
            figPath = self.myFigPath/fcDiagName
            figPath.mkdir(parents=True, exist_ok=True)
            dataPath = figPath/'data'
            dataPath.mkdir(parents=True, exist_ok=True)

            if self.nFC > 1:
              nxplots = min([self.nFC, self.MAX_FC_SUBFIGS])
              nyplots = self.nVars
              nsubplots = nxplots * nyplots
            else:
              nsubplots = self.nVars
              nxplots = np.int(np.ceil(np.sqrt(nsubplots)))
              while nsubplots%nxplots > 0 and nsubplots%nxplots / nxplots <= 0.5: nxplots += 1
              nyplots = np.int(np.ceil(np.true_divide(nsubplots, nxplots)))

            subplotLoc = deepcopy(caseLoc)

            # establish a new figure
            fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, interiorLabels)
            fig_normalized = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, interiorLabels)

            iplot = 0

            #subplot loop 1
            for (varName, varLabel) in self.varMap:
              subplotLoc['varName'] = varName

              #subplot loop 2
              for fcTDelta in self.fcTDeltas:
                subplotLoc['fcTDelta'] = fcTDelta

                expCounts = []
                expBinMethodLabelsValues = []
                expBinMethodLabelsWithMetrics = []

                #file loop 1
                for expName in self.expNames:
                  subplotLoc['expName'] = expName

                  subplotdfw = sdb.DFWrapper.fromLoc(mydfwDict['agg'], subplotLoc)

                  # aggregate across binVal for each binMethod
                  eachBinMethod = sdb.DFWrapper.fromAggStats(subplotdfw, ['binVal'])

                  #Setting to avoid over-crowding
                  if self.fcTDeltas.index(fcTDelta) > (self.MAX_FC_SUBFIGS-1): continue

                  #collect aggregated statNames, varying across fcTDelta
                  statsVals = {}
                  for statName in self.statsToPlot:
                    statsVals[statName] = {}
                    #statsVals[statName] = OrderedDict()
                    for binMethod in caseBinMethods:
                      statsVals[statName][binMethod] = np.empty(len(binStrVals), dtype=su.statDtypes[statName])

                  ptLoc = deepcopy(subplotLoc)
                  binMethodAggStats = {}
                  for binMethod in caseBinMethods:
                    ptLoc['binMethod'] = binMethod

                    for ii, binVal in enumerate(binStrVals):
                      ptLoc['binVal'] = binVal
                      for statName in statsVals.keys():
                        statsVals[statName][binMethod][ii] = subplotdfw.loc1(ptLoc, statName)

                    # extract gross statistics for each binMethod
                    binMethodAggStats[binMethod] = {}
                    for statName in su.aggregatableFileStats:
                      binMethodAggStats[binMethod][statName] = eachBinMethod.loc1(
                        {'binMethod': binMethod}, statName)

                  if createAllNumericMethod:
                    # aggregate across binMethod for each binVal
                    eachBinVal = sdb.DFWrapper.fromAggStats(subplotdfw, ['binMethod'])

                    # extract self.statsToPlot at each binVal for allNumericMethod
                    for statName in self.statsToPlot:
                      statsVals[statName][allNumericMethod] = np.empty(len(binStrVals), dtype=su.statDtypes[statName])
                      for ii, binVal in enumerate(binStrVals):
                        statsVals[statName][allNumericMethod][ii] = eachBinVal.loc1(
                          {'binVal': binVal}, statName)

                    # extract gross statistics for allNumericMethod
                    allNumeric = sdb.DFWrapper.fromAggStats(subplotdfw, ['binMethod', 'binVal'])

                    binMethodAggStats[allNumericMethod] = {}
                    for statName in su.aggregatableFileStats:
                      s = allNumeric.var(statName).to_numpy()
                      if isinstance(s, Iterable):
                        if len(s) == 1:
                          binMethodAggStats[allNumericMethod][statName] = s[0]
                      else:
                        binMethodAggStats[allNumericMethod][statName] = s

                  # extract counts for easy plotting
                  countsVals = []
                  for binMethod in binMethodLabels.keys():
                    countsVals.append(statsVals['Count'][binMethod])

                  # append Mean, STD, Skew, ExcessKurtosis to labels
                  binMethodLabelsValues = list(binMethodLabels.values())

                  binMethodLabelsWithMetrics = []
                  label = []
                  for ii, (binMethod, binMethodLabel) in enumerate(list(binMethodLabels.items())):
                    if binMethodLabel != '':
                      label += [binMethodLabel]
                    label += [r'$\mu=${:.1f}'.format(binMethodAggStats[binMethod]['Mean'])]
                    label += [r'$\sigma=${:.2f}'.format(binMethodAggStats[binMethod]['STD'])]
                    #label += [r'$\nu=${:.2f}'.format(binMethodAggStats[binMethod]['Skew'])]
                    #label += [r'$\kappa=${:.2f}'.format(binMethodAggStats[binMethod]['ExcessKurtosis'])]

                    binMethodLabelsWithMetrics.append(','.join(label))

                  # define subplot title
                  title = varLabel
                  if len(self.fcTDeltas) > 1:
                    title += ' @ '+str(float(fcTDelta.total_seconds()) / 3600.0 / 24.0)+' days'

                  expCounts += countsVals
                  for b in binMethodLabelsValues:
                    if b == '':
                      expBinMethodLabelsValues.append(expName)
                    else:
                      expBinMethodLabelsValues.append(expName+'--'+b)

                  for b in binMethodLabelsWithMetrics:
                    if b == '':
                      expBinMethodLabelsWithMetrics.append(expName)
                    else:
                      expBinMethodLabelsWithMetrics.append(expName+'--'+b)

                # end expName loop

                # perform subplot agnostic plotting (all expNames)
                # raw counts
                pdffunc(
                  fig,
                  expCounts, binNumVals,
                  expBinMethodLabelsValues,
                  title,
                  binLabel,
                  nyplots, nxplots, nsubplots, iplot,
                  lineAttribOffset = lineAttribOffset,
                  interiorLabels = interiorLabels,
                  normalized = False,
                  standardGaussian = False,
                )

                # counts normalized by bin-width and total count
                pdffunc(
                  fig_normalized,
                  expCounts, binNumVals,
                  expBinMethodLabelsWithMetrics,
                  title,
                  binLabel,
                  nyplots, nxplots, nsubplots, iplot,
                  lineAttribOffset = lineAttribOffset,
                  interiorLabels = interiorLabels,
                  normalized = True,
                  standardGaussian = standardGaussian,
                )

                iplot = iplot + 1

              # end fcTDelta loop

            # end varMap loop

            # save each figure
            filename = ('%s_%s_BinValAxis_%s-%smin_%s_%s_allExp'%(
                       binVar, binMethodCase, self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                       self.DiagSpaceName, fcDiagName))

            pu.finalize_fig(fig, str(figPath/filename), figureFileType, interiorLabels, 0.6)

            filename = ('%s_%s-normalized_BinValAxis_%s-%smin_%s_%s_allExp'%(
                       binVar, binMethodCase, self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                       self.DiagSpaceName, fcDiagName))

            pu.finalize_fig(fig_normalized, str(figPath/filename), figureFileType, interiorLabels, 0.6)



# TODO: generalize as a sub-class of MultiDimBinMethodBase
class BinValAxisStatsComposite(AnalysisBase):
    '''
    Similar to BinValAxisProfile, except
      all statistics (Count, Mean, RMS, STD) are placed on the same axis
      -  x-axis: binVal
      -    line: per statistic
      - subplot: per DiagSpace variable
      -    file: combination of FC lead time, experiment, and binMethod (if applicable)
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.maxBinVarTier = 2

        self.binVarDict = {
            vu.obsVarAlt: {
              'statsfunc': bpf.plotCompositeProfile,
              'binFilter': {
                # maximum altitude to show on all figures
                'maxvalue': 30000.,
              },
            },
            vu.obsVarImpact: {
              'statsfunc': bpf.plotCompositeProfile,
              'binFilter': {
                # maximum altitude to show on all figures
                'maxvalue': 30000.,
              },
            },
            vu.obsVarACI: {'statsfunc': bpf.plotCompositeSeries, 'binVarTier': 3},
            vu.obsVarCI: {'statsfunc': bpf.plotfitRampComposite},
            vu.obsVarCldFracX: {'statsfunc': bpf.plotCompositeSeries, 'binVarTier': 2},
            vu.obsVarCldFracY: {'statsfunc': bpf.plotCompositeSeries, 'binVarTier': 1},
            vu.obsVarLat: {'statsfunc': bpf.plotCompositeProfile},
            #vu.obsVarLogCI: {'statsfunc': bpf.plotfitRampComposite},
            vu.obsVarPrs: {'statsfunc': bpf.plotCompositeProfile},
            vu.modVarDiagPrs: {'profilefunc': bpf.plotCompositeProfile},
            # vu.modVarLat is redundant with vu.obsVarLat (both have varShort=="lat")
            #vu.modVarLat: {'statsfunc': bpf.plotCompositeProfile, 'binVarTier': 1},
            vu.modVarLev: {
              'statsfunc': bpf.plotCompositeProfile,
              'binFilter': {
                # maximum model level to show on all figures
                #'maxvalue': 40,
              },
            },
            vu.obsVarGlint: {'statsfunc': bpf.plotCompositeSeries, 'binVarTier': 3},
            vu.obsVarLandFrac: {'statsfunc': bpf.plotCompositeSeries, 'binVarTier': 3},
            vu.obsVarLT: {'statsfunc': bpf.plotCompositeSeries, 'binVarTier': 3},
            vu.obsVarSenZen: {'statsfunc': bpf.plotCompositeSeries, 'binVarTier': 3},
        }

        self.requestAggDFW = True

        # Force serial processing so that console output is contiguous
        # TODO(JJG): output to an ascii file and remove this line
        self.blocking = True

        self.subplotWidth = 1.9
        self.subplotAspect = 0.9

        self.requiredStatistics = ['Count', 'Mean', 'RMS', 'STD']

    def analyze_(self, workers = None):
        for diagnosticName, diagnosticConfig in self.diagnosticConfigs.items():
            if diagnosticName not in self.db.dfw.levels('diagName'): continue
            selectedStatistics = self.requiredStatistics
            availableStatistics = diagnosticConfig['availableStatistics']
            if not set(self.requiredStatistics).issubset(availableStatistics): continue

            diagBinVars = self.db.dfw.levels('binVar', {'diagName': diagnosticName})

            for fullBinVar, options in self.binVarDict.items():
                binVar = vu.varDictAll.get(fullBinVar, [None, fullBinVar])[1]
                if (binVar not in diagBinVars): continue
                if options.get('binVarTier', 1) > self.maxBinVarTier: continue

                binFilter = options.get('binFilter', None)

                myLoc = {}
                myLoc['diagName'] = diagnosticName
                myLoc['binVar'] = binVar
                myLoc['binVal'] = self.allBinNumVals2DasStr

                # reducing to mydfwDict speeds extractions in innerloops
                mydfwDict = {'dfw': self.db.loc(myLoc)}

                ## Get all float/int binVals associated with binVar
                binMethods = mydfwDict['dfw'].levels('binMethod')
                binUnits = mydfwDict['dfw'].uniquevals('binUnits')[0]

                # assume all bins represent same variable/units
                binLabel = binVar
                if binUnits != vu.miss_s:
                    binLabel += ' ('+binUnits+')'
                for orig, sub in self.labelReplacements.items():
                    binLabel = binLabel.replace(orig, sub)

                fcDiagName = self.fcName(diagnosticName)
                figPath = self.myFigPath/fcDiagName
                figPath.mkdir(parents=True, exist_ok=True)
                dataPath = figPath/'data'
                dataPath.mkdir(parents=True, exist_ok=True)

                # only analyze variables that have non-zero Count when sliced by myLoc
                nVarsLoc = 0
                varMapLoc = []
                for (varName, varLabel) in self.varMap:
                  if 'Count' in selectedStatistics:
                      countDF = mydfwDict['dfw'].loc({'varName': varName}, 'Count')
                      if countDF.shape[0] > 0:
                        if np.nansum(countDF.to_numpy()) > 0:
                          nVarsLoc += 1
                          varMapLoc.append((varName, varLabel))
                  else:
                      statDF = mydfwDict['dfw'].loc({'varName': varName}, list(selectedStatistics)[0])
                      if statDF.shape[0] > 0:
                        if np.isfinite(statDF.to_numpy()).sum() > 0:
                          nVarsLoc += 1
                          varMapLoc.append((varName, varLabel))

                nsubplots = nVarsLoc
                nxplots = np.int(np.ceil(np.sqrt(nsubplots)))
                while nsubplots%nxplots > 0 and nsubplots%nxplots / nxplots <= 0.5: nxplots += 1
                nyplots = np.int(np.ceil(np.true_divide(nsubplots, nxplots)))

                ptLoc = {}

                #file loop 1
                for binMethod in binMethods:
                    ptLoc['binMethod'] = binMethod

                    methodDFW = sdb.DFWrapper.fromLoc(mydfwDict['dfw'], {'binMethod': binMethod})

                    # parse binVals, which may be different for each binMethod
                    binStrVals = methodDFW.levels('binVal')

                    # bin info
                    binNumVals = []
                    for binVal in binStrVals:
                        ibin = self.allBinStrVals.index(binVal)
                        binNumVals.append(self.allBinNumVals[ibin])

                    # sort bins by numeric value
                    indices = list(range(len(binNumVals)))
                    indices.sort(key=binNumVals.__getitem__)
                    binNumVals = list(map(binNumVals.__getitem__, indices))
                    binStrVals = list(map(binStrVals.__getitem__, indices))

                    if binFilter is not None:
                      remove = np.full_like(binNumVals, False, bool)
                      minvalue = binFilter.get('minvalue', None)
                      if minvalue is not None:
                        less = bu.lessBound(np.asarray(binNumVals), minvalue)
                        remove[less] = True
                      maxvalue = binFilter.get('maxvalue', None)
                      if maxvalue is not None:
                        great = bu.greatBound(np.asarray(binNumVals), maxvalue)
                        remove[great] = True

                      binStrVals = list(np.asarray(binStrVals)[~remove])
                      binNumVals = list(np.asarray(binNumVals)[~remove])

                    nBins = len(binStrVals)
                    if nBins < 2: continue

                    self.logger.info('binVar=>'+binVar+', binMethod=>'+binMethod)

                    methodDFWAgg = sdb.DFWrapper.fromAggStats(methodDFW, ['cyDTime'])

                    #file loop 2
                    for expName in self.expNames:
                        ptLoc['expName'] = expName

                        #file loop 3
                        for (fcTDelta, fcTDelta_totmin) in self.fcMap:
                            ptLoc['fcTDelta'] = fcTDelta

                            # establish a new figure
                            fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, interiorLabels)
                            iplot = 0

                            ERRParams = {}
                            ERRParams[self.DiagSpaceName] = {}

                            #subplot loop 1
                            for (varName, varLabel) in varMapLoc:
                                ptLoc['varName'] = varName

                                #collect aggregated statNames, varying across fcTDelta
                                countsVals = np.full(nBins, 0)
                                meansVals  = np.full(nBins, np.NaN)
                                rmssVals   = np.full(nBins, np.NaN)
                                stdsVals   = np.full(nBins, np.NaN)

                                for ibin, binVal in enumerate(binStrVals):
                                    ptLoc['binVal'] = binVal
                                    c = methodDFWAgg.loc1(ptLoc,'Count')
                                    if np.isfinite(c):
                                      countsVals[ibin] = c
                                    meansVals[ibin] = methodDFWAgg.loc1(ptLoc,'Mean')
                                    rmssVals[ibin] = methodDFWAgg.loc1(ptLoc,'RMS')
                                    stdsVals[ibin] = methodDFWAgg.loc1(ptLoc,'STD')

                                # define subplot title
                                title = varLabel

                                # perform subplot agnostic plotting (all expNames)
                                FitParams = options['statsfunc'](
                                    fig,
                                    binNumVals,
                                    countsVals,
                                    meansVals,
                                    rmssVals,
                                    stdsVals,
                                    title,
                                    'STATS('+fcDiagName+')',
                                    binLabel,
                                    nyplots, nxplots, nsubplots, iplot,
                                    interiorLabels = interiorLabels)

                                paramKey = self.chlist[iplot]
                                if paramKey < 0: paramKey = varName
                                if FitParams is not None:
                                    ERRParams[self.DiagSpaceName][(paramKey, binMethod)] = FitParams

                                iplot = iplot + 1

                            YAMLParams = {}
                            self.logger.info('#For binning_params('+expName+'):')
                            for key in sorted(ERRParams[self.DiagSpaceName]):
                                self.logger.info(binVar+"ErrParams['"+self.DiagSpaceName+"']["+str(key)+"]   = "+
                                       str(ERRParams[self.DiagSpaceName][key]['bin_utils']))
                                for param, val in ERRParams[self.DiagSpaceName][key]['YAML'].items():
                                    if param not in YAMLParams: YAMLParams[param] = []
                                    YAMLParams[param] += val
                            self.logger.info('#For UFO YAML config('+expName+'):')
                            for param, val in YAMLParams.items():
                                self.logger.info('#  '+param+': '+str(val))

                            # save each figure
                            filename = ('%s%s_BinValAxis_%smin_%s_%s_%s'%(
                                       binVar, self.binMethodFile(binMethod), fcTDelta_totmin,
                                       self.DiagSpaceName, fcDiagName, expName))

                            pu.finalize_fig(fig, str(figPath/filename), figureFileType, interiorLabels, 0.6)

                    # end expName loop

                # end binMethod loop

            # end fullBinVar loop


#===========================
# Calculate gross statistics
#===========================
class GrossValues(AnalysisBase):
    '''
    Calculate gross statistics for specified category binMethods at first forecast length
      NOTE: currently only calculates statistics at self.fcTDeltas[0]
            adjust minimum forecast length in order to calculate
            for non-zero forecast lengths, assuming those lengths
            are present in db
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.requestAggDFW = True

        # Force serial processing so that console output is contiguous
        # TODO(JJG): output to an ascii file and remove this line
        self.blocking = True

        self.binVarDict = {
            (vu.obsVarQC, bu.goodQCMethod): {},
            (vu.obsVarCldFracY, bu.cloudbandsMethod): {},
        }

    def analyze_(self, workers = None):
        for diagnosticName, diagnosticConfig in self.diagnosticConfigs.items():
            if diagnosticName not in self.db.dfw.levels('diagName'): continue
            selectedStatistics = diagnosticConfig['selectedStatistics']
            availableStatistics = diagnosticConfig['availableStatistics']
            if not set(self.requiredStatistics).issubset(availableStatistics): continue

            diagLoc = {'diagName': diagnosticName}
            diagBinVars = self.db.dfw.levels('binVar', diagLoc)
            diagBinMethods = self.db.dfw.levels('binMethod', diagLoc)

            for (fullBinVar, binMethod), options in self.binVarDict.items():
                binVar = vu.varDictAll.get(fullBinVar, [None, fullBinVar])[1]
                if (binVar not in diagBinVars or
                    binMethod not in diagBinMethods): continue

                # narrow mydfwDict by binVar and binMethod to reduce run-time and memory
                myLoc = {}
                myLoc['binVar'] = binVar
                myLoc['binMethod'] = binMethod

                # reducing to mydfwDict speeds extractions in innerloops
                mydfwDict = {'dfw': self.db.loc(myLoc)}

                if self.requestAggDFW:
                    mydfwDict['agg'] = sdb.DFWrapper.fromAggStats(mydfwDict['dfw'], ['cyDTime'])
                    sdb.createORreplaceDerivedDiagnostics(mydfwDict['agg'], {diagnosticName: diagnosticConfig})

                # further narrow mydfwDict by diagName
                # NOTE: derived diagnostics may require multiple diagName values;
                # can only narrow by diagName after aggregation
                myLoc['diagName'] = diagnosticName
                for key in mydfwDict.keys():
                    mydfwDict[key] = sdb.DFWrapper.fromLoc(mydfwDict[key], myLoc)

                print(' Calculate gross statistics: binVar=>'+binVar+', binMethod=>'+binMethod)

                binValsMap = categoryBinValsAttributes(
                    mydfwDict['dfw'], fullBinVar, binMethod, options)

                print(' at FC length ', self.fcTDeltas[0])
                # Calculate gross statistics for this binVal
                statsLoc = {}
                statsLoc['fcTDelta'] = self.fcTDeltas[0]
                for binVal, binTitle in binValsMap:
                    statsLoc['binVal'] = binVal
                    GrossValues = {}
                    for varName in self.varNames:
                        statsLoc['varName'] = varName
                        for expName in self.expNames:
                            statsLoc['expName'] = expName
                            statsDFW = sdb.DFWrapper.fromLoc(mydfwDict['agg'], statsLoc)
                            for statName in selectedStatistics:
                                GrossValues[(statName, expName, varName)] = statsDFW.var(statName).to_numpy()
                    for expName in self.expNames:
                        print('Gross statistics for')
                        print('experiment=>'+expName)
                        if len(binValsMap) > 1:
                            print('binVal=>'+binVal)
                        print(' variables: ', self.varNames)

                        for statName in selectedStatistics:
                            print(statName)
                            tmp = np.asarray([])
                            for varName in self.varNames:
                                tmp = np.append(tmp, GrossValues[(statName, expName, varName)])
                            print(tmp)

AnalysisTypeDict = {
    #Derived from CategoryBinMethodBase(AnalysisBase)
    'CYAxisExpLines': CYAxisExpLines,
    'FCAxisExpLines': FCAxisExpLines,
    'FCAxisExpLinesDiffCI': FCAxisExpLinesDiffCI,
    'FCScoreCard': FCScoreCard,
    'CYAxisFCLines': CYAxisFCLines,
    'CYAxisBinValLines': CYAxisBinValLines,
    #Derived from MultiDimBinMethodBase(AnalysisBase)
    'CYandBinValAxes2D': CYandBinValAxes2D,
    'FCandBinValAxes2D': FCandBinValAxes2D,
    'BinValAxisProfile': BinValAxisProfile,
    'BinValAxisProfileDiffCI': BinValAxisProfileDiffCI,
    'BinValAxes2D': BinValAxes2D,
    #Derived from AnalysisBase
    'BinValAxisPDF': BinValAxisPDF,
    'BinValAxisPDFMultiExp': BinValAxisPDFMultiExp,
    'BinValAxisStatsComposite': BinValAxisStatsComposite,
    'GrossValues': GrossValues,
}

# NOTES:
# (1) FCAxis* types require non-zero forecast length
# (2) CYAxis* types require > 1 analysis cycle
# (3) CYAxisFCLines requires (1) and (2)
# (4) *DiffCI types require more than one experiment

def AnalysisFactory(db, analysisType, diagnosticGroupings):
    myClass = AnalysisTypeDict.get(analysisType, None)
    assert (myClass is not None and inspect.isclass(myClass)), \
        ('\n\nERROR: AnalysisFactory cannot construct ', analysisType, ' without instructions in AnalysisTypeDict')
    return myClass(db, analysisType, diagnosticGroupings)


class Analyses():
    def __init__(self, db, analysisTypes, diagnosticGroupings = {}, nproc = 1):
        self.nproc = nproc
        self.analyses = []
        for anType in analysisTypes:
            self.analyses.append(AnalysisFactory(db, anType, diagnosticGroupings))
        self.logger = logging.getLogger(__name__+'.'+db.DiagSpaceName)
        self.logger.info('Analyses Constructed')

    def analyze(self):
        self.logger.info("Entering Analyses.analyze()")

        if self.nproc > 1:
            workers = mp.Pool(self.nproc)
        else:
            workers = None

        for an in self.analyses:
            an.analyze(workers)

        if workers is not None:
            workers.close()
            workers.join()

        self.logger.info("Exiting Analyses.analyze()")

