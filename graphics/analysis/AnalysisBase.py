#!/usr/bin/env python3

import binning_utils as bu
import predefined_configs as pconf
from collections.abc import Iterable
from config import DiagSpaceConfig
from copy import deepcopy
import datetime as dt
import diag_utils as du
import logging
import numpy as np
from pathlib import Path
import plot_utils as pu
import os
import var_utils as vu

import analysis.StatisticsDatabase as sdb

def anWorkingDir(DiagSpace, analysisType):
  return DiagSpace+'_analyses'+'/'+analysisType

class AnalysisBase():
    '''
    Base class for all analysis types
    '''
    ## plot settings
    figureFileType = 'pdf' #['pdf','png']
    interiorLabels = True

    ## yaml writing options
    __dataYAMLMissingFloat = 999. # fill value when ~np.isfinite
    __dataYAMLPrecision = 4 # number of significant digits

    ## Establish default configuration
    blocking = False
    parallelism = False
    requestAggDFW = False

    subplotWidth = 2.5
    subplotAspect = 1.0

    MAX_FC_SUBFIGS = 6
    MAX_FC_LINES = 6

    requiredStatistics = []

    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict = {}):
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

        # TODO(JJG): decide if nproc is needed
        # nproc could be used to initialize workers within a single
        # analysisType to be distributed however they would be most useful.
        # That would allow for more granularity of computational effort
        # but requires applicable analysisType's to be "blocking" due to the nature
        # of multiprocessing Pool's.
        # self.nproc = nproc

        self.blankBinMethodFile = bu.identityBinMethod

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

    def dataYAMLFmtFloat(self, f):
      if np.isfinite(f):
        return float(('{:.'+str(self.__dataYAMLPrecision-1)+'e}').format(f))
      else:
        return self.__dataYAMLMissingFloat

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
