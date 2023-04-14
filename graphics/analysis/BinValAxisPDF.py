#!/usr/bin/env python3

import basic_plot_functions as bpf
import binning_utils as bu
import predefined_configs as pconf
from collections.abc import Iterable
from collections import OrderedDict
from copy import deepcopy
import numpy as np
import plot_utils as pu
import stat_utils as su
import var_utils as vu

from analysis.AnalysisBase import AnalysisBase
import analysis.StatisticsDatabase as sdb

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

    requestAggDFW = True

    subplotWidth = 1.2
    subplotAspect = 1.3

    requiredStatistics = ['Count', 'Mean', 'STD']

    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict):
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
              fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, self.interiorLabels)
              fig_normalized = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, self.interiorLabels)

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
                    interiorLabels = self.interiorLabels,
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
                    interiorLabels = self.interiorLabels,
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

              pu.finalize_fig(fig, str(figPath/filename), self.figureFileType, self.interiorLabels, 0.6)

              filename = ('%s_%s-normalized_BinValAxis_%s-%smin_%s_%s_%s'%(
                         binVar, binMethodCase, self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                         self.DiagSpaceName, fcDiagName, expName))

              pu.finalize_fig(fig_normalized, str(figPath/filename), self.figureFileType, self.interiorLabels, 0.6)

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

    requestAggDFW = True

    subplotWidth = 1.2
    subplotAspect = 1.3

    requiredStatistics = ['Count', 'Mean', 'STD']

    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict):
      super().__init__(db, analysisType, diagnosticGroupings)
      self.binVarDict = OrderedDict()

      ## vu.obsVarNormDep
      self.binVarDict[vu.obsVarNormDep] = {}

      self.binVarDict[vu.obsVarNormDep][bu.identityBinMethod] = {
        'standard gaussian': True,
        'binMethodEqualsAny': [bu.identityBinMethod],
      }

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
            fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, self.interiorLabels)
            fig_normalized = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, self.interiorLabels)

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
                  interiorLabels = self.interiorLabels,
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
                  interiorLabels = self.interiorLabels,
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

            pu.finalize_fig(fig, str(figPath/filename), self.figureFileType, self.interiorLabels, 0.6)

            filename = ('%s_%s-normalized_BinValAxis_%s-%smin_%s_%s_allExp'%(
                       binVar, binMethodCase, self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                       self.DiagSpaceName, fcDiagName))

            pu.finalize_fig(fig_normalized, str(figPath/filename), self.figureFileType, self.interiorLabels, 0.6)
