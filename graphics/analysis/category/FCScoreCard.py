#!/usr/bin/env python3

import basic_plot_functions as bpf
import binning_utils as bu
from copy import deepcopy
import diag_utils as du
import itertools as itt
import numpy as np
import plot_utils as pu
import stat_utils as su
import var_utils as vu

from analysis.category.CategoryBinMethodBase import CategoryBinMethodBase
import analysis.StatisticsDatabase as sdb

class FCScoreCard(CategoryBinMethodBase):
    '''
    Similar to FCAxisExpLinesDiffCI, except
      - instead of individual time-series plots, data is displayed in a 2D score card
        with a separate row in 2D subplot placing each FCAxisExpLinesDiffCI subplot axis
      - subplot: column by experiment, row by DiagSpace variable
      -    file: combination of binVar, binMethod, and statistic
    '''

    # used for percent ratio plots
    requestAggDFW = True

    requiredStatistics = ['Count']

    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict):
        super().__init__(db, analysisType, diagnosticGroupings)

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

            ## cloud category bins
            (vu.obsVarCldFracY, bu.cloudbandsMethod): {
              'subplot index': 'varName',
              'binVals': [
                ('clear', 'clear'),
                ('mixed-clrcld', 'mixed'),
                ('cloudy', 'cloudy'),
                ('allsky', 'all'),
              ],
            },
        }

        for key in self.binVarDict:
            if 'onlyStatNames' in self.binVarDict[key]:
                self.binVarDict[key]['onlyStatNames'] += su.bootStrapStats
            else:
                self.binVarDict[key]['onlyStatNames'] = su.bootStrapStats

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

        pu.finalize_fig(fig, str(figPath/filename), self.figureFileType, True, outerWidthFraction, outerHeightFraction)
