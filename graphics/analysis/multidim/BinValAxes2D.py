#!/usr/bin/env python3

import basic_plot_functions as bpf
import binning_utils as bu
import config as conf
import predefined_configs as pconf
from collections import OrderedDict
from copy import deepcopy
import diag_utils as du
from fit2D import fit2D, poly2DEquation
import numpy as np
import plot_utils as pu
import re
import stat_utils as su
from textwrap import indent
import var_utils as vu
import yaml

from analysis.multidim.MultiDimBinMethodBase import MultiDimBinMethodBase
import analysis.StatisticsDatabase as sdb

class BinValAxes2D(MultiDimBinMethodBase):
    '''
    Creates raster maps with binVar binVals on both x- and y-axis
      - applicable to pairs of binned diagnostics (e.g., latitude+vertical dimension, longitude+latitude)
      - this is a valid plot even for a single cycle and/or forecast length (omb)
      - subplot: column by experiment, row by DiagSpace variable
      -    file: combination of binVar2D, (binMethod), statistic, and FC lead time (if applicable)
    '''

    maxBinVarTier = 2
    requestAggDFW = True

    requiredStatistics = ['Count']

    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict):
        super().__init__(db, analysisType, diagnosticGroupings)

        # default binVars
        self.binVarDict = {
            pconf.LonLat2D: {
                'plotfunc': bpf.map2D,
                'binVarTier': 1,
                'subplotWidth': 2.0,
                'subplotAspect': {
                    'default': 0.55,
                    'abi_g16': 0.9,
                    'abi_g18': 0.9,
                    'abi_g19': 0.9,
                    'ahi_himawari8': 0.9,
                    'ahi_himawari9': 0.9,
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

    def _region_extents_for_lonlat(self, lonlat_binvar):
        if lonlat_binvar == pconf.LonLat2D:
            region_var = vu.obsRegionBinVar
            lon_var = vu.lonMeta
            lat_var = vu.latMeta
        elif lonlat_binvar == pconf.ModelLonLat2D:
            region_var = vu.modelRegionBinVar
            lon_var = vu.lonModel
            lat_var = vu.latModel
        else:
            return []

        ds_conf = conf.DiagSpaceConfig.get(self.DiagSpaceName, {})
        ds_binvars = ds_conf.get('binVarConfigs', {})
        region_methods = ds_binvars.get(region_var, [])

        region_configs = pconf.binVarConfigs.get(region_var, {})
        extents = []
        for method_name in region_methods:
            method_config = region_configs.get(method_name, pconf.nullBinMethod)
            filters = method_config.get('filters', [])

            lon_min = None
            lon_max = None
            lat_min = None
            lat_max = None

            for flt in filters:
                if not isinstance(flt, dict):
                    continue
                where = flt.get('where', None)
                variable = flt.get('variable', None)
                bound = flt.get('bounds', None)
                if not np.isscalar(bound) or isinstance(bound, str):
                    continue

                if variable == lon_var:
                    if where in [bu.lessBound, bu.lessEqualBound]:
                        lon_min = float(bound)
                    elif where in [bu.greatBound, bu.greatEqualBound]:
                        lon_max = float(bound)
                elif variable == lat_var:
                    if where in [bu.lessBound, bu.lessEqualBound]:
                        lat_min = float(bound)
                    elif where in [bu.greatBound, bu.greatEqualBound]:
                        lat_max = float(bound)

            if None in [lon_min, lon_max, lat_min, lat_max]:
                continue

            # Normalize longitudes to [-180, 180) for cartopy extents.
            # This maps 180.0 -> -180.0 so dateline-touching regions are handled consistently.
            lon_min = ((lon_min + 180.0) % 360.0) - 180.0
            lon_max = ((lon_max + 180.0) % 360.0) - 180.0

            # Skip regions that wrap the dateline; map2D extent in PlateCarree
            # expects minLon <= maxLon.
            if lon_min > lon_max:
                continue

            lon_pad = 2.0
            lat_pad = 2.0
            extent = [
                max(-180.0, lon_min - lon_pad),
                min(180.0, lon_max + lon_pad),
                max(-90.0, lat_min - lat_pad),
                min(90.0, lat_max + lat_pad),
            ]
            extents.append((method_name, extent))

        return extents

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options):

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
        xVals = np.asarray(xUnique, dtype=float)
        xSort = np.argsort(xVals)
        xVals = xVals[xSort]
        nXVals = len(xVals)
        xValsStr = list(xUnique[xSort])

        yUnique = np.array(pu.uniqueMembers(yCoords))
        yVals = np.asarray(yUnique, dtype=float)
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
            fig = pu.setup_fig(nxplots, nyplots, subplotWidth, subplotAspect, self.interiorLabels)
            zoomed_figs = {}
            if myLoc['binVar'] in [pconf.LonLat2D, pconf.ModelLonLat2D]:
                for region_name, extent in self._region_extents_for_lonlat(myLoc['binVar']):
                    zoomed_figs[region_name] = {
                        'extent': extent,
                        'fig': pu.setup_fig(nxplots, nyplots, subplotWidth, subplotAspect, self.interiorLabels),
                    }

            iplot = 0

            figureData = {}
            figureData['xVals'] = [float(f) for f in xVals]
            figureData['yVals'] = [float(f) for f in yVals]
            figureData['xLabel'] = xLabel
            figureData['yLabel'] = yLabel
            figureData['subplots'] = []

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
                cntrlAggPlaneVals = np.full((nYVals, nXVals), np.nan)
                for ibin, binVal in enumerate(binCoordsLevels):
                    cntrlAggLoc['binVal'] = binVal
                    cntrlAggPlaneVals[yIndex[ibin], xIndex[ibin]] = \
                        dfwDict['agg'].loc1(cntrlAggLoc, statName)

                # intermediate tempdfw reduces extraction time in inner loops
                tempdfw = sdb.DFWrapper.fromLoc(dfwDict['dfw'], planeLoc)

                #subplot loop 2
                dmin_relative = np.nan
                dmax_relative = np.nan
                for expName in self.expNames:
                    if useRelativeDifference:
                        title = varName
                    else:
                        title = varLabel
                    title = expName+'\n'+title
                    expFileName = re.sub(r'\.', '', re.sub(r'\s+', '-', expName))

                    bgstatDiagLabel = bgstatDiagLabel_abs
                    sciTicks = sciTicks_abs
                    logScale = logScale_abs
                    centralValue = centralValue_abs
                    dmin = dmin_abs
                    dmax = dmax_abs

                    planeVals = {}
                    for trait in su.ciTraits:
                        planeVals[trait] = np.full_like(cntrlAggPlaneVals, np.nan)

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
                        expAggPlaneVals.fill(np.nan)
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
                                        su.cimin: np.nan,
                                        su.cimax: np.nan,
                                    }}
                            else:
                                ciVals = {statName: {
                                    su.cimean: expAggPlaneVals[iy, ix],
                                    su.cimin: np.nan,
                                    su.cimax: np.nan,
                                }}

                            for trait in su.ciTraits:
                                t = float(ciVals[statName][trait])
                                # automatically generate relative difference plots for positive-semi-definite statistics
                                if useRelativeDifference:
                                    if normalizingStat != 0 and np.isfinite(normalizingStat):
                                        # divide by cntrlLoc aggregated statName
                                        t /= normalizingStat
                                        if self.relativeErrorType == 'one hundred centered':
                                            t += 1.0
                                        t *= 100.0
                                    else:
                                        t = np.nan

                                planeVals[trait][iy, ix] = t

                        # automatically generate relative difference plots for positive-semi-definite statistics
                        if useRelativeDifference:
                            # TODO: replace with same mechanisms as other DiffCI classes
                            sciTicks = False
                            logScale = False

                            val_diffs, dmin_relative, dmax_relative, centralValue, label = self.relativeErrorFunction(
                              expAggPlaneVals,
                              cntrlAggPlaneVals,
                              dmin_relative,
                              dmax_relative,
                            )
                            
                            dmin = dmin_relative
                            dmax = dmax_relative

                            bgstatDiagLabel = statName.replace('RMS','rms').replace('Mean','mean')+': '+label
                            # check to see if relative differences exceed 3%
                            if (statName == 'RMS'):
                                if not np.isnan(val_diffs).all():
                                    max_all = np.nanmax(val_diffs)
                                    min_all = np.nanmin(val_diffs)
                                    if (abs(max_all) > 3 or abs(min_all) > 3):
                                        self.logger.warning('Experiment '+expName+
                                                            ' RMS variance for ' +varName+
                                                            ' exceeds 3, max:'+str(max_all)+
                                                            ' min:'+str(min_all))

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
                        fitFig = pu.setup_fig(nxFit, nyFit, subplotWidth, subplotAspect, self.interiorLabels)
                        fplot = 0

                        # plot statName
                        bpf.plot2D(
                            fitFig,
                            xVals, yVals, expAggPlaneVals,
                            'Data, '+title, xLabel, yLabel, cLabel,
                            xConfig, yConfig,
                            sciTicks_abs, logScale_abs, centralValue_abs,
                            nyFit, nxFit, nyFit*nxFit, fplot,
                            dmin = dmin, dmax = dmax_abs,
                            interiorLabels = self.interiorLabels)

                        # plot count
                        bpf.plot2D(
                            fitFig,
                            xVals, yVals, countVals,
                            'Data, '+title, xLabel, yLabel, 'Count',
                            xConfig, yConfig,
                            True, True, None,
                            nyFit, nxFit, nyFit*nxFit, fplot+nxFit,
                            dmin = np.nan, dmax = np.nan,
                            interiorLabels = self.interiorLabels)

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
                                interiorLabels = self.interiorLabels)

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
                                interiorLabels = self.interiorLabels)

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
                                dmin = np.nan, dmax = delta/5.,
                                interiorLabels = self.interiorLabels)

                        self.logger.info('\nfit2D L2 norms: '+str(L2Norms))
                        self.logger.info('\nfit2D count-weighted L2 norms: '+str(weightedL2Norms))

                        filename = ('%s_%s_fit2D_%s_%s_%s%s_%smin'%(
                                   varName, expFileName,
                                   myLoc['binVar'], diagnosticGroup, statName,
                                   self.binMethodFile(myLoc['binMethod']),
                                   fcTDelta_totmin))

                        pu.finalize_fig(fitFig, str(figPath/filename), self.figureFileType, self.interiorLabels, 0.6)

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

                        pu.finalize_fig(LFig, str(figPath/filename), self.figureFileType, True, 0.6)

                    subplotData = {}
                    subplotData['varName'] = str(varName)
                    subplotData['expName'] = str(expName)
                    subplotData['title'] = title
                    subplotData['dataLabel'] = cLabel
                    subplotData['dmin'] = self.dataYAMLFmtFloat(dmin)
                    subplotData['dmax'] = self.dataYAMLFmtFloat(dmax)
                    subplotData['contourVals'] = {
                        trait: self.dataYAMLFmtArray(planeVals[trait])
                        for trait in su.ciTraits
                    }
                    figureData['subplots'].append(subplotData)

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
                            interiorLabels = self.interiorLabels)

                        for zoom in zoomed_figs.values():
                            options['plotfunc'](
                                zoom['fig'],
                                xVals, yVals, planeVals[su.cimean],
                                title, cLabel,
                                sciTicks, logScale, centralValue,
                                nyplots, nxplots, nsubplots, iplot,
                                contourValsMinCI = planeVals[su.cimin],
                                contourValsMaxCI = planeVals[su.cimax],
                                dmin = dmin, dmax = dmax,
                                extent = zoom['extent'],
                                interiorLabels = self.interiorLabels)

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
                            interiorLabels = self.interiorLabels)

                    iplot = iplot + 1

            # save each figure
            filename = ('%s%s_BinValAxes2D_%smin_%s_%s_%s'%(
                       myLoc['binVar'],
                       self.binMethodFile(myLoc['binMethod']),
                       fcTDelta_totmin, self.DiagSpaceName,
                       diagnosticGroup, statName))

            pu.finalize_fig(fig, str(figPath/filename), self.figureFileType, self.interiorLabels, xbuffer, ybuffer)

            # save figure data as yaml (zoomed_figs below replot the same
            # planeVals with a different map extent, so no separate sidecar
            # is written for them)
            self.write_figure_yaml(figureData, dataPath, filename)

            for region_name, zoom in zoomed_figs.items():
                zoom_filename = ('%s%s_BinValAxes2D_%smin_%s_%s_%s'%(
                               myLoc['binVar'],
                               self.binMethodFile(region_name),
                               fcTDelta_totmin, self.DiagSpaceName,
                               diagnosticGroup, statName))
                pu.finalize_fig(
                    zoom['fig'],
                    str(figPath/zoom_filename),
                    self.figureFileType,
                    self.interiorLabels,
                    xbuffer,
                    ybuffer,
                )

            if statName in twoDFittingStatistics:

                for expName, e1 in fitEquationConfigs.items():
                    expFileName = re.sub(r'\.', '', re.sub(r'\s+', '-', expName))
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
