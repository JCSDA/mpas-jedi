#!/usr/bin/env python3
import binning_utils as bu
from collections.abc import Iterable
from collections import defaultdict
from copy import deepcopy
import datetime as dt
import glob
import logging
import multiprocessing as mp
import numpy as np
import pandas as pd
import plot_utils as pu
import re
import os
import stat_utils as su
from typing import List
import var_utils as vu

class MultipleBinnedStatistics():
    def __init__(self, nrows):
        self.nrows = nrows
        self.values = {}
        self.values['expName'] = np.empty(nrows, np.chararray)
        self.values['fcTDelta'] = np.empty(nrows, dt.timedelta)
        self.values['cyDTime'] = np.empty(nrows, dt.datetime)
        for attribName in su.fileStatAttributes:
            self.values[attribName] = np.empty(nrows, np.chararray)
        for statName in su.allFileStats:
            self.values[statName] = np.empty(nrows, np.float)

    @classmethod
    def read(cls, statsFile, expName, fcTDelta, cyDTime):
        statsDict = statsFile.read()
        nrows = len(statsDict[su.fileStatAttributes[0]])
        statsDict['expName'] = np.full(nrows, expName)
        statsDict['fcTDelta'] = np.full(nrows, fcTDelta)
        statsDict['cyDTime'] = np.full(nrows, cyDTime)

        new = cls(nrows)
        for key, val in statsDict.items():
            assert key in new.values, "ERROR: MultipleBinnedStatistics.read() "+key+" not in values"
            new.values[key][:] = val[:]
        return new

    @classmethod
    def concatasync(cls, asyncresults):
        nrows = 0
        for asyncresult in asyncresults:
            nrows += asyncresult.get().nrows

        new = cls(nrows)
        srow = 0
        for asyncresult in asyncresults:
            new.insert(asyncresult.get(), srow)
            srow += asyncresult.get().nrows
        return new

    def insert(self, other, srow):
        assert srow >= 0, ("Error: can only insert MultipleBinnedStatistics rows >= 0, not ", srow)
        erow = srow + other.nrows - 1
        assert erow < self.nrows, ("Error: can only insert MultipleBinnedStatistics rows < ", self.nrows, ", not ", erow)
        for key, val in other.values.items():
            if isinstance(val, Iterable):
                assert key in self.values, key+" not in MultipleBinnedStatistics"
                self.values[key][srow:erow+1] = val[:]

    def destroy(self):
        del self.values


def dfIndexLevels(df, index):
    mi = df.index.names
    return pu.uniqueMembers(
               df.index.get_level_values(
                   mi.index(index) ).tolist() )


def dfVarVals(df, loc, var):
    return pu.uniqueMembers(df.loc[loc, var].tolist())


class StatsDB:
    '''A container class for a pandas DataFrame of
       statistics from multiple cycle and/or forecast times.
    '''
    def __init__(self, conf):
#        self.conf = conf
    ## Examples of directory structures from which this container can extract statistics.
    ## The su.BinnedStatisticsFile's are produced by DiagnoseObsStatistics.py and writediagstats_modelspace.py during
    ## cycling experiments.  The directory structure is controlled by the self.statsDirectory method.  The default
    # conventions for cycling runs (on cheyenne) are as follows:
    # (hasFCLenDir == True or self.fcTDeltas[-1] > self.fcTDeltas[0]):
    #statFile = '/glade/scratch/user/pandac/FC/3dvar/2018041500/{fcDirFormats}/diagnostic_stats/stats_omb_amsua_n19.nc'
    #            |                         |        |          |              |                |     |   |         |
    #                       ^                   ^      ^         ^              ^                 ^    ^        ^
    #                expDirectory       expLongName cyDTime  fcTDelta statsFileSubDir statsFilePrefix appIdentifier DiagSpaceName

    # (hasFCLenDir == False and self.fcTDeltas[-1] == self.fcTDeltas[0]):
    #statFile = '/glade/scratch/user/pandac/DA/3dvar/2018041500/diagnostic_stats/stats_3dvar_bumpcov_amsua_n19.nc'
    #            |                         |        |          |                |     |             |         |
    #                       ^                   ^        ^       ^                 ^              ^        ^
    #                  expDirectory       expLongName  cyDTime statsFileSubDir statsFilePrefix  appIdentifier  DiagSpaceName
        self.available = False

        # selected DiagSpace (ObsSpace name or ModelSpace name)
        self.DiagSpaceName = conf['DiagSpaceName']
        self.logger = logging.getLogger(__name__+'.'+self.DiagSpaceName)

        # cycle DateTimes
        firstCycleDTime = conf['firstCycleDTime']
        lastCycleDTime = conf['lastCycleDTime']
        cyTimeInc = conf['cyTimeInc']
        assert cyTimeInc > dt.timedelta(0), "cyTimeInc must be > 0"

        # forecast TimeDeltas
        fcTDeltaFirst = conf['fcTDeltaFirst']
        fcTDeltaLast = conf['fcTDeltaLast']
        fcTimeInc = conf['fcTimeInc']
        assert fcTimeInc > dt.timedelta(0), "fcTimeInc must be > 0"

        # experiment info
        self.expDirectory = conf['expDirectory']
        self.expLongNames = conf['expLongNames']
        self.expNames = conf['expNames']
        self.cntrlExpName = conf['cntrlExpName']
        self.cntrlExpIndex = self.expNames.index(self.cntrlExpName)
        self.noncntrlExpNames = [x for x in self.expNames if x != self.cntrlExpName]
        self.logger.info('Control Experiment: '+self.cntrlExpName)
        self.logger.info(('Non-control Experiment(s): ', self.noncntrlExpNames))

        self.appIdentifiers = conf['appIdentifiers']

        self.diagnosticConfigs = conf['diagnosticConfigs']

        self.statsFileSubDirs = conf['statsFileSubDirs']

        fcDirFormats = conf['fcDirFormats']

        self.fcTDeltas = []
        self.fcTDeltas_dir = defaultdict(list)
        self.fcTDeltas_totmin = []
        dumTimeDelta = fcTDeltaFirst
        while dumTimeDelta <= fcTDeltaLast:
            for expName, fcDirFormat in list(zip(self.expNames, fcDirFormats)):
                self.fcTDeltas_dir[expName] += [TDelta_dir(dumTimeDelta, fcDirFormat)]
            self.fcTDeltas_totmin.append(TDelta_dir(dumTimeDelta, "%m"))
            # self.fcTDeltas_totsec.append(TDelta_dir(dumTimeDelta, "%s"))

            self.fcTDeltas.append(dumTimeDelta)
            dumTimeDelta = dumTimeDelta + fcTimeInc

        # whether directory structure includes forecast length
        self.hasFCLenDir = conf['hasFCLenDir']
        if self.fcTDeltas[-1] > self.fcTDeltas[0]: self.hasFCLenDir = True

        self.cyDTimes_dir = []
        self.cyDTimes = []
        dumDateTime = firstCycleDTime
        while dumDateTime <= lastCycleDTime:
            cy_date_str = "{:04d}".format(dumDateTime.year)  \
                        + "{:02d}".format(dumDateTime.month) \
                        + "{:02d}".format(dumDateTime.day)   \
                        + "{:02d}".format(dumDateTime.hour)
            self.cyDTimes_dir.append(cy_date_str)
            self.cyDTimes.append(dumDateTime)
            dumDateTime = dumDateTime + cyTimeInc

        # Check if first statsFile is present for self.DiagSpaceName and for all experiments
        # TODO: self.available requires all cyDTimes and fcTDeltas to meet
        #       these conditions:
        #         (1) all statsFiles are present
        #         (2) all statsFiles generate DataFrames that contain the same number of rows (easier)
        #             or equivalent rows (harder)
        # TODO: add the capability to calculate the stats when files are missing/incomplete
        #       and then output to correct file name
        expsDiagSpaceNames = []
        for expName, expLongName, statsFileSubDir, appIdentifier in list(zip(
            self.expNames, self.expLongNames, self.statsFileSubDirs, self.appIdentifiers)):

            directory = self.statsDirectory(expLongName, self.cyDTimes_dir[0], self.fcTDeltas_dir[expName][0], statsFileSubDir)

            statsFile = su.BinnedStatisticsFile(
                            appIdentifier=appIdentifier,
                            DiagSpace=self.DiagSpaceName,
                            directory=directory)

            if (not os.path.exists(statsFile.fileName()+'.nc') and
                not os.path.exists(statsFile.fileName()+'.h5')):
                self.logger.warning("stats file not available at first cycle/forecast combination"+
                                    "DiagSpace = "+self.DiagSpaceName+
                                    " and expName ="+expName)
                self.logger.warning("attempted to find stats file at "+
                                    statsFile.fileName()+'.h5')

                return

        self.available = True

    def read(self, np=1):
        if not self.available: return

        self.logger.info("=====================================================")
        self.logger.info("Construct pandas dataframe from static database files")
        self.logger.info("=====================================================")

        nprocs = min(mp.cpu_count(), np)

        # Read stats for this DiagSpaceName
        self.logger.info("Reading intermediate statistics files")
        self.logger.info("with "+str(nprocs)+" out of "+str(mp.cpu_count())+" processors")
        workers = mp.Pool(nprocs)
        binnedStatsParts = []
        for cyDTime, cyDTime_dir in list(zip(self.cyDTimes, self.cyDTimes_dir)):
            self.logger.info("  Working on cycle time "+str(cyDTime))
            missingFiles = []

            for expName, expLongName, statsFileSubDir, appIdentifier in list(zip(
                self.expNames, self.expLongNames, self.statsFileSubDirs, self.appIdentifiers)):

                for fcTDelta, fcTDelta_dir in list(zip(
                    self.fcTDeltas, self.fcTDeltas_dir[expName])):

                    #Read all stats/attributes from su.BinnedStatisticsFile for ExpName, fcTDelta, cyDTime
                    directory = self.statsDirectory(expLongName, cyDTime_dir, fcTDelta_dir, statsFileSubDir)
                    statsFile = su.BinnedStatisticsFile(
                                    appIdentifier=appIdentifier,
                                    DiagSpace=self.DiagSpaceName,
                                    directory=directory)
                    if (os.path.exists(statsFile.fileName()+'.nc') or
                        os.path.exists(statsFile.fileName()+'.h5')):
                        binnedStatsParts.append(workers.apply_async(MultipleBinnedStatistics.read,
                            args = (statsFile, expName, fcTDelta, cyDTime)))
                    else:
                        missingFiles.append(statsFile.fileName())

            if len(missingFiles) > 0:
                self.logger.warning("The following files do not exist.  Matching times are excluded from the statistsics.")
                for File in missingFiles:
                    self.logger.warning(File)
        workers.close()
        workers.join()

        self.logger.info("Concatenating statistics sub-dictionaries from multiple processors")
        binnedStats = MultipleBinnedStatistics.concatasync(binnedStatsParts)

        for expName in self.expNames:
            assert expName in binnedStats.values['expName'], \
                'ERROR: no statsFiles found for expName = '+expName

        ## Convert binnedStats to DataFrame
        self.logger.info("Constructing a dataframe from statistics dictionary")
        dsDF = pd.DataFrame.from_dict(binnedStats.values)
        binnedStats.destroy()
        del binnedStatsParts

        self.logger.info("Sorting the dataframe index")

        indexNames = ['expName', 'fcTDelta', 'cyDTime', 'DiagSpaceGrp',
                      'varName', 'diagName', 'binVar', 'binVal', 'binMethod']

        dsDF.set_index(indexNames, inplace=True)
        dsDF.sort_index(inplace=True)

        self.logger.info("Extracting index values")
        ##  diagspace group
        self.DiagSpaceGrp = dsDF.index.levels[indexNames.index('DiagSpaceGrp')]

        # remove the DiagSpaceGrp dimension, because it's common across all rows and therefore extraneous
        #       expName      fcTDelta    cyDTime                     varName     diagName    binVar      binVal      binMethod
        dsLoc = (slice(None), slice(None), slice(None), self.DiagSpaceGrp[0], slice(None), slice(None), slice(None), slice(None), slice(None))
        self.dfw = DFWrapper(dsDF.xs(dsLoc))

        # drop non-required diagnostics
        requiredDiagnostics = set(list(self.diagnosticConfigs.keys()))
        for config in self.diagnosticConfigs.values():
            requiredDiagnostics = set(list(requiredDiagnostics) + list(config['requiredDiagnostics']))
        availableDiagnostics = set(self.dfw.levels('diagName'))
        keepDiagnostics = availableDiagnostics & requiredDiagnostics
        if len(keepDiagnostics) < 1:
            message = 'No remaining diagnostics! availableDiagnostics: '
            for diag in availableDiagnostics: message += diag+', '
            message += '; requiredDiagnostics: '
            for diag in requiredDiagnostics: message += diag+', '
            self.logger.error(message)
        self.dfw = DFWrapper.fromLoc(self.dfw, {'diagName': keepDiagnostics})
        # TODO: would rather drop non-required diagnostics in place to avoid memory overhead, but gives warning message
        # for diagName in availableDiagnostics:
        #     if diagName not in requiredDiagnostics:
        #         # drop unused diagName from dfw
        #         self.dfw.df.drop(diagName, level='diagName', inplace=True)

        # initialize self attributes
        self.initAttributes()

        # add non-aggregated derived diagnostics as needed
        createORreplaceDerivedDiagnostics(self.dfw, self.diagnosticConfigs)

        self.logger.info('availableDiagnostics: '+str(self.dfw.levels('diagName')))

    def initAttributes(self):
        ## diagnostics (currently unused)
        #self.containedDiagNames = self.dfw.levels('diagName')

        ##  variables
        # get varNames and sort alphabetically
        varNames = self.dfw.levels('varName')
        nVars = len(varNames)
        indices = list(range(nVars))

        # sort by channel number (int) for radiances
        chlist = [-1]*nVars
        for ivar, varName in enumerate(varNames):
            for c in list(range(len(varName))):
                sub = varName[c:]
                if pu.isint(sub):
                    chlist[ivar] = int(sub)
                    break
        if any(np.array(chlist) < 0):
            indices.sort(key=varNames.__getitem__)
        else:
            indices.sort(key=chlist.__getitem__)
        self.varNames = list(map(varNames.__getitem__, indices))
        self.chlist = list(map(chlist.__getitem__, indices))

        ## extract units for each varName from varUnits DF column
        self.varUnitss = []
        varLoc = {}
        #varLoc['fcTDelta'] = self.fcTDeltas[0]
        #varLoc['cyDTime'] = self.cyDTimes[0]

        for varName in self.varNames:
            varLoc['varName'] = varName
            units = self.dfw.uniquevals('varUnits', varLoc)
            #assert len(units) == 1, ("\n\nERROR: too many units values for varName = "+varName,
                                    #units, varLoc)
            self.varUnitss.append(units[0])

        ##  bin values --> combination of numerical and string, all stored as strings
        self.allBinStrVals = self.dfw.levels('binVal')

        # convert allBinStrVals to numeric type that can be used as axes values
        self.allBinNumVals = []
        self.allBinNumVals2DasStr = []
        for binVal in self.allBinStrVals:
            # int
            if pu.isint(binVal):
                self.allBinNumVals.append(int(binVal))
                self.allBinNumVals2DasStr.append(binVal)
            # float
            elif pu.isfloat(binVal):
                self.allBinNumVals.append(float(binVal))
                self.allBinNumVals2DasStr.append(binVal)
            else:
                self.allBinNumVals.append(vu.miss_i)
                # comma-separated lists of float/int
                if ',' in binVal:
                    binVals = binVal.split(',')
                    if all([(pu.isint(b) or pu.isfloat(b)) for b in binVals]):
                        self.allBinNumVals2DasStr.append(binVal)

    def statsDirectory(self, expLongName, cyDTime, fcTDelta, statsFileSubDir):
        dateDir = cyDTime
        if self.hasFCLenDir: dateDir += '/'+fcTDelta
        return self.expDirectory+'/'+expLongName +'/'+dateDir+'/'+statsFileSubDir

    def appendDF(self, newDiagDF):
        self.dfw.append(newDiagDF)
        self.initAttributes()

    def loc(self, locDict, var=None):
        return DFWrapper(self.dfw.loc(locDict, var))


    ## not used yet, but should work
    # def agg(self, aggovers=['cyDTime']):
    #     return DFWrapper(self.dfw.aggStats(groupby))


def createORreplaceDerivedDiagnostics(dfw, diagnosticConfigs):
    for diagName, diagnosticConfig in diagnosticConfigs.items():
        if 'DerivedDiagnostic' in diagnosticConfig:
            availableDiagnostics = dfw.levels('diagName')
            if diagName in availableDiagnostics:
                # drop derived diagName from dfw
                dfw.df.drop(diagName, level='diagName', inplace=True)

            # create then append DataFrame with derived diagName
            derivedDiagDF = diagnosticConfig['DerivedDiagnostic'].evaluate(dfw)
            dfw.append(derivedDiagDF)


class DFWrapper:
    def __init__(self, df):
        self.df = df
        self.indexNames = list(self.df.index.names)

    def __str__(self):
      with pd.option_context('display.max_rows', None,
                             'display.max_columns', None,
                             'display.precision', 3,
                             ):
        return self.df.to_string()

    @classmethod
    def fromLoc(cls, other, locDict, var=None):
        return cls(other.locdf(other.locTuple(locDict), var))

    @classmethod
    def fromAggStats(cls, other, aggovers):
        return cls(other.aggStats(aggovers))

    def append(self, otherDF = None):
        if otherDF is None: return

        #Add otherDF (DataFrame object) to self.df
        # adds new column names as needed
        # adds meaningless NaN entries in columns that do not overlap between two DF's
        # TODO: reduce memory footprint of NaN's via modifications to external data flows
        appendDF = otherDF.copy(True)

        selfColumns = list(self.df.columns)
        appendColumns = list(appendDF.columns)

        selfNRows = len(self.df.index)
        for column in appendColumns:
            if column not in selfColumns:
                self.df.insert(len(list(self.df.columns)), column, [np.NaN]*selfNRows)

        appendNRows = len(appendDF.index)
        for column in selfColumns:
            if column not in appendColumns:
                appendDF.insert(len(list(appendDF.columns)), column, [np.NaN]*appendNRows)

        self.df = self.df.append(appendDF, sort=True)

    def locTuple(self, locDict={}):
        Loc = ()
        for index in list(locDict.keys()):
            assert index in self.indexNames,(
                "\n\nERROR: index name not in the multiindex, index = "+index
                +", indexNames = ", self.indexNames)

        for index in self.indexNames:
            indL = list(Loc)
            if index not in locDict:
                indL.append(slice(None))
            elif locDict[index] is None:
                indL.append(slice(None))
            elif (isinstance(locDict[index], Iterable) and
                not isinstance(locDict[index], str)):
                indL.append(list(locDict[index]))
            else:
                indL.append(list([locDict[index]]))
            Loc = tuple(indL)
        return Loc

    def locdf(self, Loc, var=None):
        if var is None:
            return self.df.loc[Loc, :]
        else:
            return self.df.loc[Loc, var]

    def levels(self, index, locDict={}):
        newDF = self.locdf(self.locTuple(locDict))
        return dfIndexLevels(newDF, index)

    def loc(self, locDict, var=None):
        return self.locdf(self.locTuple(locDict), var)

    def loc1(self, locDict, var=None):
        s = self.loc(locDict, var).to_numpy()
        if isinstance(s, Iterable):
          if len(s) == 1:
            return s[0]
          else:
            return np.NaN
          #assert len(s) == 1, "DFWrapper::loc0, locDict/var must return single location only"
          #return s[0]
        else:
          return s

    def var(self, var):
        return self.loc({}, var=var)

    def uniquevals(self, var, locDict={}):
        return pu.uniqueMembers(self.loc(locDict, var).tolist())

    def min(self, locDict, var=None):
       return self.loc(locDict, var).dropna().min()

    def max(self, locDict, var):
        return self.loc(locDict, var).dropna().max()

    def aggStats(self, aggovers):
        groupby = deepcopy(self.indexNames)
        for aggover in aggovers:
            assert aggover in self.indexNames, (
                "\n\nERROR: aggover argument not in the multiindex, aggover = "+aggover
                +", indexNames = ", self.indexNames)
            if aggover in groupby: groupby.remove(aggover)
        return self.df.groupby(groupby).apply(su.aggStatsDF)


def TDelta_dir(tdelta, fmt):
    subs = {}
    fmts = {}
    i = '{:d}'
    i02 = '{:02d}'

    # "%D %HH:%MM:%SS"
    subs["D"] = tdelta.days
    fmts["D"] = i

    subs["HH"], hrem = divmod(tdelta.seconds, 3600)
    fmts["HH"] = i02

    subs["MM"], subs["SS"] = divmod(hrem, 60)
    fmts["MM"] = i02
    fmts["SS"] = i02

    ts = int(tdelta.total_seconds())

    # "%h"
    subs["h"], hrem = divmod(ts, 3600)
    fmts["h"] = i

    # "%MIN:%SEC"
    subs["MIN"], subs["SEC"] = divmod(ts, 60)
    fmts["MIN"] = i
    fmts["SEC"] = i02

    subs["m"] = subs["MIN"]
    fmts["m"] = fmts["MIN"]

    # "%s"
    subs["s"] = ts
    fmts["s"] = i

    out = fmt
    for key in subs.keys():
        out = out.replace("%"+key, fmts[key].format(subs[key]))

    return out
