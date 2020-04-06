from collections.abc import Iterable
import binning_utils as bu
import binning_configs as bcs
from copy import deepcopy
import datetime as dt
import glob
import numpy as np
import pandas as pd
import plot_utils as pu
import re
import os
import stat_utils as su
import var_utils as vu

singleFCLen = "singleFCLen"
multiFCLen = "multiFCLen"

def dfIndexLevels(df,index):
    mi = df.index.names
    return pu.uniqueMembers(
               df.index.get_level_values(
                   mi.index(index) ).tolist() )

def dfVarVals(df,loc,var):
    return pu.uniqueMembers(df.loc[loc,var].tolist())

class StatsDataBase:
    '''A container class for a pandas DataFrame of
       statistics from multiple cycle and/or forecast times.
    '''
    def __init__(self,conf):
#        self.conf = conf
    ## Examples of directory structures from which this container can extract statistics.
    ## The stats*.nc files are produced by writediagstats_obsspace.py during cycling experiments.
    # multiFCLen ASCII statistics file example for cycling run (on cheyenne):
    #statFile = '/glade/scratch/user/pandac/FC/3dvar/2018041500/0/diagnostic_stats/stats_omb_amsua_n19.nc'
    #            |                         |        |          | |                      |   |         |
    #                       ^                   ^        ^      ^           ^             ^        ^
    #                  expDirectory       expLongName  cyDTime fcTDelta statsFilePrefix DAMethod DiagSpaceName

    # singleFCLen ASCII statistics file example for cycling run (on cheyenne):
    #statFile = '/glade/scratch/user/pandac/DA/3dvar/2018041500/diagnostic_stats/stats_3dvar_bumpcov_amsua_n19.nc'
    #            |                         |        |          |                      |             |         |
    #                       ^                   ^        ^                ^                  ^           ^
    #                  expDirectory       expLongName  cyDTime      statsFilePrefix     DAMethod   DiagSpaceName

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
        self.DAMethods = conf['DAMethods']

        # group of selected plots (single/multi)
        self.plotGroup = conf['plotGroup']

        # selected DiagSpace (ObsSpace name or ModelSpace name)
        self.DiagSpaceName = conf['DiagSpaceName']

        self.statsFilePrefix = 'diagnostic_stats/'+su.statsFilePrefix

        self.fcTDeltas_dir = []
        self.fcTDeltas = []
        dumTimeDelta = fcTDeltaFirst
        while dumTimeDelta <= fcTDeltaLast:
            #TODO: define FC directory names in terms of seconds or d_hh-mm-ss
            #      in order to allow for increments < 1 day --> modify workflow scripts
            #fc_date_str = str(dumTimeDelta.total_seconds())
            #fc_date_str = dumTimeDelta.__str__() #[d days, h:mm:ss] --> need to parse/format to d_hh-mm-ss
            fc_date_str = str(dumTimeDelta.days)
            self.fcTDeltas_dir.append(fc_date_str)
            self.fcTDeltas.append(dumTimeDelta)
            dumTimeDelta = dumTimeDelta + fcTimeInc

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

        # Retrieve list of DiagSpaceNames from files available for all experiments
        # TODO: Only populate DiagSpaceNames if all cyDTimes and fcTDeltas meet
        #       these conditions:
        #         (1) all stats files are present
        #         (2) all stats files contain the same number of rows (easier)
        #             or equivalent rows (harder)
        # TODO: add the capability to calculate the stats when files are missing/incomplete
        #       and then output to correct file name
        expsDiagSpaceNames = []
        for expName in self.expLongNames:
            dateDir = self.cyDTimes_dir[0]
            if self.plotGroup == multiFCLen:
                dateDir = dateDir+'/'+self.fcTDeltas_dir[0]

            FILEPREFIX0 = self.expDirectory + expName +'/'+dateDir+'/' \
                          +self.statsFilePrefix+self.DAMethods[0]+"_"

            DiagSpaceNames = []
            for File in glob.glob(FILEPREFIX0+'*.nc'):
               DiagSpaceName = re.sub(".nc","",re.sub(FILEPREFIX0,"",File))
               if DiagSpaceName == self.DiagSpaceName:
                   DiagSpaceNames.append(DiagSpaceName)
            expsDiagSpaceNames.append(DiagSpaceNames)

        # Remove DiagSpaceNames that are not common to all experiments
        self.availDiagSpaceNames = deepcopy(expsDiagSpaceNames[0])
        if len(expsDiagSpaceNames) > 1:
            for expDiagSpaceNames in expsDiagSpaceNames[1:]:
                for DiagSpaceName in expDiagSpaceNames:
                    if DiagSpaceName not in expDiagSpaceNames:
                        self.availDiagSpaceNames.remove(DiagSpaceName)


    def initialize(self):
        self.available = False
        if (len(self.availDiagSpaceNames) < 1):
            print("\nWARNING: stats files not available for creating a StatsDataBase"+
                  "\nobject for the selected DiagSpace => "+self.DiagSpaceName)
            return

        assert len(self.availDiagSpaceNames) == 1, (
            "\n\nERROR: only one DiagSpaceName per object is allowed.")
        self.available = True

        print("")
        print("\n\n\n=========================================================")
        print("Initialize DataFrame for DiagSpaceName = "+self.DiagSpaceName)
        print("=========================================================")

        # Read stats and make figures for this DiagSpaceName
        print("\nReading NC intermediate files into common Pandas DataFrame...")

        #TODO: appending numpy arrays requires reallocation in every inner loop (expensive).
        #      Instead, pre-allocate the numpy arrays after determining their global lengths
        #      nRows * nExp * nFC * nCY, as long as nRows is consistent across all times
        dsDict = {}
        dsDict['expName'] = np.asarray([])
        dsDict['fcTDelta'] = np.asarray([])
        dsDict['cyDTime'] = np.asarray([])
        for attribName in su.fileStatAttributes:
            dsDict[attribName] = np.asarray([])
        for statName in su.allFileStats:
            dsDict[statName] = np.asarray([])

        for iexp, expName in enumerate(self.expNames):
            expPrefix = self.expDirectory + self.expLongNames[iexp] +'/'
            ncStatsFile = self.statsFilePrefix+self.DAMethods[iexp]+'_'+self.DiagSpaceName+'.nc'
            for ifc, fcTDelta in enumerate(self.fcTDeltas):
                fcstr=fcTDelta.__str__()
                for icy, cyDTime in enumerate(self.cyDTimes):
                    #Read all stats/attributes from NC file for DiagSpaceName, ExpName, fcTDelta, cyDTime
                    dateDir = self.cyDTimes_dir[icy]
                    if self.plotGroup == multiFCLen:
                        dateDir = dateDir+'/'+self.fcTDeltas_dir[ifc]
                    cyStatsFile = expPrefix+dateDir+'/'+ncStatsFile

                    if os.path.exists(cyStatsFile):
                        statsDict = su.read_stats_nc(cyStatsFile)
                    else:
                        print("\nWARNING: stats file does not exist: "+str(cyStatsFile)+
                              "\n    -> this time will be excluded from all statistics and figures")
                        continue
                    nrows = len(statsDict[su.fileStatAttributes[0]])
                    dsDict['expName'] = \
                        np.append(dsDict['expName'], [expName] * nrows)
                    dsDict['fcTDelta'] = \
                        np.append(dsDict['fcTDelta'], [fcTDelta] * nrows)
                    dsDict['cyDTime'] = \
                        np.append(dsDict['cyDTime'], [cyDTime] * nrows)

                    for attribName in su.fileStatAttributes:
                        dsDict[attribName] = \
                            np.append(dsDict[attribName],statsDict[attribName])
                    for statName in su.allFileStats:
                        dsDict[statName] = \
                            np.append(dsDict[statName],statsDict[statName])

        #Convert dsDict to DataFrame
        dsDF = pd.DataFrame.from_dict(dsDict)

        indexNames = ['expName','fcTDelta','cyDTime','DiagSpaceGrp',
                      'varName','diagName','binVar','binVal','binMethod']

        dsDF.set_index(indexNames,inplace=True)
        dsDF.sort_index(inplace=True)

        ##  diagspace
        self.DiagSpaceGrp = dsDF.index.levels[indexNames.index('DiagSpaceGrp')]

        # remove the DiagSpaceGrp dimension, because it's common across all rows
        #       expName      fcTDelta    cyDTime                     varName     diagName    binVar      binVal      binMethod
        dsLoc = (slice(None),slice(None),slice(None),self.DiagSpaceGrp[0],slice(None),slice(None),slice(None),slice(None),slice(None))
        self.dfw = DFWrapper(dsDF.xs(dsLoc))


        ## diagnostics
        self.allDiagNames = self.dfw.levels('diagName')


        ##  variables
        # get varNames and sort alphabetically
        varNames = self.dfw.levels('varName')
        nVars = len(varNames)
        indices = list(range(nVars))

        # sort by channel number (int) for radiances
        chlist = ['']*nVars
        for ivar, varName in enumerate(varNames):
            for c in list(range(len(varName))):
                sub = varName[c:]
                if pu.isint(sub):
                    chlist[ivar] = int(sub)
                    break
        if '' in chlist:
            indices.sort(key=varNames.__getitem__)
        else:
            indices.sort(key=chlist.__getitem__)
        self.varNames = list(map(varNames.__getitem__, indices))
        self.chlist = list(map(chlist.__getitem__, indices))

        # extract units for each varName from varUnits DF column
        self.varUnitss = []
        varLoc = {}
        varLoc['expName'] = self.expNames[0]
        varLoc['fcTDelta'] = self.fcTDeltas[0]
        varLoc['cyDTime'] = self.cyDTimes[0]
        varLoc['binVar'] = vu.obsVarQC
        varLoc['binVal'] = bcs.goodFlagName
        varLoc['binMethod'] = bu.goodQCMethod
        goodDiags = self.dfw.levels('diagName',varLoc)
        varLoc['diagName'] = goodDiags[0]

        for varName in self.varNames:
            varLoc['varName'] = varName
            units = self.dfw.uniquevals('varUnits',varLoc)
            assert len(units) == 1, ("\n\nERROR: too many units values for varName = "+varName,
                                    units, varLoc)
            self.varUnitss.append(units[0])


        ##  bin values --> combination of numerical and string, all stored as strings
        self.allBinVals = self.dfw.levels('binVal')

        # convert allBinVals to numeric type that can be used as axes values
        self.binNumVals = []
        self.binNumVals2DasStr = []
        for binVal in self.allBinVals:
            if pu.isint(binVal):
                self.binNumVals.append(int(binVal))
                self.binNumVals2DasStr.append(binVal)
            elif pu.isfloat(binVal):
                self.binNumVals.append(float(binVal))
                self.binNumVals2DasStr.append(binVal)
            else:
                self.binNumVals.append(np.NaN)

    def loc(self,locDict,var=None):
        return DFWrapper(self.dfw.loc(locDict,var))

    ## not used yet
    # def agg(self,aggovers=['cyDTime']):
    #     return DFWrapper(self.dfw.aggStats(groupby))

class DFWrapper:
    def __init__(self,df):
        self.df = df
        self.indexNames = list(df.index.names)

    @classmethod
    def fromLoc(cls,other,locDict,var=None):
        return cls(other.locdf(other.locTuple(locDict),var))

    @classmethod
    def fromAggStats(cls,other,aggovers):
        return cls(other.aggStats(aggovers))

    def locTuple(self,locDict={}):
        Loc = ()
        for index in list(locDict.keys()):
            assert index in self.indexNames,(
                "\n\nERROR: index name not in the multiindex, index = "+index
                +", indexNames = ",self.indexNames)

        for index in self.indexNames:
            indL = list(Loc)
            if index not in locDict:
                indL.append(slice(None))
            elif locDict[index] is None:
                indL.append(slice(None))
            elif (isinstance(locDict[index],Iterable) and
                not isinstance(locDict[index],str)):
                indL.append(locDict[index])
            else:
                indL.append([locDict[index]])
            Loc = tuple(indL)
        return Loc

    def locdf(self,Loc,var=None):
        if var is None:
            return self.df.loc[Loc,:]
        else:
            return self.df.loc[Loc,var]

    def levels(self,index,locDict={}):
        newDF = self.locdf(self.locTuple(locDict))
        return dfIndexLevels(newDF,index)

    def loc(self,locDict,var=None):
        return self.locdf(self.locTuple(locDict),var)

    def var(self,var):
        return self.loc({},var=var)

    def uniquevals(self,var,locDict={}):
        return pu.uniqueMembers(self.loc(locDict,var).tolist())

    def min(self,locDict,var=None):
        return self.locdf(self.locTuple(locDict),var).dropna().min()

    def max(self,locDict,var):
        return self.locdf(self.locTuple(locDict),var).dropna().max()

    def aggStats(self,aggovers):
        groupby = deepcopy(self.indexNames)
        for aggover in aggovers:
            assert aggover in self.indexNames, (
                "\n\nERROR: aggover argument not in the multiindex, aggover = "+aggover
                +", indexNames = ",self.indexNames)
            if aggover in groupby: groupby.remove(aggover)
        return self.df.groupby(groupby).apply(su.aggStatsDF)
