from collections.abc import Iterable
from copy import deepcopy
import numpy as np
import os
import plot_utils as pu
import var_utils as vu

#==========================
# names and values of bins
#==========================

## heterogeneous/named bins
# latitude bands, north to south
allNamedLatBands = {}
allNamedLatBands['values']    = ['NPol','NXTro','Tro','SXTro','SPol']
allNamedLatBands['minBounds'] = [60.0, 30.0, -30.0, -90.0, -90.0]
allNamedLatBands['maxBounds'] = [90.0, 90.0,  30.0, -30.0, -60.0]

namedLatBands = {}
namedLatBands['values'] = ['NXTro','Tro','SXTro']
namedLatBands['minBounds'] = []
namedLatBands['maxBounds'] = []

for latBand in namedLatBands['values']:
    iband = allNamedLatBands['values'].index(latBand)
    namedLatBands['minBounds'].append(allNamedLatBands['minBounds'][iband])
    namedLatBands['maxBounds'].append(allNamedLatBands['maxBounds'][iband])
latbandsMethod = 'LatBands'


## homogeneous bins
binLims = {}

binLims[vu.obsVarPrs] = {}
binLims[vu.obsVarPrs]['start']  =    0.0
binLims[vu.obsVarPrs]['finish'] = 1000.0
binLims[vu.obsVarPrs]['step']   =  100.0
binLims[vu.obsVarPrs]['format'] = '{:.0f}'

binLims[vu.obsVarAlt] = {}
binLims[vu.obsVarAlt]['start']  = 1000.0
binLims[vu.obsVarAlt]['finish'] = 50000.0
binLims[vu.obsVarAlt]['step']   = 2000.0
binLims[vu.obsVarAlt]['format'] = '{:.0f}'

binLims[vu.obsVarLat] = {}
binLims[vu.obsVarLat]['start']  = -90.0
binLims[vu.obsVarLat]['finish'] =  90.0
binLims[vu.obsVarLat]['step']   =  10.0
binLims[vu.obsVarLat]['format'] = '{:.0f}'

binLims[vu.obsVarLT] = {}
binLims[vu.obsVarLT]['start']  = 0.0
binLims[vu.obsVarLT]['finish'] = 23.0
binLims[vu.obsVarLT]['step']   = 1.0
binLims[vu.obsVarLT]['format'] = '{:.0f}'

binLims[vu.obsVarSatZen] = {}
binLims[vu.obsVarSatZen]['start']  = 0.0
binLims[vu.obsVarSatZen]['finish'] = 70.0
binLims[vu.obsVarSatZen]['step']   = 5.0
binLims[vu.obsVarSatZen]['format'] = '{:.0f}'

binLims[vu.obsVarCldFrac] = {}
binLims[vu.obsVarCldFrac]['start']  = 0.0
binLims[vu.obsVarCldFrac]['finish'] = 1.0
binLims[vu.obsVarCldFrac]['step']   = 0.05
binLims[vu.obsVarCldFrac]['format'] = '{:.2f}'

binLims[vu.obsVarSCI] = {}
binLims[vu.obsVarSCI]['start']  = 0.0
binLims[vu.obsVarSCI]['finish'] = 60.0
binLims[vu.obsVarSCI]['step']   = 1.0
binLims[vu.obsVarSCI]['format'] = '{:.0f}'

binLims[vu.obsVarNormErr] = {}
binLims[vu.obsVarNormErr]['start']  = -7.0
binLims[vu.obsVarNormErr]['finish'] =  7.0
binLims[vu.obsVarNormErr]['step']   =  0.25
binLims[vu.obsVarNormErr]['format'] = '{:.2f}'

#binLims[vu.modVarLat] = {}
#binLims[vu.modVarLat]['start']  = -90.0
#binLims[vu.modVarLat]['finish'] =  90.0
#binLims[vu.modVarLat]['step']   =  30.0
#binLims[vu.modVarLat]['format'] = '{:.0f}'

#binLims[vu.modVarAlt] = {}
#binLims[vu.modVarAlt]['start']  = 0.0
#binLims[vu.modVarAlt]['finish'] = 50000.0
#binLims[vu.modVarAlt]['step']   = 5000.0
#binLims[vu.modVarAlt]['format'] = '{:.0f}'

for binType, param in binLims.items():
    binBounds = list(np.arange(
        param['start']-0.5*np.abs(param['step']),
        param['finish']+1.5*param['step'],
        param['step']))
    binLims[binType]['minBounds'] = []
    binLims[binType]['maxBounds'] = []
    binLims[binType]['values'] = []
    for ibin in list(range(len(binBounds)-1)):
        binLims[binType]['minBounds'].append(binBounds[ibin])
        binLims[binType]['maxBounds'].append(binBounds[ibin+1])

        binVal = 0.5 * (binBounds[ibin+1] + binBounds[ibin])
        binLims[binType]['values'].append(param['format'].format(binVal))

#jet-stream bin method specifications
P_jet_min = 250.0
P_jet_max = 350.0
P_jet_val = '{:.0f}'.format(0.5 * (P_jet_min + P_jet_max))
PjetMethod = 'P='+P_jet_val+'hPa'

alt_jet_min = 9500.0
alt_jet_max = 10500.0
alt_jet_val = '{:.0f}'.format(0.5 * (alt_jet_min + alt_jet_max))
altjetMethod = 'alt='+alt_jet_val+'m'


# clear/cloudy bin method specifications
clrskyMethod = 'clear-sky'
clrskyThresh = 0.1

cldskyMethod = 'cloud-sky'
cldskyThresh = 0.9

# SCI method specifications
OkamotoMethod       = 'Okamoto'
ScaleOkamotoMethod  = 'ScaledOkamoto'
ModHarnischMethod      = 'ModHarnisch'
ScaleModHarnischMethod = 'ScaledModHarnisch'


#@EffectiveQC* values:
# pass    = 0;   // we like that one!
# missing = 1;   // missing values prevent use of observation
# preQC   = 2;   // observation rejected by pre-processing
# bounds  = 3;   // observation value out of bounds
# domain  = 4;   // observation not within domain of use
# black   = 5;   // observation black listed
# Hfailed = 6;   // H(x) computation failed
# thinned = 7;   // observation removed due to thinning
# diffref = 8;   // metadata too far from reference
# clw     = 9;   // observation removed due to cloud field
# fguess  = 10;  // observation too far from guess
# seaice  = 11;  // observation based sea ice detection, also flags land points
#
# Static list above copied on 3 Oct 2019
# see ufo/src/ufo/filters/QCflags.h for up-to-date list

goodFlag = 0
goodFlagName = 'pass'

badFlags     = [1,         2,         3,         4,
                5,         6,         7,         8,
                9,         10,        11]
badFlagNames = ['missing', 'preQC',   'bounds',  'domain',
                'black',   'Hfailed', 'thinned', 'diffref',
                'clw',     'fguess',  'seaice']


#========================
# binning where functions
#========================

# Note: NaN values have mask set to true by default

def equalBound(x, bound, nanmask=True):
    nanlocs = np.isnan(x)
    mask = np.empty_like(x,dtype=bool)
    mask[~nanlocs] = np.equal(x[~nanlocs],bound)
    mask[nanlocs] = nanmask
    return mask

def notEqualBound(x, bound, nanmask=True):
    nanlocs = np.isnan(x)
    mask = np.empty_like(x,dtype=bool)
    mask[~nanlocs] = np.not_equal(x[~nanlocs],bound)
    mask[nanlocs] = nanmask
    return mask

def lessEqualBound(x, bound, nanmask=True):
    nanlocs = np.isnan(x)
    mask = np.empty_like(x,dtype=bool)
    mask[~nanlocs] = np.less_equal(x[~nanlocs],bound)
    mask[nanlocs] = nanmask
    return mask

def lessBound(x, bound, nanmask=True):
    nanlocs = np.isnan(x)
    mask = np.empty_like(x,dtype=bool)
    mask[~nanlocs] = np.less(x[~nanlocs],bound)
    mask[nanlocs] = nanmask
    return mask

def greatEqualBound(x, bound, nanmask=True):
    nanlocs = np.isnan(x)
    mask = np.empty_like(x,dtype=bool)
    mask[~nanlocs] = np.greater_equal(x[~nanlocs],bound)
    mask[nanlocs] = nanmask
    return mask

def greatBound(x, bound, nanmask=True):
    nanlocs = np.isnan(x)
    mask = np.empty_like(x,dtype=bool)
    mask[~nanlocs] = np.greater(x[~nanlocs],bound)
    mask[nanlocs] = nanmask
    return mask

def betweenBounds(x, bound1, bound2, nanmask=True):
    nanlocs = np.isnan(x)
    belowbounds = lessEqualBound(x,bound1)
    abovebounds = greatEqualBound(x,bound2)
    mask        = np.logical_not(np.logical_or(
                      belowbounds,abovebounds ))
    mask[nanlocs] = nanmask
    return mask


#========================================================
# specific function classes
# i.e., functions of variables contained in the database
#========================================================

class LocalHour:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.dtMeta)
        self.baseVars.append(vu.lonMeta)

    def evaluate(self,dbVals,caseParams):
        TimeStr = dbVals[vu.dtMeta]
        tzOffset = np.divide(dbVals[vu.lonMeta],15.0)

        LH = np.full_like(tzOffset,0.0)
        t0 = binLims[vu.obsVarLT]['start']
        t1 = binLims[vu.obsVarLT]['finish']
        dt = binLims[vu.obsVarLT]['step']
        for ii, time in enumerate(TimeStr):
            ## Expecting time to fit YYYY-MM-DDThh:mm:ssZ
            # YYYY = float(time[0:4])
            # MMo  = float(time[5:7])
            # DD   = float(time[8:10])
            hh   = float(time[11:13])
            mmi  = float(time[14:16]) / 60.0
            ss   = float(time[17:19]) / 3600.0

            LH[ii] = hh + mmi + ss + tzOffset[ii]
            if LH[ii] < t0 - 0.5*dt: LH[ii] += 24.0
            LH[ii] = LH[ii] % 24.0
            if LH[ii] >= t1 + 0.5*dt: LH[ii] -= 24.0

        return LH

class SCIOkamoto:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.clrskyBTDiag)

    def evaluate(self,dbVals,caseParams):
        # Okamoto, et al.
        BTobs = dbVals[caseParams['base2db'][vu.selfObsValue]]
        BTdep = dbVals[caseParams['base2db'][vu.selfDepValue]]
        BTbak = np.add(BTdep,BTobs)
        BTclr = deepcopy(dbVals[caseParams['base2db'][vu.clrskyBTDiag]])
        BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
        SCI = np.multiply( 0.5,
                 np.add(np.abs(np.subtract(BTobs,BTclr)),
                        np.abs(np.subtract(BTbak,BTclr))) )
        return SCI


class ScaledSCIOkamoto:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.clrskyBTDiag)
        self.baseVars.append(vu.cldfracMeta)

    def evaluate(self,dbVals,caseParams):
        # Okamoto, et al.
        BTobs = dbVals[caseParams['base2db'][vu.selfObsValue]]
        BTdep = dbVals[caseParams['base2db'][vu.selfDepValue]]
        BTbak = np.add(BTdep,BTobs)
        BTclr = deepcopy(dbVals[caseParams['base2db'][vu.clrskyBTDiag]])
        BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
        CldFrac = dbVals[caseParams['base2db'][vu.cldfracMeta]]

        SCI = np.multiply( 0.5,
                 np.multiply(CldFrac,
                     np.add(np.abs(np.subtract(BTobs,BTclr)),
                            np.abs(np.subtract(BTbak,BTclr))) ) )
        return SCI


class SCIModHarnisch:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.clrskyBTDiag)

    def evaluate(self,dbVals,caseParams):
        # Modified Harnisch, et al.
        BTobs = dbVals[caseParams['base2db'][vu.selfObsValue]]
        BTdep = dbVals[caseParams['base2db'][vu.selfDepValue]]
        BTbak = np.add(BTdep,BTobs)
        BTclr = deepcopy(dbVals[caseParams['base2db'][vu.clrskyBTDiag]])
        BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
        zeros = np.full_like(BTbak,0.0)
        SCI = np.multiply( 0.5,
                 np.add(np.maximum(zeros,np.subtract(BTclr,BTobs)),
                        np.maximum(zeros,np.subtract(BTclr,BTbak))) )
        return SCI


class ScaledSCIModHarnisch:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.clrskyBTDiag)
        self.baseVars.append(vu.cldfracMeta)

    def evaluate(self,dbVals,caseParams):
        # Modified Harnisch, et al.
        BTobs = dbVals[caseParams['base2db'][vu.selfObsValue]]
        BTdep = dbVals[caseParams['base2db'][vu.selfDepValue]]
        BTbak = np.add(BTdep,BTobs)
        BTclr = deepcopy(dbVals[caseParams['base2db'][vu.clrskyBTDiag]])
        BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
        CldFrac = dbVals[caseParams['base2db'][vu.cldfracMeta]]

        zeros = np.full_like(BTbak,0.0)
        SCI = np.multiply( 0.5,
                 np.multiply(CldFrac,
                     np.add(np.maximum(zeros,np.subtract(BTclr,BTobs)),
                            np.maximum(zeros,np.subtract(BTclr,BTbak))) ) )
        return SCI


class NormalizedError:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.selfErrorValue)

    def evaluate(self,dbVals,caseParams):
        BTerr = dbVals[caseParams['base2db'][vu.selfErrorValue]]
        BTerr[BTerr==0.0] = np.NaN
        BTdep = dbVals[caseParams['base2db'][vu.selfDepValue]]

        return np.divide(BTdep,BTerr)

SCIERRParams = {}
SCIERRParams['abi_g16'] = {}
#This dictionary was generated automatically by
# creating AGGREGATEFC_StatsComposite figures
# with plot_stats_timeseries.py
SCIERRParams['abi_g16'][ (7, 'ModHarnisch') ]   =  {'X': [0.0, 30.61], 'ERR': [2.83, 28.13]}
SCIERRParams['abi_g16'][ (8, 'ModHarnisch') ]   =  {'X': [0.0, 29.45], 'ERR': [2.46, 27.5]}
SCIERRParams['abi_g16'][ (9, 'ModHarnisch') ]   =  {'X': [1.0, 24.16], 'ERR': [1.22, 22.92]}
SCIERRParams['abi_g16'][ (10, 'ModHarnisch') ]   =  {'X': [0.0, 33.15], 'ERR': [4.17, 21.59]}
SCIERRParams['abi_g16'][ (11, 'ModHarnisch') ]   =  {'X': [0.0, 11.03], 'ERR': [1.46, 12.98]}
SCIERRParams['abi_g16'][ (13, 'ModHarnisch') ]   =  {'X': [0.0, 14.47], 'ERR': [1.57, 15.88]}
SCIERRParams['abi_g16'][ (14, 'ModHarnisch') ]   =  {'X': [0.0, 21.49], 'ERR': [1.79, 17.77]}
SCIERRParams['abi_g16'][ (15, 'ModHarnisch') ]   =  {'X': [1.0, 28.29], 'ERR': [1.56, 25.72]}
SCIERRParams['abi_g16'][ (16, 'ModHarnisch') ]   =  {'X': [1.0, 30.01], 'ERR': [1.8, 27.81]}

SCIERRParams['abi_g16'][ (7, 'Okamoto') ]   =  {'X': [1.0, 30.05], 'ERR': [2.18, 28.22]}
SCIERRParams['abi_g16'][ (8, 'Okamoto') ]   =  {'X': [1.0, 29.09], 'ERR': [2.06, 27.64]}
SCIERRParams['abi_g16'][ (9, 'Okamoto') ]   =  {'X': [1.0, 24.19], 'ERR': [1.35, 22.93]}
SCIERRParams['abi_g16'][ (10, 'Okamoto') ]   =  {'X': [1.0, 32.35], 'ERR': [3.1, 21.59]}
SCIERRParams['abi_g16'][ (11, 'Okamoto') ]   =  {'X': [1.0, 11.59], 'ERR': [1.62, 14.68]}
SCIERRParams['abi_g16'][ (13, 'Okamoto') ]   =  {'X': [1.0, 15.35], 'ERR': [1.54, 17.25]}
SCIERRParams['abi_g16'][ (14, 'Okamoto') ]   =  {'X': [0.0, 21.21], 'ERR': [1.33, 17.77]}
SCIERRParams['abi_g16'][ (15, 'Okamoto') ]   =  {'X': [1.0, 28.32], 'ERR': [1.77, 25.72]}
SCIERRParams['abi_g16'][ (16, 'Okamoto') ]   =  {'X': [1.0, 29.77], 'ERR': [1.94, 27.92]}

SCIERRParams['abi_g16'][ (7, 'ScaledModHarnisch') ]   =  {'X': [0.0, 29.96], 'ERR': [4.79, 26.16]}
SCIERRParams['abi_g16'][ (8, 'ScaledModHarnisch') ]   =  {'X': [0.0, 29.55], 'ERR': [4.64, 25.64]}
SCIERRParams['abi_g16'][ (9, 'ScaledModHarnisch') ]   =  {'X': [0.0, 25.92], 'ERR': [3.68, 21.35]}
SCIERRParams['abi_g16'][ (10, 'ScaledModHarnisch') ]   =  {'X': [0.0, 34.42], 'ERR': [4.79, 21.59]}
SCIERRParams['abi_g16'][ (11, 'ScaledModHarnisch') ]   =  {'X': [0.0, 11.43], 'ERR': [1.87, 12.53]}
SCIERRParams['abi_g16'][ (13, 'ScaledModHarnisch') ]   =  {'X': [0.0, 16.17], 'ERR': [2.24, 15.63]}
SCIERRParams['abi_g16'][ (14, 'ScaledModHarnisch') ]   =  {'X': [0.0, 21.63], 'ERR': [2.63, 17.43]}
SCIERRParams['abi_g16'][ (15, 'ScaledModHarnisch') ]   =  {'X': [0.0, 30.94], 'ERR': [4.39, 24.42]}
SCIERRParams['abi_g16'][ (16, 'ScaledModHarnisch') ]   =  {'X': [0.0, 30.11], 'ERR': [4.64, 25.98]}

SCIERRParams['abi_g16'][ (7, 'ScaledOkamoto') ]   =  {'X': [0.0, 29.88], 'ERR': [4.71, 26.16]}
SCIERRParams['abi_g16'][ (8, 'ScaledOkamoto') ]   =  {'X': [0.0, 29.5], 'ERR': [4.58, 25.65]}
SCIERRParams['abi_g16'][ (9, 'ScaledOkamoto') ]   =  {'X': [0.0, 25.89], 'ERR': [3.65, 21.35]}
SCIERRParams['abi_g16'][ (10, 'ScaledOkamoto') ]   =  {'X': [0.0, 34.01], 'ERR': [4.43, 21.59]}
SCIERRParams['abi_g16'][ (11, 'ScaledOkamoto') ]   =  {'X': [0.0, 11.52], 'ERR': [1.68, 13.42]}
SCIERRParams['abi_g16'][ (13, 'ScaledOkamoto') ]   =  {'X': [0.0, 15.97], 'ERR': [2.12, 16.01]}
SCIERRParams['abi_g16'][ (14, 'ScaledOkamoto') ]   =  {'X': [0.0, 22.32], 'ERR': [2.62, 17.68]}
SCIERRParams['abi_g16'][ (15, 'ScaledOkamoto') ]   =  {'X': [0.0, 30.89], 'ERR': [4.32, 24.42]}
SCIERRParams['abi_g16'][ (16, 'ScaledOkamoto') ]   =  {'X': [0.0, 30.02], 'ERR': [4.56, 25.96]}


class SCINormalizedError:
    def __init__(self):
        pass

    def evaluate(self,dbVals,caseParams,SCISTDName,SCI):
        # Parameterize BTerr as a ramped step function
        # --------------------------------------------
        # STD1 |. . . ____
        #      |     /.
        #      |    / .
        #      |   /  .
        # STD0 |__/   .
        #      |  .   .
        #      |__.___.___
        #         .   .
        #        SCI0 SCI1
        #---------------------------------------------
        osName = caseParams['osName']
        if osName not in SCIERRParams:
            print("ERROR: osName not available in SCIERRParams => "+osName)
            os._exit(1)

        varName, ch = vu.splitIntSuffix(caseParams['base2db'][vu.selfDepValue])
        STD0 = SCIERRParams[osName][(int(ch),SCISTDName)]['ERR'][0]
        STD1 = SCIERRParams[osName][(int(ch),SCISTDName)]['ERR'][1]
        SCI0  = SCIERRParams[osName][(int(ch),SCISTDName)]['X'][0]
        SCI1  = SCIERRParams[osName][(int(ch),SCISTDName)]['X'][1]
        slope = (STD1 - STD0) / (SCI1 - SCI0)

        belowramp = lessEqualBound(SCI,SCI0,False)
        aboveramp = greatEqualBound(SCI,SCI1,False)
        onramp    = betweenBounds(SCI,SCI0,SCI1,False)

        BTerr = np.full_like(SCI,np.NaN)
        BTerr[belowramp] = STD0
        BTerr[onramp]    = STD0 + slope * (SCI[onramp] - SCI0)
        BTerr[aboveramp] = STD1

        BTdep = dbVals[caseParams['base2db'][vu.selfDepValue]]

        return np.divide(BTdep,BTerr)

class OkamotoNormalizedError(SCINormalizedError):
    def __init__(self):
        super().__init__()
        self.baseVars = []
        self.baseVars.append(vu.selfDepValue)
        self.SCI = SCIOkamoto()
        self.baseVars = pu.uniqueMembers(self.baseVars + self.SCI.baseVars)

    def evaluate(self,dbVals,caseParams):
        SCI = self.SCI.evaluate(dbVals,caseParams)
        return super().evaluate(dbVals,caseParams,OkamotoMethod,SCI)

class ScaledOkamotoNormalizedError(SCINormalizedError):
    def __init__(self):
        super().__init__()
        self.baseVars = []
        self.baseVars.append(vu.selfDepValue)
        self.SCI = ScaledSCIOkamoto()
        self.baseVars = pu.uniqueMembers(self.baseVars + self.SCI.baseVars)

    def evaluate(self,dbVals,caseParams):
        SCI = self.SCI.evaluate(dbVals,caseParams)
        return super().evaluate(dbVals,caseParams,ScaleOkamotoMethod,SCI)

class ModHarnischNormalizedError(SCINormalizedError):
    def __init__(self):
        super().__init__()
        self.baseVars = []
        self.baseVars.append(vu.selfDepValue)
        self.SCI = SCIModHarnisch()
        self.baseVars = pu.uniqueMembers(self.baseVars + self.SCI.baseVars)

    def evaluate(self,dbVals,caseParams):
        SCI = self.SCI.evaluate(dbVals,caseParams)
        return super().evaluate(dbVals,caseParams,ModHarnischMethod,SCI)

class ScaledModHarnischNormalizedError(SCINormalizedError):
    def __init__(self):
        super().__init__()
        self.baseVars = []
        self.baseVars.append(vu.selfDepValue)
        self.SCI = ScaledSCIModHarnisch()
        self.baseVars = pu.uniqueMembers(self.baseVars + self.SCI.baseVars)

    def evaluate(self,dbVals,caseParams):
        SCI = self.SCI.evaluate(dbVals,caseParams)
        return super().evaluate(dbVals,caseParams,ScaleModHarnischMethod,SCI)

#TODO: use shapefiles/polygons to describe geographic regions instead of lat/lon boxes, e.g.,
#def outsideRegion(dbVals,REGION_NAME):
#    Note: depending on shape file definitions, LON may need to be -180 to 180 instead of 0 to 360
#
#    shp = READ_SHAPEFILE(REGION_NAME)
#    lons = dbVals['longitdue']
#    lats = dbVals['latitude']
#    nlocs = len(lons)
#    REGIONS = isinsideregions(lons,lats,shp)
#    return REGIONS


#======================
# generic bin classes
#======================
class BaseFilterFunc:
    def __init__(self, baseVars):
        self.baseVars = deepcopy(baseVars)
        pass

    def dbVars(self,varName,outerIters):
        dbVars = []
        for base in self.baseVars:
            for outerIter in outerIters:
                dbVar = vu.base2dbVar(
                    base,varName,outerIter)
                dbVars.append(dbVar)
        return pu.uniqueMembers(dbVars)


class IdFilterFunc(BaseFilterFunc):
    def __init__(self,variable):
        super().__init__([variable])
        self.result = []

    def evaluate(self,dbVals,caseParams):
        self.result = dbVals[caseParams['base2db'][self.baseVars[0]]]


class FilterFuncWrapper(BaseFilterFunc):
    def __init__(self,function):
        self.function = function()
        super().__init__(self.function.baseVars)
        self.result = []

    def evaluate(self,dbVals,caseParams):
        self.result = self.function.evaluate(dbVals,caseParams)


class BinFilter:
    def __init__(self,config,nBins,osName):
        self.where  = config['where']
        tmp         = config['bounds']

        # allow for scalar and iterable bounds
        self.bounds = np.empty(nBins)
        if (not isinstance(tmp, Iterable) or
            len(tmp) == 1):
            self.bounds[:] = tmp
        elif len(tmp) == nBins:
            for ii in list(range(nBins)):
                self.bounds[ii] = tmp[ii]
        else:
            print("ERROR: 'bounds' need to be a scalar or an Iterable with the same length as 'values'!")
            os._exit(1)

        self.apply_to_diags = config.get('apply_to_diags',vu.allDiags)
        self.mask_value     = config.get('mask_value',np.NaN)
        #TODO: add other actions besides mask_value/blacklist

        function = config.get('function',None)
        variable = config.get('variable',None)
        if variable is not None and function is None:
            # print("VARIABLE ",variable)
            self.function = IdFilterFunc(variable)
        elif function is not None and variable is None:
            # print("FUNC ",function)
            self.function = FilterFuncWrapper(function)
        else:
            print("ERROR: either 'variable' or 'function' must be provided to BinFilter constructor, but never both")
            os._exit(1)

        self.osName = osName

#    def baseVars(self):
#        return pu.uniqueMembers(self.function.baseVars)

    def dbVars(self,varName,outerIters):
        dbVars = []
        for dbVar in self.function.dbVars(
            varName,outerIters):
            dbVars.append(dbVar)
        return pu.uniqueMembers(dbVars)

    def evaluate(self,dbVals,varName,outerIter):
        caseParams = {}
        caseParams['base2db'] = {}
        for base in self.function.baseVars:
            caseParams['base2db'][base] = vu.base2dbVar(
                    base,varName,outerIter)
        caseParams['osName'] = self.osName
        self.function.evaluate(dbVals,caseParams)

    def apply(self,array,diagName,ibin):
        # blacklist locations where the mask is True
        mask = self.where(self.function.result,(self.bounds)[ibin])

        if diagName in self.apply_to_diags:
            if len(mask) == len(array):
                array[mask] = self.mask_value
            else:
                print('\n\nERROR: BinFilter mask is incorrectly defined!')
                os._exit(1)

        return array


class BinMethod:
    def __init__(self,config,osName):
        #handle scalar/str and iterable values inputs
        tmp = config['values']
        self.values = []
        if (not isinstance(tmp, Iterable) or
            isinstance(tmp,str)):
            self.values += [tmp]
        else:
            self.values += tmp

        self.filters = []
        for filterConf in config['filters']:
            self.filters.append(
                BinFilter(filterConf,len(self.values),osName) )

        enoughBounds = False
        for Filter in self.filters:
            if len(Filter.bounds) == len(self.values):
                enoughBounds = True
        if not enoughBounds:
            print('\n\nERROR: BinMethod : at least one filter must have len(bounds) == len(values)!')
            os._exit(1)

#    def baseVars(self):
#        baseVars = []
#        for Filter in self.filters:
#            for variable in Filter.baseVars():
#                baseVars.append(variable)
#        return pu.uniqueMembers(baseVars)

    def dbVars(self,varName,outerIters=['0']):
        dbVars = []
        for Filter in self.filters:
            for dbVar in Filter.dbVars(
                varName,outerIters):
                dbVars.append(dbVar)
        return pu.uniqueMembers(dbVars)

    def evaluate(self,dbVals,varName,outerIter):
        for ii in list(range(len(self.filters))):
            self.filters[ii].evaluate(
                dbVals,varName,outerIter)

    def apply(self,array,diagName,binVal):
        ibin = self.values.index(binVal)
        masked_array = deepcopy(array)
        for Filter in self.filters:
            masked_array = Filter.apply(
                masked_array,diagName,ibin)
        return masked_array


#=================================================
# binning configurations used for all bin methods
#=================================================

# Each binVarConfig member has the following properties
# key: binVar string describes the variable to bin over
# binMethod: used to distinguish between multiple methods with the same binVar
#  (e.g., defaultBinMethod, bad, latbandsMethod, etc.)
#     filters: list of filters that will blacklist locations for each method
#     for each filter:
#         where: logical function that determines locations that are blacklisted
#         variable (optional): variable that is used to initialize the IdFilterFunc class
#         function (optional): the where test is applied to function.values. defaults to None
#         bounds: numerical value(s) used in the where test (scalar or Iterable same length as values). At least one filter must have len(bounds)==len(values).
#         apply_to_diags (optional): list of diagnostics for which a particular filter applies
#     values (string): list of values associated with each bin
#
# Note: either variable or function must be provided to each filter
#
#TODO: add "diags" selection for each binMethod (some only need to be calculated for omb/oma)
#TODO: classify each binmethod as 1D, 2D, etc. to indicate which types of figures apply to it

# standard binMethod name
defaultBinMethod = 'default'

nullBinMethod = { 'filters': [], 'values': [] }
nullBinVarConfig = { defaultBinMethod:nullBinMethod }

binVarConfigs = {
    vu.obsVarQC:{
        defaultBinMethod:{
            'filters':[
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': goodFlagName,
        },
        'bad':{
            'filters':[
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': badFlags,
             'apply_to_diags': vu.nonObsDiags},
            {'where': equalBound,
             'variable': vu.selfQCValue,
             'bounds': badFlags,
             'apply_to_diags': vu.nonObsDiags,
             'mask_value': 0.0},
            ],
            'values': badFlagNames,
        },
    },
    vu.obsVarPrs:{
        defaultBinMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.prsMeta,
             'bounds': binLims[vu.obsVarPrs]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.prsMeta,
             'bounds': binLims[vu.obsVarPrs]['maxBounds']},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarPrs]['values'],
        },
        PjetMethod:{
            'filters':[
# eliminate locations outside P_jet_min to P_jet_max
            {'where': lessBound,
             'variable': vu.prsMeta,
             'bounds': P_jet_min},
            {'where': greatEqualBound,
             'variable': vu.prsMeta,
             'bounds': P_jet_max},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': P_jet_val,
        },
    },
    vu.obsVarAlt:{
        defaultBinMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.altMeta,
             'bounds': binLims[vu.obsVarAlt]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.altMeta,
             'bounds': binLims[vu.obsVarAlt]['maxBounds']},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarAlt]['values'],
        },
        altjetMethod:{
            'filters':[
# eliminate locations outside alt_jet_min to alt_jet_max
            {'where': lessBound,
             'variable': vu.altMeta,
             'bounds': alt_jet_min},
            {'where': greatEqualBound,
             'variable': vu.altMeta,
             'bounds': alt_jet_max},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': alt_jet_val,
        },
    },
    vu.obsVarLat:{
        defaultBinMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['maxBounds']},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarLat]['values'],
        },
        latbandsMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.latMeta,
             'bounds': namedLatBands['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.latMeta,
             'bounds': namedLatBands['maxBounds']},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': namedLatBands['values'],
        },
        PjetMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['maxBounds']},
            {'where': lessBound,
             'variable': vu.prsMeta,
             'bounds': P_jet_min},
            {'where': greatEqualBound,
             'variable': vu.prsMeta,
             'bounds': P_jet_max},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarLat]['values'],
        },
        altjetMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['maxBounds']},
            {'where': lessBound,
             'variable': vu.altMeta,
             'bounds': alt_jet_min},
            {'where': greatEqualBound,
             'variable': vu.altMeta,
             'bounds': alt_jet_max},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarLat]['values'],
        },
        clrskyMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['maxBounds']},
            {'where': greatEqualBound,
             'variable': vu.cldfracMeta,
             'bounds': clrskyThresh},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarLat]['values'],
        },
        cldskyMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['maxBounds']},
            {'where': lessBound,
             'variable': vu.cldfracMeta,
             'bounds': cldskyThresh},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarLat]['values'],
        },
    },
    vu.obsVarLT:{
        defaultBinMethod:{
            'filters':[
            {'where': lessBound,
             'function': LocalHour,
             'bounds': binLims[vu.obsVarLT]['minBounds']},
            {'where': greatEqualBound,
             'function': LocalHour,
             'bounds': binLims[vu.obsVarLT]['maxBounds']},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarLT]['values'],
        },
        clrskyMethod:{
            'filters':[
            {'where': lessBound,
             'function': LocalHour,
             'bounds': binLims[vu.obsVarLT]['minBounds']},
            {'where': greatEqualBound,
             'function': LocalHour,
             'bounds': binLims[vu.obsVarLT]['maxBounds']},
            {'where': greatEqualBound,
             'variable': vu.cldfracMeta,
             'bounds': clrskyThresh},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarLT]['values'],
        },
        cldskyMethod:{
            'filters':[
            {'where': lessBound,
             'function': LocalHour,
             'bounds': binLims[vu.obsVarLT]['minBounds']},
            {'where': greatEqualBound,
             'function': LocalHour,
             'bounds': binLims[vu.obsVarLT]['maxBounds']},
            {'where': lessBound,
             'variable': vu.cldfracMeta,
             'bounds': cldskyThresh},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarLT]['values'],
        },
    },
    vu.obsVarSatZen:{
        defaultBinMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.satzenMeta,
             'bounds': binLims[vu.obsVarSatZen]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.satzenMeta,
             'bounds': binLims[vu.obsVarSatZen]['maxBounds']},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarSatZen]['values'],
        },
        clrskyMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.satzenMeta,
             'bounds': binLims[vu.obsVarSatZen]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.satzenMeta,
             'bounds': binLims[vu.obsVarSatZen]['maxBounds']},
            {'where': greatEqualBound,
             'variable': vu.cldfracMeta,
             'bounds': clrskyThresh},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarSatZen]['values'],
        },
        cldskyMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.satzenMeta,
             'bounds': binLims[vu.obsVarSatZen]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.satzenMeta,
             'bounds': binLims[vu.obsVarSatZen]['maxBounds']},
            {'where': lessBound,
             'variable': vu.cldfracMeta,
             'bounds': cldskyThresh},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarSatZen]['values'],
        },
    },
    vu.obsVarCldFrac:{
        defaultBinMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.cldfracMeta,
             'bounds': binLims[vu.obsVarCldFrac]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.cldfracMeta,
             'bounds': binLims[vu.obsVarCldFrac]['maxBounds']},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarCldFrac]['values'],
        },
        clrskyMethod:{
            'filters':[
            {'where': greatEqualBound,
             'variable': vu.cldfracMeta,
             'bounds': clrskyThresh},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': clrskyMethod,
        },
        cldskyMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.cldfracMeta,
             'bounds': cldskyThresh},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': cldskyMethod,
        },
    },
    vu.obsVarSCI:{
        OkamotoMethod:{
            'filters':[
            {'where': lessBound,
             'function': SCIOkamoto,
             'bounds': binLims[vu.obsVarSCI]['minBounds']},
            {'where': greatEqualBound,
             'function': SCIOkamoto,
             'bounds': binLims[vu.obsVarSCI]['maxBounds']},
            ],
            'values': binLims[vu.obsVarSCI]['values'],
        },
        ScaleOkamotoMethod:{
            'filters':[
            {'where': lessBound,
             'function': ScaledSCIOkamoto,
             'bounds': binLims[vu.obsVarSCI]['minBounds']},
            {'where': greatEqualBound,
             'function': ScaledSCIOkamoto,
             'bounds': binLims[vu.obsVarSCI]['maxBounds']},
            ],
            'values': binLims[vu.obsVarSCI]['values'],
        },
        ModHarnischMethod:{
            'filters':[
            {'where': lessBound,
             'function': SCIModHarnisch,
             'bounds': binLims[vu.obsVarSCI]['minBounds']},
            {'where': greatEqualBound,
             'function': SCIModHarnisch,
             'bounds': binLims[vu.obsVarSCI]['maxBounds']},
            ],
            'values': binLims[vu.obsVarSCI]['values'],
        },
        ScaleModHarnischMethod:{
            'filters':[
            {'where': lessBound,
             'function': ScaledSCIModHarnisch,
             'bounds': binLims[vu.obsVarSCI]['minBounds']},
            {'where': greatEqualBound,
             'function': ScaledSCIModHarnisch,
             'bounds': binLims[vu.obsVarSCI]['maxBounds']},
            ],
            'values': binLims[vu.obsVarSCI]['values'],
        },
    },
    vu.obsVarNormErr:{
        defaultBinMethod:{
            'filters':[
            {'where': lessBound,
             'function': NormalizedError,
             'bounds': binLims[vu.obsVarNormErr]['minBounds']},
            {'where': greatEqualBound,
             'function': NormalizedError,
             'bounds': binLims[vu.obsVarNormErr]['maxBounds']},
            ],
            'values': binLims[vu.obsVarNormErr]['values'],
        },
        OkamotoMethod:{
            'filters':[
            {'where': lessBound,
             'function': OkamotoNormalizedError,
             'bounds': binLims[vu.obsVarNormErr]['minBounds']},
            {'where': greatEqualBound,
             'function': OkamotoNormalizedError,
             'bounds': binLims[vu.obsVarNormErr]['maxBounds']},
            ],
            'values': binLims[vu.obsVarNormErr]['values'],
        },
        ScaleOkamotoMethod:{
            'filters':[
            {'where': lessBound,
             'function': ScaledOkamotoNormalizedError,
             'bounds': binLims[vu.obsVarNormErr]['minBounds']},
            {'where': greatEqualBound,
             'function': ScaledOkamotoNormalizedError,
             'bounds': binLims[vu.obsVarNormErr]['maxBounds']},
            ],
            'values': binLims[vu.obsVarNormErr]['values'],
        },
        ModHarnischMethod:{
            'filters':[
            {'where': lessBound,
             'function': ModHarnischNormalizedError,
             'bounds': binLims[vu.obsVarNormErr]['minBounds']},
            {'where': greatEqualBound,
             'function': ModHarnischNormalizedError,
             'bounds': binLims[vu.obsVarNormErr]['maxBounds']},
            ],
            'values': binLims[vu.obsVarNormErr]['values'],
        },
        ScaleModHarnischMethod:{
            'filters':[
            {'where': lessBound,
             'function': ScaledModHarnischNormalizedError,
             'bounds': binLims[vu.obsVarNormErr]['minBounds']},
            {'where': greatEqualBound,
             'function': ScaledModHarnischNormalizedError,
             'bounds': binLims[vu.obsVarNormErr]['maxBounds']},
            ],
            'values': binLims[vu.obsVarNormErr]['values'],
        },
    },
    'ObsRegion':{
        'AFRICA':    nullBinMethod,
        'ATLANTIC':  nullBinMethod,
        'AUSTRALIA': nullBinMethod,
        'CONUS':{
            'filters':[
            {'where': lessBound,
             'variable': vu.lonMeta,
             'bounds': [234.0]},
            {'where': greatBound,
             'variable': vu.lonMeta,
             'bounds': [294.0]},
            {'where': lessBound,
             'variable': vu.latMeta,
             'bounds': [ 25.0]},
            {'where': greatBound,
             'variable': vu.latMeta,
             'bounds': [ 50.0]},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': ['CONUS'],
        },
        'EUROPE':    nullBinMethod,
        'E_EUROPE':  nullBinMethod,
        'NAMERICA':  nullBinMethod,
        'PACIFIC':   nullBinMethod,
        'SAMERICA':  nullBinMethod,
        'SE_ASIA':   nullBinMethod,
        'S_ASIA':    nullBinMethod,
#TODO: use shapefiles/polygons to describe geographic regions instead of lat/lon boxes, e.g.,
#        'CONUS_POLYGON':{
#            'filters':[
#            {'where': notInBound,
#             'function': Regions,
#             'variable': [vu.lonMeta,vu.latMeta],
#             'bounds': ['CONUS']},
#            ],
#            'values': ['CONUS'],
#        },
#        'EUROPE_POLYGON':{
#            'filters':[
#            {'where': notInBound,
#             'function': Regions,
#             'variable': [vu.lonMeta,vu.latMeta],
#             'bounds': ['EUROPE']},
#            ],
#            'values': ['EUROPE'],
#        },
    },
#TODO: enable binning in MPAS Model space
#    vu.modVarPrs:{
#        defaultBinMethod:{
#            'filters':[
#            {'where': lessBound,
#             'variable': vu.modVarPrs,
#             'bounds': binLims[vu.modVarPrs]['minBounds']},
#            {'where': greatEqualBound,
#             'variable': vu.modVarPrs,
#             'bounds': binLims[vu.modVarPrs]['maxBounds']},
#            ],
#            'values': binLims[vu.modVarPrs]['values'],
#        },
#    },
#    vu.modVarAlt:{
#        defaultBinMethod:{
#            'filters':[
#            {'where': lessBound,
#             'variable': vu.modVarAlt,
#             'bounds': binLims[vu.modVarAlt]['minBounds']},
#            {'where': greatEqualBound,
#             'variable': vu.modVarAlt,
#             'bounds': binLims[vu.modVarAlt]['maxBounds']},
#            ],
#            'values': binLims[vu.modVarAlt]['values'],
#        },
#    },
#    vu.modVarLat:{
#        defaultBinMethod:{
#            'filters':[
#            {'where': lessBound,
#             'variable': vu.modVarLat,
#             'bounds': binLims[vu.modVarLat]['minBounds']},
#            {'where': greatEqualBound,
#             'variable': vu.modVarLat,
#             'bounds': binLims[vu.modVarLat]['maxBounds']},
#            ],
#            'values: binLims[vu.modVarLat]['values'],
#        },
#        latbandsMethod:{
#            'filters':[
#            {'where': lessBound,
#             'variable': vu.modVarLat,
#             'bounds': namedLatBands['minBounds']},
#            {'where': greatEqualBound,
#             'variable': vu.modVarLat,
#             'bounds': namedLatBands['maxBounds']},
#            ],
#            'values': namedLatBands['values'],
#        },
#    },
#    'ModelRegion':{
#        'AFRICA':    nullBinMethod,
#        'ATLANTIC':  nullBinMethod,
#        'AUSTRALIA': nullBinMethod,
#        'CONUS':{
#            'filters':[
#            {'where': lessBound,
#             'variable': vu.modVarLon,
#             'bounds': [234.0]},
#            {'where': greatBound,
#             'variable': vu.modVarLon,
#             'bounds': [294.0]},
#            {'where': lessBound,
#             'variable': vu.modVarLat,
#             'bounds': [ 25.0]},
#            {'where': greatBound,
#             'variable': vu.modVarLat,
#             'bounds': [ 50.0]},
#            ],
#            'values': ['CONUS'],
#        },
#        'EUROPE':    nullBinMethod,
#        'E_EUROPE':  nullBinMethod,
#        'NAMERICA':  nullBinMethod,
#        'PACIFIC':   nullBinMethod,
#        'SAMERICA':  nullBinMethod,
#        'SE_ASIA':   nullBinMethod,
#        'S_ASIA':    nullBinMethod,
#    },
}


#=============================
# Parameterized binVarConfigs
#=============================

# Add named latitude-band-specific pressure bins
for iband, latBand in enumerate(namedLatBands['values']):
    binVarConfigs[vu.obsVarPrs][latBand] = {
        'filters':[
        {'where': lessBound,
         'variable': vu.prsMeta,
         'bounds': binLims[vu.obsVarPrs]['minBounds']},
        {'where': greatEqualBound,
         'variable': vu.prsMeta,
         'bounds': binLims[vu.obsVarPrs]['maxBounds']},
        {'where': lessBound,
         'variable': vu.latMeta,
         'bounds': namedLatBands['minBounds'][iband]},
        {'where': greatEqualBound,
         'variable': vu.latMeta,
         'bounds': namedLatBands['maxBounds'][iband]},
        {'where': notEqualBound,
         'variable': vu.selfQCValue,
         'bounds': goodFlag,
         'apply_to_diags': vu.nonObsDiags},
        ],
        'values': binLims[vu.obsVarPrs]['values'],
    }


# Add named latitude-band-specific altitude bins
for iband, latBand in enumerate(namedLatBands['values']):
    binVarConfigs[vu.obsVarAlt][latBand] = {
        'filters':[
        {'where': lessBound,
         'variable': vu.altMeta,
         'bounds': binLims[vu.obsVarAlt]['minBounds']},
        {'where': greatEqualBound,
         'variable': vu.altMeta,
         'bounds': binLims[vu.obsVarAlt]['maxBounds']},
        {'where': lessBound,
         'variable': vu.latMeta,
         'bounds': namedLatBands['minBounds'][iband]},
        {'where': greatEqualBound,
         'variable': vu.latMeta,
         'bounds': namedLatBands['maxBounds'][iband]},
        {'where': notEqualBound,
         'variable': vu.selfQCValue,
         'bounds': goodFlag,
         'apply_to_diags': vu.nonObsDiags},
        ],
        'values': binLims[vu.obsVarAlt]['values'],
    }

def main():
    print ('This is not a runnable program.')
    os._exit(0)

if __name__ == '__main__': main()
