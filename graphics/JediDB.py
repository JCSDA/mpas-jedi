#!/usr/bin/env python3

from collections import defaultdict
import config as conf
from copy import deepcopy
import glob
from JediDBArgs import \
  obsFKey, geoFKey, diagFKey, default_path, fPrefixes
import h5py as h5
import logging
import multiprocessing as mp
from netCDF4 import Dataset
import numpy as np
import os
import plot_utils as pu
import var_utils as vu
import sys

_logger = logging.getLogger(__name__)

MAXINT32  = np.int32(1e9)
MAXFLOAT  = np.float32(1.e12)
MAXDOUBLE = np.float64(1.e12)

def getIODAFileRank(pathPlusFile):
    # pathPlusFile = PATH/fileName
    # fileName = prefix_suffix.fileExt
    # suffix may or may not be a numerical value indicating processing element rank
    # return the rank if any

    # fileName excludes the path
    fileName = pathPlusFile.split('/')[-1]

    # fileParts as list, removing extension (after '.') first
    fileParts = ('.'.join(fileName.split('.')[:-1])).split('_')

    rank = fileParts[-1]

    # check if suffix is numerical, indicating the rank number
    if not pu.isint(rank):
        return 'none'
    else:
        return rank


def IODAFileIsRanked(pathPlusFile):
    return pu.isint(getIODAFileRank(pathPlusFile))


def getIODAFileHandle(File, mode='r'):
  # returns a FileHandle object for File with permissions declared in mode
  ext = File.split('.')[-1]
  isNC = (ext in 'nc4')
  isHDF = (ext == 'h5')
  if isNC:
    return NCFileHandle(File, mode)
  elif isHDF:
    return HDFFileHandle(File, mode)
  else:
    _logger.error('getIODAFileHandle:unsupported file extension => '+ext)


class FileHandles():
  '''wrapper for list of FileHandle objects'''
  def __init__(self, files, mode):
    self.handles = []
    for File in files:
      self.handles.append(getIODAFileHandle(File, mode))

    self.size = len(self.handles)

    self.nlocs = 0
    for h in self.handles:
        self.nlocs += h.nlocs

    self.variables = self.handles[0].variables
    self.fileFormat = self.handles[0].fileFormat

  def close(self):
    for h in self.handles:
      h.close()

  def varsINGroup(self, group):
    # comprehensive list of variables in group
    filevars = self.handles[0].varsINGroup(group)

    # only include variables with finite data
    # ensures only assimilated channels are selected
    varlist = []
    for var in filevars:
        grpVar = vu.AllVarCtors[self.fileFormat](var, group)
        var1D = self.singleVarAtLocations(grpVar)
        if np.isfinite(var1D).sum() > 0: varlist.append(var)

        # TODO(JJG-05OCT2022): turn off above check for gnssrobndnbam to get through due to all missing values
        #varlist.append(var)

    return varlist

  def datatype(self, var):
    return self.handles[0].datatype(var)

  def singleVarAtLocations(self, grpVar):
    varName, group = vu.splitObsVarGrp(grpVar)

    # GeoVaLs and ObsDiagnostics do not have ObsGroups in the file variable names
    if group==vu.geoGroup or group==vu.diagGroup:
      fileVar = varName
    else:
      fileVar = grpVar

    dtype = self.datatype(fileVar)

    if fileVar not in self.variables:
      _logger.error('FileHandles.singleVarAtLocations: fileVar not present: '+fileVar)

    if 'byte' in dtype:
      valsOut = np.empty(self.nlocs, dtype=np.object_)
      istart = 0
      for h in self.handles:
        iend = istart + h.nlocs
        tmp = np.empty(h.nlocs, dtype=np.object_)
        for ii, bytelist in enumerate(h.singleVarAtLocations(fileVar)):
          tmp[ii] = (b''.join(bytelist)).decode('utf-8')
        valsOut[istart:iend] = tmp
        istart = iend
    else:
      # determine number of levels
      h = self.handles[0]
      field0 = np.asarray(h.singleVarAtLocations(fileVar))
      shape = field0.shape

      istart, iend = 0, -1
      if len(shape) > 1 and shape[1] > 1:
        # treat level-dependent variables
        valsOut = np.empty((self.nlocs, shape[1]), dtype=dtype)
        for h in self.handles:
          iend = istart + h.nlocs
          valsOut[istart:iend,:] = h.singleVarAtLocations(fileVar)
          istart = iend

      elif len(shape) < 3:
        # treat single-level/surface variables
        valsOut = np.empty(self.nlocs, dtype=dtype)
        for h in self.handles:
          iend = istart + h.nlocs
          if len(shape) == 1:
            valsOut[istart:iend] = h.singleVarAtLocations(fileVar)
          elif shape[1] == 1:
            # GeoVaLs and ObsDiagnostics are always stored as 2D, even for a single level
            valsOut[istart:iend] = np.transpose(np.asarray(h.singleVarAtLocations(fileVar)))[0][:]
          istart = iend

      else:
        _logger.error('FileHandles.singleVarAtLocations: unable to handle more than 2 dimensions')

    assert iend == self.nlocs, (
      'FileHandles.singleVarAtLocations: incorrect nlocs (',iend,'!=',nlocs,') for ', grpVar)

    # missing data
    if 'int32' in dtype:
      missing = np.greater(np.abs(valsOut), MAXINT32)
    elif 'float32' in dtype:
      missing = np.greater(np.abs(valsOut), MAXFLOAT)
    elif 'float64' in dtype:
      missing = np.greater(np.abs(valsOut), MAXDOUBLE)
    else:
      missing = np.full_like(valsOut, False, dtype=bool)

    if 'float' in dtype:
      valsOut[missing] = np.NaN

    #TODO: missing value handling for integers, strings, and others?

    # convert pressure from Pa to hPa if needed (kludge)
    if (vu.obsVarPrs in grpVar and np.max(valsOut) > 10000.0):
      finite = np.isfinite(valsOut)
      valsOut[finite] = np.divide(valsOut[finite], 100.0)

    # force longitude values into the range 0 to 360 (affects gnssro from ncdiag)
    if vu.obsVarLon in grpVar:
      finite = np.isfinite(valsOut)

      greater = np.full_like(finite, False, bool)
      greater[finite] = np.greater(valsOut[finite], 360.0)
      valsOut[greater] = np.subtract(valsOut[greater], 360.0)

      less = np.full_like(finite, False, bool)
      less[finite] = np.less(valsOut[finite], 0.0)
      valsOut[less] = np.add(valsOut[less], 360.0)

    return valsOut


class FileHandle():
  '''This serves as a base class for netcdf4 and hdf5 wrapper classes for IODA-formatted and UFO GeoVaLs-formatted files'''
  def __init__(self, File, mode):
    self.variables = {}
    self.File = File
    self.mode = mode

  def datatype(self, var):
    '''
    virtual method, returns the basic data type of var
    '''
    raise NotImplementedError()

  def singleVarAtLocations(self, var):
    '''
    virtual method, returns var as a 1D array
    '''
    raise NotImplementedError()

  def close(self):
    del self.variables

  def varsINGroup(self, group):
    # returns all variables in an ObsGroup with the ObsGroup label removed
    assert group!=vu.geoGroup and group!=vu.diagGroup, (
      'varsINGroup cannot handle GeoVaLs or ObsDiagnostics')
    varlist = []
    for vargrp in self.variables.keys():
        var, grp = vu.splitObsVarGrp(vargrp)
        if grp == group: varlist.append(var)
    return varlist


# NETCDF4
# ObsSpace IODA-v1, GeoVaLS, ObsDiagnostics
class NCFileHandle(FileHandle):
  def __init__(self, File, mode):
    super().__init__(File, mode)

    self.fileFormat = vu.ncFileFormat
    self.h = Dataset(File, mode)
    self.h.set_auto_mask(False)
    self.nlocs = self.h.dimensions['nlocs'].size

    self.variables = self.h.variables

  def datatype(self, var):
    return self.h.variables[var].datatype.name

  def singleVarAtLocations(self, var):
    return self.h.variables[var][:]

  def close(self):
    super().close()
    self.h.close()


# HDF5, including 2D variables
# ObsSpace IODA-v2
class HDF2DtoArray1D():
  def __init__(self, h, var2D, index):
    self.h = h
    self.var2D = var2D
    self.index = index

  def get(self):
    return self.h[self.var2D][:,self.index]

class HDFFileHandle(FileHandle):
  def __init__(self, File, mode):
    super().__init__(File, mode)

    self.fileFormat = vu.hdfFileFormat
    self.h = h5.File(File, mode)
    self.nlocs = self.h['Location'].size

    varlist = []
    for node in self.h:
      if type(self.h[node]) is h5._hl.group.Group:
        for var in self.h[node]:
          varlist += [node+'/'+var]

    self.variables = {}
    for var in varlist:
      shape = self.h[var].shape

      # retrieve dimensions
      dims = np.empty(len(shape), object)
      for ii, dim in enumerate(self.h[var].attrs['DIMENSION_LIST']):
        dims[ii] = self.h[dim[0]]
        assert len(dims[ii]) == shape[ii], ('HDFFileHandle.init: len(dims[ii]) and shape[ii] do not match, ii = ', ii)
      if dims[0] != self.h['Location']: continue

      if len(shape) == 1:
        self.variables[var] = self.h[var]
      elif len(shape) == 2:
        # assign proxy HDF2DtoArray1D objects along second dimension
        for ii, d2Value in enumerate(dims[1]):
          self.variables[vu.appendSuffix(var, d2Value)] = HDF2DtoArray1D(self.h, var, ii)
      else:
        _logger.error('HDFFileHandle.init: unable to handle more than 2 dimensions')

  def datatype(self, var):
    # first branch of logic accounts for the fact that some variables are stored in hdf5 files with a channel suffix
    if var in self.h:
      fileVarName = var
    # second branch of logic removes the channel suffix that was artifically added in self.variables
    else:
      varName, grp = vu.splitObsVarGrp(var)
      fileVarName, suf = vu.splitIntSuffix(varName)
      fileVarName = grp+'/'+fileVarName
    return self.h[fileVarName].dtype.name

  def singleVarAtLocations(self, var):
    if type(self.variables[var]) is h5._hl.dataset.Dataset:
      return self.variables[var][:]
    elif type(self.variables[var]) is HDF2DtoArray1D:
      return self.variables[var].get()
    else:
      _logger.error('HDFFileHandle.singleVarAtLocations: unsupported type')

  def close(self):
    super().close()
    self.h.close()


class JediDB:
    '''This class provides access to UFO feedback files.'''
    def __init__(self, data_path=default_path, osKeySelect=[]):
        # data_path: location of feedback files
        # fileExt: extention of feedback files
        # osKeySelect (optional): allows for the selection of particular osKey's

        supportedFileExts=['nc4','h5']

        self.filePrefixes = {}
        for key, prefix in fPrefixes.items():
            self.filePrefixes[key] = prefix

        # Group all files by experiment-ObsSpace combination
        # Two possible formats for multi-processor or single-processor/concatenated
        #  multi-processor: self.filePrefixes[fileType]+"_"+osKey+"_"+rank+"."+fileExt
        #  single-processor: self.filePrefixes[fileType]+"_"+osKey+"."+fileExt

        #  osKey includes strings associated with both the experiment and ObsSpaceName
        #  rank, when present, is the four digit processor element rank.  If a file matching
        #  the single-processor format is present, it will be selected in place of the other
        #  files.
        #
        #  examples:
        #   +  obsout_hofx_sondes_0000.nc4 OR
        #   +  obsout_hofx_sondes.nc4
        # osKey = "hofx_sondes"
        #
        #   +  obsout_variational_amsua_n19_0000.nc4 OR
        #   +  obsout_variational_amsua_n19.nc4
        # osKey = "variational_amsua_n19"
        #
        #   +  obsout_3dvar_bumpcov_bumploc_amsua_n19_0000.nc4 OR
        #   +  obsout_3dvar_bumpcov_bumploc_amsua_n19.nc4
        # osKey = "3dvar_bumpcov_bumploc_amsua_n19"
        #
        #   +  obsout_abi_g16_0000.nc4 OR
        #   +  obsout_abi_g16.nc4
        # osKey = "abi_g16"

        # note: when obsFKey (ObsSpace) files are concatenated into a single file, the
        #   optional geoval and ydiag files must also be concatenated in order to
        #   preserve locations ordering

        allFiles = {}
        for fileType, prefix in self.filePrefixes.items():
            allFiles[fileType] = defaultdict(list)
            for fileExt in supportedFileExts:
                for pathPlusFile in glob.glob(data_path+'/'+prefix+'*.'+fileExt):
                    # fileName excludes the path
                    fileName = pathPlusFile.split('/')[-1]

                    # fileParts as list, removing extension after final '.'
                    fileParts = ('.'.join(fileName.split('.')[:-1])).split('_')

                    # remove known prefix
                    fileParts.remove(prefix)

                    if IODAFileIsRanked(pathPlusFile):
                        # osKey excludes PE after the final '_'
                        osKey =  '_'.join(fileParts[:-1])
                    else:
                        # osKey includes all remaining parts
                        osKey = '_'.join(fileParts)

                    allFiles[fileType][osKey].append(pathPlusFile)

        # + transpose files for easy access by osKey
        # + remove PE-ranked files when concatenated file is present
        self.Files = {}
        for fileType, osKeys in allFiles.items():
            for osKey, files in osKeys.items():
                if osKey not in self.Files:
                  self.Files[osKey] = defaultdict(list)
                for pathPlusFile in files:
                    if not IODAFileIsRanked(pathPlusFile):
                        # for non-numerical suffix assume that this one file includes ALL locations
                        self.Files[osKey][fileType] = [pathPlusFile]
                        break
                    else:
                        # otherwise, append this file to list
                        self.Files[osKey][fileType].append(pathPlusFile)
        del allFiles

        # error checking
        self.loggers = {}
        for osKey, fileTypes in self.Files.items():
            self.loggers[osKey] = logging.getLogger(__name__+'.'+osKey)
            nObsFiles = len(fileTypes.get(obsFKey,[]))
            # force crash when there are no obsFKey files available for an osKey
            if nObsFiles < 1:
                self.loggers[osKey].error('''
                           There are no '''+obsFKey+'''
                           feedback files with prefix '''+self.filePrefixes[obsFKey]+'''
                           osKey = '''+osKey)

        # eliminate osKeys that are designated to not be processed
        self.ObsSpaceName = {}
        self.ObsSpaceGroup = {}
        for osKey in list(self.Files.keys()):
            # determine ObsSpace from osKey string
            # and only process valid ObsSpaceNames
            expt_parts = osKey.split("_")
            nstr = len(expt_parts)
            ObsSpaceInfo = conf.nullDiagSpaceInfo
            for i in range(0,nstr):
                ObsSpaceName_ = '_'.join(expt_parts[i:nstr])
                ObsSpaceInfo_ = conf.DiagSpaceConfig.get( ObsSpaceName_,conf.nullDiagSpaceInfo)
                if ObsSpaceInfo_['process']:
                    ObsSpaceName = deepcopy(ObsSpaceName_)
                    ObsSpaceInfo = deepcopy(ObsSpaceInfo_)
            if ((len(osKeySelect)>0 and osKey not in osKeySelect) or
                not ObsSpaceInfo['process'] or
                ObsSpaceInfo.get('DiagSpaceGrp',conf.model_s) == conf.model_s):
                del self.Files[osKey]
            else:
                self.ObsSpaceName[osKey] = deepcopy(ObsSpaceName)
                self.ObsSpaceGroup[osKey] = deepcopy(ObsSpaceInfo['DiagSpaceGrp'])

        # sort remaining files by rank
        for osKey, fileTypes in self.Files.items():
            for fileType, files in fileTypes.items():
                ranks = []
                for fileName in files:
                    rank = getIODAFileRank(fileName)
                    if pu.isint(rank):
                        ranks.append(int(rank))
                    elif len(files) == 1:
                        ranks = [1]
                        break
                    else:
                        self.loggers[osKey].error('too many files chosen for a concatenated ObsSpace=> '+fileType)

                indices = list(range(len(files)))
                indices.sort(key=ranks.__getitem__)
                self.Files[osKey][fileType] = \
                    list(map(files.__getitem__, indices))

        self.Handles = {}

    # TODO: add function that concatenates files together, then
    # TODO: replaces the original Files and/or Handles

    def initHandles(self, osKey):
    #initialize Handles for osKey
    # osKey (string) - experiment-ObsSpace key
        self.loggers[osKey].info('Initializing UFO file handles...')

        self.Handles[osKey] = {}

        for fileType, files in self.Files[osKey].items():
            self.loggers[osKey].info(' fileType = '+fileType)
            self.Handles[osKey][fileType] = FileHandles(files, 'r')


    def destroyHandles(self, osKey):
    #destroys file handles for osKey
    # osKey (string) - experiment-ObsSpace key
        if osKey in self.Handles:
            for fileType, h in self.Handles[osKey].items():
                h.close()
            del self.Handles[osKey]


    def fileFormat(self, osKey, fileType):
    #return the file format for the first file handle
    # osKey (string)     - experiment-ObsSpace key
    # fileType (string)  - file type key for self.Files[osKey]

        # extract file handles
        fHandles = self.Handles[osKey].get(fileType, None)
        if fHandles is None:
            self.loggers[osKey].error('no files exist => '+fileType)

        return fHandles.fileFormat


    def varList(self, osKey, fileType, selectGrp):
    #return a list of variable names that fits the input parameters
    # osKey (string)     - experiment-ObsSpace key
    # selectGrp (string) - variable group name desired from self.Files
    # fileType (string)  - file type key for self.Files[osKey]

        # extract file handles
        fHandles = self.Handles[osKey].get(fileType, None)
        if fHandles is None:
            self.loggers[osKey].error('no files exist => '+fileType)

        assert fileType == obsFKey, 'varList not implemented for '+fileType

        varlist = fHandles.varsINGroup(selectGrp)

        # sort varlist alphabetically or
        # by integer suffix for uniform dictName (e.g., channel number)
        indices = list(range(len(varlist)))
        dictName0, suf = vu.splitIntSuffix(varlist[0])
        intlist = []
        sortbyInt = True
        for var in varlist:
            dictName, suf = vu.splitIntSuffix(var)
            # do not sort by integer suffix when
            # + there is any non-integer suffix
            # + there is more than one unique dictName
            if not pu.isint(suf) or dictName != dictName0:
                sortbyInt = False
                break
            intlist.append(suf)

        if sortbyInt:
            indices.sort(key=[int(i) for i in intlist].__getitem__)
        else:
            indices.sort(key=varlist.__getitem__)
        varlist = list(map(varlist.__getitem__, indices))

        return varlist


    def readVars(self, osKey, dbVars, np=1):
    # osKey (string)  - experiment-ObsSpace key
    # dbVars (string) - list of variables desired from self.Handles[osKey]

        self.loggers[osKey].info('Reading requested variables from UFO file(s)...')

        ObsSpace = self.Handles[osKey][obsFKey]
        GeoVaLs = self.Handles[osKey].get(geoFKey, None)
        ObsDiagnostics = self.Handles[osKey].get(diagFKey, None)

        #TODO: loop over grpVar with workers when workers is not None
        nprocs = min(mp.cpu_count(), np)

        # for now, only use one pe
        nprocs = 1

        if nprocs > 1:
            workers = mp.Pool(processes = nprocs)
        else:
            workers = None

        # Construct output dictionary
        varsVals = {}
        varsAsync = {}
        for grpVar in pu.uniqueMembers(dbVars):
            if workers is None:
                varsVals[grpVar] = self.readVar(grpVar, ObsSpace, GeoVaLs, ObsDiagnostics)
            else:
                varsAsync[grpVar] = workers.apply_async(self.readVar,
                    args = (grpVar, ObsSpace, GeoVaLs, ObsDiagnostics))

        if workers is not None:
            workers.close()
            workers.join()
            for grpVar, v in varsAsync.items():
                print(grpVar, v)
                varsVals[grpVar] = v.get()

        return varsVals


    def readVar(self, grpVar, o, g, d):

        varName, grpName = vu.splitObsVarGrp(grpVar)

        if grpVar in o.variables:
            values = o.singleVarAtLocations(grpVar)

        elif vu.geoGroup in grpName and g is not None:
            values = g.singleVarAtLocations(grpVar)

        elif vu.diagGroup in grpName and d is not None:
            values = d.singleVarAtLocations(grpVar)

        else:
            self.loggers[osKey].error('grpVar not found => '+grpVar)

        return values
