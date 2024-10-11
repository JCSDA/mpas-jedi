import numpy as np
from datetime import datetime, timedelta
import sys, itertools, os, collections
import logging
import xarray as xr
import warnings
import multiprocessing
import argparse

float_missingVal = float(-999)
int_missingVal = int(-999)
IODAfloat_missing_value = float(-3.3687953e+38)
IODAint_missing_value   = int(-2147483643)
hydro_initial = float(5.4999992258578274e-15)
ioda_float_type = 'float32'
ioda_int_type = 'int32'

criteria_list_2by2 = ['Na_Ocloud-Mcloud','Nb_Oclear-Mcloud','Nc_Ocloud-Mclear','Nd_Oclear-Mclear','all','Ocloud','Oclear']
criteria_3p = ['Na_Ocloud-Mcloud','Nd_Oclear-Mclear', 'all']
criteria_list_4by4 = ['a_Olow-Mlow', 'b_Omid-Mlow', 'c_Ohigh-Mlow', 'd_Oclear-Mlow', 'e_Olow-Mmid', 'f_Omid-Mmid', 'g_Ohigh-Mmid', 'h_Oclear-Mmid', 'i_Olow-Mhigh', 'j_Omid-Mhigh', 'k_Ohigh-Mhigh', 'l_Oclear-Mhigh', 'm_Olow-Mclear', 'n_Omid-Mclear', 'o_Ohigh-Mclear', 'p_Oclear-Mclear', 'all', 'Ocloud', 'Oclear','lowOnly', 'midOnly', 'highOnly', 'Olow', 'Omid', 'Ohigh', 'Na_Ocloud-Mcloud','Nb_Oclear-Mcloud','Nc_Ocloud-Mclear','Nd_Oclear-Mclear']
criteria_only = ['lowOnly', 'midOnly', 'highOnly']

Na = ['a_Olow-Mlow', 'b_Omid-Mlow', 'c_Ohigh-Mlow', 'e_Olow-Mmid', 'f_Omid-Mmid', 'g_Ohigh-Mmid', 'i_Olow-Mhigh', 'j_Omid-Mhigh', 'k_Ohigh-Mhigh']
Nb = ['d_Oclear-Mlow', 'h_Oclear-Mmid', 'l_Oclear-Mhigh']
Nc = ['m_Olow-Mclear', 'n_Omid-Mclear', 'o_Ohigh-Mclear']
Nd = ['p_Oclear-Mclear']

cat_stats = ['FBIAS','POD','FAR','GSS','CSI','PODF','PODN','HRATE','ACC']
cat_stats_reduced = cat_stats[:5]
cont_stats = ['mean', 'count', 'rmse', 'std']

layering_list = ['all','LowClouds','MidClouds','HighClouds']

binningDict = {'CloudFraction': criteria_list_2by2,
               'CTPLayering': criteria_list_4by4}

# UPP top layer bounds (Pa) for cloud layers
PTOP_LOW_UPP  = 64200. # low for > 64200 Pa
PTOP_MID_UPP  = 35000. # mid between 35000-64200 Pa
PTOP_HIGH_UPP = 15000. # high between 15000-35000 Pa

varsList_  = ['obs', 'bkg', 'lat', 'lon', 'ocldfrac', 'fcldfrac', 'octp']

hydrosList = ['ihumr', 'imclw', 'imcli', 'imclr', 'imcls', 'imclg']
hydroDict = {'ihumr': 'water_vapor_mixing_ratio_wrt_dry_air',
             'imclw': 'mass_content_of_cloud_liquid_water_in_atmosphere_layer',
             'imcli': 'mass_content_of_cloud_ice_in_atmosphere_layer',
             'imclr': 'mass_content_of_rain_in_atmosphere_layer',
             'imcls': 'mass_content_of_snow_in_atmosphere_layer',
             'imclg': 'mass_content_of_graupel_in_atmosphere_layer'}

colorsL = ['k','r','b','g','c','m','y']

############################# FUNCTIONS #######################################

def str2bool(v):
    # https://stackoverflow.com/a/43357954/1264304
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def list_of_strings(arg):
    return arg.split(',')

def available_processors():
    host = os.environ['NCAR_HOST']
    if host == 'derecho':
      select = os.environ['PBS_SELECT']
      nnodes = int(select.split(':')[0])
      ncpus = int(select.split(':')[1].split('=')[1])
      memm = int(select.split(':')[3].split('=')[1][:-2])
      threads = int(select.split(':')[4].split('=')[1])
      nprocs = nnodes * threads * ncpus
    else:
      nprocs = multiprocessing.cpu_count()
    return nprocs

def logg():
    logFile = 'log.verification'
    try:
        os.remove(logFile)
    except OSError:
        pass
    return logging.basicConfig(filename=logFile, filemode='w', format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.INFO)

def get_list4ploting(binningList):
    if len(binningList) >= 16:
      criteList = criteria_list_4by4[:16]
      statLayer = list(itertools.product(cat_stats_reduced,layering_list))
    else:
      criteList = criteria_list_2by2[:-3]
      statLayer = list(itertools.product(cat_stats,layering_list))
    return criteList, statLayer

def get_label(data,stat,criteria=None):
    if stat == 'count':
      if criteria is not None:
        label = ' ('+str(int(np.nansum(data[criteria])))+')'
    elif stat == 'mean':
      if criteria is not None:
        label = ' ('+str(round(nanmean(data[criteria]),2))+')'
    else:
      label = ''
    return label

def varDict_binned(List,varsList):
    fillValue = np.nan
    return {var: {k: fillValue for k in List} for var in varsList}

def varDict3L(List1,List2,List3,varsList,fillValue=None):
    if fillValue == 'NaN':
      return {var: {l1: {l2: {l3: np.nan for l3 in List3} for l2 in List2} for l1 in List1} for var in varsList}
    elif fillValue == 0:
      return {var: {l1: {l2: {l3: 0 for l3 in List3} for l2 in List2} for l1 in List1} for var in varsList}
    else:
      return {var: {l1: {l2: {l3: np.array([],dtype=ioda_float_type) for l3 in List3} for l2 in List2} for l1 in List1} for var in varsList}

def varContStatDict(List1,List2,fillValue=None):
    if fillValue == 'NaN':
      return {stat: {l1: {l2: np.nan for l2 in List2} for l1 in List1} for stat in cont_stats}
    elif fillValue == 0:
      return {stat: {l1: {l2: 0 for l2 in List2} for l1 in List1} for stat in cont_stats}
    else:
      return {stat: {l1: {l2: np.array([],dtype=ioda_float_type) for l2 in List2} for l1 in List1} for stat in cont_stats}

def varContHydroStatDict(List1,List2,fillValue=None):
    if fillValue == 'NaN':
      return {var: {stat: {l1: {l2: np.nan for l2 in List2} for l1 in List1} for stat in cont_stats} for var in hydrosList}
    elif fillValue == 0:
      return {var: {stat: {l1: {l2: 0 for l2 in List2} for l1 in List1} for stat in cont_stats} for var in hydrosList}
    else:
      return {var: {stat: {l1: {l2: np.array([],dtype=ioda_float_type) for l2 in List2} for l1 in List1} for stat in cont_stats} for var in hydrosList}

def createNestedDictT(List1,List2,List3,fillValue=None):
    if fillValue == 'NaN':
      return {l1: {l2: {l3: np.nan for l3 in List3} for l2 in List2} for l1 in List1}
    elif fillValue == 0:
      return {l1: {l2: {l3: 0 for l3 in List3} for l2 in List2} for l1 in List1}
    else:
      return {l1: {l2: {l3: np.array([],dtype=ioda_int_type) for l3 in List3} for l2 in List2} for l1 in List1}

def createDictL(List,fillList):
    return {List[k]: fillList[k] for k in range(len(List))}

def m2tokm2(data):
    return data*1e-6

def get_meshSpec(meshRes):
    if meshRes == '120km':
      nCellsVMesh = 40962
    elif meshRes == '60km':
      nCellsVMesh = 163842
    elif meshRes == '30km':
      nCellsVMesh = 655362
    elif meshRes == '60-3km':
      nCellsVMesh = 835586
    mr = [int(i) for i in meshRes.split('km')[0].split('-')]
    if len(mr) > 1: # variable mesh
      meshRatio = int(mr[0]/mr[1])
    else:
      meshRatio = 1 # uniform mesh

    return nCellsVMesh, meshRatio

def open_dataset(filename, enginetype, group=None):
    if enginetype == 'h5':
      engine = 'h5netcdf'
      return xr.open_dataset(filename, group=group, engine=engine)
    elif enginetype == 'nc':
      engine = 'netcdf4'
      return xr.open_dataset(filename, engine=engine)
    else:
      raise NotImplementedError()

def getDomainIndexes(threshold,domain,metadata,staticnc):
    ocellindex = np.array(metadata.cellIndex.data-int(1),dtype=ioda_int_type) # ocellindex starts in 1 (from fortran), here it has to start from 0 (python)
    olatitude  = np.array(metadata.latitude.data,dtype=ioda_float_type)
    olongitude = np.array(metadata.longitude.data,dtype=ioda_float_type) # (longitude + 180) % 360 - 180
    areaCell = np.array(staticnc.variables['areaCell'][:],dtype=ioda_float_type) # spherical area of a Voroni cell
    areakm2 = m2tokm2(areaCell)
    matches = ocellindex; obs_ind = None

    if domain != 'fullDisk':
      if domain == 'finer':
        cellindex_ids = np.where( areakm2 <= threshold )[0]
      elif domain == 'coarser':
        cellindex_ids = np.where( areakm2 > threshold )[0]
      else:
        raise NotImplementedError()
      matches, obs_ind, model_ind = np.intersect1d(ocellindex, cellindex_ids, return_indices=True)

    return matches, obs_ind

def get_obsvalue(obsvalue,ch,ids):
    obt        = np.array(obsvalue.brightnessTemperature.isel(Channel=ch).data,dtype=ioda_float_type)
    obt[np.ma.where(obt == 0.0)] = np.nan
    if ids is not None:
      return obt[ids]
    else:
      return obt

def get_hofx(hofx,ch,ids):
    fbt        = np.array(hofx.brightnessTemperature.isel(Channel=ch).data,dtype=ioda_float_type)
    fbt[np.ma.where(fbt == 0.0)] = np.nan
    if ids is not None:
      return fbt[ids]
    else:
      return fbt

def get_cldfrac(metadata,ids):
    ocldfrac   = np.array(metadata.cloudAmount.data,dtype=ioda_float_type) *100 #Unitless to percent
    if ids is not None:
      return ocldfrac[ids]
    else:
      return ocldfrac

def get_latlon(metadata,ids):
    latitude  = np.array(metadata.latitude.data,dtype=ioda_float_type)
    longitude = np.array(metadata.longitude.data,dtype=ioda_float_type)
    #longitude = (longitude + 180) % 360 - 180
    if ids is not None:
      return latitude[ids], longitude[ids]
    else:
      return latitude, longitude

def get_geoval(geovalnc,ids):
    fcldfrac   = np.array(geovalnc.cloud_area_fraction_in_atmosphere_layer.data,dtype=ioda_float_type) *100.0 #Unitless to percent
    pressure   = np.array(geovalnc.air_pressure.data,dtype=ioda_float_type)
    psfc       = np.array(geovalnc.air_pressure_levels.data[:,-1],dtype=ioda_float_type)
    if ids is not None:
      return fcldfrac[ids], pressure[ids], psfc[ids]
    else:
      return fcldfrac, pressure, psfc

def get_ctp(obsctpnc,ids):
    ctp = np.array(obsctpnc.PRES_G16.data[0,:],dtype=ioda_float_type) *100.0 # hPa to Pa
    ctp[np.ma.where( (ctp > 0.0) & (ctp < PTOP_HIGH_UPP) )] = np.nan
    bctp_obsloc = np.full_like(ids,np.nan,dtype=ioda_float_type)
    l = [ctp[ids[i]] for i in range(len(ids))]
    bctp_obsloc = np.asarray(l)
    return bctp_obsloc

def get_hydrometeors(geovalnc,ids,variable):
    hydroName = hydroDict[variable]
    hydro_ = np.array(geovalnc[hydroName].data,dtype=ioda_float_type)
    hydro = np.nansum(hydro_,axis=1)
    hydro[np.ma.where(hydro == hydro_initial)] = np.nan
    if ids is not None:
      return hydro[ids]
    else:
      return hydro

def cleanData(binningMethod, varsDict):
    if binningMethod == 'CloudFraction':
      refvar = 'ocldfrac'
    elif binningMethod == 'CTPLayering':
      refvar = 'octp'
    else:
      raise NotImplementedError()
    validInd = np.isfinite(varsDict[refvar])
    for k in varsDict.keys():
      varsDict[k] = varsDict[k][validInd]
    return varsDict

def check(list):
    # check that all elements in a list are identical
    return all(i == list[0] for i in list)

def get_layerDefinitions_era5(psfc):
    # https://github.com/byoung-joo/cloud_vx/blob/master/bin/python_stuff.py
    PTOP_LOW = 0.8*psfc # these are arrays
    PTOP_MID = 0.45*psfc
    PTOP_HIGH = PTOP_HIGH_UPP * np.ones_like(psfc)
    return PTOP_LOW, PTOP_MID, PTOP_HIGH

def get_layerDefinitions_upp(psfc):
    # https://github.com/byoung-joo/cloud_vx/blob/master/bin/python_stuff.py
    # psfc added to get the shape
    PTOP_LOW = PTOP_LOW_UPP * np.ones_like(psfc)
    PTOP_MID = PTOP_MID_UPP * np.ones_like(psfc)
    PTOP_HIGH = PTOP_HIGH_UPP * np.ones_like(psfc)
    return PTOP_LOW, PTOP_MID, PTOP_HIGH

def get_FcstCloudFrac_low(cfr,pmid,psfc,PTOP_LOW):
    # https://github.com/byoung-joo/cloud_vx/blob/master/bin/python_stuff.py
    # cfr is cloud fraction(%), pmid is 3D pressure(Pa), psfc is surface pressure (Pa) code from UPP ./INITPOST.F
    nlocs = pmid.shape[0]
    if pmid.shape != cfr.shape:  # sanity check
      print('dimension mismatch bewteen cldfra and pressure')
      sys.exit()

    if len(psfc) != nlocs: # another sanity check
      print('dimension mismatch bewteen cldfra and surface pressure')
      sys.exit()

    idx  = [np.where( pmid[i,:] >= PTOP_LOW[i] )[0] for i in range(nlocs)]
    cldfrac_low  = np.array([ np.max(cfr[i,idx[i]]) for i in range(nlocs) if (len(idx[i]) >0) ])

    return cldfrac_low

def get_FcstCloudFrac_mid(cfr,pmid,psfc,PTOP_LOW,PTOP_MID):
    # https://github.com/byoung-joo/cloud_vx/blob/master/bin/python_stuff.py
    # cfr is cloud fraction(%), pmid is 3D pressure(Pa), psfc is surface pressure (Pa) code from UPP ./INITPOST.F

    nlocs = pmid.shape[0]
    if pmid.shape != cfr.shape:  # sanity check
      print('dimension mismatch bewteen cldfra and pressure')
      sys.exit()

    if len(psfc) != nlocs: # another sanity check
      print('dimension mismatch bewteen cldfra and surface pressure')
      sys.exit()

    idx  = [np.where( (pmid[i,:] <  PTOP_LOW[i]) & (pmid[i,:] >= PTOP_MID[i]) )[0] for i in range(nlocs)]
    cldfrac_mid  = np.array([ np.max(cfr[i,idx[i]]) for i in range(nlocs) if (len(idx[i]) >0) ])

    return cldfrac_mid

def get_FcstCloudFrac_high(cfr,pmid,psfc,PTOP_MID,PTOP_HIGH):
    # https://github.com/byoung-joo/cloud_vx/blob/master/bin/python_stuff.py
    # cfr is cloud fraction(%), pmid is 3D pressure(Pa), psfc is surface pressure (Pa) code from UPP ./INITPOST.F
    nlocs = pmid.shape[0]
    if pmid.shape != cfr.shape:  # sanity check
      print('dimension mismatch bewteen cldfra and pressure')
      sys.exit()

    if len(psfc) != nlocs: # another sanity check
      print('dimension mismatch bewteen cldfra and surface pressure')
      sys.exit()

    idx  = [np.where( (pmid[i,:] <  PTOP_MID[i]) & (pmid[i,:] >= PTOP_HIGH[i]) )[0] for i in range(nlocs)]
    cldfrac_high  = np.array([ np.max(cfr[i,idx[i]]) for i in range(nlocs) if (len(idx[i]) >0) ])

    return cldfrac_high

def count_non_nan(data):
    return data.size - np.count_nonzero(np.isnan(data))

def dateList(dateIni,dateEnd,delta):
    datei = datetime.strptime(str(dateIni), "%Y%m%d%H")
    datef = datetime.strptime(str(dateEnd), "%Y%m%d%H")
    date  = datei
    dateList = list()
    ddList = list()

    while (date <= datef):
      datestr = date.strftime("%Y%m%d%H")
      dd      = date.strftime("%d")
      dateList.append(datestr)
      ddList.append(dd)
      date = date + timedelta(hours=int(delta))
    return dateList, ddList

def calc_stat(stat,var1,refvar=None):
    if refvar is None:
      valid = np.isfinite(var1)
    else:
      valid = np.isfinite(refvar)
    if len(var1[valid]) > 0:
      if stat == 'count':
        return len(var1[valid])
      elif (stat == 'mean'):
        return nanmean(var1[valid])
      elif (stat == 'rmse'):
        return rmse(var1[valid])
      elif (stat == 'std'):
        return np.nanstd(var1[valid])
      else:
        raise NotImplementedError()
    else:
      return np.nan

def nanmean(var):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(var)

def rmse(var1):
    return np.sqrt(np.nanvar(var1))

def omb(obs, bkg, bias=None):
    if bias is None:
      bias = 0
    return (obs - bkg - bias)

def log10(arr):
    # np.where(hist_v1>0, np.log10(hist_v1), np.nan)
    # https://www.geeksforgeeks.org/how-to-fix-runtimewarning-divide-by-zero-encountered-in-log/
    # Ignore the warning for invalid logarithmic operations
    with np.errstate(divide='ignore', invalid='ignore'):
        result = np.log10(arr)
    return result

def varbinning(bins,criteria,var1):
    return var1[bins[criteria]]

def get_binsIndexes_cldfrac(cldfraThresh,cfx,cfy):

    # Initialize dictionary
    bins_ids = dict.fromkeys(criteria_list_2by2,[])

    bins_ids['Na_Ocloud-Mcloud'] = np.where( (cfy >  cldfraThresh) & (cfx >  cldfraThresh) )[0]
    bins_ids['Nb_Oclear-Mcloud'] = np.where( (cfy <= cldfraThresh) & (cfx >  cldfraThresh) )[0]
    bins_ids['Nc_Ocloud-Mclear'] = np.where( (cfy >  cldfraThresh) & (cfx <= cldfraThresh) )[0]
    bins_ids['Nd_Oclear-Mclear'] = np.where( (cfy <= cldfraThresh) & (cfx <= cldfraThresh) )[0]
    bins_ids['all'] = np.arange(len(cfy))
    bins_ids['Ocloud'] = np.where( (cfy >  cldfraThresh) )[0]
    bins_ids['Oclear'] = np.where( (cfy <= cldfraThresh) )[0]

    return bins_ids

def get_binsIndexes_ctpLayering(cldfraThresh,cfx, cfxlow, cfxmid, cfxhigh, ctp, PTOP_LOW, PTOP_MID, PTOP_HIGH):

    # Initialize dictionary
    bins_ids = dict.fromkeys(criteria_list_4by4,[])

    bins_ids['Na_Ocloud-Mcloud'] = np.where( (ctp > 0.0)       & (cfx > cldfraThresh) )[0]

    bins_ids['Nb_Oclear-Mcloud'] = np.where( (ctp == -77700.0) & (cfx >  cldfraThresh) )[0]
    bins_ids['Nc_Ocloud-Mclear'] = np.where( (ctp > 0.0)       & (cfx <= cldfraThresh) )[0]
    bins_ids['Nd_Oclear-Mclear'] = np.where( (ctp == -77700.0) & (cfx <= cldfraThresh) )[0]

    bins_ids['all']    = np.arange(len(ctp))
    bins_ids['Ocloud'] = np.where( (ctp  > 0.0)                            )[0]
    bins_ids['Oclear'] = np.where( (ctp == -77700.0)                       )[0]
    bins_ids['Olow']   = np.where(                     (ctp >= PTOP_LOW)   )[0]
    bins_ids['Omid']   = np.where( ((ctp < PTOP_LOW) & (ctp >= PTOP_MID))  )[0]
    bins_ids['Ohigh']  = np.where( ((ctp < PTOP_MID) & (ctp >= PTOP_HIGH)) )[0]

    low  =  np.where( (cfxmid > cldfraThresh) | (cfxhigh > cldfraThresh), float_missingVal, cfxlow  )
    mid  =  np.where( (cfxlow > cldfraThresh) | (cfxhigh > cldfraThresh), float_missingVal, cfxmid  )
    high =  np.where( (cfxlow > cldfraThresh) | ( cfxmid > cldfraThresh), float_missingVal, cfxhigh )
    bins_ids['lowOnly']         = np.where(  (ctp >= PTOP_LOW)                       & (low  > cldfraThresh) )[0]
    bins_ids['midOnly']         = np.where( ((ctp <  PTOP_LOW) & (ctp >= PTOP_MID))  & (mid  > cldfraThresh) )[0]
    bins_ids['highOnly']        = np.where( ((ctp <  PTOP_MID) & (ctp >= PTOP_HIGH)) & (high > cldfraThresh) )[0]

    bins_ids['a_Olow-Mlow']     = np.where((ctp >= PTOP_LOW) & (low  > cldfraThresh))[0]
    bins_ids['f_Omid-Mmid']     = np.where((ctp <  PTOP_LOW) & (ctp >= PTOP_MID)  & (mid  > cldfraThresh))[0]
    bins_ids['k_Ohigh-Mhigh']   = np.where((ctp <  PTOP_MID) & (ctp >= PTOP_HIGH) & (high > cldfraThresh))[0]
    bins_ids['b_Omid-Mlow']     = np.where((ctp <  PTOP_LOW) & (ctp >= PTOP_MID)  & (cfxlow  > cldfraThresh))[0]
    bins_ids['c_Ohigh-Mlow']    = np.where((ctp <  PTOP_MID) & (ctp >= PTOP_HIGH) & (cfxlow  > cldfraThresh))[0]
    bins_ids['d_Oclear-Mlow']   = np.where((ctp == -77700.0) & (cfxlow  > cldfraThresh))[0]
    bins_ids['e_Olow-Mmid']     = np.where((ctp >= PTOP_LOW) & (cfxmid  > cldfraThresh))[0]
    bins_ids['g_Ohigh-Mmid']    = np.where((ctp <  PTOP_MID) & (ctp >= PTOP_HIGH) & (cfxmid  > cldfraThresh))[0]
    bins_ids['h_Oclear-Mmid']   = np.where((ctp == -77700.0) & (cfxmid  > cldfraThresh))[0]
    bins_ids['i_Olow-Mhigh']    = np.where((ctp >= PTOP_LOW) & (cfxhigh > cldfraThresh))[0]
    bins_ids['j_Omid-Mhigh']    = np.where((ctp <  PTOP_LOW) & (ctp >= PTOP_MID)  & (cfxhigh > cldfraThresh))[0]
    bins_ids['l_Oclear-Mhigh']  = np.where((ctp == -77700.0) & (cfxhigh > cldfraThresh))[0]
    bins_ids['m_Olow-Mclear']   = np.where((ctp >= PTOP_LOW) & (cfx <= cldfraThresh))[0]
    bins_ids['n_Omid-Mclear']   = np.where((ctp <  PTOP_LOW) & (ctp >= PTOP_MID)  & (cfx  <=  cldfraThresh))[0]
    bins_ids['o_Ohigh-Mclear']  = np.where((ctp <  PTOP_MID) & (ctp >= PTOP_HIGH) & (cfx  <=  cldfraThresh))[0]
    bins_ids['p_Oclear-Mclear'] = np.where((ctp == -77700.0) & (cfx <= cldfraThresh))[0]

    return bins_ids

def get_binsIndexes(binningMethod, cldfraThresh, ocldfrac, fcfrmax, fcstLow, fcstMid, fcstHigh, octp, PTOP_LOW, PTOP_MID, PTOP_HIGH):
    if binningMethod == 'CloudFraction':
      return get_binsIndexes_cldfrac(cldfraThresh, fcfrmax, ocldfrac)
    elif binningMethod == 'CTPLayering':
      return get_binsIndexes_ctpLayering(cldfraThresh, fcfrmax, fcstLow, fcstMid, fcstHigh, octp, PTOP_LOW, PTOP_MID, PTOP_HIGH)
    else:
      raise NotImplementedError()

def fbias(Na,Nb,Nc,Nd,T):
    try:
      FBIAS = (Na + Nb)  / (Na + Nc)
    except ZeroDivisionError as error:
      logging.exception('ZeroDivisionError: %s', error)
      FBIAS = np.nan
    return FBIAS

def pod(Na,Nb,Nc,Nd,T):
    try:
      POD = Na / (Na + Nc)
    except ZeroDivisionError as error:
      logging.exception('ZeroDivisionError: %s', error)
      POD = np.nan
    return POD

def podf(Na,Nb,Nc,Nd,T):
    try:
      PODF = Nb / (Nb + Nd)
    except ZeroDivisionError as error:
      logging.exception('ZeroDivisionError: %s', error)
      PODF = np.nan
    return PODF

def podn(Na,Nb,Nc,Nd,T):
    try:
      PODN = Nd / (Nb + Nd)
    except ZeroDivisionError as error:
      logging.exception('ZeroDivisionError: %s', error)
      PODN = np.nan
    return PODN

def far(Na,Nb,Nc,Nd,T):
    try:
      FAR = Nb / (Na + Nb)
    except ZeroDivisionError as error:
      logging.exception('ZeroDivisionError: %s', error)
      FAR = np.nan
    return FAR

def csi(Na,Nb,Nc,Nd,T):
    try:
      CSI = Na / (Na + Nb + Nc)
    except ZeroDivisionError as error:
      logging.exception('ZeroDivisionError: %s', error)
      CSI = np.nan
    return CSI

def gss(Na,Nb,Nc,Nd,T):
    try:
      C   = (Na + Nb)*(Na + Nc) / T
      GSS = (Na - C) / (Na + Nb + Nc - C)
    except ZeroDivisionError as error:
      logging.exception('ZeroDivisionError: %s', error)
      GSS = np.nan
    return GSS

def hrate(Na,Nb,Nc,Nd,T):
    try:
      HRATE = Na / T
    except ZeroDivisionError as error:
      logging.exception('ZeroDivisionError: %s', error)
      HRATE = np.nan
    return HRATE

def acc(Na,Nb,Nc,Nd,T):
    try:
      ACC = (Na + Nd) / T
    except ZeroDivisionError as error:
      logging.exception('ZeroDivisionError: %s', error)
      ACC = np.nan
    return ACC

def contigency_table_2by2(countDict, stat):

    Na = countDict[criteria_list_2by2[0]]
    Nb = countDict[criteria_list_2by2[1]]
    Nc = countDict[criteria_list_2by2[2]]
    Nd = countDict[criteria_list_2by2[3]]

    T = Na + Nb + Nc + Nd

    if stat == 'FBIAS':
      return fbias(Na,Nb,Nc,Nd,T)
    elif stat == 'POD':
      return pod(Na,Nb,Nc,Nd,T)
    elif stat == 'PODF':
      return podf(Na,Nb,Nc,Nd,T)
    elif stat == 'PODN':
      return podn(Na,Nb,Nc,Nd,T)
    elif stat == 'FAR':
      return far(Na,Nb,Nc,Nd,T)
    elif stat == 'CSI':
      return csi(Na,Nb,Nc,Nd,T)
    elif stat == 'GSS':
      return gss(Na,Nb,Nc,Nd,T)
    elif stat == 'HRATE':
      return hrate(Na,Nb,Nc,Nd,T)
    elif stat == 'ACC':
      return acc(Na,Nb,Nc,Nd,T)
    else:
      raise NotImplementedError()

def pod_cloudlayers(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,T,Na,Nb,Nc,Nd):
    # Initialize dictionary
    POD = dict.fromkeys(layering_list,np.nan)
    try:
      POD['all']        = Na / (Na + Nc)
      POD['LowClouds']  = a / (a + e + i + m)
      POD['MidClouds']  = b / (b + f + j + n)
      POD['HighClouds'] = c / (c + g + k + o)
    except ZeroDivisionError as error:
      logging.exception('ZeroDivisionError: %s', error)
      POD['all']        = np.nan
      POD['LowClouds']  = np.nan
      POD['MidClouds']  = np.nan
      POD['HighClouds'] = np.nan
    return POD

def far_cloudlayers(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,T,Na,Nb,Nc,Nd):
    # Initialize dictionary
    FAR = dict.fromkeys(layering_list,np.nan)
    try:
      FAR['all']        = Nb / (Na + Nb)
      FAR['LowClouds']  = d / (a + b + c + d)
      FAR['MidClouds']  = h / (e + f + g + h)
      FAR['HighClouds'] = l / (i + j + k + l)
    except ZeroDivisionError as error:
      logging.exception('ZeroDivisionError: %s', error)
      FAR['all']        = np.nan
      FAR['LowClouds']  = np.nan
      FAR['MidClouds']  = np.nan
      FAR['HighClouds'] = np.nan
    return FAR

def csi_cloudlayers(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,T,Na,Nb,Nc,Nd):
    # Initialize dictionary
    CSI = dict.fromkeys(layering_list,np.nan)
    try:
      CSI['all']        = Na / (Na + Nb + Nc)
      CSI['LowClouds']  = a  / (a + b + c + d + e + i + m)
      CSI['MidClouds']  = b  / (b + f + j + n + e + g + h)
      CSI['HighClouds'] = c  / (a + b + c + d + g + k + o)
    except ZeroDivisionError as error:
      logging.exception('ZeroDivisionError: %s', error)
      CSI['all']        = np.nan
      CSI['LowClouds']  = np.nan
      CSI['MidClouds']  = np.nan
      CSI['HighClouds'] = np.nan
    return CSI

def gss_cloudlayers(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,T,Na,Nb,Nc,Nd):
    # Initialize dictionary
    GSS = dict.fromkeys(layering_list,np.nan)
    try:
      C_all    = (Na + Nb)*(Na + Nc) / T
      C_low    = (a + b + c + d)*(a + e + i + m) / T
      C_mid    = (e + f + g + h)*(b + f + j + n) / T
      C_high   = (c + g + k + o)*(a + b + c + d) / T
      GSS['all']        = (Na - C_all) / (Na + Nb + Nc - C_all)
      GSS['LowClouds']  = (a - C_low)  / (a + b + c + d + e + i + m - C_low)
      GSS['MidClouds']  = (b - C_mid)  / (b + f + j + n + e + g + h - C_mid)
      GSS['HighClouds'] = (c - C_high) / (a + b + c + d + g + k + o - C_high)
    except ZeroDivisionError as error:
      logging.exception('ZeroDivisionError: %s', error)
      GSS['all']        = np.nan
      GSS['LowClouds']  = np.nan
      GSS['MidClouds']  = np.nan
      GSS['HighClouds'] = np.nan
    return GSS

def fbias_cloudlayers(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,T,Na,Nb,Nc,Nd):
    # Initialize dictionary
    FBIAS = dict.fromkeys(layering_list,np.nan)
    try:
      FBIAS['all']        = (Na + Nb)   / (Na + Nc)
      FBIAS['LowClouds']  = (a + b + c + d) / (a + e + i + m)
      FBIAS['MidClouds']  = (e + f + g + h) / (b + f + j + n)
      FBIAS['HighClouds'] = (i + j + k + l) / (c + g + k + o)
    except ZeroDivisionError as error:
      logging.exception('ZeroDivisionError: %s', error)
      FBIAS['all']        = np.nan
      FBIAS['LowClouds']  = np.nan
      FBIAS['MidClouds']  = np.nan
      FBIAS['HighClouds'] = np.nan
    return FBIAS

def contigency_table_4by4(countDict,stat):

    a  = countDict[criteria_list_4by4[0]]
    b  = countDict[criteria_list_4by4[1]]
    c  = countDict[criteria_list_4by4[2]]
    d  = countDict[criteria_list_4by4[3]]
    e  = countDict[criteria_list_4by4[4]]
    f  = countDict[criteria_list_4by4[5]]
    g  = countDict[criteria_list_4by4[6]]
    h  = countDict[criteria_list_4by4[7]]
    i  = countDict[criteria_list_4by4[8]]
    j  = countDict[criteria_list_4by4[9]]
    k  = countDict[criteria_list_4by4[10]]
    l  = countDict[criteria_list_4by4[11]]
    m  = countDict[criteria_list_4by4[12]]
    n  = countDict[criteria_list_4by4[13]]
    o  = countDict[criteria_list_4by4[14]]
    p  = countDict[criteria_list_4by4[15]]
    T = a + b + c + d + e + f + g + h + i + j + k + l + m + n + o + p

    Na = a + b + c + e + f + g + i + j + k
    Nb = d + h + l
    Nc = m + n + o
    Nd = p

    if stat == 'POD':
      return pod_cloudlayers(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,T,Na,Nb,Nc,Nd)
    elif stat == 'FAR':
      return far_cloudlayers(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,T,Na,Nb,Nc,Nd)
    elif stat == 'CSI':
      return csi_cloudlayers(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,T,Na,Nb,Nc,Nd)
    elif stat == 'GSS':
      return gss_cloudlayers(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,T,Na,Nb,Nc,Nd)
    elif stat == 'FBIAS':
      return fbias_cloudlayers(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,T,Na,Nb,Nc,Nd)
    else:
      raise NotImplementedError()
