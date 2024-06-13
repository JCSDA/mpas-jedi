import numpy as np
from datetime import datetime, timedelta
import sys, itertools, os, collections
import logging
import xarray as xr
import warnings
import multiprocessing

float_missingVal = float(-999)
int_missingVal = int(-999)
IODAfloat_missing_value = float(-3.3687953e+38)
IODAint_missing_value   = int(-2147483643)
hydro_initial = float(5.4999992258578274e-15)
ioda_float_type = 'float32'

cldfraThresh = 20.0

criteria_list_2by2 = ['Na_Ocloud-Mcloud','Nb_Oclear-Mcloud','Nc_Ocloud-Mclear','Nd_Oclear-Mclear','all','Ocloud','Oclear']
criteria_3p = ['Na_Ocloud-Mcloud','Nd_Oclear-Mclear', 'all']
criteria_list_4by4 = ['a_Olow-Mlow', 'b_Omid-Mlow', 'c_Ohigh-Mlow', 'd_Oclear-Mlow', 'e_Olow-Mmid', 'f_Omid-Mmid', 'g_Ohigh-Mmid', 'h_Oclear-Mmid', 'i_Olow-Mhigh', 'j_Omid-Mhigh', 'k_Ohigh-Mhigh', 'l_Oclear-Mhigh', 'm_Olow-Mclear', 'n_Omid-Mclear', 'o_Ohigh-Mclear', 'p_Oclear-Mclear', 'all', 'Ocloud', 'Oclear','lowOnly', 'midOnly', 'highOnly', 'Olow', 'Omid', 'Ohigh', 'Na_Ocloud-Mcloud','Nb_Oclear-Mcloud','Nc_Ocloud-Mclear','Nd_Oclear-Mclear']
criteria_only = ['lowOnly', 'midOnly', 'highOnly']
Na = ['a_Olow-Mlow', 'b_Omid-Mlow', 'c_Ohigh-Mlow', 'e_Olow-Mmid', 'f_Omid-Mmid', 'g_Ohigh-Mmid', 'i_Olow-Mhigh', 'j_Omid-Mhigh', 'k_Ohigh-Mhigh']
Nb = ['d_Oclear-Mlow', 'h_Oclear-Mmid', 'l_Oclear-Mhigh']
Nc = ['m_Olow-Mclear', 'n_Omid-Mclear', 'o_Ohigh-Mclear']
Nd = ['p_Oclear-Mclear']

cat_stats = ['FBIAS','POD','FAR','CSI','GSS','POFD','PODN','HRATE','ACC']
cat_stats_reduced = cat_stats[:5]
cont_stats = ['mean', 'count', 'rmse', 'std']

layering_list = ['all','LowClouds','MidClouds','HighClouds']

binningDict = {'CloudFraction': criteria_list_2by2,
               'CTPLayering': criteria_list_4by4}

# UPP top layer bounds (Pa) for cloud layers
PTOP_LOW_UPP  = 64200. # low for > 64200 Pa
PTOP_MID_UPP  = 35000. # mid between 35000-64200 Pa
PTOP_HIGH_UPP = 15000. # high between 15000-35000 Pa

varsList = ['obs', 'bkg', 'ihumr', 'imclw', 'imcli', 'imclr', 'imcls', 'imclg']
hydrosList = varsList[2:]

colorsL = ['k','r','b','g','c','m','y']

############################# FUNCTIONS #######################################
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

def varDict_binned(List):
    fillValue = np.nan
    return {var: {k: fillValue for k in List} for var in varsList}

def varDict3L(List1,List2,List3,fillValue=None):
    if fillValue == 'NaN':
      return {var: {l1: {l2: {l3: np.nan for l3 in List3} for l2 in List2} for l1 in List1} for var in varsList}
    elif fillValue == 0:
      return {var: {l1: {l2: {l3: 0 for l3 in List3} for l2 in List2} for l1 in List1} for var in varsList}
    else:
      return {var: {l1: {l2: {l3: np.array([]) for l3 in List3} for l2 in List2} for l1 in List1} for var in varsList}

def varContStatDict(List1,List2,fillValue=None):
    if fillValue == 'NaN':
      return {stat: {l1: {l2: np.nan for l2 in List2} for l1 in List1} for stat in cont_stats}
    elif fillValue == 0:
      return {stat: {l1: {l2: 0 for l2 in List2} for l1 in List1} for stat in cont_stats}
    else:
      return {stat: {l1: {l2: np.array([]) for l2 in List2} for l1 in List1} for stat in cont_stats}

def varContHydroStatDict(List1,List2,fillValue=None):
    if fillValue == 'NaN':
      return {var: {stat: {l1: {l2: np.nan for l2 in List2} for l1 in List1} for stat in cont_stats} for var in hydrosList}
    elif fillValue == 0:
      return {var: {stat: {l1: {l2: 0 for l2 in List2} for l1 in List1} for stat in cont_stats} for var in hydrosList}
    else:
      return {var: {stat: {l1: {l2: np.array([]) for l2 in List2} for l1 in List1} for stat in cont_stats} for var in hydrosList}

def createDict(List,fillValue=None):
    if fillValue == 'NaN':
      return {k: np.nan for k in List}
    elif fillValue == 0:
      return {k: 0 for k in List}
    else:
      return {k: [] for k in List}

def createNestedDictT(List1,List2,List3,fillValue=None):
    if fillValue == 'NaN':
      return {l1: {l2: {l3: np.nan for l3 in List3} for l2 in List2} for l1 in List1}
    elif fillValue == 0:
      return {l1: {l2: {l3: 0 for l3 in List3} for l2 in List2} for l1 in List1}
    else:
      return {l1: {l2: {l3: np.array([]) for l3 in List3} for l2 in List2} for l1 in List1}

def createDictL(List,fillList):
    return {List[k]: fillList[k] for k in range(len(List))}

def m2tokm2(data):
    return data*1e-6

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
    ocellindex = np.array(metadata.cellIndex.data-int(1),dtype=int) # ocellindex starts in 1 (from fortran), here it has to start from 0 (python)
    olatitude  = metadata.latitude.data
    olongitude = metadata.longitude.data # (longitude + 180) % 360 - 180
    areaCell = np.array( staticnc.variables['areaCell'][:] )
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
    obt        = obsvalue.brightnessTemperature.isel(Channel=ch).data
    obt[np.ma.where(obt == 0.0)] = np.nan
    if ids is not None:
      return obt[ids]
    else:
      return obt

def get_hofx(hofx,ch,ids):
    fbt        = hofx.brightnessTemperature.isel(Channel=ch).data
    fbt[np.ma.where(fbt == 0.0)] = np.nan
    if ids is not None:
      return fbt[ids]
    else:
      return fbt

def get_cldfrac(metadata,ids):
    ocldfrac   = metadata.cloudAmount.data *100 #Unitless to percent
    if ids is not None:
      return ocldfrac[ids]
    else:
      return ocldfrac

def get_geoval(geovalnc,ids):
    fcldfrac   = geovalnc.cloud_area_fraction_in_atmosphere_layer.data *100.0 #Unitless to percent
    pressure   = geovalnc.air_pressure.data
    psfc       = geovalnc.air_pressure_levels.data[:,-1]
    if ids is not None:
      return fcldfrac[ids], pressure[ids], psfc[ids]
    else:
      return fcldfrac, pressure, psfc

def get_ctp(obsctpnc,ids):
    ctp = obsctpnc.PRES_G16.data[0,:] *100.0 # hPa to Pa
    ctp[np.ma.where( (ctp > 0.0) & (ctp < PTOP_HIGH_UPP) )] = np.nan
    bctp_obsloc = np.full_like(ids,np.nan,dtype='float32')
    l = [ctp[ids[i]] for i in range(len(ids))]
    bctp_obsloc = np.asarray(l)
    return bctp_obsloc

def get_hydrometeors(geovalnc,ids):
    humr  = geovalnc.humidity_mixing_ratio.data
    ihumr = np.sum(humr,axis=1)
    ihumr[np.ma.where(ihumr == hydro_initial)] = np.nan

    mclw  = geovalnc.mass_content_of_cloud_liquid_water_in_atmosphere_layer.data
    imclw = np.sum(mclw,axis=1)
    imclw[np.ma.where(imclw == hydro_initial)] = np.nan

    mcli  = geovalnc.mass_content_of_cloud_ice_in_atmosphere_layer.data
    imcli = np.sum(mcli,axis=1)
    imcli[np.ma.where(imcli == hydro_initial)] = np.nan

    mclr  = geovalnc.mass_content_of_rain_in_atmosphere_layer.data
    imclr = np.sum(mclr,axis=1)
    imclr[np.ma.where(imclr == hydro_initial)] = np.nan

    mcls  = geovalnc.mass_content_of_snow_in_atmosphere_layer.data
    imcls = np.sum(mcls,axis=1)
    imcls[np.ma.where(imcls == hydro_initial)] = np.nan

    mclg  = geovalnc.mass_content_of_graupel_in_atmosphere_layer.data
    imclg = np.sum(mclg,axis=1)
    imclg[np.ma.where(imclg == hydro_initial)] = np.nan

    if ids is not None:
      return ihumr[ids], imclw[ids], imcli[ids], imclr[ids], imcls[ids], imclg[ids]
    else:
      return ihumr, imclw, imcli, imclr, imcls, imclg

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

    del idx
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

    del idx
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

    del idx
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
      p = np.isfinite(var1)
    else:
      p = np.isfinite(refvar)
    nonnan = count_non_nan(var1[p])
    if ((nonnan > 0) & (len(var1[p]) > 0)):
      if stat == 'count':
        return len(var1[p])
      elif (stat == 'mean'):
        return nanmean(var1[p])
      elif (stat == 'rmse'):
        return rmse(var1[p])
      elif (stat == 'std'):
        return np.nanstd(var1[p])
      else:
        raise NotImplementedError()
    else:
      return np.nan

def nanmean(var):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(var)

def rmse(var1):
    return np.sqrt(np.var(var1))

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

def get_binsIndexes_cldfrac(cfx,ncfy):
    # Use only valid locations
    nonan = np.isfinite(ncfy)
    cfy = np.array(ncfy[nonan], dtype=ioda_float_type)
    cfx = np.array(cfx[nonan], dtype=ioda_float_type)

    # Initialize dictionary
    bins_ids = createDict(criteria_list_2by2)

    if count_non_nan(cfy) == count_non_nan(cfx):  # these should all match
      bins_ids['Na_Ocloud-Mcloud'] = np.where( (cfy >= cldfraThresh) & (cfx >= cldfraThresh) )[0]
      bins_ids['Nb_Oclear-Mcloud'] = np.where( (cfy <  cldfraThresh) & (cfx >= cldfraThresh) )[0]
      bins_ids['Nc_Ocloud-Mclear'] = np.where( (cfy >= cldfraThresh) & (cfx <  cldfraThresh) )[0]
      bins_ids['Nd_Oclear-Mclear'] = np.where( (cfy <  cldfraThresh) & (cfx <  cldfraThresh) )[0]
      bins_ids['all'] = np.arange(count_non_nan(cfy))
      bins_ids['Ocloud'] = np.where( cfy >= cldfraThresh )[0]
      bins_ids['Oclear'] = np.where( cfy < cldfraThresh )[0]

    return bins_ids

def get_binsIndexes_ctpLayering(fcfrmax, fcstLow, fcstMid, fcstHigh, octp, fPTOP_LOW, fPTOP_MID, fPTOP_HIGH):
    # Use only valid locations
    nonan = np.isfinite(octp)
    ctp = np.array(octp[nonan], dtype=ioda_float_type)
    cfx = np.array(fcfrmax[nonan], dtype=ioda_float_type)
    cfxlow = np.array(fcstLow[nonan], dtype=ioda_float_type)
    cfxmid = np.array(fcstMid[nonan], dtype=ioda_float_type)
    cfxhigh = np.array(fcstHigh[nonan], dtype=ioda_float_type)
    PTOP_LOW = np.array(fPTOP_LOW[nonan], dtype=ioda_float_type)
    PTOP_MID = np.array(fPTOP_MID[nonan], dtype=ioda_float_type)
    PTOP_HIGH = np.array(fPTOP_HIGH[nonan], dtype=ioda_float_type)

    # Initialize dictionary
    bins_ids = createDict(criteria_list_4by4)

    if count_non_nan(ctp) == count_non_nan(cfx) == count_non_nan(cfxlow) == count_non_nan(cfxmid) == count_non_nan(cfxhigh) == count_non_nan(PTOP_LOW) == count_non_nan(PTOP_LOW) == count_non_nan(PTOP_MID) == count_non_nan(PTOP_HIGH):  # these should all match

      bins_ids['Na_Ocloud-Mcloud'] = np.where( (ctp > 0.0)       & (cfx >= cldfraThresh) )[0]
      bins_ids['Nb_Oclear-Mcloud'] = np.where( (ctp == -77700.0) & (cfx >= cldfraThresh) )[0]
      bins_ids['Nc_Ocloud-Mclear'] = np.where( (ctp > 0.0)       & (cfx <  cldfraThresh) )[0]
      bins_ids['Nd_Oclear-Mclear'] = np.where( (ctp == -77700.0) & (cfx <  cldfraThresh) )[0]

      bins_ids['all']    = np.arange(count_non_nan(ctp))
      bins_ids['Ocloud'] = np.where( ctp  > 0.0 )[0]
      bins_ids['Oclear'] = np.where( ctp == -77700.0 )[0]
      bins_ids['Olow']   = np.where(                     (ctp >= PTOP_LOW)   )[0]
      bins_ids['Omid']   = np.where( ((ctp < PTOP_LOW) & (ctp >= PTOP_MID))  )[0]
      bins_ids['Ohigh']  = np.where( ((ctp < PTOP_MID) & (ctp >= PTOP_HIGH)) )[0]

      low  =  np.where( (cfxmid >= cldfraThresh) | (cfxhigh >= cldfraThresh), float_missingVal, cfxlow  )
      mid  =  np.where( (cfxlow >= cldfraThresh) | (cfxhigh >= cldfraThresh), float_missingVal, cfxmid  )
      high =  np.where( (cfxlow >= cldfraThresh) | ( cfxmid >= cldfraThresh), float_missingVal, cfxhigh )
      bins_ids['lowOnly']         = np.where(  (ctp >= PTOP_LOW)                       & (low  >= cldfraThresh) )[0]
      bins_ids['midOnly']         = np.where( ((ctp <  PTOP_LOW) & (ctp >= PTOP_MID))  & (mid  >= cldfraThresh) )[0]
      bins_ids['highOnly']        = np.where( ((ctp <  PTOP_MID) & (ctp >= PTOP_HIGH)) & (high >= cldfraThresh) )[0]

      remaining = []
      for i in range(len(ctp)):
        if (ctp[i] >= PTOP_LOW[i])                                & (low[i]  >= cldfraThresh):
          bins_ids['a_Olow-Mlow'].append(i)

        elif ((ctp[i] <  PTOP_LOW[i]) & (ctp[i] >= PTOP_MID[i]))  & (mid[i]  >= cldfraThresh):
          bins_ids['f_Omid-Mmid'].append(i)

        elif ((ctp[i] <  PTOP_MID[i]) & (ctp[i] >= PTOP_HIGH[i])) & (high[i] >= cldfraThresh):
          bins_ids['k_Ohigh-Mhigh'].append(i)

        elif ((ctp[i] <  PTOP_LOW[i]) & (ctp[i] >= PTOP_MID[i]))  & (cfxlow[i]  >= cldfraThresh):
          bins_ids['b_Omid-Mlow'].append(i)

        elif ((ctp[i] <  PTOP_MID[i]) & (ctp[i] >= PTOP_HIGH[i])) & (cfxlow[i]  >= cldfraThresh):
           bins_ids['c_Ohigh-Mlow'].append(i)

        elif (ctp[i] == -77700.0)                                 & (cfxlow[i]  >= cldfraThresh):
          bins_ids['d_Oclear-Mlow'].append(i)

        elif (ctp[i] >= PTOP_LOW[i])                              & (cfxmid[i]  >= cldfraThresh):
          bins_ids['e_Olow-Mmid'].append(i)

        elif ((ctp[i] <  PTOP_MID[i]) & (ctp[i] >= PTOP_HIGH[i])) & (cfxmid[i]  >= cldfraThresh):
          bins_ids['g_Ohigh-Mmid'].append(i)

        elif (ctp[i] == -77700.0)                                 & (cfxmid[i]  >= cldfraThresh):
          bins_ids['h_Oclear-Mmid'].append(i)

        elif (ctp[i] >= PTOP_LOW[i])                              & (cfxhigh[i] >= cldfraThresh):
          bins_ids['i_Olow-Mhigh'].append(i)

        elif ((ctp[i] <  PTOP_LOW[i]) & (ctp[i] >= PTOP_MID[i]))  & (cfxhigh[i] >= cldfraThresh):
          bins_ids['j_Omid-Mhigh'].append(i)

        elif (ctp[i] == -77700.0)                                 & (cfxhigh[i] >= cldfraThresh):
          bins_ids['l_Oclear-Mhigh'].append(i)

        elif (ctp[i] >= PTOP_LOW[i])                              & (cfx[i]  <  cldfraThresh):
          bins_ids['m_Olow-Mclear'].append(i)

        elif ((ctp[i] <  PTOP_LOW[i]) & (ctp[i] >= PTOP_MID[i]))  & (cfx[i]  <  cldfraThresh):
          bins_ids['n_Omid-Mclear'].append(i)

        elif ((ctp[i] <  PTOP_MID[i]) & (ctp[i] >= PTOP_HIGH[i])) & (cfx[i]  <  cldfraThresh):
          bins_ids['o_Ohigh-Mclear'].append(i)

        elif (ctp[i] == -77700.0)                                 & (cfx[i]  <  cldfraThresh):
          bins_ids['p_Oclear-Mclear'].append(i)
        else:
          remaining.append(i)

      try:
        if len(remaining) != 0:  # sanity check
          raise Exception('nObs: {} outside conditions'.format(str(len(remaining))))
      except Exception as e:
          logging.error(str(e))
          raise   # propagate the exception to the log file

    return bins_ids

def get_binsIndexes(binningMethod, ocldfrac, fcfrmax, fcstLow, fcstMid, fcstHigh, octp, PTOP_LOW, PTOP_MID, PTOP_HIGH):
    if binningMethod == 'CloudFraction':
      return get_binsIndexes_cldfrac(fcfrmax, ocldfrac)
    elif binningMethod == 'CTPLayering':
      return get_binsIndexes_ctpLayering(fcfrmax, fcstLow, fcstMid, fcstHigh, octp, PTOP_LOW, PTOP_MID, PTOP_HIGH)
    else:
      raise NotImplementedError()

def contigency_table_2by2(countDict):

    Na = countDict[criteria_list_2by2[0]]
    Nb = countDict[criteria_list_2by2[1]]
    Nc = countDict[criteria_list_2by2[2]]
    Nd = countDict[criteria_list_2by2[3]]

    T = Na + Nb + Nc + Nd
    try:
      FBIAS = (Na + Nb) / (Na + Nc)
      POD   = Na / (Na + Nc)
      POFD  = Nb / (Nb + Nd)
      PODN  = Nd / (Nb + Nd)
      FAR   = Nb / (Nb + Na)
      CSI   = Na / (Na + Nb + Nc)
      C     = (Na + Nb)*(Na + Nc) / T
      GSS   = (Na - C) / (Na + Nb + Nc - C)
      RATE  = Na / T
      ACC   = (Na + Nd) / T
    except ZeroDivisionError as error:
      logging.exception('ZeroDivisionError: %s', error)
      FBIAS = np.nan
      POD   = np.nan
      POFD  = np.nan
      PODN  = np.nan
      FAR   = np.nan
      CSI   = np.nan
      C     = np.nan
      GSS   = np.nan
      RATE  = np.nan
      ACC   = np.nan
    del Na, Nb, Nc, Nd
    return FBIAS, POD, POFD, PODN, FAR, CSI, GSS, RATE, ACC

def pod_cloudlayers(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,T,Na,Nb,Nc,Nd):
    # Initialize dictionary
    POD = createDict(layering_list)

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
    FAR = createDict(layering_list)

    try:
      FAR['all']        = Nb / (Nb + Na)
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
    CSI = createDict(layering_list)

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
    GSS = createDict(layering_list)

    try:
      C_all    = (Na + Nb)*(Na + Nc) / T
      C_low    = (a + b + c + d)*(a + e + i + m) / T
      C_mid    = (e + f + g + h)*(b + f + j + n) / T
      C_high   = (c + g + k + o)*(a + b + c + d) / T
      GSS['all']        = (Na - C_all) / (Na + Nb + Nc - C_all)
      GSS['LowClouds']  = (a - C_low)     / (a + b + c + d + e + i + m - C_low)
      GSS['MidClouds']  = (b - C_mid)     / (b + f + j + n + e + g + h - C_mid)
      GSS['HighClouds'] = (c - C_high)    / (a + b + c + d + g + k + o - C_high)
    except ZeroDivisionError as error:
      logging.exception('ZeroDivisionError: %s', error)
      GSS['all']        = np.nan
      GSS['LowClouds']  = np.nan
      GSS['MidClouds']  = np.nan
      GSS['HighClouds'] = np.nan
    return GSS

def fbias_cloudlayers(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,T,Na,Nb,Nc,Nd):
    # Initialize dictionary
    FBIAS = createDict(layering_list)

    try:
      FBIAS['all']        = (Na + Nb)   / (Na + Nc)
      FBIAS['LowClouds']  = (a + b + c) / (a + e + i)
      FBIAS['MidClouds']  = (e + f + g) / (b + f + j)
      FBIAS['HighClouds'] = (i + j + k) / (c + g + k)
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
