import os
import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('pdf')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime as dt 
from datetime import *
import ast
import var_utils as vu

#This script is for plotting MPAS background, analysis, increment, MPAS background and analysis compare to GFS analysis  
#python plot_inc.py $DATE method      variable   level_interval   USE_GFSANA   $GFSANA_DIR
#python plot_inc.py $DATE 3denvar_bump theta        5              True        /glade/work/${USER}/pandac/120km_GFSANA
#python plot_inc.py $DATE 3denvar_bump theta        5              False

'''
Directory Structure:
test/
 ├── mpas.3denvar_bump.2018-04-15_00.00.00.nc
 ├── restart.2018-04-15_00.00.00.nc
 ├── graphics/
 │   ├── plot_inc.py
 │   ├── ...
'''

DATE = str(sys.argv[1])
DA_METHOD = str(sys.argv[2])
VAR_NAME = str(sys.argv[3])
INTERVAL = str(sys.argv[4])

USE_GFSANA = ast.literal_eval(str(sys.argv[5])) #convert str to boolean
if (USE_GFSANA):
    GFSANA_DIR = str(sys.argv[6])

def readdata():
     datein = datetime.strptime(DATE+'00','%Y%m%d%H%M%S')
     print('datein=',datein)
     datefile = datein.strftime('%Y-%m-%d_%H.%M.%S')

     bakfile = '../restart.'+datefile+'.nc'
     baknc = Dataset(bakfile, "r")

     anafile = '../mpas.'+DA_METHOD+'.'+datefile+'.nc'
     ananc = Dataset(anafile, "r")

     #use gfsana
     if (USE_GFSANA):
         gfsfile = GFSANA_DIR+'/x1.40962.init.'+datefile+'.nc'
         gfsnc = Dataset(gfsfile, "r")

     lats = np.array( baknc.variables['latCell'][:] ) * 180.0 / np.pi
     lons = np.array( baknc.variables['lonCell'][:] ) * 180.0 / np.pi
     xtime = np.array( baknc.variables['xtime'][:] )
     nVertLSize = baknc.dimensions['nVertLevels'].size
     if ( VAR_NAME == "surface_pressure" ): nVertLSize=1

     for lvl in range(0,nVertLSize,int(INTERVAL)):
         print('lvl=',lvl)
         if ( VAR_NAME != "surface_pressure"):
             bakmpas = np.array( baknc.variables[VAR_NAME][0,:,lvl] )
             anampas = np.array( ananc.variables[VAR_NAME][0,:,lvl] )
         else:
             bakmpas = np.array( baknc.variables[VAR_NAME][0,:] )
             anampas = np.array( ananc.variables[VAR_NAME][0,:] )
         ambmpas = anampas - bakmpas
         plot(lats,lons,anampas,DATE,lvl+1,'MPASANA')
         plot(lats,lons,bakmpas,DATE,lvl+1,'MPASBAK')
         plot(lats,lons,ambmpas,DATE,lvl+1,'MPASAMB')

         if (USE_GFSANA):
             if ( VAR_NAME != "surface_pressure"):
                 gfs     = np.array( gfsnc.variables[VAR_NAME][0,:,lvl] )
             else:
                 gfs     = np.array( gfsnc.variables[VAR_NAME][0,:] )
             anamgfs = anampas - gfs
             bakmgfs = bakmpas - gfs

            #compare with gfsana: plot()
             plot(lats,lons,anamgfs,DATE,lvl+1,'MPASANA-GFSANA')
             plot(lats,lons,bakmgfs,DATE,lvl+1,'MPASBAK-GFSANA')

def plot(lats,lons,data,yyyymmddhh,lvl,source):
    fig,ax=plt.subplots(figsize=(8,8))

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_global()
    ax.gridlines(draw_labels=True, xlocs=np.arange(-180,180,60),linestyle='--')

    cont=plt.tricontourf(lons, lats, data, 60,
             transform=ccrs.PlateCarree())
    plt.title( source+'  '+VAR_NAME+ ' (' +vu.varDictModel[VAR_NAME][0]+'), lev=' + str(lvl)+",  "+yyyymmddhh)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom",size="5%", pad=0.2,axes_class=plt.Axes)
    plt.colorbar(cont,cax=cax,orientation='horizontal') #,cax=cax,ax=ax,orientation='horizontal')
    plt.savefig('distr_%s_%s_%s_%s.png'%(VAR_NAME,lvl,yyyymmddhh,source),dpi=200,bbox_inches='tight')
    plt.close()

def main():
    readdata()

if __name__ == '__main__': main()


