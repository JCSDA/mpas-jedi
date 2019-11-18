import os
import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('pdf')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime as dt 
from datetime import *
import ast

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
DA_METHOD1 = str(sys.argv[2])
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

     anafile = '../mpas.'+DA_METHOD1+'.'+datefile+'.nc'
     ananc = Dataset(anafile, "r")

     #use gfsana
     if (USE_GFSANA):
         gfsfile = GFSANA_DIR+'/x1.40962.init.'+datefile+'.nc'
         gfsnc = Dataset(gfsfile, "r")

     lats = np.array( baknc.variables['latCell'][:] ) * 180.0 / np.pi
     lons = np.array( baknc.variables['lonCell'][:] ) * 180.0 / np.pi
     xtime = np.array( baknc.variables['xtime'][:] )
     nVertL = baknc.dimensions['nVertLevels']

     for lvl in range(0,nVertL.size,int(INTERVAL)):
         print('lvl=',lvl)
         bakmpas = np.array( baknc.variables[VAR_NAME][0,:,lvl-1] )
         anampas = np.array( ananc.variables[VAR_NAME][0,:,lvl-1] )
         ambmpas = anampas - bakmpas
         plot(lats,lons,anampas,DATE,lvl,'MPASANA')
         plot(lats,lons,bakmpas,DATE,lvl,'MPASBAK')
         plot(lats,lons,ambmpas,DATE,lvl,'MPASAMB')

         if (USE_GFSANA):
             gfs     = np.array( gfsnc.variables[VAR_NAME][0,:,lvl-1] )
             anamgfs = anampas - gfs
             bakmgfs = bakmpas - gfs

            #compare with gfsana: plot()
             plot(lats,lons,anamgfs,DATE,lvl,'MPASANA-GFSANA')
             plot(lats,lons,bakmgfs,DATE,lvl,'MPASBAK-GFSANA')

def plot(lats,lons,data,yyyymmddhh,lvl,source):
    fig,ax=plt.subplots(figsize=(8,8))
    map=Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=0,urcrnrlon=360,resolution='c')
    map.drawcoastlines()
    map.drawcountries()
    map.drawparallels(np.arange(-90,90,30),labels=[1,1,0,1])
    map.drawmeridians(np.arange(0,360,60),labels=[1,1,0,1])

    plt.title( source+'  '+VAR_NAME+', lvl=' + str(lvl)+",  "+yyyymmddhh)
    #cont=plt.tricontourf(lons,lats,data,30,vmin=-(max(abs(data))),vmax=max(abs(data)),extend='both',cmap=cm.bwr)
    cont=plt.tricontourf(lons,lats,data,30,cmap=cm.bwr) 
#create axes for colorbar======================================================
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom",
        size="3%", #width of the colorbar
        pad=0.3)   #space between colorbar and graph
    plt.colorbar(cont,cax=cax,ax=ax,orientation='horizontal')
    cax.xaxis.set_ticks_position('bottom')  #set the orientation of the ticks of the colorbar

#save figures 
    plt.savefig('distr_%s_%s_%s_%s.png'%(VAR_NAME,lvl,yyyymmddhh,source),dpi=200,bbox_inches='tight')
    plt.close()

def main():
    readdata()

if __name__ == '__main__': main()


