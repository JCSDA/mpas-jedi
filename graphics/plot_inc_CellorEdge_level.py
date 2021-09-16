import os
import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('pdf')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime as dt 
from datetime import *
import ast
import var_utils as vu

#This script is for plotting increments in 2D. 
#x-axis:Cells/edges, y-axis: verfical levels
#How to use it:
#python plot_inc_celloredge_level.py 2018041500 theta False

DATE = str(sys.argv[1])
VAR_NAME = str(sys.argv[2])

USE_GFSANA = ast.literal_eval(str(sys.argv[3])) #convert str to boolean
if (USE_GFSANA):
    GFSANA_DIR = str(sys.argv[4])

#nums used in 120km:
if (VAR_NAME =='u'):
   nums = 122880
else:
   nums = 40962

def readdata():
    xarray = np.arange(0,nums)
    binsfory = np.arange(0,55)
     
    datein = datetime.strptime(DATE+'00','%Y%m%d%H%M%S')
    datefile = datein.strftime('%Y-%m-%d_%H.%M.%S')

    bakfile = '../restart.'+datefile+'.nc_orig'
    baknc = Dataset(bakfile, "r")

    anafile = '../analysis.'+datefile+'.nc'
    ananc = Dataset(anafile, "r")

    #use gfsana
    if (USE_GFSANA):
        gfsfile = GFSANA_DIR+'/x1.40962.init.'+datefile+'.nc'
        gfsnc = Dataset(gfsfile, "r")

    lats = np.array( baknc.variables['latCell'][:] ) * 180.0 / np.pi
    lons = np.array( baknc.variables['lonCell'][:] ) * 180.0 / np.pi
    xtime = np.array( baknc.variables['xtime'][:] )
    nVertL = baknc.dimensions['nVertLevels']

    bakmpas = np.array( baknc.variables[VAR_NAME][0,:,:] )
    if (VAR_NAME == 'pressure_p'):
        pressure = np.array( ananc.variables['pressure'][0,:,:] )
        pressure_base = np.array( ananc.variables['pressure_base'][0,:,:] )
        anampas = pressure - pressure_base
    else:
        anampas = np.array( ananc.variables[VAR_NAME][0,:,:] )
    ambmpas = anampas - bakmpas
    #plot(xarray,binsfory,anampas,DATE,'MPASANA')
    #plot(xarray,binsfory,bakmpas,DATE,'MPASBAK')
    plot(xarray,binsfory,ambmpas,DATE,'MPASAMB')
    if (USE_GFSANA):
        gfs     = np.array( gfsnc.variables[VAR_NAME][0,:,:] )
        anamgfs = anampas - gfs
        bakmgfs = bakmpas - gfs

        #compare with gfsana: plot()
        plot(xarray,binsfory,anamgfs,DATE,'MPASANA-GFSANA')
        plot(xarray,binsfory,bakmgfs,DATE,'MPASBAK-GFSANA')
def plot(xarray,binsfory,model,yyyymmddhh,source):
    zgrid = np.loadtxt("/glade/work/jban/pandac/fix_input/graphics/zgrid_v55.txt")

    fig, ax1 = plt.subplots()
    cmap = 'coolwarm' #'bwr' #'jet'# 'coolwarm'
    model = model.reshape(len(xarray),len(binsfory)).T
    valuemin = np.amin(model) 
    valuemax = np.amax(model)
    norm = matplotlib.colors.DivergingNorm(vmin=valuemin, vcenter=0, vmax=valuemax)
    plt.contourf(xarray,binsfory,model,10,vmin=valuemin,vmax=valuemax,norm=norm,cmap=cmap)
    plt.title( source+' '+vu.varDictModel[VAR_NAME][1]+'('+ vu.varDictModel[VAR_NAME][0]+') max=' +str(round(valuemax,4))+' min='+str(round(valuemin,4)))
    major_ticks = np.arange(0, 56, 5)
    ax1.set_yticks(major_ticks)
    ax1.set_ylim([0,54])
    ax1.set_ylabel('Vertical level',fontsize=15)

    ax2 = ax1.twinx()
    ax2.set_yticks(major_ticks-1)
    ax2.set_yticklabels((zgrid[::5]).astype(int))

    ax2.set_ylabel('Height (m)',fontsize=13)

    if ( VAR_NAME == 'u'):
        ax1.set_xlabel( 'Edge',fontsize=15)
    else:
        ax1.set_xlabel( 'Cell',fontsize=15)

    plt.colorbar(extend='both',orientation="horizontal")
    plt.savefig('%s_%s_%s.png'%(vu.varDictModel[VAR_NAME][1],yyyymmddhh,source),dpi=200,bbox_inches='tight')

def main():
    readdata()

if __name__ == '__main__': main()
