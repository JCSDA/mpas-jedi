#import os
#import sys
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
#from mpl_toolkits.basemap import Basemap

#load background file & read lat/lon/time
file = 'restart.2018-04-15_00.00.00.nc'
f = Dataset(file, "r", format="NETCDF4") #Dataset is the class behavior to open the file
lats = np.array( f.variables['latCell'][:] )   # extract/copy the data
lons = np.array( f.variables['lonCell'][:] ) 
zgrid = np.array( f.variables['zgrid'][:] )
#print(zgrid.shape)  # (40962, 56)

nz=f.dimensions['nVertLevels'].size
nzp1=f.dimensions['nVertLevelsP1'].size
print( str(nz) + ', ' + str(nzp1) )

z_havg = np.zeros(nz)
z_havgP1 = np.average( zgrid , axis=0 )
for k in range(z_havg.shape[0]):
   z_havg[k] = 0.5 * ( z_havgP1[k] + z_havgP1[k+1] ) / 1000.0  # [km] unit


#for current exp
file = 'mpas.cor_rv.2018-04-15_00.00.00.nc'
f_exp = Dataset(file, "r", format="NETCDF4")

#for comparison
file = '../../20220111/CMAT_00/mpas.cor_rv.2018-04-15_00.00.00.nc'
f_ref = Dataset(file, "r", format="NETCDF4")

list_variables=['stream_function', 'velocity_potential', 'temperature', 'spechum'] #, 'surface_pressure']

for var in list_variables:
   print('Processing for variable ' + var)   # [Time, nCells, nVertLevels]

   if var == 'surface_pressure' :
      lvl=0
      data_exp = np.array( f_exp.variables[var][0,:] ) # (nl0, nc0)
      data_ref = np.array( f_ref.variables[var][0,:] ) # (nl0, nc0)
      print( '*** For '+var+': exp / ref = ', np.mean(data_exp, axis=0), np.mean(data_ref, axis=0) ) # [lvls] unit
   else:
      data_exp = np.array( f_exp.variables[var][0,:,:] ) # (nl0, nc0)
      data_ref = np.array( f_ref.variables[var][0,:,:] ) # (nl0, nc0)

   #set map
   if var == 'stream_function' :
      vdata1_exp=np.mean(data_exp, axis=0) # [lvls] unit
      vdata1_ref=np.mean(data_ref, axis=0) # [lvls] unit
   if var == 'velocity_potential' :
      vdata2_exp=np.mean(data_exp, axis=0) # [lvls] unit
      vdata2_ref=np.mean(data_ref, axis=0) # [lvls] unit
   if var == 'temperature' :
      vdata3_exp=np.mean(data_exp, axis=0) # [lvls] unit
      vdata3_ref=np.mean(data_ref, axis=0) # [lvls] unit
   if var == 'spechum' :
      vdata4_exp=np.mean(data_exp, axis=0) # [lvls] unit
      vdata4_ref=np.mean(data_ref, axis=0) # [lvls] unit
   np.delete(data_exp,0)
   np.delete(data_ref,0)

fig, (ax1,ax2,ax3,ax4) = plt.subplots(1,4,figsize=(10,2.5))
plt.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.9, wspace=None, hspace=None)
ax1.set(ylim=(0,max(z_havg)))
ax1.plot(vdata1_exp, z_havg)
ax1.plot(vdata1_ref, z_havg)
ax1.set(xlabel='lvls', ylabel='Vert. Coord [km]', title='psi')
ax1.grid()
ax2.set(ylim=(0,max(z_havg)))
ax2.plot(vdata2_exp, z_havg)
ax2.plot(vdata2_ref, z_havg)
ax2.set(xlabel='lvls',title='chi_u')
ax2.grid()
ax3.set(ylim=(0,max(z_havg)))
ax3.plot(vdata3_exp, z_havg)
ax3.plot(vdata3_ref, z_havg)
ax3.set(xlabel='lvls', title='T_u')
ax3.grid()
ax4.set(ylim=(0,max(z_havg)))
ax4.plot(vdata4_exp, z_havg)
ax4.plot(vdata4_ref, z_havg)
ax4.set(xlabel='lvls', title='Q')
ax4.grid()
#ax4.legend(['anisotropic','develop'])
ax4.legend(['exp','ref'])

plt.savefig(fname="zz_cmpr_CMAT_havg_rv.png")

f.close()
