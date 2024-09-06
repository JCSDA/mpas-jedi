import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as tck
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

import utils as ut

import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature

def line(stat=str,clist=list(),fhours=int,experiments=None,data1=None,data2=None,labelsL=list(),savename=str):

    all_panels = len(clist)
    nrows = int(np.sqrt(all_panels))
    ncols = all_panels // nrows
    d = plt.rcParams.get('figure.figsize') # [6.4, 4.8]
    figsize_ = d
    if ncols == 2:
      figsize_ = tuple(np.multiply(d,2.5))
    elif ncols == 4:
      figsize_ = (np.multiply(d,3.75))
    elif ncols == 3 and nrows == 1:
      figsize_ = (int(np.multiply(d,3)[0]),d[1])#2.8125
    else:
      figsize_ = d

    save_fig = True
    fig = plt.figure(figsize=figsize_)
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=list(np.ones(ncols,int)), height_ratios=list(np.ones(nrows,int)), wspace=0.25, hspace=0.3)

    panel_names = clist
    labels = ut.createDictL(experiments,labelsL)
    colors = ut.createDictL(experiments,ut.colorsL)

    for i in range(all_panels):
      axes = fig.add_subplot(gs[i])
      axes.xaxis.grid(linestyle=":", alpha=0.1, color='grey')
      axes.yaxis.grid(linestyle=":", alpha=0.1, color='grey')
      title = panel_names[i]
      if len(title.split('_')) > 1:
        title = title.split('_')[1]
      elif title in ut.cat_stats:
        title = ''
      axes.set_title(title, loc='center', fontsize=16,  fontweight='bold')
      axes.xaxis.set_tick_params(labelsize=12)
      axes.yaxis.set_tick_params(labelsize=12)
      axes.xaxis.set_minor_locator(tck.AutoMinorLocator(3))
      axes.yaxis.set_minor_locator(tck.AutoMinorLocator())

      ya = axes.get_yaxis()
      ya.set_major_locator(tck.MaxNLocator(integer=True))

      x = range(fhours+1)
      axes.set_xticks(x[::3])
      axes.set_xlabel('Forecast hour', fontsize=14)
      axes.set_ylabel(stat, fontsize=14)

      if stat == 'count':
        axes.ticklabel_format(axis='y',style='sci', scilimits=(0,4))

      for e in range(len(experiments)):
        if stat in ut.cont_stats:
          outputFname = 'continous'

          lStat = ut.get_label(data1[experiments[e]],stat,panel_names[i])
          axes.plot(data1[experiments[e]][panel_names[i]],label=labels[experiments[e]]+lStat,c=colors[experiments[e]],ls='-',lw=2)

          if data2 != None and e == 0: # add obs values only for obs
            if   'Ocloud-Mcloud' in panel_names[i]:
              lObs = ut.get_label(data2[experiments[e]],stat,'Ocloud')
              axes.plot(data2[experiments[e]]['Ocloud'],label='Obs cloudy'+lObs,c='grey',ls='--',lw = 3)
            elif 'Oclear-Mclear' in panel_names[i]:
              lObs = ut.get_label(data2[experiments[e]],stat,'Oclear')
              axes.plot(data2[experiments[e]]['Oclear'],label='Obs clear'+lObs,c='grey',ls='--',lw=3)
            elif 'Olow-Mlow' in panel_names[i]:
              lObs = ut.get_label(data2[experiments[e]],stat,'Olow')
              axes.plot(data2[experiments[e]]['Olow'],label='Obs low'+lObs,c='grey',ls='--',lw=3)
            elif 'Omid-Mmid' in panel_names[i]:
              lObs = ut.get_label(data2[experiments[e]],stat,'Omid')
              axes.plot(data2[experiments[e]]['Omid'],label='Obs mid'+lObs,c='grey',ls='--',lw=3)
            elif 'Ohigh-Mhigh' in panel_names[i]:
              lObs = ut.get_label(data2[experiments[e]],stat,'Ohigh')
              axes.plot(data2[experiments[e]]['Ohigh'],label='Obs high'+lObs,c='grey',ls='--',lw=3)
            elif 'all' in panel_names[i]:
              lObs = ut.get_label(data2[experiments[e]],stat,'all')
              axes.plot(data2[experiments[e]]['all'],label='Obs all'+lObs,c='grey',ls='--',lw=3)
        else:
          outputFname = 'categorical'
          layer = clist[0]
          if len(data1[experiments[e]][stat][layer]) > 0:
            axes.plot(data1[experiments[e]][stat][layer],label=labels[experiments[e]],c=colors[experiments[e]],ls='-',lw=2)
          else:
            save_fig = False

      if save_fig:
        axes.legend()

    outputfolder = './'+outputFname
    if not os.path.exists(outputfolder):
      os.makedirs(outputfolder)

    if save_fig:
      plt.savefig(outputfolder+'/'+savename+'.png',dpi=300,bbox_inches='tight')
    plt.close()


def pdf(data1=None,data2=None,title=None,savename=None):

    hist_v1, bins_v1, patches_v1 = plt.hist(data1, bins='auto', density=True)
    hist_v2, bins_v2, patches_v2 = plt.hist(data2, bins='auto', density=True)
    res_v1 = ut.log10(hist_v1)
    res_v2 = ut.log10(hist_v2)

    fig = plt.figure()
    gs = GridSpec(ncols=1, nrows=1, height_ratios=[1], wspace=0, hspace=0)
    ax = fig.add_subplot(gs[0])
    ax.plot(bins_v1[1:],res_v1, label='Forecast ('+str(len(data1))+')', c='k')
    ax.plot(bins_v2[1:],res_v2, label='Observations ('+str(len(data2))+')', c='r')
    ax.set_xlabel('Brightness temperature [K]')
    ax.set_ylabel('log10(PDF)')
    ax.set_title(title, loc='center', fontsize=16,  fontweight='bold')

    plt.legend()

    outputfolder = './pdf'
    if not os.path.exists(outputfolder):
      os.makedirs(outputfolder)

    plt.savefig(outputfolder+'/'+savename+'.png',dpi=300,bbox_inches='tight')
    plt.close()

def scatter(lon,lat,data,title,colormap,savename):
    proj = ccrs.PlateCarree(central_longitude=180)
    lon = (lon + 180) % 360 - 180
    extent = [-140,-10,-65,65]
    fig = plt.figure(figsize=(8,8))
    ax = plt.axes(projection=proj)
    background(ax, extent)
    vmin = np.nanmin(data)
    vmax = np.nanmax(data)
    cmap = colormap

    cntr = ax.scatter(lon,lat, c=data, s=0.5, vmin=vmin, vmax=vmax, cmap=cmap, transform=ccrs.PlateCarree())
    plt.title( title+', nlocs: '+str(len(data)) )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom",size="3%", pad=0.2,axes_class=plt.Axes)
    cbar = plt.colorbar(cntr,cax=cax,orientation='horizontal',extend='both')

    outputfolder = './scatter'
    if not os.path.exists(outputfolder):
      os.makedirs(outputfolder)

    plt.savefig(outputfolder+'/'+savename+'.png',dpi=300,bbox_inches='tight')
    plt.close()


def background(ax,extent):
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='black', alpha=0.5, linestyle='dotted')
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = True
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}
    return ax
