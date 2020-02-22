import datetime as dt
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import matplotlib
matplotlib.use('AGG')
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import BoundaryNorm
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes
import numpy as np
import plot_utils as pu
import var_utils as vu
import os

#This script includes basic plotting functions.

def plotDistri(lats,lons,values,ObsType,VarName,var_unit,out_name,nstation,levbin):
#================================================================
#INPUTS:
# lats     - latitude
# lons     - longitude
# values   - values will be plotted
# ObsType - observation type
# VarName - variable name
# var_unit - variable units
# out_name - will be included in output file name. It can be experiment name.
# nstation - station numbers for sondes.
# levbin   - plot all levels together (levbin=all); or plot every level.
#================================================================

#set map=======================================================================
    fig,ax=plt.subplots(figsize=(8,8))
    m=Basemap(projection='cyl',     #map projection
        llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90, #the lat lon of leftlow and upper right corner
#       llcrnrlon=-20,llcrnrlat=10,urcrnrlon=60,urcrnrlat=80,   #zoom in
        resolution='l') #c: crude; l:low; i:intermediate, h:high, f:full; sometimes can only >=h
    m.drawcoastlines()  #draw coastline
#   m.drawmapboundary(fill_color='lightblue') #the whole map is filled with the specified color
#   m.fillcontinents(color='wheat',lake_color='lightblue') #the continents are filled

#set the lat lon grid and the ticks of x y axies ==============================
    parallels=np.arange(-90,90,30)      #lat
    meridians=np.arange(-180,180,60)    #lon
    m.drawparallels(parallels,
        labels=[True,True,True,True])   #Turn the ticks on or off
    m.drawmeridians(meridians,
        labels=[True,True,True,True])

#set title  ===================================================================
    if nstation == 0:
        plt.text(0.5, 1.25, '%s   variable: %s %s nlocs:%s' \
            %(ObsType,VarName,var_unit,len(values)),    \
            horizontalalignment='center', \
            fontsize=12, transform = ax.transAxes)
    else:
        plt.text(0.5, 1.25, '%s   variable: %s %s nlocs:%s nstation:%s' \
            %(ObsType,VarName,var_unit,len(values),nstation),    \
            horizontalalignment='center', \
            fontsize=12, transform = ax.transAxes)

#draw points onto map =========================================================
    cm=plt.cm.get_cmap('rainbow')
    sc=m.scatter(lons[:],lats[:],c=values[:],s=6,cmap=cm,
        zorder=10) #10 ; zorder determines the order of the layer. If not set,
                        #the point on continent will be blocked

#create axes for colorbar======================================================
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top",
        size="3%", #width of the colorbar
        pad=0.3)   #space between colorbar and graph
    plt.colorbar(sc,cax=cax,ax=ax,orientation='horizontal')
    cax.xaxis.set_ticks_position('top')  #set the orientation of the ticks of the colorbar

    plt.savefig('distri_%s_%s_%s.png'%(VarName,out_name,levbin),dpi=200,bbox_inches='tight')
    plt.close()

def plotTimeserial2D(Stats,xlabeltime,ylevels,VarName):
#================================================================
#INPUTS:
# Stats      - statistics
# xlabeltime - time labels for x-axis 
# ylevels    - vertical levels for y-axis
# VarName    - variable name
#================================================================
    zgrid = np.loadtxt("/glade/work/jban/pandac/fix_input/graphics/zgrid_v55.txt")

    fig, ax1 = plt.subplots()

    xarray = range(len(xlabeltime))
    valuemin = np.amin(Stats)
    valuemax = np.amax(Stats)

    if (valuemin > 0 or valuemax < 0):
        color = 'rainbow'
        plt.contourf(xarray,ylevels,Stats,40,vmin=valuemin, vmax=valuemax,cmap=color)
    else:
        cmap = 'coolwarm'
        norm = matplotlib.colors.DivergingNorm(vmin=valuemin, vcenter=0, vmax=valuemax)
        plt.contourf(xarray,ylevels,Stats,40,vmin=valuemin, vmax=valuemax,norm=norm,cmap=cmap)

    xarray = range(len(xlabeltime))
    major_ticks = np.arange(0, 56, 5)
    ax1.set_yticks(major_ticks)
    ax1.set_ylim([0,54])
    ax1.set_ylabel('Vertical level',fontsize=15)

    ax2 = ax1.twinx()
    ax2.set_yticks(major_ticks-1)
    ax2.set_yticklabels((zgrid[::5]).astype(int))

    ax2.set_ylabel('Height (m)',fontsize=13)

    FCDay = ''.join(VarName.split("_")[1:][:-3])
    if (FCDay == 'day0.0'):
        ax1.set_xlabel('Analysis Time',fontsize=15)
        ax1.set_xticks(xarray[::4])
        ax1.set_xticklabels(xlabeltime[::4],rotation=90)
    elif (FCDay == 'day0.25'):
        ax1.set_xlabel( '6h Forecast',fontsize=15)
        ax1.set_xticks(xarray[::4])
        ax1.set_xticklabels(xlabeltime[::4],rotation=90)
    else:
        ax1.set_xlabel( 'Lead Time',fontsize=15)

    plt.colorbar(extend='both',orientation="horizontal",pad=0.2)
    ax1.grid(True)
    region = ''.join(VarName.split("_")[2:][:-2])
    var    = ''.join(VarName.split("_")[3:][:-1])
    stats  = ''.join(VarName.split("_")[4:])
    plt.title(stats+'  variable:'+vu.varDictModel[var][1]+'('+ vu.varDictModel[var][0]+')  '+region, fontsize = 12)
    plt.savefig(VarName+'_TS_2d.png',dpi=200,bbox_inches='tight')
    plt.close()


###############################################################################
lenWarnSer = 0
nanWarnSer = 0
def plotSeries(fig, \
               linesVals, xVals, \
               linesLabel, \
               title="", dataLabel="y", \
               sciticks=False, signdef=False, \
               indepLabel="x", invert_ind_axis=False, \
               ny=1, nx=1, nplots=1, iplot=0, \
               linesValsMinCI=None, linesValsMaxCI=None, \
               dmin=np.NaN, dmax=np.NaN, \
               lineAttribOffset=0, \
               legend_inside=True,
               interiorLabels=True):

# ARGUMENTS
# fig              - matplotlib figure object
# linesVals        - dependent variable (list of arrays)
# xVals            - independent variable on x-axis (array)
# linesLabel       - legend label for linesVals (list)

# title            - subplot title, optional
# dataLabel        - label for linesVals, optional
# sciticks         - whether linesVals needs scientific formatting for ticks, optional
# signdef          - whether linesVals is positive/negative definite, optional
# indepLabel       - label for xVals, optional
# invert_ind_axis  - whether to invert x-axis orientation, optional

# ny, nx           - number of subplots in x/y direction, optional
# nplots           - total number of subplots, optional
# iplot            - this subplot index (starting at 0), optional

# linesValsMinCI   - minimum error bound for linesVals (list of arrays), optional
# linesValsMaxCI   - maximum error bound for linesVals (list of arrays), optional
# Note: linesValsMinCI and linesValsMaxCI must be specified together

# lineAttribOffset - offset for selecting line attributes, optional
# dmin, dmax       - min/max values of linesVals, optional
# legend_inside    - whether legend should be placed inside the subplot, optional

    ax = fig.add_subplot(ny, nx, iplot+1)

    #title
    ax.set_title(title,fontsize=5)

    #add lines
    plotVals = []
    nLines = 0
    for iline, lineVals in enumerate(linesVals):
        if all(np.isnan(lineVals)):
            global nanWarnSer
            if nanWarnSer==0:
                print("\nWARNING: skipping all-NaN data")
                print(title,indepLabel,linesLabel[iline])
            nanWarnSer=nanWarnSer+1
            continue
        if len(lineVals)!=len(xVals):
            global lenWarnSer
            if lenWarnSer==0:
                print("\nWARNING: skipping data where len(x)!=len(y)")
                print(title,indepLabel,linesLabel[iline])
            lenWarnSer=lenWarnSer+1
            continue

        # Plot line for each lineVals that has non-missing data
        pColor = pu.plotColors[iline+lineAttribOffset]

        ax.plot(xVals, lineVals, \
                color=pColor, \
                label=linesLabel[iline], \
                ls=pu.plotLineStyles[iline+lineAttribOffset], \
                linewidth=0.5)
        nLines = nLines + 1
        plotVals.append(lineVals)

        # Add shaded error regions if specified
        if linesValsMinCI is not None and \
           linesValsMaxCI is not None:

            # test statistical significance versus zero
            if signdef:
                significant = np.empty(len(lineVals))
                significant[:] = np.NaN
            else:
                significant = np.multiply(linesValsMinCI[iline], linesValsMaxCI[iline])
            significant = np.array([x if not np.isnan(x) else -1.0 for x in significant])

            siginds = np.array([i for i,x in enumerate(significant) if x > 0.0],dtype=int)
            if len(siginds) > 0:
                ax.plot(np.array(xVals)[siginds], np.array(lineVals)[siginds], \
                        color=pColor, \
                        ls='', \
                        marker=pu.plotMarkers[iline+lineAttribOffset], \
                        markersize=1.5)
            ax.plot(xVals, linesValsMinCI[iline], \
                    color=pColor, \
                    alpha=0.4, \
                    ls='-', \
                    linewidth=0.5)
            ax.plot(xVals, linesValsMaxCI[iline], \
                    color=pColor, \
                    alpha=0.4, \
                    ls='-', \
                    linewidth=0.5)
            ax.fill_between(xVals, linesValsMinCI[iline], linesValsMaxCI[iline], \
                            color=pColor, \
                            edgecolor=pColor, \
                            linewidth=0.0, alpha = 0.1)
            ax.fill_between(xVals, linesValsMinCI[iline], linesValsMaxCI[iline], \
                            where=significant > 0.0, \
                            color=pColor, \
                            edgecolor=pColor, \
                            linewidth=0.2, alpha = 0.3)

    if nLines == 0:
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    #axes settings
    ax.xaxis.set_tick_params(labelsize=3)
    ax.yaxis.set_tick_params(labelsize=3)

    if sciticks:
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.yaxis.get_offset_text().set_fontsize(3)

    # add vertical zero line for unbounded quantities
    if not signdef:
        ax.plot([xVals[0], xVals[-1]], [0., 0.], ls="--", c=".3", \
            linewidth=0.7,markersize=0)

    # standardize x-limits
    mindval, maxdval = pu.get_clean_ax_limits(dmin,dmax,plotVals,signdef)
    if (not np.isnan(mindval) and
        not np.isnan(maxdval) and
        not np.isinf(mindval) and
        not np.isinf(maxdval)):
        ax.set_ylim(mindval,maxdval)
        if maxdval-mindval < 1.0 or \
           maxdval-mindval > 100.0:
            ax.tick_params(axis='y',rotation=-35)

    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if not interiorLabels \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if interiorLabels or ix == 0:
        ax.set_xlabel(indepLabel,fontsize=4)
        ax.set_ylabel(dataLabel,fontsize=4)

    #legend
    if legend_inside:
        #INSIDE AXES
        lh = ax.legend(loc='best',fontsize=3,frameon=True,\
                       framealpha=0.4,ncol=1)
        lh.get_frame().set_linewidth(0.0)
    elif ix==nx-1 or iplot==nplots-1:
        #OUTSIDE AXES
        ax.legend(loc='upper left',fontsize=3,frameon=False, \
                  bbox_to_anchor=(1.02, 1), borderaxespad=0)

    if invert_ind_axis:
        ax.invert_xaxis()

    ax.grid()

    return

###############################################################################
lenWarnProf = 0
nanWarnProf = 0
def plotProfile(fig, \
                linesVals, yVals, \
                linesLabel, \
                title="", dataLabel="x", \
                sciticks=False, signdef=False, \
                indepLabel="y", invert_ind_axis=False, \
                ny=1, nx=1, nplots=1, iplot=0, \
                linesValsMinCI=None, linesValsMaxCI=None, \
                dmin=np.NaN, dmax=np.NaN, \
                lineAttribOffset=0, \
                legend_inside=True,
                interiorLabels=True):

# ARGUMENTS
# fig              - matplotlib figure object
# linesVals        - dependent variable (list of arrays)
# yVals            - independent variable on y-axis (array)
# linesLabel       - legend label for linesVals (list)

# title            - subplot title, optional
# dataLabel           - label for linesVals, optional
# sciticks         - whether linesVals needs scientific formatting for ticks, optional
# signdef          - whether linesVals is positive/negative definite, optional
# indepLabel       - label for yVals, optional
# invert_ind_axis  - whether to invert y-axis orientation, optional

# ny, nx           - number of subplots in x/y direction, optional
# nplots           - total number of subplots, optional
# iplot            - this subplot index (starting at 0), optional

# linesValsMinCI   - minimum error bound for linesVals (list of arrays), optional
# linesValsMaxCI   - maximum error bound for linesVals (list of arrays), optional
# Note: linesValsMinCI and linesValsMaxCI must be specified together

# lineAttribOffset - offset for selecting line attributes, optional
# dmin, dmax       - min/max values of linesVals, optional
# legend_inside    - whether legend should be placed inside the subplot, optional

    ax = fig.add_subplot(ny, nx, iplot+1)

    #title
    ax.set_title(title,fontsize=5)

    #add lines
    plotVals = []
    nLines = 0
    for iline, lineVals in enumerate(linesVals):
        if all(np.isnan(lineVals)):
            global nanWarnProf
            if nanWarnProf==0:
                print("\nWARNING: skipping all-NaN data")
                print(title,dataLabel,linesLabel[iline])
            nanWarnProf=nanWarnProf+1
            continue
        if len(lineVals)!=len(yVals):
            global lenWarnProf
            if lenWarnProf==0:
                print("\nWARNING: skipping data where len(x)!=len(y)")
                print(title,dataLabel,linesLabel[iline])
            lenWarnProf=lenWarnProf+1
            continue

        # Plot line for each lineVals that has non-missing data
        pColor = pu.plotColors[iline+lineAttribOffset]

        ax.plot(lineVals, yVals, \
                color=pColor, \
                label=linesLabel[iline], \
                ls=pu.plotLineStyles[iline+lineAttribOffset], \
                linewidth=0.5)
        nLines = nLines + 1
        plotVals.append(lineVals)

        # Add shaded error regions if specified
        if linesValsMinCI is not None and \
           linesValsMaxCI is not None:

            # test statistical significance versus zero
            if signdef:
                significant = np.empty(len(lineVals))
                significant[:] = np.NaN
            else:
                significant = np.multiply(linesValsMinCI[iline], linesValsMaxCI[iline])
            significant = np.array([x if not np.isnan(x) else -1.0 for x in significant])

            siginds = np.array([i for i,x in enumerate(significant) if x > 0.0],dtype=int)
            if len(siginds) > 0:
                ax.plot(np.array(lineVals)[siginds], np.array(yVals)[siginds], \
                        color=pColor, \
                        ls='', \
                        marker=pu.plotMarkers[iline+lineAttribOffset], \
                        markersize=1.5)
            ax.plot(linesValsMinCI[iline], yVals, \
                    color=pColor, \
                    alpha=0.4, \
                    ls='-', \
                    linewidth=0.5)
            ax.plot(linesValsMaxCI[iline], yVals, \
                    color=pColor, \
                    alpha=0.4, \
                    ls='-', \
                    linewidth=0.5)
            ax.fill_betweenx(yVals, linesValsMinCI[iline], linesValsMaxCI[iline], \
                            color=pColor, \
                            edgecolor=pColor, \
                            linewidth=0.0, alpha = 0.1)
            ax.fill_betweenx(yVals, linesValsMinCI[iline], linesValsMaxCI[iline], \
                            where=significant > 0.0, \
                            color=pColor, \
                            edgecolor=pColor, \
                            linewidth=0.2, alpha = 0.3)

    if nLines == 0:
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    #axes settings
    ax.xaxis.set_tick_params(labelsize=3)
    ax.yaxis.set_tick_params(labelsize=3)

    if sciticks:
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax.xaxis.get_offset_text().set_fontsize(3)

    # add vertical zero line for unbounded quantities
    if not signdef:
        ax.plot([0., 0.], [yVals[0], yVals[-1]], ls="--", c=".3", \
            linewidth=0.7,markersize=0)

    # standardize x-limits
    mindval, maxdval = pu.get_clean_ax_limits(dmin,dmax,plotVals,signdef)
    if (not np.isnan(mindval) and
        not np.isnan(maxdval) and
        not np.isinf(mindval) and
        not np.isinf(maxdval)):
        ax.set_xlim(mindval,maxdval)
        if maxdval-mindval < 1.0 or \
           maxdval-mindval > 100.0:
            ax.tick_params(axis='x',rotation=-35)

    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if not interiorLabels \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if interiorLabels or ix == 0:
        ax.set_xlabel(dataLabel,fontsize=4)
        ax.set_ylabel(indepLabel,fontsize=4)

    #legend
    if legend_inside:
        #INSIDE AXES
        lh = ax.legend(loc='best',fontsize=3,frameon=True,\
                       framealpha=0.4,ncol=1)
        lh.get_frame().set_linewidth(0.0)
    elif ix==nx-1 or iplot==nplots-1:
        #OUTSIDE AXES
        ax.legend(loc='upper left',fontsize=3,frameon=False, \
                  bbox_to_anchor=(1.02, 1), borderaxespad=0)

    if invert_ind_axis:
        ax.invert_yaxis()

    ax.grid()

    return


###############################################################################
lenWarnTS=0
nanWarnTS=0
def plotTimeSeries(fig, \
                   xsDates, linesVals, \
                   linesLabel, \
                   title="", ylabel="", \
                   sciticks=False, signdef=False, \
                   ny=1, nx=1, nplots=1, iplot=0, \
                   linesValsMinCI=None, linesValsMaxCI=None, \
                   dmin=np.NaN, dmax=np.NaN, \
                   lineAttribOffset=0, \
                   legend_inside=True,
                   interiorLabels=True):

# ARGUMENTS
# fig              - matplotlib figure object
# xsDates          - date x-values (list/array or list of lists/arrays
#                                   of float seconds, dt.timedelta, dt.datetime)
# linesVals        - dependent variable (list of arrays)
# linesLabel       - legend label for linesVals (list)

# title            - subplot title, optional
# ylabel           - label for linesVals, optional
# sciticks         - whether linesVals needs scientific formatting for ticks, optional
# signdef          - whether linesVals is positive/negative definite, optional

# ny, nx           - number of subplots in x/y direction, optional
# nplots           - total number of subplots, optional
# iplot            - this subplot index (starting at 0), optional

# linesValsMinCI   - minimum error bound for linesVals (list of arrays), optional
# linesValsMaxCI   - maximum error bound for linesVals (list of arrays), optional
# Note: linesValsMinCI and linesValsMaxCI must be specified together

# lineAttribOffset - offset for selecting line attributes, optional
# dmin, dmax       - min/max values of linesVals, optional
# legend_inside    - whether legend should be placed inside the subplot, optional

    ax = fig.add_subplot(ny, nx, iplot+1)

    #title
    ax.set_title(title,fontsize=5)

    #add lines
    plotVals = []
    nLines = 0
    for iline, lineVals in enumerate(linesVals):
        if all(np.isnan(lineVals)):
            global nanWarnTS
            if nanWarnTS==0:
                print("\nWARNING: skipping all-NaN data")
                print(title,ylabel,linesLabel[iline])
            nanWarnTS=nanWarnTS+1
            continue

        #float xVals
        if isinstance(xsDates[0],(list,np.ndarray)):
            xVals = pu.TDeltas2Seconds(xsDates[min([iline,len(xsDates)-1])])
        else:
            xVals = pu.TDeltas2Seconds(xsDates)

        if len(lineVals)!=len(xVals):
            global lenWarnTS
            if lenWarnTS==0:
                print("\nWARNING: skipping data where len(x)!=len(y)")
                print(title,ylabel,linesLabel[iline])
            lenWarnTS=lenWarnTS+1
            continue

        if iline == 0:
            minX = xVals[0]
            maxX = xVals[-1]
        else:
            minX = min([xVals[0], minX])
            maxX = max([xVals[-1], maxX])

        # Plot line for each lineVals that has non-missing data
        pColor = pu.plotColors[iline+lineAttribOffset]

        ax.plot(xVals, lineVals, \
                label=linesLabel[iline], \
                color=pColor, \
                ls=pu.plotLineStyles[iline+lineAttribOffset], \
                linewidth=0.5)
        nLines = nLines + 1
        plotVals.append(lineVals)

        # Add shaded CI regions if specified
        if linesValsMinCI is not None and \
           linesValsMaxCI is not None:

            # test statistical significance versus zero
            if signdef:
                significant = np.empty(len(lineVals))
                significant[:] = np.NaN
            else:
                significant = np.multiply(linesValsMinCI[iline], linesValsMaxCI[iline])
            significant = np.array([x if not np.isnan(x) else -1.0 for x in significant])
            siginds = np.array([i for i,x in enumerate(significant) if x > 0.0],dtype=int)
            if len(siginds) > 0:
                ax.plot(np.array(xVals)[siginds], np.array(lineVals)[siginds], \
                        color=pColor, \
                        ls='', \
                        marker=pu.plotMarkers[iline+lineAttribOffset], \
                        markersize=1.5)
            ax.plot(xVals, linesValsMinCI[iline], \
                    color=pColor, \
                    alpha=0.4, \
                    ls='-', \
                    linewidth=0.5)
            ax.plot(xVals, linesValsMaxCI[iline], \
                    color=pColor, \
                    alpha=0.4, \
                    ls='-', \
                    linewidth=0.5)
            ax.fill_between(xVals, linesValsMinCI[iline], linesValsMaxCI[iline], \
                            color=pColor, \
                            edgecolor=pColor, \
                            linewidth=0.0, alpha = 0.1)
            ax.fill_between(xVals, linesValsMinCI[iline], linesValsMaxCI[iline], \
                            where=significant > 0.0, \
                            color=pColor, \
                            edgecolor=pColor, \
                            linewidth=0.2, alpha = 0.3)

    if nLines == 0:
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    #axes settings
    if isinstance(xsDates[0],(list,np.ndarray)):
        pu.format_x_for_dates(ax, xsDates[0])
    else:
        pu.format_x_for_dates(ax, xsDates)
    ax.xaxis.set_tick_params(labelsize=3)
    ax.yaxis.set_tick_params(labelsize=3)

    if sciticks:
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.yaxis.get_offset_text().set_fontsize(3)

    # add horizontal zero line for unbounded quantities
    if not signdef:
        ax.plot([minX, maxX], [0., 0.], ls="--", c=".3", \
            linewidth=0.7,markersize=0)

    # standardize y-limits
    mindval, maxdval = pu.get_clean_ax_limits(dmin,dmax,plotVals,signdef)
    if (not np.isnan(mindval) and
        not np.isnan(maxdval) and
        not np.isinf(mindval) and
        not np.isinf(maxdval)):
        ax.set_ylim(mindval,maxdval)

    ax.grid()

    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if not interiorLabels \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if interiorLabels or ix == 0:
        ax.set_ylabel(ylabel,fontsize=4)

    #legend
    if legend_inside:
        #INSIDE AXES
        nlcol = np.int(np.ceil(np.sqrt(nLines)))
        lh = ax.legend(loc='best',fontsize=3,frameon=True,\
                       framealpha=0.4,ncol=nlcol)
        lh.get_frame().set_linewidth(0.0)
    elif ix==nx-1 or iplot==nplots-1:
        #OUTSIDE AXES
        ax.legend(loc='upper left',fontsize=3,frameon=False, \
                  bbox_to_anchor=(1.02, 1), borderaxespad=0)


    return


###############################################################################
def plotTimeSeries2D(fig, \
                     xDates, yVals, contourVals, \
                     title="", clabel="", \
                     sciticks=False, signdef=False, \
                     ylabel="y", invert_ind_axis=False, \
                     ny=1, nx=1, nplots=1, iplot=0, \
                     dmin=np.NaN, dmax=np.NaN,
                     interiorLabels=True):

# ARGUMENTS
# fig           - matplotlib figure object
# xDates        - date x-values (array of float seconds, dt.timedelta, dt.datetime)
# yVals         - second independent variable
# contourVals   - dependent variable (2d array)

# title         - subplot title, optional
# clabel        - label for dependent variable, optional
# sciticks      - whether contourVals needs scientific formatting for ticks, optional
# signdef       - whether contourVals is positive/negative definite, optional
# ylabel        - label for yVals, optional
# invert_ind_axis - whether to invert y-axis orientation, optional

# ny, nx        - number of subplots in x/y direction, optional
# nplots        - total number of subplots, optional
# iplot         - this subplot index (starting at 0), optional

# dmin, dmax    - min/max values of contourVals, optional

    ax = fig.add_subplot(ny, nx, iplot+1)

    if (np.isnan(contourVals)).all():
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    xVals = pu.TDeltas2Seconds(xDates)

    # standardize c-limits
    mindval, maxdval = pu.get_clean_ax_limits(dmin,dmax,contourVals,signdef)
    if signdef:
        cmapName = 'BuPu'
        nlevs = 18
    else:
        cmapName = 'seismic'
        nlevs = 28

    # plot contour
    # option 1: smoothed contours
    #cp = ax.contourf(xVals, yVals, contourVals, nlevs, cmap=cmapName, extend='both', \
    #     vmin=mindval, vmax=maxdval)

    # option 2: pixel contours
    levels = MaxNLocator(nbins=nlevs).tick_values(mindval,maxdval)
    cmap = plt.get_cmap(cmapName)
    cmap.set_bad(color = 'k', alpha = 1.0)
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    xVals_pcolor, yVals_pcolor = transformXY_for_pcolor(xVals,yVals)
    cp = ax.pcolormesh(xVals_pcolor, yVals_pcolor, contourVals, cmap=cmap, norm = norm)

    #title
    ax.set_title(title,fontsize=5)

    #axes settings
    pu.format_x_for_dates(ax, xDates)
    ax.xaxis.set_tick_params(labelsize=3)
    ax.yaxis.set_tick_params(labelsize=3)

    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if not interiorLabels \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if interiorLabels or ix == 0:
        ax.set_ylabel(ylabel,fontsize=4)
    if interiorLabels or ix == nx-1:
        #colorbar
        m = plt.cm.ScalarMappable(cmap=cmap)
        m.set_array(contourVals)
        if (not np.isnan(mindval) and
            not np.isnan(maxdval) and
            not np.isinf(mindval) and
            not np.isinf(maxdval)):
            m.set_clim(mindval,maxdval)
        cb = plt.colorbar(m, ax=ax)
        #scientific formatting
        if sciticks:
            cb.ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            cb.ax.yaxis.get_offset_text().set_fontsize(3)

        cb.ax.tick_params(labelsize=3)
        cb.set_label(clabel,fontsize=5)

    if invert_ind_axis:
        ax.invert_yaxis()

    # optionally add a grid
    #ax.grid()

    return


###############################################################################
def transformXY_for_pcolor(xs,ys):
    # adjust centered x and y values to edges to work with pcolormesh 
    # note: works best for regularly spaced data
    xs_diff = xs[1] - xs[0]
    # extend xs by 2
    # fill in first endpoint
    xs_extend = [xs[0]-xs_diff]
    # fill in internal values
    for x in xs: xs_extend.append(x)
    # fill in last endpoint
    xs_extend.append(xs_extend[-1]+(xs[-1]-xs[-2]))
    # calculate the midpoints
    xs_pcolormesh_midpoints = []
    for ii, x in enumerate(xs_extend[:-1]):
        xs_pcolormesh_midpoints.append(x+0.5*(xs_extend[ii+1] - xs_extend[ii]))

    ys_diff = ys[1] - ys[0]
    # extend ys by 2
    # fill in first endpoint
    ys_extend = [ys[0]-ys_diff]
    # fill in internal values
    for y in ys: ys_extend.append(y)
    # fill in last endpoint
    ys_extend.append(ys_extend[-1]+(ys[-1]-ys[-2]))
    # calculate the midpoints
    ys_pcolormesh_midpoints = []
    for ii, y in enumerate(ys_extend[:-1]):
        ys_pcolormesh_midpoints.append(y+0.5*(ys_extend[ii+1] - ys_extend[ii]))

    return xs_pcolormesh_midpoints, ys_pcolormesh_midpoints

def main():
    print ('This is not a runnable program.')
    os._exit(0)

if __name__ == '__main__': main()
