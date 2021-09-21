#!/usr/bin/env python3

#import datetime as dt
#import getopt
from copy import deepcopy
import matplotlib
import matplotlib.pyplot as plt
#from matplotlib.dates import DateFormatter, AutoDateLocator
import numpy as np
import os
import pandas as pd
#import sys

import plot_utils as pu
import stat_utils as su


def bootstrap_test():

    NOMINAL_GROUP_SIZE = 500

    NOMINAL_NUMBER_OF_GROUPS = 60


    N_SUPER_SAMPLE = NOMINAL_GROUP_SIZE * NOMINAL_NUMBER_OF_GROUPS

    # number of bootstrap samples
    nBootSamples = np.array([  np.around(10.0 ** x).astype(int) \
                               for x in list(np.arange(1.5,4.5,0.25))])

#    nBootSamples = np.array([  np.around(10.0 ** x).astype(int) \
#                               for x in list(np.arange(1.5,2.75,0.25))])


    # general settings for plotting
    fullColor = pu.plotColors[1]
    clustColor = pu.plotColors[2]

    # establish f(g(X),g(Y)) for full population sampling
    #g(x)
    vecFuncs = [su.identityFunc, np.square, su.rmsFunc]
    bootFuncs = [su.bootStrapVector,su.bootStrapVector,su.bootStrapVectorRMSFunc]
    # vecFuncs = [np.square]
    nFuncs = len(vecFuncs)

    #f(x)
    statFunc = np.subtract

    # draw data from N(XMEAN, XVAR):
    randSample1 = np.random.randn(N_SUPER_SAMPLE)
    X1MEAN = 2.0   ; X1STD = 1.5
    X1_SAMPLE = X1MEAN + X1STD*randSample1

    randSample2 = np.random.randn(N_SUPER_SAMPLE)
    X2MEAN = 1.85 ; X2STD = 1.65
    X2_SAMPLE = X2MEAN + X2STD*randSample2

    print("\n\nX1_SAMPLE stats:")
    print(su.calcStats(X1_SAMPLE))

    print("\n\nX2_SAMPLE stats:")
    print(su.calcStats(X2_SAMPLE))

    print("\n\nstatFunc(vecFunc(X2),vecFunc(X1)) stats:")
    sampleStats = []
    for ifunc, vecFunc in enumerate(vecFuncs):
        if vecFunc is su.rmsFunc:
            sampleStats.append({'Mean': statFunc(vecFunc(X2_SAMPLE),vecFunc(X1_SAMPLE))})

        else:
            sampleStats.append(su.calcStats(statFunc(vecFunc(X2_SAMPLE),vecFunc(X1_SAMPLE))))

    print(sampleStats)

    # Bootstrap on mean of point-by-point statistics for the full population
    # MEAN(statFunc(vecFunc(X2),vecFunc(X1)))

    # Bootstrap on full population
    statsCIFullSamples = su.bootStrapVectorFunc( \
                             X2_SAMPLE, X1_SAMPLE, 
                             n_samples = nBootSamples, \
                             vecFuncs = vecFuncs, \
                             bootFuncs = bootFuncs, \
                             statFunc = statFunc)

    # Loop over variations in cluster sizes
#    SIGMANS = [5, 15, 50, 150, 400]
    SIGMANS = np.around(np.multiply(NOMINAL_GROUP_SIZE, \
                  [0.005, 0.015, 0.050, 0.150, 0.400])).astype(int)


    for SIGMAN in SIGMANS:

        #Create DataFrames containing cluster statistics for both experiments
        # establish clusters
        clustsN = []
        clustsStart = []
        clustsEnd = []

## alternative way to initialize cluster sizes
#        clustsStart.append(0)
#        nremain = N_SUPER_SAMPLE
#        while nremain > 0:
#            clustN = int(max([np.around(NOMINAL_GROUP_SIZE + np.random.randn(1)*SIGMAN), 1.0]))
#            if np.sum(clustsN) + clustN > N_SUPER_SAMPLE - NOMINAL_GROUP_SIZE/4:
#                clustN = N_SUPER_SAMPLE - np.sum(clustsN)
#            clustsN.append( clustN )
#            nremain = nremain - clustN
#            ntotal = np.sum(clustsN)
#            clustsEnd.append(ntotal)
#            clustsStart.append(clustsEnd[-1])
#

        clustsN = np.multiply(NOMINAL_GROUP_SIZE,np.ones(NOMINAL_NUMBER_OF_GROUPS,dtype=int))
        for i in list(range(0,int(NOMINAL_NUMBER_OF_GROUPS/2 - 1))):
            delta = NOMINAL_GROUP_SIZE - \
                    int(max([np.around(NOMINAL_GROUP_SIZE + np.random.randn(1)*SIGMAN), 1.0]))

            delta = np.abs(delta)

            ### note: always causes non-negligible difference in CI
            clustsN[i]      = np.max([clustsN[i]      + delta, 1])
            clustsN[-(i+1)] = np.max([clustsN[-(i+1)] - delta, 1])


            # clustsN[i]      = np.max([clustsN[i]      + delta*2, 1])
            # clustsN[i+1]    = np.max([clustsN[i+1]    - delta,   1])
            # clustsN[-(i+1)] = np.max([clustsN[-(i+1)] - delta,   1])

        clustsN = np.random.choice(clustsN, NOMINAL_NUMBER_OF_GROUPS, replace = False )

        clustsStart.append(0)

        for i in list(range(0,NOMINAL_NUMBER_OF_GROUPS)):
            clustsEnd.append(np.sum(clustsN[0:i+1]))
            clustsStart.append(clustsEnd[-1])

        print("\n\nclustsN = ",clustsN)


        #X1_SAMPLE
        X1ClustsValues = []
        X1ClustsStats = {}
        for stat in su.aggregatableFileStats:
            X1ClustsStats[stat] = []

        for ic, clusts in enumerate(clustsN):
            X1ClustsValues.append(X1_SAMPLE[clustsStart[ic]:clustsEnd[ic]])
            clustStats = su.calcStats(X1ClustsValues[ic]) 
            for stat in su.aggregatableFileStats:
                X1ClustsStats[stat].append(clustStats[stat])

        X1ClustsStatsDF = pd.DataFrame.from_dict(X1ClustsStats)
        X1ClustsStatsDF.sort_index()

        # print("\n\nX1 cluster stats:")
        # print(X1ClustsStatsDF)


        print("\n\nX1 aggregated stats:")
        aggClustStats = su.aggStatsDF(X1ClustsStatsDF)
        print(aggClustStats)


        #X2_SAMPLE
        X2ClustsValues = []
        X2ClustsStats = {}
        for stat in su.aggregatableFileStats:
            X2ClustsStats[stat] = []

        for ic, clusts in enumerate(clustsN):
            X2ClustsValues.append(X2_SAMPLE[clustsStart[ic]:clustsEnd[ic]])
            clustStats = su.calcStats(X2ClustsValues[ic]) 
            for stat in su.aggregatableFileStats:
                X2ClustsStats[stat].append(clustStats[stat])

        X2ClustsStatsDF = pd.DataFrame.from_dict(X2ClustsStats)
        X2ClustsStatsDF.sort_index()


        # print("\n\nX2 cluster stats:")
        # print(X2ClustsStatsDF)


        print("\n\nX2 aggregated stats:")
        aggClustStats = su.aggStatsDF(X2ClustsStatsDF)
        print(aggClustStats)


        print("\n\nPerforming bootstrap for sigmaN = ",SIGMAN)

        # Bootstrap on cluster subpopulations
        statsCIClustSamples = su.bootStrapClusterFunc( \
                                  X2ClustsStatsDF, X1ClustsStatsDF, \
                                  n_samples = nBootSamples, \
                                  statFunc = statFunc)

        #
        # Results for point-by-point
        #
        xVals = nBootSamples
        ny = nFuncs
        nx = 2

        fig = pu.setup_fig(nx, ny, inch_width=2.2)

        for ifunc, vecFunc in enumerate(vecFuncs):
            bootFunc = bootFuncs[ifunc]
            if vecFunc is su.identityFunc:
                sampleAggStat = "Mean"
            elif vecFunc is np.square:
                sampleAggStat = "MS"
            elif vecFunc is su.rmsFunc:
                sampleAggStat = "RMS"
            else:
                print("ERROR: vecFunc has no equivalent in aggregated cluster bootstrap")
                os._exit(1)
    
    
            ax1 = fig.add_subplot(ny, nx, nx*ifunc+1)

            plotVals = []

            #Full Sampling
            lineVals = deepcopy(statsCIFullSamples[ifunc]['VALUE'])
            lineValsMin = deepcopy(statsCIFullSamples[ifunc]['LO'])
            lineValsMax = deepcopy(statsCIFullSamples[ifunc]['HI'])
            plotVals.append(lineVals)
            plotVals.append(lineValsMin)
            plotVals.append(lineValsMax)


            ax1.plot(xVals, lineVals, \
                    label='std bootstrap', \
                    color=fullColor, \
                    ls=pu.plotLineStyles[1], \
                    linewidth=0.7)

            ax1.plot(xVals, lineValsMin, \
                    label=None, \
                    color=fullColor, \
                    alpha=0.7, \
                    ls='-', \
                    linewidth=0.7)
            ax1.plot(xVals, lineValsMax, \
                    label=None, \
                    color=fullColor, \
                    alpha=0.7, \
                    ls='-', \
                    linewidth=0.7)

            ax1.plot([xVals[0], xVals[-1]], np.multiply(lineValsMin[-1],[1., 1.]), \
                    ls="--", c=fullColor, \
                    alpha=0.4, \
                    linewidth=0.5,markersize=0)
            ax1.plot([xVals[0], xVals[-1]], np.multiply(lineValsMax[-1],[1., 1.]), \
                    ls="--", c=fullColor, \
                    alpha=0.4, \
                    linewidth=0.5,markersize=0)


            #Cluster Sampling
            lineVals = deepcopy(statsCIClustSamples[sampleAggStat]['VALUE'])
            lineValsMin = deepcopy(statsCIClustSamples[sampleAggStat]['LO'])
            lineValsMax = deepcopy(statsCIClustSamples[sampleAggStat]['HI'])
            plotVals.append(lineVals)
            plotVals.append(lineValsMin)
            plotVals.append(lineValsMax)

            ax1.plot(xVals, lineVals, \
                    label='cluster bootstrap', \
                    color=clustColor, \
                    ls=pu.plotLineStyles[1], \
                    linewidth=0.7)

            ax1.plot(xVals, lineValsMin, \
                    label=None, \
                    color=clustColor, \
                    alpha=0.6, \
                    ls='-', \
                    linewidth=0.7)
            ax1.plot(xVals, lineValsMax, \
                    label=None, \
                    color=clustColor, \
                    alpha=0.6, \
                    ls='-', \
                    linewidth=0.7)

            ax1.plot([xVals[0], xVals[-1]], np.multiply(lineValsMin[-1],[1., 1.]), \
                    ls="--", c=clustColor, \
                    alpha=0.4, \
                    linewidth=0.5,markersize=0)
            ax1.plot([xVals[0], xVals[-1]], np.multiply(lineValsMax[-1],[1., 1.]), \
                    ls="--", c=clustColor, \
                    alpha=0.4, \
                    linewidth=0.5,markersize=0)


            # Population Statistics
            ax1.plot([xVals[0], xVals[-1]], np.multiply(sampleStats[ifunc]['Mean'],[1., 1.]), \
                    ls="--", c=".2", \
                linewidth=0.9,markersize=0,label='sample pop.')

            lh = ax1.legend(loc='best',fontsize=3,frameon=True,\
                           framealpha=0.4,ncol=2)
            lh.get_frame().set_linewidth(0.0)

            ax1.set_ylabel(sampleAggStat+"_2 - "+sampleAggStat+"_1",fontsize=4)
            ax1.set_xlabel('# bootstrap samples',fontsize=4)
            ax1.set_xscale('log')
            ax1.xaxis.set_tick_params(labelsize=3)
            ax1.yaxis.set_tick_params(labelsize=3)

            minyval, maxyval = pu.get_clean_ax_limits(plotVals=plotVals,symmetric=False)
            ax1.set_ylim(minyval,maxyval)
            #ax1.grid(True)


            ax2 = fig.add_subplot(ny, nx, nx*ifunc+2)
            ax2.hist(clustsN,bins=int(1.2*np.floor(np.sqrt(NOMINAL_NUMBER_OF_GROUPS))))
            ax2.set_ylabel("Count",fontsize=4)
            ax2.set_xlabel('Cluster Size',fontsize=4)
            ax2.grid(True)

        filename = 'cluster_bootstrap_test/comparison' \
                   +'_Nmean='+str(NOMINAL_GROUP_SIZE) \
                   +'_Nstd='+str(SIGMAN) \
                   +'_Nt='+str(N_SUPER_SAMPLE)

        pu.finalize_fig(fig, filename, 'png', True)


def main():
    bootstrap_test()

if __name__ == '__main__': main()
