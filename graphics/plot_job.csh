#!/bin/csh
#PBS -N plot_stats_timeseries
#PBS -A NMMM0015
#PBS -q premium
#PBS -l select=6:ncpus=10:mpiprocs=10
#PBS -l walltime=5:00:00
#PBS -m ae
#PBS -k eod
#PBS -o plot.log.job.out 
#PBS -e plot.log.job.err

date

#
# set environment:
# =============================================
module load python/3.7.5
source /glade/u/apps/ch/opt/usr/bin/npl/ncar_pylib.csh

set PLOTSCRIPT="plot_stats_timeseries.py"
set NUMPROC=`cat $PBS_NODEFILE | wc -l`

#
# set environment:
# =============================================

#MULTIPLE PROCESSORS
python ${PLOTSCRIPT} -n ${NUMPROC} >& plot.log

#SINGLE PROCESSOR
#python ${PLOTSCRIPT} >& plot.log

deactivate

date

exit
