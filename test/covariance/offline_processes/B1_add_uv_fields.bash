#!/bin/bash
#PBS -N mpasinit.2018041512
#PBS -A NMMM0015
#PBS -l walltime=06:00:00
#PBS -j oe 
#PBS -q share
#PBS -l select=1:ncpus=1
#PBS -V 

date=2018030100
date=2018041212
lastdate=2018060100

module load nco

while [ $date -le $lastdate ]
do
  echo "running for "$date
  vyyyy=`echo $date | cut -c1-4`
  vmm=`echo $date | cut -c5-6`
  vdd=`echo $date | cut -c7-8`
  vhh=`echo $date | cut -c9-10`

  ncks -A -v uReconstructZonal,uReconstructMeridional gfs_f12/${date}/x1.40962.init.${vyyyy}-${vmm}-${vdd}_${vhh}.00.00.nc ${date}/FULL_f12.nc
  ncks -A -v uReconstructZonal,uReconstructMeridional gfs_f24/${date}/x1.40962.init.${vyyyy}-${vmm}-${vdd}_${vhh}.00.00.nc ${date}/FULL_f24.nc
  ncks -A -v uReconstructZonal,uReconstructMeridional gfs_f48/${date}/x1.40962.init.${vyyyy}-${vmm}-${vdd}_${vhh}.00.00.nc ${date}/FULL_f48.nc


  date=`/glade/work/bjung/panda-c/testdata_prepare/da_advance_time.exe $date +6h`
done #--- for date

