#!/bin/bash
#PBS -N mpasinit.2018041512
#PBS -A NMMM0015
#PBS -l walltime=06:00:00
#PBS -j oe 
#PBS -q share
#PBS -l select=1:ncpus=1
#PBS -V 

#date=2018030100
date=2018040100
lastdate=2018060100


while [ $date -le $lastdate ]
do
  echo "running for "$date
  vyyyy=`echo $date | cut -c1-4`
  vmm=`echo $date | cut -c5-6`
  vdd=`echo $date | cut -c7-8`
  vhh=`echo $date | cut -c9-10`

  cd ${date}
  ncdiff -O FULL_f24.nc FULL_f12.nc PTB_f24mf12.nc
  ncdiff -O FULL_f48.nc FULL_f24.nc PTB_f48mf24.nc
  cd ..

  date=`/glade/work/bjung/panda-c/testdata_prepare/da_advance_time.exe $date +6h`
done #--- for date

