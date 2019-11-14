#!/bin/csh

#
#set environment:
# =============================================
source ./setup.csh
source /glade/u/apps/ch/opt/usr/bin/npl/ncar_pylib.csh

echo "CYLC_TASK_CYCLE_POINT=" $CYLC_TASK_CYCLE_POINT

set yymmdd = `echo ${CYLC_TASK_CYCLE_POINT} | cut -c 1-8`
set hh = `echo ${CYLC_TASK_CYCLE_POINT} | cut -c 10-11`
set DATE = ${yymmdd}${hh}

setenv plotdir /glade/scratch/${USER}/pandac/DA_$CYLC_SUITE_NAME/$DATE
mkdir -p $plotdir/graphics
cp  $GRAPHICS_DIR/*.py  $plotdir/graphics/
cd $plotdir/graphics

python plot_cost_grad.py

#python plot_diag_omaomb.py

#python plot_obs_nc_loc.py cycling $DATE 

#python plot_inc.py $DATE method      variable   level_interval                           
python plot_inc.py $DATE 3denvar_bump theta        5       True   $GFSANA_DIR 
#python plot_inc.py $DATE 3denvar_bump theta        5        False
#python plot_inc.py $DATE 3denvar_bump qv           5       True   $GFSANA_DIR 
#python plot_inc.py $DATE 3denvar_bump uReconstructZonal        5  True   $GFSANA_DIR 
#python plot_inc.py $DATE 3denvar_bump uReconstructMeridional   5  True   $GFSANA_DIR 

#
# write diagnostics
# =============================================
cd ..
mkdir diagnostic_stats
cd diagnostic_stats
set DIAGSTATSCRIPT = "write_diagnostic_stats.py"
ln ../graphics/plot_utils.py .
ln ../graphics/write_diagnostic_stats.py .

module load parallel
#set NUMPROC=`cat $PBS_NODEFILE | wc -l`
set NUMPROC=32     #temporary solution. how cylc define $PBS_NODEFILE?
parallel -j${NUMPROC} --plus "python ${DIAGSTATSCRIPT} -n {##} -i {#} >& diags{#}.log" ::: `seq ${NUMPROC}`

#TODO:
if ( ${CYLC_TASK_CYCLE_POINT} == ${CYLC_SUITE_FINAL_CYCLE_POINT} ) then
   echo "cd directory for time serial plot"
   echo "starting time serial plot"
endif

exit
