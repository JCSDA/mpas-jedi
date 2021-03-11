#!/bin/csh

echo "CYLC_TASK_CYCLE_POINT=" $CYLC_TASK_CYCLE_POINT

set yymmdd = `echo ${CYLC_TASK_CYCLE_POINT} | cut -c 1-8`
set hh = `echo ${CYLC_TASK_CYCLE_POINT} | cut -c 10-11`
set DATE = ${yymmdd}${hh}

setenv plotdir /glade/scratch/${USER}/pandac/FC2_$CYLC_SUITE_NAME/$DATE
echo $plotdir

echo "fcdiag. Will be updated soon." 

#
# write diagnostics
# =============================================
source /glade/u/apps/ch/opt/usr/bin/npl/ncar_pylib.csh

#TODO: Write diagnostics for model space. Will be updated soon. 

exit
