#!/bin/csh

#
# Set up environment:
# =============================================
source ./setup.csh

echo "Run FC: CYLC_TASK_CYCLE_POINT=" $CYLC_TASK_CYCLE_POINT

#
# Time info for namelist, yaml etc:
# =============================================
set yymmdd = `echo ${CYLC_TASK_CYCLE_POINT} | cut -c 1-8`
set hh = `echo ${CYLC_TASK_CYCLE_POINT} | cut -c 10-11`

set DATE = ${yymmdd}${hh}

set yy1 = `echo ${DATE} | cut -c 1-4`
set mm1 = `echo ${DATE} | cut -c 5-6`
set dd1 = `echo ${DATE} | cut -c 7-8`
set hh1 = `echo ${DATE} | cut -c 9-10`
set FILE_DATE  = ${yy1}-${mm1}-${dd1}_${hh1}.00.00
set NAMELIST_DATE  = ${yy1}-${mm1}-${dd1}_${hh1}:00:00

#
# cd working directory:
# =============================================
#delete old results
rm -fr  $FC1_WORK_DIR/$DATE
#create new dir
mkdir -p $FC1_WORK_DIR/$DATE
cd $FC1_WORK_DIR/$DATE

#
# Copy/link files: 
# =============================================
cp $FC_NML_DIR/* .
cp namelist.atmosphere orig_namelist.atmosphere
ln -fs $GRAPHINFO_DIR/x1.40962.graph.info* .
ln -fs ${CODE_DIR}/libs/build/MPAS_gnu-openmpi_debug=0/src/core_atmosphere/physics/physics_wrf/files/* .
#link background:
ln -fs $DA_WORK_DIR/$DATE/restart.${FILE_DATE}.nc  .
#link analysis:
ln -fs $DA_WORK_DIR/$DATE/mpas.3denvar_bump.${FILE_DATE}.nc .

#
# Update analyzed variables:
# =============================================
ncks -A -v theta,rho,u,qv,uReconstructZonal,uReconstructMeridional mpas.3denvar_bump.${FILE_DATE}.nc restart.${FILE_DATE}.nc

#
# Revise time info in namelist
# =============================================
cat >! newnamelist << EOF
  /config_start_time /c\
   config_start_time      = '${NAMELIST_DATE}'
EOF
sed -f newnamelist orig_namelist.atmosphere >! namelist.atmosphere
rm newnamelist
#
# Run the executable:
# =============================================
mpiexec ${CODE_DIR}/libs/build/MPAS_gnu-openmpi_debug=0/atmosphere_model

#
# write diagnostics
# =============================================
source /glade/u/apps/ch/opt/usr/bin/npl/ncar_pylib.csh

#TODO: Write diagnostics for model space. Will be updated soon.

exit
