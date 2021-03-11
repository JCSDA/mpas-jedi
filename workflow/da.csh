#!/bin/csh

#
#set environment:
# =============================================
source ./setup.csh

echo "CYLC_TASK_CYCLE_POINT=" $CYLC_TASK_CYCLE_POINT
#
# Time info for namelist, yaml etc:
# =============================================
set yymmdd = `echo ${CYLC_TASK_CYCLE_POINT} | cut -c 1-8`
set hh = `echo ${CYLC_TASK_CYCLE_POINT} | cut -c 10-11`

set DATE = ${yymmdd}${hh}
set PREV_DATE = `${BIN_DIR}/advance_cymdh ${DATE} -${CYCLE_PERIOD}`

set yy1 = `echo ${DATE} | cut -c 1-4`
set mm1 = `echo ${DATE} | cut -c 5-6`
set dd1 = `echo ${DATE} | cut -c 7-8`
set hh1 = `echo ${DATE} | cut -c 9-10`

set FILE_DATE  = ${yy1}-${mm1}-${dd1}_${hh1}.00.00
set NAMELIST_DATE  = ${yy1}-${mm1}-${dd1}_${hh1}:00:00
set YAML_DATE  = ${yy1}-${mm1}-${dd1}T${hh1}:00:00Z

set WIND_DATE = `${BIN_DIR}/advance_cymdh ${DATE} -3`
set yy2 = `echo ${WIND_DATE} | cut -c 1-4`
set mm2 = `echo ${WIND_DATE} | cut -c 5-6`
set dd2 = `echo ${WIND_DATE} | cut -c 7-8`
set hh2 = `echo ${WIND_DATE} | cut -c 9-10`
set WIND_DATE_NAME = ${yy2}-${mm2}-${dd2}T${hh2}:00:00Z

#
# Link restart files for inital time
# =============================================
if ($CYLC_TASK_CYCLE_POINT == $CYLC_SUITE_INITIAL_CYCLE_POINT) then
   mkdir -p ${FC1_WORK_DIR}/
   cd ${FC1_WORK_DIR}/
   ln -fs  ${GFSANA6HFC_DIR}/$PREV_DATE   .
endif

#
# cd working directory:
# =============================================
mkdir -p ${DA_WORK_DIR}/${DATE}
cd ${DA_WORK_DIR}/${DATE}

#
# Copy/link files: namelist, yaml, graph.info, TBL/DBL, background, obs data
# =============================================
cp $DA_NML_DIR/* .
cp 3denvar_bumploc.yaml  orig_3denvar_bumploc.yaml
cp namelist.atmosphere orig_namelist.atmosphere
ln -fs $GRAPHINFO_DIR/x1.40962.graph.info* .
ln -fs ${CODE_DIR}/libs/build/MPAS_gnu-openmpi_debug=0/src/core_atmosphere/physics/physics_wrf/files/* .
ln -sf $FC1_WORK_DIR/$PREV_DATE/restart.$FILE_DATE.nc restart.$FILE_DATE.nc_orig
cp $FC1_WORK_DIR/$PREV_DATE/restart.$FILE_DATE.nc .

mkdir Data
ln -fs $OBS_DIR/${DATE}/* Data/

#
# Revise time info in namelist and yaml 
# =============================================
cat >! newnamelist << EOF
  /config_start_time /c\
   config_start_time      = '${NAMELIST_DATE}'
EOF
sed -f newnamelist orig_namelist.atmosphere >! namelist.atmosphere_${DATE}
rm newnamelist

sed 's/restart.2018-04-15_00.00.00.nc/restart.'${FILE_DATE}'.nc/g; s/2018041500/'${DATE}'/g; s/2018-04-15T00:00:00Z/'${YAML_DATE}'/g'  orig_3denvar_bumploc.yaml  > new0.yaml
sed 's/x1.40962.init.2018-04-15_00.00.00.nc/x1.40962.init.'${FILE_DATE}'.nc/g' new0.yaml > new1.yaml

cat >! new2.yaml << EOF
  /window_begin: /c\
  window_begin: '${WIND_DATE_NAME}'
  /datadir: /c\
          datadir: ${filesbump}
EOF

sed -f new2.yaml new1.yaml >! 3denvar_bumploc.yaml
rm new0.yaml new1.yaml new2.yaml

cp namelist.atmosphere_${DATE}  namelist.atmosphere

#
# Update sst,xice:
# =============================================
if ($UPDATESST == true) then
   #delete, then, append
   #delete sst 
   setenv SST_FILE ${GFSANA_DIR}/${DATE}/x1.40962.sfc_update.${FILE_DATE}.nc
   ncks -a -x -v sst,xice restart.${FILE_DATE}.nc restart.${FILE_DATE}_nosstice.nc

   #append sst,xice 
   ncks -A -v sst,xice ${SST_FILE} restart.${FILE_DATE}_nosstice.nc
   mv  restart.${FILE_DATE}_nosstice.nc  restart.${FILE_DATE}.nc
endif

#
# Run the executable:
# =============================================
mpiexec ${CODE_DIR}/build/mpas-bundle/bin/mpas_variational.x  ./3denvar_bumploc.yaml ./3denvar_bumploc.run
#
# Check status:
# =============================================
grep "Finished running the atmosphere core" log.atmosphere.0000.out
if ( $status != 0 ) then
    touch ./FAIL
    echo "ERROR in $0 : mpas oops failed" >> ./FAIL
    exit 1
endif

#
# Link log file to testoutput dir:
# =============================================
mkdir testoutput
cd testoutput
ln -fs ../3denvar_bumploc.run .

exit
