#!/bin/csh

#
#set environment:
# =============================================
source ./setup.csh
source /glade/u/apps/ch/opt/usr/bin/npl/ncar_pylib.csh

echo "PARAM=" $1
echo "CYLC_TASK_CYCLE_POINT=" ${CYLC_TASK_CYCLE_POINT}
#
# Time info for namelist, yaml etc:
# =============================================
set yymmdd = `echo ${CYLC_TASK_CYCLE_POINT} | cut -c 1-8`
set hh = `echo ${CYLC_TASK_CYCLE_POINT} | cut -c 10-11`

set DADATE = ${yymmdd}${hh}
set PARAM = $1
set fc_num = $1
set ntimefc = `expr $PARAM \* 24 `

set timefc  = `${HOME}/bin/advance_cymdh $DADATE  $ntimefc`

echo "timefc=" ${timefc}
echo "DATE=" ${DADATE}

setenv DATE       $DADATE
setenv VTIME      $timefc

set yy1 = `echo ${VTIME} | cut -c 1-4`
set mm1 = `echo ${VTIME} | cut -c 5-6`
set dd1 = `echo ${VTIME} | cut -c 7-8`
set hh1 = `echo ${VTIME} | cut -c 9-10`

set FILE_DATE  = ${yy1}-${mm1}-${dd1}_${hh1}.00.00
set NAMELIST_DATE  = ${yy1}-${mm1}-${dd1}_${hh1}:00:00
set YAML_DATE  = ${yy1}-${mm1}-${dd1}T${hh1}:00:00Z

set WIND_DATE = `${BIN_DIR}/advance_cymdh ${VTIME} -3`
set yy2 = `echo ${WIND_DATE} | cut -c 1-4`
set mm2 = `echo ${WIND_DATE} | cut -c 5-6`
set dd2 = `echo ${WIND_DATE} | cut -c 7-8`
set hh2 = `echo ${WIND_DATE} | cut -c 9-10`
set WIND_DATE_NAME = ${yy2}-${mm2}-${dd2}T${hh2}:00:00Z

#
# cd working directory: 
# ==================================================
mkdir -p ${OMF_WORK_DIR}/${DATE}/${fc_num}
cd ${OMF_WORK_DIR}/${DATE}/${fc_num}

#
# Copy/link files: namelist, yaml, graph.info, TBL/DBL, background, obs data
# =============================================
cp $DA_NML_DIR/* .
cp v.yaml  orig_v.yaml
cp namelist.atmosphere orig_namelist.atmosphere
ln -fs $GRAPHINFO_DIR/x1.40962.graph.info* .
ln -fs ${CODE_DIR}/libs/build/MPAS_gnu-openmpi_debug=0/src/core_atmosphere/physics/physics_wrf/files/* .
ln -sf $FC2_WORK_DIR/$DATE/restart.$FILE_DATE.nc restart.$FILE_DATE.nc_orig
cp $FC2_WORK_DIR/$DATE/restart.$FILE_DATE.nc .

mkdir Data
ln -fs $OBS_DIR/${VTIME}/* Data/

#
# Revise time info in namelist and yaml 
# =============================================
cat >! newnamelist << EOF
  /config_start_time /c\
   config_start_time      = '${NAMELIST_DATE}'
EOF
sed -f newnamelist orig_namelist.atmosphere >! namelist.atmosphere_${VTIME}
rm newnamelist

sed 's/restart.2018-04-15_00.00.00.nc/restart.'${FILE_DATE}'.nc/g; s/2018041500/'${VTIME}'/g; s/2018-04-15T00:00:00Z/'${YAML_DATE}'/g'  orig_v.yaml  > new1.yaml

cat >! new2.yaml << EOF
  /window_begin: /c\
  window_begin: '${WIND_DATE_NAME}'
EOF
sed -f new2.yaml new1.yaml >! v.yaml
rm new1.yaml new2.yaml
    
cp namelist.atmosphere_${VTIME}  namelist.atmosphere

mpiexec ${CODE_DIR}/build/mpas-bundle/bin/mpas_variational.x  ./v.yaml ./verify.run

#
# check status
# =============================================
grep "Finished running the atmosphere core" log.atmosphere.0000.out
if ( $status != 0 ) then
    touch ./FAIL
    echo "ERROR in $0 : mpas oops failed" >> ./FAIL
    exit 1
endif

#
# write diagnostics
# =============================================
source /glade/u/apps/ch/opt/usr/bin/npl/ncar_pylib.csh " "
mkdir diagnostic_stats
cd diagnostic_stats
set DIAGSTATSCRIPT = "write_diagnostic_stats.py"
cp $GRAPHICS_DIR/plot_utils.py .
cp $GRAPHICS_DIR/write_diagnostic_stats.py .

module load parallel
#set NUMPROC=`cat $PBS_NODEFILE | wc -l`
set NUMPROC=36     #temporary solution. how cylc define $PBS_NODEFILE?
parallel -j${NUMPROC} --plus "python ${DIAGSTATSCRIPT} -n {##} -i {#} >& diags{#}.log" ::: `seq ${NUMPROC}`

exit
