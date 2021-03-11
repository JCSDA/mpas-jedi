#!/bin/csh -f

#
# Set up environment:
# =============================================
source /etc/profile.d/modules.csh
setenv OPT /glade/work/miesch/modules
module use $OPT/modulefiles/core
module purge
module load jedi/gnu-openmpi
#module list
module load nco
module load python
limit stacksize unlimited
#setenv OOPS_TRACE 1
setenv GFORTRAN_CONVERT_UNIT 'native;big_endian:101-200'
setenv F_UFMTENDIAN 'big:101-200'

#
# Settings
# =============================================
setenv CYCLE_PERIOD      6
setenv FC2_FCST_LENGTH   10_00:00:00   # 10day fc
setenv FC2_FCST_INTERVAL 1_00:00:00
setenv UPDATESST         True #False 
setenv TOP_DIR           /gpfs/fs1/scratch/${USER}/pandac
setenv BIN_DIR           /glade/u/home/${USER}/bin  # path for advance_cymdh.f90
setenv CODE_DIR          /glade/work/${USER}/pandac # path for mpas-bundle and libs
setenv GFSANA6HFC_DIR    ${TOP_DIR}/FC_GFSANA6HFC   # restart file path for initial time DA input
setenv GFSANA_DIR        ${TOP_DIR}/120km_GFSANA    # path for gfs analysis and sst files
setenv GRAPHINFO_DIR     ${TOP_DIR}/120km_graph     # path for graph info
setenv DA_NML_DIR        ${TOP_DIR}/120km_DA_NML    # path for yaml, stream list, namelist
setenv FC_NML_DIR        ${TOP_DIR}/120km_FC_NML    # path for stream list, namelist
setenv OBS_DIR           ${TOP_DIR}/1mon_obs        # path for 1-month ioda obs
setenv filesbump         ${TOP_DIR}/120km_filesbump # path for bump files
setenv GRAPHICS_DIR      ${TOP_DIR}/graphics        # path for post-processing scripts 
setenv currdir  `basename "$PWD"`

setenv DA_WORK_DIR       ${TOP_DIR}/DA_$currdir     # path for DA output
setenv FC1_WORK_DIR      ${TOP_DIR}/FC1_$currdir    # path for 6-h FC output
setenv FC2_WORK_DIR      ${TOP_DIR}/FC2_$currdir    # path for 10-DAY FC output
setenv OMF_WORK_DIR      ${TOP_DIR}/OMF_$currdir    # path for obs minus FC

setenv SCRIPT_DIR        /gpfs/fs1/p/mmm/parc/${USER}/test/pandac/work/$currdir

