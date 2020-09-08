#!/bin/bash
#PBS -N mpasinit.2018041512
#PBS -A NMMM0015
#PBS -l walltime=06:00:00
#PBS -j oe 
#PBS -q share
#PBS -l select=1:ncpus=1
#PBS -V 

#date=2018041500 #2018041500
#lastdate=2018051500
date=2018030100
#lastdate=2018041418
#date=2018051506
lastdate=2018060100


while [ $date -le $lastdate ]
do
  echo "running for "$date
  vyyyy=`echo $date | cut -c1-4`
  vmm=`echo $date | cut -c5-6`
  vdd=`echo $date | cut -c7-8`
  vhh=`echo $date | cut -c9-10`
  mkdir -p ${date}

cat > tmp.ncl << EOF
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
;-----------------------------------------------------------

begin
;-----------------------------------------------------------------------------
; user defined

  FILE_IN  = "./gfs_f48/${date}/x1.40962.init.${vyyyy}-${vmm}-${vdd}_${vhh}.00.00.nc"   ; input filename
  FILE_WGT = "./ESMF_weights/MPAS_x1.40962_to_latlon_1p0_bilinear.nc"  ; weight filename
  FILE_WGT2= "./ESMF_weights/latlon_1p0_to_MPAS_x1.40962_bilinear.nc"  ; weight filename


; end of user defined
;-----------------------------------------------------------------------------
  start_time = systemfunc("date")
  setfileoption("nc","Format","LargeFile")

  f_in = addfile(FILE_IN, "r")

  u_cell = transpose( f_in->uReconstructZonal(0,:,:) )      ; Originally [Time, nCells, nVertLevels] 
  v_cell = transpose( f_in->uReconstructMeridional(0,:,:) ) ; Grap [ nCells, nVerLevels] then transpose into [nVertLevels, nCells]

  dims_cell = dimsizes(u_cell)

  print( " == Interpolate MPAS mesh into Lat/Lon with ESMF " +systemfunc("date"))
  Opt                = True
  Opt@PrintTimings   = True
  ;Opt@Debug          = True
  u_ll = ESMF_regrid_with_weights(u_cell,FILE_WGT, Opt)
  v_ll = ESMF_regrid_with_weights(v_cell,FILE_WGT, Opt)
  print( " == DONE " +systemfunc("date"))

  dims = dimsizes(u_ll)
  nZ = dims(0)
  nY = dims(1)
  nX = dims(2)

  u = new( (/nZ,nY,nX/), double )
  v = new( (/nZ,nY,nX/), double )
  sf = new( (/nZ,nY,nX/), double )
  vp = new( (/nZ,nY,nX/), double )

  u(:,:,:) = u_ll(:,:,:)
  v(:,:,:) = v_ll(:,:,:)

  print (" == Call uv2sfvpf ================== " +systemfunc("date"))
  uv2sfvpf (u, v, sf, vp )
  print (" == Done =========================== " +systemfunc("date"))

;------------------------------------------------
  sf_cell4write = f_in->theta(:,:,:)
  vp_cell4write = f_in->theta(:,:,:)

  sf_cell = ESMF_regrid_with_weights(sf,FILE_WGT2, Opt)
  vp_cell = ESMF_regrid_with_weights(vp,FILE_WGT2, Opt)

  dims_new = dimsizes(sf_cell)

  dims_tmp = dimsizes(sf_cell4write)

  ratio=6371229.0d0/6371220.0d0
  sf_cell_transpose = transpose(sf_cell(:,:) * ratio )
  vp_cell_transpose = transpose( -1.0d0 * vp_cell(:,:) * ratio )
  sf_cell4write(0,:,:)= (/ sf_cell_transpose(:,:) /)
  vp_cell4write(0,:,:)= (/ vp_cell_transpose(:,:) /)

  ;-- attribute
  sf_cell4write@units = "m^2 s^(-2)"
  sf_cell4write@long_name = "stream function"
  vp_cell4write@units = "m^2 s^(-2)"
  vp_cell4write@long_name = "velocity potential"

;------------------------------------------------
 fdir = "./${date}"
 system("cp ./template_PTB.nc "+fdir+"/FULL_f48.nc")
 f_out = addfile(fdir+"/FULL_f48.nc","rw")
 f_out->stream_function    = sf_cell4write
 f_out->velocity_potential = vp_cell4write
 delete(f_out)


end

EOF

ncl tmp.ncl
rm tmp.ncl


  date=`/glade/work/bjung/panda-c/testdata_prepare/da_advance_time.exe $date +6h`
done #--- for date

