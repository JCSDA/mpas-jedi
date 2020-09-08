#!/bin/bash

#date=2018041500
#lastdate=2018051500 #2018042018
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

cat > tmpPTB_alt_step1.ncl << EOF
begin

EnsSamp1=new( (/1,40962,55,1/), "double" ) ;stream_function   [ time, cell, level, member ]
EnsSamp2=new( (/1,40962,55,1/), "double" ) ;velocity_potential
EnsSamp3=new( (/1,40962,55,1/), "double" ) ;temperature
EnsSamp4=new( (/1,40962,55,1/), "double" ) ;spechum
EnsSamp5=new( (/1,40962,1/),    "double" )    ;surface_pressure

;======== Read Ensemble
do k=1,1
 FILE_IN  = "./gfs_f48/${date}/x1.40962.init.${vyyyy}-${vmm}-${vdd}_${vhh}.00.00.nc"   ; input filename
 f = addfile(FILE_IN, "r")
 fdir = "./${date}"
 f_inout = addfile(fdir+"/FULL_f48.nc","rw")

;T = TH * ( P / P0 ) ^ (2./7.)
;sh = mixing_ratio / ( 1.0 + mixing_ratio )

 ;EnsSamp1(:,:,:,k-1)=f_inout->stream_function(:,:,:)
 ;EnsSamp2(:,:,:,k-1)=f_inout->velocity_potential(:,:,:)
 EnsSamp3(:,:,:,k-1)=f->theta(:,:,:) * ( (f->pressure_p(:,:,:)+f->pressure_base(:,:,:))/100000.0d0 ) ^ (2.0d0/7.0d0)  ;temperature(:,:,:)
 EnsSamp4(:,:,:,k-1)=f->qv(:,:,:) / (1.0d0 + f->qv(:,:,:))  ;spechum(:,:,:)
 EnsSamp5(:,:,k-1)  =f->surface_pressure(:,:)
end do

;======== Calc Ensemble PTB
do k=1,1

 ;f_inout->stream_function   =(/ EnsSamp1(:,:,:,k-1) /)
 ;f_inout->velocity_potential=(/ EnsSamp2(:,:,:,k-1) /)
 f_inout->temperature       =(/ EnsSamp3(:,:,:,k-1) /)
 f_inout->spechum           =(/ EnsSamp4(:,:,:,k-1) /)
 f_inout->surface_pressure  =(/ EnsSamp5(:,:,k-1)   /)
end do

delete(f)
delete(f_inout)

;======== Save as file
end

EOF

ncl tmpPTB_alt_step1.ncl
rm tmpPTB_alt_step1.ncl

  date=`/glade/work/bjung/panda-c/testdata_prepare/da_advance_time.exe $date +6h`
done #--- for date

