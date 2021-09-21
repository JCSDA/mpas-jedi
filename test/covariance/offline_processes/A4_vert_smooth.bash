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

  cp ${date}/PTB_f24mf12.nc ${date}/PTB_f24mf12_vertsmooth.nc

cat > tmpPTB_alt_step1.ncl << EOF
begin

EnsSamp1=new( (/1,40962,55,1/), "double" ) ;stream_function   [ time, cell, level, member ]
EnsSamp2=new( (/1,40962,55,1/), "double" ) ;velocity_potential
EnsSamp3=new( (/1,40962,55,1/), "double" ) ;temperature
EnsSamp4=new( (/1,40962,55,1/), "double" ) ;spechum
EnsSamp1s=new( (/1,40962,55,1/), "double" ) ;stream_function   [ time, cell, level, member ]
EnsSamp2s=new( (/1,40962,55,1/), "double" ) ;velocity_potential
EnsSamp3s=new( (/1,40962,55,1/), "double" ) ;temperature
EnsSamp4s=new( (/1,40962,55,1/), "double" ) ;spechum

;======== Read Ensemble
do k=1,1
 fdir = "./${date}"
 f_inout= addfile(fdir+"/PTB_f24mf12_vertsmooth.nc","rw")

 EnsSamp1(:,:,:,k-1)=f_inout->stream_function(:,:,:)
 EnsSamp2(:,:,:,k-1)=f_inout->velocity_potential(:,:,:)
 EnsSamp3(:,:,:,k-1)=f_inout->temperature(:,:,:)
 EnsSamp4(:,:,:,k-1)=f_inout->spechum(:,:,:)
end do

;======== Vertical Smoothing
do k=1,1
  do kk=2, 54 ;for interior levels
    EnsSamp1s(:,:,kk-1,k-1) = 0.5d0 * EnsSamp1(:,:,kk-1,k-1) + 0.25d0 * ( EnsSamp1(:,:,kk-2,k-1) + EnsSamp1(:,:,kk,k-1) )
    EnsSamp2s(:,:,kk-1,k-1) = 0.5d0 * EnsSamp2(:,:,kk-1,k-1) + 0.25d0 * ( EnsSamp2(:,:,kk-2,k-1) + EnsSamp2(:,:,kk,k-1) )
    EnsSamp3s(:,:,kk-1,k-1) = 0.5d0 * EnsSamp3(:,:,kk-1,k-1) + 0.25d0 * ( EnsSamp3(:,:,kk-2,k-1) + EnsSamp3(:,:,kk,k-1) )
    EnsSamp4s(:,:,kk-1,k-1) = 0.5d0 * EnsSamp4(:,:,kk-1,k-1) + 0.25d0 * ( EnsSamp4(:,:,kk-2,k-1) + EnsSamp4(:,:,kk,k-1) )
  end do
  ;for sfc and top
  EnsSamp1s(:,:,0,k-1) = 0.5d0 * ( EnsSamp1(:,:,0,k-1) + EnsSamp1(:,:,1,k-1) )
  EnsSamp2s(:,:,0,k-1) = 0.5d0 * ( EnsSamp2(:,:,0,k-1) + EnsSamp2(:,:,1,k-1) )
  EnsSamp3s(:,:,0,k-1) = 0.5d0 * ( EnsSamp3(:,:,0,k-1) + EnsSamp3(:,:,1,k-1) )
  EnsSamp4s(:,:,0,k-1) = 0.5d0 * ( EnsSamp4(:,:,0,k-1) + EnsSamp4(:,:,1,k-1) )
  EnsSamp1s(:,:,54,k-1) = 0.5d0 * ( EnsSamp1(:,:,54,k-1) + EnsSamp1(:,:,53,k-1) )
  EnsSamp2s(:,:,54,k-1) = 0.5d0 * ( EnsSamp2(:,:,54,k-1) + EnsSamp2(:,:,53,k-1) )
  EnsSamp3s(:,:,54,k-1) = 0.5d0 * ( EnsSamp3(:,:,54,k-1) + EnsSamp3(:,:,53,k-1) )
  EnsSamp4s(:,:,54,k-1) = 0.5d0 * ( EnsSamp4(:,:,54,k-1) + EnsSamp4(:,:,53,k-1) )

end do
;======== Write 
do k=1,1
 f_inout->stream_function   =(/ EnsSamp1s(:,:,:,k-1) /)
 f_inout->velocity_potential=(/ EnsSamp2s(:,:,:,k-1) /)
 f_inout->temperature       =(/ EnsSamp3s(:,:,:,k-1) /)
 f_inout->spechum           =(/ EnsSamp4s(:,:,:,k-1) /)
end do

delete(f_inout)

end

EOF

ncl tmpPTB_alt_step1.ncl
rm tmpPTB_alt_step1.ncl

  date=`/glade/work/bjung/panda-c/testdata_prepare/da_advance_time.exe $date +6h`
done #--- for date

