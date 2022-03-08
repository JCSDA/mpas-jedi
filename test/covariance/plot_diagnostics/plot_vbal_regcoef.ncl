;netcdf mpas_vbal_vbal {
;dimensions:
;	nc2 = 200 ;
;	nl0_1 = 55 ;
;	nl0_2 = 55 ;
;group: stream_function-temperature {
;  variables:
;  	double reg_c2(nl0_2, nl0_1, nc2) ;
;  		reg_c2:_FillValue = -3.33476705790481e+38 ;
;  	double cov_c2(nl0_2, nl0_1, nc2) ;
;  		cov_c2:_FillValue = -3.33476705790481e+38 ;
;  	double explained_var_c2(nl0_1, nc2) ;
;  		explained_var_c2:_FillValue = -3.33476705790481e+38 ;
;  } // group stream_function-temperature


begin
  DATADIR="../VBAL_00/"
  f=addfile(DATADIR+"mpas_vbal_vbal.nc","r")
  fsamp=addfile(DATADIR+"mpas_vbal_sampling.nc","r")


nz=55          ; number of vertical levels
idx_nc2=150    ; nc2 index for plotting

  RAD2DEG = 180.0d/(atan(1)*4.0d)  ; Radian to Degree
  clat=fsamp->lat_c2(idx_nc2)*RAD2DEG
  clon=fsamp->lon_c2(idx_nc2)*RAD2DEG
  print( "=== lat / lon = " + clat + " / " + clon )

  var="stream_function-temperature"
  grp=f=>$var$ 
  dum= grp->reg_c2(:,:,idx_nc2)  ;(/55,55,nc2/)) 
  print( "min/max of reg13 = "+min(dum)+" / "+max(dum) )
  var="stream_function-velocity_potential"
  grp=f=>$var$  
  dum2= grp->reg_c2(:,:,idx_nc2)  ;(/55,55,nc2/)) 
  dum3=dum2(0,:)
  do k=0,55-1
    dum3(k)=(/ dum2(k,k) /)
  end do
  print( "min/max of reg12 = "+min(dum3)+" / "+max(dum3) )
  var="stream_function-surface_pressure"
  grp=f=>$var$  
  dum4= grp->reg_c2(:,:,idx_nc2)   ;(/55,55,nc2/)) 
  dum5=dum2(0,:)
  do k=0,55-1
    dum5(k)=(/ dum4(k,0) /)
  end do
  print( "min/max of reg15 = "+min(dum5)+" / "+max(dum5) )


;------ plot wks -------
  plot = new(3,graphic)
  wks  = gsn_open_wks("pdf","plot_vbal_regcoef_"+idx_nc2)
  
  res                      = True
  res@gsnDraw              = False
  res@gsnFrame             = False

  res@cnFillOn             = True         ; turn on color
  res@cnLineLabelsOn       = False        ; no contour line labels
  res@cnInfoLabelOn        = False        ; no info label
  res@lbOrientation        = "vertical"
  res@tiXAxisString        = "psi"        ; x-axis label
  res@tiYAxisString        = "T"          ; y-axis label
  aa=dum(:,:)
  symMinMaxPlt (aa(0:40,0:40),20,False,res)
  res@tiMainString = "psi-T reg.coef. @lat/lon = "+sprintf("%6.2f",clat)+" / "+sprintf("%6.2f",clon)
  plot(0)=gsn_csm_contour(wks,aa,res)           ; plot zonal ave

  res@trYMinF = 1
  res@trYMaxF = 55

  x=ispan(1,55,1)
  res@tiYAxisString     = "Vertical levels"
  res@tiMainString      = "chi-psi"
  plot(1) = gsn_csm_xy(wks,dum3,x,res)

  res@tiMainString      = "ps-psi [* 1e-4]"
  plot(2) = gsn_csm_xy(wks,1e4*dum5,x,res)


;************************************************
; create panel
;************************************************
  pres                    = True              ; mod panel plot
  pres@gsnPaperOrientation = "landscape"
  pres@gsnMaximize         = True        
  gsn_panel(wks,plot,(/1,3/),pres)             ; now draw as one plot

end
