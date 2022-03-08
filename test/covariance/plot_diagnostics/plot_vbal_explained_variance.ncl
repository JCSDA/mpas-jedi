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

  RAD2DEG = 180.0d/(atan(1)*4.0d)  ; Radian to Degree
  x_axis = fsamp->lat_c2 * RAD2DEG


  var="stream_function-temperature"
  grp=f=>$var$  
  dum1= grp->explained_var_c2(:,:)   ; (/55,nc2/)

  var="stream_function-velocity_potential"
  grp=f=>$var$  
  dum2= grp->explained_var_c2(:,:)     ; (/55,nc2/)) 
  print( "min/max of reg12 = "+min(dum2)+" / "+max(dum2) )

  var="stream_function-surface_pressure"
  grp=f=>$var$  
  dum3= grp->explained_var_c2(:,:)     ; (/55,nc2/)) 
  print( "min/max of reg15 = "+min(dum3)+" / "+max(dum3) )


;------ plot wks -------
  plot = new(3,graphic)
  wks  = gsn_open_wks("pdf","plot_vbal_explained")
  
  res                      = True
  res@gsnDraw              = False
  res@gsnFrame             = False

  res@cnFillOn             = True         ; turn on color
  res@cnLineLabelsOn       = False        ; no contour line labels
  res@cnInfoLabelOn        = False        ; no info label
  res@lbOrientation        = "vertical"
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = 0.05
  res@cnMaxLevelValF       = 1.0
  res@cnLevelSpacingF      = 0.05

  res@tiMainString         = " Explained T with psi"
  res@tiXAxisString        = "Latitude(nc2, south-to-north)"
  res@tiYAxisString        = "Explained Variance"
  plot(0) = gsn_csm_contour(wks,dum1,res)

  res@tiMainString = " Explained chi with psi"
  plot(1) = gsn_csm_contour(wks,dum2,res)

  res@trXMinF = 1
  res@trXMaxF = 200  ; nc2
  x=ispan(1,200,1)
  res@tiMainString  = " Explained Ps with psi"
  plot(2) = gsn_csm_xy(wks,x,dum3(0,:),res)

;************************************************
; create panel
;************************************************
  pres                    = True              ; mod panel plot
  pres@gsnPaperOrientation = "landscape"
  pres@gsnMaximize         = True        
  gsn_panel(wks,plot,(/1,3/),pres)             ; now draw as one plot

end
