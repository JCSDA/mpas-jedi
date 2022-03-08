;group: stream_function-stream_function {
;  dimensions:
;  	nsa = 1083 ;
;  	nsb = 1483 ;
;  	nsc = 4468 ;
;  variables:
;  	double lon_sa(nsa) ;
;  		lon_sa:_FillValue = -3.33476705790481e+38 ;
;  	double lat_sa(nsa) ;
;  		lat_sa:_FillValue = -3.33476705790481e+38 ;
;  	int lev_sa(nsa) ;
;  		lev_sa:_FillValue = -999 ;
;  	double lon_sb(nsb) ;
;  		lon_sb:_FillValue = -3.33476705790481e+38 ;
;  	double lat_sb(nsb) ;
;  		lat_sb:_FillValue = -3.33476705790481e+38 ;
;  	int lev_sb(nsb) ;
;  		lev_sb:_FillValue = -999 ;
;  	double lon_sc(nsc) ;
;  		lon_sc:_FillValue = -3.33476705790481e+38 ;
;  	double lat_sc(nsc) ;
;  		lat_sc:_FillValue = -3.33476705790481e+38 ;
;  	int lev_sc(nsc) ;
;  		lev_sc:_FillValue = -999 ;
;  } // group stream_function-stream_function

begin
  DATADIR="../CMAT_00.tune/"
  full = addfile(DATADIR+"restart.2018-04-15_00.00.00.nc","r")
  RAD2DEG = 180.0d/(atan(1)*4.0d)  ; Radian to Degree
  lonCell = full->lonCell * RAD2DEG
  latCell = full->latCell * RAD2DEG

  f=addfile(DATADIR+"mpas_parametersbump_cov_nicas_grids_local_000128-000049.nc","r")
  GROUP="stream_function-stream_function"
  ;GROUP="temperature-temperature"
  GRP=f=>$GROUP$
  lon_sa=GRP->lon_sa
  lat_sa=GRP->lat_sa
  lev_sa=GRP->lev_sa
  lon_sb=GRP->lon_sb
  lat_sb=GRP->lat_sb
  lev_sb=GRP->lev_sb
  lon_sc=GRP->lon_sc
  lat_sc=GRP->lat_sc
  lev_sc=GRP->lev_sc

  lvl4sort = lev_sc
  lev_scMed = dim_median(lvl4sort)  ; median value
  qsort(lvl4sort)
  idx_lev_closest = closest_val(10,lvl4sort) ; closest level to "10"
  lev4plot=lvl4sort(idx_lev_closest)
  print("lev4plot="+lev4plot)

  lat_sa=where( lev_sa.eq.lev4plot, lat_sa, lat_sa@_FillValue)  ; mask
  lon_sa=where( lev_sa.eq.lev4plot, lon_sa, lon_sa@_FillValue)  ; mask
  lat_sb=where( lev_sb.eq.lev4plot, lat_sb, lat_sb@_FillValue)  ; mask
  lon_sb=where( lev_sb.eq.lev4plot, lon_sb, lon_sb@_FillValue)  ; mask
  lat_sc=where( lev_sc.eq.lev4plot, lat_sc, lat_sc@_FillValue)  ; mask
  lon_sc=where( lev_sc.eq.lev4plot, lon_sc, lon_sc@_FillValue)  ; mask

;==========PLOT
  wks = gsn_open_wks("pdf","plot_nicas_halos_"+GROUP+"_"+sprinti("%0.2i",lev4plot))

  res                             = True
  res@gsnFrame                    = False   ; So we can draw markers
  res@gsnMaximize                 = True    ; affects ps, eps, pdf only
  res@gsnPaperOrientation         = "portrait"  ; force portrait 

  res@tiMainString                = "Full & halo grids, lev = "+sprinti("%i",lev4plot)

  res@mpMinLatF                   = min(lat_sc)-7.    ; range to zoom in on
  res@mpMaxLatF                   = max(lat_sc)+7. 
  res@mpMinLonF                   = min(lon_sc)-7.
  res@mpMaxLonF                   = max(lon_sc)+7.
;  res@mpCenterLonF                = (lonL+lonR)*0.5

  res@mpFillOn = False
  map = gsn_csm_map(wks,res)
;
; Draw markers on the plot in the lat/lon locations.
;
  mkres               = True
  mkres@gsMarkerIndex = 16     ; Filled circle
  mkres@gsMarkerSizeF = 0.005
  mkres@gsMarkerColor = "Green"
  map@locations = gsn_add_polymarker(wks,map,lon_sc,lat_sc,mkres)
  draw(map)
  mkres@gsMarkerIndex = 16     ; Filled circle
  mkres@gsMarkerSizeF = 0.005
  mkres@gsMarkerColor = "Orange"
  map@locations = gsn_add_polymarker(wks,map,lon_sb,lat_sb,mkres)
  draw(map)
  mkres@gsMarkerIndex = 16     ; Filled circle
  mkres@gsMarkerSizeF = 0.005
  mkres@gsMarkerColor = "Red"
  map@locations = gsn_add_polymarker(wks,map,lon_sa,lat_sa,mkres)
  draw(map)
  mkres@gsMarkerIndex = 16     ; Filled circle
  mkres@gsMarkerSizeF = 0.0001
  mkres@gsMarkerColor = "Black"
  map@locations = gsn_add_polymarker(wks,map,lonCell,latCell,mkres)
  draw(map)

  txres                       = True
  txres@txPerimOn             = True
  ;txres@txBackgroundFillColor = "Salmon"
  txres@txFontHeightF         = 0.02
  gsn_text(wks,map,GROUP,-110.,63.,txres)
  draw(map)

  frame(wks)   ; Now advance the frame.

end
