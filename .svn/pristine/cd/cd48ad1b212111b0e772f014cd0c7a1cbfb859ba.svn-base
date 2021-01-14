function read_arr2d, filename, varname, t1, t2
;; function that read input file and return 2d array with monthly timecounter 
nyear = (t2-t1+1)/12
arr2d = ncdf_lec(filename, VAR=varname)
arr2d = arr2d[t1:t2]
arr2d = reform(arr2d,12,nyear) ; put in 2D array
;arr2d = total(arr2d,2)/nyear ; total over 2th dimension (i.e.years)
arr2d = arr2d[*, nyear-1] ;  select last year

return, arr2d
end 

;; here start procedure that use function read_arr2d
pro std_ts_ICE_FRAM, masknp, s_iodir_data,  POSTSCRIPT = postscript, _extra = ex

  compile_opt idl2, strictarrsubs
  
@common
@std_common

; get exp1 info
  vICE1     = getenv('VAR1_ICE')     &   prefix = getenv('V1ICE_PREF')  &   suffix = getenv('V1ICE_SUFF')
  v1_Ithick = getenv('VAR1_Ithick')  &   prefix = getenv('V1It_PREF')   &   suffix = getenv('V1It_SUFF')
  v1_IV     = getenv('VAR1_IvelV')   &   prefix = getenv('V1IvV_PREF')  &   suffix = getenv('V1IvV_SUFF')
; get exp2 info
  vICE2     = getenv('VAR2_ICE')     &   prefix2 = getenv('V2ICE_PREF')   &   suffix2 = getenv('V2ICE_SUFF')
  v2_Ithick = getenv('VAR2_Ithick')  &   prefix2 = getenv('V2It_PREF')    &   suffix2 = getenv('V2It_SUFF')
  v2_IV     = getenv('VAR2_IvelV')   &   prefix2 = getenv('V2IvV_PREF')   &   suffix2 = getenv('V2IvV_SUFF')

  cdti3 = string(cnt, format = '(i3.3)')
  print, cdti3 + ') ' + blabla
  filename = cdti3 + '_ts_ICE_FRAM_'+prefix
  if prefix NE prefix2 then filename = filename + '_'+prefix2
  if KEYWORD_SET(postscript) then openps, filename+'.ps', portrait = 1
;
  d1_d2 = '('+strtrim(date1, 1)+' - '+strtrim(date2, 1)+')'
  d2 = '('+strtrim(date2, 1)+')'
;
  iodir = std_iodir_data
  ; ICE Area(=Surface) in FRAM Strait
  ;; FRAM Strait points, computed with ncview meshmask for ORCA2 grid (133-1434-135-136, 137)
  indx0 = 133 
  indx0_last = 136
  indy0 = 137
  indy0_last = 137
  ;; ORI domdef, 133, 137, 136, 136, /xindex, /yindex,/memeindices
  indx1= (indx0 - ixminmesh + key_shift)mod(jpi)
  indx2= (indx0_last - ixminmesh + key_shift)mod(jpi)
  indy1= indy0 - iyminmesh
  indy2= indy1 
 ;
 ;OBSERVATIONS : mean seasonal cycle/month
 ; vol_obs = [0.261625, 0.230750, 0.325375, 0.252000, 0.172500, 0.0805000, 0.0805000, 0.0805000, 0.0805000, 0.176500, 0.148500, 0.235000]
  vol_obs = [261.625, 230.750, 325.375, 252.000, 172.500, 80.5000, 80.5000, 80.5000, 80.5000, 176.500, 148.500, 235.000]
  area_obs = [0.103292, 0.0997500, 0.107625, 0.0944167, 0.0612083, 0.0262500, 0.0262500, 0.0262500, 0.0262500, 0.0843750, 0.0914583, 0.104083]
 ;
  domdef, indx1, indx2, indy1, indy2, /xindex, /yindex,/memeindices
  ICE = rseries_ncdf(vICE1, date1, date2, prefix, suffix, FIRSTONLY = 1 - allrec, /nostruct)
  ICE_thick = rseries_ncdf(v1_Ithick, date1, date2, prefix, suffix, FIRSTONLY = 1 - allrec, /nostruct)
 ; domdef for V-Point in j-1
 ; REALLY NOT NECESSARY, BECAUSE FLUX CAN BE COMPUTED IN J POINT, is the same
  domdef, indx1, indx2, indy1-1, indy2-1, /xindex, /yindex,/memeindices
  VN = rseries_ncdf(v1_IV, date1, date2, prefix, suffix, FIRSTONLY = 1 - allrec,/nostruct) ;!! warning positive northward

  ;; Area export
  ICE_area_export = (-1) * ICE * VN * ((e1v[firstxv:lastxv, firstyv:lastyv])[*]#replicate(1., jpt))
  ICE_area_export = total(reform(ICE_area_export),1) ; in m2/s -> need to change the unit?

  ;; Volume export
  ICE_vol_export = (-1) * ICE * ICE_thick * VN * ((e1v[firstxv:lastxv, firstyv:lastyv])[*]#replicate(1., jpt))
  ICE_vol_export = total(reform(ICE_vol_export),1) ;! in m3/s -> need to change the unit?

  ;   needed for seasonal cycle :
  if jpt mod 12 ne 0 then stop
  nyr=jpt/12.
  ;; AREA
  ICE_area_export = reform(ICE_area_export, 12, nyr)
  ;ICE_area_export = total(ICE_area_export,2)/nyr  ; old version monthly mean average  over all years
  ICE_area_export = ICE_area_export[*,nyr-1]  ; dim= 12, index 0-11; last year choosen
  ; ICE_area_export = {arr:ICE_area_export * 1.e-12 * 86400 * 365 , unit : '10^6 Km^2/year'}  ; annual mean
  ICE_area_export = {arr:ICE_area_export * 1.e-12 * 86400 * 30 , unit : '10^6 Km^2/month'}    ; monthly mean
  
  ;
  ICE_vol_export = reform(ICE_vol_export, 12, nyr)
  ;ICE_vol_export = total(ICE_vol_export,2)/nyr
  ICE_vol_export = ICE_vol_export[*,nyr-1]
  ; ICE_vol_export = {arr:ICE_vol_export * 1.e-9 * 86400 * 365 , unit : '10^3 Km^3/year'}      ; annual mean
  ICE_vol_export = {arr:ICE_vol_export * 1.e-9 * 86400 * 30 , unit : '10^3 Km^3/month'}      ; monthly mean

  ;
  ;title = 'Fram Strait Areal Export: LAST YEAR'+'!C'+prefix+' '+d1_d2
  title = 'Fram Strait Areal Export: LAST YEAR'+'!C'+prefix+' '+d2
  jpt=12
  time=julday(1,15,1900)+30*lindgen(12)
  pltt, ICE_area_export, 't', /REMPLI, /PORTRAIT, MIN = 0., MAX = .5 , XGRIDSTYLE = 1 $
        , small = [1, 2, 1],YTITLE = '10^6 Km^2/month',  TITLE = title, DATE_FORMAT = '%M', _extra = ex    
  pltt, area_obs, 't', /REMPLI, /NOERASE, psym = 4, THICK = 4  $    ; light blue
         , /ov1d, COLOR = 100, small = [1, 2, 2], YTITLE = '10^6 Km^2/month', TITLE = title, DATE_FORMAT = '%M', _extra = ex
  ;
  tot_area_expo = total(ICE_area_export.arr)
  ;
  xyouts, julday(5,15,1900), 0.49, 'Tot. Annual Export OBS = 0.851 million Km2', ALIGN = 0, CHARTHICK = 2, CHARSIZE=1, COLOR=2
  xyouts, julday(5,15,1900), 0.46, 'Tot. Annual Export Model = '+strtrim(tot_area_expo, 1)+' million Km2', ALIGN = 0, CHARTHICK = 2, CHARSIZE=1, COLOR=2
  xyouts, julday(7,15,1900), 0.41, 'Data from Kwok et al.(2004), 1979-2002 ', ALIGN = 0, CHARTHICK = 2, CHARSIZE=0.8, COLOR=2
  ;
  ;title = 'Fram Strait Volume Export LAST YEAR'+'!C'+prefix+' '+d1_d2
  title = 'Fram Strait Volume Export LAST YEAR'+'!C'+prefix+' '+d2
  pltt, ICE_vol_export, 't',  /REMPLI, MIN = 60., MAX = 500. , /NOERASE, XGRIDSTYLE = 1 $
       , small = [1, 2, 2], YTITLE = '10^3 Km^3/month', TITLE = title, DATE_FORMAT = '%M', _extra = ex
  pltt, vol_obs, 't', /REMPLI, /NOERASE, psym = 4, THICK = 4  $    ; light blue
         , /ov1d, COLOR = 100, small = [1, 2, 2], YTITLE = '10^6 Km^2/month',DATE_FORMAT = '%M', TITLE = title, _extra = ex
  ;
  tot_vol_expo = total(ICE_vol_export.arr)
  tot_vol_expo_Sv = tot_vol_expo * 1.e06 * 1/86400 * 1/365 ; annual mean in Sverdrup
  xyouts, julday(5,15,1900), 490, 'Tot. Annual Export OBS = 2124 10^3 Km3/year', ALIGN = 0, CHARTHICK = 2, CHARSIZE=1, COLOR=2
  xyouts, julday(5,15,1900), 460, 'Tot. Annual Export Model = '+strtrim(tot_vol_expo, 1)+' 10^3 Km3/year', ALIGN = 0, CHARTHICK = 2, CHARSIZE=1, COLOR=2
  xyouts, julday(9,15,1900), 440, 'in Sv = '+strtrim(tot_vol_expo_Sv, 1)+' Sv', ALIGN = 0, CHARTHICK = 2, CHARSIZE=1, COLOR=2
  xyouts, julday(7,15,1900), 360, 'Data from Kwok et al.(2004), 1992-1998 ', ALIGN = 0, CHARTHICK = 2, CHARSIZE=0.8, COLOR=2

  htmltxt = [ htmltxt, '<hr>'+blabla, '<br><img width="80%" src='+filename+'.png  />  ' ]
  if KEYWORD_SET(postscript) then closeps

  if prefix NE prefix2 then BEGIN

    d1_d2_2 = '('+strtrim(date1_2, 1)+' - '+strtrim(date2_2, 1)+')'
    d2_2 = '('+strtrim(date2_2, 1)+')'
    tsave = time
    domdef, indx1, indx2, indy1, indy2, /xindex, /yindex,/memeindices
    ICE_2= rseries_ncdf(vICE2, date1_2, date2_2, prefix2, suffix2, FIRSTONLY = 1 - allrec, /nostruct)
    ICE_thick_2 = rseries_ncdf(v2_Ithick, date1_2, date2_2, prefix2, suffix2, FIRSTONLY = 1 - allrec, /nostruct)
    ; domdef for V-Point in j-1
    ; REALLY NOT NECESSARY, BECAUSE FLUX CAN BE COMPUTED IN J POINT, is the same
    domdef, indx1, indx2, indy1-1, indy2-1, /xindex, /yindex,/memeindices
    VN_2 = rseries_ncdf(v2_IV, date1_2, date2_2, prefix2, suffix2, FIRSTONLY = 1 - allrec,/nostruct) ;!! warning positive northward

    ;; Area export
    ICE_area_export_2 = (-1) * ICE_2 * VN_2 * ((e1v[firstxv:lastxv, firstyv:lastyv])[*]#replicate(1., jpt))
    ICE_area_export_2 = total(reform(ICE_area_export_2),1) ; in m2/s -> need to change the unit?

    ;; Volume export
    ICE_vol_export_2 = (-1) * ICE_2 * ICE_thick_2 * VN_2 * ((e1v[firstxv:lastxv, firstyv:lastyv])[*]#replicate(1., jpt))
    ICE_vol_export_2 = total(reform(ICE_vol_export_2),1) ;! in m3/s -> need to change the unit?

    ;   needed for seasonal cycle :
    if jpt mod 12 ne 0 then stop
    nyr=jpt/12.
    ;; AREA
    ICE_area_export_2 = reform(ICE_area_export_2, 12, nyr)
    ;ICE_area_export_2 = total(ICE_area_export_2,2)/nyr
    ICE_area_export_2 = ICE_area_export_2[*,nyr-1]
    ; ICE_area_export_2 = {arr:ICE_area_export_2 * 1.e-12 * 86400 * 365 , unit : '10^6 Km^2/year'}   ; annual mean
    ICE_area_export_2 = {arr:ICE_area_export_2 * 1.e-12 * 86400 * 30 , unit : '10^6 Km^2/month'}     ; monthly mean
    ;
    ICE_vol_export_2 = reform(ICE_vol_export_2, 12, nyr)
    ;ICE_vol_export_2 = total(ICE_vol_export_2,2)/nyr
    ICE_vol_export_2 = ICE_vol_export_2[*,nyr-1]
    ; ICE_vol_export_2 = {arr:ICE_vol_export_2 * 1.e-12 * 86400 * 365 , unit : '10^3 Km^3/year'}     ; annual mean
    ICE_vol_export_2 = {arr:ICE_vol_export_2 * 1.e-9 * 86400 * 30 , unit : '10^3 Km^3/month'}       ; monthly mean
    ;
    ;
    if KEYWORD_SET(postscript) then openps, filename+'_2.ps', portrait = 1

    ;title = 'Fram Strait Areal Export'+'!C'+prefix+' (BLACK) - '+prefix2+' (RED) '+d1_d2_2
    title = 'Fram Strait Areal Export LAST YEAR'+'!C'+prefix+' (BLACK) - '+prefix2+' (RED) '+d2_2
    jpt=12
    time=julday(1,15,1900)+30*lindgen(12)
    pltt, ICE_area_export, 't', /REMPLI, /PORTRAIT, MIN = 0., MAX = .5, XGRIDSTYLE = 1, window = 2 $
        , small = [1, 2, 1], YTITLE = '10^6 Km^2/month',  TITLE = title, DATE_FORMAT = '%M', _extra = ex    
    pltt, ICE_area_export_2 ,'t', /REMPLI, /PORTRAIT, /NOERASE $
        , /ov1d, COLOR = 250, small = [1, 2, 1],YTITLE = '10^6 Km^2/month',  TITLE = title, DATE_FORMAT = '%M', _extra = ex
    pltt, area_obs, 't', /REMPLI, /NOERASE, psym = 4, THICK = 4  $    ; light blue
         , /ov1d, COLOR = 100, small = [1, 2, 1], YTITLE = '10^6 Km^2/month', TITLE = title, DATE_FORMAT = '%M', _extra = ex
    ;
    tot_area_expo = total(ICE_area_export.arr)
    tot_area_expo_2 = total(ICE_area_export_2.arr)
    ;
    xyouts, julday(5,15,1900), 0.49, 'Tot. Annual Export OBS = 0.851 million Km2', ALIGN = 0, CHARTHICK = 2, CHARSIZE=1, COLOR=2
    xyouts, julday(5,15,1900), 0.46, 'Tot. Annual Export Model 1= '+strtrim(tot_area_expo, 1)+' million Km2', ALIGN = 0, CHARTHICK = 2, CHARSIZE=1, COLOR=2
    xyouts, julday(5,15,1900), 0.43, 'Tot. Annual Export Model 2= '+strtrim(tot_area_expo_2, 1)+' million Km2', ALIGN = 0, CHARTHICK = 2, CHARSIZE=1, COLOR=2
    xyouts, julday(7,15,1900), 0.41, 'Data from Kwok et al.(2004), 1979-2002 ', ALIGN = 0, CHARTHICK = 2, CHARSIZE=0.8, COLOR=2
    ;
    ;title = 'Fram Strait Volume Export'+'!C'+prefix+' (BLACK) - '+prefix2+' (RED) '+d1_d2_2
    title = 'Fram Strait Volume Export LAST YEAR'+'!C'+prefix+' (BLACK) - '+prefix2+' (RED) '+d2_2
    pltt, ICE_vol_export, 't', /REMPLI,  MIN = 60., MAX = 500., /NOERASE, XGRIDSTYLE = 1 $
        , small = [1, 2, 2], YTITLE = '10^3 Km^3/month', TITLE = title, DATE_FORMAT = '%M', _extra = ex
    pltt, ICE_vol_export_2, 't', /REMPLI, /NOERASE  $
        , /ov1d, COLOR = 250, small = [1, 2, 2], YTITLE = '10^3 Km^3/month', TITLE = title, DATE_FORMAT = '%M',  _extra = ex
    pltt, vol_obs, 't', /REMPLI, /NOERASE, psym = 4, THICK = 4  $    ; light blue
         , /ov1d, COLOR = 100, small = [1, 2, 2], YTITLE = '10^6 Km^2/month', TITLE = title, DATE_FORMAT = '%M', _extra = ex
    ;
    tot_vol_expo = total(ICE_vol_export.arr)
    tot_vol_expo_2 = total(ICE_vol_export_2.arr)
    tot_vol_expo_Sv = tot_vol_expo * 1.e06 * 1/86400 * 1/365 ; annual mean in Sverdrup
    tot_vol_expo_2_Sv = tot_vol_expo_2 * 1.e06 * 1/86400 * 1/365 ; annual mean in Sverdrup
    xyouts, julday(5,15,1900), 490, 'Tot. Annual Export OBS = 2124 10^3 Km3/year', ALIGN = 0, CHARTHICK = 2, CHARSIZE=1, COLOR=2
    xyouts, julday(5,15,1900), 460, 'Tot. Annual Export Model 1 = '+strtrim(tot_vol_expo, 1)+' 10^3 Km3/year', ALIGN = 0, CHARTHICK = 2, CHARSIZE=1, COLOR=2
    xyouts, julday(9,15,1900), 440, 'in Sv = '+strtrim(tot_vol_expo_Sv, 1)+' Sv', ALIGN = 0, CHARTHICK = 2, CHARSIZE=1, COLOR=2
    xyouts, julday(5,15,1900), 410, 'Tot. Annual Export Model 2 = '+strtrim(tot_vol_expo_2, 1)+' 10^3 Km3/year', ALIGN = 0, CHARTHICK = 2, CHARSIZE=1, COLOR=2
    xyouts, julday(9,15,1900), 390, 'in Sv = '+strtrim(tot_vol_expo_2_Sv, 1)+' Sv', ALIGN = 0, CHARTHICK = 2, CHARSIZE=1, COLOR=2
    xyouts, julday(7,15,1900), 360, 'Data from Kwok et al.(2004), 1992-1998 ', ALIGN = 0, CHARTHICK = 2, CHARSIZE=0.8, COLOR=2

    htmltxt = [ htmltxt, '<hr>'+blabla, '<br><img width="80%" src='+filename+'_2.png  />  ' ]
    if KEYWORD_SET(postscript) then closeps

  endif

  domdef
  

  return
end
