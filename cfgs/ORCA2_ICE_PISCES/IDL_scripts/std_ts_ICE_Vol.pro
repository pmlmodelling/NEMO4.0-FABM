function read_arr2d, filename, varname, t1, t2
;; function that read input file and return 2d array with monthly timecounter 
nyear = (t2-t1+1)/12
arr2d = ncdf_lec(filename, VAR=varname)
arr2d = arr2d[t1:t2]
arr2d = reform(arr2d,12,nyear) ; put in 2D array
arr2d = total(arr2d,2)/nyear ; total over 2th dimension (i.e.years)

return, arr2d
end 

;
;
pro std_ts_ICE_Vol, masknp, s_iodir_data,  POSTSCRIPT = postscript, _extra = ex

  compile_opt idl2, strictarrsubs
  
@common
@std_common

; get exp1 info
  vICE1     = getenv('VAR1_ICE')     &   prefix = getenv('V1ICE_PREF')  &   suffix = getenv('V1ICE_SUFF')
  v1_Ithick = getenv('VAR1_Ithick')  &   prefix = getenv('V1It_PREF')   &   suffix = getenv('V1It_SUFF')
; get exp2 info
  vICE2     = getenv('VAR2_ICE')     &   prefix2 = getenv('V2ICE_PREF')   &   suffix2 = getenv('V2ICE_SUFF')
  v2_Ithick = getenv('VAR2_Ithick')  &   prefix2 = getenv('V2It_PREF')    &   suffix2 = getenv('V2It_SUFF')

  cdti3 = string(cnt, format = '(i3.3)')
  print, cdti3 + ') ' + blabla
  filename = cdti3 + '_ts_ICE_Vol_'+prefix
  if prefix NE prefix2 then filename = filename + '_'+prefix2
  if KEYWORD_SET(postscript) then openps, filename+'.ps', portrait = 1
;
  d1_d2 = '('+strtrim(date1, 1)+' - '+strtrim(date2, 1)+')'
;
  iodir = std_iodir_data
  ; ICE Area(=Surface) in NORTH Hemisphere
  domdef, 0, jpi-1, 30, 90, /xindex
  ICE_N = rseries_ncdf(vICE1, date1, date2, prefix, suffix, FIRSTONLY = 1 - allrec)
  ICE_thick = rseries_ncdf(v1_Ithick, date1, date2, prefix, suffix, FIRSTONLY = 1 - allrec)
  ; Volume = Area(=Surface) * Thickness
  ICE_vol_N = (ICE_N.arr < 1.e10 ) * ( ICE_thick.arr < 1.e10)  ; limited mask value of 1.e20, because 1.e20 * 1.e20 = inf for idl
  ICE_vol_N = grossemoyenne(ICE_vol_N, 'xy', /integration, mask2d = masknp)

  ;
  if jpt mod 12 ne 0 then stop
  nyr=jpt/12.
  ICE_vol_N = reform(ICE_vol_N, 12, nyr)
  ICE_vol_N = total(ICE_vol_N,2)/nyr
  ICE_vol_N = {arr:ICE_vol_N * 1.e-9, unit : '10^9 Km^3'}
  ;
  ;ICE Area(=Surface) in SOUTH Hemisphere
  domdef, 0, jpi-1, -90, -30, /xindex 
  ICE_S = rseries_ncdf(vICE1, date1, date2, prefix, suffix, FIRSTONLY = 1 - allrec)
  ICE_thick = rseries_ncdf(v1_Ithick, date1, date2, prefix, suffix, FIRSTONLY = 1 - allrec)
  ; Volume = Area(=Surface) * Thickness
  ICE_vol_S = (ICE_S.arr < 1.e10 ) * ( ICE_thick.arr < 1.e10)  ; limited mask value of 1.e20, because 1.e20 * 1.e20 = inf for idl
  ICE_vol_S = grossemoyenne(ICE_vol_S, 'xy', /integration, mask2d = masknp)
  if jpt mod 12 ne 0 then stop
  nyr=jpt/12.
  ICE_vol_S = reform(ICE_vol_S, 12, nyr)
  ICE_vol_S = total(ICE_vol_S,2)/nyr
  ICE_vol_S = {arr:ICE_vol_S * 1.e-9, unit : '10^9 Km^3'}
  ;
  title = 'Northern Hemisphere'+'!C'+prefix+' '+d1_d2+'!C'+'Global Annual Mean Ice Volume (BLACK) '
  jpt=12
  time=julday(1,15,1900)+30*lindgen(12)
  pltt, ICE_vol_N, 't', MIN = 0., MAX = 40000., /REMPLI, /PORTRAIT, XGRIDSTYLE = 1, DATE_FORMAT = '%M' $
        , small = [1, 2, 1], YTITLE = '10^9 Km^3 ', TITLE = title, _extra = ex
;
  title ='Southern Hemisphere' +'!C'+prefix+' '+d1_d2+' - '+'!C'+'Global Annual Mean Ice Volume (BLACK)'
  pltt, ICE_vol_S, 't', MIN = 0., MAX = 12000., /REMPLI, /NOERASE, XGRIDSTYLE = 1 , DATE_FORMAT = '%M' $
        , small = [1, 2, 2], YTITLE = '10^9 Km^3 ', TITLE = title, _extra = ex
;
  htmltxt = [ htmltxt, '<hr>'+blabla, '<br><img width="80%" src='+filename+'.png  />  ' ]
  if KEYWORD_SET(postscript) then closeps

  if prefix NE prefix2 then BEGIN

    d1_d2_2 = '('+strtrim(date1_2, 1)+' - '+strtrim(date2_2, 1)+')'
    tsave = time
    domdef, 0, jpi-1, 30, 90, /xindex
    ;ICE Area(=Surface) in NORTH Hemisphere
    ICE_N2 = rseries_ncdf(vICE2, date1_2, date2_2, prefix2, suffix2, FIRSTONLY = 1 - allrec)
    ICE_thick2 = rseries_ncdf(v2_Ithick, date1_2, date2_2, prefix2, suffix2, FIRSTONLY = 1 - allrec)
    ; Volume = Area(=Surface) * Thickness
    ICE_vol_N2 = (ICE_N2.arr < 1.e10 ) * ( ICE_thick2.arr < 1.e10) ; limited mask value of 1.e20, because 1.e20 * 1.e20 = inf for idl
    ICE_vol_N2 = grossemoyenne(ICE_vol_N2, 'xy', /integration, mask2d = masknp)
    if jpt mod 12 ne 0 then stop
    nyr=jpt/12.
    ICE_vol_N2 = reform(ICE_vol_N2, 12, nyr)
    ICE_vol_N2 = total(ICE_vol_N2,2)/nyr
    ICE_vol_N2 = {arr:ICE_vol_N2 * 1.e-9, unit : '10^3 Km^3'}

    ;ICE Area(=Surface) in SOUTH Hemisphere
    domdef, 0, jpi-1, -90, -30, /xindex 
    ICE_S2 = rseries_ncdf(vICE2, date1_2, date2_2, prefix2, suffix2, FIRSTONLY = 1 - allrec)
    ICE_thick2 = rseries_ncdf(v2_Ithick, date1_2, date2_2, prefix2, suffix2, FIRSTONLY = 1 - allrec)
    ; Volume = Area(=Surface) * Thickness
    ICE_vol_S2 = (ICE_S2.arr < 1.e10 ) * ( ICE_thick.arr < 1.e10) ; limited mask value of 1.e20, because 1.e20 * 1.e20 = inf for idl
    ICE_vol_S2 = grossemoyenne(ICE_vol_S2, 'xy', /integration, mask2d = masknp)
    if jpt mod 12 ne 0 then stop
    nyr=jpt/12.
    ICE_vol_S2 = reform(ICE_vol_S2, 12, nyr)
    ICE_vol_S2 = total(ICE_vol_S2,2)/nyr
    ICE_vol_S2 = {arr:ICE_vol_S2 * 1.e-9, unit : '10^3 Km^3'}

   ; time = tsave   &   IF n_elements(time) NE jpt THEN stop

    if KEYWORD_SET(postscript) then openps, filename+'_2.ps', portrait = 1

    title = 'Northern Hemisphere'+'!C'+prefix+' (BLACK) - '+prefix2+' (RED) '+d1_d2_2 +'!C'+'Global Annual Mean Ice Volume (BLACK) '
    jpt=12
    time=julday(1,15,1900)+30*lindgen(12)
    pltt, ICE_vol_N, 't', MIN = 0., MAX = 40000., /REMPLI, /PORTRAIT, XGRIDSTYLE = 1, window = 2, DATE_FORMAT = '%M' $
          , small = [1, 2, 1], YTITLE = '10^9 Km^3 ', TITLE = title, _extra = ex
    pltt, ICE_vol_N2 , 't', /REMPLI, /PORTRAIT $ ; linee tratteggiate LINESTYLE=2  $
        , /ov1d, COLOR = 250, small = [1, 2, 1], YTITLE = '10^9 Km^3 ', TITLE = title, _extra = ex
    ;
    title = 'Southern Hemisphere'+'!C'+prefix+' (BLACK) - '+prefix2+' (RED) '+d1_d2_2+'!C'+'Global Annual Mean Ice Volume'
    pltt, ICE_vol_S, 't', MIN = 0., MAX = 12000., /REMPLI, /NOERASE, XGRIDSTYLE = 1, DATE_FORMAT = '%M' $
          , small = [1, 2, 2], YTITLE = '10^9 Km^3 ', TITLE = title, _extra = ex   
    pltt, ICE_vol_S2, 't',  /REMPLI, /NOERASE $
          , /ov1d, COLOR = 250, small = [1, 2, 2], YTITLE = '10^9 Km^3 ', TITLE = title, _extra = ex   

    htmltxt = [ htmltxt, '<hr>'+blabla, '<br><img width="80%" src='+filename+'_2.png  />  ' ]
    if KEYWORD_SET(postscript) then closeps

  endif

  domdef
  

  return
end
