function read_arr2d, filename, varname, t1, t2
;; function that read input file and return 2d array with monthly timecounter 
nyear = (t2-t1+1)/12
arr2d = ncdf_lec(filename, VAR=varname)
arr2d = arr2d[t1:t2]
arr2d = reform(arr2d,12,nyear) ; put in 2D array
arr2d = total(arr2d,2)/nyear ; total over 2th dimension (i.e.years)

return, arr2d
end 

;; here start procedure that use function read_arr2d  
pro std_ts_ICE, masknp, s_iodir_data,  POSTSCRIPT = postscript, _extra = ex

  compile_opt idl2, strictarrsubs
  
@common
@std_common

; get exp1 info
  vICE1 = getenv('VAR1_ICE')   &   prefix = getenv('V1ICE_PREF')    &   suffix = getenv('V1ICE_SUFF')
; get exp2 info
  vICE2 = getenv('VAR2_ICE')   &   prefix2 = getenv('V2ICE_PREF')   &   suffix2 = getenv('V2ICE_SUFF')
; get ice climatology info
  std_file_ice =  isafile(getenv('FILE_ICE'), title = 'ICE Extent Climatology', iodir = std_iodir_climato)
;
  time_ice = ncdf_lec( std_file_ice, VAR='time' )
  time_ice = (time_ice - FLOOR(time_ice) ) * 12
  time_ice = (round(time_ice) + 11) mod 12; round to nearest integer
  t1 = where(time_ice eq 0)
  t1 = t1[0] ;  january
  t2 = where(time_ice eq 11, count)
  t2 = t2[count-1] ; last day of december

  vICE_ext_NH = read_arr2d(std_file_ice, getenv('VAR_ICE_EXT_NH'), t1, t2 )
  vICE_ext_SH = read_arr2d(std_file_ice, getenv('VAR_ICE_EXT_SH'), t1, t2 )
;
  vICE_area_NH = read_arr2d(std_file_ice, getenv('VAR_ICE_area_NH'), t1, t2 )
  vICE_area_SH = read_arr2d(std_file_ice, getenv('VAR_ICE_area_SH'), t1, t2 )
;  
  cdti3 = string(cnt, format = '(i3.3)')
  print, cdti3 + ') ' + blabla
  filename = cdti3 + '_ts_ICE_'+prefix
  if prefix NE prefix2 then filename = filename + '_'+prefix2
  if KEYWORD_SET(postscript) then openps, filename+'.ps', portrait = 1
;
  d1_d2 = '('+strtrim(date1, 1)+' - '+strtrim(date2, 1)+')'
;
  iodir = std_iodir_data
  ; ICE Area in NORTH Hemisphere
  domdef, 0, jpi-1, 30, 90, /xindex
  ICE_N = rseries_ncdf(vICE1, date1, date2, prefix, suffix, FIRSTONLY = 1 - allrec)
  ICE_N = grossemoyenne(ICE_N.arr, 'xy', /integration, mask2d = masknp)

  if jpt mod 12 ne 0 then stop
  nyr=jpt/12.
  ICE_N = reform(ice_n, 12,nyr)
  ICE_n = total(ice_n,2)/nyr
  ICE_N = {arr:ICE_N * 1.e-12, unit : '10^12 m^2'}

  ; ICE EXTENT (Area minus 15%) in NORTH Hemisphere
  ICE_N_15 = rseries_ncdf(vICE1, date1, date2, prefix, suffix, FIRSTONLY = 1 - allrec)
  msk = ICE_N_15.arr gt 0.15 ; remove 0.15% for observations
  ICE_N_15 = grossemoyenne( msk, 'xy', /integration, mask2d = masknp)
  if jpt mod 12 ne 0 then stop
  nyr=jpt/12.
  ICE_N_15 = reform(ice_n_15, 12,nyr)
  ICE_n_15 = total(ice_n_15,2)/nyr
  ICE_N_15 = {arr:ICE_N_15 * 1.e-12, unit : '10^12 m^2'}
  ;
  ;ICE Area in SOUTH Hemisphere
  domdef, 0, jpi-1, -90, -30, /xindex 
  ICE_S = rseries_ncdf(vICE1, date1, date2, prefix, suffix, FIRSTONLY = 1 - allrec)
  ICE_S = grossemoyenne(ICE_S.arr, 'xy', /integration, mask2d = masknp)
  if jpt mod 12 ne 0 then stop
  nyr=jpt/12.
  ICE_S = reform(ice_S, 12,nyr)
  ICE_S = total(ice_S,2)/nyr
  ICE_S = {arr:ICE_S * 1.e-12, unit : '10^12 m^2'}
  ; ICE EXTENT (Area minus 15%) in SOUTH Hemisphere
  ICE_S_15 = rseries_ncdf(vICE1, date1, date2, prefix, suffix, FIRSTONLY = 1 - allrec)
  msk = ICE_S_15.arr gt 0.15 ; remove 0.15% for observations
  ICE_S_15 = grossemoyenne(msk, 'xy', /integration, mask2d = masknp)
  if jpt mod 12 ne 0 then stop
  nyr=jpt/12.
  ICE_S_15 = reform(ice_S_15, 12,nyr)
  ICE_S_15 = total(ice_S_15,2)/nyr
  ICE_S_15 = {arr:ICE_S_15 * 1.e-12, unit : '10^12 m^2'}
  ;
  ;;title = 'Northern Hemisphere'+'!C'+prefix+' (BLACK) '+d1_d2+'!C'+'OBSERVATION (light blue) '+'!C'+' Global Annual Mean Ice Area (CONTINUOUS) '+'!C'+ 'and Extend minus 15% (DASHED)'
  title = 'Northern Hemisphere'+'!C'+prefix+' (BLACK) '+d1_d2+'!C'+'OBSERVATION (light blue) '+'!C'+' Global Annual Mean Ice Area (DASHED) '+'!C'+ 'and Extend minus 15% (CONTINUOUS)'
  jpt=12
  time=julday(1,15,1900)+30*lindgen(12)
  pltt, ICE_N, 't', MIN = 4., MAX = 16., /REMPLI, /PORTRAIT, LINESTYLE=2, XGRIDSTYLE = 1, DATE_FORMAT = '%M' $
       , COLOR = 000 , small = [1, 2, 1], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex
  pltt, ICE_N_15, 't', /REMPLI, /PORTRAIT $       ;;; dashed lines is LINESTYLE=2  $
        , /ov1d, COLOR = 000, small = [1, 2, 1], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex 
  pltt, vICE_area_NH, 't',  /REMPLI, /PORTRAIT, LINESTYLE=2  $ 
        , /ov1d, COLOR = 100, small = [1, 2, 1], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex
  pltt, vICE_ext_NH, 't', /REMPLI, /PORTRAIT $   ;;; dashed lines is LINESTYLE=2  $
        , /ov1d, COLOR = 100, small = [1, 2, 1], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex
;
  title = 'Southern Hemisphere'+'!C'+prefix+' (BLACK) '+d1_d2+'!C'+'OBSERVATION (light blue) '+'!C'+' Global Annual Mean Ice Area (DASHED) '+'!C'+ 'and Extend minus 15% (CONTINUOUS)'
  pltt, ICE_S, 't', MIN = 0., MAX = 20., /REMPLI, LINESTYLE=2, /NOERASE , XGRIDSTYLE = 1 , DATE_FORMAT = '%M' $
        ,COLOR = 000, small = [1, 2, 2], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex
  pltt, ICE_S_15, 't', /REMPLI, /PORTRAIT $ 
        , /ov1d, COLOR = 000, small = [1, 2, 1], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex 
  pltt,  vICE_area_SH, 't', /REMPLI, /PORTRAIT, LINESTYLE=2  $ 
         , /ov1d, COLOR = 100, small = [1, 2, 1], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex
  pltt,  vICE_ext_SH, 't', /REMPLI, /PORTRAIT $ 
        , /ov1d, COLOR = 100, small = [1, 2, 1], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex
;
  htmltxt = [ htmltxt, '<hr>'+blabla, '<br><img width="80%" src='+filename+'.png  />  ' ]
  if KEYWORD_SET(postscript) then closeps

  if prefix NE prefix2 then BEGIN

    d1_d2_2 = '('+strtrim(date1_2, 1)+' - '+strtrim(date2_2, 1)+')'
    tsave = time
    domdef, 0, jpi-1, 30, 90, /xindex
    ;ICE Extent in NORTH Hemisphere
    ICE_N2 = rseries_ncdf(vICE2, date1_2, date2_2, prefix2, suffix2, FIRSTONLY = 1 - allrec)
    ICE_N2 = grossemoyenne(ICE_N2.arr, 'xy', /integration, mask2d = masknp)
    if jpt mod 12 ne 0 then stop
    nyr=jpt/12.
    ICE_N2 = reform(ICE_N2, 12,nyr)
    ICE_N2 = total(ICE_N2,2)/nyr
    ICE_N2 = {arr:ICE_N2 * 1.e-12, unit : '10^12 m^2'}
    ;ICE Extent minus 15% in NORTH Hemisphere
    ICE_N2_15 = rseries_ncdf(vICE2, date1_2, date2_2, prefix2, suffix2, FIRSTONLY = 1 - allrec)
    msk = ICE_N2_15.arr gt 0.15  ; remove 0.15% for observations
    ICE_N2_15 = grossemoyenne( msk, 'xy', /integration, mask2d = masknp)
    if jpt mod 12 ne 0 then stop
    nyr=jpt/12.
    ICE_N2_15 = reform(ICE_N2_15, 12,nyr)
    ICE_N2_15 = total(ICE_N2_15,2)/nyr
    ICE_N2_15 = {arr:ICE_N2_15 * 1.e-12, unit : '10^12 m^2'}
    ;ICE Extent in SOUTH Hemisphere
    domdef, 0, jpi-1, -90, -30, /xindex 
    ICE_S2 = rseries_ncdf(vICE2, date1_2, date2_2, prefix2, suffix2, FIRSTONLY = 1 - allrec)
    ICE_S2 = grossemoyenne(ICE_S2.arr, 'xy', /integration, mask2d = masknp)
    if jpt mod 12 ne 0 then stop
    nyr=jpt/12.
    ICE_S2 = reform(ICE_S2, 12,nyr)
    ICE_S2 = total(ICE_S2,2)/nyr
    ICE_S2 = {arr:ICE_S2 * 1.e-12, unit : '10^12 m^2'}
    ;ICE Extent minus 15% in SOUTH Hemisphere
    ICE_S2_15 = rseries_ncdf(vICE2, date1_2, date2_2, prefix2, suffix2, FIRSTONLY = 1 - allrec)
    msk = ICE_S2_15.arr gt 0.15  ; remove 0.15% for observations
    ICE_S2_15 = grossemoyenne(msk, 'xy', /integration, mask2d = masknp)
    if jpt mod 12 ne 0 then stop
    nyr=jpt/12.
    ICE_S2_15 = reform(ICE_S2_15, 12,nyr)
    ICE_S2_15 = total(ICE_S2_15,2)/nyr
    ICE_S2_15 = {arr:ICE_S2_15 * 1.e-12, unit : '10^12 m^2'}
  ;
 ;   time = tsave   &   IF n_elements(time) NE jpt THEN stop

    if KEYWORD_SET(postscript) then openps, filename+'_2.ps', portrait = 1


  ;;title = 'Northern Hemisphere'+'!C'+prefix+' (BLACK) - '+prefix2+' (RED) '+d1_d2_2+'!C'+'OBSERVATION (light blue) '+'!C'+' Global Annual Mean Ice Area (CONTINUOUS) '+'!C'+ 'and Extend minus 15% (DASHED)'
  title = 'Northern Hemisphere'+'!C'+prefix+' (BLACK) - '+prefix2+' (RED) '+d1_d2_2+'!C'+'OBSERVATION (light blue) '+'!C'+' Global Annual Mean Ice Area (DASHED) '+'!C'+ 'and Extend minus 15% (CONTINUOUS)'
  jpt=12
  time=julday(1,15,1900)+30*lindgen(12)
  pltt, ICE_N, 't', MIN = 4, MAX = 16,  /REMPLI, /PORTRAIT, LINESTYLE=2, XGRIDSTYLE = 1, window = 2, DATE_FORMAT = '%M' $
        , COLOR = 000, small = [1, 2, 1], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex ; BLACK
  pltt, ICE_N2, 't', /REMPLI, /PORTRAIT , LINESTYLE=2 $
        , /ov1d, COLOR = 250, small = [1, 2, 1], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex  ; RED
  pltt, ICE_N_15, 't', /REMPLI, /PORTRAIT $ ; linee tratteggiate LINESTYLE=2  $
        , /ov1d, COLOR = 000, small = [1, 2, 1], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex 
  pltt, ICE_N2_15, 't', /REMPLI, /PORTRAIT $ ; linee tratteggiate LINESTYLE=2  $
        , /ov1d, COLOR = 250, small = [1, 2, 1], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex
  pltt, vICE_area_NH, 't', /REMPLI, /PORTRAIT, LINESTYLE=2  $ 
        , /ov1d, COLOR = 100, small = [1, 2, 1], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex  ; light blue
  pltt, vICE_ext_NH, 't', /REMPLI, /PORTRAIT $ 
        , /ov1d, COLOR = 100, small = [1, 2, 1], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex  ; blu scuro
;
  title ='Southern Hemisphere'+'!C'+prefix+' (BLACK) - '+prefix2+' (RED) '+d1_d2_2+'!C'+'OBSERVATION (light blue) '+'!C'+'Global Annual Mean Ice Area (DASHED)'+'!C'+ 'and Extend minus 15% (CONTINUOUS)'
;  title ='Southern Hemisphere'+'!C'
  pltt, ICE_S, 't', MIN = 0., MAX = 20., /REMPLI, LINESTYLE=2, /NOERASE, XGRIDSTYLE = 1, DATE_FORMAT = '%M' $
         , COLOR = 000, small = [1, 2, 2], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex
  pltt, ICE_S2, 't', /REMPLI, /NOERASE, LINESTYLE=2 $
        , /ov1d, COLOR = 250, small = [1, 2, 2], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex
  pltt, ICE_S_15 , 't', /REMPLI, /PORTRAIT $ ; linee tratteggiate LINESTYLE=2  $
        , /ov1d, COLOR = 000, small = [1, 2, 2], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex
  pltt, ICE_S2_15, 't', /REMPLI, /PORTRAIT $ ; linee tratteggiate LINESTYLE=2  $
        , /ov1d, COLOR = 250, small = [1, 2, 2], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex
  pltt,  vICE_area_SH, 't', /REMPLI, /PORTRAIT, LINESTYLE=2 $ 
          , /ov1d, COLOR = 100, small = [1, 2, 2], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex
  pltt,  vICE_ext_SH, 't', /REMPLI, /PORTRAIT $ 
        , /ov1d, COLOR = 100, small = [1, 2, 2], YTITLE = '10^12 m^2 ', TITLE = title, _extra = ex
;

    htmltxt = [ htmltxt, '<hr>'+blabla, '<br><img width="80%" src='+filename+'_2.png  />  ' ]
    if KEYWORD_SET(postscript) then closeps

  endif

  domdef
  

  return
end
