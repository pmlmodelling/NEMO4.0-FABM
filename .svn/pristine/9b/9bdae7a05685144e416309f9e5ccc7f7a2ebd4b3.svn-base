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
pro std_ts_ICE_Vel, masknp, s_iodir_data,  POSTSCRIPT = postscript, _extra = ex

  compile_opt idl2, strictarrsubs
  
@common
@std_common
; get exp1 info
  vICE1 = getenv('VAR1_ICE')   &   prefix = getenv('V1ICE_PREF')    &   suffix = getenv('V1ICE_SUFF')
; get exp2 info
  vICE2 = getenv('VAR2_ICE')   &   prefix2 = getenv('V2ICE_PREF')   &   suffix2 = getenv('V2ICE_SUFF')
; get exp1 info
  vICE_vel_1 = getenv('VAR1_Ivel')   &   prefix = getenv('V1Iv_PREF')  &   suffix = getenv('V1Iv_SUFF')
; get exp2 info
  vICE_vel_2 = getenv('VAR2_Ivel')   &   prefix2 = getenv('V2Iv_PREF')  &  suffix2 = getenv('V2Iv_SUFF')

  cdti3 = string(cnt, format = '(i3.3)')
  print, cdti3 + ') ' + blabla
  filename = cdti3 + '_ts_ICE_Vel_'+prefix
  if prefix NE prefix2 then filename = filename + '_'+prefix2
  if KEYWORD_SET(postscript) then openps, filename+'.ps', portrait = 1
;
  d1_d2 = '('+strtrim(date1, 1)+' - '+strtrim(date2, 1)+')'
;
  iodir = std_iodir_data
  ; ICE Velocity in NORTH Hemisphere
  domdef, 0, jpi-1, 30, 90, /xindex
  Velo_N = rseries_ncdf(vICE_vel_1, date1, date2, prefix, suffix, FIRSTONLY = 1 - allrec, /nostruct) ;!! warning positive northward
  ICE_N_15 = rseries_ncdf(vICE1, date1, date2, prefix, suffix, FIRSTONLY = 1 - allrec,  /nostruct) 
  print, 'N15', max(ICE_N_15)
  
  
  ICE_N_15[where(ICE_N_15 lt 0.15)] = 0.
  
  ICE_velo_N = grossemoyenne( (Velo_N < 1.e10) * (ICE_N_15 < 1.e10), 'xy',/integration, mask2d = masknp)
  ICE_N_15 = grossemoyenne(ICE_N_15, 'xy',/integration, mask2d = masknp)
  
  ICE_velo_N = ICE_velo_N / ICE_N_15
  
  if jpt mod 12 ne 0 then stop
  nyr=jpt/12.
  ICE_velo_N = reform(ICE_velo_N, 12, nyr)
  ICE_velo_N = total(ICE_velo_N,2)/nyr
  ICE_velo_N = {arr:ICE_velo_N , unit : 'm/s'}  
   
 
  ;ICE Velocity in SOUTH Hemisphere
  domdef, 0, jpi-1, -90, -30, /xindex
  Velo_S = rseries_ncdf(vICE_vel_1, date1, date2, prefix, suffix, FIRSTONLY = 1 - allrec,/nostruct) ;!! warning positive northward

   
    
  title = 'Northern Hemisphere'+'!C'+prefix+' '+d1_d2+'!C'+'Ice Velocity (Black SOLID simulation)'
  jpt=12
  time=julday(1,15,1900)+30*lindgen(12)


  pltt, ICE_velo_N, 't', /REMPLI, /PORTRAIT, XGRIDSTYLE = 1 $
        , small = [1, 2, 1], YTITLE = varunit, TITLE = title, _extra = ex    
  ;
  title ='Southern Hemisphere' +'!C'+prefix+' '+d1_d2+' - '+'!C'+'Ice Velocity (Black SOLID simulation)'
  pltt, Velo_S, 't',  /REMPLI, /NOERASE, XGRIDSTYLE = 1 $
       , small = [1, 2, 2], YTITLE = varunit, TITLE = title, _extra = ex
  ;

  htmltxt = [ htmltxt, '<hr>'+blabla, '<br><img width="80%" src='+filename+'.png  />  ' ]
  if KEYWORD_SET(postscript) then closeps

  if prefix NE prefix2 then BEGIN

    d1_d2_2 = '('+strtrim(date1_2, 1)+' - '+strtrim(date2_2, 1)+')'
    tsave = time
    ; ICE Velocity in NORTH Hemisphere
    domdef, 0, jpi-1, 30, 90, /xindex
    Velo_N_2 = rseries_ncdf(vICE_vel_2, date1_2, date2_2, prefix2, suffix2, FIRSTONLY = 1 - allrec,/nostruct)
    
    ;;;;;; if grossempoyenne is not present dimension of array is not ok: 54   Array[180, 54, 120] 
    ;;;;;; with the domain dimensions [jpi/nx, jpj/ny, jpk/nz, jpt] = [180/180, 149/50, 31/31, 12]
  
    ICE_velo_N_2 = grossemoyenne(Velo_N_2, 'xy', /integration, mask2d = masknp)
    print, 'max ice velo', max(ICE_velo_N_2)
  
    if jpt mod 12 ne 0 then stop
    nyr=jpt/12.
    ICE_velo_N_2 = reform(ICE_velo_N_2, 12, nyr)
    ICE_velo_N_2 = total(ICE_velo_N_2,2)/nyr
    ICE_velo_N_2 = {arr:ICE_velo_N_2 * 86400 * 365, unit : 'm/year'}  
    print, 'max ice velo', max(ICE_velo_N_2.arr)
    
  
    ;ICE Velocity in SOUTH Hemisphere
    domdef, 0, jpi-1, -90, -30, /xindex
    Velo_S_2 = rseries_ncdf(vICE_vel_2, date1_2, date2_2, prefix2, suffix2, FIRSTONLY = 1 - allrec,/nostruct)
      
    if KEYWORD_SET(postscript) then openps, filename+'_2.ps', portrait = 1

    jpt=12
    time=julday(1,15,1900)+30*lindgen(12)

    title = 'Northern Hemisphere'+'!C'+prefix+' (BLACK) - '+prefix2+' (RED) '+d1_d2_2 +'!C'+'Ice Velocity (BLACK) '
    pltt, ICE_velo_N , 't', /REMPLI, XGRIDSTYLE = 1, window = 2  $
          , small = [1, 2, 1], YTITLE = varunit, TITLE = title, /noerase, _extra = ex
    pltt, ICE_velo_N_2, 't', /REMPLI $
          , /ov1d, color = 250, small = [1, 2, 1], YTITLE = varunit, TITLE = title, /noerase, _extra = ex
 
    title = 'Southern Hemisphere'+'!C'+prefix+' (BLACK) - '+prefix2+' (RED) '+d1_d2_2 +'!C'+'Ice Velocity (BLACK) '
    pltt, Velo_S , 't', /REMPLI, XGRIDSTYLE = 1 $
          , small = [1, 2, 2], YTITLE = varunit, TITLE = title, /noerase, _extra = ex
    pltt, Velo_S_2, 't', /REMPLI $
         , /ov1d, color = 250, small = [1, 2, 2], YTITLE = varunit, TITLE = title, /noerase, _extra = ex

    htmltxt = [ htmltxt, '<hr>'+blabla, '<br><img width="80%" src='+filename+'_2.png  />  ' ]
    if KEYWORD_SET(postscript) then closeps

  endif

  domdef
  

  return
end
