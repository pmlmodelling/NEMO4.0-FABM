pro std_plot_MLD_ortho, MLD1, MLD2, ARC = arc, ANT = ant, POSTSCRIPT = postscript, _extra = ex

  compile_opt idl2, strictarrsubs

@common  
@std_common

  cdti3 = string(cnt, format = '(i3.3)')
  print, cdti3 + ') ' + blabla

  var = 'MLD'
  IF keyword_set(arc) THEN var = var+'_Arc'
  IF keyword_set(ant) THEN var = var+'_Ant'

  filename = cdti3 + '_'+var+'_'+std_file1_T
  if std_file1_T NE std_file2_T then filename = filename + '_'+std_file2_T
  if KEYWORD_SET(postscript) then openps, filename+'.ps', portrait = 1
  
  IF keyword_set(arc) THEN BEGIN 
    domdef, 20, 380, 40, 90
    map = [90, 0, 0]
  ENDIF
  IF keyword_set(ant) THEN BEGIN 
    domdef, 20, 380, -90, -40
    map = [-90, 0, 0]
  ENDIF
;
;  varunit = .unit
;
  title = var+'!C'+std_file1_T
;                                ;
  plt, MLD1.arr, label=4, cb_label=[0, 15, 25, 50, 75, 100, 125, 150, 250, 400, 600, 850], format = '(I3)' $
       , small = [1, 2, 1], COAST_THICK = 2, TITLE = title $
       , /ORTHO, MAP = map, /PORTRAIT, _extra = ex
  if std_file1_T NE std_file2_T then begin            
    title = title + ' - '+ std_file2_T
    plt, MLD1.arr - MLD2.arr,  MIN = -80., MAX = 80., INTER = 10., format = '(I3)' $
       , small = [1, 2, 2], COAST_THICK = 2, TITLE = title $
       , /ORTHO, MAP = map, /NOERASE, _extra = ex
  endif

  domdef

  htmltxt = [ htmltxt, '<hr>'+blabla, '<br><img width="80%" src='+filename+'.png  />  ' ]
  if KEYWORD_SET(postscript) then closeps

  return
end
