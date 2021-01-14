pro std_plot_IceVel, IvelU1, IvelU2, IvelV1, IvelV2, Ivelo1, Ivelo2, ARC = arc, ANT = ant, FEBR = febr, MARCH = march, SEPT = sept, POSTSCRIPT = postscript, _extra = ex


@common  
@std_common

  cdti3 = string(cnt, format = '(i3.3)')
  print, cdti3 + ') ' + blabla

  var = 'IceVel'
  IF keyword_set(arc) THEN var = var+'_Arc_'
  IF keyword_set(ant) THEN var = var+'_Ant_'
  IF keyword_set(febr) THEN var = var+'Febr'
  IF keyword_set(march) THEN var = var+'March'
  IF keyword_set(sept) THEN var = var+'Sept'

  filename = cdti3 + '_'+var+'_'+std_file1_I
  if std_file1_I NE std_file2_I then filename = filename + '_'+std_file2_I
  if KEYWORD_SET(postscript) then openps, filename+'.ps', portrait = 1
 
  IF keyword_set(arc) THEN BEGIN
     domdef, 20, 380, 50, 90
     map = [90, 0, 0]
  ENDIF
  IF keyword_set(ant) THEN BEGIN
     domdef, 20, 380, -90, -50
     map = [-90, 0, 0]
  ENDIF
  ;
  title = var+'!C'+std_file1_I
   
  Ivelo1 = {arr: Ivelo1.arr - 1.E-04, g: 'T'}
  ; Ivelo1 = {arr: Ivelo1.arr, g: 'T'}
  
  plt, Ivelo1, vecteur={u:IvelU1, v:IvelV1}, unvectsur=[3,3], normeref = 0.5, cmref = 1.  $
        , small = [1, 2, 1], TITLE = title $
        , /ORTHO, MAP = map, /PORTRAIT, _extra = ex
 ;                                
 
  if std_file1_I NE std_file2_I then begin
     title = title + std_file2_I  
  ;   
     Ivelo = {arr: Ivelo1.arr - Ivelo2.arr, g: 'T'}
     ; Ivelo1 = {arr: Ivelo1.arr, g: 'T'}
     ; Ivelo2 = {arr: Ivelo2.arr - 1.E-04, g: 'T'}
     Ivelo2 = {arr: Ivelo2.arr , g: 'T'}
    
     plt, Ivelo, vecteur={u: IvelU2, v: IvelV2}, unvectsur=[3,3], normeref = 0.5, cmref = 1.  $
          , small = [1, 2, 2], TITLE = title $
          , /ORTHO, MAP = map, /NOERASE, _extra = ex
 endif
 
   domdef
 
   htmltxt = [ htmltxt, '<hr>'+blabla, '<br><img width="80%" src='+filename+'.png  />  ' ]
   if KEYWORD_SET(postscript) then closeps
 
   return
 end
 