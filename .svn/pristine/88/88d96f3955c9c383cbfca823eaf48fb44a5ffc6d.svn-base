pro std_plot_max_mld, MLD1, MLD2, MARCH = march, SEPT = sept, POSTSCRIPT = postscript, _extra = ex

  compile_opt idl2, strictarrsubs

@common  
@std_common

;
  cdti3 = string(cnt, format = '(i3.3)')
  print, cdti3 + ') ' + blabla
  
  var = 'max_MLD'
  IF keyword_set(march) THEN var = var+'_March'
  IF keyword_set(sept) THEN var = var+'_Sept'
  
  filename = cdti3 + '_'+var+'_'+std_file1_T
  IF std_file1_T NE std_file2_T then filename = filename + '_'+std_file2_T
  IF KEYWORD_SET(postscript) then openps, filename+'.ps', portrait = 1
;
  varunit = MLD1.unit
  titleorg =  var+'!C'
;
;;   IF MLD2.arr[0] NE -1  THEN BEGIN
;;     title = titleorg+std_file1_T+ ' - '+std_file2_T
;;     plt, MLD1.arr - MLD2.arr, MIN = -80., MAX = 80., INTER = 10., FORMAT = '(I3)'  $
;;          , small = [1, 2, 1], COAST_THICK = 2, TITLE = title, /NOCONTOUR, /PORTRAIT, _extra = ex
;;   ENDIF ELSE BEGIN 
;;     title = titleorg+std_file1_T
;;     plt, MLD1,label=4, cb_label=[0, 15, 25, 50, 75, 100, 125, 150, 250, 400, 600, 850],/NOCONTOUR,FORMAT = '(I3)' $
;;        , small = [1, 2, 1],COAST_THICK = 2, TITLE = title, /PORTRAIT, _extra = ex
;;  ;
;;  ;ORI SF    plt, MLD1, MIN = 0., MAX = 500., INTER = 25., /NOCONTOUR, FORMAT = '(I3)' $
;;  ;ORI SF         , small = [1, 2, 1], COAST_THICK = 2, TITLE = title, /PORTRAIT, _extra = ex
;;  ;
;;    ENDELSE 
;

   IF MLD2.arr[0] NE -1  THEN BEGIN
     title = titleorg+ std_file2_T
     ; plt, MLD2.arr, MIN = 0., MAX = 800., INTER = 50., FORMAT = '(I3)'  $
     plt, MLD2.arr, label=4, cb_label=[0, 15, 25, 50, 75, 100, 125, 150, 250, 400, 600, 850],FORMAT = '(I3)' $
          , small = [1, 2, 1], COAST_THICK = 2, TITLE = title, /NOCONTOUR, /PORTRAIT, _extra = ex   
     ;
     title = titleorg+std_file1_T+ ' - '+std_file2_T
     plt, MLD1.arr - MLD2.arr, MIN = -80., MAX = 80., INTER = 10., FORMAT = '(I3)'  $
          , small = [1, 2, 2], COAST_THICK = 2, TITLE = title, /NOCONTOUR, /PORTRAIT, _extra = ex
   ENDIF ELSE BEGIN 
     title = titleorg+std_file1_T
     plt, MLD1,label=4, cb_label=[0, 15, 25, 50, 75, 100, 125, 150, 250, 400, 600, 850],/NOCONTOUR,FORMAT = '(I3)' $
        , small = [1, 2, 1],COAST_THICK = 2, TITLE = title, /PORTRAIT, _extra = ex
   ENDELSE


  
;;  IF MLD2.arr[0] NE -1  THEN BEGIN
;;    title = titleorg+std_file2_T+ ' - DeBoyer'
;;    tmp = MLD2.arr - MLD.arr
;;  ENDIF ELSE BEGIN 
;;    title = titleorg+std_file1_T+ ' - DeBoyer'
;;    tmp = MLD1.arr - MLD.arr
;;  ENDELSE
;;  plt, temporary(tmp), MIN = -80., MAX = 80., INTER = 10., FORMAT = '(I3)'  $
;;       , small = [1, 2, 2], COAST_THICK = 2, TITLE = title, /NOCONTOUR, /NOERASE, _extra = ex
;  
  htmltxt = [ htmltxt, '<hr>'+blabla, '<br><img width="80%" src='+filename+'.png  />  ' ]
  if KEYWORD_SET(postscript) then closeps
;
  return
end

