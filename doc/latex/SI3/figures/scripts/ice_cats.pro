
;
;==============================================================================
; 1) User parameters
;==============================================================================
;

out_dir=''
file_out   = 'ThicknessSpace'

numplot_x = 1       ; number of horizontal plots
numplot_y = 1       ; number of vertical plots
ct = 13             ; colortable
cs = 2.0            ; size of fonts on plots (1=normal)
th = 2.             ; thickness of curves
device = 'PS'       ; 'PS' or 'X'

;
;==============================================================================
; 2) Initialize graphics
;==============================================================================
;

dx = 12.
dy = 8.
figuresize_x = numplot_x * dx; figure size on x direction
figuresize_y = numplot_y * dy; figure size on y direction
SET_PLOT, DEVICE
;!P.FONT=1
TVLCT, [0,255,0,0], [0,0,255,0], [0,0,0,255]
IF ( device EQ 'PS' ) THEN BEGIN
;   DEVICE, /COLOR, /LANDSCAPE, filename=out_dir+file_out+'.ps', $
;           XSIZE=figuresize_x,YSIZE=figuresize_y,FONT_SIZE=12.0, SET_FONT='Helvetica', BITS = 16
   DEVICE, /COLOR, /LANDSCAPE, filename=out_dir+file_out+'.ps', $
           XSIZE=figuresize_x,YSIZE=figuresize_y,FONT_SIZE=9.0
ENDIF

IF ( device EQ 'X' ) THEN BEGIN
   xsize = 1200
   ysize = 800
   colorkey = 'rd'
   init_graphics_x, xsize, ysize, colorkey
ENDIF

!P.MULTI=[0,numplot_x,numplot_y]
!X.MARGIN   = [7,3]
!Y.MARGIN   = [4,3]

;
;==============================================================================
; 2) Initialize graphics
;==============================================================================
;

jpl = 5
jpl1 = 10

;-----------------------
; original CICE version
;-----------------------

hi_max_cice_5 = FLTARR(jpl+1)
hi_max_cice_5(0) = 0
zc1 =  3. / FLOAT(jpl)
zc2 = 10. * zc1 
zc3 =  3.

FOR jl = 1, jpl  DO BEGIN
   zx1 = FLOAT( jl-1 ) / FLOAT( jpl )
   hi_max_cice_5(jl) = hi_max_cice_5(jl-1) + zc1 + zc2 * (1. + TANH( zc3 * (zx1 - 1. ) ) ) 
ENDFOR

hi_max_cice_10 = FLTARR(jpl1+1)
hi_max_cice_10(0) = 0
zc1 =  3. / FLOAT(jpl1)
zc2 = 10. * zc1 
zc3 =  3.

FOR jl = 1, jpl1  DO BEGIN
   zx1 = FLOAT( jl-1 ) / FLOAT( jpl1 )
   hi_max_cice_10(jl) = hi_max_cice_10(jl-1) + zc1 + zc2 * (1. + TANH( zc3 * (zx1 - 1. ) ) ) 
ENDFOR

;-------------------
; LIM3 distribution
;-------------------
alpha = 1./20.

; 5 categories, Arctic ice
hmean = 2.5
hL = 3.*hmean

hi_max_1 = FLTARR(jpl + 1)
hi_max_1(0) = 0.

FOR jl = 1, jpl DO BEGIN
   znum = jpl * ( hL+1 )^alpha
   zden = ( jpl - jl ) * ( hL+1 )^alpha + jl
   hi_max_1(jl) = ( znum / zden )^(1./alpha) - 1
ENDFOR

; 10 categories, Arctic ice
hi_max_2 = FLTARR(jpl1 + 1)
hi_max_2(0) = 0.
FOR jl = 1, jpl1 DO BEGIN
   znum = jpl1 * ( hL+1 )^alpha
   zden = ( jpl1 - jl ) * ( hL+1 )^alpha + jl
   hi_max_2(jl) = ( znum / zden )^(1./alpha) - 1
ENDFOR

; 5 categories, Spitzberg
hmean = 1.0
hL = 3.*hmean

hi_max_3 = FLTARR(jpl + 1)
hi_max_3(0) = 0.

FOR jl = 1, jpl DO BEGIN
   znum = jpl * ( hL+1 )^alpha
   zden = ( jpl - jl ) * ( hL+1 )^alpha + jl
   hi_max_3(jl) = ( znum / zden )^(1./alpha) - 1
ENDFOR

; 10 categories, Spitzberg
hmean = 1.0
hi_max_4 = FLTARR(jpl1 + 1)
hi_max_4(0) = 0.

FOR jl = 1, jpl1 DO BEGIN
   znum = jpl1 * ( hL+1 )^alpha
   zden = ( jpl1 - jl ) * ( hL+1 )^alpha + jl
   hi_max_4(jl) = ( znum / zden )^(1./alpha) - 1
ENDFOR

!P.MULTI = [ 0, numplot_x, numplot_y ]

PLOT, [0.1, 10.], [0., 7.], /NODATA, charsize = cs, XSTYLE = 1, YSTYLE = 1
OPLOT, hi_max_cice_5, FINDGEN(jpl+1)*0.+1, psym = 4, symsize = 1, thick = 2
OPLOT, hi_max_cice_10, FINDGEN(jpl1+1)*0.+2, psym = 4, symsize = 1, thick = 2, color =150
OPLOT, hi_max_1, FINDGEN(jpl+1)*0.+3   , psym = 1, symsize = 1, thick = 2, color = 0
OPLOT, hi_max_2, FINDGEN(jpl1+1)*0.+4  , psym = 1, symsize = 1, thick = 2, color = 150
OPLOT, hi_max_3, FINDGEN(jpl+1)*0.+5   , psym = 2, symsize = 1, thick = 2, color = 0
OPLOT, hi_max_4, FINDGEN(jpl1+1)*0.+6   , psym = 2, symsize = 1, thick = 2, color = 150

PLOT, FINDGEN(jpl+1)/jpl, hi_max_cice_5, psym = 4, symsize = 1, thick = 2, charsize = cs
OPLOT, FINDGEN(jpl1+1)/jpl1, hi_max_cice_10, psym = 4, symsize = 1, thick = 2, color = 150
OPLOT, FINDGEN(jpl+1)/jpl, hi_max_1,   psym = 1, symsize = 1, thick = 2
OPLOT, FINDGEN(jpl1+1)/jpl1, hi_max_2, psym = 1, symsize = 1, thick = 2, color = 150
OPLOT, FINDGEN(jpl+1)/jpl, hi_max_3,   psym = 2, symsize = 1, thick = 2 
OPLOT, FINDGEN(jpl1+1)/jpl1, hi_max_4, psym = 2, symsize = 1, thick = 2, color = 150


DEVICE, /CLOSE
SET_PLOT, 'X'
!P.MULTI=[0,1,1]
END

