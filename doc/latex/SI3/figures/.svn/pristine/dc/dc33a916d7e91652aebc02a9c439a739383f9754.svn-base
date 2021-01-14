;
;------------------------------------------------------------------------------
;
; Plot LIM thermal properties
; 
;------------------------------------------------------------------------------

; Model parameters
; Sample experiments
N_spl = 3 & Ts_spl = FLTARR(N_spl)
Ts_spl[*] = [ -15., -10., -5. ] & Tw_spl = 2.
colors = [ 50., 100., 150. ]

;
;==============================================================================
; 1) User parameters
;==============================================================================
;

out_dir=''
file_out   = 'Thermal_properties'

numplot_x = 2       ; number of horizontal plots
numplot_y = 2       ; number of vertical plots
ct = 13             ; colortable
cs = 1.0            ; size of fonts on plots (1=normal)
th = 2.             ; thickness of curves
device = 'PS'       ; 'PS' or 'X'

;
;==============================================================================
; 2) Initialize graphics
;==============================================================================
;

dx = 6 
dy = 6.
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
!X.MARGIN   = [6,3]
!Y.MARGIN   = [4,3]

;==============================================================================
; 3) Computations and plots
;==============================================================================

; arrays
N_T = 200
N_S = 3
miss_val = 9.99e32

;--- Physical constants
c0 = 2067.d0
cw = 3991.d0
L  = 334000.d0
mu = 0.054d0
rhoi = 917.d0
rhow = 1026.d0
rhos = 330.d0

Sw = 34.d0

beta1 = 0.09; pringle coefficients
beta2 = 0.011;

;--- Ice salinity
;S = ( rhoi - rhos ) / rhoi * Sw
S = FLTARR(N_S)
Smin = 5.   
Smax = 15.
dS = ( Smax - Smin ) / FLOAT(N_S-1)
S[0:N_S-1] = Smin + FINDGEN(N_S) * dS

Tf = -mu * S

;--- Temperatures
T = FLTARR(N_T)
Tmin = -20.
Tf   = -mu * Smax
Tmax = 0. 
dT = ( Tmax - Tmin ) / FLOAT(N_T-1)
T[0:N_T-1] = Tmin + FINDGEN(N_T) * dT

;--- Grids
S_grid = FLTARR(N_S,N_T)
T_grid = FLTARR(N_S,N_T)

FOR j = 0, N_T - 1 DO BEGIN
   S_grid(0:N_S-1,j) = S(0:N_S-1)
ENDFOR

FOR i = 0, N_S - 1 DO BEGIN
   T_grid(i,0:N_T-1) = T(0:N_T-1)
ENDFOR

;--- Brine fraction
phi = -mu * S_grid / T_grid

;--- specific enthalpy
q_i = c0 * ( T_grid + mu * S_grid ) - L * ( 1 + mu * S_grid / T_grid ) - cw*mu*S_grid

;--- thermal conductivity (pringle)
kpri = 2.11 + beta1 * S_grid / T_grid - beta2 * T_grid; 

;--- specific heat
cp = c0 + L*mu*S_grid/(T_grid*T_grid)

;--- PLOTS
!P.MULTI=[0, numplot_x, numplot_y]
LOADCT, 13

PLOT, [Tmin, Tmax], [0., 50.], /NODATA, XTITLE='T (C)', YTITLE='phi (%)', CHARSIZE = cs
FOR i = 0, N_S - 1 DO BEGIN
   zaddr = where(phi LT 100) 
   OPLOT, T_grid(i,zaddr), phi(i,zaddr)*100., thick = th, color = i* 50
ENDFOR

PLOT, [Tmin, Tmax], [-4., 0.], /NODATA, XTITLE='T (C)', YTITLE='q_m (10^5 J/kg)', CHARSIZE = cs
FOR i = 0, N_S - 1 DO BEGIN
   zaddr = where(phi LT 100) 
   OPLOT, T_grid(i,zaddr), q_i(i,zaddr)/1.0e5, thick = th, color = i* 50
ENDFOR

PLOT, [Tmin, Tmax], [0., 3.], /NODATA, XTITLE='T (C)', YTITLE='k_i (W/m/K)', CHARSIZE = cs
FOR i = 0, N_S - 1 DO BEGIN
   zaddr = where(phi LT 100) 
   OPLOT, T_grid(i,zaddr), kpri(i,zaddr), thick = th, color = i* 50
ENDFOR

PLOT, [Tmin, Tmax], [0., 20.], /NODATA, XTITLE='T (C)', YTITLE='c_i (10^3 J/kg/K)', CHARSIZE = cs
FOR i = 0, N_S - 1 DO BEGIN
   zaddr = where(phi LT 100) 
   OPLOT, T_grid(i,zaddr), cp(i,zaddr)/1000., thick = th, color = i* 50
ENDFOR


;T[0:N_T-1] = FINDGEN(N_T) * dT + Tmin
;dT = ( Tmax - Tf) / FLOAT(N/2-1)
;T[N/2:N-1] = FINDGEN(N/2) * dT + Tf

;;PLOT, T, TITLE = "T", CHARSIZE = cs
;
;; Specific heat
;ci = FLTARR(N) & ci(*) = miss_val
;cis = FLTARR(N) & cis(*) = miss_val
;ci(0:N/2) = c0 + L*mu*S/T(0:N/2)/T(0:N/2)
;cis(0:N/2) = mu*(c0 - L/T(0:N/2))
;;PLOT, T(0:N/2), ci(0:N/2)/L, CHARSIZE = cs
;;OPLOT, T(0:N/2), cis(0:N/2)/L, LINESTYLE = 1
;
;;--- Brine volume
;e_i = FLTARR(N) & e_i(*) = 1.
;i_addr = FINDGEN(N/2+1)
;e_i(i_addr) = -mu * S / T(i_addr)
;PLOT, T, e_i, CHARSIZE = cs, TITLE = "e"
;OPLOT, T, REPLICATE(0.5, N), linestyle = 1
;T_50 = -mu * S / 0.50 ; temperature at which brine volume is 50%
;
;; Notz (2005), Phd thesis, page 44
;i_addr = WHERE(T LT Tf)
;Sbr = - 21.4 * T(i_addr) - 0.886 * T(i_addr) * T(i_addr) - 0.0170 * T(i_addr) * T(i_addr) * T(i_addr)
;zrS = S / Sbr(i_addr)
;beta_ocs = 0.81
;
;zrhoi      = 916.8 - 0.1403 * T(i_addr)              ; pure ice density Pounder(1965)
;zrhol      = rhow + beta_ocs * ( Sbr(i_addr) - Sw )
;zrR        = zrhoi / zrhol
;znum       = 1 - zrS
;zden       = 1 + zrS * ( zrR - 1 )
;e_i2 = FLTARR(N) & e_i2(*) = 1.
;e_i2(i_addr) =  1. - znum / zden
;;OPLOT, T, e_i2, linestyle = 3
;
;;OPLOT, [ T_50, T_50], [ 0., 1.], linestyle = 1
;
;;--- Enthalpies
;Ei = FLTARR(N) & Ei(*) = miss_val
;Ei2= FLTARR(N) & Ei2(*) = miss_val
;Es = FLTARR(N) & Es(*) = miss_val
;Es_wgt = FLTARR(N) & Es(*) = miss_val
;Ew = FLTARR(N) & Ew(*) = miss_val
;
;i_addr = WHERE( T LE 0. )
;Es(i_addr) = c0 * T(i_addr) - L
;Es_wgt(i_addr) = Es(i_addr) * rhos / rhoi
;i_addr = WHERE( T GT 0. )
;Es(i_addr) = cw * T(i_addr)
;Es_wgt(i_addr) = Es(i_addr) * ( rhoi - rhos ) / rhoi
;
;i_addr = WHERE( T GE Tf )
;Ew(i_addr) = cw * T(i_addr)
;i_addr = FINDGEN(N/2+1)
;
;Ei(i_addr) = c0 * ( T(i_addr) + mu * S ) - L * ( 1 + mu * S / T(i_addr) ) - cw*mu*S
;Ei2(i_addr) = - ( 1 - e_i(i_addr) ) * L
;
;Emax = MAX(Ew(WHERE(Ew NE miss_val)))
;
;PLOT, [Tmin, Tmax], [ MIN(Ei)/1000., Emax/1000. ], /NODATA, TITLE = "Enthalpy", CHARSIZE = cs, XTITLE = "degC", YTITLE = "kJ/kg", YMINOR = 2
;OPLOT, T, Ei/1000., MAX=0., THICK = 3.;, solid black
;OPLOT, T, Ei2/1000., MAX=0., COLOR = 150, THICK = 3.; solid grey
;;OPLOT, T, Es/1000., LINESTYLE = 1, THICK = 3.
;OPLOT, T, Es_wgt/1000., COLOR = 100, THICK = 3
;i_addr = WHERE( T GE 0. )
;OPLOT, T(i_addr), Es_wgt(i_addr)/1000., COLOR = 200, THICK = 3
;;OPLOT, T, Ew/1000., MAX=Emax/1000., LINESTYLE = 0, THICK = 3.
;OPLOT, [ T_50, T_50], [ -400., 100.], COLOR = 100
;OPLOT, [ Tmin, Tmax], [0., 0.], COLOR = 100
;
;;--- Sample experiment
;Es_spl = c0 * Ts_spl - L
;Ew_spl = cw * Tw_spl
;
;LOADCT, ct
;FOR i = 0, N_spl - 1 DO BEGIN
;   OPLOT, [ Ts_spl(i) ], [ rhos / rhoi * Es_spl(i) ]/1000., PSYM = 1, SYMSIZE = 2, THICK = 3, COLOR = COLORS(i)
;ENDFOR
;LOADCT, 0
;OPLOT, [ Tw_spl ], [ (rhoi - rhos) / rhoi * Ew_spl ]/1000., PSYM = 1, SYMSIZE = 2, THICK = 3
;
;Ei_spl = Es_spl * rhos / rhoi + Ew_spl * ( rhoi - rhos ) / rhoi
;
;; Retrieve T
;zA = ( L *( rhos/rhoi ) - (rhos/rhoi) * c0 * Ts_spl + (rhos-rhoi)/rhoi * cw * Tw_spl ) * rhoi
;
;aaa = c0
;bbb = mu * S * ( c0 - cw ) - L + zA / rhoi
;ccc = - L * mu * S
;
;discr = SQRT ( bbb * bbb - 4 * aaa * ccc )
;
;Ti_spl = ( - bbb - discr ) / 2.0 / aaa
;
;Ei_test = c0 * ( Ti_spl + mu * S ) - L * ( 1 + mu * S / Ti_spl ) - cw*mu*S
;PRINT, "Ei_test : ", Ei_test
;PRINT, "Ei_spl  : ", Ei_spl
;
;; Retrieve brine volume
;e_i_spl = -mu*S/Ti_spl
;PRINT, "e_i_spl = ", e_i_spl
;
;LOADCT, ct
;FOR i = 0, N_spl - 1 DO BEGIN
;   OPLOT, [ Ti_spl(i) ], [ Ei_spl(i) ]/1000., PSYM = 1, SYMSIZE = 2, COLOR = COLORS(i), THICK = 3
;ENDFOR
;LOADCT, 0

;--------------------------------------------------------------------------------------------------
DEVICE, /CLOSE
SET_PLOT, 'X'
END
