PRO std_main_last_year_5D
@initenv
;
  type = getenv('PLOTTYPE')
  CASE type OF
    '':print, 'The environment variable PLOTTYPE is not defined. We stop'
    'ts':std_ts_all, /postscript
    'plot':std_plot_all_last_year_5D, /postscript
    ELSE:print, 'Wrong definition of the environment variable PLOTTYPE. We stop'
  ENDCASE
;
END
