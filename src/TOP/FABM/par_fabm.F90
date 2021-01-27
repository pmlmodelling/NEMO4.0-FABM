MODULE par_fabm

   USE fabm

   IMPLICIT NONE

   TYPE (type_model) :: model !FABM model instance

   INTEGER, PUBLIC :: jp_fabm0, jp_fabm1, jp_fabm, &
                      jp_fabm_surface, jp_fabm_bottom, &
                      jp_fabm_m1

#if defined key_fabm
   !!---------------------------------------------------------------------
   !!   'key_fabm'                     FABM tracers
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC ::   ln_fabm     = .TRUE.   !: FABM flag 
   LOGICAL, PUBLIC, ALLOCATABLE, DIMENSION(:) ::   lk_rad_fabm !: FABM negativity correction flag array 
#else
   !!---------------------------------------------------------------------
   !!   Default                           No user defined tracers (FABM)
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   ln_fabm     = .FALSE.  !: FABM flag 
#endif

   !!======================================================================
END MODULE par_fabm
