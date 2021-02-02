MODULE par_fabm

#if defined key_fabm
#  include "fabm_version.h"
#  if _FABM_API_VERSION_ < 1
#    error You need FABM 1.0 or later
#  endif
   USE fabm
#endif

   IMPLICIT NONE

   INTEGER, PUBLIC :: jp_fabm0, jp_fabm1, jp_fabm, &
                      jp_fabm_surface, jp_fabm_bottom, &
                      jp_fabm_m1

   LOGICAL, PUBLIC, ALLOCATABLE, DIMENSION(:) ::   lk_rad_fabm !: FABM negativity correction flag array

#if defined key_fabm
   CLASS (type_fabm_model), POINTER :: model !FABM model instance

   !!---------------------------------------------------------------------
   !!   'key_fabm'                     FABM tracers
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_fabm     = .TRUE.   !: FABM flag 
#else
   !!---------------------------------------------------------------------
   !!   Default                           No user defined tracers (FABM)
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_fabm     = .FALSE.  !: FABM flag 
#endif

   !!======================================================================
END MODULE par_fabm
