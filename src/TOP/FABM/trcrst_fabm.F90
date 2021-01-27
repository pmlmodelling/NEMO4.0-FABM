MODULE trcrst_fabm
   !!======================================================================
   !!                      ***  MODULE trcrst_fabm  ***
   !! Read and write additional restart fields used by FABM
   !!======================================================================
   !! History :
   !!----------------------------------------------------------------------
#if defined key_fabm
   !!----------------------------------------------------------------------
   !!   'key_fabm'   :                                       FABM model
   !!----------------------------------------------------------------------
   !! trc_nam_fabm      : FABM initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE iom

   USE par_fabm
   USE trcsms_fabm
   USE st2D_fabm

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_rst_read_fabm   ! called by trcrst.F90 module
   PUBLIC   trc_rst_wri_fabm    ! called by trcrst.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_rst_read_fabm
      INTEGER :: jn

      DO jn=1,jp_fabm_surface
         CALL iom_get( numrtr, jpdom_autoglo, 'fabm_st2Db'//TRIM(model%surface_state_variables(jn)%name), fabm_st2Db(:,:,jn) )
      END DO
      DO jn=1,jp_fabm_bottom
         CALL iom_get( numrtr, jpdom_autoglo, 'fabm_st2Db'//TRIM(model%bottom_state_variables(jn)%name), fabm_st2Db(:,:,jp_fabm_surface+jn) )
      END DO

      DO jn=1,jp_fabm_surface
         CALL iom_get( numrtr, jpdom_autoglo, 'fabm_st2Dn'//TRIM(model%surface_state_variables(jn)%name), fabm_st2Dn(:,:,jn) )
      END DO
      DO jn=1,jp_fabm_bottom
         CALL iom_get( numrtr, jpdom_autoglo, 'fabm_st2Dn'//TRIM(model%bottom_state_variables(jn)%name), fabm_st2Dn(:,:,jp_fabm_surface+jn) )
      END DO
   END SUBROUTINE trc_rst_read_fabm

   SUBROUTINE trc_rst_wri_fabm(kt)
      INTEGER, INTENT( in ) ::   kt    ! ocean time-step index

      INTEGER :: jn

      DO jn=1,jp_fabm_surface
         CALL iom_rstput( kt, nitrst, numrtw, 'fabm_st2Db'//TRIM(model%surface_state_variables(jn)%name), fabm_st2Db(:,:,jn) )
      END DO
      DO jn=1,jp_fabm_bottom
         CALL iom_rstput( kt, nitrst, numrtw, 'fabm_st2Db'//TRIM(model%bottom_state_variables(jn)%name), fabm_st2Db(:,:,jp_fabm_surface+jn) )
      END DO

      DO jn=1,jp_fabm_surface
         CALL iom_rstput( kt, nitrst, numrtw, 'fabm_st2Dn'//TRIM(model%surface_state_variables(jn)%name), fabm_st2Dn(:,:,jn) )
      END DO
      DO jn=1,jp_fabm_bottom
         CALL iom_rstput( kt, nitrst, numrtw, 'fabm_st2Dn'//TRIM(model%bottom_state_variables(jn)%name), fabm_st2Dn(:,:,jp_fabm_surface+jn) )
      END DO

   END SUBROUTINE trc_rst_wri_fabm

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                             No FABM
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_rst_read_fabm
   END  SUBROUTINE trc_rst_read_fabm

   SUBROUTINE trc_rst_wri_fabm(kt)
      INTEGER, INTENT( in ) ::   kt    ! ocean time-step index
   END SUBROUTINE trc_rst_wri_fabm
#endif

   !!======================================================================
END MODULE trcrst_fabm
