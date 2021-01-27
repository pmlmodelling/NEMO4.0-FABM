MODULE trcwri_fabm
   !!======================================================================
   !!                       *** MODULE trcwri_fabm ***
   !!    fabm :   Output of FABM tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
#if defined key_top && key_fabm && defined key_iomput
   !!----------------------------------------------------------------------
   !!   'key_fabm'                                           FABM model
   !!----------------------------------------------------------------------
   !! trc_wri_fabm   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE trc         ! passive tracers common variables 
   USE iom         ! I/O manager
   USE trcsms_fabm, only: trc_sms_fabm_check_mass
   USE par_fabm
   USE st2d_fabm
   USE fabm, only: fabm_get_bulk_diagnostic_data, fabm_get_horizontal_diagnostic_data

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_fabm 

!jth #  include "top_substitute.h90"
CONTAINS

   SUBROUTINE trc_wri_fabm
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      INTEGER              :: jn
      !!---------------------------------------------------------------------
 
      ! write the tracer concentrations in the file
      ! ---------------------------------------
      DO jn = 1, jp_fabm
         CALL iom_put( model%state_variables(jn)%name, trn(:,:,:,jp_fabm0+jn-1) )
      END DO
      DO jn = 1, jp_fabm_surface
         CALL iom_put( model%surface_state_variables(jn)%name, fabm_st2dn(:,:,jn) )
      END DO
      DO jn = 1, jp_fabm_bottom
         CALL iom_put( model%bottom_state_variables(jn)%name, fabm_st2dn(:,:,jp_fabm_surface+jn) )
      END DO


      ! write 3D diagnostics in the file
      ! ---------------------------------------
      DO jn = 1, size(model%diagnostic_variables)
         IF (model%diagnostic_variables(jn)%save) &
             CALL iom_put( model%diagnostic_variables(jn)%name, fabm_get_bulk_diagnostic_data(model,jn))
      END DO

      ! write 2D diagnostics in the file
      ! ---------------------------------------
      DO jn = 1, size(model%horizontal_diagnostic_variables)
         IF (model%horizontal_diagnostic_variables(jn)%save) &
             CALL iom_put( model%horizontal_diagnostic_variables(jn)%name, fabm_get_horizontal_diagnostic_data(model,jn))
      END DO
      !

      CALL trc_sms_fabm_check_mass
   END SUBROUTINE trc_wri_fabm

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   PUBLIC trc_wri_fabm
CONTAINS
   SUBROUTINE trc_wri_fabm                     ! Empty routine  
   END SUBROUTINE trc_wri_fabm
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcwri_fabm.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcwri_fabm
