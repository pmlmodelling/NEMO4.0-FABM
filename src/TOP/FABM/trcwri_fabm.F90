MODULE trcwri_fabm
   !!======================================================================
   !!                       *** MODULE trcwri_fabm ***
   !!    fabm :   Output of FABM tracers
   !!======================================================================
   !! History :   1.0  !  2015-04  (PML) Original code
   !! History :   1.1  !  2020-06  (PML) Update to FABM 1.0, improved performance
   !!----------------------------------------------------------------------
#if defined key_top && key_fabm && defined key_iomput
   !!----------------------------------------------------------------------
   !!   'key_fabm'                                           FABM model
   !!----------------------------------------------------------------------
   !! trc_wri_fabm   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE trc         ! passive tracers common variables 
   USE iom         ! I/O manager
   USE trdtrc_oce
   USE trcsms_fabm, only: trc_sms_fabm_check_mass
   USE par_fabm
   USE st2d_fabm

   IMPLICIT NONE
   PRIVATE

#if defined key_tracer_budget
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:), SAVE :: tr_temp
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: fabm_st2d_temp
#endif

   INTERFACE trc_wri_fabm
       MODULE PROCEDURE wri_fabm,wri_fabm_fl
   END INTERFACE trc_wri_fabm

   PUBLIC trc_wri_fabm 

#  include "top_substitute.h90"
CONTAINS

   SUBROUTINE wri_fabm_fl(kt,fl)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in )               :: fl
      INTEGER, INTENT( in )               :: kt

#if defined key_tracer_budget
      INTEGER              :: jn
      CHARACTER (len=20)   :: cltra
      REAL(wp), DIMENSION(jpi,jpj,jpk)    :: trpool !temporary storage pool 3D
      REAL(wp), DIMENSION(jpi,jpj)    :: st2dpool !temporary storage pool 2D
      !!---------------------------------------------------------------------
 
      ! write the tracer concentrations in the file
      ! ---------------------------------------
! depth integrated
! for strict budgetting write this out at end of timestep as an average between 'now' and 'after' at kt
      DO jn = 1, jp_fabm
        IF(ln_trdtrc (jp_fabm_m1+jn))THEN
         trpool(:,:,:) = 0.5 * ( trn(:,:,:,jp_fabm_m1+jn)*fse3t_a(:,:,:) + &
                             tr_temp(:,:,:,jn)*fse3t(:,:,:) )
         cltra = TRIM( model%interior_state_variables(jn)%name )//"_e3t"     ! depth integrated output
         IF( kt == nittrc000 ) write(6,*)'output pool ',cltra
         CALL iom_put( cltra, trpool)
        ENDIF
      END DO
#else
      CONTINUE
#endif

   END SUBROUTINE wri_fabm_fl

   SUBROUTINE wri_fabm(kt)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in )               :: kt
      INTEGER              :: jn, jk
      REAL(wp), DIMENSION(jpi,jpj)    :: vint

#if defined key_tracer_budget
      IF( kt == nittrc000 ) THEN
         ALLOCATE(tr_temp(jpi,jpj,jpk,jp_fabm),fabm_st2d_temp(jpi,jpj,jp_fabm_surface+jp_fabm_bottom))
      ENDIF
      tr_temp(:,:,:,:)=trn(:,:,:,jp_fabm0:jp_fabm1) ! slwa save for tracer budget (unfiltered trn)
      fabm_st2d_temp(:,:,:)=fabm_st2dn(:,:,:)
#endif
      DO jn = 1, jp_fabm
         ! Save 3D field
         CALL iom_put(model%interior_state_variables(jn)%name, trn(:,:,:,jp_fabm_m1+jn))

         ! Save depth integral if selected for output in XIOS
         IF (iom_use(TRIM(model%interior_state_variables(jn)%name)//'_VINT')) THEN
            vint = 0._wp
            DO jk = 1, jpkm1
               vint = vint + trn(:,:,jk,jp_fabm_m1+jn) * fse3t(:,:,jk) * tmask(:,:,jk)
            END DO
            CALL iom_put(TRIM(model%interior_state_variables(jn)%name)//'_VINT', vint)
         END IF
      END DO
      DO jn = 1, jp_fabm_surface
         CALL iom_put( model%surface_state_variables(jn)%name, fabm_st2dn(:,:,jn) )
      END DO
      DO jn = 1, jp_fabm_bottom
         CALL iom_put( model%bottom_state_variables(jn)%name, fabm_st2dn(:,:,jp_fabm_surface+jn) )
      END DO

      CALL trc_sms_fabm_check_mass

   END SUBROUTINE wri_fabm

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   INTERFACE trc_wri_fabm
       MODULE PROCEDURE wri_fabm,wri_fabm_fl
   END INTERFACE trc_wri_fabm

   PUBLIC trc_wri_fabm

   CONTAINS

   SUBROUTINE wri_fabm_fl(kt,fl)
      INTEGER, INTENT( in )               :: fl
      INTEGER, INTENT( in )               :: kt
   END SUBROUTINE wri_fabm_fl

   SUBROUTINE wri_fabm(kt)                 ! Empty routine  
      INTEGER, INTENT( in )               :: kt
   END SUBROUTINE wri_fabm

#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcwri_fabm.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcwri_fabm
