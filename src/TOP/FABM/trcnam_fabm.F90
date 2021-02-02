MODULE trcnam_fabm
   !!======================================================================
   !!                      ***  MODULE trcnam_fabm  ***
   !! TOP :   initialisation of some run parameters for FABM bio-model
   !!======================================================================
   !! History :   1.0  !  2015-04  (PML) Original code
   !! History :   1.1  !  2020-06  (PML) Update to FABM 1.0, improved performance
   !!----------------------------------------------------------------------
   USE trc             ! TOP variables
#if defined key_fabm
   !!----------------------------------------------------------------------
   !!   'key_fabm'   :                                       FABM model
   !!----------------------------------------------------------------------
   !! trc_nam_fabm      : FABM initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters

   USE par_fabm
   USE trcsms_fabm

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_fabm   ! called by trcnam.F90 module
   PUBLIC   trc_nam_fabm_override ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_fabm
      LOGICAL :: l_ext
      INTEGER :: nmlunit, ios
      NAMELIST/namfabm/ nn_adv

      ! Read NEMO-FABM coupler settings from namfabm
      nn_adv = 3
      INQUIRE( FILE='namelist_fabm_ref', EXIST=l_ext )
      IF (l_ext) then
         CALL ctl_opn( nmlunit, 'namelist_fabm_ref', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE.)
         READ(nmlunit, nml=namfabm, iostat=ios)
         IF( ios /= 0 ) CALL ctl_nam ( ios , 'namfabm in namelist_fabm_ref', .TRUE. )
      END IF
      INQUIRE( FILE='namelist_fabm_cfg', EXIST=l_ext )
      IF (l_ext) then
         CALL ctl_opn( nmlunit, 'namelist_fabm_cfg', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE.)
         READ(nmlunit, nml=namfabm, iostat=ios)
         IF( ios /= 0 ) CALL ctl_nam ( ios , 'namfabm in namelist_fabm_cfg', .TRUE. )
      END IF
      IF (nn_adv /= 1 .and. nn_adv /= 3) CALL ctl_stop( 'STOP', 'trc_ini_fabm: nn_adv must be 1 or 3.' )
   END SUBROUTINE trc_nam_fabm

   SUBROUTINE trc_nam_fabm_override(sn_tracer)
      TYPE(PTRACER), DIMENSION(jpmaxtrc), INTENT(INOUT) :: sn_tracer

      INTEGER :: jn
      CHARACTER(LEN=3) :: index

      DO jn=1,jp_fabm
         IF (sn_tracer(jn)%clsname /= 'NONAME' .AND. sn_tracer(jn)%clsname /= model%interior_state_variables(jn)%name) THEN
            WRITE (index,'(i0)') jn
            CALL ctl_stop('Tracer name mismatch in namtrc: '//TRIM(sn_tracer(jn)%clsname)//' found at sn_tracer('//TRIM(index)//') where '//TRIM(model%interior_state_variables(jn)%name)//' was expected.')
         END IF
         sn_tracer(jn)%clsname = TRIM(model%interior_state_variables(jn)%name)
         sn_tracer(jn)%cllname = TRIM(model%interior_state_variables(jn)%long_name)
         sn_tracer(jn)%clunit = TRIM(model%interior_state_variables(jn)%units)
         sn_tracer(jn)%llinit = .FALSE.
      END DO
   END SUBROUTINE trc_nam_fabm_override

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                             No FABM
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_fabm                      ! Empty routine
   END  SUBROUTINE  trc_nam_fabm

   SUBROUTINE trc_nam_fabm_override (dummy)
       TYPE(PTRACER), DIMENSION(jpmaxtrc), INTENT(INOUT), optional :: dummy
   END SUBROUTINE trc_nam_fabm_override
#endif  

   !!======================================================================
END MODULE trcnam_fabm
