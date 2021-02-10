MODULE trcwri
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    TOP :   Output of passive tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
#if defined key_top && defined key_iomput
   !!----------------------------------------------------------------------
   !!   'key_top'                                           TOP models
   !!----------------------------------------------------------------------
   !! trc_wri_trc   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE dom_oce     ! ocean space and time domain variables
   USE oce_trc     ! shared variables between ocean and passive tracers
   USE trc         ! passive tracers common variables 
   USE iom         ! I/O manager
   USE dianam      ! Output file name
   USE trcwri_pisces
   USE trcwri_cfc
   USE trcwri_c14
   USE trcwri_age
   USE trcwri_my_trc
   ! +++>>> FABM
   USE trcwri_fabm
   ! FABM <<<+++

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri      

CONTAINS

#if defined key_tracer_budget
   SUBROUTINE trc_wri( kt , fl)  !slwa
#else
   SUBROUTINE trc_wri( kt )
#endif
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri  ***
      !! 
      !! ** Purpose :   output passive tracers fields and dynamical trends
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in )     :: kt
      ! +++>>>FABM
#if defined key_tracer_budget
      INTEGER, INTENT( in ), OPTIONAL     :: fl  ! slwa
#endif
      ! FABM <<<+++

      !
      INTEGER                   :: jn
      CHARACTER (len=20)        :: cltra
      CHARACTER (len=40)        :: clhstnam
      INTEGER ::   inum = 11            ! temporary logical unit
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_wri')
      !
      IF( l_offline ) THEN    ! WRITE root name in date.file for use by postpro
         IF(  kt == nittrc000 .AND. lwp ) THEN    ! WRITE root name in date.file for use by postpro
           CALL dia_nam( clhstnam, nn_writetrc,' ' )
           CALL ctl_opn( inum, 'date.file', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
           WRITE(inum,*) clhstnam
           CLOSE(inum)
        ENDIF
        ! Output of initial vertical scale factor
        CALL iom_put("e3t_0", e3t_0(:,:,:) )
        CALL iom_put("e3u_0", e3u_0(:,:,:) )
        CALL iom_put("e3v_0", e3v_0(:,:,:) )
        !
        CALL iom_put( "e3t" , e3t_n(:,:,:) )
        CALL iom_put( "e3u" , e3u_n(:,:,:) )
        CALL iom_put( "e3v" , e3v_n(:,:,:) )
        !
      ENDIF
      ! write the tracer concentrations in the file
      ! ---------------------------------------
      IF( ln_pisces  )   CALL trc_wri_pisces     ! PISCES 
      IF( ll_cfc     )   CALL trc_wri_cfc        ! surface fluxes of CFC
      IF( ln_c14     )   CALL trc_wri_c14        ! surface fluxes of C14
      IF( ln_age     )   CALL trc_wri_age        ! AGE tracer
#if defined key_tracer_budget
      ! +++>>>FABM
      IF( PRESENT(fl) ) THEN
         IF( ln_fabm    )   CALL trc_wri_fabm (kt, fl) ! MY_TRC  tracers for budget
!         IF( ln_my_trc ) CALL trc_wri_my_trc (kt, fl)    ! MY_TRC  tracers for budget
      ELSE
         IF( ln_fabm    )   CALL trc_wri_fabm (kt) ! FABM  tracers for budget
!         IF( ln_my_trc  )   CALL trc_wri_my_trc (kt) ! MY_TRC  tracers
      ENDIF
#else
      IF( ln_fabm  )   CALL trc_wri_fabm (kt)     ! FABM  tracers
!      IF( ln_my_trc  )   CALL trc_wri_my_trc(kt)     ! MY_TRC  tracers
      ! FABM <<<+++
      !
#endif
      IF( ln_timing )   CALL timing_stop('trc_wri')
      !
   END SUBROUTINE trc_wri

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   PUBLIC trc_wri
CONTAINS
   SUBROUTINE trc_wri( kt )                     ! Empty routine   
   INTEGER, INTENT(in) :: kt
   END SUBROUTINE trc_wri
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcwri.F90 12280 2019-12-21 10:42:44Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!======================================================================
END MODULE trcwri
