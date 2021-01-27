MODULE trcsms_fabm
   !!======================================================================
   !!                         ***  MODULE trcsms_fabm  ***
   !! TOP :   Main module of the FABM tracers
   !!======================================================================
   !! History :   1.0  !  2015-04  (PML) Original code
   !!----------------------------------------------------------------------
#if defined key_fabm
   !!----------------------------------------------------------------------
   !!   'key_fabm'                                               FABM tracers
   !!----------------------------------------------------------------------
   !! trc_sms_fabm       : FABM model main routine
   !! trc_sms_fabm_alloc : allocate arrays specific to FABM sms
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trcbc
   USE trd_oce
   USE trdtrc
#if defined key_trdtrc && defined key_iomput
   USE trdtrc_oce
#endif

   USE oce, only: tsn  ! Needed?
   USE sbc_oce, only: lk_oasis
   USE dom_oce
   USE zdf_oce
   !USE iom
   USE xios
   USE cpl_oasis3
   USE st2D_fabm
   USE inputs_fabm
   USE vertical_movement_fabm

   !USE fldread         !  time interpolation

   USE fabm

   IMPLICIT NONE

!jth#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"

   PRIVATE

   PUBLIC   trc_sms_fabm       ! called by trcsms.F90 module
   PUBLIC   trc_sms_fabm_alloc ! called by trcini_fabm.F90 module
   PUBLIC   trc_sms_fabm_check_mass
   PUBLIC   st2d_fabm_nxt ! 2D state intergration
   PUBLIC   compute_fabm ! Compute FABM sources, sinks and diagnostics

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: flux    ! Cross-interface flux of pelagic variables (# m-2 s-1)

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)   :: ext     ! Light extinction coefficient (m-1)

   ! Work array for mass aggregation
   REAL(wp), ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)   :: current_total


   ! Arrays for environmental variables
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:,:) :: prn,rho
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:) :: taubot

   ! repair counters
   INTEGER :: repair_interior_count,repair_surface_count,repair_bottom_count

   ! state check type
   TYPE type_state
      LOGICAL             :: valid
      LOGICAL             :: repaired
   END TYPE

   REAL(wp), PUBLIC :: daynumber_in_year

   TYPE (type_bulk_variable_id),SAVE :: swr_id

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_fabm( kt )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_fabm  ***
      !!
      !! ** Purpose :   main routine of FABM model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER :: jn
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztrfabm

!!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sms_fabm')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,'(a,i0,a,i4.4,a,i2.2,a,i2.2,a,i5,a)') &
          ' trc_sms_fabm:  FABM model, iteration ',kt,' ', &
          nyear,'-',nmonth,'-',nday,' ',nsec_day," secs"
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      CALL update_inputs( kt )

      CALL compute_fabm

      CALL compute_vertical_movement( kt )

      CALL st2d_fabm_nxt( kt )

      IF( l_trdtrc )  CALL wrk_alloc( jpi, jpj, jpk, ztrfabm )
!jth
!      CALL trc_bc_read  ( kt )       ! tracers: surface and lateral Boundary Conditions
       CALL trc_bc       ( kt )       ! tracers: surface and lateral Boundary Conditions
      CALL trc_rnf_fabm ( kt ) ! River forcings

      IF( l_trdtrc ) THEN      ! Save the trends in the mixed layer
          DO jn = jp_fabm0, jp_fabm1
            ztrfabm(:,:,:) = tra(:,:,:,jn)
            CALL trd_trc( ztrfabm, jn, jptra_sms, kt )   ! save trends
          END DO
          CALL wrk_dealloc( jpi, jpj, jpk, ztrfabm )
      END IF
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_sms_fabm')

   END SUBROUTINE trc_sms_fabm

   SUBROUTINE compute_fabm()
      INTEGER :: ji,jj,jk,jn
      TYPE(type_state) :: valid_state
      REAL(wp) :: zalfg,zztmpx,zztmpy

      ! Validate current model state (setting argument to .TRUE. enables repair=clipping)
      valid_state = check_state(.TRUE.)
      IF (.NOT.valid_state%valid) THEN
         WRITE(numout,*) "Invalid value in FABM encountered in area ",narea,"!!!"
#if defined key_iomput
         CALL xios_finalize                ! end mpp communications with xios
         IF( lk_oasis ) CALL cpl_finalize    ! end coupling and mpp communications with OASIS
#else
         IF( lk_oasis ) THEN
            CALL cpl_finalize              ! end coupling and mpp communications with OASIS
         ELSE
            IF( lk_mpp )   CALL mppstop    ! end mpp communications
         ENDIF
#endif
      END IF
      IF (valid_state%repaired) THEN
         WRITE(numout,*) "Total interior repairs up to now on process",narea,":",repair_interior_count
         WRITE(numout,*) "Total surface repairs up to now on process",narea,":",repair_surface_count
         WRITE(numout,*) "Total bottom repairs up to now on process",narea,":",repair_bottom_count
      ENDIF

      ! Compute the now hydrostatic pressure
      ! copied from istate.F90
      ! ------------------------------------

      zalfg = 0.5e-4 * grav ! FABM wants dbar, convert from Pa

      rho = rau0 * ( 1. + rhd )

      prn(:,:,1) = 10.1325 + zalfg * e3t_n(:,:,1) * rho(:,:,1)

      daynumber_in_year=(fjulday-fjulstartyear+1)*1._wp

      DO jk = 2, jpk                                              ! Vertical integration from the surface
         prn(:,:,jk) = prn(:,:,jk-1) + zalfg * ( &
                     e3t_n(:,:,jk-1) * rho(:,:,jk-1)  &
                     + e3t_n(:,:,jk) * rho(:,:,jk) )
      END DO

      ! Bottom stress
      taubot(:,:) = 0._wp
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
               zztmpx = (  bfrua(ji  ,jj) * un(ji  ,jj,mbku(ji  ,jj))  &
                      &  + bfrua(ji-1,jj) * un(ji-1,jj,mbku(ji-1,jj))  )
               zztmpy = (  bfrva(ji,  jj) * vn(ji,jj  ,mbkv(ji,jj  ))  &
                      &  + bfrva(ji,jj-1) * vn(ji,jj-1,mbkv(ji,jj-1))  )
               taubot(ji,jj) = 0.5_wp * rau0 * SQRT( zztmpx * zztmpx + zztmpy * zztmpy ) * tmask(ji,jj,1)
               !
         ENDDO
      ENDDO
      ! Compute light extinction
      DO jk=1,jpk
          DO jj=1,jpj
            call fabm_get_light_extinction(model,1,jpi,jj,jk,ext)
         END DO
      END DO

      ! Compute light field (stored among FABM's internal diagnostics)
      DO jj=1,jpj
          DO ji=1,jpi
            call fabm_get_light(model,1,jpk,ji,jj)
         END DO
      END DO

      ! TODO: retrieve 3D shortwave and store in etot3

      ! Zero rate array of interface-attached state variables
      fabm_st2Da = 0._wp

      ! Compute interfacial source terms and fluxes
      DO jj=1,jpj
         ! Process bottom (fabm_do_bottom increments rather than sets, so zero flux array first)
         flux = 0._wp
         CALL fabm_do_bottom(model,1,jpi,jj,flux,fabm_st2Da(:,jj,jp_fabm_surface+1:))
         DO jn=1,jp_fabm
             DO ji=1,jpi
                 ! Divide bottom fluxes by height of bottom layer and add to source terms.
                 ! TODO: is there perhaps an existing variable for e3t_n(ji,jj,mbkt(ji,jj))??
                 tra(ji,jj,mbkt(ji,jj),jp_fabm_m1+jn) = tra(ji,jj,mbkt(ji,jj),jp_fabm_m1+jn) + flux(ji,jn)/e3t_n(ji,jj,mbkt(ji,jj))
             END DO
         END DO

         ! Process surface (fabm_do_surface increments rather than sets, so zero flux array first)
         flux = 0._wp
         CALL fabm_do_surface(model,1,jpi,jj,flux,fabm_st2Da(:,jj,1:jp_fabm_surface))
         DO jn=1,jp_fabm
             ! Divide surface fluxes by height of surface layer and add to source terms.
             tra(:,jj,1,jp_fabm_m1+jn) = tra(:,jj,1,jp_fabm_m1+jn) + flux(:,jn)/e3t_n(:,jj,1)
         END DO
      END DO

      ! Compute interior source terms (NB fabm_do increments rather than sets)
      DO jk=1,jpk
          DO jj=1,jpj
              CALL fabm_do(model,1,jpi,jj,jk,tra(:,jj,jk,jp_fabm0:jp_fabm1))
          END DO
      END DO
   END SUBROUTINE compute_fabm

   FUNCTION check_state(repair) RESULT(exit_state)
      LOGICAL, INTENT(IN) :: repair
      TYPE(type_state) :: exit_state

      INTEGER             :: jj,jk
      LOGICAL             :: valid_int,valid_sf,valid_bt

      exit_state%valid = .TRUE.
      exit_state%repaired =.FALSE.
      DO jk=1,jpk
         DO jj=1,jpj
            CALL fabm_check_state(model,1,jpi,jj,jk,repair,valid_int)
            IF (repair.AND..NOT.valid_int) THEN
               repair_interior_count = repair_interior_count + 1
               exit_state%repaired = .TRUE.
            END IF
            IF (.NOT.(valid_int.OR.repair)) exit_state%valid = .FALSE.
         END DO
      END DO
      DO jj=1,jpj
         CALL fabm_check_surface_state(model,1,jpi,jj,repair,valid_sf)
         IF (repair.AND..NOT.valid_sf) THEN
            repair_surface_count = repair_surface_count + 1
            exit_state%repaired = .TRUE.
         END IF
         IF (.NOT.(valid_sf.AND.valid_bt).AND..NOT.repair) exit_state%valid = .FALSE.
         CALL fabm_check_bottom_state(model,1,jpi,jj,repair,valid_bt)
         IF (repair.AND..NOT.valid_bt) THEN
            repair_bottom_count = repair_bottom_count + 1
            exit_state%repaired = .TRUE.
         END IF
         IF (.NOT.(valid_sf.AND.valid_bt).AND..NOT.repair) exit_state%valid = .FALSE.
      END DO
   END FUNCTION

   SUBROUTINE trc_sms_fabm_check_mass()
      REAL(wp) :: total(SIZE(model%conserved_quantities))
      INTEGER :: jk,jj,jn

      total = 0._wp

      DO jk=1,jpk
         DO jj=1,jpj
            CALL fabm_get_conserved_quantities(model,1,jpi,jj,jk,current_total)
            DO jn=1,SIZE(model%conserved_quantities)
               total(jn) = total(jn) + SUM(cvol(:,jj,jk)*current_total(:,jn)*tmask_i(:,jj))
            END DO
         END DO
      END DO

      DO jj=1,jpj
         CALL fabm_get_horizontal_conserved_quantities(model,1,jpi,jj,current_total)
         DO jn=1,SIZE(model%conserved_quantities)
            total(jn) = total(jn) + SUM(e1e2t(:,jj)*current_total(:,jn)*tmask_i(:,jj))
         END DO
      END DO

      IF( lk_mpp ) CALL mpp_sum(total,SIZE(model%conserved_quantities))

      DO jn=1,SIZE(model%conserved_quantities)
         IF(lwp) WRITE(numout,*) 'FABM '//TRIM(model%conserved_quantities(jn)%name),total(jn),TRIM(model%conserved_quantities(jn)%units)//'*m3'
      END DO

   END SUBROUTINE trc_sms_fabm_check_mass

   SUBROUTINE st2d_fabm_nxt( kt )
      !!----------------------------------------------------------------------
      !!                     ***  st2d_fabm_nxt  ***
      !!
      !! ** Purpose :   routine to integrate 2d states in time
      !!
      !! ** Method  :   based on integration of 3D passive tracer fields
      !!                implemented in TOP_SRC/TRP/trcnxt.F90, plus
      !!                tra_nxt_fix in OPA_SRC/TRA/tranxt.F90. Similar to
      !!                time integration of sea surface height in
      !!                OPA_SRC/DYN/sshwzv.F90.
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      REAL(wp) :: z2dt
      INTEGER :: jn

!!----------------------------------------------------------------------
      !
      IF( neuler == 0 .AND. kt == nittrc000 ) THEN
          z2dt = rdt                  ! set time step size (Euler)
      ELSE
          z2dt = 2._wp * rdt          ! set time step size (Leapfrog)
      ENDIF

      ! Forward Euler time step to compute "now"
      DO jn=1,jp_fabm_surface+jp_fabm_bottom
         fabm_st2Da(:,:,jn) = (fabm_st2db(:,:,jn) + z2dt * fabm_st2da(:,:,jn)) * tmask(:,:,1)
      ENDDO

      IF( neuler == 0 .AND. kt == nittrc000 ) THEN        ! Euler time-stepping at first time-step
         !                                                ! (only swap)
         fabm_st2Dn(:,:,:) = fabm_st2Da(:,:,:)
         !
      ELSE
         ! Update now state + Asselin filter time stepping
         fabm_st2Db(:,:,:) = (1._wp - 2._wp*atfp) * fabm_st2Dn(:,:,:) + &
             atfp * ( fabm_st2Db(:,:,:) + fabm_st2Da(:,:,:) )
         fabm_st2Dn(:,:,:) = fabm_st2Da(:,:,:)
      ENDIF

   END SUBROUTINE st2d_fabm_nxt

   INTEGER FUNCTION trc_sms_fabm_alloc()
      INTEGER :: jj,jk,jn
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_fabm_alloc  ***
      !!----------------------------------------------------------------------
      !
      ! ALLOCATE here the arrays specific to FABM
      ALLOCATE( lk_rad_fabm(jp_fabm))
      ALLOCATE( prn(jpi, jpj, jpk))
      ALLOCATE( rho(jpi, jpj, jpk))
      ALLOCATE( taubot(jpi, jpj))
      ! ALLOCATE( tab(...) , STAT=trc_sms_fabm_alloc )

      ! Allocate arrays to hold state for surface-attached and bottom-attached state variables
      ALLOCATE(fabm_st2Dn(jpi, jpj, jp_fabm_surface+jp_fabm_bottom))
      ALLOCATE(fabm_st2Da(jpi, jpj, jp_fabm_surface+jp_fabm_bottom))
      ALLOCATE(fabm_st2Db(jpi, jpj, jp_fabm_surface+jp_fabm_bottom))

      ! Work array to hold surface and bottom fluxes
      ALLOCATE(flux(jpi,jp_fabm))

      ! Work array to hold extinction coefficients
      ALLOCATE(ext(jpi))
      ext=0._wp

      ! Allocate work arrays for vertical movement
      ALLOCATE(w_ct(jpi,jpk,jp_fabm))
      ALLOCATE(w_if(jpk,jp_fabm))
      ALLOCATE(zwgt_if(jpk,jp_fabm))
      ALLOCATE(flux_if(jpk,jp_fabm))
      ALLOCATE(current_total(jpi,SIZE(model%conserved_quantities)))
#if defined key_trdtrc && defined key_iomput
      IF( lk_trdtrc ) ALLOCATE(tr_vmv(jpi,jpj,jpk,jp_fabm))
      IF( lk_trdtrc ) ALLOCATE(tr_inp(jpi,jpj,jpk))
#endif

      trc_sms_fabm_alloc = 0      ! set to zero if no array to be allocated
      !
      IF( trc_sms_fabm_alloc /= 0 ) CALL ctl_warn('trc_sms_fabm_alloc : failed to allocate arrays')
      !

      ! Make FABM aware of diagnostics that are not needed [not included in output]
      DO jn=1,size(model%diagnostic_variables)
          !model%diagnostic_variables(jn)%save = iom_use(model%diagnostic_variables(jn)%name)
      END DO
      DO jn=1,size(model%horizontal_diagnostic_variables)
          !model%horizontal_diagnostic_variables(jn)%save = iom_use(model%horizontal_diagnostic_variables(jn)%name)
      END DO

      ! Provide FABM with domain extents [after this, the save attribute of diagnostic variables can no longe change!]
      call fabm_set_domain(model,jpi, jpj, jpk)

      ! Provide FABM with the vertical indices of the surface and bottom, and the land-sea mask.
      call model%set_bottom_index(mbkt)  ! NB mbkt extents should match dimension lengths provided to fabm_set_domain
      call model%set_surface_index(1)
      call fabm_set_mask(model,tmask,tmask(:,:,1)) ! NB tmask extents should match dimension lengths provided to fabm_set_domain

      ! Send pointers to state data to FABM
      do jn=1,jp_fabm
         call fabm_link_bulk_state_data(model,jn,trn(:,:,:,jp_fabm_m1+jn))
      end do
      DO jn=1,jp_fabm_surface
         CALL fabm_link_surface_state_data(model,jn,fabm_st2Dn(:,:,jn))
      END DO
      DO jn=1,jp_fabm_bottom
         CALL fabm_link_bottom_state_data(model,jn,fabm_st2Dn(:,:,jp_fabm_surface+jn))
      END DO

      ! Send pointers to environmental data to FABM
      call fabm_link_bulk_data(model,standard_variables%temperature,tsn(:,:,:,jp_tem))
      call fabm_link_bulk_data(model,standard_variables%practical_salinity,tsn(:,:,:,jp_sal))
      call fabm_link_bulk_data(model,standard_variables%density,rho(:,:,:))
      call fabm_link_bulk_data(model,standard_variables%pressure,prn)
      call fabm_link_horizontal_data(model,standard_variables%bottom_stress,taubot(:,:))
      ! correct target for cell thickness depends on NEMO configuration:
!jth not used any more
!#ifdef key_vvl
      call fabm_link_bulk_data(model,standard_variables%cell_thickness,e3t_n)
!#else
!      call fabm_link_bulk_data(model,standard_variables%cell_thickness,e3t_0)
!#endif
      call fabm_link_horizontal_data(model,standard_variables%latitude,gphit)
      call fabm_link_horizontal_data(model,standard_variables%longitude,glamt)
      call fabm_link_scalar_data(model,standard_variables%number_of_days_since_start_of_the_year,daynumber_in_year)
      call fabm_link_horizontal_data(model,standard_variables%wind_speed,wndm(:,:))
      call fabm_link_horizontal_data(model,standard_variables%surface_downwelling_shortwave_flux,qsr(:,:))
      call fabm_link_horizontal_data(model,standard_variables%bottom_depth_below_geoid,ht_0(:,:))

      swr_id = model%get_bulk_variable_id(standard_variables%downwelling_shortwave_flux)

      ! Obtain user-specified input variables (read from NetCDF file)
      call link_inputs
      call update_inputs( nit000, .false. )

      ! Check whether FABM has all required data
      call fabm_check_ready(model)

      ! Initialize state
      DO jj=1,jpj
         CALL fabm_initialize_surface_state(model,1,jpi,jj)
         CALL fabm_initialize_bottom_state(model,1,jpi,jj)
      END DO
      DO jk=1,jpk
         DO jj=1,jpj
            CALL fabm_initialize_state(model,1,jpi,jj,jk)
         END DO
      END DO

      ! Set mask for negativity corrections to the relevant states
      DO jn=1,jp_fabm
        IF (model%state_variables(jn)%minimum.ge.0) THEN
          lk_rad_fabm(jn)=.TRUE.
          IF(lwp) WRITE(numout,*) 'FABM clipping for '//TRIM(model%state_variables(jn)%name)//' activated.'
        END IF
      END DO

      ! Mask land points in states with zeros, not nice, but coherent
      ! with NEMO "convention":
      DO jn=jp_fabm0,jp_fabm1
        WHERE (tmask==0._wp)
          trn(:,:,:,jn)=0._wp
        END WHERE
      END DO
      DO jn=1,jp_fabm_surface+jp_fabm_bottom
        WHERE (tmask(:,:,1)==0._wp)
          fabm_st2Dn(:,:,jn)=0._wp
        END WHERE
      END DO

      ! Copy initial condition for interface-attached state variables to "previous" state field
      ! NB NEMO does this itself for pelagic state variables (trb) in TOP_SRC/trcini.F90.
      fabm_st2Db = fabm_st2Dn

      ! Initialise repair counters
      repair_interior_count = 0
      repair_surface_count = 0
      repair_bottom_count = 0

   END FUNCTION trc_sms_fabm_alloc

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No FABM model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_fabm( kt )             ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms_fabm: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_fabm
#endif

   !!======================================================================
END MODULE trcsms_fabm
