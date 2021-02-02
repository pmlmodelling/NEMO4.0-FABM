MODULE trcsms_fabm
   !!======================================================================
   !!                         ***  MODULE trcsms_fabm  ***
   !! TOP :   Main module of the FABM tracers
   !!======================================================================
   !! History :   1.0  !  2015-04  (PML) Original code
   !! History :   1.1  !  2020-06  (PML) Update to FABM 1.0, improved performance
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
   USE sbc_oce, only: lk_oasis,fr_i
   USE dom_oce
   USE zdf_oce
   USE iom
   USE xios
   USE cpl_oasis3
   USE st2D_fabm
   USE inputs_fabm
   USE vertical_movement_fabm

   !USE fldread         !  time interpolation

   IMPLICIT NONE

#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"

   PRIVATE

   PUBLIC   trc_sms_fabm            ! called by trcsms.F90 module
   PUBLIC   trc_sms_fabm_alloc      ! called by trcini_fabm.F90 module
   PUBLIC   trc_sms_fabm_check_mass ! called by trcwri_fabm.F90

   ! Work arrays
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: flux          ! Cross-interface flux of pelagic variables (# m-2 s-1)
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: current_total ! Totals of conserved quantities

   ! Arrays for environmental variables that are computed by the coupler
   REAL(wp), PUBLIC, TARGET, ALLOCATABLE, DIMENSION(:,:,:) :: prn,rho
   REAL(wp), PUBLIC, TARGET, ALLOCATABLE, DIMENSION(:,:) :: taubot
   REAL(wp), PUBLIC, TARGET :: daynumber_in_year

   ! State repair counters
   INTEGER, SAVE :: repair_interior_count = 0
   INTEGER, SAVE :: repair_surface_count  = 0
   INTEGER, SAVE :: repair_bottom_count   = 0

   ! Coupler parameters
   INTEGER, PUBLIC :: nn_adv  ! Vertical advection scheme for sinking/floating/movement
                              ! (1: 1st order upwind, 3: 3rd order TVD)

   ! Flag indicating whether model%start has been called (will be done on-demand)
   LOGICAL, SAVE :: started = .false.

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
      INTEGER :: jn, jk
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztrfabm, pdat
      REAL(wp), DIMENSION(jpi,jpj)    :: vint

!!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sms_fabm')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,'(a,i0,a,i4.4,a,i2.2,a,i2.2,a,i5,a)') &
          ' trc_sms_fabm:  FABM model, iteration ',kt,' ', &
          nyear,'-',nmonth,'-',nday,' ',nsec_day," secs"
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      IF (.NOT. started) CALL nemo_fabm_start

      CALL update_inputs( kt )

      CALL compute_fabm( kt )

      CALL compute_vertical_movement( kt, nn_adv )

      CALL st2d_fabm_nxt( kt )

      IF( l_trdtrc )  CALL wrk_alloc( jpi, jpj, jpk, ztrfabm )

      CALL trc_bc_read  ( kt )       ! tracers: surface and lateral Boundary Conditions
      CALL trc_rnf_fabm ( kt ) ! River forcings

      ! Send 3D diagnostics to output (these apply to time "n")
      DO jn = 1, size(model%interior_diagnostic_variables)
         IF (model%interior_diagnostic_variables(jn)%save) THEN
            ! Save 3D field
            pdat => model%get_interior_diagnostic_data(jn)
            CALL iom_put(model%interior_diagnostic_variables(jn)%name, pdat)

            ! Save depth integral if selected for output in XIOS
            IF (iom_use(TRIM(model%interior_diagnostic_variables(jn)%name)//'_VINT')) THEN
               vint = 0._wp
               DO jk = 1, jpkm1
                  vint = vint + pdat(:,:,jk) * fse3t(:,:,jk) * tmask(:,:,jk)
               END DO
               CALL iom_put(TRIM(model%interior_diagnostic_variables(jn)%name)//'_VINT', vint)
            END IF
         END IF
      END DO

      ! Send 2D diagnostics to output (these apply to time "n")
      DO jn = 1, size(model%horizontal_diagnostic_variables)
         IF (model%horizontal_diagnostic_variables(jn)%save) &
             CALL iom_put( model%horizontal_diagnostic_variables(jn)%name, model%get_horizontal_diagnostic_data(jn))
      END DO

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

   SUBROUTINE compute_fabm( kt )
      INTEGER, INTENT(in) :: kt   ! ocean time-step index

      INTEGER :: ji,jj,jk,jn
      LOGICAL :: valid, repaired
      REAL(wp) :: zalfg,zztmpx,zztmpy

      ! Validate current model state (setting argument to .TRUE. enables repair=clipping)
      CALL check_state(.TRUE., valid, repaired)
      IF (.NOT. valid) THEN
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
      IF (repaired) THEN
         WRITE(numout,*) "Total interior repairs up to now on process",narea,":",repair_interior_count
         WRITE(numout,*) "Total surface repairs up to now on process",narea,":",repair_surface_count
         WRITE(numout,*) "Total bottom repairs up to now on process",narea,":",repair_bottom_count
      ENDIF

      daynumber_in_year = fjulday - fjulstartyear + 1

      ! Compute the now hydrostatic pressure
      ! copied from istate.F90
      ! ------------------------------------

      IF (ALLOCATED(rho)) rho = rau0 * ( 1._wp + rhd )

      IF (ALLOCATED(prn)) THEN
         zalfg = 0.5e-4_wp * grav ! FABM wants dbar, convert from Pa (and multiply with 0.5 to average 2 cell thicknesses below)
         prn(:,:,1) = 10.1325_wp + zalfg * fse3t(:,:,1) * rho(:,:,1)
         DO jk = 2, jpkm1                                              ! Vertical integration from the surface
            prn(:,:,jk) = prn(:,:,jk-1) + zalfg * ( &
                        fse3t(:,:,jk-1) * rho(:,:,jk-1)  &
                        + fse3t(:,:,jk) * rho(:,:,jk) )
         END DO
      END IF

      ! Compute the bottom stress
      ! copied from diawri.F90
      ! ------------------------------------

      IF (ALLOCATED(taubot)) THEN
         taubot(:,:) = 0._wp
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
                  zztmpx = (  bfrua(ji  ,jj) * un(ji  ,jj,mbku(ji  ,jj))  &
                        &  + bfrua(ji-1,jj) * un(ji-1,jj,mbku(ji-1,jj))  )
                  zztmpy = (  bfrva(ji,  jj) * vn(ji,jj  ,mbkv(ji,jj  ))  &
                        &  + bfrva(ji,jj-1) * vn(ji,jj-1,mbkv(ji,jj-1))  )
                  taubot(ji,jj) = 0.5_wp * rau0 * SQRT( zztmpx * zztmpx + zztmpy * zztmpy ) * tmask(ji,jj,1)
                  !
            END DO
         END DO
      END IF

      CALL model%prepare_inputs(real(kt, wp),nyear,nmonth,nday,REAL(nsec_day,wp))

      ! TODO: retrieve 3D shortwave and store in etot3

      ! Zero rate array of interface-attached state variables
      fabm_st2Da = 0._wp

      ! Compute interfacial source terms and fluxes
      DO jj=2,jpjm1
         ! Process bottom (get_bottom_sources increments rather than sets, so zero flux array first)
         flux = 0._wp
         CALL model%get_bottom_sources(fs_2,fs_jpim1,jj,flux,fabm_st2Da(fs_2:fs_jpim1,jj,jp_fabm_surface+1:))
         DO jn=1,jp_fabm
            ! Divide bottom fluxes by height of bottom layer and add to source terms.
            DO ji=fs_2,fs_jpim1
               tra(ji,jj,mbkt(ji,jj),jp_fabm_m1+jn) = tra(ji,jj,mbkt(ji,jj),jp_fabm_m1+jn) + flux(ji,jn)/fse3t(ji,jj,mbkt(ji,jj))
            END DO
         END DO

         ! Process surface (get_surface_sources increments rather than sets, so zero flux array first)
         flux = 0._wp
         CALL model%get_surface_sources(fs_2,fs_jpim1,jj,flux,fabm_st2Da(fs_2:fs_jpim1,jj,1:jp_fabm_surface))
         ! Divide surface fluxes by height of surface layer and add to source terms.
         DO jn=1,jp_fabm
            DO ji=fs_2,fs_jpim1
               tra(ji,jj,1,jp_fabm_m1+jn) = tra(ji,jj,1,jp_fabm_m1+jn) + flux(ji,jn)/fse3t(ji,jj,1)
            END DO
         END DO
      END DO

      ! Compute interior source terms (NB get_interior_sources increments rather than sets)
      DO jk=1,jpkm1
         DO jj=2,jpjm1
            CALL model%get_interior_sources(fs_2,fs_jpim1,jj,jk,tra(fs_2:fs_jpim1,jj,jk,jp_fabm0:jp_fabm1))
         END DO
      END DO

      CALL model%finalize_outputs()
   END SUBROUTINE compute_fabm

   SUBROUTINE check_state(repair, valid, repaired)
      LOGICAL, INTENT(IN)  :: repair
      LOGICAL, INTENT(OUT) :: valid, repaired

      INTEGER :: jj, jk
      LOGICAL :: valid_int, valid_sf, valid_bt

      valid = .TRUE.     ! Whether the model state is valid after this subroutine returns
      repaired = .FALSE. ! Whether the model state has been repaired by this subroutine
      DO jk=1,jpkm1
         DO jj=2,jpjm1
            CALL model%check_interior_state(fs_2, fs_jpim1, jj, jk, repair, valid_int)
            IF (repair .AND. .NOT. valid_int) THEN
               repair_interior_count = repair_interior_count + 1
               repaired = .TRUE.
            END IF
            IF (.NOT. (valid_int .OR. repair)) valid = .FALSE.
         END DO
      END DO
      DO jj=2,jpjm1
         CALL model%check_surface_state(fs_2, fs_jpim1, jj, repair, valid_sf)
         IF (repair .AND. .NOT. valid_sf) THEN
            repair_surface_count = repair_surface_count + 1
            repaired = .TRUE.
         END IF
         IF (.NOT. (valid_sf .AND. valid_bt) .AND. .NOT. repair) valid = .FALSE.
         CALL model%check_bottom_state(fs_2, fs_jpim1, jj, repair, valid_bt)
         IF (repair .AND. .NOT. valid_bt) THEN
            repair_bottom_count = repair_bottom_count + 1
            repaired = .TRUE.
         END IF
         IF (.NOT. (valid_sf .AND. valid_bt) .AND. .NOT. repair) valid = .FALSE.
      END DO
   END SUBROUTINE

   SUBROUTINE trc_sms_fabm_check_mass()
      REAL(wp) :: total(SIZE(model%conserved_quantities))
      INTEGER :: ji,jk,jj,jn

      total = 0._wp

      IF (.NOT. started) CALL nemo_fabm_start

      DO jk=1,jpkm1
         DO jj=2,jpjm1
            CALL model%get_interior_conserved_quantities(fs_2,fs_jpim1,jj,jk,current_total)
            DO jn=1,SIZE(model%conserved_quantities)
               DO ji=fs_2,fs_jpim1
                  total(jn) = total(jn) + cvol(ji,jj,jk) * current_total(ji,jn) * tmask_i(ji,jj)
               END DO
            END DO
         END DO
      END DO

      DO jj=2,jpjm1
         CALL model%get_horizontal_conserved_quantities(fs_2,fs_jpim1,jj,current_total)
         DO jn=1,SIZE(model%conserved_quantities)
            DO ji=fs_2,fs_jpim1
               total(jn) = total(jn) + e1e2t(ji,jj) * current_total(ji,jn) * tmask_i(ji,jj)
            END DO
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
      INTEGER :: jn
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_fabm_alloc  ***
      !!----------------------------------------------------------------------
      !
      ! ALLOCATE here the arrays specific to FABM
      ALLOCATE( lk_rad_fabm(jp_fabm))
      IF (model%variable_needs_values(fabm_standard_variables%pressure)) ALLOCATE(prn(jpi, jpj, jpk))
      IF (ALLOCATED(prn) .or. model%variable_needs_values(fabm_standard_variables%density)) ALLOCATE(rho(jpi, jpj, jpk))
      IF (model%variable_needs_values(fabm_standard_variables%bottom_stress)) ALLOCATE(taubot(jpi, jpj))
      ! ALLOCATE( tab(...) , STAT=trc_sms_fabm_alloc )

      ! Allocate arrays to hold state for surface-attached and bottom-attached state variables
      ALLOCATE(fabm_st2Dn(jpi, jpj, jp_fabm_surface+jp_fabm_bottom))
      ALLOCATE(fabm_st2Da(jpi, jpj, jp_fabm_surface+jp_fabm_bottom))
      ALLOCATE(fabm_st2Db(jpi, jpj, jp_fabm_surface+jp_fabm_bottom))

      ! Work array to hold surface and bottom fluxes
      ALLOCATE(flux(fs_2:fs_jpim1,jp_fabm))

      ! Allocate work arrays for vertical movement
      ALLOCATE(w_ct(fs_2:fs_jpim1,1:jpkm1,jp_fabm))
      ALLOCATE(current_total(fs_2:fs_jpim1,SIZE(model%conserved_quantities)))
#if defined key_trdtrc && defined key_iomput
      IF( lk_trdtrc ) ALLOCATE(tr_vmv(jpi,jpj,jpk,jp_fabm))
      IF( lk_trdtrc ) ALLOCATE(tr_inp(jpi,jpj,jpk))
#endif

      trc_sms_fabm_alloc = 0      ! set to zero if no array to be allocated
      !
      IF( trc_sms_fabm_alloc /= 0 ) CALL ctl_warn('trc_sms_fabm_alloc : failed to allocate arrays')
      !

      ! Provide FABM with domain extents
      CALL model%set_domain(jpi, jpj, jpk, rdt)
      CALL model%set_domain_start(fs_2, 2, 1)
      CALL model%set_domain_stop(fs_jpim1, jpjm1, jpkm1)

      ! Provide FABM with the vertical indices of the bottom, and the land-sea mask.
      CALL model%set_bottom_index(mbkt)  ! NB mbkt extents should match dimension lengths provided to set_domain
      CALL model%set_mask(tmask,tmask(:,:,1)) ! NB tmask extents should match dimension lengths provided to set_domain

      ! Initialize state and send pointers to state data to FABM
      ! We mask land points in states with zeros, as per with NEMO "convention"
      ! NB we cannot call model%initialize_*_state at this point, because model%start has not been called yet.
      DO jn=1,jp_fabm
         trn(:,:,:,jp_fabm_m1+jn) = model%interior_state_variables(jn)%initial_value * tmask
         CALL model%link_interior_state_data(jn,trn(:,:,:,jp_fabm_m1+jn))
      END DO
      DO jn=1,jp_fabm_surface
         fabm_st2Dn(:,:,jn) = model%surface_state_variables(jn)%initial_value * tmask(:,:,1)
         CALL model%link_surface_state_data(jn,fabm_st2Dn(:,:,jn))
      END DO
      DO jn=1,jp_fabm_bottom
         fabm_st2Dn(:,:,jp_fabm_surface+jn) = model%bottom_state_variables(jn)%initial_value * tmask(:,:,1)
         CALL model%link_bottom_state_data(jn,fabm_st2Dn(:,:,jp_fabm_surface+jn))
      END DO

      ! Send pointers to environmental data to FABM
      CALL model%link_interior_data(fabm_standard_variables%depth, fsdept(:,:,:))
      CALL model%link_interior_data(fabm_standard_variables%temperature, tsn(:,:,:,jp_tem))
      CALL model%link_interior_data(fabm_standard_variables%practical_salinity, tsn(:,:,:,jp_sal))
      IF (ALLOCATED(rho)) CALL model%link_interior_data(fabm_standard_variables%density, rho(:,:,:))
      IF (ALLOCATED(prn)) CALL model%link_interior_data(fabm_standard_variables%pressure, prn)
      IF (ALLOCATED(taubot)) CALL model%link_horizontal_data(fabm_standard_variables%bottom_stress, taubot(:,:))
      CALL model%link_interior_data(fabm_standard_variables%cell_thickness, fse3t(:,:,:))
      CALL model%link_horizontal_data(fabm_standard_variables%latitude, gphit)
      CALL model%link_horizontal_data(fabm_standard_variables%longitude, glamt)
      CALL model%link_scalar(fabm_standard_variables%number_of_days_since_start_of_the_year, daynumber_in_year)
      CALL model%link_horizontal_data(fabm_standard_variables%wind_speed, wndm(:,:))
      CALL model%link_horizontal_data(fabm_standard_variables%surface_downwelling_shortwave_flux, qsr(:,:))
      CALL model%link_horizontal_data(fabm_standard_variables%bottom_depth_below_geoid, bathy(:,:))
      CALL model%link_horizontal_data(fabm_standard_variables%ice_area_fraction, fr_i(:,:))

      ! Obtain user-specified input variables (read from NetCDF file)
      CALL link_inputs
      CALL update_inputs(nit000, .FALSE.)

      ! Set mask for negativity corrections to the relevant states
      lk_rad_fabm(:) = .FALSE.
      DO jn=1,jp_fabm
        IF (model%interior_state_variables(jn)%minimum >= 0._wp) THEN
          lk_rad_fabm(jn) = .TRUE.
          IF(lwp) WRITE(numout,*) 'FABM clipping for '//TRIM(model%interior_state_variables(jn)%name)//' activated.'
        END IF
      END DO


      ! Copy initial condition for interface-attached state variables to "previous" state field
      ! NB NEMO does this itself for pelagic state variables (trb) in TOP_SRC/trcini.F90.
      fabm_st2Db = fabm_st2Dn

   END FUNCTION trc_sms_fabm_alloc

   SUBROUTINE nemo_fabm_start()
      INTEGER :: jn

      ! Make FABM aware of diagnostics that are not needed [not included in output]
      ! This works only after iom has completely initialised, because it depends on iom_use
      DO jn=1,size(model%interior_diagnostic_variables)
         model%interior_diagnostic_variables(jn)%save = iom_use(model%interior_diagnostic_variables(jn)%name) &
            .or. iom_use(TRIM(model%interior_diagnostic_variables(jn)%name)//'_VINT')
      END DO
      DO jn=1,size(model%horizontal_diagnostic_variables)
         model%horizontal_diagnostic_variables(jn)%save = iom_use(model%horizontal_diagnostic_variables(jn)%name)
      END DO

      ! Check whether FABM has all required data
      ! [after this, the save attribute of diagnostic variables can no longer change!]
      CALL model%start()

      started = .TRUE.
   END SUBROUTINE

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
