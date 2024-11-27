MODULE vertical_movement_fabm
   !!======================================================================
   !!                         ***  MODULE vertical_movement_fabm  ***
   !! TOP :   Module for the vertical movement of the FABM tracers
   !!======================================================================

#if defined key_fabm
   !!----------------------------------------------------------------------
   !!   'key_fabm'                                               FABM tracers
   !!----------------------------------------------------------------------
   !! compute_vertical_movement : compute vertical movement of FABM fields
   !!----------------------------------------------------------------------
   USE par_trc
   USE oce_trc
   USE trc
   USE fabm
   USE par_fabm
   USE dom_oce
#if defined key_trdtrc && defined key_iomput
   USE iom
   USE trdtrc_oce
#endif

   IMPLICIT NONE
   INTEGER, PUBLIC :: nn_sink_lbc    !: Type of boundary conditons for sinking ( ln_sink_slg )

#  include "vectopt_loop_substitute.h90"

   PRIVATE

   PUBLIC compute_vertical_movement

   ! Work arrays for vertical advection (residual movement/sinking/floating)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:,:) :: w_ct
#if defined key_trdtrc && defined key_iomput
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:,:,:) :: tr_vmv
#endif

   CONTAINS

   SUBROUTINE compute_vertical_movement( kt, method )
      !!----------------------------------------------------------------------
      !!                     ***  compute_vertical_movement  ***
      !!
      !! ** Purpose : compute vertical movement of FABM tracers through the water
      !!              (sinking/floating/active movement)
      !!
      !! ** Method  : Retrieves additional vertical velocity field and applies
      !!              advection scheme.
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   method ! advection method (1: 1st order upstream, 3: 3rd order TVD with QUICKEST limiter)

      INTEGER :: ji,jj,jk,jn,k_floor,k
      REAL(wp) :: zwgt_if(1:jpkm1-1), dc(1:jpkm1), w_if(1:jpkm1-1), z2dt, h(1:jpkm1), hb(1:jpkm1)
#if defined key_trdtrc
      CHARACTER (len=20) :: cltra
#endif

#if defined key_trdtrc && defined key_iomput
      IF( lk_trdtrc ) tr_vmv = 0.0_wp
#endif

      IF( neuler == 0 .AND. kt == nittrc000 ) THEN
          z2dt = rdt                  ! set time step size (Euler)
      ELSE
          z2dt = 2._wp * rdt          ! set time step size (Leapfrog)
      ENDIF
      ! Compute interior vertical velocities and include them in source array.
      DO jj=2,jpjm1 ! j-loop
         ! Get vertical velocities at layer centres (entire i-k slice).
         DO jk=1,jpkm1
            CALL model%get_vertical_movement(fs_2,fs_jpim1,jj,jk,w_ct(:,jk,:))
         END DO
         DO ji=fs_2,fs_jpim1 ! i-loop
            ! Only process this horizontal point (ji,jj) if number of layers exceeds 1
            k_floor = mbkt(ji,jj)
            IF (k_floor > 1) THEN ! Level check
               ! Linearly interpolate to velocities at the interfaces between layers
               ! Note:
               !    - interface k sits between cell centre k and k+1 (k=0 for surface)
               !    - k [1,jpkm1] increases downwards
               !    - upward velocity is positive, downward velocity is negative
               h(1:k_floor) = e3t_n(ji,jj,1:k_floor)
               hb(1:k_floor) = e3t_b(ji,jj,1:k_floor)
               zwgt_if(1:k_floor-1) = h(2:k_floor) / (h(1:k_floor-1) + h(2:k_floor))

               ! Advect:
               DO jn=1,jp_fabm ! State loop
                  IF (ALL(w_ct(ji,1:k_floor,jn) == 0._wp)) CYCLE

                  ! Compute velocities at interfaces
                  w_if(1:k_floor-1) = zwgt_if(1:k_floor-1) * w_ct(ji,1:k_floor-1,jn) + (1._wp - zwgt_if(1:k_floor-1)) * w_ct(ji,2:k_floor,jn)
                  

                WRITE(numout,*) 'Vertical movement computed using method = ',method
                  ! Compute change (per volume) due to vertical movement per layer
                  IF (method == 1) THEN
                     CALL advect_1(k_floor, trn(ji,jj,1:k_floor,jp_fabm_m1+jn), w_if(1:k_floor-1), h(1:k_floor), z2dt, dc(1:k_floor))
                 ELSE IF (method == 2) THEN
                    CALL semi_lagrangian_sedimentation(k_floor, tra(ji,jj,1:k_floor,jp_fabm_m1+jn), w_if(1:k_floor-1), h(1:k_floor), z2dt, gdepw_n(ji,jj,1:k_floor), tmask(ji,jj,k_floor), dc(1:k_floor))
                 ELSE IF (method == 3) THEN
                     CALL advect_3(k_floor, trb(ji,jj,1:k_floor,jp_fabm_m1+jn), w_if(1:k_floor-1), h(1:k_floor), z2dt, dc(1:k_floor))
                  END IF

                  ! Incorporate change due to vertical movement in sources-sinks
                  tra(ji,jj,1:k_floor,jp_fabm_m1+jn) = tra(ji,jj,1:k_floor,jp_fabm_m1+jn) + dc(1:k_floor)
                  !print tra(ji,jj,k_floor,jp_fabm_m1+jn)

#if defined key_trdtrc && defined key_iomput
                  ! Store change due to vertical movement as diagnostic
                  IF( lk_trdtrc .AND. ln_trdtrc( jp_fabm_m1+jn)) tr_vmv(ji,jj,1:k_floor,jn) = dc(1:k_floor)
#endif
              END DO ! State loop
            END IF ! Level check
         END DO ! i-loop
      END DO ! j-loop
#if defined key_trdtrc && defined key_iomput
      DO jn=1,jp_fabm ! State loop
        IF( lk_trdtrc .AND. ln_trdtrc(jp_fabm_m1+jn) ) THEN
          cltra = 'VMV_'//TRIM(ctrcnm(jp_fabm_m1+jn))
          CALL iom_put( cltra,  tr_vmv(:,:,:,jn) )
        END IF
      ENDDO
#endif

   END SUBROUTINE compute_vertical_movement

   SUBROUTINE semi_lagrangian_sedimentation(nk, c_old, w, h, dt, gdepw1, tmask1, trend)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nk                   ! Number of vertical levels
      REAL(wp), INTENT(IN) :: c_old(1:nk)         ! Old concentration
      REAL(wp), INTENT(IN) :: w(1:nk-1)           ! Sinking velocities (between layers)
      REAL(wp), INTENT(IN) :: h(1:nk)             ! Layer thicknesses
      REAL(wp), INTENT(IN) :: dt                  ! Time step
      REAL(wp), INTENT(IN) :: gdepw1(1:nk)        ! Depth of each grid point
      REAL(wp), INTENT(IN) :: tmask1(1:nk)        ! Mask indicating water presence
      REAL(wp), INTENT(OUT) :: trend(1:nk)        ! Trend/output flux due to sinking

      ! Local variables
      REAL(wp) :: zqR(nk), zqL(nk), zWR(nk), zWL(nk), zFC(nk+1)
      REAL(wp) :: zdltR, zdltL, zcff, zHz_inv2, zHz_inv3, zcu, zflx
      REAL(wp) :: zcffL, zcffR
      INTEGER :: jk, ik, ksource(nk)

      ! Initialize variables
      zqR = 0.0_wp
      zqL = 0.0_wp
      zWR = 0.0_wp
      zWL = 0.0_wp
      zFC = 0.0_wp
      ksource = 0

      ! Semi-Lagrangian flux computation
      DO jk = 2, nk
        zHz_inv2 = 1.0_wp / (h(jk) + h(jk-1))
        zFC(jk) = (c_old(jk-1) - c_old(jk)) * zHz_inv2
      END DO

      ! Apply PPM and WENO constraints
      DO jk = 2, nk-1
        zdltR = h(jk) * zFC(jk)
        zdltL = h(jk) * zFC(jk+1)
        zcff = h(jk+1) + 2.0_wp * h(jk) + h(jk-1)
        zcffR = zcff * zFC(jk)
        zcffL = zcff * zFC(jk+1)

        ! PPM monotonicity constraint
        IF (zdltR * zdltL <= 0.0_wp) THEN
          zdltR = 0.0_wp
          zdltL = 0.0_wp
        ELSE
          zdltR = MIN(ABS(zdltR), ABS(zcffL)) * SIGN(1.0_wp, zdltR)
          zdltL = MIN(ABS(zdltL), ABS(zcffR)) * SIGN(1.0_wp, zdltL)
        END IF

        ! Reconstruct right (zqR) and left (zqL) sides
        zHz_inv3 = 1.0_wp / (h(jk) + h(jk-1) + h(jk+1))
        zcff = (zdltR - zdltL) * zHz_inv3
        zdltR = zdltR - zcff * h(jk-1)
        zdltL = zdltL + zcff * h(jk+1)
        zqR(jk) = c_old(jk) + zdltR
        zqL(jk) = c_old(jk) - zdltL
        zWR(jk) = (2.0_wp * zdltR - zdltL) ** 2
        zWL(jk) = (zdltR - 2.0_wp * zdltL) ** 2
      END DO

      ! Reconciliation of parabolic profiles using WENO procedure
      DO jk = 2, nk-2
        zdltL = MAX(1.0e-14_wp, zWL(jk))
        zdltR = MAX(1.0e-14_wp, zWR(jk-1))
        zqR(jk) = (zdltR * zqR(jk) + zdltL * zqL(jk-1)) / (zdltR + zdltL)
        zqL(jk-1) = zqR(jk)
      END DO

      ! Boundary conditions
      zFC(1) = 0.0_wp
      zqL(1) = zqR(2)
      zqR(1) = 2.0_wp * c_old(1) - zqL(1)
      zqR(nk) = zqL(nk-1)
      zqL(nk) = 2.0_wp * c_old(nk) - zqR(nk)

      ! Reapply monotonicity constraint
      DO jk = 1, nk
        zdltR = zqR(jk) - c_old(jk)
        zdltL = c_old(jk) - zqL(jk)
        zcffR = 2.0_wp * zdltR
        zcffL = 2.0_wp * zdltL

        IF (zdltR * zdltL < 0.0_wp) THEN
          zdltR = 0.0_wp
          zdltL = 0.0_wp
        ELSEIF (ABS(zdltR) > ABS(zcffL)) THEN
          zdltR = zcffL
        ELSEIF (ABS(zdltL) > ABS(zcffR)) THEN
          zdltL = zcffR
        END IF

        zqR(jk) = c_old(jk) + zdltR
        zqL(jk) = c_old(jk) - zdltL
      END DO

      ! Compute the semi-Lagrangian advective flux
      DO jk = 1, nk-1
        zcff = dt * ABS(w(jk)) / rday * tmask1(jk)
        zFC(jk+1) = 0.0_wp
        zWL(jk) = -gdepw1(jk+1) + zcff
        zWR(jk) = h(jk) * c_old(jk)
        ksource(jk) = jk
      END DO

      DO jk = 1, nk
        DO ik = 2, jk
          IF (zWL(jk) > -gdepw1(ik)) THEN
            ksource(jk) = ik - 1
            zFC(jk+1) = zFC(jk+1) + zWR(ik)
          END IF
        END DO
      END DO

      ! Finalize flux computation
      DO jk = 1, nk-1
        ik = ksource(jk)
        zHz_inv2 = 1.0_wp / h(jk)
        zcu = MIN(1.0_wp, (zWL(jk) + gdepw1(ik+1)) * zHz_inv2)
        zFC(jk+1) = zFC(jk+1) + h(ik) * zcu * (zqL(ik) + zcu * (0.5_wp * (zqR(ik) - zqL(ik)) - &
                            (1.5_wp - zcu) * (zqR(ik) + zqL(ik) - 2.0_wp * c_old(ik))))
      END DO

      ! Update tracer concentration based on fluxes
      DO jk = 1, nk-1
        zHz_inv2 = 1.0_wp / h(jk)
        zflx = (zFC(jk) - zFC(jk+1)) * zHz_inv2
        trend(jk) = zflx
      END DO
      trend(nk) = -SUM(trend)

      trend(:) = trend(:) / dt


    END SUBROUTINE semi_lagrangian_sedimentation



   SUBROUTINE trc_sink2_slg(nk, c_old, w, h, dt, gdepw1, tmask1, trend)
        !!---------------------------------------------------------------------
        !!                     ***  ROUTINE trc_sink2_slg  ***
        !!
        !! ** Purpose :   Compute the sedimentation terms for the various sinking particles.
        !!                The scheme used to compute the trends is based on
        !!                a semi-Lagrangian advective flux algorithm
        !!
        !!---------------------------------------------------------------------

        INTEGER,  INTENT(IN)  :: nk               ! Number of vertical levels
        REAL(wp), INTENT(IN)  :: c_old(1:nk)      ! Old concentration
        REAL(wp), INTENT(IN)  :: w(1:nk-1)        ! Sinking velocities (between layers)
        REAL(wp), INTENT(IN)  :: h(1:nk)          ! Layer thicknesses
        REAL(wp), INTENT(IN)  :: dt               ! Time step
        REAL(wp), INTENT(IN) :: gdepw1(1:nk)      ! Trend/output flux due to sinking
        REAL(wp), INTENT(IN) :: tmask1(1:nk)      ! Trend/output flux due to sinking
        REAL(wp), INTENT(OUT) :: trend(1:nk)      ! Trend/output flux due to sinking
        !
        INTEGER  :: k
        REAL(wp) :: zcff, zcu, zdltL, zdltR, zflx, zcffR, zcffL, ik
        REAL(wp) :: zWR, zWL, zHz_inv, zHz_inv2, zHz_inv3
        REAL(wp) :: zFC(nk), zWL_arr(nk), zWR_arr(nk), zqR(nk), zqL(nk)
        INTEGER :: ksource(nk)
        !---------------------------------------------------------------------
        ! Initialize the trend and intermediate variables
        trend(:) = 0._wp
        zFC(:) = 0._wp
        zWL_arr(:) = 0._wp
        zWR_arr(:) = 0._wp
        ksource(:) = 0

        !-----------------------------------------------------------------------
        !  Compute semi-Lagrangian flux due to vertical sinking.
        !-----------------------------------------------------------------------

        ! Step 1: Compute flux difference (zFC) between vertical layers
        DO k = 2, nk-1
            zHz_inv2 = 1._wp / (h(k) + h(k-1))
            zFC(k) = (c_old(k-1) - c_old(k)) * zHz_inv2
        END DO

        ! Step 2: Apply parabolic reconstruction (PPM) to ensure monotonicity
        DO k = 2, nk-1
            zdltR = h(k) * zFC(k)
            zdltL = h(k) * zFC(k+1)
            zcff  = h(k+1) + 2. * h(k) + h(k-1)
            zcffR = zcff * zFC(k)
            zcffL = zcff * zFC(k+1)


            IF (zdltR * zdltL <= 0._wp) THEN
                zdltR = 0._wp
                zdltL = 0._wp
            ELSE IF (ABS(zdltR) >= zcffL) THEN
                zdltR = zdltL
            ELSE IF (ABS(zdltL) > ABS(zcffR)) THEN
                zdltL = zcffR
            END IF

            zHz_inv3 = 1._wp / (h(k) + h(k-1) + h(k+1))
            zcff = (zdltR - zdltL) * zHz_inv3
            zdltR = zdltR - zcff * h(k-1)
            zdltL = zdltL + zcff * h(k+1)

            zqR(k) = c_old(k) + zdltR
            zqL(k) = c_old(k) - zdltL
            zWR = (2._wp * zdltR - zdltL)**2
            zWL = (zdltR - 2._wp * zdltL)**2

            zWR_arr(k) = zWR
            zWL_arr(k) = zWL
        END DO


        ! Step 3: Reconcile parabolic segments for monotonicity
        zcff = 1.e-14
        DO k = 2, nk-2
            zdltL = MAX(zcff, zWL_arr(k))
            zdltR = MAX(zcff, zWR_arr(k-1))
            zqR(k) = (zdltR * zqR(k) + zdltL * zqL(k-1)) / (zdltR + zdltL)
            zqL(k-1) = zqR(k)
        END DO

        DO k = 1, nk
            zdltR = zqR(k) - c_old(k)
            zdltL = c_old(k) - zqL(k)
            zcffR = 2._wp * zdltR
            zcffL = 2._wp * zdltL
            IF( zdltR * zdltL < 0._wp ) THEN
               zdltR = 0._wp
               zdltL = 0._wp
            ELSE IF( ABS( zdltR ) > ABS( zcffL ) ) THEN
               zdltR = zcffL
            ELSE IF( ABS( zdltL ) > ABS( zcffR ) ) THEN
               zdltL = zcffR
            ENDIF
            zqR(k) = c_old(k)+ zdltR
            zqL(k) = c_old(k)- zdltL
        END DO 

        ! Step 4: Finalize flux computation using semi-Lagrangian approach
        DO k = 2, nk-1
            zcff = dt * ABS(w(k)) / rday * tmask1(k)
            zFC(k+1)   = 0._wp
            zWL_arr(k) = -gdepw1(k+1) + zcff
            zWR_arr(k) = h(k) * c_old(k)
            ksource(k) = k
        END DO
        DO k = 1, nk
            IF( zWL_arr(k) > -gdepw1(k) ) THEN
               ksource(k) = k - 1
               zFC(k+1) = zFC(k+1) + zWR_arr(k)
            ENDIF
        END DO

        ! Step 5: Compute final trend (sinking flux)
        DO k = 2, nk
            ik = ksource(k)
            zHz_inv = 1._wp / h(k)
            zcu = MIN(1._wp, (zWL_arr(k) + gdepw1(ik)) * zHz_inv)
            zFC(k) = zFC(k) + h(k) * zcu * (zqL(ik) + zcu * (0.5_wp * (zqR(ik) - zqL(ik)) - (1.5_wp - zcu) * (zqR(ik) + zqL(ik) - 2._wp * c_old(ik))))
        END DO

        ! Step 6: Compute the final trend and update the output array
        DO k = 1, nk
            zHz_inv = 1._wp / h(k)
            zflx = (zFC(k) - zFC(k+1)) * zHz_inv
            trend(k) = trend(k) + zflx
            write(numout,*),'trend(',k,') = ', trend(k)
        END DO
    END SUBROUTINE trc_sink2_slg

   SUBROUTINE advect_1(nk, c, w, h, dt, trend)
      INTEGER,  INTENT(IN)  :: nk
      REAL(wp), INTENT(IN)  :: c(1:nk)
      REAL(wp), INTENT(IN)  :: w(1:nk-1)
      REAL(wp), INTENT(IN)  :: h(1:nk)
      REAL(wp), INTENT(IN)  :: dt
      REAL(wp), INTENT(OUT) :: trend(1:nk)

      REAL(wp) :: flux(0:nk)
      INTEGER  :: jk
      ! Compute fluxes (per surface area) over at interfaces (remember: positive for upwards)
      flux(0) = 0._wp
      DO jk=1,nk-1 ! k-loop
         IF (w(jk) > 0) THEN
            ! Upward movement (source layer is jk+1)
            flux(jk) = min(w(jk), h(jk+1)/dt) * c(jk+1)
         ELSE
            ! Downward movement (source layer is jk)
            flux(jk) = max(w(jk), -h(jk)/dt) * c(jk)
         END IF
      END DO
      flux(nk) = 0._wp
      trend = (flux(1:nk) - flux(0:nk-1)) / h
   END SUBROUTINE

   SUBROUTINE advect_3(nk, c_old, w, h, dt, trend)
      INTEGER,  INTENT(IN)  :: nk
      REAL(wp), INTENT(IN)  :: c_old(1:nk)
      REAL(wp), INTENT(IN)  :: w(1:nk-1)
      REAL(wp), INTENT(IN)  :: h(1:nk)
      REAL(wp), INTENT(IN)  :: dt
      REAL(wp), INTENT(OUT) :: trend(1:nk)

      INTEGER, PARAMETER :: n_itermax=100
      REAL(wp) :: cmax_no
      REAL(wp) :: cfl(1:nk-1)
      INTEGER  :: n_iter, n_count, jk
      REAL(wp) :: c(1:nk)
      REAL(wp) :: tr_u(1:nk-1)
      REAL(wp) :: tr_c(1:nk-1)
      REAL(wp) :: tr_d(1:nk-1)
      REAL(wp) :: delta_tr_u(1:nk-1)
      REAL(wp) :: delta_tr(1:nk-1)
      REAL(wp) :: ratio(1:nk-1)
      REAL(wp) :: x_fac(1:nk-1)
      REAL(wp) :: phi_lim(1:nk-1)
      REAL(wp) :: limiter(1:nk-1)
      REAL(wp) :: flux_if(1:nk-1)

      c(:) = c_old(:)

      ! get maximum Courant number:
      cfl = ABS(w) * dt / (0.5_wp * (h(2:nk) + h(1:nk-1)))
      cmax_no = MAXVAL(cfl)

      ! number of iterations:
      n_iter = MIN(n_itermax, INT(cmax_no) + 1)
      IF (ln_ctl.AND.(n_iter .gt. 1)) THEN
         WRITE(numout,*) 'compute_vertical_movement::advect_3():'
         WRITE(numout,*) '   Maximum Courant number is ',cmax_no,'.'
         WRITE(numout,*) '   ',n_iter,' iterations used for vertical advection.'
      ENDIF

      ! effective Courant number:
      cfl = cfl/n_iter

      DO n_count=1,n_iter ! Iterative loop
         ! Determine tracer concentration at 1.5 upstream (tr_u), 0.5 upstream (tr_c), 0.5 downstream (tr_d) from interface
         IF (nk.gt.2) THEN
            ! More than 2 vertical wet points
            IF (nk.gt.3) THEN
               WHERE (w(2:nk-2).ge.0._wp)
                  !upward movement
                  tr_u(2:nk-2)=c(4:nk)
                  tr_c(2:nk-2)=c(3:nk-1)
                  tr_d(2:nk-2)=c(2:nk-2)
               ELSEWHERE
                  ! downward movement
                  tr_u(2:nk-2)=c(1:nk-3)
                  tr_c(2:nk-2)=c(2:nk-2)
                  tr_d(2:nk-2)=c(3:nk-1)
               ENDWHERE
            ENDIF

            ! Interface between surface layer and the next
            IF (w(1).ge.0._wp) THEN
               ! upward movement
               tr_u(1)=c(3)
               tr_c(1)=c(2)
               tr_d(1)=c(1)
            ELSE
               ! downward movement
               tr_u(1)=c(1)
               tr_c(1)=c(1)
               tr_d(1)=c(2)
            ENDIF

            ! Interface between bottom layer and the previous
            IF (w(nk-1).ge.0._wp) THEN
               ! upward movement
               tr_u(nk-1)=c(nk)
               tr_c(nk-1)=c(nk)
               tr_d(nk-1)=c(nk-1)
            ELSE
               ! downward movement
               tr_u(nk-1)=c(nk-2)
               tr_c(nk-1)=c(nk-1)
               tr_d(nk-1)=c(nk)
            ENDIF
         ELSE
            ! only 2 vertical wet points, i.e. only 1 interface
            IF (w(1).ge.0._wp) THEN
               ! upward movement
               tr_u(1)=c(2)
               tr_c(1)=c(2)
               tr_d(1)=c(1)
            ELSE
               ! downward movement
               tr_u(1)=c(1)
               tr_c(1)=c(1)
               tr_d(1)=c(2)
            ENDIF
         ENDIF

         delta_tr_u = tr_c - tr_u
         delta_tr = tr_d - tr_c
         WHERE (delta_tr * delta_tr_u > 0._wp)
            ! Monotonic function over tr_u, tr_c, r_d

            ! Compute slope ratio
            ratio = delta_tr_u / delta_tr

            ! QUICKEST flux limiter
            x_fac = (1._wp - 2._wp * cfl) / 6._wp
            phi_lim = (0.5_wp + x_fac) + (0.5_wp - x_fac) * ratio
            limiter = MIN(phi_lim, 2._wp / (1._wp - cfl), 2._wp * ratio / (cfl + 1.e-10_wp))

            ! Compute limited flux
            flux_if = w * (tr_c + 0.5_wp * limiter * (1._wp - cfl) * delta_tr)
         ELSEWHERE
            ! Non-monotonic, use 1st order upstream
            flux_if = w * tr_c
         ENDWHERE

         ! Compute pseudo update for trend aggregation:
         c(1:nk-1) = c(1:nk-1) + dt / real(n_iter, wp) / h(1:nk-1) * flux_if
         c(2:nk)   = c(2:nk)   - dt / real(n_iter, wp) / h(2:nk)   * flux_if

      ENDDO ! Iterative loop

      ! Estimate rate of change from pseudo state updates (source splitting):
      trend = (c - c_old) / dt
   END SUBROUTINE

   SUBROUTINE semi_lagrangian(nk, c_old, w, h, dt, trend)
        INTEGER,  INTENT(IN)  :: nk
        REAL(wp), INTENT(IN)  :: c_old(1:nk)
        REAL(wp), INTENT(IN)  :: w(1:nk-1)
        REAL(wp), INTENT(IN)  :: h(1:nk)
        REAL(wp), INTENT(IN)  :: dt
        REAL(wp), INTENT(OUT) :: trend(1:nk)

        REAL(wp) :: zcff, zcu, zcffL, zcffR, zdltL, zdltR
        REAL(wp) :: zHz_inv, zHz_inv2, zHz_inv3
        REAL(wp), DIMENSION(1:nk) :: zFC, zqR, zqL, zWR, zWL
        INTEGER, DIMENSION(1:nk) :: ksource
        INTEGER :: ik, jk

        ! Initialize variables
        zFC = 0._wp
        zqR = 0._wp
        zqL = 0._wp
        zWR = 0._wp
        zWL = 0._wp
         
        !-----------------------------------------------------------------------
        ! Vertical sinking of particle concentration
        !-----------------------------------------------------------------------

        ! Loop over the vertical levels
        DO jk = 2, nk-1
            ! Compute semi-Lagrangian flux due to sinking
            zHz_inv2 = 1._wp / (h(jk) + h(jk-1))
            zFC(jk) = (c_old(jk-1) - c_old(jk)) * zHz_inv2

            ! Apply PPM monotonicity constraint
            zdltR = h(jk) * zFC(jk)
            zdltL = h(jk) * zFC(jk+1)
            zcff = h(jk+1) + 2._wp * h(jk) + h(jk-1)
            zcffR = zcff * zFC(jk)
            zcffL = zcff * zFC(jk+1)
            !  Apply PPM monotonicity constraint to prevent oscillations within the grid box.
            IF (zdltR * zdltL <= 0._wp) THEN
                zdltR = 0._wp
                zdltL = 0._wp
            ELSE IF (ABS(zdltR) >= zcffL) THEN
                zdltR = zcffL
            ELSE IF (ABS(zdltL) > ABS(zcffR)) THEN
                zdltL = zcffR
            ENDIF
            !
            !  Compute right and left side values (qR,qL) of parabolic segments
            !  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
            !
            !  NOTE: Although each parabolic segment is monotonic within its grid
            !        box, monotonicity of the whole profile is not guaranteed,
            !        because qL(k+1)-qR(k) may still have different sign than
            !        tr(k+1)-tr(k).  This possibility is excluded, after qL and qR
            !        are reconciled using WENO procedure.
            !
            zHz_inv3 = 1._wp / (h(jk) + h(jk-1) + h(jk+1))
            zcff = (zdltR - zdltL) * zHz_inv3
            zdltR = zdltR - zcff * h(jk-1)
            zdltL = zdltL + zcff * h(jk+1)
            zqR(jk) = c_old(jk) + zdltR
            zqL(jk) = c_old(jk) - zdltL
            zWR(jk) = (2._wp * zdltR - zdltL)**2
            zWL(jk) = (zdltR - 2._wp * zdltL)**2
        END DO
        zHz_inv2 = 1._wp / (h(nk) + h(nk-1))
        zFC(nk) = (c_old(nk-1) - c_old(nk)) * zHz_inv2

        ! Boundary conditions
        zFC(1) = 0._wp
        zFC(nk) = 0._wp
        zqR(1) = c_old(1)
        zqL(1) = c_old(1)
        zqR(nk) = c_old(nk)
        zqL(nk) = c_old(nk)

        !  Apply monotonicity constraint again, since the reconciled interfacial
        !  values may cause a non-monotonic behavior of the parabolic segments
        !  inside the grid box.
        DO jk = 2, nk-1
            zdltR = zqR(jk) - c_old(jk)
            zdltL = c_old(jk) - zqL(jk)
            zcffR = 2._wp * zdltR
            zcffL = 2._wp * zdltL
            IF (zdltR * zdltL < 0._wp) THEN
                zdltR = 0._wp
                zdltL = 0._wp
            ELSE IF (ABS(zdltR) > ABS(zcffL)) THEN
                zdltR = zcffL
            ELSE IF (ABS(zdltL) > ABS(zcffR)) THEN
                zdltL = zcffR
            ENDIF
            zqR(jk) = c_old(jk) + zdltR
            zqL(jk) = c_old(jk) - zdltL
        END DO
        !  After this moment reconstruction is considered complete. The next
        !  stage is to compute vertical advective fluxes, FC. It is expected
        !  that sinking may occurs relatively fast, the algorithm is designed
        !  to be free of CFL criterion, which is achieved by allowing
        !  integration bounds for semi-Lagrangian advective flux to use as
        !  many grid boxes in upstream direction as necessary.

        !  In the two code segments below, WL is the z-coordinate of the
        !  departure point for grid box interface z_w with the same indices;
        !  FC is the finite volume flux; ksource(:,k) is index of vertical
        !  grid box which contains the departure point (restricted by N(ng)).
        !  During the search: also add in content of whole grid boxes
        !  participating in FC.
        DO jk = 2, nk-1
            zcff = dt * ABS(w(jk)) / h(jk)
            zFC(jk) = 0._wp
            zWL(jk) = -zcff
            zWR(jk) = h(jk) * c_old(jk)
            ksource(jk) = jk
        END DO

        DO jk = 1, nk
            DO ik = 2, jk
                IF (zWL(jk) > -zcff) THEN
                    ksource(jk) = ik - 1
                    zFC(jk+1) = zFC(jk+1) + zWR(ik)
                ENDIF
            END DO
        END DO

        ! Finalize computation of flux and calculate trend
        DO jk = 2, nk-1
            zHz_inv = 1._wp / h(jk)
            zcu = MIN(1._wp, zWL(jk) * zHz_inv)
            zFC(jk+1) = zFC(jk+1) + h(jk) * zcu * (zqL(jk) + zcu * (0.5_wp * (zqR(jk) - zqL(jk)) - (1.5_wp - zcu) * (zqR(jk) + zqL(jk) - 2._wp * c_old(jk))))
        END DO

        DO jk = 2, nk-1
            zHz_inv = 1._wp / h(jk)
            trend(jk) = (zFC(jk) - zFC(jk+1)) * zHz_inv
        END DO
    END SUBROUTINE semi_lagrangian


#endif
END MODULE
