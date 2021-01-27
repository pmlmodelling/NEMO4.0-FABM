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
   USE par_fabm
   USE fabm
   USE dom_oce
#if defined key_trdtrc && defined key_iomput
   USE iom
   USE trdtrc_oce
#endif

   IMPLICIT NONE

!jth#  include "domzgr_substitute.h90"

   PRIVATE

   PUBLIC compute_vertical_movement

   ! Work arrays for vertical advection (residual movement/sinking/floating)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:,:) :: w_ct
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)   :: w_if
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)   :: zwgt_if
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)   :: flux_if
#if defined key_trdtrc && defined key_iomput
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:,:,:) :: tr_vmv
#endif

   CONTAINS

   SUBROUTINE compute_vertical_movement( kt )
      !!----------------------------------------------------------------------
      !!                     ***  compute_vertical_movement  ***
      !!
      !! ** Purpose :   compute vertical movement of FABM tracers
      !!
      !! ** Method  : Sets additional vertical velocity field and computes
      !!              resulting advection using a conservative 3rd upwind
      !!              scheme with QUICKEST TVD limiter, based on GOTM
      !!              module adv_center.F90 (www.gotm.net). Currently assuming
      !!              zero flux at sea surface and sea floor.
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER :: ji,jj,jk,jn,k_floor,n_iter,n_count
      INTEGER,PARAMETER :: n_itermax=100
      REAL(wp) :: cmax_no,z2dt
      REAL(wp),DIMENSION(jpk) :: tr_it,tr_u,tr_d,tr_c,tr_slope,c_no,flux_lim
      REAL(wp),DIMENSION(jpk) :: phi_lim,x_fac
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
      DO jj=1,jpj ! j-loop
         ! Get vertical velocities at layer centres (entire 1:jpi,1:jpk slice).
         DO jk=1,jpk
            CALL fabm_get_vertical_movement(model,1,jpi,jj,jk,w_ct(:,jk,:))
         END DO

         DO ji=1,jpi ! i-loop
            ! Only process this horizontal point (ji,jj) if number of layers exceeds 1
            IF (mbkt(ji,jj)>1) THEN ! Level check
               k_floor=mbkt(ji,jj)
               ! Linearly interpolate to velocities at the interfaces between layers
               ! Note:
               !    - interface k sits between cell centre k and k-1,
               !    - k [1,jpk] increases downwards
               !    - upward velocity is positive, downward velocity is negative
               zwgt_if(1,:)=0._wp ! surface
               w_if(1,:)=0._wp ! surface
               zwgt_if(2:k_floor,:)=spread(&
                   e3t_n(ji,jj,2:k_floor)/ (e3t_n(ji,jj,1:k_floor-1)+e3t_n(ji,jj,2:k_floor))&
                   ,2,jp_fabm)
               w_if(2:k_floor,:) = zwgt_if(2:k_floor,:)*w_ct(ji,1:k_floor-1,:)&
                  +(1._wp-zwgt_if(1:k_floor-1,:))*w_ct(ji,2:k_floor,:)
               zwgt_if(k_floor+1:,:)=0._wp ! sea floor and below
               w_if(k_floor+1:,:)=0._wp ! sea floor and below

               ! Advect:
               DO jn=1,jp_fabm ! State loop
                  ! get maximum Courant number:
                  c_no(2:k_floor)=abs(w_if(2:k_floor,jn))*z2dt/ &
                                ( 0.5_wp*(e3t_n(ji,jj,2:k_floor) + &
                                e3t_n(ji,ji,1:k_floor-1)) )
                  cmax_no=MAXVAL(c_no(2:k_floor))

                  ! number of iterations:
                  n_iter=min(n_itermax,int(cmax_no)+1)
                  IF (ln_ctl.AND.(n_iter .gt. 1)) THEN
                      WRITE(numout,*) 'vertical_movement_fabm():'
                      WRITE(numout,*) '   Maximum Courant number is ',cmax_no,'.'
                      WRITE(numout,*) '   ',n_iter,' iterations used for vertical advection.'
                  ENDIF

                  ! effective Courant number:
                  c_no=c_no/n_iter

                  tr_it(1:k_floor)=trb(ji,jj,1:k_floor,jp_fabm_m1+jn)
                  DO n_count=1,n_iter ! Iterative loop
                     !Compute slope ratio
                     IF (k_floor.gt.2) THEN !More than 2 vertical wet points
                        IF (k_floor.gt.3) THEN
                          WHERE (w_if(3:k_floor-1,jn).ge.0._wp) !upward movement
                           tr_u(3:k_floor-1)=tr_it(4:k_floor)
                           tr_c(3:k_floor-1)=tr_it(3:k_floor-1)
                           tr_d(3:k_floor-1)=tr_it(2:k_floor-2)
                          ELSEWHERE !downward movement
                           tr_u(3:k_floor-1)=tr_it(1:k_floor-3)
                           tr_c(3:k_floor-1)=tr_it(2:k_floor-2)
                           tr_d(3:k_floor-1)=tr_it(3:k_floor-1)
                          ENDWHERE
                        ENDIF
                        IF (w_if(2,jn).ge.0._wp) THEN
                           tr_u(2)=tr_it(3)
                           tr_c(2)=tr_it(2)
                           tr_d(2)=tr_it(1)
                        ELSE
                           tr_u(2)=tr_it(1)
                           tr_c(2)=tr_it(1)
                           tr_d(2)=tr_it(2)
                        ENDIF
                        IF (w_if(k_floor,jn).ge.0._wp) THEN
                           tr_u(k_floor)=tr_it(k_floor)
                           tr_c(k_floor)=tr_it(k_floor)
                           tr_d(k_floor)=tr_it(k_floor-1)
                        ELSE
                           tr_u(k_floor)=tr_it(k_floor-2)
                           tr_c(k_floor)=tr_it(k_floor-1)
                           tr_d(k_floor)=tr_it(k_floor)
                        ENDIF
                     ELSE !only 2 vertical wet points, i.e. only 1 interface
                        IF (w_if(k_floor,jn).ge.0._wp) THEN
                           tr_u(2)=tr_it(2)
                           tr_c(2)=tr_it(2)
                           tr_d(2)=tr_it(1)
                        ELSE
                           tr_u(2)=tr_it(1)
                           tr_c(2)=tr_it(1)
                           tr_d(2)=tr_it(2)
                        ENDIF
                     ENDIF
                     WHERE (abs(tr_d(2:k_floor)-tr_c(2:k_floor)).gt.1.e-10_wp)
                        tr_slope(2:k_floor)= &
                           (tr_c(2:k_floor)-tr_u(2:k_floor))/ &
                           (tr_d(2:k_floor)-tr_c(2:k_floor))
                     ELSEWHERE
                        tr_slope(2:k_floor)=SIGN(1._wp,w_if(2:k_floor,jn))* &
                              (tr_c(2:k_floor)-tr_u(2:k_floor))*1.e10_wp
                     ENDWHERE

                     !QUICKEST flux limiter:
                     x_fac(2:k_floor)=(1._wp-2._wp*c_no(2:k_floor))/6._wp
                     phi_lim(2:k_floor)=(0.5_wp+x_fac(2:k_floor)) + &
                        (0.5_wp-x_Fac(2:k_floor))*tr_slope(2:k_floor)
                     flux_lim(2:k_floor)=max( 0._wp, &
                       min( phi_lim(2:k_floor),2._wp/(1._wp-c_no(2:k_floor)), &
                         2._wp*tr_slope(2:k_floor)/(c_no(2:k_floor)+1.e-10_wp)) )

                     ! Compute limited flux:
                     flux_if(2:k_floor,jn) = w_if(2:k_floor,jn)* &
                        ( tr_c(2:k_floor) + &
                        0.5_wp*flux_lim(2:k_floor)*(1._wp-c_no(2:k_floor))* &
                        (tr_d(2:k_floor)-tr_c(2:k_floor)) )

                     ! Compute pseudo update for trend aggregation:
                     tr_it(1:k_floor-1) = tr_it(1:k_floor-1) + &
                        z2dt/float(n_iter)/e3t_n(ji,jj,1:k_floor-1)* &
                        flux_if(2:k_floor,jn)
                     tr_it(2:k_floor) = tr_it(2:k_floor) - &
                        z2dt/float(n_iter)/e3t_n(ji,jj,2:k_floor)* &
                        flux_if(2:k_floor,jn)

                  ENDDO ! Iterative loop

                  ! Estimate rate of change from pseudo state updates (source
                  ! splitting):
                  tra(ji,jj,1:k_floor,jp_fabm_m1+jn) = &
                     tra(ji,jj,1:k_floor,jp_fabm_m1+jn) + &
                     (tr_it(1:k_floor) - trb(ji,jj,1:k_floor,jp_fabm_m1+jn))/z2dt
#if defined key_trdtrc && defined key_iomput
                  IF( lk_trdtrc .AND. ln_trdtrc( jp_fabm_m1+jn ) ) THEN
                    tr_vmv(ji,jj,1:k_floor,jn)=(tr_it(1:k_floor) - trb(ji,jj,1:k_floor,jn))/z2dt
                  END IF
#endif
               ENDDO ! State loop
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

#endif
END MODULE
