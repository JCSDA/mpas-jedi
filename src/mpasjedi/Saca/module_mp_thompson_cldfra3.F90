!=================================================================================================================
 module module_mp_thompson_cldfra3_saca

!module_mp_thompson_cldfra3 contains the subroutine cal_cldfra3 which calculates the cloud fraction as a function
!of relative humidity. The subroutine cal_cldfra3 was tested in WRF 3.8.1 with the Thompson cloud microphysics
!and should not be used with other cloud microphysics schemes.

!subroutine cal_cldfra3 was originally copied from ./phys/module_radiation_driver.F from WRF version 3.8.1.
!Laura D. Fowler (laura@ucar.edu) / 2016-09-22.

! add-ons and modifications to sourcecode:
! ----------------------------------------
! * in subroutine find_cloudLayers, changed the line k = k_m12C+2 to k = min(k_m12C,k_m12C+2) to avoid k greater
!   than the model-top index.
!   Laura D. Fowler (laura@ucar.edu)/2016-09-23.

 use mpas_kind_types, only : RKIND
 use mpas_atmphys_functions, only: rslf,rsif
 use mpas_constants, only : rgas, rv

 implicit none
 private
 public:: saca_param_type, cal_cldfra3

!---------------------------------------------------------------------
!> Fortran derived type to hold MPAS field
 type :: saca_param_type
   private

   logical, public :: l_build_madwrf
   logical, public :: l_build_gsdcloud
   logical, public :: l_saturate_qv
   logical, public :: l_conserve_thetaV
   REAL(KIND=RKIND), public :: cldfra_def
   REAL(KIND=RKIND), public :: cldfra_thresh
   REAL(KIND=RKIND), public :: cldmask_thresh
   REAL(KIND=RKIND), public :: cld_bld_hgt

 end type saca_param_type

!---------------------------------------------------------------------

 contains

!=================================================================================================================
!+---+-----------------------------------------------------------------+
!..Cloud fraction scheme by G. Thompson (NCAR-RAL), not intended for
!.. combining with any cumulus or shallow cumulus parameterization
!.. scheme cloud fractions.  This is intended as a stand-alone for
!.. cloud fraction and is relatively good at getting widespread stratus
!.. and stratoCu without caring whether any deep/shallow Cu param schemes
!.. is making sub-grid-spacing clouds/precip.  Under the hood, this
!.. scheme follows Mocko and Cotton (1995) in applicaiton of the
!.. Sundqvist et al (1989) scheme but using a grid-scale dependent
!.. RH threshold, one each for land v. ocean points based on
!.. experiences with HWRF testing.
!+---+-----------------------------------------------------------------+
!
!+---+-----------------------------------------------------------------+

      SUBROUTINE cal_cldfra3(CLDFRA, qv, qc, qi, qs, ni,                &
     &                 dz,p,t,rho, XLAND, gridkm,                       &
!    &                 rand_perturb_on, kme_stoch, rand_pert,           &
     &                 modify_qvapor,                                   &
     &                 ids,ide, jds,jde, kds,kde,                       &
     &                 ims,ime, jms,jme, kms,kme,                       &
     &                 its,ite, jts,jte, kts,kte,                       &
     &                 k_tropo, cldmask, cldtopz,                       &
     &                 saca_params, diag, diag_cldfra)
!
!     USE module_mp_thompson   , ONLY : rsif, rslf
!     IMPLICIT NONE
!
      INTEGER, INTENT(IN):: ids,ide, jds,jde, kds,kde,                  &
     &                      ims,ime, jms,jme, kms,kme,                  &
!    &                      kme_stoch,                                  &
     &                      its,ite, jts,jte, kts,kte

!     INTEGER, INTENT(IN):: rand_perturb_on
      LOGICAL, INTENT(IN):: modify_qvapor
      REAL(KIND=RKIND), DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN):: p,rho,dz
      REAL(KIND=RKIND), DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT):: t,qv,qc,qi,qs,ni
!     REAL(KIND=RKIND), DIMENSION(ims:ime,kms:kme_stoch,jms:jme), INTENT(IN):: rand_pert
      REAL(KIND=RKIND), DIMENSION(ims:ime,jms:jme), INTENT(IN):: XLAND

      REAL(KIND=RKIND), DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT):: cldfra
      REAL(KIND=RKIND), DIMENSION(ims:ime,jms:jme), INTENT(IN):: gridkm

      integer, dimension(ims:ime,jms:jme), intent(in), optional :: k_tropo
      REAL(KIND=RKIND), dimension(ims:ime,jms:jme), intent(in), optional :: cldmask, cldtopz
      type(saca_param_type), intent(in) :: saca_params
      REAL(KIND=RKIND), DIMENSION(ims:ime,kms:kme,jms:jme),             &
     &                  INTENT(INOUT), optional :: diag, diag_cldfra

!..Local vars.
      REAL(KIND=RKIND)::  RH_00, RHI_max, entrmnt
      REAL(KIND=RKIND), DIMENSION(its:ite,jts:jte):: RH_00L, RH_00O
      REAL(KIND=RKIND), DIMENSION(ims:ime,kms:kme,jms:jme):: qvsat, rh
      INTEGER:: i,j,k
      REAL(KIND=RKIND):: TK, TC, qvsi, qvsw, RHUM, xx, yy
      REAL(KIND=RKIND), DIMENSION(kms:kme):: qvs1d, cfr1d, T1d,                     &
     &                           P1d, R1d, qv1d, qc1d, qi1d, qs1d

      character*512 dbg_msg
      LOGICAL:: debug_flag
      logical :: is_tropo_init, impose_cldtopz
      logical, dimension(ims:ime,jms:jme) :: keep_clouds_below_lcl

      integer, dimension(ims:ime,jms:jme) :: cldtopk
      REAL(KIND=RKIND), PARAMETER :: CLDTHICK_DEF = -999.9  ! Default thickness [m] of new clouds
      REAL(KIND=RKIND)            :: CLDFRA_DEF             ! cloud fraction to insert in new clouds
      REAL(KIND=RKIND)            :: CLDFRA_THRESH          ! cloud fraction threshold to determine presence of clouds
      REAL(KIND=RKIND), PARAMETER :: RH_NOCLOUD = 0.3
      INTEGER :: cldthick_bot_k(ims:ime,jms:jme), cldfra_thresh_k(ims:ime,jms:jme)
!+---+
      REAL(KIND=RKIND)            :: cldmask_thresh         ! threshold to distinguish clear/cloudy for "observed" cldmask (or cloud fraction)
      REAL(KIND=RKIND), DIMENSION(ims:ime,kms:kme,jms:jme) :: cldfra_bk
      logical :: l_build_madwrf, l_build_gsdcloud
      logical :: l_saturate_qv, l_conserve_thetaV

!---- gsdcloud parameters
      REAL(KIND=RKIND), PARAMETER :: sat_cloud_pthick_p = 30.,  & !from cloudCover_NESDIS.f90
                                     cloud_q_qvis_rat_p = 0.05    !
      REAL(KIND=RKIND)            :: cld_bld_hgt
      REAL(KIND=RKIND), PARAMETER :: auto_conver = 0.0002,      & !from cloudLWC.f90
                                     temp_qvis1 = 268.15,       & !
                                     temp_qvis2 = 263.15          !
      REAL(KIND=RKIND), PARAMETER :: rh_clear_p = 0.8,         & !from cloud_saturation.f90
                                     qv_max_inc = 0.005,       & ! kg/kg
                                     miter = 3
      REAL(KIND=RKIND) :: temp, qtemp, cloudqvis, q_bk, constantTv
      integer :: nnn

      REAL(KIND=RKIND) :: watwgt, qavail, cloud_q_qvis_ratio
      integer :: k_cld_bld_hgt


!---SACA 0: set up parameters
!!! added for cloud init !!!
      is_tropo_init = .true.

      keep_clouds_below_lcl = .false.

      impose_cldtopz = .true.

      cldtopk = -999
      cldthick_bot_k = -999
      cldfra_thresh_k = -999

      l_build_gsdcloud  = saca_params % l_build_gsdcloud
      l_build_madwrf    = saca_params % l_build_madwrf
      l_saturate_qv     = saca_params % l_saturate_qv
      l_conserve_thetaV = saca_params % l_conserve_thetaV
      CLDFRA_DEF        = saca_params % cldfra_def
      CLDFRA_THRESH     = saca_params % cldfra_thresh
      cldmask_thresh    = saca_params % cldmask_thresh
      cld_bld_hgt       = saca_params % cld_bld_hgt
      write(0,*) " l_build_gsdcloud = ",l_build_gsdcloud
      write(0,*) "   l_build_madwrf = ",l_build_madwrf
      write(0,*) "    l_saturate_qv = ",l_saturate_qv
      write(0,*) "l_conserve_thetaV = ",l_conserve_thetaV
      write(0,*) "       cldfra_def = ",CLDFRA_DEF
      write(0,*) "    cldfra_thresh = ",CLDFRA_THRESH
      write(0,*) "   cldmask_thresh = ",cldmask_thresh
      write(0,*) "      cld_bld_hgt = ",cld_bld_hgt

      if(present(diag)) diag(:,1,:) = real( k_tropo(:,:) ) !BJJ, diag, k-index for tropopause (bg)


!---SACA 1: Calculate the model cloud fraction
!..First cut scale-aware. Higher resolution should require closer to
!.. saturated grid box for higher cloud fraction.  Simple functions
!.. chosen based on Mocko and Cotton (1995) starting point and desire
!.. to get near 100% RH as grid spacing moves toward 1.0km, but higher
!.. RH over ocean required as compared to over land.

      do j = jts,jte
      do i = its,ite
         RH_00L(i,j) = 0.781 + SQRT(1./(35.0+gridkm(i,j)*gridkm(i,j)*gridkm(i,j)*0.5))
         RH_00O(i,j) = 0.831 + SQRT(1./(70.0+gridkm(i,j)*gridkm(i,j)*gridkm(i,j)*0.5))
      enddo
      enddo

      DO j = jts,jte
      DO k = kts,kte
      DO i = its,ite

         CLDFRA(I,K,J) = 0.0

         if (qc(i,k,j).gt.1.E-6 .or. qi(i,k,j).ge.1.E-7 .or. qs(i,k,j).gt.1.E-5) then
            CLDFRA(I,K,J) = 1.0
            qvsat(i,k,j) = qv(i,k,j)
         else
            TK   = t(i,k,j)
            TC   = TK - 273.16

            qvsw = rslf(P(i,k,j), TK)
            qvsi = rsif(P(i,k,j), TK)

            if (tc .ge. -12.0) then
               qvsat(i,k,j) = qvsw
            elseif (tc .lt. -20.0) then
               qvsat(i,k,j) = qvsi
            else
               qvsat(i,k,j) = qvsw - (qvsw-qvsi)*(-12.0-tc)/(-12.0+20.)
            endif

            if (modify_qvapor) then
               if (qc(i,k,j).gt.1.E-8) then
                  qv(i,k,j) = MAX(qv(i,k,j), qvsw)
                  qvsat(i,k,j) = qvsw
               endif
               if (qc(i,k,j).le.1.E-8 .and. qi(i,k,j).ge.1.E-9) then
                  qv(i,k,j) = MAX(qv(i,k,j), qvsi)
                  qvsat(i,k,j) = qvsi
               endif
            endif

            RHUM = MAX(0.01, MIN(qv(i,k,j)/qvsat(i,k,j), 0.9999))
            rh(i,k,j) = MAX(0.01, qv(i,k,j)/qvsat(i,k,j))

            IF ((XLAND(I,J)-1.5).GT.0.) THEN                             !--- Ocean
               RH_00 = RH_00O(i,j)
            ELSE                                                         !--- Land
               RH_00 = RH_00L(i,j)
            ENDIF

            if (tc .ge. -12.0) then
               RHUM = MIN(0.999, RHUM)
               CLDFRA(I,K,J) = MAX(0.0, 1.0-SQRT((1.0-RHUM)/(1.-RH_00)))
            elseif (tc.lt.-12..and.tc.gt.-70. .and. RHUM.gt.RH_00O(i,j)) then
               RHUM = MAX(0.01, MIN(qv(i,k,j)/qvsat(i,k,j), qvsw/qvsi - 1.E-6))
               RHI_max = MAX(RHUM+1.E-6, qvsw/qvsi)
               CLDFRA(I,K,J) = MAX(0., 1.0-SQRT((RHI_max-RHUM)/(RHI_max-RH_00O(i,j))))
            endif
            CLDFRA(I,K,J) = MIN(0.90, CLDFRA(I,K,J))

         endif
      ENDDO
      ENDDO
      ENDDO

      cldfra_bk = CLDFRA !-BJJ
      if(present(diag_cldfra)) diag_cldfra = CLDFRA

!---SACA 2: Cloud clearing
!---CF/QX Remove---
      ! Remove hydrometeors if necessary and
      ! calculate height related variables

      if (impose_cldtopz) then
         do j = jts,jte
         do i = its,ite
            if (cldmask(i,j) <= cldmask_thresh) then
               do k = kts, kte
                  cldfra(i,k,j) = 0.0
                  qc(i,k,j) = 0.0
                  qi(i,k,j) = 0.0
                  qs(i,k,j) = 0.0
               end do
               if(present(diag)) diag(i,11,j) = 1.0 !BJJ, diag, remove entire cloud profile at clear observed cloud mask

            else if (cldmask(i,j) > cldmask_thresh .and. cldmask(i,j) <= 1.0 .and. cldtopz(i,j) > 0.0) then
               ! Zero out hydrometeors above cldtopz
               call find_k_from_z_agl(cldtopz(i,j), cldtopk(i,j), dz(i,kts:kte,j), kts, kte)
               if(present(diag)) diag(:,2,:) = real( cldtopk(:,:) ) !BJJ, diag, k-index for observed cloud top (ABI)
               if (cldtopk(i,j) > 0 .and. cldtopk(i,j) < kte) then
                  do k = cldtopk(i,j) + 1, kte
                     cldfra(i,k,j) = 0.0
                     qc(i,k,j) = 0.0
                     qi(i,k,j) = 0.0
                     qs(i,k,j) = 0.0
                  end do
               end if
               if(present(diag)) diag(i,12,j) = 1.0 !BJJ, diag, remove cloud above obs_cloud_top

               ! Calculate the cloud bottom k if the specified thickness were to be applied
               ! If that would be below ground, then set the k level to 1
               if (CLDTHICK_DEF > 0.0) then
                  if (cldtopz(i,j) - CLDTHICK_DEF > 0.0) then
                     call find_k_from_z_agl(cldtopz(i,j) - CLDTHICK_DEF, cldthick_bot_k(i,j), dz(i,kts:kte,j), kts, kte)
                  else
                     cldthick_bot_k(i,j) = 1
                  end if
               end if
            end if
         end do
         end do
      end if


!---SACA 3a: Cloud building w/ gsdcloud algorithm
!---GSDCLOUD CF/QX Add---
      if( l_build_gsdcloud ) then
         do j = jts,jte
         do i = its,ite
            if ( cldtopk(i,j) > 0 .and. cldtopk(i,j) < kte) then

               call find_k_from_z_agl( cld_bld_hgt, k_cld_bld_hgt, dz(i,kts:kte,j), kts, kte)

               if ( cldtopk(i,j) .le. k_cld_bld_hgt ) then

                  do k=cldtopk(i,j), 1, -1

                     if ( p(i,k,j) - p(i, cldtopk(i,j), j) .le. sat_cloud_pthick_p ) then

                        cldfra(i, k ,j) = max( cldfra(i, k ,j), CLDFRA_DEF )
                        if(present(diag)) diag(i,21,j) = 1.0 !BJJ, diag, build cloud w/ gsdcloud

                        temp = t(i,k,j)
                        if ( temp >= temp_qvis2 ) then
                           watwgt = 1.
                           cloud_q_qvis_ratio = watwgt*cloud_q_qvis_rat_p
                           qavail = min(0.25*auto_conver, cloud_q_qvis_ratio*qvsat(i,k,j))
                        else
                           watwgt = 0.
                           cloud_q_qvis_ratio = 0.1*cloud_q_qvis_rat_p
                           qavail = min(0.1*auto_conver, cloud_q_qvis_ratio*qvsat(i,k,j))
                        end if

                        qc(i,k,j) = max( qc(i,k,j), watwgt*qavail )
                        qi(i,k,j) = max( qi(i,k,j), (1.-watwgt)*qavail )

                     end if
                  end do !- k
               end if
            end if
         end do
         end do
      end if !- l_build_gsdcloud


!---SACA 3b: Cloud building w/ madwrf algorithm
!---MADWRF CF/QX ADD---
      if( l_build_madwrf ) then
         if (impose_cldtopz) then
            do j = jts,jte
            do i = its,ite
               if (cldtopk(i,j) > 0 .and. cldtopk(i,j) < kte) then

                  keep_clouds_below_lcl(i,j) = .true.
                  ! Find the first cloudy level
                  call find_thresh_k_downward(cldfra(i,kts:kte,j), CLDFRA_THRESH, cldfra_thresh_k(i,j), cldtopk(i,j), kts, kts, kte)
                  if(present(diag)) diag(:,3,:) = real( cldfra_thresh_k(:,:) ) !BJJ, diag, k-index for bg cloud top, after cloud clearing
                  if (cldfra_thresh_k(i,j) > 0) then
                     ! Extend cloud until cldtopz
                     if (cldfra_thresh_k(i,j) < cldtopk(i,j)) cldfra(i,cldfra_thresh_k(i,j) + 1:cldtopk(i,j),j) = CLDFRA_DEF
                     if(present(diag) .and. cldfra_thresh_k(i,j) < cldtopk(i,j)) diag(i,31,j) = 1.0 !BJJ, diag, build cloud up-to obs_cloud_top w/ madwrf
                  else
                     ! No clouds present
                     if (CLDTHICK_DEF > 0.0) then
                        cldfra(i,cldthick_bot_k(i,j):cldtopk(i,j),j) = CLDFRA_DEF
                     else
                        if (rh(i,cldtopk(i,j),j) > RH_NOCLOUD) cldfra(i,cldtopk(i,j),j) = CLDFRA_DEF
                        if(present(diag) .and. rh(i,cldtopk(i,j),j) > RH_NOCLOUD) diag(i,32,j) = 1.0 !BJJ diag, build cloud over no bg cloud w/ madwrf
                     end if
                  end if
!        call mpas_log_write('qc444 $i $i $i $i $r $r $r $r', intArgs=(/i,j,cldthick_bot_k(i,j),cldfra_thresh_k(i,j)/), realArgs=(/rh(i,cldtopk(i,j),j),cldfra(i,cldfra_thresh_k(i,j) + 1:cldtopk(i,j),j), cldfra(i,cldthick_bot_k(i,j):cldtopk(i,j),j),cldfra(i,cldtopk(i,j),j)/))
               end if
            end do
            end do
         end if

!..Prepare for a 1-D column to find various cloud layers.

         DO j = jts,jte
         DO i = its,ite
!        if (i.gt.10.and.i.le.20 .and. j.gt.10.and.j.le.20) then
!          debug_flag = .true.
!        else
!           debug_flag = .false.
!        endif

!        if (rand_perturb_on .eq. 1) then
!           entrmnt = MAX(0.01, MIN(0.99, 0.5 + rand_pert(i,1,j)*0.5))
!        else
            entrmnt = 0.5
!        endif

            DO k = kts,kte
               qvs1d(k) = qvsat(i,k,j)
               cfr1d(k) = cldfra(i,k,j)
               T1d(k) = t(i,k,j)
               P1d(k) = p(i,k,j)
               R1d(k) = rho(i,k,j)
               qc1d(k) = qc(i,k,j)
               qi1d(k) = qi(i,k,j)
               qs1d(k) = qs(i,k,j)
               qv1d(k) = qv(i,k,j)
            ENDDO

!     if (debug_flag) then
!       WRITE (dbg_msg,*) 'DEBUG-GT: finding cloud layers at point  (', i, ', ', j, ')'
!       CALL wrf_debug (150, dbg_msg)
!     endif

            if (is_tropo_init) then
               call find_cloudLayers(qvs1d, cfr1d, T1d, P1d, R1d, entrmnt,    &
     &                               debug_flag, qc1d, qi1d, qs1d, kts,kte, keep_clouds_below_lcl(i,j), k_tropo(i,j))
            else
               call find_cloudLayers(qvs1d, cfr1d, T1d, P1d, R1d, entrmnt,    &
     &                               debug_flag, qc1d, qi1d, qs1d, kts,kte, keep_clouds_below_lcl(i,j))
            end if

            DO k = kts,kte
               cldfra(i,k,j) = cfr1d(k)
               qc(i,k,j) = qc1d(k)
               qi(i,k,j) = qi1d(k)
               qs(i,k,j) = qs1d(k)
            ENDDO

         end do
         end do
      end if ! l_build_madwrf


!---SACA 4: Update the cloud number concentrations
!           madwrf literature does not modify these quantities.
!           gsdcloud literature modifies (Nc, Ni), based on aerosol-aware Thompson MP
!           In MPAS, only Ni (and Nr) is available as a prognostic variable.
!           So here we only modify Ni.
      do j = jts, jte
      do k = kts, kte
      do i = its, ite
         !- If there is a change in CF
         if( cldfra(i,k,j) .ne. cldfra_bk(i,k,j) ) then

            !- If CF is added, set constant Ni
            if( cldfra(i,k,j) .eq. CLDFRA_DEF ) then
               ni(i,k,j) = max( ni(i,k,j), 1.0e6_RKIND )
            end if

            !- If CF is removed, remove the background Ni
            if( cldfra(i,k,j) .lt. cldfra_bk(i,k,j) .and. &
                cldfra(i,k,j) .eq. 0.0                    ) then
               ni(i,k,j) = 0.0_RKIND
            end if

         end if
      enddo
      enddo
      enddo

!---SACA 5a: apply saturation w/ madwrf algorithm
!          : also update the cloud number concentrations (Nc, Ni)
      if (l_saturate_qv) then
         do j = jts, jte
         do k = kts, kte
         do i = its, ite
            if (cldfra(i,k,j).gt.0.20 .and. cldfra(i,k,j).lt.1.0) then
               qv(i,k,j) = MAX( qv(i,k,j), qvsat(i,k,j) )
            endif
         enddo
         enddo
         enddo
      endif


!---SACA 5b: apply saturation & sub-saturation w/ gsdcloud algorithm
!          : also update the cloud number concentrations (Nc, Ni)
      if (l_conserve_thetaV) then

         DO j = jts,jte
         DO i = its,ite
         DO k = kts, kte

            !- If there is a change in CF
            if( cldfra(i,k,j) .ne. cldfra_bk(i,k,j) ) then

               !- set variables
               temp = t(i,k,j)
               cloudqvis = qvsat(i,k,j)
               q_bk = min( qv(i,k,j), cloudqvis )
               constantTv= temp * ( 1.0 + (rv/rgas-1.0) * q_bk )

               !- If CF is added, make it saturated.
               if( cldfra(i,k,j) .eq. CLDFRA_DEF ) then
                  qtemp = cloudqvis * 1.00
                  if( q_bk < cloudqvis * 1.00 ) then
                     do nnn=1,miter
                        temp=constantTv/(1.0 + (rv/rgas-1.0) * qtemp)

                        !cloudqvis = ....
                        TK   = temp
                        TC   = TK - 273.16
                        qvsw = rslf(P(i,k,j), TK)
                        qvsi = rsif(P(i,k,j), TK)
                        if (tc .ge. -12.0) then
                           cloudqvis = qvsw
                        elseif (tc .lt. -20.0) then
                           cloudqvis = qvsi
                        else
                           cloudqvis = qvsw - (qvsw-qvsi)*(-12.0-tc)/(-12.0+20.)
                        endif

                        qtemp = cloudqvis * 1.00
                     enddo
                  endif
                  !limit increment
                  q_bk = min(qtemp, q_bk+qv_max_inc)
               endif

               !- If CF is removed, make it sub-saturated. w/ at least 100mb above sfc condition.
               if( cldfra(i,k,j) .lt. cldfra_bk(i,k,j) .and. &
                   cldfra(i,k,j) .eq. 0.0              .and. &
                   ( p(i,1,j) - p(i,k,j) ) .gt. 100.0        ) then
                  if( q_bk > cloudqvis * rh_clear_p ) then
                     qtemp = cloudqvis * rh_clear_p
                     do nnn=1,miter
                        temp=constantTv/(1.0 + (rv/rgas-1.0) * qtemp)

                        !cloudqvis = ....
                        TK   = temp
                        TC   = TK - 273.16
                        qvsw = rslf(P(i,k,j), TK)
                        qvsi = rsif(P(i,k,j), TK)
                        if (tc .ge. -12.0) then
                           cloudqvis = qvsw
                        elseif (tc .lt. -20.0) then
                           cloudqvis = qvsi
                        else
                           cloudqvis = qvsw - (qvsw-qvsi)*(-12.0-tc)/(-12.0+20.)
                        endif

                        qtemp = cloudqvis * rh_clear_p
                     enddo
                  endif
                  !limit increment
                  q_bk = min(qtemp, q_bk+qv_max_inc)
               endif

               !- update vars
               t(i,k,j)  = temp
               qv(i,k,j) = q_bk
            endif

         ENDDO
         ENDDO
         ENDDO
      endif !- l_conserve_thetaV


      END SUBROUTINE cal_cldfra3

!+---+-----------------------------------------------------------------+

      SUBROUTINE find_thresh_k_downward(var, var_thresh, k_lev, k_top, k_bot, kts, kte)

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: kts, kte
      INTEGER, INTENT(IN)  :: k_top, k_bot
      REAL(KIND=RKIND), DIMENSION(kts:kte), INTENT(IN) :: var
      REAL(KIND=RKIND), INTENT(IN)     :: var_thresh
      INTEGER, INTENT(OUT) :: k_lev
      INTEGER :: k

      !! Identify the k-level where the input variable first exceeds a threshold, searching downward over a given range
      do k = k_top, k_bot, -1
         if (var(k) > var_thresh) then
            k_lev = k
            exit
         end if
      end do

      END SUBROUTINE find_thresh_k_downward

!+---+-----------------------------------------------------------------+

      SUBROUTINE find_k_from_z_agl(z_agl, k_lev, dz, kts, kte)

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: kts, kte
      REAL(KIND=RKIND), INTENT(IN)     :: z_agl
      INTEGER, INTENT(OUT) :: k_lev
      REAL(KIND=RKIND), DIMENSION(kts:kte), INTENT(IN) :: dz
      INTEGER :: k
      REAL(KIND=RKIND)    :: z_this, z_next, z_full_lev

      !! Identify the k-level of a given height AGL, starting from ground level
      z_full_lev = 0.0
      do k = kts, kte
         z_this = z_full_lev
         z_next = z_this + dz(k)
         if (z_agl > 0.0) then
            if (z_agl > z_this .and. z_agl <= z_next) then
               k_lev = k
            end if
         end if
         z_full_lev = z_next
      end do

      END SUBROUTINE find_k_from_z_agl

!+---+-----------------------------------------------------------------+
!..From cloud fraction array, find clouds of multi-level depth and compute
!.. a reasonable value of LWP or IWP that might be contained in that depth,
!.. unless existing LWC/IWC is already there.

      SUBROUTINE find_cloudLayers(qvs1d, cfr1d, T1d, P1d, R1d, entrmnt, &
     &                            debugfl, qc1d, qi1d, qs1d, kts,kte, keep_clouds_below_lcl, ktropo)
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: kts, kte
      LOGICAL, INTENT(IN):: debugfl, keep_clouds_below_lcl
      REAL(KIND=RKIND), INTENT(IN):: entrmnt
      REAL(KIND=RKIND), DIMENSION(kts:kte), INTENT(IN):: qvs1d,T1d,P1d,R1d
      REAL(KIND=RKIND), DIMENSION(kts:kte), INTENT(INOUT):: cfr1d
      REAL(KIND=RKIND), DIMENSION(kts:kte), INTENT(INOUT):: qc1d, qi1d, qs1d
      integer, intent(in), optional :: ktropo

!..Local vars.
      REAL(KIND=RKIND), DIMENSION(kts:kte):: theta, dz
      REAL(KIND=RKIND):: Z1, Z2, theta1, theta2, ht1, ht2
      INTEGER:: k, k2, k_tropo, k_m12C, k_m40C, k_cldb, k_cldt, kbot
      LOGICAL:: in_cloud
      character*512 dbg_msg

      logical :: is_tropo_init
!+---+

      if (present(ktropo)) then
        is_tropo_init = .true.
      else
        is_tropo_init = .false.
      end if

      k_m12C = 0
      k_m40C = 0
      DO k = kte, kts, -1
         theta(k) = T1d(k)*((100000.0/P1d(k))**(287.05/1004.))
         if (T1d(k)-273.16 .gt. -40.0) k_m40C = MAX(k_m40C, k)
         if (T1d(k)-273.16 .gt. -12.0) k_m12C = MAX(k_m12C, k)
      ENDDO
      if (k_m40C .le. kts) k_m40C = kts
      if (k_m12C .le. kts) k_m12C = kts

      Z2 = 44307.692 * (1.0 - (P1d(kte)/101325.)**0.190)
      DO k = kte-1, kts, -1
         Z1 = 44307.692 * (1.0 - (P1d(k)/101325.)**0.190)
         dz(k+1) = Z2 - Z1
         Z2 = Z1
      ENDDO
      dz(kts) = dz(kts+1)

!..Find tropopause height, best surrogate, because we would not really
!.. wish to put fake clouds into the stratosphere.  The 10/1500 ratio
!.. d(Theta)/d(Z) approximates a vertical line on typical SkewT chart
!.. near typical (mid-latitude) tropopause height.  Since messy data
!.. could give us a false signal of such a transition, do the check over
!.. three K-level change, not just a level-to-level check.  This method
!.. has potential failure in arctic-like conditions with extremely low
!.. tropopause height, as would any other diagnostic, so ensure resulting
!.. k_tropo level is above 4km.

      DO k = kte-3, kts, -1
         theta1 = theta(k)
         theta2 = theta(k+2)
         ht1 = 44307.692 * (1.0 - (P1d(k)/101325.)**0.190)
         ht2 = 44307.692 * (1.0 - (P1d(k+2)/101325.)**0.190)
         if ( (((theta2-theta1)/(ht2-ht1)) .lt. 10./1500. ) .AND.       &
     &                       (ht1.lt.19000.) .and. (ht1.gt.4000.) ) then
            goto 86
         endif
      ENDDO
 86   continue
!      k_tropo = MAX(kts+2, k+2)

      if (is_tropo_init) then
        k_tropo = ktropo
      else
        k_tropo = MAX(kts+2, k+2)
      end if

!     if (debugfl) then
!     print*, ' FOUND TROPOPAUSE ', k_tropo, ' near ', ht2, ' m'
!       WRITE (dbg_msg,*) 'DEBUG-GT: FOUND TROPOPAUSE ', k_tropo, ' near ', ht2, ' m'
!       CALL wrf_debug (150, dbg_msg)
!     endif

!..Eliminate possible fractional clouds above supposed tropopause.
      DO k = k_tropo+1, kte
         if (cfr1d(k).gt.0.0 .and. cfr1d(k).lt.0.999) then
            cfr1d(k) = 0.
         endif
      ENDDO

!..We would like to prevent fractional clouds below LCL in idealized
!.. situation with deep well-mixed convective PBL, that otherwise is
!.. likely to get clouds in more realistic capping inversion layer.

      kbot = kts+2
      DO k = kbot, k_m12C
         if ( (theta(k)-theta(k-1)) .gt. 0.05E-3*dz(k)) EXIT
      ENDDO
      kbot = MAX(kts+1, k-2)

!      DO k = kts, kbot
!         if (cfr1d(k).gt.0.0 .and. cfr1d(k).lt.0.999) cfr1d(k) = 0.
!      ENDDO

      if (.not. keep_clouds_below_lcl) then
        DO k = kts, kbot
!           print *, 'Reducing cloud below LCL', k
           if (cfr1d(k).gt.0.0 .and. cfr1d(k).lt.1.0) cfr1d(k) = 0.5*cfr1d(k)
        ENDDO
      else
        kbot = kts + 1
      end if


!..Starting below tropo height, if cloud fraction greater than 1 percent,
!.. compute an approximate total layer depth of cloud, determine a total
!.. liquid water/ice path (LWP/IWP), then reduce that amount with tuning
!.. parameter to represent entrainment factor, then divide up LWP/IWP
!.. into delta-Z weighted amounts for individual levels per cloud layer.

      k_cldb = k_tropo
      in_cloud = .false.
      k = k_tropo

      DO WHILE (.not. in_cloud .AND. k.gt.k_m12C)
         k_cldt = 0
         if (cfr1d(k).ge.0.01) then
            in_cloud = .true.
            k_cldt = MAX(k_cldt, k)
         endif
         if (in_cloud) then
            DO k2 = k_cldt-1, k_m12C, -1
               if (cfr1d(k2).lt.0.01 .or. k2.eq.k_m12C) then
                  k_cldb = k2+1
                  goto 87
               endif
            ENDDO
 87         continue
            in_cloud = .false.
         endif
         if ((k_cldt - k_cldb + 1) .ge. 2) then
!     if (debugfl) then
!           print*, 'An ice cloud layer is found between ', k_cldt, k_cldb, P1d(k_cldt)*0.01, P1d(k_cldb)*0.01
!       WRITE (dbg_msg,*) 'DEBUG-GT: An ice cloud layer is found between ', k_cldt, k_cldb, P1d(k_cldt)*0.01, P1d(k_cldb)*0.01
!       CALL wrf_debug (150, dbg_msg)
!     endif
            call adjust_cloudIce(cfr1d, qi1d, qs1d, qvs1d, T1d,R1d,dz,  &
     &                           entrmnt, k_cldb,k_cldt,kts,kte)
            k = k_cldb
         else
            if (cfr1d(k_cldb).gt.0.and.qi1d(k_cldb).lt.1.E-6)           &
     &               qi1d(k_cldb)=1.E-5*cfr1d(k_cldb)
         endif
         k = k - 1
      ENDDO


      k_cldb = k_tropo
      in_cloud = .false.

!     k = k_m12C + 2
      k = min(k_m12C,k_m12C+2)
      DO WHILE (.not. in_cloud .AND. k.gt.kbot)
         k_cldt = 0
         if (cfr1d(k).ge.0.01) then
            in_cloud = .true.
            k_cldt = MAX(k_cldt, k)
         endif
         if (in_cloud) then
            DO k2 = k_cldt-1, kbot, -1
               if (cfr1d(k2).lt.0.01 .or. k2.eq.kbot) then
                  k_cldb = k2+1
                  goto 88
               endif
            ENDDO
 88         continue
            in_cloud = .false.
         endif
         if ((k_cldt - k_cldb + 1) .ge. 2) then
!     if (debugfl) then
!           print*, 'A water cloud layer is found between ', k_cldt, k_cldb, P1d(k_cldt)*0.01, P1d(k_cldb)*0.01
!       WRITE (dbg_msg,*) 'DEBUG-GT: A water cloud layer is found between ', k_cldt, k_cldb, P1d(k_cldt)*0.01, P1d(k_cldb)*0.01
!       CALL wrf_debug (150, dbg_msg)
!     endif
            call adjust_cloudH2O(cfr1d, qc1d, qvs1d, T1d,R1d,dz,        &
     &                           entrmnt, k_cldb,k_cldt,kts,kte)
            k = k_cldb
         else
            if (cfr1d(k_cldb).gt.0.and.qc1d(k_cldb).lt.1.E-6)           &
     &               qc1d(k_cldb)=1.E-5*cfr1d(k_cldb)
         endif
         k = k - 1
      ENDDO

!..Do a final total column adjustment since we may have added more than 1mm
!.. LWP/IWP for multiple cloud decks.

      call adjust_cloudFinal(cfr1d, qc1d, qi1d, R1d,dz, kts,kte,k_tropo)

!     if (debugfl) then
!     print*, ' Made-up fake profile of clouds'
!     do k = kte, kts, -1
!        write(*,'(i3, 2x, f8.2, 2x, f9.2, 2x, f6.2, 2x,  f15.7, 2x, f15.7)') &
!    &        K, T1d(k)-273.15, P1d(k)*0.01, cfr1d(k)*100., qc1d(k)*1000.,qi1d(k)*1000.
!     enddo
!       WRITE (dbg_msg,*) 'DEBUG-GT:  Made-up fake profile of clouds'
!       CALL wrf_debug (150, dbg_msg)
!       do k = kte, kts, -1
!          write(dbg_msg,'(f8.2, 2x, f9.2, 2x, f6.2, 2x,  f15.7, 2x, f15.7)') &
!    &          T1d(k)-273.15, P1d(k)*0.01, cfr1d(k)*100., qc1d(k)*1000.,qi1d(k)*1000.
!          CALL wrf_debug (150, dbg_msg)
!       enddo
!     endif


      END SUBROUTINE find_cloudLayers

!+---+-----------------------------------------------------------------+

      SUBROUTINE adjust_cloudIce(cfr,qi,qs,qvs, T,Rho,dz, entr, k1,k2,kts,kte)
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: k1,k2, kts,kte
      REAL(KIND=RKIND), INTENT(IN):: entr
      REAL(KIND=RKIND), DIMENSION(kts:kte), INTENT(IN):: cfr, qvs, T, Rho, dz
      REAL(KIND=RKIND), DIMENSION(kts:kte), INTENT(INOUT):: qi, qs
      REAL(KIND=RKIND):: iwc, max_iwc, tdz, this_iwc, this_dz, iwp_exists
      INTEGER:: k, kmid

      tdz = 0.
      do k = k1, k2
         tdz = tdz + dz(k)
      enddo
      kmid = NINT(0.5*(k1+k2))
      max_iwc = ABS(qvs(k2-1)-qvs(k1))
!     print*, ' max_iwc = ', max_iwc, ' over DZ=',tdz

      iwp_exists = 0.
      do k = k1, k2
         iwp_exists = iwp_exists + (qi(k)+qs(k))*Rho(k)*dz(k)
      enddo
      if (iwp_exists .gt. 1.0) RETURN

      this_dz = 0.0
      do k = k1, k2
         if (k.eq.k1) then
            this_dz = this_dz + 0.5*dz(k)
         else
            this_dz = this_dz + dz(k)
         endif
         this_iwc = max_iwc*this_dz/tdz
         iwc = MAX(1.E-6, this_iwc*(1.-entr))
         if (cfr(k).gt.0.01.and.cfr(k).lt.0.99.and.T(k).ge.203.16) then
            qi(k) = qi(k) + 0.1*cfr(k)*iwc
         elseif (qi(k).lt.1.E-5.and.cfr(k).ge.0.99.and.T(k).ge.203.16) then
            qi(k) = qi(k) + 0.01*iwc
         endif
      enddo

      END SUBROUTINE adjust_cloudIce

!+---+-----------------------------------------------------------------+

      SUBROUTINE adjust_cloudH2O(cfr, qc, qvs, T,Rho,dz, entr, k1,k2,kts,kte)
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: k1,k2, kts,kte
      REAL(KIND=RKIND), INTENT(IN):: entr
      REAL(KIND=RKIND), DIMENSION(kts:kte):: cfr, qc, qvs, T, Rho, dz
      REAL(KIND=RKIND):: lwc, max_lwc, tdz, this_lwc, this_dz, lwp_exists
      INTEGER:: k, kmid

      tdz = 0.
      do k = k1, k2
         tdz = tdz + dz(k)
      enddo
      kmid = NINT(0.5*(k1+k2))
      max_lwc = ABS(qvs(k2-1)-qvs(k1))
!     print*, ' max_lwc = ', max_lwc, ' over DZ=',tdz

      lwp_exists = 0.
      do k = k1, k2
         lwp_exists = lwp_exists + qc(k)*Rho(k)*dz(k)
      enddo
      if (lwp_exists .gt. 1.0) RETURN

      this_dz = 0.0
      do k = k1, k2
         if (k.eq.k1) then
            this_dz = this_dz + 0.5*dz(k)
         else
            this_dz = this_dz + dz(k)
         endif
         this_lwc = max_lwc*this_dz/tdz
         lwc = MAX(1.E-6, this_lwc*(1.-entr))
         if (cfr(k).gt.0.01.and.cfr(k).lt.0.99.and.T(k).lt.298.16.and.T(k).ge.253.16) then
            qc(k) = qc(k) + cfr(k)*cfr(k)*lwc
         elseif (cfr(k).ge.0.99.and.qc(k).lt.1.E-5.and.T(k).lt.298.16.and.T(k).ge.253.16) then
            qc(k) = qc(k) + 0.1*lwc
         endif
      enddo

      END SUBROUTINE adjust_cloudH2O

!+---+-----------------------------------------------------------------+

!..Do not alter any grid-explicitly resolved hydrometeors, rather only
!.. the supposed amounts due to the cloud fraction scheme.

      SUBROUTINE adjust_cloudFinal(cfr, qc, qi, Rho,dz, kts,kte,k_tropo)
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: kts,kte,k_tropo
      REAL(KIND=RKIND), DIMENSION(kts:kte), INTENT(IN):: cfr, Rho, dz
      REAL(KIND=RKIND), DIMENSION(kts:kte), INTENT(INOUT):: qc, qi
      REAL(KIND=RKIND):: lwp, iwp, xfac
      INTEGER:: k

      lwp = 0.
      do k = kts, k_tropo
         if (cfr(k).gt.0.01 .and. cfr(k).lt.0.99) then
            lwp = lwp + qc(k)*Rho(k)*dz(k)
         endif
      enddo

      iwp = 0.
      do k = kts, k_tropo
         if (cfr(k).gt.0.01 .and. cfr(k).lt.0.99) then
            iwp = iwp + qi(k)*Rho(k)*dz(k)
         endif
      enddo

      if (lwp .gt. 1.0) then
         xfac = 1./lwp
         do k = kts, k_tropo
            if (cfr(k).gt.0.01 .and. cfr(k).lt.0.99) then
               qc(k) = qc(k)*xfac
            endif
         enddo
      endif

      if (iwp .gt. 1.0) then
         xfac = 1./iwp
         do k = kts, k_tropo
            if (cfr(k).gt.0.01 .and. cfr(k).lt.0.99) then
               qi(k) = qi(k)*xfac
            endif
         enddo
      endif

      END SUBROUTINE adjust_cloudFinal

!=================================================================================================================
  end module module_mp_thompson_cldfra3_saca
!=================================================================================================================
