module mpas_saca_interface_mod

!MPAS-Model
use mpas_kind_types
use mpas_pool_routines
use mpas_derived_types

!saca
use module_mp_thompson_cldfra3_saca, only: cal_cldfra3

!mpas-jedi
use mpas_fields_mod

implicit none
private

public :: update_cloud_fields

contains

subroutine update_cloud_fields ( state, obs )
   implicit none

   class(mpas_fields), intent(inout) :: state !< state
   class(mpas_fields), intent(in   ) :: obs   !< state

   !-- of MPAS pool
   real(kind=RKIND),dimension(:),pointer:: meshDensity, xland, ter
   real(kind=RKIND),dimension(:),pointer:: cldmask, brtemp
   real(kind=RKIND),dimension(:,:),pointer:: qv, qc, qi, qs, cldfrac
   real(kind=RKIND),dimension(:,:),pointer:: rho, t, p
   real(kind=RKIND),pointer:: len_disp

   type (mpas_pool_type), pointer :: meshPool
   type (block_type), pointer :: block_ptr

   !-- for physics interface
   integer:: ids,ide,jds,jde,kds,kde
   integer:: ims,ime,jms,jme,kms,kme
   integer:: its,ite,jts,jte,kts,kte

   real(kind=RKIND),dimension(:,:),allocatable :: dx_p, xland_p, ter_p
   real(kind=RKIND),dimension(:,:),allocatable :: cldmask_p, brtemp_p
   real(kind=RKIND),dimension(:,:,:),allocatable :: qv_p, qc_p, qi_p, qs_p, cldfrac_p
   real(kind=RKIND),dimension(:,:,:),allocatable :: rho_p, t_p, p_p
   real(kind=RKIND),dimension(:,:,:),allocatable :: dz_p
   integer,dimension(:,:),allocatable :: k_tropo_p
   real(kind=RKIND),dimension(:,:),allocatable :: tropoz_p, cldtopz_p, cldtopz_p0

   integer :: i, j, k
   integer :: k_tropo
   logical :: debug_flag, modify_qvapor
   real(kind=RKIND) :: tropo_z, this_height

!some configurable variable
   modify_qvapor = .true.

!set dimensions
!originally from "subroutine physics_run_init" of "mpas_atmphys_manager.F"
!WS put these "subroutine init_atm_setup_case" of "mpas_init_atm_cases.F"
   ims=1   ; ime=state % geom % nCellsSolve
   jms=1   ; jme=1
   kms=1   ; kme=state % geom % nVertLevels + 1
   ids=ims ; ide=ime + 1
   jds=jms ; jde=jme + 1
   kds=kms ; kde=kme 
   its=ims ; ite=ime 
   jts=jms ; jte=jme
   kts=kms ; kte=kme - 1

!allocate
!originally from "subroutine allocate_forall_physics" of "mpas_atmphys_interface.F"
   allocate(dx_p(ims:ime,jms:jme))
   allocate(xland_p(ims:ime,jms:jme))
   allocate(ter_p(ims:ime,jms:jme))
   allocate(cldmask_p(ims:ime,jms:jme))
   allocate(brtemp_p(ims:ime,jms:jme))
   allocate(qv_p(ims:ime,kms:kme,jms:jme))
   allocate(qc_p(ims:ime,kms:kme,jms:jme))
   allocate(qi_p(ims:ime,kms:kme,jms:jme))
   allocate(qs_p(ims:ime,kms:kme,jms:jme))
   allocate(cldfrac_p(ims:ime,kms:kme,jms:jme))
   allocate(rho_p(ims:ime,kms:kme,jms:jme))
   allocate(t_p(ims:ime,kms:kme,jms:jme))
   allocate(p_p(ims:ime,kms:kme,jms:jme))
   allocate(dz_p(ims:ime,kms:kme,jms:jme))
   allocate(k_tropo_p(ims:ime,jms:jme))
   allocate(tropoz_p(ims:ime,jms:jme))
   allocate(cldtopz_p(ims:ime,jms:jme))
   allocate(cldtopz_p0(ims:ime,jms:jme))

!pass the pool variables to physics variables
!originally from "subroutine MPAS_to_physics" of "mpas_atmphys_interface.F"
   block_ptr => state % geom % domain % blocklist
   call mpas_pool_get_subpool(block_ptr % structs,'mesh',meshPool)
   call mpas_pool_get_config(block_ptr % configs,'config_len_disp',len_disp)
   call mpas_pool_get_array(meshPool,'meshDensity',meshDensity)
   call state%get('xland',xland)
   call obs%get('cldmask',cldmask)
   call obs%get('brtemp',brtemp)

   do j = jts,jte
      do i = its,ite
         dx_p(i,j)      = len_disp / meshDensity(i)**0.25
         xland_p(i,j)   = xland(i)
         ter_p(i,j)     = state % geom % zgrid(1,i)  !ter(i) BJJ
         cldmask_p(i,j) = cldmask(i)
         brtemp_p(i,j)  = brtemp(i)
      enddo
   enddo

   call state%get('qv', qv)
   call state%get('qc', qc)
   call state%get('qi', qi)
   call state%get('qs', qs)
   call state%get('cldfrac'    , cldfrac)
   call state%get('rho'        , rho)
   call state%get('temperature', t)
   call state%get('pressure'   , p)
   do j = jts, jte
   do k = kts, kte
   do i = its, ite
      !water vapor and moist arrays:
      qv_p(i,k,j) = max(0.,qv(k,i))
      qc_p(i,k,j) = max(0.,qc(k,i))
      qi_p(i,k,j) = max(0.,qi(k,i))
      qs_p(i,k,j) = max(0.,qs(k,i))

      cldfrac_p(i,k,j) = 0._RKIND

      rho_p(i,k,j) = rho(k,i)*(1._RKIND + qv_p(i,k,j))
      t_p(i,k,j)   = t(k,i)
      p_p(i,k,j)   = p(k,i)
      dz_p(i,k,j)   = state % geom % zgrid(k+1,i) - state % geom % zgrid(k,i)
   enddo
   enddo
   enddo

!determine the tropopause height from Model
   debug_flag = .false.
   do j = jts,jte
   do i = its,ite
      call Calc_tropo_height (t_p(i,kts:kte-1,j), p_p(i,kts:kte-1,j), dz_p(i,kts:kte-1,j), kts, kte-1, debug_flag, k_tropo, tropo_z)
      k_tropo_p(i,j) = k_tropo
   end do
   end do

!determine the cloud top height from Obs
   debug_flag = .false.
   do j = jts,jte
   do i = its,ite
      !..Staring at model top, go downwards first to the level of the tropopause,
      !.. then keep going until the model temperature is greater (warmer) than the
      !.. incoming GOES satellite longwave IR tempeature (channel13 on GOES-R series).
      !.. As version1.0, this will be the highest altitude of potential cloud based on
      !.. IR temperature whereas an improvement is possible for inversions in which the
      !.. cloud top could be lower altitude that might be detected using the RH field.

      tropoz_p(i,j) = tropo_z + ter_p(i,j) 

      cldtopz_p0(i,j) = -999.
      if (cldmask_p(i,j) .gt. 0.) then
         this_height = SUM(dz_p(i,kts:kte-1,j)) + ter_p(i,j)
         do k = kte-2, kts, -1
            this_height = this_height - dz_p(i,k+1,j)
            if (this_height .gt. tropoz_p(i,j)) CYCLE
            if (t_p(i,k,j) .gt. brtemp_p(i,j)) then
               cldtopz_p0(i,j) = this_height + dz_p(i,k,j)  &
                  * (brtemp_p(i,j)-t_p(i,k,j)) / (t_p(i,k+1,j)-t_p(i,k,j))

               if ( (cldtopz_p0(i,j).lt.this_height) .OR. &
                    (cldtopz_p0(i,j).gt.this_height+dz_p(i,k,j)) ) then
                  cldtopz_p0(i,j) = tropoz_p(i,j)
               end if
               EXIT
            endif
         enddo
      endif
   end do
   end do

   do j = jts,jte
   do i = its,ite
      if (cldtopz_p0(i,j) > 0.0) then
         cldtopz_p(i,j) = cldtopz_p0(i, j) - ter_p(i, j)
      else
         cldtopz_p(i,j) = cldtopz_p0(i, j)
      end if
   end do
   end do

!call cal_cldfra3
   call cal_cldfra3( &
         cldfra = cldfrac_p , qv     = qv_p    , qc = qc_p , qi  = qi_p    ,       &
         qs     = qs_p   ,    dz     = dz_p    , p = p_p   , t   = t_p     , rho = rho_p   ,    &
         xland  = xland_p   , gridkm = dx_p    , modify_qvapor = modify_qvapor ,   &
         ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,   &
         ims = ims , ime = ime , jms = jms , jme = jme , kms = kds , kme = kme ,   &
         its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte ,   &
         k_tropo = k_tropo_p, cldmask = cldmask_p, cldtopz = cldtopz_p             &
         )

!update the pool variables
   call state%get('cldfrac',cldfrac)
   do j = jts,jte
   do k = kts,kte
   do i = its,ite
      cldfrac(k,i) = cldfrac_p(i,k,j)
      qc(k,i) = qc_p(i,k,j)
      qi(k,i) = qi_p(i,k,j)
      qs(k,i) = qs_p(i,k,j)
      if (modify_qvapor) then
         qv(k,i) = qv_p(i,k,j)
      end if
   end do
   end do
   end do

!deallocate
!originally from "subroutine deallocate_forall_physics" of "mpas_atmphys_interface.F"
   deallocate(dx_p)
   deallocate(xland_p)
   deallocate(ter_p)
   deallocate(cldmask_p)
   deallocate(brtemp_p)
   deallocate(qv_p)
   deallocate(qc_p)
   deallocate(qi_p)
   deallocate(qs_p)
   deallocate(cldfrac_p)
   deallocate(rho_p)
   deallocate(t_p)
   deallocate(p_p)
   deallocate(dz_p)
   deallocate(k_tropo_p)
   deallocate(tropoz_p)
   deallocate(cldtopz_p)
   deallocate(cldtopz_p0)

end subroutine update_cloud_fields

subroutine Calc_tropo_height (T1d, P1d, dz1d, kts, kte, debugfl, k_tropo, tropo_z)

   !..Find tropopause height, best surrogate, because we would not really
   !.. wish to put fake clouds into the stratosphere.  The 10/1500 ratio
   !.. d(Theta)/d(Z) approximates a vertical line on typical SkewT chart
   !.. near typical (mid-latitude) tropopause height.  Since messy data
   !.. could give us a false signal of such a transition, do the check over 
   !.. three K-level change, not just a level-to-level check.  This method
   !.. has potential failure in arctic-like conditions with extremely low
   !.. tropopause height, as would any other diagnostic, so ensure resulting
   !.. k_tropo level is above 700hPa.

   implicit none

   REAL(kind=RKIND), DIMENSION(kts:kte), INTENT(IN) :: T1d, P1d, Dz1d
   integer, intent(in) :: kts, kte
   LOGICAL, INTENT(IN):: debugfl
   integer, intent(out) :: k_tropo
   real(kind=RKIND), intent (out) :: tropo_z

     ! Local vars
   integer :: k, k_p200
   real(kind=RKIND) :: theta1, theta2, delz
   character*512 dbg_msg
   REAL(kind=RKIND), DIMENSION(kts:kte):: theta

   k_p200 = 0
   DO k = kte, kts, -1
      theta(k) = T1d(k)*((100000.0/P1d(k))**(287.05/1004.))
      if (P1d(k).gt.19999.0 .and. k_p200.eq.0) k_p200 = k
   END DO

   if ( (kte-k_p200) .lt. 3) k_p200 = kte-3
   DO k = k_p200-2, kts, -1
      theta1 = theta(k)
      theta2 = theta(k+2)
      delz = 0.5*dz1d(k) + dz1d(k+1) + 0.5*dz1d(k+2)
      if ( (((theta2-theta1)/delz).lt.10./1500.) .OR. P1d(k).gt.70000.) EXIT
   ENDDO
   k_tropo = MAX(kts+2, MIN(k+2, kte-1))

   if (k_tropo .gt. k_p200) then
      DO k = kte-3, k_p200-2, -1
         theta1 = theta(k)
         theta2 = theta(k+2)
         delz = 0.5*dz1d(k) + dz1d(k+1) + 0.5*dz1d(k+2)
         if ( (((theta2-theta1)/delz).lt.10./1500.) .AND. P1d(k).gt.9500.) EXIT
      ENDDO
      k_tropo = MAX(k_p200-1, MIN(k+2, kte-1))
   endif
   tropo_z = SUM(dz1d(kts:k_tropo))

   if (k_tropo.gt.kte-2) then
      WRITE (dbg_msg,*) 'DEBUG-GT: CAUTION, tropopause appears to be very high up: ', k_tropo
      call mpas_log_write(dbg_msg)
      do k = kte, kts, -1
         WRITE (dbg_msg,*) 'DEBUG-GT,   P, T : ', k,P1d(k)*0.01,T1d(k)-273.16
         call mpas_log_write(dbg_msg)
      enddo
   elseif (debugfl) then
      WRITE (dbg_msg,*) 'DEBUG-GT: FOUND TROPOPAUSE k=', k_tropo
      call mpas_log_write(dbg_msg)
   endif

end subroutine Calc_tropo_height

end module mpas_saca_interface_mod
