   !***********************************************************************
   ! Copyright (c) 2018, National Atmospheric for Atmospheric Research (NCAR).
   !
   ! Unless noted otherwise source code is licensed under the BSD license.
   ! Additional copyright and license information can be found in the LICENSE file
   ! distributed with this code, or at http://mpas-dev.github.com/license.html
   !
   !-----------------------------------------------------------------------

   module mpas2ufo_vars_mod

   !***********************************************************************
   !
   !  Module mpas2ufo_vars_mod to convert/extract variables from MPAS  
   !  needed for Data assimilation purpose.
   !> \author  Gael Descombes/Mickael Duda NCAR/MMM
   !> \date    Februray 2018
   !
   !-----------------------------------------------------------------------

!oops
use kinds, only : kind_real

!ufo
use gnssro_mod_transform, only: geometric2geop
use ufo_vars_mod

!MPAS-Model
use atm_core
use mpas_abort, only : mpas_dmpar_global_abort
use mpas_constants, only : gravity, rgas, rv, cp
use mpas_dmpar
use mpas_derived_types
use mpas_field_routines
use mpas_pool_routines

!MPAS-JEDI
use mpas_constants_mod
use mpas_geom_mod
use mpas_kinds, only : kind_double

private

public :: convert_type_soil, convert_type_veg
public :: uv_to_wdir
public :: w_to_q, q_to_w
public :: theta_to_temp, temp_to_theta
public :: twp_to_rho, hydrostatic_balance

public :: update_mpas_field
public :: convert_mpas_field2ufo,   &
          convert_mpas_field2ufoTL, &
          convert_mpas_field2ufoAD

contains

!-------------------------------------------------------------------------------------------

subroutine update_mpas_field(domain, pool_a)

   implicit none

   type (domain_type), pointer, intent(in) :: domain
   type (mpas_pool_type), pointer, intent(inout) :: pool_a

   type (mpas_pool_type), pointer :: state, diag, mesh, pool_b
   type (mpas_pool_iterator_type) :: poolItr
   integer, parameter :: maxfield = 1
   character (len=22) :: fieldname(1:maxfield)
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   integer :: ii

   call mpas_pool_begin_iteration(pool_a)
   call mpas_pool_get_subpool(domain % blocklist % structs,'state',state)
   call mpas_pool_get_subpool(domain % blocklist % structs, 'diag', diag)
   call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh) 
   call atm_compute_output_diagnostics(state, 1, diag, mesh)

   fieldname(1) = 'pressure' 

   do while ( mpas_pool_get_next_member(pool_a, poolItr) )

     ! Pools may in general contain dimensions, namelist options, fields, or other pools,
     ! so we select only those members of the pool that are fields
     if (poolItr % memberType == MPAS_POOL_FIELD) then
       ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
       if (poolItr % dataType == MPAS_POOL_REAL) then
!          write(*,*)'Looking for ',trim(poolItr % memberName)
          do ii=1, maxfield
             if ( trim(fieldname(ii)).eq.trim(poolItr % memberName)) then
               if (poolItr % nDims == 0) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r0d_ptr_a)
                  call mpas_pool_get_array(diag, trim(poolItr % memberName), r0d_ptr_b)
                  r0d_ptr_a = r0d_ptr_a + r0d_ptr_b
               else if (poolItr % nDims == 1) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
                  call mpas_pool_get_array(diag, trim(poolItr % memberName), r1d_ptr_b)
                  r1d_ptr_a = r1d_ptr_b
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_a)
                  call mpas_pool_get_array(diag, trim(poolItr % memberName), r2d_ptr_b)
                  r2d_ptr_a = r2d_ptr_b
!                  write(0,*)'Update MIN/MAX: ',trim(fieldname(ii)),minval(r2d_ptr_a),maxval(r2d_ptr_a)
               end if
             end if
          end do
       end if
     end if

   end do

end subroutine update_mpas_field

!-------------------------------------------------------------------------------------------

   !-- from WRFDA/da_crtm.f90
   integer function convert_type_soil(type_in)
   implicit none
   integer, intent(in) :: type_in
   integer, parameter :: n_soil_type = 16  ! wrf num_soil_cat
   integer, parameter :: wrf_to_crtm_soil(n_soil_type) = &
      (/ 1, 1, 4, 2, 2, 8, 7, 2, 6, 5, 2, 3, 8, 1, 6, 9 /)

   convert_type_soil = max( 1, wrf_to_crtm_soil(type_in) )

   end function convert_type_soil

   integer function convert_type_veg(type_in)
   integer, intent(in) :: type_in
   integer, parameter :: USGS_n_type = 24 
   integer, parameter :: IGBP_n_type = 20 
   ! vegetation type mapping for GFS classification scheme
   ! REL-2.1.3.CRTM_User_Guide.pdf table 4.16
   integer, parameter :: usgs_to_crtm_mw(USGS_n_type) = &
      (/  7, 12, 12, 12, 12, 12,  7,  9,  8,  6, &
          2,  5,  1,  4,  3,  0,  8,  8, 11, 10, &
         10, 10, 11, 13 /)  
   integer, parameter :: igbp_to_crtm_mw(IGBP_n_type) = &
      (/  4,  1,  5,  2,  3,  8,  9,  6,  6,  7, &
          8, 12,  7, 12, 13, 11,  0, 10, 10, 11 /)

   !TODO: make this general: consider both dataset, usgs & igbp
   convert_type_veg = max( 1, usgs_to_crtm_mw(type_in) )

   end function convert_type_veg

!  !-- from GSI/crtm_interface.f90
!  integer(i_kind), parameter :: USGS_N_TYPES = 24
!  integer(i_kind), parameter :: IGBP_N_TYPES = 20
!  integer(i_kind), parameter :: NAM_SOIL_N_TYPES = 16
!  integer(i_kind), parameter, dimension(1:IGBP_N_TYPES) :: igbp_to_gfs=(/4, &
!    1, 5, 2, 3, 8, 9, 6, 6, 7, 8, 12, 7, 12, 13, 11, 0, 10, 10, 11/)
!  integer(i_kind), parameter, dimension(1:USGS_N_TYPES) :: usgs_to_gfs=(/7, &
!    12, 12, 12, 12, 12, 7, 9, 8, 6, 2, 5, 1, 4, 3, 0, 8, 8, 11, 10, 10, &
!    10, 11, 13/)
!  integer(i_kind), parameter, dimension(1:NAM_SOIL_N_TYPES) :: nmm_soil_to_crtm=(/1, &
!    1, 4, 2, 2, 8, 7, 2, 6, 5, 2, 3, 8, 1, 6, 9/)

!-------------------------------------------------------------------------------------------

!-from subroutine call_crtm in GSI/crtm_interface.f90
subroutine uv_to_wdir(uu5, vv5, wind10_direction)

   implicit none

   real (kind=kind_real), intent(in)  :: uu5, vv5
   real (kind=kind_real), intent(out) :: wind10_direction
   real (kind=kind_real)              :: windratio, windangle
   integer                            :: iquadrant  
   real(kind=kind_real),parameter:: windscale = 999999.0_kind_real
   real(kind=kind_real),parameter:: windlimit = 0.0001_kind_real
   real(kind=kind_real),parameter:: quadcof  (4, 2  ) = &
      reshape((/MPAS_JEDI_ZERO_kr,  MPAS_JEDI_ONE_kr,  MPAS_JEDI_ONE_kr,  MPAS_JEDI_TWO_kr, &
                MPAS_JEDI_ONE_kr,  -MPAS_JEDI_ONE_kr,  MPAS_JEDI_ONE_kr, -MPAS_JEDI_ONE_kr/), (/4, 2/))

   if (uu5 >= MPAS_JEDI_ZERO_kr .and. vv5 >= MPAS_JEDI_ZERO_kr) iquadrant = 1
   if (uu5 >= MPAS_JEDI_ZERO_kr .and. vv5 <  MPAS_JEDI_ZERO_kr) iquadrant = 2
   if (uu5 <  MPAS_JEDI_ZERO_kr .and. vv5 >= MPAS_JEDI_ZERO_kr) iquadrant = 4
   if (uu5 <  MPAS_JEDI_ZERO_kr .and. vv5 <  MPAS_JEDI_ZERO_kr) iquadrant = 3
   if (abs(vv5) >= windlimit) then
      windratio = uu5 / vv5
   else
      windratio = MPAS_JEDI_ZERO_kr
      if (abs(uu5) > windlimit) then
         windratio = windscale * uu5
      endif
   endif
   windangle        = atan(abs(windratio))   ! wind azimuth is in radians
   wind10_direction = ( quadcof(iquadrant, 1) * MPAS_JEDI_PII_kr + windangle * quadcof(iquadrant, 2) )

end subroutine uv_to_wdir

!-------------------------------------------------------------------------------------------
elemental subroutine w_to_q(mixing_ratio, specific_humidity)
   implicit none
   real (kind=kind_real), intent(in)  :: mixing_ratio
   real (kind=kind_real), intent(out) :: specific_humidity

   specific_humidity = mixing_ratio / (MPAS_JEDI_ONE_kr + mixing_ratio)
end subroutine w_to_q
!-------------------------------------------------------------------------------------------
elemental subroutine q_to_w(specific_humidity, mixing_ratio)
   implicit none
   real (kind=kind_real), intent(in)  :: specific_humidity
   real (kind=kind_real), intent(out) :: mixing_ratio

   mixing_ratio = specific_humidity / (MPAS_JEDI_ONE_kr - specific_humidity)
end subroutine q_to_w
!-------------------------------------------------------------------------------------------
elemental subroutine q_to_w_tl(specific_humidity_tl, sh_traj, mixing_ratio_tl)
   implicit none
   real (kind=kind_real), intent(in)  :: specific_humidity_tl
   real (kind=kind_real), intent(in)  :: sh_traj
   real (kind=kind_real), intent(out) :: mixing_ratio_tl

   mixing_ratio_tl = specific_humidity_tl / (MPAS_JEDI_ONE_kr - sh_traj)**2
end subroutine q_to_w_tl
!-------------------------------------------------------------------------------------------
elemental subroutine q_to_w_ad(specific_humidity_ad, sh_traj, mixing_ratio_ad)
   implicit none
   real (kind=kind_real), intent(inout) :: specific_humidity_ad
   real (kind=kind_real), intent(in)    :: sh_traj
   real (kind=kind_real), intent(in)    :: mixing_ratio_ad

   specific_humidity_ad = specific_humidity_ad + &
                 MPAS_JEDI_ONE_kr / ( MPAS_JEDI_ONE_kr - sh_traj)**2 * mixing_ratio_ad
end subroutine q_to_w_ad
!-------------------------------------------------------------------------------------------
elemental subroutine tw_to_tv(temperature,mixing_ratio,virtual_temperature)
   implicit none
   real (kind=kind_real), intent(in)  :: temperature
   real (kind=kind_real), intent(in)  :: mixing_ratio
   real (kind=kind_real), intent(out) :: virtual_temperature

   virtual_temperature = temperature * &
                  ( MPAS_JEDI_ONE_kr + (rv/rgas - MPAS_JEDI_ONE_kr)*mixing_ratio )
end subroutine tw_to_tv
!-------------------------------------------------------------------------------------------
elemental subroutine tw_to_tv_tl(temperature_tl,mixing_ratio_tl,t_traj,m_traj,virtual_temperature_tl)
   implicit none
   real (kind=kind_real), intent(in)  :: temperature_tl
   real (kind=kind_real), intent(in)  :: mixing_ratio_tl
   real (kind=kind_real), intent(in)  :: t_traj
   real (kind=kind_real), intent(in)  :: m_traj
   real (kind=kind_real), intent(out) :: virtual_temperature_tl

   virtual_temperature_tl = temperature_tl * &
                  ( MPAS_JEDI_ONE_kr + (rv/rgas - MPAS_JEDI_ONE_kr)*m_traj )  + &
                  t_traj * (rv/rgas - MPAS_JEDI_ONE_kr)*mixing_ratio_tl
end subroutine tw_to_tv_tl
!-------------------------------------------------------------------------------------------
elemental subroutine tw_to_tv_ad(temperature_ad,mixing_ratio_ad,t_traj,m_traj,virtual_temperature_ad)
   implicit none
   real (kind=kind_real), intent(inout) :: temperature_ad
   real (kind=kind_real), intent(inout) :: mixing_ratio_ad
   real (kind=kind_real), intent(in)    :: t_traj
   real (kind=kind_real), intent(in)    :: m_traj
   real (kind=kind_real), intent(in)    :: virtual_temperature_ad

   temperature_ad = temperature_ad + virtual_temperature_ad * &
                  ( MPAS_JEDI_ONE_kr + (rv/rgas - MPAS_JEDI_ONE_kr)*m_traj )
   mixing_ratio_ad = mixing_ratio_ad + virtual_temperature_ad * &
                  t_traj * (rv/rgas - MPAS_JEDI_ONE_kr)
end subroutine tw_to_tv_ad
!-------------------------------------------------------------------------------------------
elemental subroutine theta_to_temp(theta,pressure,temperature)
   implicit none
   real (kind=kind_real), intent(in)  :: theta
   real (kind=kind_real), intent(in)  :: pressure
   real (kind=kind_real), intent(out) :: temperature
   temperature = theta / &
             ( MPAS_JEDI_P0_kr / pressure ) ** ( rgas / cp )
   !TODO: Following formula would give the same result with the formular above,
   !      but gives a slightly different precision. Need to check.
   !temperature = theta * &
   !          ( pressure / MPAS_JEDI_P0_kr ) ** ( rgas / cp )
end subroutine theta_to_temp
!-------------------------------------------------------------------------------------------
elemental subroutine temp_to_theta(temperature,pressure,theta)
   implicit none
   real (kind=kind_real), intent(in)  :: temperature
   real (kind=kind_real), intent(in)  :: pressure
   real (kind=kind_real), intent(out) :: theta
   theta = temperature * &
             ( MPAS_JEDI_P0_kr / pressure ) ** ( rgas / cp )
end subroutine temp_to_theta
!-------------------------------------------------------------------------------------------
elemental subroutine twp_to_rho(temperature,mixing_ratio,pressure,rho)
   implicit none
   real (kind=kind_real), intent(in)  :: temperature
   real (kind=kind_real), intent(in)  :: mixing_ratio
   real (kind=kind_real), intent(in)  :: pressure
   real (kind=kind_real), intent(out) :: rho
   rho = pressure / ( rgas * temperature * &
                                ( MPAS_JEDI_ONE_kr + (rv/rgas) * mixing_ratio ) )
end subroutine twp_to_rho
!-------------------------------------------------------------------------------------------
subroutine pressure_half_to_full(pressure, zgrid, surface_pressure, nC, nV, pressure_f)
   implicit none
   real (kind=kind_real), dimension(nV,nC), intent(in) :: pressure
   real (kind=kind_real), dimension(nV+1,nC), intent(in) :: zgrid
   real (kind=kind_real), dimension(nC), intent(in) :: surface_pressure
   integer, intent(in) :: nC, nV
   real (kind=kind_real), dimension(nV+1,nC), intent(out) :: pressure_f

   real (kind=kind_real), dimension(nC,nV) :: fzm_p, fzp_p 
   real (kind=kind_real) :: tem1, z0, z1, z2, w1, w2
   integer :: i, k, its, ite, kts, kte

        !-- ~/libs/MPAS-Release/src/core_atmosphere/physics/mpas_atmphys_manager.F   >> dimension Line 644.
        !-- ~/libs/MPAS-Release/src/core_atmosphere/physics/mpas_atmphys_interface.F >> formula   Line 365.
        !-- ~/libs/MPAS-Release/src/core_atmosphere/physics/mpas_atmphys_vars.F      >> declarations
        ! This routine needs to access GEOM (dimension, zgrid )
        ! TODO: Check: ite = nCells??  nCellsSolve???  MPI consideration ???
        !              Seems to be nCellsSolve
        its=1 ; ite = nC
        kts=1 ; kte = nV

        do k = kts+1,kte
        do i = its,ite
          tem1 = MPAS_JEDI_ONE_kr/(zgrid(k+1,i)- zgrid(k-1,i))
          fzm_p(i,k) = ( zgrid(k,i)- zgrid(k-1,i)) * tem1
          fzp_p(i,k) = ( zgrid(k+1,i)- zgrid(k,i)) * tem1
          pressure_f(k,i) = fzm_p(i,k)*pressure(k,i) + fzp_p(i,k)*pressure(k-1,i)
        enddo
        enddo
        k = kte+1
        do i = its,ite
          z0 = zgrid(k,i)
          z1 = MPAS_JEDI_HALF_kr*(zgrid(k,i)+zgrid(k-1,i))
          z2 = MPAS_JEDI_HALF_kr*(zgrid(k-1,i)+zgrid(k-2,i))
          w1 = (z0-z2)/(z1-z2)
          w2 = MPAS_JEDI_ONE_kr-w1
          !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
          pressure_f(k,i) = exp( w1*log(pressure(k-1,i)) + w2*log(pressure(k-1,i)) )
        enddo
        k = kts
        do i = its,ite
          pressure_f(k,i) = surface_pressure(i)
        enddo

end subroutine pressure_half_to_full
!-------------------------------------------------------------------------------------------
subroutine hydrostatic_balance(ncells, nlevels, zw, t, qv, ps, p, rho, theta)
!---------------------------------------------------------------------------------------
! Purpose: 
!    Derive 3D pressure, dry air density, and dry potential temperature
!    from analyzed surface pressure, temperature, and water vapor mixing 
!    ratio by applying hypsometric equation, equation of state, and
!    Poisson's equation from bottom to top.
! Method:  
!    1. Vertical integral of pressure from half level k-1 to k is split 
!       into two steps: step-1 does integral from half level k-1 to full
!       level k, followed by a step-2 integral from full level k to half level k.
!       Full-level Tv is derived using bilinear interpolation of half level Tv.
!    2. Layer's virtual T: Tv = T * (1 + 0.608*Qv)
!    3. Hypsometric Eq.: P(z2) = P(z1) * exp[-g*(z2-z1)/(Rd*Tv)]
!    4. Eq. of State: rho_m = P/(Rd*Tv); rho_d = rho_m/(1+Qv). Add Qc/Qi/Qs/Qg?
!    5. Poisson's Eq.: theta_d = T * (P0/P)^(Rd/cp)
! References: 
!    MPAS-A model user's guide version 7.0, Appendix C.2 for vertical grid.
!--------------------------------------------------------------------------
   implicit none
   integer,                                            intent(in)  :: ncells, nlevels
   real (kind=kind_real), dimension(nlevels+1,ncells), intent(in)  :: zw    ! physical height m at w levels
   real (kind=kind_real), dimension(nlevels,ncells),   intent(in)  :: t     ! temperature, K
   real (kind=kind_real), dimension(nlevels,ncells),   intent(in)  :: qv    ! mixing ratio, kg/kg
   real (kind=kind_real), dimension(ncells),           intent(in)  :: ps    ! surface P, Pa
   real (kind=kind_real), dimension(nlevels,ncells),   intent(out) :: p     ! 3D P, Pa
   real (kind=kind_real), dimension(nlevels,ncells),   intent(out) :: rho   ! dry air density, kg/m^3
   real (kind=kind_real), dimension(nlevels,ncells),   intent(out) :: theta ! dry potential T, K
   
   integer                :: icell, k
   real (kind=kind_real)  :: rvordm1   ! rv/rd - 1. = 461.6/287-1 ~= 0.60836
   real (kind=kind_real), dimension(nlevels)   :: tv_h  ! half level virtual T
   real (kind=kind_real), dimension(nlevels+1) :: pf    ! full level pressure
   real (kind=kind_real), dimension(nlevels)   :: zu    ! physical height at u level
   real (kind=kind_real)  :: tv_f, tv, w
   
      rvordm1 = rv/rgas - MPAS_JEDI_ONE_kr

      do icell=1,ncells

         k = 1 ! 1st half level
         pf(k) = ps(icell) ! first full level P is at surface
         zu(k) = MPAS_JEDI_HALF_kr * ( zw(k,icell) + zw(k+1,icell) )
         tv_h(k) = t(k,icell) * ( MPAS_JEDI_ONE_kr + rvordm1*qv(k,icell) )

       ! integral from full level k to half level k
         p(k,icell) = pf(k) * exp( -gravity * (zu(k)-zw(k,icell))/(rgas*tv_h(k)) )
         rho(k,icell) = p(k,icell)/( rgas*tv_h(k)*(MPAS_JEDI_ONE_kr+qv(k,icell)) )
         theta(k,icell) = t(k,icell) * ( MPAS_JEDI_P0_kr/p(k,icell) )**(rgas/cp)

       do k=2,nlevels ! loop over half levels

       ! half level physical height, zu(k-1) -> zw(k) -> zu(k)
         zu(k) = MPAS_JEDI_HALF_kr * ( zw(k,icell) + zw(k+1,icell) )
   
       ! bilinear interpolation weight for value at the half level below the full level
         w = ( zu(k) - zw(k,icell) )/( zu(k) - zu(k-1) )

       ! half level Tv
         tv_h(k)   = t(k,  icell) * ( MPAS_JEDI_ONE_kr + rvordm1*qv(k,  icell) )

       ! interpolate two half levels Tv to a full level Tv
         tv_f = w * tv_h(k-1) + (MPAS_JEDI_ONE_kr - w) * tv_h(k)

       ! two-step integral from half level k-1 to k
       !-----------------------------------------------------------

       ! step-1: tv used for integral from half level k-1 to full level k
         tv = MPAS_JEDI_HALF_kr * ( tv_h(k-1) + tv_f )

       ! step-1: integral from half level k-1 to full level k, hypsometric Eq.
         pf(k) = p(k-1,icell) * exp( -gravity * (zw(k,icell)-zu(k-1))/(rgas*tv) )

       ! step-2: tv used for integral from full level k to half level k
         tv = MPAS_JEDI_HALF_kr * ( tv_h(k) + tv_f )

       ! step-2: integral from full level k to half level k
         p(k,icell) = pf(k) * exp( -gravity * (zu(k)-zw(k,icell))/(rgas*tv) )

       !------------------------------------------------------------
       ! dry air density at half level, equation of state
         rho(k,icell) = p(k,icell)/( rgas*tv_h(k)*(MPAS_JEDI_ONE_kr+qv(k,icell)) )

       ! dry potential T at half level, Poisson's equation
         theta(k,icell) = t(k,icell) * ( MPAS_JEDI_P0_kr/p(k,icell) )**(rgas/cp)

       end do
      end do

end subroutine hydrostatic_balance
!-------------------------------------------------------------------------------------------
subroutine geometricZ_full_to_half(zgrid_f, nC, nV, zgrid)
   implicit none
   real (kind=kind_real), dimension(nV+1,nC), intent(in) :: zgrid_f
   integer, intent(in) :: nC, nV
   real (kind=kind_real), dimension(nV,nC), intent(out) :: zgrid
   integer :: i, k

!  calculate midpoint geometricZ:
   do i=1,nC
      do k=1,nV
         zgrid(k,i) = ( zgrid_f(k,i) + zgrid_f(k+1,i) ) * MPAS_JEDI_HALF_kr
      enddo
   enddo
end subroutine geometricZ_full_to_half
!-------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------
subroutine index_q_fields_forward(indexName, geovalName, subFields, convFields, pressure_levels, ngrid, nVertLevels)

   implicit none

   character (len=*),                     intent(in)    :: indexName
   character (len=*),                     intent(in)    :: geovalName
   type (mpas_pool_type), pointer,        intent(in)    :: subFields    !< self % subFields
   type (mpas_pool_type), pointer,        intent(inout) :: convFields   !< pool with geovals variable
   real (kind=kind_real), dimension(:,:), intent(in)    :: pressure_levels
   integer,                               intent(in)    :: ngrid        !< number of grid cells
   integer,                               intent(in)    :: nVertLevels

   type (field2DReal), pointer :: converted_field, index_field_src
   real (kind=kind_real) :: kgkg_kgm2
   integer :: i, k

   call mpas_pool_get_field(subFields, indexName, index_field_src) !- [kg/kg]
   call mpas_duplicate_field(index_field_src, converted_field)
   do i=1,ngrid
      do k=1,nVertLevels
         kgkg_kgm2=( pressure_levels(k,i)-pressure_levels(k+1,i) ) / gravity !- Still bottom-to-top
         converted_field % array(k,i) = index_field_src%array(k,i) * kgkg_kgm2
      enddo
   enddo
   
   ! Ensure positive-definite hydrometeor mixing ratios
   !  with respect to precision of crtm::CRTM_Parameters::ZERO.
   !  Note: replacing MPAS_JEDI_GREATERZERO_kr with MPAS_JEDI_ZERO_kr
   !   will cause some profiles to fail in CRTM.
   where(converted_field % array(:,1:ngrid) < MPAS_JEDI_GREATERZERO_kr)
      converted_field % array(:,1:ngrid) = MPAS_JEDI_GREATERZERO_kr
   end where

   converted_field % fieldName = geovalName
   call mpas_pool_add_field(convFields, geovalName, converted_field)
   !        write(*,*) "end-of ",geovalName
end subroutine index_q_fields_forward
!-------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------
subroutine index_q_fields_TL(indexName, geovalName, subFields, convFields, pressure_levels, ngrid, nVertLevels)

   implicit none

   character (len=*),                     intent(in)    :: indexName
   character (len=*),                     intent(in)    :: geovalName
   type (mpas_pool_type), pointer,        intent(in)    :: subFields    !< self % subFields
   type (mpas_pool_type), pointer,        intent(inout) :: convFields   !< pool with geovals variable
   real (kind=kind_real), dimension(:,:), intent(in)    :: pressure_levels
   integer,                               intent(in)    :: ngrid        !< number of grid cells
   integer,                               intent(in)    :: nVertLevels

   type (field2DReal), pointer :: converted_field, index_field_src
   real (kind=kind_real) :: kgkg_kgm2
   integer :: i, k

   call mpas_pool_get_field(subFields, indexName, index_field_src) !- [kg/kg]
   call mpas_duplicate_field(index_field_src, converted_field)
   do i=1,ngrid
      do k=1,nVertLevels
         kgkg_kgm2=( pressure_levels(k,i)-pressure_levels(k+1,i) ) / gravity !- Still bottom-to-top
         converted_field % array(k,i) = index_field_src%array(k,i) * kgkg_kgm2
      enddo
   enddo

   converted_field % fieldName = geovalName
   call mpas_pool_add_field(convFields, geovalName, converted_field)
   !        write(*,*) "end-of ",geovalName
end subroutine index_q_fields_TL
!-------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------
subroutine index_q_fields_AD(geom, trajFields, indexName, geovalName, subFields, convFields, ngrid)

   implicit none

   type (mpas_geom),               intent(in)    :: geom          !< geometry
   type (mpas_pool_type), pointer, intent(in)    :: trajFields    !< linearization state for conversion
   character (len=*),              intent(in)    :: indexName
   character (len=*),              intent(in)    :: geovalName
   type (mpas_pool_type), pointer, intent(inout) :: subFields    !< self % subFields
   type (mpas_pool_type), pointer, intent(in)    :: convFields   !< pool with geovals variable
   integer,                        intent(in )   :: ngrid        !< number of grid cells
   
   real (kind=kind_real), dimension(:,:), allocatable :: pressure_levels
   real (kind=kind_real), dimension(:,:), pointer :: index_array
   type (field2DReal), pointer :: index_increment, geoval_field_src
   real (kind=kind_real), dimension(:,:), pointer :: pressure_array
   real (kind=kind_real), dimension(:), pointer :: surface_pressure_array
   real (kind=kind_real) :: kgkg_kgm2
   integer :: i, k

   ! Calcluate air_pressure_levels
   call mpas_pool_get_array(trajFields, 'pressure', pressure_array)
   call mpas_pool_get_array(trajFields, 'surface_pressure', surface_pressure_array)
   allocate (pressure_levels(geom%nVertLevels+1,ngrid))
   call pressure_half_to_full(pressure_array(:,1:ngrid), geom%zgrid(:,1:ngrid), &
                              surface_pressure_array(1:ngrid), ngrid, geom%nVertLevels, &
                              pressure_levels)

   call mpas_pool_get_array(subFields, indexName, index_array) !- [kg/kg]
   call mpas_pool_get_field(convFields, geovalName, geoval_field_src)
   call mpas_duplicate_field(geoval_field_src, index_increment)
   do i=1,ngrid
      do k=1,geom%nVertLevels
         kgkg_kgm2=( pressure_levels(k,i)-pressure_levels(k+1,i) ) / gravity !- Still bottom-to-top
         index_increment % array(k,i) = geoval_field_src%array(k,i) * kgkg_kgm2
      enddo
   enddo
   index_array(:,1:ngrid) = index_array(:,1:ngrid) + index_increment % array(:,1:ngrid)
   call mpas_deallocate_field(index_increment) ! not used
   deallocate (pressure_levels)

end subroutine index_q_fields_AD
!-------------------------------------------------------------------------------------------
real (kind=kind_real) function wgamma(y)
implicit none
real (kind=kind_real), intent(in) :: y

wgamma = exp(gammln(y))

end function wgamma
!-------------------------------------------------------------------------------------------
real (kind=kind_real) function gammln(xx)
implicit none
real (kind=kind_real), intent(in)  :: xx
real (kind=kind_double), parameter :: stp = 2.5066282746310005_kind_double
real (kind=kind_double), parameter :: &
   cof(6) = (/   76.18009172947146_kind_double,    -86.50532032941677_kind_double, &
                 24.01409824083091_kind_double,    -1.231739572450155_kind_double, &
              0.001208650973866179_kind_double, -0.000005395239384953_kind_double/)
real (kind=kind_double) :: ser,tmp,x,y
integer :: j

x=real(xx,kind_double)
y=x
tmp=x+5.5_kind_double
tmp=(x+0.5_kind_double)*log(tmp)-tmp
ser=1.000000000190015_kind_double
do j=1,6
   y=y+1.0_kind_double
   ser=ser+cof(j)/y
end do
gammln=real(tmp+log(stp*ser/x),kind_real)
end function gammln

!-------------------------------------------------------------------------------------------
subroutine effectRad_rainwater (qr, rho, nr, re_qr, mp_scheme, ngrid, nVertLevels)
!-----------------------------------------------------------------------
!  Compute radiation effective radii of rainwater for WSM6/Thompson microphysics.
!  These are consistent with microphysics equations.
!  References: WRF/phys/module_mp_wsm6.F
!              Lin, Y. L., Farley, R. D., & Orville, H. D. (1983). Bulk parameterization of the snow field in a cloud model. Journal of climate and applied meteorology, 22(6), 1065-1092.
!              WRF/phys/module_mp_thompson.F
!              Thompson, G., Rasmussen, R. M., & Manning, K. (2004). Explicit forecasts of winter precipitation using an improved bulk microphysics scheme. Part I: Description and sensitivity analysis. Monthly Weather Review, 132(2), 519-542.
!              Thompson, G., Field, P. R., Rasmussen, R. M., & Hall, W. D. (2008). Explicit forecasts of winter precipitation using an improved bulk microphysics scheme. Part II: Implementation of a new snow parameterization. Monthly Weather Review, 136(12), 5095-5115.
!-----------------------------------------------------------------------

implicit none

real(kind=kind_real), dimension( nVertLevels, ngrid ), intent(in)  :: qr, rho
real(kind=kind_real), dimension( nVertLevels, ngrid ), intent(in)  :: nr
integer,                                               intent(in)  :: ngrid, nVertLevels
character(len=StrKIND),                                intent(in)  :: mp_scheme
real(kind=kind_real), dimension( nVertLevels, ngrid ), intent(out) :: re_qr

!Local variables
! constants
real(kind=kind_real), parameter :: denr = MPAS_JEDI_THOUSAND_kr, n0r = 8.e6
real(kind=kind_real), parameter :: R1 = MPAS_JEDI_ONE_kr / &
                                        MPAS_JEDI_MILLION_kr / &
                                        MPAS_JEDI_MILLION_kr, &
                                   R2 = MPAS_JEDI_ONE_kr / &
                                        MPAS_JEDI_MILLION_kr
real(kind=kind_real), parameter :: mu_r = MPAS_JEDI_ZERO_kr
real(kind=kind_real), parameter :: am_r = MPAS_JEDI_PII_kr*denr/6.0_kind_real
real(kind=kind_real), parameter :: bm_r = MPAS_JEDI_THREE_kr

real(kind=kind_double) :: lamdar
integer                :: i, k

real(kind=kind_real), dimension( nVertLevels, ngrid ) :: rqr, nr_rho
!For Thompson scheme
real(kind=kind_real) :: cre2,cre3,crg2,crg3,org2,ombr,obmr
!Generalized gamma distributions for rain
! N(D) = N_0 * D**mu * exp(-lamda*D);  mu=0 is exponential.

!-----------------------------------------------------------------------
cre2 = mu_r + MPAS_JEDI_ONE_kr
cre3 = bm_r + mu_r + MPAS_JEDI_ONE_kr
crg2 = wgamma(cre2)
crg3 = wgamma(cre3)
org2 = MPAS_JEDI_ONE_kr/crg2
obmr = MPAS_JEDI_ONE_kr/bm_r

do i = 1, ngrid
   do k = 1, nVertLevels
      rqr(k,i) = max(R1, qr(k,i)*rho(k,i))
   enddo
enddo

if (any(rqr > R1)) then
   do i = 1, ngrid
      do k = 1, nVertLevels
         re_qr(k,i) = 99.e-6
         if (rqr(k,i).le.R1) CYCLE
         select case (trim(mp_scheme))
         case ('mp_wsm6')
            lamdar = sqrt(sqrt(MPAS_JEDI_PII_kr*denr*n0r/rqr(k,i)))
            re_qr(k,i) =  max(99.9D-6,min(1.5_kind_double/lamdar,1999.D-6))
         case ('mp_thompson')
            nr_rho(k,i) = max(R2, nr(k,i)*rho(k,i))
            lamdar = (am_r*crg3*org2*nr_rho(k,i)/rqr(k,i))**obmr
            re_qr(k,i) = max(99.9e-6, min(real(MPAS_JEDI_HALF_kr*(MPAS_JEDI_THREE_kr+mu_r),kind_double)/lamdar, 1999.e-6))
         case default
            re_qr(k,i) = 999.e-6
         end select
      enddo
   enddo
endif

end subroutine effectRad_rainwater

!-------------------------------------------------------------------------------------------
subroutine effectRad_graupel (qg, rho, re_qg, mp_scheme, ngrid, nVertLevels)

!-----------------------------------------------------------------------
!  Compute radiation effective radii of graupel for WSM6/Thompson microphysics.
!  These are consistent with microphysics equations.
!  References: WRF/phys/module_mp_wsm6.F
!              Lin, Y. L., Farley, R. D., & Orville, H. D. (1983). Bulk parameterization of the snow field in a cloud model. Journal of climate and applied meteorology, 22(6), 1065-1092.
!              WRF/phys/module_mp_thompson.F
!              Thompson, G., Rasmussen, R. M., & Manning, K. (2004). Explicit forecasts of winter precipitation using an improved bulk microphysics scheme. Part I: Description and sensitivity analysis. Monthly Weather Review, 132(2), 519-542.
!              Thompson, G., Field, P. R., Rasmussen, R. M., & Hall, W. D. (2008). Explicit forecasts of winter precipitation using an improved bulk microphysics scheme. Part II: Implementation of a new snow parameterization. Monthly Weather Review, 136(12), 5095-5115.
!-----------------------------------------------------------------------

implicit none

real(kind=kind_real), dimension( nVertLevels, ngrid ), intent(in) :: qg, rho
integer,                                               intent(in) :: ngrid, nVertLevels
character (len=StrKIND),                               intent(in) :: mp_scheme
real(kind=kind_real), dimension( nVertLevels, ngrid ), intent(out):: re_qg

!Local variables
integer                :: i, k, k_0
real(kind=kind_double) :: lamdag, lam_exp, N0_exp
real(kind=kind_real)   :: n0g, deng
real(kind=kind_real), parameter :: R1 = MPAS_JEDI_ONE_kr / &
                                        MPAS_JEDI_MILLION_kr / &
                                        MPAS_JEDI_MILLION_kr
real(kind=kind_real), dimension( nVertLevels, ngrid ):: rqg
!MPAS model set it as 0, WRF model set it through namelist
integer:: hail_opt = 0
! for Thompson scheme
real(kind=kind_real) :: mu_g = MPAS_JEDI_ZERO_kr
real(kind=kind_real) :: obmg, cge1,cgg1,oge1,cge3,cgg3,ogg1,cge2,cgg2,ogg2
real(kind=kind_real) :: xslw1, ygra1, zans1
real(kind=kind_real), parameter :: am_g = MPAS_JEDI_PII_kr*500.0_kind_real/6.0_kind_real
real(kind=kind_real), parameter :: bm_g = MPAS_JEDI_THREE_kr

!-----------------------------------------------------------------------
obmg = MPAS_JEDI_ONE_kr/bm_g
cge1 = bm_g + MPAS_JEDI_ONE_kr
cgg1 = wgamma(cge1)
oge1 = MPAS_JEDI_ONE_kr/cge1
cge3 = bm_g + mu_g + MPAS_JEDI_ONE_kr
cgg3 = wgamma(cge3)
ogg1 = MPAS_JEDI_ONE_kr/cgg1
cge2 = mu_g + MPAS_JEDI_ONE_kr
cgg2 = wgamma(cge2)
ogg2 = MPAS_JEDI_ONE_kr/cgg2

if (hail_opt .eq. 1) then
   n0g  = 4.e4
   deng = 700.0_kind_real
else
   n0g  = 4.e6
   deng = 500.0_kind_real
endif

do i = 1, ngrid
   do k = 1, nVertLevels
      rqg(k,i) = max(R1, qg(k,i)*rho(k,i))
   enddo
enddo

if (any( rqg > R1 )) then
   select case (trim(mp_scheme))
   case ('mp_wsm6')
      re_qg = 49.7e-6
      do i = 1, ngrid
         do k = 1, nVertLevels
            if (rqg(k,i).le.R1) CYCLE
            lamdag = sqrt(sqrt(MPAS_JEDI_PII_kr*deng*n0g/rqg(k,i)))
            re_qg(k,i) = max(50.D-6,min(1.5_kind_double/lamdag,9999.D-6))
         end do
      end do
   case ('mp_thompson')
      re_qg = 99.5e-6
      do i = 1, ngrid
         do k = nVertLevels, 1, -1
            if (rqg(k,i).le.R1) CYCLE
            ygra1 = alog10(sngl(max(1.e-9, rqg(k,i))))
            zans1 = (2.5_kind_real + 2.5_kind_real/7.0_kind_real * (ygra1+7.0_kind_real))
            zans1 = max(MPAS_JEDI_TWO_kr, min(zans1, 7.0_kind_real)) ! new in WRF V4.2
            N0_exp = 10.0_kind_real**(zans1)
            lam_exp = (N0_exp*am_g*cgg1/rqg(k,i))**oge1
            lamdag = lam_exp * (cgg3*ogg2*ogg1)**obmg ! we can simplify this without considering rainwater
            re_qg(k,i) = max(99.9e-6, min(real(MPAS_JEDI_HALF_kr*(MPAS_JEDI_THREE_kr+mu_g),kind_double)/lamdag, 9999.e-6))
         enddo
      enddo
   case default
      re_qg = 600.e-6
   end select
endif

end subroutine effectRad_graupel
!-------------------------------------------------------------------------------------------

!--- variables can be found in subFields
subroutine convert_mpas_field2ufo(geom, subFields, convFields, fieldname, nfield, ngrid)

   implicit none

   type(mpas_geom),                intent(in)  :: geom         !< geometry
   type (mpas_pool_type), pointer, intent(in ) :: subFields    !< self % subFields
   type (mpas_pool_type), pointer, intent(out) :: convFields   !< pool with geovals variable
   character (len=*),              intent(in ) :: fieldname(:) !< list of variables for geovals
   integer,                        intent(in ) :: nfield       !< number of variables
   integer,                        intent(in ) :: ngrid        !< number of grid cells

   type (mpas_pool_iterator_type) :: poolItr
   type (mpas_pool_type), pointer :: allFields
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b, r2d_ptr_c
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   integer, dimension(:), pointer :: i1d_ptr_a, i1d_ptr_b

   type (field1DInteger), pointer :: field1di, field1di_src
   type (field1DReal), pointer :: field1d, field1d_src
   type (field2DReal), pointer :: field2d, field2d_a, field2d_src
   integer :: ivar, i, k

   real (kind=kind_real) :: lat
   type (field2DReal), pointer :: field2d_nr, field2d_qr, field2d_qg, field2d_rho
   character(len=StrKIND),  pointer :: config_microp_scheme, &
                                       config_radt_cld_scheme
   logical,  pointer :: config_microp_re

   !--- create new pool for geovals
   call mpas_pool_create_pool(convFields)

   do ivar=1, nfield
!     write(*,*) 'convert_mpas_field2ufo  :inside do/select case, &
!               & ivar, trim(fieldname(ivar))=',ivar,trim(fieldname(ivar))

     select case (trim(fieldname(ivar)))

     case ( var_tv ) !-virtual_temperature
        call mpas_pool_get_field(subFields, 'temperature', field2d_src) !< get temperature
        call mpas_pool_get_array(subFields,     'spechum', r2d_ptr_b)   !< get specific_humidity
!        write(*,*) 'MIN/MAX of temperature=',minval(field2d_src%array),maxval(field2d_src%array)
!        write(*,*) 'MIN/MAX of     spechum=',minval(r2d_ptr_b),maxval(r2d_ptr_b)

        call mpas_duplicate_field(field2d_src, field2d)   ! for virtual_temperature
        call mpas_duplicate_field(field2d_src, field2d_a) ! for mixing_ratio, intermediate variable

        call q_to_w( r2d_ptr_b(:,1:ngrid) , field2d_a % array(:,1:ngrid) )
        call tw_to_tv( field2d_src%array(:,1:ngrid), field2d_a % array(:,1:ngrid), field2d % array(:,1:ngrid) )
        call mpas_deallocate_field(field2d_a) ! not used
!        write(*,*) 'MIN/MAX of Tv=',minval(field2d % array(:,1:ngrid)),maxval(field2d % array(:,1:ngrid))

        field2d % fieldName = var_tv

        call mpas_pool_add_field(convFields, var_tv, field2d)

!        write(*,*) "end-of ",var_tv

     case ( var_ts, var_t ) !-air_temperature, temperature
        call mpas_pool_get_field(subFields, 'temperature', field2d_src) !< get temperature
        call mpas_duplicate_field(field2d_src, field2d)! for air_temperature
        field2d % fieldName = trim(fieldname(ivar))
        call mpas_pool_add_field(convFields, trim(fieldname(ivar)), field2d)

     case ( var_u ) !-eastward_wind
        call mpas_pool_get_field(subFields, 'uReconstructZonal', field2d_src) !< get zonal wind
        call mpas_duplicate_field(field2d_src, field2d)!  for eastward_wind
!        write(*,*) 'MIN/MAX of zonal wind=',minval(field2d % array),maxval(field2d % array)
        field2d % fieldName = var_u
        call mpas_pool_add_field(convFields, var_u, field2d)
!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ( var_v ) !-northward_wind
        call mpas_pool_get_field(subFields, 'uReconstructMeridional', field2d_src) !< get meridional wind
        call mpas_duplicate_field(field2d_src, field2d)!  for northward_wind
!        write(*,*) 'MIN/MAX of meridional wind=',minval(field2d % array),maxval(field2d % array)
        field2d % fieldName = var_v
        call mpas_pool_add_field(convFields, var_v, field2d)
!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ( var_mixr ) !-humidity_mixing_ratio
        call mpas_pool_get_field(subFields, 'spechum', field2d_src) !< get specific_humidity
        call mpas_duplicate_field(field2d_src, field2d)! for humidity_mixing_ratio
        call q_to_w(field2d_src % array(:,1:ngrid) , field2d % array(:,1:ngrid))
        field2d % array(:,1:ngrid) = max(MPAS_JEDI_ZERO_kr, field2d % array(:,1:ngrid)) * MPAS_JEDI_THOUSAND_kr ! [kg/kg] -> [g/kg]
        field2d % fieldName = var_mixr
        call mpas_pool_add_field(convFields, var_mixr, field2d)
!        write(*,*) "end-of ",var_mixr

     case( var_q ) !-specific_humidity
        call mpas_pool_get_field(subFields, "spechum", field2d_src)
        call mpas_duplicate_field(field2d_src, field2d)
        field2d % fieldName = var_q
        call mpas_pool_add_field(convFields, var_q, field2d)
!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ( var_prs ) !-air_pressure
        call mpas_pool_get_field(subFields, 'pressure', field2d_src) !< get pressure
        call mpas_duplicate_field(field2d_src, field2d)
        field2d % array(:,1:ngrid) = field2d_src%array(:,1:ngrid)
        field2d % fieldName = var_prs
        call mpas_pool_add_field(convFields, var_prs, field2d)
!        write(*,*) "end-of ",var_prs

     case ( var_ps ) !-surface_pressure
        call mpas_pool_get_field(subFields, 'surface_pressure', field1d_src) !< get surface_pressure
        call mpas_duplicate_field(field1d_src, field1d)
        field1d % array (1:ngrid) = field1d_src%array(1:ngrid)
        field1d % fieldName = trim(fieldname(ivar))
        call mpas_pool_add_field(convFields, trim(fieldname(ivar)), field1d)

     case ( var_prsi ) !-air_pressure_levels
        call mpas_pool_get_array(subFields, "pressure", r2d_ptr_a)
        call mpas_pool_get_array(subFields, "surface_pressure", r1d_ptr_a)
        call mpas_pool_get_field(subFields, 'w', field2d_src) ! as a dummy array
        call mpas_duplicate_field(field2d_src, field2d)

        call pressure_half_to_full(r2d_ptr_a(:,1:ngrid), geom%zgrid(:,1:ngrid), r1d_ptr_a(1:ngrid), &
                                   ngrid, geom % nVertLevels, field2d%array(:,1:ngrid))
!        write(*,*) 'MIN/MAX of prsi=',minval(field2d % array),maxval(field2d % array)
!        write(*,*) 'test prs       =',r2d_ptr_a(:,1)
!        write(*,*) 'test prsi      =',field2d % array(:,1)

        field2d % fieldName = var_prsi
        call mpas_pool_add_field(convFields, var_prsi, field2d)
!        write(*,*) "end-of ",var_prsi

     case ( var_oz ) !-mole_fraction_of_ozone_in_air :TODO: not directly available from MPAS
!        call mpas_pool_get_array(subFields, "o3", r2d_ptr_a)
        call mpas_pool_get_field(subFields, 'theta', field2d_src) ! as a dummy array
        call mpas_duplicate_field(field2d_src, field2d)
        field2d % array(:,1:ngrid) = MPAS_JEDI_ZERO_kr !r2d_ptr_a(:,1:ngrid) ! convert ??
        field2d % fieldName = var_oz
        call mpas_pool_add_field(convFields, var_oz, field2d)
!        write(*,*) "end-of ",var_oz

     case ( var_co2 ) !-mole_fraction_of_carbon_dioxide_in_air :TODO: not directly available from MPAS
!        call mpas_pool_get_array(subFields, "co2", r2d_ptr_a)
        call mpas_pool_get_field(subFields, 'theta', field2d_src) ! as a dummy array
        call mpas_duplicate_field(field2d_src, field2d)
        field2d % array(:,1:ngrid) = MPAS_JEDI_ZERO_kr !r2d_ptr_a(:,1:ngrid) ! convert ??
        field2d % fieldName = var_co2
        call mpas_pool_add_field(convFields, var_co2, field2d)
!        write(*,*) "end-of ",var_co2

     case ( var_clw ) !-mass_content_of_cloud_liquid_water_in_atmosphere_layer
      !--TODO: Trial: Should already have "var_prsi"
      call mpas_pool_get_array(convFields, "air_pressure_levels", r2d_ptr_b) !- [hPa]
      call index_q_fields_forward('index_qc', var_clw, subFields, convFields, r2d_ptr_b, ngrid, geom % nVertLevels)

     case ( var_cli ) !-mass_content_of_cloud_ice_in_atmosphere_layer
         !--TODO: Trial: Should already have "var_prsi"
         call mpas_pool_get_array(convFields, "air_pressure_levels", r2d_ptr_b) !- [hPa]
         call index_q_fields_forward('index_qi', var_cli, subFields, convFields, r2d_ptr_b, ngrid, geom % nVertLevels)

      case ( var_clr ) !-mass_content_of_rain_in_atmosphere_layer
         !--TODO: Trial: Should already have "var_prsi"
         call mpas_pool_get_array(convFields, "air_pressure_levels", r2d_ptr_b) !- [hPa]
         call index_q_fields_forward('index_qr', var_clr, subFields, convFields, r2d_ptr_b, ngrid, geom % nVertLevels)

      case ( var_cls ) !-mass_content_of_snow_in_atmosphere_layer
         !--TODO: Trial: Should already have "var_prsi"
         call mpas_pool_get_array(convFields, "air_pressure_levels", r2d_ptr_b) !- [hPa]
         call index_q_fields_forward('index_qs', var_cls, subFields, convFields, r2d_ptr_b, ngrid, geom % nVertLevels)

      case ( var_clg ) !-mass_content_of_graupel_in_atmosphere_layer
         !--TODO: Trial: Should already have "var_prsi"
         call mpas_pool_get_array(convFields, "air_pressure_levels", r2d_ptr_b) !- [hPa]
         call index_q_fields_forward('index_qg', var_clg, subFields, convFields, r2d_ptr_b, ngrid, geom % nVertLevels)

      case ( var_clh ) !-mass_content_of_hail_in_atmosphere_layer
         !--TODO: Trial: Should already have "var_prsi"
         call mpas_pool_get_array(convFields, "air_pressure_levels", r2d_ptr_b) !- [hPa]
         call index_q_fields_forward('index_qh', var_clh, subFields, convFields, r2d_ptr_b, ngrid, geom % nVertLevels)

     case ( var_clwefr ) !-effective_radius_of_cloud water particle
        call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_re', config_microp_re)
        if (config_microp_re) then
           call mpas_pool_get_field(subFields, 're_cloud', field2d_src) !- [m]
           call mpas_duplicate_field(field2d_src, field2d)
           field2d % array(:,1:ngrid) = field2d_src%array(:,1:ngrid) * MPAS_JEDI_MILLION_kr ! [m] -> [micron]
        else
           call mpas_pool_get_field(subFields, 'index_qc', field2d_src) ! We do not really need 're_cloud' unless it is calculated in MPAS model microp
           call mpas_duplicate_field(field2d_src, field2d)
           field2d % array(:,1:ngrid) = 10.0_kind_real ![micron]  ! Or, we can refer to rre_cloud which is default for MPAS radiance scheme
        end if
        field2d % fieldName = var_clwefr
        call mpas_pool_add_field(convFields, var_clwefr, field2d)
!        write(*,*) "end-of ",var_clwefr

     case ( var_cliefr ) !-effective_radius_of_cloud ice particle
        call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_re', config_microp_re)
        if (config_microp_re) then
           call mpas_pool_get_field(subFields, 're_ice', field2d_src) !- [m]
           call mpas_duplicate_field(field2d_src, field2d)
           field2d % array(:,1:ngrid) = field2d_src%array(:,1:ngrid) * MPAS_JEDI_MILLION_kr ! [m] -> [micron]
        else
           call mpas_pool_get_field(subFields, 'index_qi', field2d_src) ! We do not really need 're_ice' unless it is calculated in MPAS model microp
           call mpas_duplicate_field(field2d_src, field2d)
           field2d % array(:,1:ngrid) = 30.0_kind_real ! [micron]
        end if
        field2d % fieldName = var_cliefr
        call mpas_pool_add_field(convFields, var_cliefr, field2d)
!        write(*,*) "end-of ",var_cliefr

      case ( var_clrefr )  !-effective_radius_of_rain water_particle
         call mpas_pool_get_field(subFields, 'index_qr', field2d_qr) ! effective_radius of rain_water is not calculated in MPAS model physics
         call mpas_duplicate_field(field2d_qr, field2d)
         call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_scheme', config_microp_scheme)
         call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_re', config_microp_re)
         if (config_microp_re) then
            call mpas_duplicate_field(field2d_qr, field2d_a)
            if (trim(config_microp_scheme) == 'mp_thompson') then
               call mpas_pool_get_field(subFields, 'index_nr', field2d_nr) !- [nb kg^{-1}]: MPAS output for 2-moment MP scheme
            else
               call mpas_duplicate_field(field2d_qr, field2d_nr)
               field2d_nr % array(:,1:ngrid) = MPAS_JEDI_ONE_kr
            end if
            call mpas_pool_get_field(subFields, 'rho', field2d_rho) !- [kg m^{-3}]: Dry air density
            call effectRad_rainwater(field2d_qr%array(:,1:ngrid), field2d_rho%array(:,1:ngrid),&
                                     field2d_nr%array(:,1:ngrid), field2d_a%array(:,1:ngrid), config_microp_scheme, &
                                     ngrid, geom%nVertLevels)
            field2d % array(:,1:ngrid) = field2d_a % array(:,1:ngrid) * MPAS_JEDI_MILLION_kr ! [m] -> [micron]
         else
            field2d % array(:,1:ngrid) = 999.0_kind_real ! [micron]
         end if
         field2d % fieldName = var_clrefr
         call mpas_pool_add_field(convFields, var_clrefr, field2d)
 !        write(*,*) "end-of ",var_clrefr

      case ( var_clsefr ) !-effective_radius_of_snow particle
         call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_re', config_microp_re)
         if (config_microp_re) then
            call mpas_pool_get_field(subFields, 're_snow', field2d_src) !- [m]
            call mpas_duplicate_field(field2d_src, field2d)
            field2d % array(:,1:ngrid) = field2d_src%array(:,1:ngrid) * MPAS_JEDI_MILLION_kr ! [m] -> [micron]
         else
            call mpas_pool_get_field(subFields, 'index_qs', field2d_src) ! We do not really need 're_snow' unless it is calculated in MPAS model microp
            call mpas_duplicate_field(field2d_src, field2d)
            field2d % array(:,1:ngrid) = 600.0_kind_real ! [micron]
         end if
         field2d % fieldName = var_clsefr
         call mpas_pool_add_field(convFields, var_clsefr, field2d)
 !        write(*,*) "end-of ",var_clsefr

      case ( var_clgefr ) !-effective_radius_of_graupel_particle
         call mpas_pool_get_field(subFields, 'index_qg', field2d_qg) !effective_radius_of_graupel is not calculated in MPAS model physics
         call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_scheme', config_microp_scheme)
         call mpas_pool_get_config(geom % domain % blocklist % configs, 'config_microp_re', config_microp_re)
         call mpas_duplicate_field(field2d_qg, field2d)
         if (config_microp_re) then
            call mpas_duplicate_field(field2d_qg, field2d_a)
            call mpas_pool_get_field(subFields, 'rho', field2d_rho) !- [kg m^{-3}]: Dry air density
            call effectRad_graupel(field2d_qg%array(:,1:ngrid), field2d_rho%array(:,1:ngrid), &
                                   field2d_a% array(:,1:ngrid), config_microp_scheme,         &
                                   ngrid, geom%nVertLevels)
            field2d % array(:,1:ngrid) = field2d_a % array(:,1:ngrid) * MPAS_JEDI_MILLION_kr ! [m] -> [micron]
         else
            field2d % array(:,1:ngrid) = 600.0_kind_real ! [micron]
         end if
         field2d % fieldName = var_clgefr
         call mpas_pool_add_field(convFields, var_clgefr, field2d)
 !        write(*,*) "end-of ",var_clgefr

      case ( var_clhefr ) !-effective_radius_of_hail_particle :TODO: currently filled w/ default value
         call mpas_pool_get_field(subFields, 'index_qg', field2d_src) !The current MP schemes do not include hail (wsm7 has)
         call mpas_duplicate_field(field2d_src, field2d)
         field2d % array(:,1:ngrid) = 600.0_kind_real !field2d_src%array(:,1:ngrid) * MPAS_JEDI_MILLION_kr ! [m] -> [micron]
         field2d % fieldName = var_clhefr
         call mpas_pool_add_field(convFields, var_clhefr, field2d)
 !        write(*,*) "end-of ",var_clhefr

      case ( var_cldfrac ) !-cloud_area_fraction_in_atmosphere_layer
         call mpas_pool_get_field(subFields, 'temperature', field2d_a)
         call mpas_duplicate_field(field2d_a, field2d)

         call mpas_pool_get_config(geom % domain % blocklist % configs, &
                                   'config_radt_cld_scheme', config_radt_cld_scheme)
         call mpas_pool_get_config(geom % domain % blocklist % configs, &
                                   'config_microp_scheme', config_microp_scheme)

         if ( trim(config_radt_cld_scheme) == MPAS_JEDI_OFF .or. &
              trim(config_microp_scheme) == MPAS_JEDI_OFF ) then
            field2d % array(:,1:ngrid) = MPAS_JEDI_LESSONE_kr
         else
            call mpas_pool_get_field(subFields, 'cldfrac', field2d_src)
            field2d % array(:,1:ngrid) = field2d_src % array(:,1:ngrid)
            where(field2d % array(:,1:ngrid) > MPAS_JEDI_LESSONE_kr)
               field2d % array(:,1:ngrid) = MPAS_JEDI_LESSONE_kr
            end where
            where(field2d % array(:,1:ngrid) < MPAS_JEDI_ZERO_kr)
               field2d % array(:,1:ngrid) = MPAS_JEDI_ZERO_kr
            end where

            ! Assume clouds fill entire column when subgrid fraction is zero everywhere
            if ( all(field2d % array(:,1:ngrid) <= MPAS_JEDI_ZERO_kr) ) then
               field2d % array(:,1:ngrid) = MPAS_JEDI_LESSONE_kr
            end if
         end if
         field2d % fieldName = var_cldfrac
         call mpas_pool_add_field(convFields, var_cldfrac, field2d)

     case ( var_sfc_wtmp, var_sfc_ltmp, var_sfc_itmp, var_sfc_stmp ) !-surface_temperature_where_sea, surface_temperature_where_land, surface_temperature_where_ice, surface_temperature_where_snow
        call mpas_pool_get_field(subFields, 'u10', field1d_src) ! as a dummy array

        !-NOTE: Currently assign "skintemp" for all temperature
        !-TODO: More proper variable for each surface temperature ?
        call mpas_pool_get_array(subFields, "skintemp", r1d_ptr_a) !"ground or water surface temperature"
!        write(*,*) 'MIN/MAX of skintemp=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
        call mpas_duplicate_field(field1d_src, field1d)
        field1d % array(1:ngrid) = r1d_ptr_a(1:ngrid) ! quantity and unit might change
        field1d % fieldName = trim(fieldname(ivar))
        call mpas_pool_add_field(convFields, trim(fieldname(ivar)), field1d)
!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ( var_sfc_sdepth ) !-surface_snow_thickness
        call mpas_pool_get_array(subFields, "snowh", r1d_ptr_a)
!        write(*,*) 'MIN/MAX of snowh=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
        call mpas_pool_get_field(subFields, 'u10', field1d_src) ! as a dummy array
        call mpas_duplicate_field(field1d_src, field1d)
        field1d % array(1:ngrid) = r1d_ptr_a(1:ngrid) * MPAS_JEDI_THOUSAND_kr ! [m] -> [mm]
        field1d % fieldName = var_sfc_sdepth
        call mpas_pool_add_field(convFields, var_sfc_sdepth, field1d)
!        write(*,*) "end-of ",var_sfc_sdepth

     case ( var_sfc_vegfrac ) !-vegetation_area_fraction
        call mpas_pool_get_array(subFields, "vegfra", r1d_ptr_a)
!        write(*,*) 'MIN/MAX of vegfra=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
        call mpas_pool_get_field(subFields, 'u10', field1d_src) ! as a dummy array
        call mpas_duplicate_field(field1d_src, field1d)
        field1d % array(1:ngrid) = r1d_ptr_a(1:ngrid) / 100.0_kind_real ! [unitless, 0~100] = [%, 0~1]
        field1d % fieldName = var_sfc_vegfrac
        call mpas_pool_add_field(convFields, var_sfc_vegfrac, field1d)
!        write(*,*) "end-of ",var_sfc_vegfrac

     case ( var_sfc_lai ) !-leaf_area_index :TODO, has a value for restart file, not for init file
        call mpas_pool_get_array(subFields, "lai", r1d_ptr_a)
!        write(*,*) 'MIN/MAX of lai=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
        call mpas_pool_get_field(subFields, 'u10', field1d_src) ! as a dummy array
        call mpas_duplicate_field(field1d_src, field1d)
        field1d % array(1:ngrid) = 3.5_kind_real !r1d_ptr_a(:) ! convert ?
        field1d % fieldName = var_sfc_lai
        call mpas_pool_add_field(convFields, var_sfc_lai, field1d)
!        write(*,*) "end-of ",var_sfc_lai

     case ( var_sfc_soilm ) !-volume_fraction_of_condensed_water_in_soil : NOTE: use 1st level
        call mpas_pool_get_array(subFields, "smois", r2d_ptr_a)
!        write(*,*) 'MIN/MAX of smois=',minval(r2d_ptr_a(1,:)),maxval(r2d_ptr_a(1,:))
        call mpas_pool_get_field(subFields, 'u10', field1d_src) ! as a dummy array
        call mpas_duplicate_field(field1d_src, field1d)
        field1d % array(1:ngrid) = 0.05_kind_real !r2d_ptr_a(1,:) ! [m3/m3] -> [g/cm3] :TODO ~ Good enough ?? range [0,1]
        field1d % fieldName = var_sfc_soilm
        call mpas_pool_add_field(convFields, var_sfc_soilm, field1d)
!        write(*,*) "end-of ",var_sfc_soilm

     case ( var_sfc_soilt ) !-soil_temperature : NOTE: use 1st level
        call mpas_pool_get_array(subFields, "tslb", r2d_ptr_a)
!        write(*,*) 'MIN/MAX of tslb=',minval(r2d_ptr_a(1,:)),maxval(r2d_ptr_a(1,:))
        call mpas_pool_get_field(subFields, 'u10', field1d_src) ! as a dummy array
        call mpas_duplicate_field(field1d_src, field1d)
        field1d % array(1:ngrid) = r2d_ptr_a(1,1:ngrid) ! [K] -> [K]
        field1d % fieldName = var_sfc_soilt
        call mpas_pool_add_field(convFields, var_sfc_soilt, field1d)
!        write(*,*) "end-of ",var_sfc_soilt

     case ( var_z )  !-geopotential_height, geopotential heights at midpoint
        call mpas_pool_get_field(subFields, 'theta', field2d_src) ! as a dummy array
        call mpas_duplicate_field(field2d_src, field2d)
        call mpas_duplicate_field(field2d_src, field2d_a)

!       calculate midpoint geometricZ (unit: m):
        call geometricZ_full_to_half(geom%zgrid(:,1:ngrid), ngrid, &
                                     geom % nVertLevels,field2d_a%array(:,1:ngrid))

        do i=1,ngrid
           lat = geom%latCell(i) * MPAS_JEDI_RAD2DEG_kr !- to Degrees
           do k=1,geom % nVertLevels
              call geometric2geop(lat,field2d_a%array(k,i), field2d%array(k,i))
           enddo
        enddo

        field2d % fieldName = var_z
        call mpas_pool_add_field(convFields, var_z, field2d)
        call mpas_deallocate_field(field2d_a) ! not used

     case ( var_sfc_z )  !-surface_geopotential_height
        call mpas_pool_get_field(subFields, 'u10', field1d_src) ! as a dummy array
        call mpas_duplicate_field(field1d_src, field1d)

        do i=1,ngrid
           lat = geom%latCell(i) * MPAS_JEDI_RAD2DEG_kr !- to Degrees
           call geometric2geop(lat,geom%zgrid(1,i), field1d%array(i))
        enddo
        field1d % fieldName = var_sfc_z
        call mpas_pool_add_field(convFields, var_sfc_z, field1d)

     case ( var_geomz )  !-height
        call mpas_pool_get_field(subFields, 'theta', field2d_src) ! as a dummy array
        call mpas_duplicate_field(field2d_src, field2d)

!       calculate midpoint geometricZ (unit: m):
        call geometricZ_full_to_half(geom%zgrid(:,1:ngrid), ngrid, &
                                     geom % nVertLevels,field2d%array(:,1:ngrid))

        field2d % fieldName = var_geomz
        call mpas_pool_add_field(convFields, var_geomz, field2d)

     case ( var_sfc_geomz )  !-surface_altitude
        call mpas_pool_get_field(subFields, 'u10', field1d_src) ! as a dummy array
        call mpas_duplicate_field(field1d_src, field1d)

        field1d % fieldName = var_sfc_geomz
        field1d % array (1:ngrid) = geom%zgrid(1,1:ngrid)

        call mpas_pool_add_field(convFields, var_sfc_geomz, field1d)

     case default
        write(*,*) 'Not processed in sub. convert_mpas_field2ufo: ',trim(fieldname(ivar))
        !- TODO: Abort processing when we get here. (Breaks hofx ctests)
        !call abor1_ftn("In convert_mpas_field2ufo: Unimplemented GeoVaLs name: "//trim(fieldname(ivar)))

     end select

   end do !ivar
  
end subroutine convert_mpas_field2ufo

!-------------------------------------------------------------------------------------------

subroutine convert_mpas_field2ufoTL(geom, trajFields, subFields_tl, convFields_tl, fieldname, nfield, ngrid)

   implicit none

   type (mpas_geom),               intent(in)  :: geom          !< geometry
   type (mpas_pool_type), pointer, intent(in ) :: trajFields    !< linearization state for conversion
   type (mpas_pool_type), pointer, intent(in ) :: subFields_tl  !< self % subFields
   type (mpas_pool_type), pointer, intent(out) :: convFields_tl !< pool with geovals variable
   character (len=*),              intent(in ) :: fieldname(:)  !< list of variables for geovals
   integer,                        intent(in ) :: nfield        !< number of variables
   integer,                        intent(in ) :: ngrid         !< number of grid cells

   type (mpas_pool_iterator_type) :: poolItr
   type (mpas_pool_type), pointer :: allFields
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b, r2d_ptr_c
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: traj_r2d_a, traj_r2d_b, traj_r2d_c !BJJ test
   real (kind=kind_real), dimension(:), pointer :: traj_r1d_a

   type (field2DReal), pointer :: field2d, field2d_src, field2d_a, traj_field2d_a
   type (field1DReal), pointer :: field1d, field1d_src
   integer :: ivar, i, k
   real (kind=kind_real), dimension (:,:), allocatable :: pressure_f

   !--- create new pull for ufo_vars
   call mpas_pool_create_pool(convFields_tl)

   do ivar=1, nfield
!     write(*,*) 'convert_mpas_field2ufoTL:inside do/select case, &
!               & ivar, trim(fieldname(ivar))=',ivar,trim(fieldname(ivar))

     select case (trim(fieldname(ivar)))

     case ( var_tv ) !-virtual_temperature
        !get TL variables
        call mpas_pool_get_field(subFields_tl, 'temperature', field2d_src)
        call mpas_pool_get_array(subFields_tl,     'spechum', r2d_ptr_b)
!        write(*,*) 'MIN/MAX of TL temperature(in)=',minval(field2d_src%array),maxval(field2d_src%array)
!        write(*,*) 'MIN/MAX of TL     spechum(in)=',minval(r2d_ptr_b),maxval(r2d_ptr_b)

        !get linearization state
        call mpas_pool_get_array(trajFields, 'temperature', traj_r2d_a)
        call mpas_pool_get_array(trajFields,     'spechum', traj_r2d_b)
!        write(*,*) 'MIN/MAX of TRAJ temperature=',minval(traj_r2d_a),maxval(traj_r2d_a)
!        write(*,*) 'MIN/MAX of TRAJ     spechum=',minval(traj_r2d_b),maxval(traj_r2d_b)

        call mpas_duplicate_field(field2d_src, field2d)        ! for TL of virtual_temperature
        call mpas_duplicate_field(field2d_src, field2d_a)      ! for TL of mixing_ratio, intermediate variable
        call mpas_duplicate_field(field2d_src, traj_field2d_a) ! for NL of mixing_ratio, intermediate variable

        call q_to_w_tl( r2d_ptr_b(:,1:ngrid), traj_r2d_b(:,1:ngrid), field2d_a % array(:,1:ngrid) )
        call q_to_w( traj_r2d_b(:,1:ngrid) , traj_field2d_a % array(:,1:ngrid) ) !NL coeff.
        call tw_to_tv_tl( field2d_src % array(:,1:ngrid), field2d_a % array(:,1:ngrid), &
                          traj_r2d_a(:,1:ngrid), traj_field2d_a % array(:,1:ngrid), &
                          field2d % array(:,1:ngrid) )
        call mpas_deallocate_field(field2d_a)      ! not used
        call mpas_deallocate_field(traj_field2d_a) ! not used 
!        write(*,*) 'MIN/MAX of TL Tv(out)=',minval(field2d % array(:,1:ngrid)),maxval(field2d % array(:,1:ngrid))

        field2d % fieldName = var_tv

        call mpas_pool_add_field(convFields_tl, var_tv, field2d)

!        write(*,*) "end-of ",var_tv

     case ( var_ts, var_t ) !-air_temperature, temperature
        call mpas_pool_get_field(subFields_tl, 'temperature', field2d_src) !< get temperature
        call mpas_duplicate_field(field2d_src, field2d)!  as a dummy array 
        field2d % fieldName = trim(fieldname(ivar))
        call mpas_pool_add_field(convFields_tl, trim(fieldname(ivar)), field2d)

     case ( var_u ) !-eastward_wind
        !get TL variable
        call mpas_pool_get_field(subFields_tl, 'uReconstructZonal', field2d_src)
        call mpas_duplicate_field(field2d_src, field2d)!  as a dummy array
!        write(*,*) 'MIN/MAX of TL zonal wind(in/out)=',minval(field2d % array),maxval(field2d % array)

        !NL: field2d % array = field2d % array
        !TL: field2d % array = field2d % array

        field2d % fieldName = var_u

        call mpas_pool_add_field(convFields_tl, var_u, field2d)

!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ( var_v ) !-northward_wind
        !get TL variable
        call mpas_pool_get_field(subFields_tl, 'uReconstructMeridional', field2d_src)
        call mpas_duplicate_field(field2d_src, field2d)!  as a dummy array
!        write(*,*) 'MIN/MAX of TL meridional wind(in/out)=',minval(field2d % array),maxval(field2d % array)

        !NL: field2d % array = field2d % array
        !TL: field2d % array = field2d % array

        field2d % fieldName = var_v

        call mpas_pool_add_field(convFields_tl, var_v, field2d)

!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ( var_mixr ) !-humidity_mixing_ratio
        call mpas_pool_get_field(subFields_tl, 'spechum', field2d_src) !< get TL of specific_huuumidity
        call mpas_pool_get_array(trajFields,   'spechum', traj_r2d_a)  !< get linearization state
        call mpas_duplicate_field(field2d_src, field2d)

        call q_to_w_tl(field2d_src % array(:,1:ngrid), traj_r2d_a(:,1:ngrid), field2d % array(:,1:ngrid)) 
        where (traj_r2d_a(:,1:ngrid) <= MPAS_JEDI_ZERO_kr)
          field2d % array(:,1:ngrid) = MPAS_JEDI_ZERO_kr
        end where
        field2d % array(:,1:ngrid) = field2d % array(:,1:ngrid) * MPAS_JEDI_THOUSAND_kr
        field2d % fieldName = var_mixr
        call mpas_pool_add_field(convFields_tl, var_mixr, field2d)
!        write(*,*) "end-of ",var_mixr

     case( var_q ) !-specific_humidity
        call mpas_pool_get_field(subFields_tl, "spechum", field2d_src)
        call mpas_duplicate_field(field2d_src, field2d)
        field2d % fieldName = var_q
        call mpas_pool_add_field(convFields_tl, var_q, field2d)
!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ( var_prs ) !-air_pressure
!        call mpas_pool_get_field(subFields_tl, 'pressure', field2d)
!        field2d % array(:,1:ngrid) = field2d_src%array(:,1:ngrid)
!        field2d % fieldName = var_prs
!        call mpas_pool_add_field(convFields_tl, var_prs, field2d)
!        write(*,*) "end-of ",var_prs

     case ( var_ps ) !-surface_pressure
        call mpas_pool_get_field(subFields_tl, 'surface_pressure', field1d_src) !< get surface_pressure
        call mpas_duplicate_field(field1d_src, field1d)!  as a dummy array
        field1d % fieldName = trim(fieldname(ivar))
        call mpas_pool_add_field(convFields_tl, trim(fieldname(ivar)), field1d)
     case ( var_prsi ) !-air_pressure_levels
     case ( var_oz )   !-mole_fraction_of_ozone_in_air
     case ( var_co2 )  !-mole_fraction_of_carbon_dioxide_in_air

     case ( var_clw ) !-mass_content_of_cloud_liquid_water_in_atmosphere_layer
        ! Calcluate air_pressure_levels
        call mpas_pool_get_array(trajFields, 'pressure', traj_r2d_b)
        call mpas_pool_get_array(trajFields, 'surface_pressure', traj_r1d_a)
        allocate (pressure_f(geom%nVertLevels+1,ngrid))
        call pressure_half_to_full(traj_r2d_b(:,1:ngrid), geom%zgrid(:,1:ngrid), traj_r1d_a(1:ngrid), &
                                   ngrid, geom%nVertLevels, pressure_f)
        call index_q_fields_TL('index_qc', var_clw, subFields_tl, convFields_tl, pressure_f, ngrid, geom % nVertLevels)
        deallocate (pressure_f)

     case ( var_cli ) !-mass_content_of_cloud_ice_in_atmosphere_layer
        ! Calcluate air_pressure_levels
        call mpas_pool_get_array(trajFields, 'pressure', traj_r2d_b)
        call mpas_pool_get_array(trajFields, 'surface_pressure', traj_r1d_a)
        allocate (pressure_f(geom%nVertLevels+1,ngrid))
        call pressure_half_to_full(traj_r2d_b(:,1:ngrid), geom%zgrid(:,1:ngrid), traj_r1d_a(1:ngrid), &
                                   ngrid, geom%nVertLevels, pressure_f)
        call index_q_fields_TL('index_qi', var_cli, subFields_tl, convFields_tl, pressure_f, ngrid, geom % nVertLevels)
        deallocate (pressure_f)

      case ( var_clr ) !-mass_content_of_rain_in_atmosphere_layer
         ! Calcluate air_pressure_levels
         call mpas_pool_get_array(trajFields, 'pressure', traj_r2d_b)
         call mpas_pool_get_array(trajFields, 'surface_pressure', traj_r1d_a)
         allocate (pressure_f(geom%nVertLevels+1,ngrid))
         call pressure_half_to_full(traj_r2d_b(:,1:ngrid), geom%zgrid(:,1:ngrid), traj_r1d_a(1:ngrid), &
                                    ngrid, geom%nVertLevels, pressure_f)
         call index_q_fields_TL('index_qr', var_clr, subFields_tl, convFields_tl, pressure_f, ngrid, geom % nVertLevels)
         deallocate (pressure_f)

      case ( var_cls ) !-mass_content_of_snow_in_atmosphere_layer
         ! Calcluate air_pressure_levels
         call mpas_pool_get_array(trajFields, 'pressure', traj_r2d_b)
         call mpas_pool_get_array(trajFields, 'surface_pressure', traj_r1d_a)
         allocate (pressure_f(geom%nVertLevels+1,ngrid))
         call pressure_half_to_full(traj_r2d_b(:,1:ngrid), geom%zgrid(:,1:ngrid), traj_r1d_a(1:ngrid), &
                                    ngrid, geom%nVertLevels, pressure_f)
         call index_q_fields_TL('index_qs', var_cls, subFields_tl, convFields_tl, pressure_f, ngrid, geom % nVertLevels)
         deallocate (pressure_f)

      case ( var_clg ) !-mass_content_of_graupel_in_atmosphere_layer
         ! Calcluate air_pressure_levels
         call mpas_pool_get_array(trajFields, 'pressure', traj_r2d_b)
         call mpas_pool_get_array(trajFields, 'surface_pressure', traj_r1d_a)
         allocate (pressure_f(geom%nVertLevels+1,ngrid))
         call pressure_half_to_full(traj_r2d_b(:,1:ngrid), geom%zgrid(:,1:ngrid), traj_r1d_a(1:ngrid), &
                                    ngrid, geom%nVertLevels, pressure_f)
         call index_q_fields_TL('index_qg', var_clg, subFields_tl, convFields_tl, pressure_f, ngrid, geom % nVertLevels)
         deallocate (pressure_f)

      case ( var_clh ) !-mass_content_of_hail_in_atmosphere_layer
         ! Calcluate air_pressure_levels
         call mpas_pool_get_array(trajFields, 'pressure', traj_r2d_b)
         call mpas_pool_get_array(trajFields, 'surface_pressure', traj_r1d_a)
         allocate (pressure_f(geom%nVertLevels+1,ngrid))
         call pressure_half_to_full(traj_r2d_b(:,1:ngrid), geom%zgrid(:,1:ngrid), traj_r1d_a(1:ngrid), &
                                    ngrid, geom%nVertLevels, pressure_f)
         call index_q_fields_TL('index_qh', var_clh, subFields_tl, convFields_tl, pressure_f, ngrid, geom % nVertLevels)
         deallocate (pressure_f)

     case ( var_clwefr ) !-effective_radius_of_cloud_liquid_water_particle
     case ( var_cliefr ) !-effective_radius_of_cloud_ice_particle

     case ( var_sfc_wfrac ) !-water_area_fraction
!JJG: All non-TL trajectory model variables should be provided in trajFields
!JJG: All TL model variables should be provided in subFields_tl
!        call mpas_pool_get_array(subFields_tl, "", r1d_ptr_a)
     case ( var_sfc_lfrac ) !-land_area_fraction
!        call mpas_pool_get_array(subFields_tl, "", r1d_ptr_a)
     case ( var_sfc_ifrac ) !-ice_area_fraction
!        call mpas_pool_get_array(subFields_tl, "", r1d_ptr_a)
     case ( var_sfc_sfrac ) !-surface_snow_area_fraction
!        call mpas_pool_get_array(subFields_tl, "", r1d_ptr_a)
     case ( var_sfc_wtmp) !-surface_temperature_where_sea
!        call mpas_pool_get_array(subFields_tl, "skintemp", r1d_ptr_a)
     case ( var_sfc_ltmp ) !-surface_temperature_where_land
!        call mpas_pool_get_array(subFields_tl, "skintemp", r1d_ptr_a)
     case ( var_sfc_itmp ) !-surface_temperature_where_ice
!        call mpas_pool_get_array(subFields_tl, "skintemp", r1d_ptr_a)
     case ( var_sfc_stmp ) !-surface_temperature_where_snow
!        call mpas_pool_get_array(subFields_tl, "skintemp", r1d_ptr_a)

     case ( var_sfc_sdepth ) !-surface_snow_thickness
!        call mpas_pool_get_array(subFields_tl, "snowh", r1d_ptr_a)
     case ( var_sfc_vegfrac ) !-vegetation_area_fraction
!        call mpas_pool_get_array(subFields_tl, "vegfra", r1d_ptr_a)
     case (var_sfc_wspeed, var_sfc_wdir) !-surface_wind_speed, surface_wind_from_direction
!        call mpas_pool_get_array(subFields_tl, "u10", r1d_ptr_a)
!        call mpas_pool_get_array(subFields_tl, "v10", r1d_ptr_b)
     case ( var_sfc_lai ) !-leaf_area_index
!        call mpas_pool_get_array(subFields_tl, "lai", r1d_ptr_a)
     case ( var_sfc_soilm ) !-volume_fraction_of_condensed_water_in_soil
!        call mpas_pool_get_array(subFields_tl, "smois", r1d_ptr_a)
     case ( var_sfc_soilt ) !-soil_temperature
!        call mpas_pool_get_array(subFields_tl, "tslb", r1d_ptr_a)
     case ( var_sfc_landtyp ) !-land_type_index
!        call mpas_pool_get_array(subFields_tl, "", r1d_ptr_a)
     case ( var_sfc_vegtyp ) !-vegetation_type_index
!        call mpas_pool_get_array(subFields_tl, "ivgtyp", r1d_ptr_a)
     case ( var_sfc_soiltyp ) !-soil_type
!        call mpas_pool_get_array(subFields_tl, "isltyp", r1d_ptr_a)

     end select

   end do !ivar

  
end subroutine convert_mpas_field2ufoTL

!-------------------------------------------------------------------------------------------

subroutine convert_mpas_field2ufoAD(geom, trajFields, subFields_ad, convFields_ad, fieldname, nfield, ngrid)

   implicit none

   type (mpas_geom),               intent(in)    :: geom          !< geometry
   type (mpas_pool_type), pointer, intent(in   ) :: trajFields    !< linearization state for conversion
   type (mpas_pool_type), pointer, intent(inout) :: subFields_ad  !< self % subFields
   type (mpas_pool_type), pointer, intent(in   ) :: convFields_ad !< pool with geovals variable
   character (len=*),              intent(in   ) :: fieldname(:)  !< list of variables for geovals
   integer,                        intent(in   ) :: nfield        !< number of variables
   integer,                        intent(in   ) :: ngrid         !< number of grid cells

   type (mpas_pool_iterator_type) :: poolItr
   type (mpas_pool_type), pointer :: allFields, clone_subFields_ad
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b, r2d_ptr_c
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: traj_r2d_a, traj_r2d_b, traj_r2d_c !BJJ test
   real (kind=kind_real), dimension(:), pointer :: traj_r1d_a

   type (field2DReal), pointer :: field2d, field2d_src, field2d_a, traj_field2d_a
   type (field1DReal), pointer :: field1d, field1d_src
   integer :: ivar, i, k
   real (kind=kind_real), dimension (:,:), allocatable :: pressure_f
   real (kind=kind_real) :: kgkg_kgm2 !-- for var_clw, var_cli

   do ivar=1, nfield
!     write(*,*) 'convert_mpas_field2ufoAD: inside do/select case, &
!               & ivar, trim(fieldname(ivar))=',ivar,trim(fieldname(ivar))

     select case (trim(fieldname(ivar)))

     case ( var_tv ) !-virtual_temperature
        !get AD variables
        call mpas_pool_get_array(subFields_ad, 'temperature', r2d_ptr_a)
        call mpas_pool_get_array(subFields_ad,     'spechum', r2d_ptr_b)
!        write(*,*) 'MIN/MAX of AD temperature(in)=',minval(r2d_ptr_a),maxval(r2d_ptr_a)
!        write(*,*) 'MIN/MAX of AD     spechum(in)=',minval(r2d_ptr_b),maxval(r2d_ptr_b)

        !get linearization state
        call mpas_pool_get_array(trajFields, 'temperature', traj_r2d_a)
        call mpas_pool_get_array(trajFields,     'spechum', traj_r2d_b)
!        write(*,*) 'MIN/MAX of TRAJ temperature=',minval(traj_r2d_a),maxval(traj_r2d_a)
!        write(*,*) 'MIN/MAX of TRAJ     spechum=',minval(traj_r2d_b),maxval(traj_r2d_b)

        call mpas_pool_get_field(convFields_ad, var_tv, field2d)
!        write(*,*) 'MIN/MAX of AD Tv(in)=',minval(field2d % array),maxval(field2d % array)

        call mpas_duplicate_field(field2d, field2d_a)      ! for AD of mixing_ratio, intermediate variable
        call mpas_duplicate_field(field2d, traj_field2d_a) ! for NL of mixing_ratio, intermediate variable

        call q_to_w( traj_r2d_b(:,1:ngrid) , traj_field2d_a % array(:,1:ngrid) ) !NL coeff.
        field2d_a % array(:,1:ngrid) = MPAS_JEDI_ZERO_kr !initialize local var.
        call tw_to_tv_ad(r2d_ptr_a(:,1:ngrid), field2d_a % array(:,1:ngrid), &
                         traj_r2d_a(:,1:ngrid), traj_field2d_a % array(:,1:ngrid), &
                         field2d % array(:,1:ngrid) )
        call q_to_w_ad( r2d_ptr_b(:,1:ngrid), traj_r2d_b(:,1:ngrid), field2d_a % array(:,1:ngrid) )
        call mpas_deallocate_field(field2d_a) ! not used
        call mpas_deallocate_field(traj_field2d_a) ! not used
!        write(*,*) 'MIN/MAX of AD temperature(out)=',minval(r2d_ptr_a),maxval(r2d_ptr_a)
!        write(*,*) 'MIN/MAX of AD     spechum(out)=',minval(r2d_ptr_b),maxval(r2d_ptr_b)

!        write(*,*) "end-of ",var_tv

     case ( var_ts, var_t ) !-air_temperature, temperature
        call mpas_pool_get_array(subFields_ad, 'temperature', r2d_ptr_a) !< get temperature
        call mpas_pool_get_field(convFields_ad, trim(fieldname(ivar)), field2d)
        r2d_ptr_a(:,1:ngrid) = r2d_ptr_a(:,1:ngrid) + &
                         field2d % array(:,1:ngrid)

     case ( var_u ) !-eastward_wind
        !get AD variable
        call mpas_pool_get_array(subFields_ad, 'uReconstructZonal', r2d_ptr_a) !< get zonal wind
!        write(*,*) 'MIN/MAX of AD zonal wind(in)  =',minval(r2d_ptr_a),maxval(r2d_ptr_a)

        !NL: field2d % array = field2d % array
        !TL: field2d % array = field2d % array
        call mpas_pool_get_field(convFields_ad, var_u, field2d)

        r2d_ptr_a(:,1:ngrid) = r2d_ptr_a(:,1:ngrid) + &
                         field2d % array(:,1:ngrid)

!        write(*,*) 'MIN/MAX of AD zonal wind(out) =',minval(r2d_ptr_a),maxval(r2d_ptr_a)

!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ( var_v ) !-northward_wind
        !get AD variable
        call mpas_pool_get_array(subFields_ad, 'uReconstructMeridional', r2d_ptr_a) !< get meridional wind
!        write(*,*) 'MIN/MAX of AD meridional wind(in)  =',minval(r2d_ptr_a),maxval(r2d_ptr_a)

        !NL: field2d % array = field2d % array
        !TL: field2d % array = field2d % array
        call mpas_pool_get_field(convFields_ad, var_v, field2d)

        r2d_ptr_a(:,1:ngrid) = r2d_ptr_a(:,1:ngrid) + &
                         field2d % array(:,1:ngrid)

!        write(*,*) 'MIN/MAX of AD meridional wind(out) =',minval(r2d_ptr_a),maxval(r2d_ptr_a)

!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ( var_mixr ) !-humidity_mixing_ratio
        call mpas_pool_get_array(subFields_ad, "spechum", r2d_ptr_a) !< get sphechum_ad
        call mpas_pool_get_array(trajFields,   'spechum', traj_r2d_a)  !< get linearization state
!        write(*,*) 'MIN/MAX of AD spechum(in) =',minval(r2d_ptr_a(:,1:ngrid)),maxval(r2d_ptr_a(:,1:ngrid))

        call mpas_pool_get_field(convFields_ad, var_mixr, field2d)
!        write(*,*) 'MIN/MAX of AD var_mixr =',minval(field2d % array(:,1:ngrid)),maxval(field2d % array(:,1:ngrid))

        field2d % array(:,1:ngrid) = field2d % array(:,1:ngrid) * MPAS_JEDI_THOUSAND_kr
        call q_to_w_ad(r2d_ptr_a(:,1:ngrid), traj_r2d_a(:,1:ngrid), field2d % array(:,1:ngrid))
        where (traj_r2d_a(:,1:ngrid) <= MPAS_JEDI_ZERO_kr)
          r2d_ptr_a(:,1:ngrid) = MPAS_JEDI_ZERO_kr
        end where
!        write(*,*) 'MIN/MAX of AD spechum(out) =',minval(r2d_ptr_a(:,1:ngrid)),maxval(r2d_ptr_a(:,1:ngrid))

!        write(*,*) "end-of ",var_mixr

     case( var_q ) !-specific_humidity
        call mpas_pool_get_array(subFields_ad, "spechum", r2d_ptr_a)
        call mpas_pool_get_field(convFields_ad, var_q, field2d)
        r2d_ptr_a(:,1:ngrid) = r2d_ptr_a(:,1:ngrid) + field2d % array(:,1:ngrid)
!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ( var_prs ) !-air_pressure
!        call mpas_pool_get_array(subFields_ad, "pressure", r2d_ptr_a)
!        call mpas_pool_get_field(clone_subFields_ad, 'pressure', field2d) ! as a dummy array
!        field2d % array(:,1:ngrid) = r2d_ptr_a(:,1:ngrid)
!        field2d % fieldName = var_prs
!        call mpas_pool_add_field(convFields_ad, var_prs, field2d)
!        write(*,*) "end-of ",var_prs

     case ( var_ps ) !-surface_pressure
        call mpas_pool_get_array(subFields_ad, 'surface_pressure', r1d_ptr_a) !< get psfc
        call mpas_pool_get_field(convFields_ad, trim(fieldname(ivar)), field1d)
        r1d_ptr_a(1:ngrid) = r1d_ptr_a(1:ngrid) + &
                         field1d % array(1:ngrid)

     case ( var_prsi ) !-air_pressure_levels
     case ( var_oz )   !-mole_fraction_of_ozone_in_air
     case ( var_co2 )  !-mole_fraction_of_carbon_dioxide_in_air

     case ( var_clw ) !-mass_content_of_cloud_liquid_water_in_atmosphere_layer
        call index_q_fields_AD(geom, trajFields, 'index_qc', var_clw, subFields_ad, convFields_ad, ngrid)

      case ( var_cli ) !-mass_content_of_cloud_ice_in_atmosphere_layer
         call index_q_fields_AD(geom, trajFields, 'index_qi', var_cli, subFields_ad, convFields_ad, ngrid)

      case ( var_clr ) !-mass_content_of_rain_in_atmosphere_layer
         call index_q_fields_AD(geom, trajFields, 'index_qr', var_clr, subFields_ad, convFields_ad, ngrid)

      case ( var_cls ) !-mass_content_of_snow_in_atmosphere_layer
         call index_q_fields_AD(geom, trajFields, 'index_qs', var_cls, subFields_ad, convFields_ad, ngrid)

      case ( var_clg ) !-mass_content_of_graupel_in_atmosphere_layer
         call index_q_fields_AD(geom, trajFields, 'index_qg', var_clg, subFields_ad, convFields_ad, ngrid)

      case ( var_clh ) !-mass_content_of_hail_in_atmosphere_layer
         call index_q_fields_AD(geom, trajFields, 'index_qh', var_clh, subFields_ad, convFields_ad, ngrid)

     case ( var_clwefr ) !-effective_radius_of_cloud_liquid_water_particle
     case ( var_cliefr ) !-effective_radius_of_cloud_ice_particle

     case ( var_sfc_wfrac ) !-water_area_fraction
!JJG: All non-AD trajectory model variables should be provided in trajFields
!JJG: All AD model variables should be provided in subFields_ad
!        call mpas_pool_get_array(subFields_ad, "", r1d_ptr_a)
     case ( var_sfc_lfrac ) !-land_area_fraction
!        call mpas_pool_get_array(subFields_ad, "", r1d_ptr_a)
     case ( var_sfc_ifrac ) !-ice_area_fraction
!        call mpas_pool_get_array(subFields_ad, "", r1d_ptr_a)
     case ( var_sfc_sfrac ) !-surface_snow_area_fraction
!        call mpas_pool_get_array(subFields_ad, "", r1d_ptr_a)
     case ( var_sfc_wtmp) !-surface_temperature_where_sea
!        call mpas_pool_get_array(subFields_ad, "skintemp", r1d_ptr_a)
     case ( var_sfc_ltmp ) !-surface_temperature_where_land
!        call mpas_pool_get_array(subFields_ad, "skintemp", r1d_ptr_a)
     case ( var_sfc_itmp ) !-surface_temperature_where_ice
!        call mpas_pool_get_array(subFields_ad, "skintemp", r1d_ptr_a)
     case ( var_sfc_stmp ) !-surface_temperature_where_snow
!        call mpas_pool_get_array(subFields_ad, "skintemp", r1d_ptr_a)
     case ( var_sfc_sdepth ) !-surface_snow_thickness
!        call mpas_pool_get_array(subFields_ad, "snowh", r1d_ptr_a)
     case ( var_sfc_vegfrac ) !-vegetation_area_fraction
!        call mpas_pool_get_array(subFields_ad, "vegfra", r1d_ptr_a)
     case (var_sfc_wspeed, var_sfc_wdir) !-surface_wind_speed, surface_wind_from_direction
!        call mpas_pool_get_array(subFields_ad, "u10", r1d_ptr_a)
!        call mpas_pool_get_array(subFields_ad, "v10", r1d_ptr_b)
     case ( var_sfc_lai ) !-leaf_area_index
!        call mpas_pool_get_array(subFields_ad, "lai", r1d_ptr_a)
     case ( var_sfc_soilm ) !-volume_fraction_of_condensed_water_in_soil
!        call mpas_pool_get_array(subFields_ad, "smois", r1d_ptr_a)
     case ( var_sfc_soilt ) !-soil_temperature
!        call mpas_pool_get_array(subFields_ad, "tslb", r1d_ptr_a)
     case ( var_sfc_landtyp ) !-land_type_index
!        call mpas_pool_get_array(subFields_ad, "", r1d_ptr_a)
     case ( var_sfc_vegtyp ) !-vegetation_type_index
!        call mpas_pool_get_array(subFields_ad, "ivgtyp", r1d_ptr_a)
     case ( var_sfc_soiltyp ) !-soil_type
!        call mpas_pool_get_array(subFields_ad, "isltyp", r1d_ptr_a)

     end select

   end do !ivar

   !call mpas_pool_destroy_pool(clone_subFields_ad)
  
end subroutine convert_mpas_field2ufoAD

!-----------------------------------------------------------------------------------

   end module mpas2ufo_vars_mod
