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
use mpas_constants, only : gravity, rgas, rv, cp, pii
use mpas_dmpar
use mpas_derived_types
use mpas_field_routines
use mpas_pool_routines

private

public :: convert_type_soil, convert_type_veg
public :: uv_to_wdir
public :: w_to_q, q_to_w
public :: theta_to_temp, temp_to_theta
public :: twp_to_rho

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
   real(kind=kind_real),parameter:: quadcof  (4, 2  ) =      &
      reshape((/0.0_kind_real, 1.0_kind_real, 1.0_kind_real, 2.0_kind_real, 1.0_kind_real, &
               -1.0_kind_real, 1.0_kind_real, -1.0_kind_real/), (/4, 2/))

   if (uu5 >= 0.0_kind_real .and. vv5 >= 0.0_kind_real) iquadrant = 1
   if (uu5 >= 0.0_kind_real .and. vv5 <  0.0_kind_real) iquadrant = 2
   if (uu5 <  0.0_kind_real .and. vv5 >= 0.0_kind_real) iquadrant = 4
   if (uu5 <  0.0_kind_real .and. vv5 <  0.0_kind_real) iquadrant = 3
   if (abs(vv5) >= windlimit) then
      windratio = uu5 / vv5
   else
      windratio = 0.0_kind_real
      if (abs(uu5) > windlimit) then
         windratio = windscale * uu5
      endif
   endif
   windangle        = atan(abs(windratio))   ! wind azimuth is in radians
   wind10_direction = ( quadcof(iquadrant, 1) * pii + windangle * quadcof(iquadrant, 2) )

end subroutine uv_to_wdir

!-------------------------------------------------------------------------------------------
elemental subroutine w_to_q(mixing_ratio, specific_humidity)
   implicit none
   real (kind=kind_real), intent(in)  :: mixing_ratio
   real (kind=kind_real), intent(out) :: specific_humidity

   specific_humidity = mixing_ratio / (1.0_kind_real + mixing_ratio)
end subroutine w_to_q
!-------------------------------------------------------------------------------------------
elemental subroutine q_to_w(specific_humidity, mixing_ratio)
   implicit none
   real (kind=kind_real), intent(in)  :: specific_humidity
   real (kind=kind_real), intent(out) :: mixing_ratio

   mixing_ratio = specific_humidity / (1.0_kind_real - specific_humidity)
end subroutine q_to_w
!-------------------------------------------------------------------------------------------
elemental subroutine q_to_w_tl(specific_humidity_tl, sh_traj, mixing_ratio_tl)
   implicit none
   real (kind=kind_real), intent(in)  :: specific_humidity_tl
   real (kind=kind_real), intent(in)  :: sh_traj
   real (kind=kind_real), intent(out) :: mixing_ratio_tl

   mixing_ratio_tl = specific_humidity_tl / (1.0_kind_real - sh_traj)**2
end subroutine q_to_w_tl
!-------------------------------------------------------------------------------------------
elemental subroutine q_to_w_ad(specific_humidity_ad, sh_traj, mixing_ratio_ad)
   implicit none
   real (kind=kind_real), intent(inout) :: specific_humidity_ad
   real (kind=kind_real), intent(in)    :: sh_traj
   real (kind=kind_real), intent(in)    :: mixing_ratio_ad

   specific_humidity_ad = specific_humidity_ad + &
                 1.0_kind_real / ( 1.0_kind_real - sh_traj)**2 * mixing_ratio_ad
end subroutine q_to_w_ad
!-------------------------------------------------------------------------------------------
elemental subroutine tw_to_tv(temperature,mixing_ratio,virtual_temperature)
   implicit none
   real (kind=kind_real), intent(in)  :: temperature
   real (kind=kind_real), intent(in)  :: mixing_ratio
   real (kind=kind_real), intent(out) :: virtual_temperature

   virtual_temperature = temperature * &
                  ( 1.0_kind_real + (rv/rgas - 1.0_kind_real)*mixing_ratio )
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
                  ( 1.0_kind_real + (rv/rgas - 1.0_kind_real)*m_traj )  + &
                  t_traj * (rv/rgas - 1.0_kind_real)*mixing_ratio_tl
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
                  ( 1.0_kind_real + (rv/rgas - 1.0_kind_real)*m_traj )
   mixing_ratio_ad = mixing_ratio_ad + virtual_temperature_ad * &
                  t_traj * (rv/rgas - 1.0_kind_real)
end subroutine tw_to_tv_ad
!-------------------------------------------------------------------------------------------
elemental subroutine theta_to_temp(theta,pressure,temperature)
   implicit none
   real (kind=kind_real), intent(in)  :: theta
   real (kind=kind_real), intent(in)  :: pressure
   real (kind=kind_real), intent(out) :: temperature
   temperature = theta / &
             ( 100000.0_kind_real / pressure ) ** ( rgas / cp )
   !TODO: Following formula would give the same result with the formular above,
   !      but gives a slightly different precision. Need to check.
   !temperature = theta * &
   !          ( pressure / 100000.0_kind_real ) ** ( rgas / cp )
end subroutine theta_to_temp
!-------------------------------------------------------------------------------------------
elemental subroutine temp_to_theta(temperature,pressure,theta)
   implicit none
   real (kind=kind_real), intent(in)  :: temperature
   real (kind=kind_real), intent(in)  :: pressure
   real (kind=kind_real), intent(out) :: theta
   theta = temperature * &
             ( 100000.0_kind_real / pressure ) ** ( rgas / cp )
end subroutine temp_to_theta
!-------------------------------------------------------------------------------------------
elemental subroutine twp_to_rho(temperature,mixing_ratio,pressure,rho)
   implicit none
   real (kind=kind_real), intent(in)  :: temperature
   real (kind=kind_real), intent(in)  :: mixing_ratio
   real (kind=kind_real), intent(in)  :: pressure
   real (kind=kind_real), intent(out) :: rho
   rho = pressure / ( rgas * temperature * &
                                ( 1.0_kind_real + (rv/rgas - 1.0_kind_real) * mixing_ratio ) )
end subroutine twp_to_rho
!-------------------------------------------------------------------------------------------
subroutine pressure_half_to_full(pressure, zgrid, nC, nV, pressure_f)
   implicit none
   real (kind=kind_real), dimension(nV,nC), intent(in) :: pressure
   real (kind=kind_real), dimension(nV+1,nC), intent(in) :: zgrid
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
          tem1 = 1.0_kind_real/(zgrid(k+1,i)- zgrid(k-1,i))
          fzm_p(i,k) = ( zgrid(k,i)- zgrid(k-1,i)) * tem1
          fzp_p(i,k) = ( zgrid(k+1,i)- zgrid(k,i)) * tem1
          pressure_f(k,i) = fzm_p(i,k)*pressure(k,i) + fzp_p(i,k)*pressure(k-1,i)
        enddo
        enddo
        k = kte+1
        do i = its,ite
          z0 = zgrid(k,i)
          z1 = 0.5_kind_real*(zgrid(k,i)+zgrid(k-1,i))
          z2 = 0.5_kind_real*(zgrid(k-1,i)+zgrid(k-2,i))
          w1 = (z0-z2)/(z1-z2)
          w2 = 1.0_kind_real-w1
          !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
          pressure_f(k,i) = exp( w1*log(pressure(k-1,i)) + w2*log(pressure(k-1,i)) )
        enddo
        k = kts
        do i = its,ite
          z0 = zgrid(k,i)
          z1 = 0.5_kind_real*(zgrid(k,i)+zgrid(k+1,i))
          z2 = 0.5_kind_real*(zgrid(k+1,i)+zgrid(k+2,i))
          w1 = (z0-z2)/(z1-z2)
          w2 = 1.0_kind_real-w1
          pressure_f(k,i) = w1*pressure(k,i) + w2*pressure(k+1,i)
        enddo

end subroutine pressure_half_to_full
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
         zgrid(k,i) = ( zgrid_f(k,i) + zgrid_f(k+1,i) ) * 0.5_kind_real
      enddo
   enddo
end subroutine geometricZ_full_to_half
!-------------------------------------------------------------------------------------------

   !--- variables can be found in subFields
subroutine convert_mpas_field2ufo(geom, subFields, convFields, fieldname, nfield, ngrid)

   use mpas_geom_mod

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

   real (kind=kind_real) :: kgkg_kgm2 !-- for var_clw, var_cli

   real (kind=kind_real), parameter :: deg2rad = pii/180.0_kind_real
   real (kind=kind_real) :: lat

   !--- create new pool for geovals
   call mpas_pool_create_pool(convFields)

   do ivar=1, nfield
     write(*,*) 'convert_mpas_field2ufo  :inside do/select case, &
               & ivar, trim(fieldname(ivar))=',ivar,trim(fieldname(ivar))

     select case (trim(fieldname(ivar)))

     case ( "virtual_temperature" ) !-var_tv 
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

     case ( "air_temperature", "temperature" ) !-var_ts
        call mpas_pool_get_field(subFields, 'temperature', field2d_src) !< get temperature
        call mpas_duplicate_field(field2d_src, field2d)! for air_temperature
        field2d % fieldName = trim(fieldname(ivar))
        call mpas_pool_add_field(convFields, trim(fieldname(ivar)), field2d)

     case ( "eastward_wind" ) !-var_??? eastward_wind

        call mpas_pool_get_field(subFields, 'uReconstructZonal', field2d_src) !< get zonal wind
        call mpas_duplicate_field(field2d_src, field2d)!  for eastward_wind
!        write(*,*) 'MIN/MAX of zonal wind=',minval(field2d % array),maxval(field2d % array)
        field2d % fieldName = trim(fieldname(ivar))
        call mpas_pool_add_field(convFields, trim(fieldname(ivar)), field2d)
!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ( "northward_wind" ) !-var_??? northward_wind
        call mpas_pool_get_field(subFields, 'uReconstructMeridional', field2d_src) !< get meridional wind
        call mpas_duplicate_field(field2d_src, field2d)!  for northward_wind
!        write(*,*) 'MIN/MAX of meridional wind=',minval(field2d % array),maxval(field2d % array)
        field2d % fieldName = trim(fieldname(ivar))
        call mpas_pool_add_field(convFields, trim(fieldname(ivar)), field2d)
!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ("atmosphere_ln_pressure_coordinate") !-var_prsl
        call mpas_pool_get_field(subFields, 'pressure', field2d_src) !< get pressure
        call mpas_duplicate_field(field2d_src, field2d)

        field2d % array(:,1:ngrid) = log( field2d_src%array(:,1:ngrid) / 100.0_kind_real / 10.0_kind_real ) !< unit: Pa -> hPa ->cb
!        write(*,*) 'MIN/MAX of ln_p=',minval(field2d % array(:,1:ngrid)),maxval(field2d % array(:,1:ngrid))

        field2d % fieldName = var_prsl

        call mpas_pool_add_field(convFields, var_prsl, field2d)

!        write(*,*) "end-of ",var_prsl

     case ("humidity_mixing_ratio") !-var_mixr
        call mpas_pool_get_field(subFields, 'spechum', field2d_src) !< get specific_humidity
        call mpas_duplicate_field(field2d_src, field2d)! for humidity_mixing_ratio
        call q_to_w(field2d_src % array(:,1:ngrid) , field2d % array(:,1:ngrid))
        field2d % array(:,1:ngrid) = max(0.0_kind_real, field2d % array(:,1:ngrid)) * 1000.0_kind_real ! [kg/kg] -> [g/kg]
        field2d % fieldName = var_mixr
        call mpas_pool_add_field(convFields, var_mixr, field2d)
!        write(*,*) "end-of ",var_mixr

     case("specific_humidity")
        call mpas_pool_get_field(subFields, "spechum", field2d_src)
        call mpas_duplicate_field(field2d_src, field2d)
        field2d % fieldName = trim(fieldname(ivar))
        call mpas_pool_add_field(convFields, trim(fieldname(ivar)), field2d)
!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ("air_pressure") !-var_prs
        call mpas_pool_get_field(subFields, 'pressure', field2d_src) !< get pressure
        call mpas_duplicate_field(field2d_src, field2d)
        field2d % array(:,1:ngrid) = field2d_src%array(:,1:ngrid) / 100.0_kind_real ! [Pa] -> [hPa]
        field2d % fieldName = var_prs
        call mpas_pool_add_field(convFields, var_prs, field2d)
!        write(*,*) "end-of ",var_prs

     case ("air_pressure_levels") !-var_prsi
        call mpas_pool_get_array(subFields, "pressure", r2d_ptr_a)
        call mpas_pool_get_field(subFields, 'w', field2d_src) ! as a dummy array
        call mpas_duplicate_field(field2d_src, field2d)

        call pressure_half_to_full(r2d_ptr_a(:,1:ngrid), geom%zgrid(:,1:ngrid), ngrid, geom % nVertLevels, field2d%array(:,1:ngrid))
        field2d % array = field2d % array / 100.0_kind_real ! [Pa] -> [hPa]
!        write(*,*) 'MIN/MAX of prsi=',minval(field2d % array),maxval(field2d % array)
!        write(*,*) 'test prs       =',r2d_ptr_a(:,1)
!        write(*,*) 'test prsi      =',field2d % array(:,1)

        field2d % fieldName = var_prsi
        call mpas_pool_add_field(convFields, var_prsi, field2d)
!        write(*,*) "end-of ",var_prsi

     case ("mass_concentration_of_ozone_in_air") !-var_oz :TODO: not directly available from MPAS
!        call mpas_pool_get_array(subFields, "o3", r2d_ptr_a)
        call mpas_pool_get_field(subFields, 'theta', field2d_src) ! as a dummy array
        call mpas_duplicate_field(field2d_src, field2d)
        field2d % array(:,1:ngrid) = 0.0_kind_real !r2d_ptr_a(:,1:ngrid) ! convert ??
        field2d % fieldName = var_oz
        call mpas_pool_add_field(convFields, var_oz, field2d)
!        write(*,*) "end-of ",var_oz

     case ("mass_concentration_of_carbon_dioxide_in_air") !-var_co2 :TODO: not directly available from MPAS
!        call mpas_pool_get_array(subFields, "co2", r2d_ptr_a)
        call mpas_pool_get_field(subFields, 'theta', field2d_src) ! as a dummy array
        call mpas_duplicate_field(field2d_src, field2d)
        field2d % array(:,1:ngrid) = 0.0_kind_real !r2d_ptr_a(:,1:ngrid) ! convert ??
        field2d % fieldName = var_co2
        call mpas_pool_add_field(convFields, var_co2, field2d)
!        write(*,*) "end-of ",var_co2

     case ("atmosphere_mass_content_of_cloud_liquid_water") !-var_clw 
        call mpas_pool_get_field(subFields, 'index_qc', field2d_src) !- [kg/kg]
        call mpas_duplicate_field(field2d_src, field2d)
        !--TODO: Trial: Should already have "var_prsi"
        call mpas_pool_get_array(convFields, "air_pressure_levels", r2d_ptr_b) !- [hPa]
        do i=1,ngrid
        do k=1,geom % nVertLevels
          kgkg_kgm2=( r2d_ptr_b(k,i)-r2d_ptr_b(k+1,i) ) * 100.0_kind_real / gravity !- Still bottom-to-top
          field2d % array(k,i) = field2d_src%array(k,i) * kgkg_kgm2 
        enddo
        enddo
!        write(*,*) 'MIN/MAX of index_qc.converted=',minval(field2d % array),maxval(field2d % array)
        !field2d % array(:,1:ngrid) = r2d_ptr_a(:,1:ngrid) ! TODO: [kg/kg] -> [kg/m2]
                                              ! multiply kgkg_kgm2=(atmosphere(1)%level_pressure(k)-atmosphere(1)%level_pressure(k-1))*r100/grav
                                              ! see gsi/crtm_interface.f90 or wrfda/da_get_innov_vector_crtm.inc
        field2d % fieldName = var_clw
        call mpas_pool_add_field(convFields, var_clw, field2d)
!        write(*,*) "end-of ",var_clw

     case ("atmosphere_mass_content_of_cloud_ice") !-var_cli 
        call mpas_pool_get_field(subFields, 'index_qi', field2d_src) !- [kg/kg]
        call mpas_duplicate_field(field2d_src, field2d)
        !--TODO: Trial: Should already have "var_prsi"
        call mpas_pool_get_array(convFields, "air_pressure_levels", r2d_ptr_b) !- [hPa]
        do i=1,ngrid
        do k=1,geom % nVertLevels
          kgkg_kgm2=( r2d_ptr_b(k,i)-r2d_ptr_b(k+1,i) ) * 100.0_kind_real / gravity !- Still bottom-to-top
          field2d % array(k,i) = field2d_src%array(k,i) * kgkg_kgm2 
        enddo
        enddo
!        write(*,*) 'MIN/MAX of index_qi.converted=',minval(field2d % array),maxval(field2d % array)
        !field2d % array(:,1:ngrid) = r2d_ptr_a(:,1:ngrid) ! TODO: [kg/kg] -> [kg/m2]
        !                                      ! same as var_clw
        field2d % fieldName = var_cli
        call mpas_pool_add_field(convFields, var_cli, field2d)
!        write(*,*) "end-of ",var_cli

     case ("effective_radius_of_cloud_liquid_water_particle") !-var_clwefr :TODO: currently filled w/ default value
        call mpas_pool_get_field(subFields, 're_cloud', field2d_src) !- [m]
        call mpas_duplicate_field(field2d_src, field2d)
        field2d % array(:,1:ngrid) = 10.0_kind_real !field2d_src%array(:,1:ngrid) * 1.0e-6 ! [m] -> [micron]
        field2d % fieldName = var_clwefr
        call mpas_pool_add_field(convFields, var_clwefr, field2d)
!        write(*,*) "end-of ",var_clwefr

     case ("effective_radius_of_cloud_ice_particle") !-var_cliefr :TODO: currently filled w/ default value
        call mpas_pool_get_field(subFields, 're_ice', field2d_src) !- [m]
        call mpas_duplicate_field(field2d_src, field2d)
        field2d % array(:,1:ngrid) = 30.0_kind_real !field2d_src%array(:,1:ngrid) * 1.0e-6 ! [m] -> [micron]
        field2d % fieldName = var_cliefr
        call mpas_pool_add_field(convFields, var_cliefr, field2d)
!        write(*,*) "end-of ",var_cliefr

     case ("Water_Temperature", "Land_Temperature", "Ice_Temperature", "Snow_Temperature" )

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

     case ("Snow_Depth")
        call mpas_pool_get_array(subFields, "snowh", r1d_ptr_a)
!        write(*,*) 'MIN/MAX of snowh=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
        call mpas_pool_get_field(subFields, 'u10', field1d_src) ! as a dummy array
        call mpas_duplicate_field(field1d_src, field1d)
        field1d % array(1:ngrid) = r1d_ptr_a(1:ngrid) * 1000.0_kind_real ! [m] -> [mm]
        field1d % fieldName = var_sfc_sdepth
        call mpas_pool_add_field(convFields, var_sfc_sdepth, field1d)
!        write(*,*) "end-of ",var_sfc_sdepth

     case ("Vegetation_Fraction")
        call mpas_pool_get_array(subFields, "vegfra", r1d_ptr_a)
!        write(*,*) 'MIN/MAX of vegfra=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
        call mpas_pool_get_field(subFields, 'u10', field1d_src) ! as a dummy array
        call mpas_duplicate_field(field1d_src, field1d)
        field1d % array(1:ngrid) = r1d_ptr_a(1:ngrid) / 100.0_kind_real ! [unitless, 0~100] = [%, 0~1]
        field1d % fieldName = var_sfc_vegfrac
        call mpas_pool_add_field(convFields, var_sfc_vegfrac, field1d)
!        write(*,*) "end-of ",var_sfc_vegfrac

     case ("Lai") !-var_sfc_lai :TODO, has a value for restart file, not for init file
        call mpas_pool_get_array(subFields, "lai", r1d_ptr_a)
!        write(*,*) 'MIN/MAX of lai=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
        call mpas_pool_get_field(subFields, 'u10', field1d_src) ! as a dummy array
        call mpas_duplicate_field(field1d_src, field1d)
        field1d % array(1:ngrid) = 3.5_kind_real !r1d_ptr_a(:) ! convert ?
        field1d % fieldName = var_sfc_lai
        call mpas_pool_add_field(convFields, var_sfc_lai, field1d)
!        write(*,*) "end-of ",var_sfc_lai

     case ("Soil_Moisture") !-var_sfc_soilm : NOTE: use 1st level
        call mpas_pool_get_array(subFields, "smois", r2d_ptr_a)
!        write(*,*) 'MIN/MAX of smois=',minval(r2d_ptr_a(1,:)),maxval(r2d_ptr_a(1,:))
        call mpas_pool_get_field(subFields, 'u10', field1d_src) ! as a dummy array
        call mpas_duplicate_field(field1d_src, field1d)
        field1d % array(1:ngrid) = 0.05_kind_real !r2d_ptr_a(1,:) ! [m3/m3] -> [g/cm3] :TODO ~ Good enough ?? range [0,1]
        field1d % fieldName = var_sfc_soilm
        call mpas_pool_add_field(convFields, var_sfc_soilm, field1d)
!        write(*,*) "end-of ",var_sfc_soilm

     case ("Soil_Temperature") !-var_sfc_soilt : NOTE: use 1st level
        call mpas_pool_get_array(subFields, "tslb", r2d_ptr_a)
!        write(*,*) 'MIN/MAX of tslb=',minval(r2d_ptr_a(1,:)),maxval(r2d_ptr_a(1,:))
        call mpas_pool_get_field(subFields, 'u10', field1d_src) ! as a dummy array
        call mpas_duplicate_field(field1d_src, field1d)
        field1d % array(1:ngrid) = r2d_ptr_a(1,1:ngrid) ! [K] -> [K]
        field1d % fieldName = var_sfc_soilt
        call mpas_pool_add_field(convFields, var_sfc_soilt, field1d)
!        write(*,*) "end-of ",var_sfc_soilt

     case ("geopotential_height")  !-var_z  geopotential heights at midpoint

        call mpas_pool_get_field(subFields, 'theta', field2d_src) ! as a dummy array
        call mpas_duplicate_field(field2d_src, field2d)
        call mpas_duplicate_field(field2d_src, field2d_a)

!       calculate midpoint geometricZ (unit: m):
        call geometricZ_full_to_half(geom%zgrid(:,1:ngrid), ngrid, &
                                     geom % nVertLevels,field2d_a%array(:,1:ngrid))

        do i=1,ngrid
           lat = geom%latCell(i) / deg2rad !- to Degrees
           do k=1,geom % nVertLevels
              call geometric2geop(lat,field2d_a%array(k,i), field2d%array(k,i))
           enddo
        enddo

        field2d % fieldName = var_z
        call mpas_pool_add_field(convFields, var_z, field2d)
        call mpas_deallocate_field(field2d_a) ! not used

     case ("sfc_geopotential_height")  !-var_sfc_z
        call mpas_pool_get_field(subFields, 'u10', field1d_src) ! as a dummy array
        call mpas_duplicate_field(field1d_src, field1d)

        do i=1,ngrid
           lat = geom%latCell(i) / deg2rad !- to Degrees
           call geometric2geop(lat,geom%zgrid(1,i), field1d%array(i))
        enddo
        field1d % fieldName = var_sfc_z
        call mpas_pool_add_field(convFields, var_sfc_z, field1d)

     case default
        write(*,*) 'Not processed in sub. convert_mpas_field2ufo: ',trim(fieldname(ivar))

     end select

   end do !ivar
  
end subroutine convert_mpas_field2ufo

!-------------------------------------------------------------------------------------------

subroutine convert_mpas_field2ufoTL(trajFields, subFields_tl, convFields_tl, fieldname, nfield, ngrid)

   implicit none

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

   type (field2DReal), pointer :: field2d, field2d_src, field2d_a, traj_field2d_a
   integer :: ivar

   !--- create new pull for ufo_vars
   call mpas_pool_create_pool(convFields_tl)


   do ivar=1, nfield
     write(*,*) 'convert_mpas_field2ufoTL:inside do/select case, &
               & ivar, trim(fieldname(ivar))=',ivar,trim(fieldname(ivar))

     select case (trim(fieldname(ivar)))

     case ( "virtual_temperature" ) !-var_tv 

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

     case ( "air_temperature", "temperature" ) !-var_ts
        call mpas_pool_get_field(subFields_tl, 'temperature', field2d_src) !< get temperature
        call mpas_duplicate_field(field2d_src, field2d)!  as a dummy array 
        field2d % fieldName = trim(fieldname(ivar))
        call mpas_pool_add_field(convFields_tl, trim(fieldname(ivar)), field2d)

     case ( "eastward_wind" ) !-var_??? eastward_wind

        !get TL variable
        call mpas_pool_get_field(subFields_tl, 'uReconstructZonal', field2d_src)
        call mpas_duplicate_field(field2d_src, field2d)!  as a dummy array
!        write(*,*) 'MIN/MAX of TL zonal wind(in/out)=',minval(field2d % array),maxval(field2d % array)

        !NL: field2d % array = field2d % array
        !TL: field2d % array = field2d % array

        field2d % fieldName = trim(fieldname(ivar))

        call mpas_pool_add_field(convFields_tl, trim(fieldname(ivar)), field2d)

!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ( "northward_wind" ) !-var_??? northward_wind

        !get TL variable
        call mpas_pool_get_field(subFields_tl, 'uReconstructMeridional', field2d_src)
        call mpas_duplicate_field(field2d_src, field2d)!  as a dummy array
!        write(*,*) 'MIN/MAX of TL meridional wind(in/out)=',minval(field2d % array),maxval(field2d % array)

        !NL: field2d % array = field2d % array
        !TL: field2d % array = field2d % array

        field2d % fieldName = trim(fieldname(ivar))

        call mpas_pool_add_field(convFields_tl, trim(fieldname(ivar)), field2d)

!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ("atmosphere_ln_pressure_coordinate") !-var_prsl

        !BJ: DO WE NEED THIS? "Currently" UFO do not access to tl of var_prsl. but for general purpose ?!?!
        !BJ: Without this array, 3DVAR gives the same results.

     case ("humidity_mixing_ratio") !-var_mixr
        call mpas_pool_get_field(subFields_tl, 'spechum', field2d_src) !< get TL of specific_huuumidity
        call mpas_pool_get_array(trajFields,   'spechum', traj_r2d_a)  !< get linearization state
        call mpas_duplicate_field(field2d_src, field2d)

        call q_to_w_tl(field2d_src % array(:,1:ngrid), traj_r2d_a(:,1:ngrid), field2d % array(:,1:ngrid)) 
        where (traj_r2d_a(:,1:ngrid) <= 0.0_kind_real)
          field2d % array(:,1:ngrid) = 0.0_kind_real
        end where
        field2d % array(:,1:ngrid) = field2d % array(:,1:ngrid) * 1000.0_kind_real
        field2d % fieldName = var_mixr
        call mpas_pool_add_field(convFields_tl, var_mixr, field2d)
!        write(*,*) "end-of ",var_mixr

     case("specific_humidity")
        call mpas_pool_get_field(subFields_tl, "spechum", field2d_src)
        call mpas_duplicate_field(field2d_src, field2d)
        field2d % fieldName = trim(fieldname(ivar))
        call mpas_pool_add_field(convFields_tl, trim(fieldname(ivar)), field2d)
!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ("air_pressure") !-var_prs
!        call mpas_pool_get_field(subFields_tl, 'pressure', field2d)
!        field2d % array(:,1:ngrid) = field2d_src%array(:,1:ngrid)
!        field2d % fieldName = var_prs
!        call mpas_pool_add_field(convFields_tl, var_prs, field2d)
!        write(*,*) "end-of ",var_prs

     case ("air_pressure_levels")
     case ("mass_concentration_of_ozone_in_air")
     case ("mass_concentration_of_carbon_dioxide_in_air")
     case ("atmosphere_mass_content_of_cloud_liquid_water")
     case ("atmosphere_mass_content_of_cloud_ice")
     case ("effective_radius_of_cloud_liquid_water_particle")
     case ("effective_radius_of_cloud_ice_particle")

     case ("Water_Fraction")
!JJG: All non-TL trajectory model variables should be provided in trajFields
!JJG: All TL model variables should be provided in subFields_tl
!        call mpas_pool_get_array(subFields_tl, "", r1d_ptr_a)
     case ("Land_Fraction")
!        call mpas_pool_get_array(subFields_tl, "", r1d_ptr_a)
     case ("Ice_Fraction")
!        call mpas_pool_get_array(subFields_tl, "", r1d_ptr_a)
     case ("Snow_Fraction")
!        call mpas_pool_get_array(subFields_tl, "", r1d_ptr_a)
     case ("Water_Temperature")
!        call mpas_pool_get_array(subFields_tl, "skintemp", r1d_ptr_a)
     case ("Land_Temperature")
!        call mpas_pool_get_array(subFields_tl, "skintemp", r1d_ptr_a)
     case ("Ice_Temperature")
!        call mpas_pool_get_array(subFields_tl, "skintemp", r1d_ptr_a)
     case ("Snow_Temperature")
!        call mpas_pool_get_array(subFields_tl, "skintemp", r1d_ptr_a)

     case ("Snow_Depth")
!        call mpas_pool_get_array(subFields_tl, "snowh", r1d_ptr_a)
     case ("Vegetation_Fraction")
!        call mpas_pool_get_array(subFields_tl, "vegfra", r1d_ptr_a)
     case ("Sfc_Wind_Speed", "Sfc_Wind_Direction")
!        call mpas_pool_get_array(subFields_tl, "u10", r1d_ptr_a)
!        call mpas_pool_get_array(subFields_tl, "v10", r1d_ptr_b)
     case ("Lai")
!        call mpas_pool_get_array(subFields_tl, "lai", r1d_ptr_a)
     case ("Soil_Moisture")
!        call mpas_pool_get_array(subFields_tl, "smois", r1d_ptr_a)
     case ("Soil_Temperature")
!        call mpas_pool_get_array(subFields_tl, "tslb", r1d_ptr_a)
     case ("Land_Type_Index")
!        call mpas_pool_get_array(subFields_tl, "", r1d_ptr_a)
     case ("Vegetation_Type")
!        call mpas_pool_get_array(subFields_tl, "ivgtyp", r1d_ptr_a)
     case ("Soil_Type")
!        call mpas_pool_get_array(subFields_tl, "isltyp", r1d_ptr_a)

     end select

   end do !ivar

  
end subroutine convert_mpas_field2ufoTL

!-------------------------------------------------------------------------------------------

subroutine convert_mpas_field2ufoAD(trajFields, subFields_ad, convFields_ad, fieldname, nfield, ngrid)

   implicit none

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

   type (field2DReal), pointer :: field2d, field2d_src, field2d_a, traj_field2d_a
   integer :: ivar

   do ivar=1, nfield
     write(*,*) 'convert_mpas_field2ufoAD: inside do/select case, &
               & ivar, trim(fieldname(ivar))=',ivar,trim(fieldname(ivar))

     select case (trim(fieldname(ivar)))

     case ( "virtual_temperature" ) !-var_tv 

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
        field2d_a % array(:,1:ngrid) = 0.0_kind_real !initialize local var.
        call tw_to_tv_ad(r2d_ptr_a(:,1:ngrid), field2d_a % array(:,1:ngrid), &
                         traj_r2d_a(:,1:ngrid), traj_field2d_a % array(:,1:ngrid), &
                         field2d % array(:,1:ngrid) )
        call q_to_w_ad( r2d_ptr_b(:,1:ngrid), traj_r2d_b(:,1:ngrid), field2d_a % array(:,1:ngrid) )
        call mpas_deallocate_field(field2d_a) ! not used
        call mpas_deallocate_field(traj_field2d_a) ! not used
!        write(*,*) 'MIN/MAX of AD temperature(out)=',minval(r2d_ptr_a),maxval(r2d_ptr_a)
!        write(*,*) 'MIN/MAX of AD     spechum(out)=',minval(r2d_ptr_b),maxval(r2d_ptr_b)

!        write(*,*) "end-of ",var_tv

     case ( "air_temperature", "temperature" ) !-var_ts
        call mpas_pool_get_array(subFields_ad, 'temperature', r2d_ptr_a) !< get temperature
        call mpas_pool_get_field(convFields_ad, trim(fieldname(ivar)), field2d)
        r2d_ptr_a(:,1:ngrid) = r2d_ptr_a(:,1:ngrid) + &
                         field2d % array(:,1:ngrid)

     case ( "eastward_wind" ) !-var_??? eastward_wind

        !get AD variable
        call mpas_pool_get_array(subFields_ad, 'uReconstructZonal', r2d_ptr_a) !< get zonal wind
!        write(*,*) 'MIN/MAX of AD zonal wind(in)  =',minval(r2d_ptr_a),maxval(r2d_ptr_a)

        !NL: field2d % array = field2d % array
        !TL: field2d % array = field2d % array
        call mpas_pool_get_field(convFields_ad, trim(fieldname(ivar)), field2d)

        r2d_ptr_a(:,1:ngrid) = r2d_ptr_a(:,1:ngrid) + &
                         field2d % array(:,1:ngrid)

!        write(*,*) 'MIN/MAX of AD zonal wind(out) =',minval(r2d_ptr_a),maxval(r2d_ptr_a)

!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ( "northward_wind" ) !-var_??? northward_wind

        !get AD variable
        call mpas_pool_get_array(subFields_ad, 'uReconstructMeridional', r2d_ptr_a) !< get meridional wind
!        write(*,*) 'MIN/MAX of AD meridional wind(in)  =',minval(r2d_ptr_a),maxval(r2d_ptr_a)

        !NL: field2d % array = field2d % array
        !TL: field2d % array = field2d % array
        call mpas_pool_get_field(convFields_ad, trim(fieldname(ivar)), field2d)

        r2d_ptr_a(:,1:ngrid) = r2d_ptr_a(:,1:ngrid) + &
                         field2d % array(:,1:ngrid)

!        write(*,*) 'MIN/MAX of AD meridional wind(out) =',minval(r2d_ptr_a),maxval(r2d_ptr_a)

!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ("atmosphere_ln_pressure_coordinate") !-var_prsl
        !BJ: DO WE NEED THIS? "Currently" UFO do not access to ad of var_prsl. but for general purpose ?!?!
        !BJ: Without this array, 3DVAR gives the same results.

     case ("humidity_mixing_ratio") !-var_mixr
        call mpas_pool_get_array(subFields_ad, "spechum", r2d_ptr_a) !< get sphechum_ad
        call mpas_pool_get_array(trajFields,   'spechum', traj_r2d_a)  !< get linearization state
!        write(*,*) 'MIN/MAX of AD spechum(in) =',minval(r2d_ptr_a(:,1:ngrid)),maxval(r2d_ptr_a(:,1:ngrid))

        call mpas_pool_get_field(convFields_ad, var_mixr, field2d)
!        write(*,*) 'MIN/MAX of AD var_mixr =',minval(field2d % array(:,1:ngrid)),maxval(field2d % array(:,1:ngrid))

        field2d % array(:,1:ngrid) = field2d % array(:,1:ngrid) * 1000.0_kind_real
        call q_to_w_ad(r2d_ptr_a(:,1:ngrid), traj_r2d_a(:,1:ngrid), field2d % array(:,1:ngrid))
        where (traj_r2d_a(:,1:ngrid) <= 0.0_kind_real)
          r2d_ptr_a(:,1:ngrid) = 0.0_kind_real
        end where
!        write(*,*) 'MIN/MAX of AD spechum(out) =',minval(r2d_ptr_a(:,1:ngrid)),maxval(r2d_ptr_a(:,1:ngrid))

!        write(*,*) "end-of ",var_mixr

     case("specific_humidity")
        call mpas_pool_get_array(subFields_ad, "spechum", r2d_ptr_a)
        call mpas_pool_get_field(convFields_ad, trim(fieldname(ivar)), field2d)
        r2d_ptr_a(:,1:ngrid) = r2d_ptr_a(:,1:ngrid) + field2d % array(:,1:ngrid)
!        write(*,*) "end-of ",trim(fieldname(ivar))

     case ("air_pressure") !-var_prs
!        call mpas_pool_get_array(subFields_ad, "pressure", r2d_ptr_a)
!        call mpas_pool_get_field(clone_subFields_ad, 'pressure', field2d) ! as a dummy array
!        field2d % array(:,1:ngrid) = r2d_ptr_a(:,1:ngrid)
!        field2d % fieldName = var_prs
!        call mpas_pool_add_field(convFields_ad, var_prs, field2d)
!        write(*,*) "end-of ",var_prs

     case ("air_pressure_levels")
     case ("mass_concentration_of_ozone_in_air")
     case ("mass_concentration_of_carbon_dioxide_in_air")
     case ("atmosphere_mass_content_of_cloud_liquid_water")
     case ("atmosphere_mass_content_of_cloud_ice")
     case ("effective_radius_of_cloud_liquid_water_particle")
     case ("effective_radius_of_cloud_ice_particle")

     case ("Water_Fraction")
!JJG: All non-AD trajectory model variables should be provided in trajFields
!JJG: All AD model variables should be provided in subFields_ad
!        call mpas_pool_get_array(subFields_ad, "", r1d_ptr_a)
     case ("Land_Fraction")
!        call mpas_pool_get_array(subFields_ad, "", r1d_ptr_a)
     case ("Ice_Fraction")
!        call mpas_pool_get_array(subFields_ad, "", r1d_ptr_a)
     case ("Snow_Fraction")
!        call mpas_pool_get_array(subFields_ad, "", r1d_ptr_a)
     case ("Water_Temperature")
!        call mpas_pool_get_array(subFields_ad, "skintemp", r1d_ptr_a)
     case ("Land_Temperature")
!        call mpas_pool_get_array(subFields_ad, "skintemp", r1d_ptr_a)
     case ("Ice_Temperature")
!        call mpas_pool_get_array(subFields_ad, "skintemp", r1d_ptr_a)
     case ("Snow_Temperature")
!        call mpas_pool_get_array(subFields_ad, "skintemp", r1d_ptr_a)

     case ("Snow_Depth")
!        call mpas_pool_get_array(subFields_ad, "snowh", r1d_ptr_a)
     case ("Vegetation_Fraction")
!        call mpas_pool_get_array(subFields_ad, "vegfra", r1d_ptr_a)
     case ("Sfc_Wind_Speed", "Sfc_Wind_Direction")
!        call mpas_pool_get_array(subFields_ad, "u10", r1d_ptr_a)
!        call mpas_pool_get_array(subFields_ad, "v10", r1d_ptr_b)
     case ("Lai")
!        call mpas_pool_get_array(subFields_ad, "lai", r1d_ptr_a)
     case ("Soil_Moisture")
!        call mpas_pool_get_array(subFields_ad, "smois", r1d_ptr_a)
     case ("Soil_Temperature")
!        call mpas_pool_get_array(subFields_ad, "tslb", r1d_ptr_a)
     case ("Land_Type_Index")
!        call mpas_pool_get_array(subFields_ad, "", r1d_ptr_a)
     case ("Vegetation_Type")
!        call mpas_pool_get_array(subFields_ad, "ivgtyp", r1d_ptr_a)
     case ("Soil_Type")
!        call mpas_pool_get_array(subFields_ad, "isltyp", r1d_ptr_a)

     end select

   end do !ivar

   !call mpas_pool_destroy_pool(clone_subFields_ad)
  
end subroutine convert_mpas_field2ufoAD

!-----------------------------------------------------------------------------------

   end module mpas2ufo_vars_mod
