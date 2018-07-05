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

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_abort, only : mpas_dmpar_global_abort
   use atm_core
   use mpi 
   use mpas_constants
   use mpas_kinds, only : kind_real


   use ufo_vars_mod !, only : var_tv, var_prsl

   private

   public :: update_mpas_field
   public :: convert_mpas_field2ufo,   &
             convert_mpas_field2ufoTL, &
             convert_mpas_field2ufoAD

   public :: convert_type_soil, convert_type_veg
   public :: uv_to_wdir

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
          write(*,*)'Looking for ',trim(poolItr % memberName)
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
                  write(0,*)'Update MIN/MAX: ',trim(fieldname(ii)),minval(r2d_ptr_a),maxval(r2d_ptr_a)
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

   !---- pool_a : self % subFields
   !---- pool_b : self % auxFields
   !---- pool_c : pool with UFO vars
   subroutine convert_mpas_field2ufo(geom, pool_a, pool_b, pool_c, fieldname, nfield)

   use mpas_geom_mod
   use mpas_constants, only : gravity

   implicit none

   type(mpas_geom), intent(in)  :: geom
   type (mpas_pool_type), pointer, intent(in) :: pool_a, pool_b ! subFields, auxFields
   type (mpas_pool_type), pointer, intent(out) :: pool_c ! pool with UFO vars
   integer, intent(in) :: nfield
   character (len=*), intent(in) :: fieldname(:) ! ufo

   type (mpas_pool_iterator_type) :: poolItr
   type (mpas_pool_type), pointer :: allFields
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b, r2d_ptr_c
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   integer, dimension(:), pointer :: i1d_ptr_a, i1d_ptr_b

   type (field1DInteger), pointer :: field1di, field1di_src
   type (field1DReal), pointer :: field1d, field1d_src
   type (field2DReal), pointer :: field2d, field2d_a, field2d_b, field2d_c, field2d_src
   integer :: ii, ivar

   !-- for var_prsi
   integer :: i, k, its, ite, kts, kte
   real (kind=kind_real) :: tem1, z0, z1, z2, w1, w2
   real (kind=kind_real), dimension(:,:), allocatable :: fzm_p, fzp_p

   real (kind=kind_real) :: kgkg_kgm2 !-- for var_clw, var_cli

   logical, SAVE :: l_detsfctemp = .false.

   !--- create new pull for ufo_vars
   call mpas_pool_create_pool(pool_c, nfield)

   !--- ufo_vars can be found in subFields or auxFields

   do ivar=1, nfield
     write(*,*) 'convert_mpas_field2ufo  :inside do/select case, ivar, trim(fieldname(ivar))=',ivar,trim(fieldname(ivar))

     select case (trim(fieldname(ivar)))

     case ( "virtual_temperature" ) !-var_tv 
        !get    theta from subFields
        !get       qv from subfields
        !get pressure from auxFields
        !do the conversion : need 'nl', 'tl', 'ad' ?!

        call mpas_pool_get_array(pool_a,    'theta', r2d_ptr_a)
        call mpas_pool_get_array(pool_a, 'index_qv', r2d_ptr_b)
        call mpas_pool_get_array(pool_b, 'pressure', r2d_ptr_c)
        write(*,*) 'MIN/MAX of    theta=',minval(r2d_ptr_a),maxval(r2d_ptr_a)
        write(*,*) 'MIN/MAX of index_qv=',minval(r2d_ptr_b),maxval(r2d_ptr_b)
        write(*,*) 'MIN/MAX of pressure=',minval(r2d_ptr_c),maxval(r2d_ptr_c)

        call mpas_pool_get_field(pool_a, 'theta', field2d_src)
        call mpas_duplicate_field(field2d_src, field2d)!  as a dummy array 

!%%%  NL 
!%%%  TODO: use clean formula
        field2d % array(:,:) = r2d_ptr_a(:,:) / ( (100000.0/r2d_ptr_c(:,:))**(287.05/1005.7) ) * &
                               (1.0 + (461.50/287.05 - 1.0)*r2d_ptr_b(:,:))
             !Tv = T * ( 1.0 + (rv/rd – 1)*qv), rv=461.50 , rd=287.05
                  !T = theta / ( (p0/pressure)**rd_over_cp) : to sensible temperature	
        write(*,*) 'MIN/MAX of Tv=',minval(field2d % array(:,:)),maxval(field2d % array(:,:))

        field2d % fieldName = var_tv

        call mpas_pool_add_field(pool_c, var_tv, field2d)

        write(*,*) "end-of ",var_tv

     case ("atmosphere_ln_pressure_coordinate") !-var_prsl
        !get pressure from subField
        !do the conversion : need 'nl', 'tl', 'ad' ?!

        call mpas_pool_get_array(pool_b, 'pressure', r2d_ptr_a)
        write(*,*) 'MIN/MAX of pressure=',minval(r2d_ptr_a),maxval(r2d_ptr_a)

        call mpas_pool_get_field(pool_a, 'rho', field2d_src)
        call mpas_duplicate_field(field2d_src, field2d)

        field2d % array(:,:) = log( r2d_ptr_a(:,:) / 100.0_kind_real / 10.0_kind_real ) !- Pa -> hPa ->cb
        write(*,*) 'MIN/MAX of ln_p=',minval(field2d % array(:,:)),maxval(field2d % array(:,:))

        field2d % fieldName = var_prsl

        call mpas_pool_add_field(pool_c, var_prsl, field2d)

        write(*,*) "end-of ",var_prsl

     case ("humidity_mixing_ratio") !-var_mixr
        call mpas_pool_get_array(pool_a, "index_qv", r2d_ptr_a)
        call mpas_pool_get_field(pool_a, 'theta', field2d_src) ! as a dummy array
        call mpas_duplicate_field(field2d_src, field2d)
        field2d % array(:,:) = r2d_ptr_a(:,:) * 1000.0_kind_real ! [kg/kg] -> [g/kg]
        field2d % fieldName = var_mixr
        call mpas_pool_add_field(pool_c, var_mixr, field2d)
        write(*,*) "end-of ",var_mixr

     case ("air_pressure") !-var_prs
        call mpas_pool_get_array(pool_b, "pressure", r2d_ptr_a)
        call mpas_pool_get_field(pool_a, 'theta', field2d_src) ! as a dummy array
        call mpas_duplicate_field(field2d_src, field2d)
        field2d % array(:,:) = r2d_ptr_a(:,:) / 100.0_kind_real ! [Pa] -> [hPa]
        field2d % fieldName = var_prs
        call mpas_pool_add_field(pool_c, var_prs, field2d)
        write(*,*) "end-of ",var_prs

     case ("air_pressure_levels") !-var_prsi
        call mpas_pool_get_array(pool_b, "pressure", r2d_ptr_a)
        call mpas_pool_get_field(pool_b, 'w', field2d_src) ! as a dummy array
        call mpas_duplicate_field(field2d_src, field2d)

        !-- ~/libs/MPAS-Release/src/core_atmosphere/physics/mpas_atmphys_manager.F   >> dimension Line 644.
        !-- ~/libs/MPAS-Release/src/core_atmosphere/physics/mpas_atmphys_interface.F >> formula   Line 365.
        !-- ~/libs/MPAS-Release/src/core_atmosphere/physics/mpas_atmphys_vars.F      >> declarations
        ! This routine needs to access GEOM (dimension, zgrid )
        ! TODO: Check: ite = nCells??  nCellsSolve???  MPI consideration ??? 
        its=1 ; ite = geom % nCells
        kts=1 ; kte = geom % nVertLevels
        allocate( fzm_p(its:ite,kts:kte), fzp_p(its:ite,kts:kte) )

        do k = kts+1,kte
        do i = its,ite
          tem1 = 1.0_kind_real/(geom % zgrid(k+1,i)-geom % zgrid(k-1,i))
          fzm_p(i,k) = (geom % zgrid(k,i)-geom % zgrid(k-1,i)) * tem1
          fzp_p(i,k) = (geom % zgrid(k+1,i)-geom % zgrid(k,i)) * tem1
          field2d % array(k,i) = fzm_p(i,k)*r2d_ptr_a(k,i) + fzp_p(i,k)*r2d_ptr_a(k-1,i)
        enddo
        enddo
        k = kte+1
        do i = its,ite
          z0 = geom % zgrid(k,i)
          z1 = 0.5_kind_real*(geom % zgrid(k,i)+geom % zgrid(k-1,i))
          z2 = 0.5_kind_real*(geom % zgrid(k-1,i)+geom % zgrid(k-2,i))
          w1 = (z0-z2)/(z1-z2)
          w2 = 1.0_kind_real-w1
          !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
          field2d % array(k,i) = exp( w1*log(r2d_ptr_a(k-1,i)) + w2*log(r2d_ptr_a(k-1,i)) )
        enddo
        k = kts
        do i = its,ite
          z0 = geom % zgrid(k,i)
          z1 = 0.5_kind_real*(geom % zgrid(k,i)+geom % zgrid(k+1,i))
          z2 = 0.5_kind_real*(geom % zgrid(k+1,i)+geom % zgrid(k+2,i))
          w1 = (z0-z2)/(z1-z2)
          w2 = 1.0_kind_real-w1
          field2d % array(k,i) = w1*r2d_ptr_a(k,i) + w2*r2d_ptr_a(k+1,i)
        enddo

        deallocate( fzm_p, fzp_p )

        field2d % array = field2d % array / 100.0_kind_real ! [Pa] -> [hPa]
        write(*,*) 'MIN/MAX of prsi=',minval(field2d % array),maxval(field2d % array)
        write(*,*) 'test prs       =',r2d_ptr_a(:,1)
        write(*,*) 'test prsi      =',field2d % array(:,1)

        field2d % fieldName = var_prsi
        call mpas_pool_add_field(pool_c, var_prsi, field2d)
        write(*,*) "end-of ",var_prsi

     case ("mass_concentration_of_ozone_in_air") !-var_oz :TODO: not directly available from MPAS
        call mpas_pool_get_array(pool_a, "theta", r2d_ptr_a)
        call mpas_pool_get_field(pool_a, 'theta', field2d_src) ! as a dummy array
        call mpas_duplicate_field(field2d_src, field2d)
        field2d % array(:,:) = 0.0_kind_real !r2d_ptr_a(:,:) ! convert ??
        field2d % fieldName = var_oz
        call mpas_pool_add_field(pool_c, var_oz, field2d)
        write(*,*) "end-of ",var_oz

     case ("mass_concentration_of_carbon_dioxide_in_air") !-var_co2 :TODO: not directly available from MPAS
        call mpas_pool_get_array(pool_a, "theta", r2d_ptr_a)
        call mpas_pool_get_field(pool_a, 'theta', field2d_src) ! as a dummy array
        call mpas_duplicate_field(field2d_src, field2d)
        field2d % array(:,:) = 0.0_kind_real !r2d_ptr_a(:,:) ! convert ??
        field2d % fieldName = var_co2
        call mpas_pool_add_field(pool_c, var_co2, field2d)
        write(*,*) "end-of ",var_co2

     case ("atmosphere_mass_content_of_cloud_liquid_water") !-var_clw 
        call mpas_pool_get_array(pool_b, "index_qc", r2d_ptr_a) !- [kg/kg] 
        write(*,*) 'MIN/MAX of index_qc=',minval(r2d_ptr_a),maxval(r2d_ptr_a)
        call mpas_pool_get_field(pool_a, 'theta', field2d_src) ! as a dummy array
        call mpas_duplicate_field(field2d_src, field2d)
        !--TODO: Trial: Should already have "var_prsi"
        call mpas_pool_get_array(pool_c, "air_pressure_levels", r2d_ptr_b) !- [hPa]
        do i=1,geom % nCells
        do k=1,geom % nVertLevels
          kgkg_kgm2=( r2d_ptr_b(k,i)-r2d_ptr_b(k+1,i) ) * 100.0_kind_real / gravity !- Still bottom-to-top
          field2d % array(k,i) = r2d_ptr_a(k,i) * kgkg_kgm2 
        enddo
        enddo
        write(*,*) 'MIN/MAX of index_qc.converted=',minval(field2d % array),maxval(field2d % array)
        !field2d % array(:,:) = r2d_ptr_a(:,:) ! TODO: [kg/kg] -> [kg/m2]
                                              ! multiply kgkg_kgm2=(atmosphere(1)%level_pressure(k)-atmosphere(1)%level_pressure(k-1))*r100/grav
                                              ! see gsi/crtm_interface.f90 or wrfda/da_get_innov_vector_crtm.inc
        field2d % fieldName = var_clw
        call mpas_pool_add_field(pool_c, var_clw, field2d)
        write(*,*) "end-of ",var_clw

     case ("atmosphere_mass_content_of_cloud_ice") !-var_cli 
        call mpas_pool_get_array(pool_b, "index_qi", r2d_ptr_a) !- [kg/kg] 
        write(*,*) 'MIN/MAX of index_qi=',minval(r2d_ptr_a),maxval(r2d_ptr_a)
        call mpas_pool_get_field(pool_a, 'theta', field2d_src) ! as a dummy array
        call mpas_duplicate_field(field2d_src, field2d)
        !--TODO: Trial: Should already have "var_prsi"
        call mpas_pool_get_array(pool_c, "air_pressure_levels", r2d_ptr_b) !- [hPa]
        do i=1,geom % nCells
        do k=1,geom % nVertLevels
          kgkg_kgm2=( r2d_ptr_b(k,i)-r2d_ptr_b(k+1,i) ) * 100.0_kind_real / gravity !- Still bottom-to-top
          field2d % array(k,i) = r2d_ptr_a(k,i) * kgkg_kgm2 
        enddo
        enddo
        write(*,*) 'MIN/MAX of index_qi.converted=',minval(field2d % array),maxval(field2d % array)
        !field2d % array(:,:) = r2d_ptr_a(:,:) ! TODO: [kg/kg] -> [kg/m2]
        !                                      ! same as var_clw
        field2d % fieldName = var_cli
        call mpas_pool_add_field(pool_c, var_cli, field2d)
        write(*,*) "end-of ",var_cli

     case ("effective_radius_of_cloud_liquid_water_particle") !-var_clwefr :TODO: currently filled w/ default value
        call mpas_pool_get_array(pool_b, "re_cloud", r2d_ptr_a) !- [m]
        write(*,*) 'MIN/MAX of re_cloud=',minval(r2d_ptr_a),maxval(r2d_ptr_a)
        call mpas_pool_get_field(pool_a, 'theta', field2d_src) ! as a dummy array
        call mpas_duplicate_field(field2d_src, field2d)
        field2d % array(:,:) = 10.0_kind_real !r2d_ptr_a(:,:) * 1.0e-6 ! [m] -> [micron]
        field2d % fieldName = var_clwefr
        call mpas_pool_add_field(pool_c, var_clwefr, field2d)
        write(*,*) "end-of ",var_clwefr

     case ("effective_radius_of_cloud_ice_particle") !-var_cliefr :TODO: currently filled w/ default value
        call mpas_pool_get_array(pool_b, "re_ice", r2d_ptr_a) !- [m]
        write(*,*) 'MIN/MAX of re_ice=',minval(r2d_ptr_a),maxval(r2d_ptr_a)
        call mpas_pool_get_field(pool_a, 'theta', field2d_src) ! as a dummy array
        call mpas_duplicate_field(field2d_src, field2d)
        field2d % array(:,:) = 30.0_kind_real !r2d_ptr_a(:,:) * 1.0e-6 ! [m] -> [micron]
        field2d % fieldName = var_cliefr
        call mpas_pool_add_field(pool_c, var_cliefr, field2d)
        write(*,*) "end-of ",var_cliefr

     case ("Water_Temperature", "Land_Temperature", "Ice_Temperature", "Snow_Temperature" )

        if( .not. l_detsfctemp) then
          call mpas_pool_get_field(pool_b, 'u10', field1d_src) ! as a dummy array

          !--- currently assign "skintemp" for all temperature, TODO: more proper variable?
          call mpas_pool_get_array(pool_b, "skintemp", r1d_ptr_a) !"ground or water surface temperature"
          write(*,*) 'MIN/MAX of skintemp=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
          call mpas_duplicate_field(field1d_src, field1d)
          field1d % array(:) = r1d_ptr_a(:) ! quantity and unit might change
          field1d % fieldName = var_sfc_wtmp
          call mpas_pool_add_field(pool_c, var_sfc_wtmp, field1d)

          call mpas_duplicate_field(field1d_src, field1d)
          field1d % array(:) = r1d_ptr_a(:) ! quantity and unit might change
          field1d % fieldName = var_sfc_ltmp
          call mpas_pool_add_field(pool_c, var_sfc_ltmp, field1d)

          call mpas_duplicate_field(field1d_src, field1d)
          field1d % array(:) = r1d_ptr_a(:) ! quantity and unit might change
          field1d % fieldName = var_sfc_itmp
          call mpas_pool_add_field(pool_c, var_sfc_itmp, field1d)

          call mpas_duplicate_field(field1d_src, field1d)
          field1d % array(:) = r1d_ptr_a(:) ! quantity and unit might change
          field1d % fieldName = var_sfc_stmp
          call mpas_pool_add_field(pool_c, var_sfc_stmp, field1d)

          l_detsfctemp = .true.
        endif

     case ("Snow_Depth")
        call mpas_pool_get_array(pool_b, "snowh", r1d_ptr_a)
        write(*,*) 'MIN/MAX of snowh=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
        call mpas_pool_get_field(pool_b, 'u10', field1d_src) ! as a dummy array
        call mpas_duplicate_field(field1d_src, field1d)
        field1d % array(:) = r1d_ptr_a(:) * 1000.0_kind_real ! [m] -> [mm]
        field1d % fieldName = var_sfc_sdepth
        call mpas_pool_add_field(pool_c, var_sfc_sdepth, field1d)
        write(*,*) "end-of ",var_sfc_sdepth

     case ("Vegetation_Fraction")
        call mpas_pool_get_array(pool_b, "vegfra", r1d_ptr_a)
        write(*,*) 'MIN/MAX of vegfra=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
        call mpas_pool_get_field(pool_b, 'u10', field1d_src) ! as a dummy array
        call mpas_duplicate_field(field1d_src, field1d)
        field1d % array(:) = r1d_ptr_a(:) / 100.0_kind_real ! [unitless, 0~100] = [%, 0~1]
        field1d % fieldName = var_sfc_vegfrac
        call mpas_pool_add_field(pool_c, var_sfc_vegfrac, field1d)
        write(*,*) "end-of ",var_sfc_vegfrac

     case ("Lai") !-var_sfc_lai :TODO, has a value for restart file, not for init file
        call mpas_pool_get_array(pool_b, "lai", r1d_ptr_a)
        write(*,*) 'MIN/MAX of lai=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
        call mpas_pool_get_field(pool_b, 'u10', field1d_src) ! as a dummy array
        call mpas_duplicate_field(field1d_src, field1d)
        field1d % array(:) = 3.5_kind_real !r1d_ptr_a(:) ! convert ?
        field1d % fieldName = var_sfc_lai
        call mpas_pool_add_field(pool_c, var_sfc_lai, field1d)
        write(*,*) "end-of ",var_sfc_lai

     case ("Soil_Moisture") !-var_sfc_soilm : NOTE: use 1st level
        call mpas_pool_get_array(pool_b, "smois", r2d_ptr_a)
        write(*,*) 'MIN/MAX of smois=',minval(r2d_ptr_a(1,:)),maxval(r2d_ptr_a(1,:))
        call mpas_pool_get_field(pool_b, 'u10', field1d_src) ! as a dummy array
        call mpas_duplicate_field(field1d_src, field1d)
        field1d % array(:) = 0.05_kind_real !r2d_ptr_a(1,:) ! [m3/m3] -> [g/cm3] :TODO ~ Good enough ?? range [0,1]
        field1d % fieldName = var_sfc_soilm
        call mpas_pool_add_field(pool_c, var_sfc_soilm, field1d)
        write(*,*) "end-of ",var_sfc_soilm

     case ("Soil_Temperature") !-var_sfc_soilt : NOTE: use 1st level
        call mpas_pool_get_array(pool_b, "tslb", r2d_ptr_a)
        write(*,*) 'MIN/MAX of tslb=',minval(r2d_ptr_a(1,:)),maxval(r2d_ptr_a(1,:))
        call mpas_pool_get_field(pool_b, 'u10', field1d_src) ! as a dummy array
        call mpas_duplicate_field(field1d_src, field1d)
        field1d % array(:) = r2d_ptr_a(1,:) ! [K] -> [K]
        field1d % fieldName = var_sfc_soilt
        call mpas_pool_add_field(pool_c, var_sfc_soilt, field1d)
        write(*,*) "end-of ",var_sfc_soilt

     case default
        write(*,*) 'Not processed in sub. convert_mpas_field2ufo: ',trim(fieldname(ivar))

     end select

   end do !ivar
  
   end subroutine convert_mpas_field2ufo

!-------------------------------------------------------------------------------------------

   !subroutine convert_mpas_field2ufoTL(pool_a, pool_b, pool_c, fieldname, nfield)
   !subroutine convert_mpas_field2ufoTL(traj_pool_a, traj_pool_b, pool_a, pool_b, pool_c, fieldname, nfield)
   subroutine convert_mpas_field2ufoTL(pool_traj, pool_a, pool_b, pool_c, fieldname, nfield)

   implicit none

   type (mpas_pool_type), pointer, intent(in) :: pool_traj !traj_pool_a, traj_pool_b ! subFields, auxFields
   type (mpas_pool_type), pointer, intent(in) :: pool_a, pool_b ! subFields, auxFields
   type (mpas_pool_type), pointer, intent(out) :: pool_c
   integer, intent(in) :: nfield
   character (len=*), intent(in) :: fieldname(:) ! ufo

   type (mpas_pool_iterator_type) :: poolItr
   type (mpas_pool_type), pointer :: allFields
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b, r2d_ptr_c
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: traj_r2d_a, traj_r2d_b, traj_r2d_c !BJJ test

   type (field2DReal), pointer :: field2d, field2d_a, field2d_b, field2d_c, field2d_src
   integer :: ii, ivar

   !--- create new pull for ufo_vars
   call mpas_pool_create_pool(pool_c, nfield)

   !--- ufo_vars can be found in subFields or auxFields

   do ivar=1, nfield
     write(*,*) 'convert_mpas_field2ufoTL:inside do/select case, ivar, trim(fieldname(ivar))=',ivar,trim(fieldname(ivar))

     select case (trim(fieldname(ivar)))

     case ( "virtual_temperature" ) !-var_tv 
        !get    theta from subFields
        !get       qv from subfields
        !get pressure from auxFields
        !do the conversion : need 'nl', 'tl', 'ad' ?!

        call mpas_pool_get_array(pool_a,    'theta', r2d_ptr_a)
        call mpas_pool_get_array(pool_a, 'index_qv', r2d_ptr_b)
        call mpas_pool_get_array(pool_b, 'pressure', r2d_ptr_c)
        write(*,*) 'MIN/MAX of    theta=',minval(r2d_ptr_a),maxval(r2d_ptr_a)
        write(*,*) 'MIN/MAX of index_qv=',minval(r2d_ptr_b),maxval(r2d_ptr_b)
        write(*,*) 'MIN/MAX of pressure=',minval(r2d_ptr_c),maxval(r2d_ptr_c)

        call mpas_pool_get_array(pool_traj,    'theta', traj_r2d_a)
        call mpas_pool_get_array(pool_traj, 'index_qv', traj_r2d_b)
        call mpas_pool_get_array(pool_traj, 'pressure', traj_r2d_c)
        write(*,*) 'MIN/MAX of TRAJ    theta=',minval(traj_r2d_a),maxval(traj_r2d_a)
        write(*,*) 'MIN/MAX of TRAJ index_qv=',minval(traj_r2d_b),maxval(traj_r2d_b)
        write(*,*) 'MIN/MAX of TRAJ pressure=',minval(traj_r2d_c),maxval(traj_r2d_c)

        call mpas_pool_get_field(pool_a, 'theta', field2d_src)
        call mpas_duplicate_field(field2d_src, field2d)

!%%%  TL with trajectory!
!%%%  TODO: use clean formula
        !field2d % array(:,:) = r2d_ptr_a(:,:) / ( (100000.0/r2d_ptr_c(:,:))**(287.05/1005.7) ) * &
        !                       (1.0 + (461.50/287.05 - 1.0)*r2d_ptr_b(:,:))
        !     !Tv = T * ( 1.0 + (rv/rd – 1)*qv), rv=461.50 , rd=287.05
        !          !T = theta / ( (p0/pressure)**rd_over_cp) : to sensible temperature	
        field2d % array(:,:) = ( ( 1.0 + (461.50/287.05 - 1.0)*traj_r2d_b(:,:) ) * r2d_ptr_a(:,:) &
                                + traj_r2d_a(:,:) * (461.50/287.05 - 1.0) * r2d_ptr_b(:,:) ) &
                               / ( (100000.0/traj_r2d_c(:,:))**(287.05/1005.7) )
        write(*,*) 'MIN/MAX of Tv=',minval(field2d % array(:,:)),maxval(field2d % array(:,:))

        field2d % fieldName = var_tv

        call mpas_pool_add_field(pool_c, var_tv, field2d)

        write(*,*) "end-of ",var_tv

     case ("atmosphere_ln_pressure_coordinate") !-var_prsl
        !get pressure from subField
        !do the conversion : need 'nl', 'tl', 'ad' ?!

        call mpas_pool_get_array(pool_b, 'pressure', r2d_ptr_a)
        write(*,*) 'MIN/MAX of pressure=',minval(r2d_ptr_a),maxval(r2d_ptr_a)

        call mpas_pool_get_field(pool_a, 'rho', field2d_src)
        call mpas_duplicate_field(field2d_src, field2d)  ! as a dummy array

        field2d % array(:,:) = log( r2d_ptr_a(:,:)/100./10. ) !- Pa -> hPa -> kPa
        write(*,*) 'MIN/MAX of ln_p=',minval(field2d % array(:,:)),maxval(field2d % array(:,:))

        field2d % fieldName = var_prsl

        call mpas_pool_add_field(pool_c, var_prsl, field2d)

        write(*,*) "end-of ",var_prsl

     case ("humidity_mixing_ratio") !-var_mixr
        !call mpas_pool_get_array(pool_a, "index_qv", r2d_ptr_a)
        !call mpas_pool_get_field(clone_pool_a, 'pressure', field2d) ! as a dummy array
        !field2d % array(:,:) = r2d_ptr_a(:,:)
        !field2d % fieldName = var_mixr
        !call mpas_pool_add_field(pool_c, var_mixr, field2d)
        !write(*,*) "end-of ",var_mixr

     case ("air_pressure") !-var_prs
        !call mpas_pool_get_array(pool_b, "pressure", r2d_ptr_a)
        !call mpas_pool_get_field(clone_pool_a, 'pressure', field2d) ! as a dummy array
        !field2d % array(:,:) = r2d_ptr_a(:,:)
        !field2d % fieldName = var_prs
        !call mpas_pool_add_field(pool_c, var_prs, field2d)
        !write(*,*) "end-of ",var_prs

     case ("air_pressure_levels")
     case ("mass_concentration_of_ozone_in_air")
     case ("mass_concentration_of_carbon_dioxide_in_air")
     case ("atmosphere_mass_content_of_cloud_liquid_water")
     case ("atmosphere_mass_content_of_cloud_ice")
     case ("effective_radius_of_cloud_liquid_water_particle")
     case ("effective_radius_of_cloud_ice_particle")

     case ("Water_Fraction")
        call mpas_pool_get_array(pool_b, "", r1d_ptr_a)
     case ("Land_Fraction")
        call mpas_pool_get_array(pool_b, "", r1d_ptr_a)
     case ("Ice_Fraction")
        call mpas_pool_get_array(pool_b, "", r1d_ptr_a)
     case ("Snow_Fraction")
        call mpas_pool_get_array(pool_b, "", r1d_ptr_a)
     case ("Water_Temperature")
        call mpas_pool_get_array(pool_b, "skintemp", r1d_ptr_a)
     case ("Land_Temperature")
        call mpas_pool_get_array(pool_b, "skintemp", r1d_ptr_a)
     case ("Ice_Temperature")
        call mpas_pool_get_array(pool_b, "skintemp", r1d_ptr_a)
     case ("Snow_Temperature")
        call mpas_pool_get_array(pool_b, "skintemp", r1d_ptr_a)

     case ("Snow_Depth")
        call mpas_pool_get_array(pool_b, "snowh", r1d_ptr_a)
     case ("Vegetation_Fraction")
        call mpas_pool_get_array(pool_b, "vegfra", r1d_ptr_a)
     case ("Sfc_Wind_Speed", "Sfc_Wind_Direction")
        call mpas_pool_get_array(pool_a, "u10", r1d_ptr_a)
        call mpas_pool_get_array(pool_a, "v10", r1d_ptr_b)
     case ("Lai")
        call mpas_pool_get_array(pool_b, "lai", r1d_ptr_a)
     case ("Soil_Moisture")
        call mpas_pool_get_array(pool_b, "smois", r1d_ptr_a)
     case ("Soil_Temperature")
        call mpas_pool_get_array(pool_b, "tslb", r1d_ptr_a)
     case ("Land_Type_Index")
        call mpas_pool_get_array(pool_b, "", r1d_ptr_a)
     case ("Vegetation_Type")
        call mpas_pool_get_array(pool_b, "ivgtyp", r1d_ptr_a)
     case ("Soil_Type")
        call mpas_pool_get_array(pool_b, "isltyp", r1d_ptr_a)

     end select

   end do !ivar

  
   end subroutine convert_mpas_field2ufoTL

!-------------------------------------------------------------------------------------------

   !subroutine convert_mpas_field2ufoAD(pool_a, pool_b, pool_c, fieldname, nfield)
   !subroutine convert_mpas_field2ufoAD(traj_pool_a, traj_pool_b, pool_a, pool_b, pool_c, fieldname, nfield)
   subroutine convert_mpas_field2ufoAD(pool_traj, pool_a, pool_b, pool_c, fieldname, nfield)

   implicit none

   type (mpas_pool_type), pointer, intent(in   ) :: pool_traj !traj_pool_a, traj_pool_b ! subFields, auxFields
   type (mpas_pool_type), pointer, intent(inout) :: pool_a, pool_b ! subFields, auxFields
   type (mpas_pool_type), pointer, intent(inout) :: pool_c
   integer, intent(in) :: nfield
   character (len=*), intent(in) :: fieldname(:) ! ufo

   type (mpas_pool_iterator_type) :: poolItr
   type (mpas_pool_type), pointer :: allFields, clone_pool_a
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b, r2d_ptr_c
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: traj_r2d_a, traj_r2d_b, traj_r2d_c !BJJ test

   type (field2DReal), pointer :: field2d, field2d_a, field2d_b, field2d_c
   integer :: ii, ivar
   !--- test
   !call mpas_pool_create_pool(clone_pool_a)
   !call mpas_pool_clone_pool(pool_a,clone_pool_a)

   !--- create new pull for ufo_vars
   !call mpas_pool_create_pool(pool_c, nfield)

   !--- ufo_vars can be found in subFields or auxFields

   do ivar=1, nfield
     write(*,*) 'convert_mpas_field2ufoAD: inside do/select case, ivar, trim(fieldname(ivar))=',ivar,trim(fieldname(ivar))

     select case (trim(fieldname(ivar)))

     case ( "virtual_temperature" ) !-var_tv 
        !get    theta from subFields
        !get       qv from subfields
        !get pressure from auxFields
        !do the conversion : need 'nl', 'tl', 'ad' ?!

        call mpas_pool_get_array(pool_a,    'theta', r2d_ptr_a)
        call mpas_pool_get_array(pool_a, 'index_qv', r2d_ptr_b)
        call mpas_pool_get_array(pool_b, 'pressure', r2d_ptr_c)
        write(*,*) 'MIN/MAX of    theta=',minval(r2d_ptr_a),maxval(r2d_ptr_a)
        write(*,*) 'MIN/MAX of index_qv=',minval(r2d_ptr_b),maxval(r2d_ptr_b)
        write(*,*) 'MIN/MAX of pressure=',minval(r2d_ptr_c),maxval(r2d_ptr_c)

        call mpas_pool_get_array(pool_traj,    'theta', traj_r2d_a)
        call mpas_pool_get_array(pool_traj, 'index_qv', traj_r2d_b)
        call mpas_pool_get_array(pool_traj, 'pressure', traj_r2d_c)
        write(*,*) 'MIN/MAX of TRAJ    theta=',minval(traj_r2d_a),maxval(traj_r2d_a)
        write(*,*) 'MIN/MAX of TRAJ index_qv=',minval(traj_r2d_b),maxval(traj_r2d_b)
        write(*,*) 'MIN/MAX of TRAJ pressure=',minval(traj_r2d_c),maxval(traj_r2d_c)

        !call mpas_pool_get_field(clone_pool_a, 'theta', field2d) ! as a dummy array

!%%%  AD with trajectory!
!%%%  TODO: use clean formula
        !field2d % array(:,:) = r2d_ptr_a(:,:) / ( (100000.0/r2d_ptr_c(:,:))**(287.05/1005.7) ) * &
        !                       (1.0 + (461.50/287.05 - 1.0)*r2d_ptr_b(:,:))
        !     !Tv = T * ( 1.0 + (rv/rd ?~@~S 1)*qv), rv=461.50 , rd=287.05
        !          !T = theta / ( (p0/pressure)**rd_over_cp) : to sensible temperature   
        !field2d % array(:,:) = ( ( 1.0 + (461.50/287.05 - 1.0)*traj_r2d_b(:,:) ) * r2d_ptr_a(:,:) &
        !                        + traj_r2d_a(:,:) * (461.50/287.05 - 1.0) * r2d_ptr_b(:,:) ) &
        !                       / ( (100000.0/traj_r2d_c(:,:))**(287.05/1005.7) )
        call mpas_pool_get_field(pool_c, 'virtual_temperature', field2d)
        write(*,*) 'MIN/MAX of Tv=',minval(field2d % array(:,:)),maxval(field2d % array(:,:))
        r2d_ptr_a(:,:)=0.0_kind_real
        r2d_ptr_b(:,:)=0.0_kind_real
        r2d_ptr_a(:,:) = r2d_ptr_a(:,:) + ( 1.0 + (461.50/287.05 - 1.0)*traj_r2d_b(:,:) ) / &
                         ( (100000.0/traj_r2d_c(:,:))**(287.05/1005.7) ) * field2d % array(:,:)
        r2d_ptr_b(:,:) = r2d_ptr_b(:,:) + traj_r2d_a(:,:) * (461.50/287.05 - 1.0) / &
                         ( (100000.0/traj_r2d_c(:,:))**(287.05/1005.7) ) * field2d % array(:,:)
        write(*,*) 'MIN/MAX of theta=',minval(r2d_ptr_a(:,:)),maxval(r2d_ptr_a(:,:))
        write(*,*) 'MIN/MAX of index_qv=',minval(r2d_ptr_b(:,:)),maxval(r2d_ptr_b(:,:))

        write(*,*) "end-of ",var_tv

     case ("atmosphere_ln_pressure_coordinate") !-var_prsl
        !get pressure from subField
        !do the conversion : need 'nl', 'tl', 'ad' ?!

!        call mpas_pool_get_array(pool_b, 'pressure', r2d_ptr_a)
!        write(*,*) 'MIN/MAX of pressure=',minval(r2d_ptr_a),maxval(r2d_ptr_a)
!
!        call mpas_pool_get_field(clone_pool_a, 'rho', field2d) ! as a dummy array
!
!        field2d % array(:,:) = log( r2d_ptr_a(:,:)/100./10. ) !- Pa -> hPa -> kPa
!        write(*,*) 'MIN/MAX of ln_p=',minval(field2d % array(:,:)),maxval(field2d % array(:,:))
!
!        field2d % fieldName = var_prsl
!
!        call mpas_pool_add_field(pool_c, var_prsl, field2d)
!
        write(*,*) "end-of ",var_prsl

     case ("humidity_mixing_ratio") !-var_mixr
        !call mpas_pool_get_array(pool_a, "index_qv", r2d_ptr_a)
        !call mpas_pool_get_field(clone_pool_a, 'pressure', field2d) ! as a dummy array
        !field2d % array(:,:) = r2d_ptr_a(:,:)
        !field2d % fieldName = var_mixr
        !call mpas_pool_add_field(pool_c, var_mixr, field2d)
        !write(*,*) "end-of ",var_mixr

     case ("air_pressure") !-var_prs
        !call mpas_pool_get_array(pool_b, "pressure", r2d_ptr_a)
        !call mpas_pool_get_field(clone_pool_a, 'pressure', field2d) ! as a dummy array
        !field2d % array(:,:) = r2d_ptr_a(:,:)
        !field2d % fieldName = var_prs
        !call mpas_pool_add_field(pool_c, var_prs, field2d)
        !write(*,*) "end-of ",var_prs

     case ("air_pressure_levels")
     case ("mass_concentration_of_ozone_in_air")
     case ("mass_concentration_of_carbon_dioxide_in_air")
     case ("atmosphere_mass_content_of_cloud_liquid_water")
     case ("atmosphere_mass_content_of_cloud_ice")
     case ("effective_radius_of_cloud_liquid_water_particle")
     case ("effective_radius_of_cloud_ice_particle")

     case ("Water_Fraction")
        call mpas_pool_get_array(pool_b, "", r1d_ptr_a)
     case ("Land_Fraction")
        call mpas_pool_get_array(pool_b, "", r1d_ptr_a)
     case ("Ice_Fraction")
        call mpas_pool_get_array(pool_b, "", r1d_ptr_a)
     case ("Snow_Fraction")
        call mpas_pool_get_array(pool_b, "", r1d_ptr_a)
     case ("Water_Temperature")
        call mpas_pool_get_array(pool_b, "skintemp", r1d_ptr_a)
     case ("Land_Temperature")
        call mpas_pool_get_array(pool_b, "skintemp", r1d_ptr_a)
     case ("Ice_Temperature")
        call mpas_pool_get_array(pool_b, "skintemp", r1d_ptr_a)
     case ("Snow_Temperature")
        call mpas_pool_get_array(pool_b, "skintemp", r1d_ptr_a)

     case ("Snow_Depth")
        call mpas_pool_get_array(pool_b, "snowh", r1d_ptr_a)
     case ("Vegetation_Fraction")
        call mpas_pool_get_array(pool_b, "vegfra", r1d_ptr_a)
     case ("Sfc_Wind_Speed", "Sfc_Wind_Direction")
        call mpas_pool_get_array(pool_a, "u10", r1d_ptr_a)
        call mpas_pool_get_array(pool_a, "v10", r1d_ptr_b)
     case ("Lai")
        call mpas_pool_get_array(pool_b, "lai", r1d_ptr_a)
     case ("Soil_Moisture")
        call mpas_pool_get_array(pool_b, "smois", r1d_ptr_a)
     case ("Soil_Temperature")
        call mpas_pool_get_array(pool_b, "tslb", r1d_ptr_a)
     case ("Land_Type_Index")
        call mpas_pool_get_array(pool_b, "", r1d_ptr_a)
     case ("Vegetation_Type")
        call mpas_pool_get_array(pool_b, "ivgtyp", r1d_ptr_a)
     case ("Soil_Type")
        call mpas_pool_get_array(pool_b, "isltyp", r1d_ptr_a)

     end select

   end do !ivar


   !call mpas_pool_empty_pool(clone_pool_a)
   !call mpas_pool_destroy_pool(clone_pool_a)
  
   end subroutine convert_mpas_field2ufoAD

!-----------------------------------------------------------------------------------

   end module mpas2ufo_vars_mod
