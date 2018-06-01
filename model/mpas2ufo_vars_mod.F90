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
   use mpas_kinds, only : kind_real


   use ufo_vars_mod !, only : var_tv, var_prsl

!   private

   public :: update_mpas_field
   public :: convert_mpas_field2ufo,   &
             convert_mpas_field2ufoTL, &
             convert_mpas_field2ufoAD

   public :: wrf_to_crtm_soil, usgs_to_crtm_mw

   !-- from WRFDA/da_crtm.f90
   integer, parameter :: n_soil_type = 16  ! wrf num_soil_cat
   integer, parameter :: USGS_n_type = 24 
   integer, parameter :: IGBP_n_type = 20 
   integer, parameter :: wrf_to_crtm_soil(n_soil_type) = &
      (/ 1, 1, 4, 2, 2, 8, 7, 2, 6, 5, 2, 3, 8, 1, 6, 9 /)
   ! vegetation type mapping for GFS classification scheme
   ! REL-2.1.3.CRTM_User_Guide.pdf table 4.16
   integer, parameter :: usgs_to_crtm_mw(USGS_n_type) = &
      (/  7, 12, 12, 12, 12, 12,  7,  9,  8,  6, &
          2,  5,  1,  4,  3,  0,  8,  8, 11, 10, &
         10, 10, 11, 13 /)  
   integer, parameter :: igbp_to_crtm_mw(IGBP_n_type) = &
      (/  4,  1,  5,  2,  3,  8,  9,  6,  6,  7, &
          8, 12,  7, 12, 13, 11,  0, 10, 10, 11 /)
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
   integer :: i, j, k, its, ite, jts, jte, kts, kte
   real (kind=kind_real) :: tem1, z0, z1, z2, w1, w2
   real (kind=kind_real), dimension(:,:,:), allocatable :: fzm_p, fzp_p

   real (kind=kind_real) :: kgkg_kgm2 !-- for var_clw, var_cli

   logical, SAVE :: l_detsfctyp = .false.

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
        ! This routine is trying to access GEOM (dimension, zgrid )
        ! TODO: Check: ite = nCells??  nCellsSolve???  MPI consideration ??? 
        its=1 ; ite = geom % nCells
        jts=1 ; jte = 1
        kts=1 ; kte = geom % nVertLevels
        allocate( fzm_p(its:ite,jts:jte,kts:kte), fzp_p(its:ite,jts:jte,kts:kte) )

        do j = jts,jte !dummy index to mimic a rectangular grid, NOT clean
        do k = kts+1,kte
        do i = its,ite
          tem1 = 1.0_kind_real/(geom % zgrid(k+1,i)-geom % zgrid(k-1,i))
          fzm_p(i,k,j) = (geom % zgrid(k,i)-geom % zgrid(k-1,i)) * tem1
          fzp_p(i,k,j) = (geom % zgrid(k+1,i)-geom % zgrid(k,i)) * tem1
          !pres2_p(i,k,j) = fzm_p(i,k,j)*pres_p(i,k,j) + fzp_p(i,k,j)*pres_p(i,k-1,j)
          field2d % array(k,i) = fzm_p(i,k,j)*r2d_ptr_a(k,i) + fzp_p(i,k,j)*r2d_ptr_a(k-1,i)
        enddo
        enddo
        enddo
        k = kte+1
        do j = jts,jte
        do i = its,ite
          z0 = geom % zgrid(k,i)
          z1 = 0.5_kind_real*(geom % zgrid(k,i)+geom % zgrid(k-1,i))
          z2 = 0.5_kind_real*(geom % zgrid(k-1,i)+geom % zgrid(k-2,i))
          w1 = (z0-z2)/(z1-z2)
          w2 = 1.0_kind_real-w1
          !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
          !pres2_p(i,k,j) = exp(w1*log(pres_p(i,k-1,j))+w2*log(pres_p(i,k-2,j)))
          field2d % array(k,i) = exp( w1*log(r2d_ptr_a(k-1,i)) + w2*log(r2d_ptr_a(k-1,i)) )
        enddo
        enddo
        k = kts
        do j = jts,jte
        do i = its,ite
          z0 = geom % zgrid(k,i)
          z1 = 0.5_kind_real*(geom % zgrid(k,i)+geom % zgrid(k+1,i))
          z2 = 0.5_kind_real*(geom % zgrid(k+1,i)+geom % zgrid(k+2,i))
          w1 = (z0-z2)/(z1-z2)
          w2 = 1.0_kind_real-w1
          !pres2_p(i,k,j) = w1*pres_p(i,k,j)+w2*pres_p(i,k+1,j)
          field2d % array(k,i) = w1*r2d_ptr_a(k,i) + w2*r2d_ptr_a(k+1,i)
        enddo
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


!     case ("Water_Fraction", "Land_Fraction", "Ice_Fraction", "Snow_Fraction", &
!           "Water_Temperature", "Land_Temperature", "Ice_Temperature", "Snow_Temperature", &
!           "Land_Type_Index", "Vegetation_Type", "Soil_Type" )
     case ("Water_Temperature", "Land_Temperature", "Ice_Temperature", "Snow_Temperature" )

        if( .not. l_detsfctyp) then
          call mpas_pool_get_field(pool_b, 'u10', field1d_src) ! as a dummy array

!----- These variables are moved into subroutine interp of mpas_fields_mod.F90
!        !--- fraction: TODO: split fractions... Need interp. weights
!                       ! see gsi/deter_sfc_mod.f90. sub deter_sfc (simple), sub deter_sfc_fov (complex)
!                       ! see wrfda/da_detsurtyp.inc
!          call mpas_pool_get_array(pool_b, "landmask", i1d_ptr_a) !"land-ocean mask (1=land ; 0=ocean)"
!          write(*,*) 'MIN/MAX of landmask=',minval(i1d_ptr_a),maxval(i1d_ptr_a)
!          call mpas_duplicate_field(field1d_src, field1d)
!          field1d % array(:) = 1.0_kind_real !real(i1d_ptr_a(:)) ! quantity and unit might change
!          field1d % fieldName = var_sfc_wfrac
!          call mpas_pool_add_field(pool_c, var_sfc_wfrac, field1d)
!
!          call mpas_duplicate_field(field1d_src, field1d)
!          field1d % array(:) = 0.0_kind_real !real(i1d_ptr_a(:)) ! quantity and unit might change
!          field1d % fieldName = var_sfc_lfrac
!          call mpas_pool_add_field(pool_c, var_sfc_lfrac, field1d)
!
!          call mpas_pool_get_array(pool_b, "xice", r1d_ptr_a) !"fractional area coverage of sea-ice"
!          write(*,*) 'MIN/MAX of xice=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
!          call mpas_duplicate_field(field1d_src, field1d)
!          field1d % array(:) = 0.0_kind_real !r1d_ptr_a(:) ! quantity and unit might change
!          field1d % fieldName = var_sfc_ifrac
!          call mpas_pool_add_field(pool_c, var_sfc_ifrac, field1d)
!
!          call mpas_pool_get_array(pool_b, "snowc", r1d_ptr_a) !"flag for snow on ground (=0 no snow; =1,otherwise"
!          write(*,*) 'MIN/MAX of snowc=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
!          call mpas_duplicate_field(field1d_src, field1d)
!          field1d % array(:) = 0.0_kind_real !r1d_ptr_a(:) ! quantity and unit might change
!          field1d % fieldName = var_sfc_sfrac
!          call mpas_pool_add_field(pool_c, var_sfc_sfrac, field1d)

        !--- temperature
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

!----- These variables are moved into subroutine interp of mpas_fields_mod.F90
!        !--- type: TODO: How to do "nearest neighbor" interp. for these integer variables??
!                      !: Choosing the grid value of maximum weights ??
!                      !: more complex one in future: Consider FOV ??
!          call mpas_pool_get_field(pool_b, 'landmask', field1di_src) ! as a dummy array
!
!          call mpas_pool_get_array(pool_b, "ivgtyp", i1d_ptr_a)
!          write(*,*) 'MIN/MAX of ivgtyp=',minval(i1d_ptr_a),maxval(i1d_ptr_a)
!          call mpas_duplicate_field(field1di_src, field1di)
!          field1di % array(:) = 1 !i1d_ptr_a(:)
!          field1di % fieldName = var_sfc_landtyp
!          call mpas_pool_add_field(pool_c, var_sfc_landtyp, field1di)
!
!          call mpas_duplicate_field(field1di_src, field1di)
!          field1di % array(:) = 1 !max(1,usgs_to_crtm_mw(i1d_ptr_a(:))  !chage category, TMP: hardcoded as 1
!          field1di % fieldName = var_sfc_vegtyp
!          call mpas_pool_add_field(pool_c, var_sfc_vegtyp, field1di)
!
!          call mpas_pool_get_array(pool_b, "isltyp", i1d_ptr_a)
!          write(*,*) 'MIN/MAX of isltyp=',minval(i1d_ptr_a),maxval(i1d_ptr_a)
!          call mpas_duplicate_field(field1di_src, field1di)
!          field1di % array(:) = 1 !max(1,wrf_to_crtm_soil(i1d_ptr_a(:)) ! chage category, TMP: hardcoded as 1
!          field1di % fieldName = var_sfc_soiltyp
!          call mpas_pool_add_field(pool_c, var_sfc_soiltyp, field1di)

          l_detsfctyp = .true.
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

!----- These two variables are moved into subroutine interp of mpas_fields_mod.F90
!     case ("Sfc_Wind_Speed")
!        call mpas_pool_get_array(pool_b, "u10", r1d_ptr_a)
!        call mpas_pool_get_array(pool_b, "v10", r1d_ptr_b)
!        write(*,*) 'MIN/MAX of u10=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
!        write(*,*) 'MIN/MAX of v10=',minval(r1d_ptr_b),maxval(r1d_ptr_b)
!        call mpas_pool_get_field(pool_b, 'skintemp', field1d_src) ! as a dummy array
!        call mpas_duplicate_field(field1d_src, field1d)
!        field1d % array(:) = sqrt( r1d_ptr_a(:)**2 + r1d_ptr_b(:)**2 ) ! ws = sqrt(u**2+v**2) [m/s]
!        field1d % fieldName = var_sfc_wspeed
!        call mpas_pool_add_field(pool_c, var_sfc_wspeed, field1d)
!        write(*,*) "end-of ",var_sfc_wspeed
!
!     case ("Sfc_Wind_Direction") !-var_sfc_wdir 
!        call mpas_pool_get_array(pool_b, "u10", r1d_ptr_a)
!        call mpas_pool_get_array(pool_b, "v10", r1d_ptr_b)
!        write(*,*) 'MIN/MAX of u10=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
!        write(*,*) 'MIN/MAX of v10=',minval(r1d_ptr_b),maxval(r1d_ptr_b)
!        call mpas_pool_get_field(pool_b, 'skintemp', field1d_src) ! as a dummy array
!        call mpas_duplicate_field(field1d_src, field1d)
!        field1d % array(:) = 0.0_kind_real !NOTE: Do this here? or after interpolated to obs location ??
!        field1d % fieldName = var_sfc_wdir
!        call mpas_pool_add_field(pool_c, var_sfc_wdir, field1d)
!        write(*,*) "end-of ",var_sfc_wdir
!
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
   subroutine convert_mpas_field2ufoTL(traj_pool_a, traj_pool_b, pool_a, pool_b, pool_c, fieldname, nfield)

   implicit none

   type (mpas_pool_type), pointer, intent(in) :: traj_pool_a, traj_pool_b ! subFields, auxFields
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

        call mpas_pool_get_array(traj_pool_a,    'theta', traj_r2d_a)
        call mpas_pool_get_array(traj_pool_a, 'index_qv', traj_r2d_b)
        call mpas_pool_get_array(traj_pool_b, 'pressure', traj_r2d_c)
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
   subroutine convert_mpas_field2ufoAD(traj_pool_a, traj_pool_b, pool_a, pool_b, pool_c, fieldname, nfield)

   implicit none

   type (mpas_pool_type), pointer, intent(inout) :: traj_pool_a, traj_pool_b ! subFields, auxFields
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

        call mpas_pool_get_array(traj_pool_a,    'theta', traj_r2d_a)
        call mpas_pool_get_array(traj_pool_a, 'index_qv', traj_r2d_b)
        call mpas_pool_get_array(traj_pool_b, 'pressure', traj_r2d_c)
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
