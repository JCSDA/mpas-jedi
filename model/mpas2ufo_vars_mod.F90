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

   use ufo_vars_mod !, only : var_tv, var_prsl

   contains


   subroutine update_mpas_field(domain, pool_a)

   implicit none

   type (domain_type), pointer, intent(in) :: domain
   type (mpas_pool_type), pointer, intent(inout) :: pool_a

   type (mpas_pool_type), pointer :: state, diag, mesh, pool_b
   type (mpas_pool_iterator_type) :: poolItr
   integer, parameter :: maxfield = 1
   character (len=22) :: fieldname(1:maxfield)
   real (kind=RKIND), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   real (kind=RKIND), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
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
   subroutine convert_mpas_field2ufo(pool_a, pool_b, pool_c, fieldname, nfield)

   implicit none

   type (mpas_pool_type), pointer, intent(in) :: pool_a, pool_b ! subFields, auxFields
   type (mpas_pool_type), pointer, intent(out) :: pool_c
   integer, intent(in) :: nfield
   character (len=*), intent(in) :: fieldname(:) ! ufo

   type (mpas_pool_iterator_type) :: poolItr
   type (mpas_pool_type), pointer :: allFields
   real (kind=RKIND), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b, r2d_ptr_c
   real (kind=RKIND), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b

   type (field2DReal), pointer :: field2d, field2d_a, field2d_b, field2d_c, field2d_src
   integer :: ii, ivar

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

#define READY_TRAJ
#ifdef READY_TRAJ
!%%%
!%%%  NL 
!%%%  TODO: use clean formula
        field2d % array(:,:) = r2d_ptr_a(:,:) / ( (100000.0/r2d_ptr_c(:,:))**(287.05/1005.7) ) * &
                               (1.0 + (461.50/287.05 - 1.0)*r2d_ptr_b(:,:))
             !Tv = T * ( 1.0 + (rv/rd – 1)*qv), rv=461.50 , rd=287.05
                  !T = theta / ( (p0/pressure)**rd_over_cp) : to sensible temperature	
#else
        field2d % array(:,:) = r2d_ptr_a(:,:) !BJJ test: pass theta as T_v
#endif
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

        field2d % array(:,:) = log( r2d_ptr_a(:,:)/100./10. ) !- Pa -> hPa ->cb
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
     case ("ice_concentration")
        call mpas_pool_get_array(pool_b, "seaice", r1d_ptr_a)

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
   real (kind=RKIND), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b, r2d_ptr_c
   real (kind=RKIND), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   real (kind=RKIND), dimension(:,:), pointer :: traj_r2d_a, traj_r2d_b, traj_r2d_c !BJJ test

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

#ifdef READY_TRAJ
!%%%
!%%%  NL currently. --> TL with trajectory!
!%%%  TODO: use clean formula
        !field2d % array(:,:) = r2d_ptr_a(:,:) / ( (100000.0/r2d_ptr_c(:,:))**(287.05/1005.7) ) * &
        !                       (1.0 + (461.50/287.05 - 1.0)*r2d_ptr_b(:,:))
        !     !Tv = T * ( 1.0 + (rv/rd – 1)*qv), rv=461.50 , rd=287.05
        !          !T = theta / ( (p0/pressure)**rd_over_cp) : to sensible temperature	
        field2d % array(:,:) = ( ( 1.0 + (461.50/287.05 - 1.0)*traj_r2d_b(:,:) ) * r2d_ptr_a(:,:) &
                                + traj_r2d_a(:,:) * (461.50/287.05 - 1.0) * r2d_ptr_b(:,:) ) &
                               / ( (100000.0/traj_r2d_c(:,:))**(287.05/1005.7) )
        write(*,*) 'MIN/MAX of theta=',minval(r2d_ptr_a(:,:)),maxval(r2d_ptr_a(:,:))
#else
        field2d % array(:,:) = r2d_ptr_a(:,:) !BJJ test: pass theta as T_v
#endif
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
     case ("ice_concentration")
        call mpas_pool_get_array(pool_b, "seaice", r1d_ptr_a)

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
   real (kind=RKIND), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b, r2d_ptr_c
   real (kind=RKIND), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   real (kind=RKIND), dimension(:,:), pointer :: traj_r2d_a, traj_r2d_b, traj_r2d_c !BJJ test

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

#ifdef READY_TRAJ
!%%%
!%%%  NL currently. --> AD with trajectory!
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
        r2d_ptr_a(:,:)=0.0
        r2d_ptr_b(:,:)=0.0
        r2d_ptr_a(:,:) = r2d_ptr_a(:,:) + ( 1.0 + (461.50/287.05 - 1.0)*traj_r2d_b(:,:) ) / &
                         ( (100000.0/traj_r2d_c(:,:))**(287.05/1005.7) ) * field2d % array(:,:)
        r2d_ptr_b(:,:) = r2d_ptr_b(:,:) + traj_r2d_a(:,:) * (461.50/287.05 - 1.0) / &
                         ( (100000.0/traj_r2d_c(:,:))**(287.05/1005.7) ) * field2d % array(:,:)
        write(*,*) 'MIN/MAX of theta=',minval(r2d_ptr_a(:,:)),maxval(r2d_ptr_a(:,:))
        write(*,*) 'MIN/MAX of index_qv=',minval(r2d_ptr_b(:,:)),maxval(r2d_ptr_b(:,:))
#else
        call mpas_pool_get_field(pool_c, 'virtual_temperature', field2d)
        write(*,*) 'MIN/MAX of Tv=',minval(field2d % array(:,:)),maxval(field2d % array(:,:))
        r2d_ptr_a(:,:)=0.0
        r2d_ptr_a(:,:) = r2d_ptr_a(:,:) + field2d % array(:,:) !BJJ test: pass theta as T_v
        write(*,*) 'MIN/MAX of theta=',minval(r2d_ptr_a(:,:)),maxval(r2d_ptr_a(:,:))
#endif

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
     case ("ice_concentration")
        call mpas_pool_get_array(pool_b, "seaice", r1d_ptr_a)

     end select

   end do !ivar


   !call mpas_pool_empty_pool(clone_pool_a)
   !call mpas_pool_destroy_pool(clone_pool_a)
  
   end subroutine convert_mpas_field2ufoAD

!-----------------------------------------------------------------------------------

   end module mpas2ufo_vars_mod
