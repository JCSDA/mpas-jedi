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
   !---- other "var_struct", such as "diag", "diag_physics", or "allFields" can be found inside domain
   subroutine convert_mpas_field2ufo(domain, pool_a, pool_b, fieldname, nfield)

   implicit none

   type (domain_type), pointer, intent(in) :: domain
   type (mpas_pool_type), pointer, intent(in) :: pool_a
   type (mpas_pool_type), pointer, intent(out) :: pool_b
   integer, intent(in) :: nfield
   character (len=*), intent(in) :: fieldname(:) ! ufo

   type (mpas_pool_iterator_type) :: poolItr
   type (mpas_pool_type), pointer :: allFields, clone_pool_a
   real (kind=RKIND), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b, r2d_ptr_c
   real (kind=RKIND), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b

   type (field2DReal), pointer :: field2d, field2d_a, field2d_b, field2d_c
   integer :: ii, ivar
   !--- test
   call mpas_pool_create_pool(clone_pool_a)
   call mpas_pool_clone_pool(pool_a,clone_pool_a)

   !--- create new pull for ufo_vars
   call mpas_pool_create_pool(pool_b, nfield)

   !--- some ufo_vars can be found outside-of-subFields.
   allFields => domain % blocklist % allFields


   do ivar=1, nfield
   write(*,*) '---------inside do/select case, ivar, trim(fieldname(ivar))=',ivar,trim(fieldname(ivar))

     select case (trim(fieldname(ivar)))

     case ( "virtual_temperature" ) !-var_tv 
        !get theta from subField
        !get pressure from subField
        !get qv from subfield
        !do the conversion : need 'nl', 'tl', 'ad' ?!
   call mpas_pool_begin_iteration(pool_a)
   do while ( mpas_pool_get_next_member(pool_a, poolItr) )
     if (poolItr % memberType == MPAS_POOL_FIELD) then
       if (poolItr % dataType == MPAS_POOL_REAL) then
          write(*,*)'Looking for ',trim(poolItr % memberName)
               if (poolItr % nDims == 0) then
                  !call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r0d_ptr_a)
               else if (poolItr % nDims == 1) then
                  !call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
               else if (poolItr % nDims == 2) then
                  if( trim(poolItr % memberName) .eq. 'theta' ) &
                       call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_a)
                  if( trim(poolItr % memberName) .eq. 'pressure' ) &
                       call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_b)
                  if( trim(poolItr % memberName) .eq. 'index_qv' ) &
                       call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_c)
               end if
       end if
     end if
   end do

        call mpas_pool_get_array(domain % blocklist % allFields, 'pressure', r2d_ptr_b)

        call mpas_pool_get_field(clone_pool_a, 'theta', field2d) ! as a dummy array
        write(*,*) 'MIN/MAX of theta=',minval(r2d_ptr_a),maxval(r2d_ptr_a)
        write(*,*) 'MIN/MAX of pressure=',minval(r2d_ptr_b),maxval(r2d_ptr_b)
        write(*,*) 'MIN/MAX of qv=',minval(r2d_ptr_c),maxval(r2d_ptr_c)
        field2d % array(:,:) = r2d_ptr_a(:,:) / ( (100000.0/r2d_ptr_b(:,:))**(287.05/1005.7) ) * &
                               (1.0 + (461.50/287.05 - 1.0)*r2d_ptr_c(:,:))
             !Tv = T * ( 1.0 + (rv/rd – 1)*qv), rv=461.50 , rd=287.05
                  !T = theta / ( (p0/pressure)**rd_over_cp) : to sensible temperature	
        write(*,*) 'MIN/MAX of Tv=',minval(field2d % array(:,:)),maxval(field2d % array(:,:))
        field2d % fieldName = var_tv
        call mpas_pool_add_field(pool_b, var_tv, field2d)
        write(*,*) "end-of ",var_tv

     case ("atmosphere_ln_pressure_coordinate") !-var_prsl
        !get pressure from subField
        !do the conversion : need 'nl', 'tl', 'ad' ?!
   call mpas_pool_begin_iteration(pool_a)
   do while ( mpas_pool_get_next_member(pool_a, poolItr) )
     if (poolItr % memberType == MPAS_POOL_FIELD) then
       if (poolItr % dataType == MPAS_POOL_REAL) then
               if (poolItr % nDims == 0) then
                  !call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r0d_ptr_a)
               else if (poolItr % nDims == 1) then
                  !call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
               else if (poolItr % nDims == 2) then
                  if( trim(poolItr % memberName) .eq. 'pressure' ) &
                       call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_a)
               end if
       end if
     end if
   end do

        write(*,*) 'step 1: test: get from domain % blocklist % allFields'
        call mpas_pool_get_array(domain % blocklist % allFields, 'pressure', r2d_ptr_a)

        call mpas_pool_get_field(clone_pool_a, 'theta', field2d) ! as a dummy array
        field2d % array(:,:) = log( r2d_ptr_a(:,:)/100./10. ) !- Pa -> hPa -> kPa
write(*,*) 'ln_p, ln_p, ln_p =', field2d % array(:,1)
        write(*,*) 'MIN/MAX of ln_p=',minval(field2d % array(:,:)),maxval(field2d % array(:,:))
        field2d % fieldName = var_prsl
        call mpas_pool_add_field(pool_b, var_prsl, field2d)
        write(*,*) "end-of ",var_prsl

     case ("humidity_mixing_ratio") !-var_mixr
        call mpas_pool_get_array(pool_a, "index_qv", r2d_ptr_a)
        call mpas_pool_get_field(clone_pool_a, 'pressure', field2d) ! as a dummy array
        field2d % array(:,:) = r2d_ptr_a(:,:)
        field2d % fieldName = var_mixr
        call mpas_pool_add_field(pool_b, var_mixr, field2d)
        write(*,*) "end-of ",var_mixr

     case ("air_pressure") !-var_prs
        call mpas_pool_get_array(pool_a, "pressure", r2d_ptr_a)
        call mpas_pool_get_field(clone_pool_a, 'pressure', field2d) ! as a dummy array
        field2d % array(:,:) = r2d_ptr_a(:,:)
        field2d % fieldName = var_prs
        call mpas_pool_add_field(pool_b, var_prs, field2d)
        write(*,*) "end-of ",var_prs

     case ("air_pressure_levels")
     case ("mass_concentration_of_ozone_in_air")
     case ("mass_concentration_of_carbon_dioxide_in_air")
     case ("atmosphere_mass_content_of_cloud_liquid_water")
     case ("atmosphere_mass_content_of_cloud_ice")
     case ("effective_radius_of_cloud_liquid_water_particle")
     case ("effective_radius_of_cloud_ice_particle")

     case ("Water_Fraction")
        call mpas_pool_get_array(allFields, "", r1d_ptr_a)
     case ("Land_Fraction")
        call mpas_pool_get_array(allFields, "", r1d_ptr_a)
     case ("Ice_Fraction")
        call mpas_pool_get_array(allFields, "", r1d_ptr_a)
     case ("Snow_Fraction")
        call mpas_pool_get_array(allFields, "", r1d_ptr_a)
     case ("Water_Temperature")
        call mpas_pool_get_array(allFields, "skintemp", r1d_ptr_a)
     case ("Land_Temperature")
        call mpas_pool_get_array(allFields, "skintemp", r1d_ptr_a)
     case ("Ice_Temperature")
        call mpas_pool_get_array(allFields, "skintemp", r1d_ptr_a)
     case ("Snow_Temperature")
        call mpas_pool_get_array(allFields, "skintemp", r1d_ptr_a)

     case ("Snow_Depth")
        call mpas_pool_get_array(allFields, "snowh", r1d_ptr_a)
     case ("Vegetation_Fraction")
        call mpas_pool_get_array(allFields, "vegfra", r1d_ptr_a)
     case ("Sfc_Wind_Speed", "Sfc_Wind_Direction")
        call mpas_pool_get_array(pool_a, "u10", r1d_ptr_a)
        call mpas_pool_get_array(pool_a, "v10", r1d_ptr_b)
     case ("Lai")
        call mpas_pool_get_array(allFields, "lai", r1d_ptr_a)
     case ("Soil_Moisture")
        call mpas_pool_get_array(allFields, "smois", r1d_ptr_a)
     case ("Soil_Temperature")
        call mpas_pool_get_array(allFields, "tslb", r1d_ptr_a)
     case ("Land_Type_Index")
        call mpas_pool_get_array(allFields, "", r1d_ptr_a)
     case ("Vegetation_Type")
        call mpas_pool_get_array(allFields, "ivgtyp", r1d_ptr_a)
     case ("Soil_Type")
        call mpas_pool_get_array(allFields, "isltyp", r1d_ptr_a)
     case ("ice_concentration")
        call mpas_pool_get_array(allFields, "seaice", r1d_ptr_a)



     end select

   end do !ivar


   call mpas_pool_empty_pool(clone_pool_a)
   call mpas_pool_destroy_pool(clone_pool_a)
  
   end subroutine convert_mpas_field2ufo

!-----------------------------------------------------------------------------------

#ifdef TL_READY
   !---- pool_a : self % subFields
   !---- other "var_struct", such as "diag", "diag_physics", or "allFields" can be found inside domain
   subroutine convert_mpas_field2ufoTL(domain, pool_a, pool_b, fieldname, nfield)

   implicit none

   type (domain_type), pointer, intent(in) :: domain
   type (mpas_pool_type), pointer, intent(in) :: pool_a
   type (mpas_pool_type), pointer, intent(out) :: pool_b
   integer, intent(in) :: nfield
   character (len=*), intent(in) :: fieldname(:) ! ufo

   type (mpas_pool_iterator_type) :: poolItr
   type (mpas_pool_type), pointer :: allFields, clone_pool_a
   real (kind=RKIND), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b, r2d_ptr_c
   real (kind=RKIND), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b

   type (field2DReal), pointer :: field2d, field2d_a, field2d_b, field2d_c
   integer :: ii, ivar
  
   end subroutine convert_mpas_field2ufoTL
#endif

!-----------------------------------------------------------------------------------

   end module mpas2ufo_vars_mod
