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

   contains


   subroutine update_mpas_field(domain, pool_a)

   implicit none

   type (domain_type), pointer, intent(in) :: domain
   type (mpas_pool_type), pointer, intent(inout) :: pool_a

   type (mpas_pool_type), pointer :: state, diag, mesh, pool_b
   type (mpas_pool_iterator_type) :: poolItr
   integer, parameter :: maxfield = 1
   character (len=20) :: fieldname(1:maxfield)
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
   character (len=20), intent(in) :: fieldname(1:nfield) ! ufo

   type (mpas_pool_iterator_type) :: poolItr
   type (mpas_pool_type), pointer :: state, diag, mesh, allFields
   real (kind=RKIND), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   real (kind=RKIND), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   integer :: ii

   !--- create new pull for ufo_vars
   call mpas_pool_create_pool(pool_b, nfield)

   !--- some ufo_vars can be found outside-of-subFields.
   allFields => domain % blocklist % allFields


   call mpas_pool_begin_iteration(pool_a)

   do while ( mpas_pool_get_next_member(pool_a, poolItr) )

     ! Pools may in general contain dimensions, namelist options, fields, or other pools,
     ! so we select only those members of the pool that are fields
     if (poolItr % memberType == MPAS_POOL_FIELD) then
       ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
       if (poolItr % dataType == MPAS_POOL_REAL) then
          write(*,*)'Looking for ',trim(poolItr % memberName)
          do ii=1, nfield
             if ( trim(fieldname(ii)).eq.trim(poolItr % memberName)) then
               if (poolItr % nDims == 0) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r0d_ptr_a)
                  !call mpas_pool_get_array(diag, trim(poolItr % memberName), r0d_ptr_b)
                  !r0d_ptr_a = r0d_ptr_a + r0d_ptr_b
               else if (poolItr % nDims == 1) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
                  !call mpas_pool_get_array(diag, trim(poolItr % memberName), r1d_ptr_b)
                  r1d_ptr_a = r1d_ptr_b
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_a)
                  !call mpas_pool_get_array(diag, trim(poolItr % memberName), r2d_ptr_b)
                  !r2d_ptr_a = r2d_ptr_b
                  write(0,*)'Update MIN/MAX: ',trim(fieldname(ii)),minval(r2d_ptr_a),maxval(r2d_ptr_a)
               end if
             end if
          end do
       end if
     end if

   end do
   
  
   end subroutine convert_mpas_field2ufo


   end module mpas2ufo_vars_mod
