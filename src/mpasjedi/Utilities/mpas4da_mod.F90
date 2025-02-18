! Copyright (c) 2018, National Atmospheric for Atmospheric Research (NCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html

module mpas4da_mod

   !***********************************************************************
   !
   !  Module mpas4da_mod to encapsulate operations needed for
   !  Data assimilation purpose.
   !  It can be used from /somewhere/MPAS/src/operators
   !  or from /somewhere/mpas-bundle/mpas/model (OOPS)
   !> \author  Gael Descombes/Mickael Duda NCAR/MMM
   !> \date    January 2018
   !
   !-----------------------------------------------------------------------

!fckit
use fckit_log_module, only: fckit_log

!oops
use kinds, only: kind_real
use random_mod, only: normal_distribution

!ufo
use ufo_vars_mod

!MPAS-Model
use mpas_abort, only: mpas_dmpar_global_abort
use mpas_constants
use mpas_derived_types
use mpas_dmpar
use mpas_field_routines
use mpas_kind_types, only: ShortStrKIND
use mpas_pool_routines

!mpas-jedi
use mpas_constants_mod
use mpas_geom_mod, only: mpas_geom, pool_has_field, getSolveDimSizes

private

public :: &
   da_operator_addition, &
   da_copy_all2sub_fields, &
   da_copy_sub2all_fields, &
   da_template_pool, &
   !mpas_pool_template_field, &
   da_random, &
   da_operator, &
   da_self_mult, &
   da_constant, &
   da_posdef, &
   da_setval, &
   da_axpy, &
   da_gpnorm, &
   da_fldrms, &
   da_dot_product, &
   cvt_oopsmpas_date, &
   uv_cell_to_edges, &
   r3_normalize

character(len=1024) :: message

contains

   !***********************************************************************
   !
   !  function field_is_scalar
   !
   !> \brief   Test for a form of water vapor or hydrometeor of interest
   !> \details
   !>  At various places in this module we wish to test for the case where
   !>  a string is one of several 'q?' values. Rather than repeat that
   !>  logic many times, it is encapsulated here.
   !
   !-----------------------------------------------------------------------
   pure function field_is_scalar(fieldName)

      implicit none

      character (len=*), intent(in) :: fieldName
      logical :: field_is_scalar
      field_is_scalar = any(trim(fieldName) == &
                    (/'qv', 'qc', 'qi', 'qr', 'qs', 'qg', 'qh', 'nc', 'ni', 'nr', 'ns', 'ng', 'nh'/))

   end function

   !***********************************************************************
   !
   !  function match_scalar
   !
   !> \brief   Test for a form of water vapor or hydrometeor of interest
   !> \author  Steven Vahl
   !> \date    11 July 2019
   !> \details
   !>  Test for the case where one string is 'scalars' and another string
   !>   matches with available scalar variables.
   !
   !-----------------------------------------------------------------------
   pure function match_scalar(scalarName, fieldName)

      implicit none

      character (len=*), intent(in) :: scalarName
      character (len=*), intent(in) :: fieldName
      logical :: match_scalar

      match_scalar = ( &
         scalarName == 'scalars' .and. &
         field_is_scalar(fieldName) )

   end function

   !***********************************************************************
   !
   !  subroutine mpas_pool_demo
   !
   !> \brief   Demonstrate basic usage of MPAS pools
   !> \author  Michael Duda
   !> \date    20 December 2017
   !> \details
   !>  This routine provides a simple demonstration of how to construct a new
   !>  pool at runtime, add members (fields) to the pool, and to perform generic
   !>  operations on that pool.
   !
   !-----------------------------------------------------------------------
   subroutine mpas_pool_demo(block)

      implicit none

      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: structs
      type (mpas_pool_type), pointer :: allFields
      type (mpas_pool_type), pointer :: da_state
      type (mpas_pool_type), pointer :: da_state_incr

      type (field2DReal), pointer :: field

      write(0,*) '****** Begin pool demo routine ******'

      structs => block % structs
      allFields => block % allFields

      !
      ! Create a new pool
      !
      call mpas_pool_create_pool(da_state)

      !
      ! Get pointers to several fields from the allFields pool, and add
      ! those fields to the da_state pool as well
      !
      call mpas_pool_get_field(allFields, 'theta', field)
      call mpas_pool_add_field(da_state, 'theta', field)
      write(0,*) 'Now, max value of theta is ', maxval(field % array),minval(field % array)
      field % array(:,:) = MPAS_JEDI_ONE_kr
      write(0,*)'Dimensions Field: ',field % dimSizes(:)

      call mpas_pool_get_field(allFields, 'rho', field)
      call mpas_pool_add_field(da_state, 'rho', field)
      write(0,*) 'Now, max value of rho is ', maxval(field % array),minval(field % array)
      field % array(:,:) = MPAS_JEDI_ONE_kr

      !
      ! Create another pool
      !
      call mpas_pool_create_pool(da_state_incr)

      !
      ! Duplicate the members of da_state into da_state_incr, and do a deep
      ! copy of the fields from da_state to da_state_incr
      !
      call mpas_pool_clone_pool(da_state, da_state_incr)

      !
      ! Call example algebra routine to compute A = A + B for all fields in
      ! the da_state and da_state_inc pools
      !
      call da_operator_addition(da_state, da_state_incr)

      call mpas_pool_get_field(da_state_incr, 'rho', field)
      write(0,*) 'Now, max value of rho_incr is ', maxval(field % array)

      call mpas_pool_get_field(da_state, 'rho', field)
      write(0,*) 'Now, max value of rho is ', maxval(field % array)

      !
      ! Before destroying a pool, we should remove any fields that are
      ! still referenced by other active pools to avoid deallocating them
      !
      call mpas_pool_empty_pool(da_state)

      !
      ! Destroy the now-empty da_state pool
      !
      call mpas_pool_destroy_pool(da_state)

      !
      ! Destroy the da_state_incr pool, deallocating all of its
      ! fields in the process (because this pool was not emptied)
      !
      call mpas_pool_destroy_pool(da_state_incr)

      write(0,*) '****** End pool demo routine ******'

   end subroutine mpas_pool_demo

   !***********************************************************************
   !
   !  subroutine da_operator_addition
   !
   !> \brief   Performs A = A + B for pools A and B
   !> \author  Michael Duda
   !> \date    20 December 2017
   !> \details
   !>  Given two pools, A and B, where the fields in B are a subset of
   !>  the fields in A, this routine adds the fields in B to fields in A
   !>  with the same name. When A and B contain identical fields, this
   !>  is equivalent to A = A + B.
   !
   !-----------------------------------------------------------------------
   subroutine da_operator_addition(pool_a, pool_b)

      implicit none

      type (mpas_pool_type), pointer :: pool_a, pool_b

      type (mpas_pool_iterator_type) :: poolItr
      real (kind=RKIND), pointer :: r0d_ptr_a, r0d_ptr_b
      real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
      real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
      real (kind=RKIND), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b

      !
      ! Iterate over all fields in pool_b, adding them to fields of the same
      ! name in pool_a
      !
      call mpas_pool_begin_iteration(pool_b)

      do while ( mpas_pool_get_next_member(pool_b, poolItr) )

         ! Pools may in general contain dimensions, namelist options, fields, or other pools,
         ! so we select only those members of the pool that are fields
         if (poolItr % memberType == MPAS_POOL_FIELD) then

            ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
            if (poolItr % dataType == MPAS_POOL_REAL) then

               ! Depending on the dimensionality of the field, we need to set pointers of
               ! the correct type
               if (poolItr % nDims == 0) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r0d_ptr_a)
                  call mpas_pool_get_array(pool_b, trim(poolItr % memberName), r0d_ptr_b)
                  r0d_ptr_a = r0d_ptr_a + r0d_ptr_b
               else if (poolItr % nDims == 1) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
                  call mpas_pool_get_array(pool_b, trim(poolItr % memberName), r1d_ptr_b)
                  r1d_ptr_a = r1d_ptr_a + r1d_ptr_b
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_a)
                  call mpas_pool_get_array(pool_b, trim(poolItr % memberName), r2d_ptr_b)
                  r2d_ptr_a = r2d_ptr_a + r2d_ptr_b
                  write(message,*) 'Operator add MIN/MAX: ',minval(r2d_ptr_a),maxval(r2d_ptr_a)
                  call fckit_log%debug(message)
               else if (poolItr % nDims == 3) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r3d_ptr_a)
                  call mpas_pool_get_array(pool_b, trim(poolItr % memberName), r3d_ptr_b)
                  r3d_ptr_a = r3d_ptr_a + r3d_ptr_b
               end if

            end if
         end if
      end do

   end subroutine da_operator_addition


   !***********************************************************************
   !
   !  subroutine da_copy_all2sub_fields
   !
   !> \brief   Performs a copy of allfield to a sub pool A
   !> \author  Gael Desccombes
   !> \date    5 February 2018
   !> \details
   !>  Given two pools, allfields and A, where the fields in A are a subset of
   !>  the fields in allfields, this routine copy the fields allfields to fields in A
   !>  with the same name.
   !
   !-----------------------------------------------------------------------
   subroutine da_copy_all2sub_fields(geom, pool_a)

      implicit none

      type (mpas_geom), pointer, intent(in) :: geom
      type (mpas_pool_type), pointer, intent(inout) :: pool_a
      type (mpas_pool_type), pointer :: pool_b, state

      type (mpas_pool_iterator_type) :: poolItr_a, poolItr_b
      real (kind=RKIND), pointer :: r0d_ptr_a, r0d_ptr_b
      real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
      real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
      integer, pointer :: index_scalar
      character(len=MAXVARLEN) :: targetName, ioName
      real(kind=RKIND) :: scaling_factor

      type (field2DReal), pointer :: field2d
      type (field3DReal), pointer :: field3d

      pool_b => geom % domain % blocklist % allFields
      call mpas_pool_get_subpool(geom % domain % blocklist % structs,'state',state)
      !
      ! Iterate over all fields in pool_b, adding them to fields of the same
      ! name in pool_a
      !
      call mpas_pool_begin_iteration(pool_b)

      do while ( mpas_pool_get_next_member(pool_b, poolItr_b) )

         ! Pools may in general contain dimensions, namelist options, fields, or other pools,
         ! so we select only those members of the pool that are fields
         if (poolItr_b % memberType == MPAS_POOL_FIELD) then

            ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
            if (poolItr_b % dataType == MPAS_POOL_REAL) then

             call mpas_pool_begin_iteration(pool_a)
             do while ( mpas_pool_get_next_member(pool_a, poolItr_a) )

               !1. start from pool field name
               targetName = trim(poolItr_a % memberName)
               !2. If exists io_name, update with that
               if ( geom%has_io_name(targetName) ) then
                   ioName = geom%io_name(targetName)
               else
                   ioName = 'none'
               end if
               !3. If exists io_scaling_factor, 
               if ( geom%has_io_scaling_factor(targetName) ) then
                   scaling_factor = real(geom%io_scaling_factor(targetName),RKIND)
               else
                   scaling_factor = 0.0_RKIND
               end if

               if ( ( trim(poolItr_b % memberName).eq.trim(targetName) ) .or. &
                    ( trim(poolItr_b % memberName).eq.trim(ioName)     ) ) then
                  ! Depending on the dimensionality of the field, we need to set pointers of
                  ! the correct type
                  if (poolItr_b % nDims == 0) then
                     call mpas_pool_get_array(pool_a, trim(poolItr_a % memberName), r0d_ptr_a)
                     call mpas_pool_get_array(pool_b, trim(poolItr_b % memberName), r0d_ptr_b)
                     r0d_ptr_a = r0d_ptr_b
                  else if (poolItr_b % nDims == 1) then
                     call mpas_pool_get_array(pool_a, trim(poolItr_a % memberName), r1d_ptr_a)
                     call mpas_pool_get_array(pool_b, trim(poolItr_b % memberName), r1d_ptr_b)
                     if ( scaling_factor /= 0.0_RKIND ) then
                        r1d_ptr_a = r1d_ptr_b * scaling_factor
                     else
                        r1d_ptr_a = r1d_ptr_b
                     end if
                  else if (poolItr_b % nDims == 2) then
                     write(message,*) 'poolItr_b % memberName=',trim(poolItr_b % memberName)
                     call fckit_log%debug(message)
                     call mpas_pool_get_array(pool_a, trim(poolItr_a % memberName), r2d_ptr_a)
                     call mpas_pool_get_array(pool_b, trim(poolItr_b % memberName), r2d_ptr_b)
                     if ( scaling_factor /= 0.0_RKIND ) then
                        r2d_ptr_a = r2d_ptr_b * scaling_factor
                     else
                        r2d_ptr_a = r2d_ptr_b
                     end if
                     write(message,*) 'Copy all2sub field MIN/MAX: ',trim(poolItr_b % memberName), &
                                      minval(r2d_ptr_a),maxval(r2d_ptr_a)
                     call fckit_log%debug(message)
                  end if

               else if ( match_scalar(trim(poolItr_b % memberName), trim(ioName)) ) then ! Here we know qx or nx variables are all defined as ioName
                  write(message,*) 'Copy all2sub field: Looking at SCALARS now',trim(poolItr_a % memberName)
                  call fckit_log%debug(message)
                  call mpas_pool_get_dimension(state, 'index_'//trim(ioName), index_scalar)
                  if (index_scalar .gt. 0) then
                     call mpas_pool_get_field(pool_a, trim(poolItr_a % memberName), field2d)
                     call mpas_pool_get_field(pool_b, trim(poolItr_b % memberName), field3d)
                     if ( scaling_factor /= 0.0_RKIND ) then
                        field2d % array(:,:) = field3d % array(index_scalar,:,:) * scaling_factor
                     else
                        field2d % array(:,:) = field3d % array(index_scalar,:,:)
                     end if
                     write(message,*) 'Copy all2sub field MIN/MAX: ',trim(poolItr_a % memberName), &
                                      minval(field2d % array), maxval(field2d % array)
                     call fckit_log%debug(message)
                  else
                     write(message,*) 'WARNING in Copy all2sub field; ',trim(poolItr_a % memberName), &
                                      'not available from MPAS'
                     call fckit_log%debug(message)
                  end if
               end if
            end do
           end if
         end if
      end do

   end subroutine da_copy_all2sub_fields


   !***********************************************************************
   !
   !  subroutine da_copy_sub2all_fields
   !
   !> \brief   Performs a copy of a sub pool A to allfields
   !> \author  Gael Desccombes
   !> \date    5 February 2018
   !> \details
   !>  Given two pools, allfields and A, where the fields in A are a subset of
   !>  the fields in allfields, this routine copy the subfields to allfields
   !>  with the same name.
   !
   !-----------------------------------------------------------------------
   subroutine da_copy_sub2all_fields(geom, pool_a)

      implicit none

      type (mpas_geom), pointer, intent(in) :: geom
      type (mpas_pool_type), pointer, intent(in) :: pool_a
      type (mpas_pool_type), pointer :: pool_b, state

      type (mpas_pool_iterator_type) :: poolItr_a, poolItr_b
      real (kind=RKIND), pointer :: r0d_ptr_a, r0d_ptr_b
      real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
      real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
      integer, pointer :: index_scalar
      character(len=MAXVARLEN) :: targetName, ioName
      real(kind=RKIND) :: scaling_factor

      type (field2DReal), pointer :: field2d
      type (field3DReal), pointer :: field3d

      pool_b => geom % domain % blocklist % allFields
      call mpas_pool_get_subpool(geom % domain % blocklist % structs,'state',state)
      !
      ! Iterate over all fields in pool_b, adding them to fields of the same
      ! name in pool_a
      !
      call mpas_pool_begin_iteration(pool_b)

      do while ( mpas_pool_get_next_member(pool_b, poolItr_b) )

         ! Pools may in general contain dimensions, namelist options, fields, or other pools,
         ! so we select only those members of the pool that are fields
         if (poolItr_b % memberType == MPAS_POOL_FIELD) then

            ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
            if (poolItr_b % dataType == MPAS_POOL_REAL) then

             call mpas_pool_begin_iteration(pool_a)
             do while ( mpas_pool_get_next_member(pool_a, poolItr_a) )

               !1. start from pool field name
               targetName = trim(poolItr_a % memberName)
               !2. If exists io_name, update with that
               if ( geom%has_io_name(targetName) ) then
                   ioName = geom%io_name(targetName)
               else
                   ioName = 'none'
               end if
               !3. If exists io_scaling_factor, 
               if ( geom%has_io_scaling_factor(targetName) ) then
                   scaling_factor = real(geom%io_scaling_factor(targetName),RKIND)
               else
                   scaling_factor = 0.0_RKIND
               end if

               if ( ( trim(poolItr_b % memberName).eq.trim(targetName) ) .or. &
                    ( trim(poolItr_b % memberName).eq.trim(ioName)     ) ) then
                  ! Depending on the dimensionality of the field, we need to set pointers of
                  ! the correct type
                  if (poolItr_b % nDims == 0) then
                     call mpas_pool_get_array(pool_a, trim(poolItr_a % memberName), r0d_ptr_a)
                     call mpas_pool_get_array(pool_b, trim(poolItr_b % memberName), r0d_ptr_b)
                     r0d_ptr_b = r0d_ptr_a
                  else if (poolItr_b % nDims == 1) then
                     call mpas_pool_get_array(pool_a, trim(poolItr_a % memberName), r1d_ptr_a)
                     call mpas_pool_get_array(pool_b, trim(poolItr_b % memberName), r1d_ptr_b)
                     if ( scaling_factor /= 0.0_RKIND ) then
                        r1d_ptr_b = r1d_ptr_a / scaling_factor
                     else
                        r1d_ptr_b = r1d_ptr_a
                     end if
                     write(message,*) 'Copy sub2all field MIN/MAX: ',trim(poolItr_b % memberName), &
                                      minval(r1d_ptr_a),maxval(r1d_ptr_a)
                     call fckit_log%debug(message)
                  else if (poolItr_b % nDims == 2) then
                     call mpas_pool_get_array(pool_a, trim(poolItr_a % memberName), r2d_ptr_a)
                     call mpas_pool_get_array(pool_b, trim(poolItr_b % memberName), r2d_ptr_b)
                     if ( scaling_factor /= 0.0_RKIND ) then
                        r2d_ptr_b = r2d_ptr_a / scaling_factor
                     else
                        r2d_ptr_b = r2d_ptr_a
                     end if
                     write(message,*) 'Copy sub2all field MIN/MAX: ',trim(poolItr_b % memberName), &
                                      minval(r2d_ptr_a),maxval(r2d_ptr_a)
                     call fckit_log%debug(message)
                  end if

               else if ( match_scalar(trim(poolItr_b % memberName), trim(ioName)) ) then ! Here we know qx or nx variables are all defined as ioName
                  write(message,*) 'Copy sub2all field: Looking at SCALARS now',trim(poolItr_a % memberName)
                  call fckit_log%debug(message)
                  call mpas_pool_get_dimension(state, 'index_'//trim(ioName), index_scalar)
                  if (index_scalar .gt. 0) then
                     call mpas_pool_get_field(pool_a, trim(poolItr_a % memberName), field2d)
                     call mpas_pool_get_field(pool_b, trim(poolItr_b % memberName), field3d)
                     if ( scaling_factor /= 0.0_RKIND ) then
                        field3d % array(index_scalar,:,:) = field2d % array(:,:) / scaling_factor
                     else
                        field3d % array(index_scalar,:,:) = field2d % array(:,:)
                     end if
                     write(message,*) 'Copy sub2all field MIN/MAX: ',trim(poolItr_a % memberName), &
                                      minval(field2d % array), maxval(field2d % array)
                     call fckit_log%debug(message)
                  else
                     write(message,*) 'WARNING in Copy sub2all field; ',trim(poolItr_a % memberName), &
                                      'not available from MPAS'
                     call fckit_log%debug(message)
                  end if
!                       end if
               end if
            end do
           end if
         end if
      end do

   end subroutine da_copy_sub2all_fields


   !***********************************************************************
   !
   !  subroutine da_template_pool
   !
   !> \brief   Subset a pool from fields described in geom
   !> \details
   !>  Given an mpas_geom object, create templatePool containing a subset
   !>  of the fields in the geom % domain's allFields pool and the
   !>  geom % templated_fields
   !
   !-----------------------------------------------------------------------
   subroutine da_template_pool(geom, templatePool, nf, fieldnames)

      implicit none

      ! Arguments
      type (mpas_geom), intent(in) :: geom
      type (mpas_pool_type), pointer, intent(out) :: templatePool
      integer, intent(in) :: nf
      character (len=*), intent(in) :: fieldnames(nf)

      ! Local variables
      type (mpas_pool_type), pointer :: allFields, state
      character(len=MAXVARLEN) :: fieldname, template
      !type (field2DReal), pointer :: qField
      !type (field3DReal), pointer :: scalars

      integer, pointer :: index_scalar, dim0d
      integer :: ii
      integer, parameter :: ndims=10
      character(len=ShortStrKIND) :: dimnames(ndims)

      call mpas_pool_create_pool(templatePool)

      dimnames ( 1) = 'nCellsSolve'
      dimnames ( 2) = 'nEdgesSolve'
      dimnames ( 3) = 'nVerticesSolve'
      dimnames ( 4) = 'nVertLevels'
      dimnames ( 5) = 'nVertLevelsP1'
      dimnames ( 6) = 'nSoilLevels'
      dimnames ( 7) = 'nCells'
      dimnames ( 8) = 'nEdges'
      dimnames ( 9) = 'nVertices'
      dimnames (10) = 'vertexDegree'

      do ii = 1, ndims
         call mpas_pool_get_dimension(geom % domain % blocklist % dimensions, trim(dimnames(ii)), dim0d)
         call mpas_pool_add_dimension(templatePool, trim(dimnames(ii)), dim0d)
      end do

      call mpas_pool_get_subpool(geom % domain % blocklist % structs, 'state',state)
      allFields => geom % domain % blocklist % allFields

      do ii=1, nf
         fieldname = fieldnames(ii)

         if (field_is_scalar(fieldname)) then
            ! Check if this scalar is activated
            call mpas_pool_get_dimension(state, 'index_'//trim(fieldname), index_scalar)
            if (index_scalar .gt. 0) then
               call mpas_pool_template_field('theta', allFields, fieldname, templatePool)
            else
               write(message,*)'--> da_template_pool: ',trim(fieldname), &
                           ' not available in MPAS domain'
               call abor1_ftn(message)
            end if
         else
            if (pool_has_field(allFields, fieldname)) then
               template = fieldname
            else if (geom % is_templated(fieldname)) then
               template = geom%template(fieldname)
               if (.not.pool_has_field(allFields, template)) then
                  write(message,*)'--> da_template_pool: ',trim(template), &
                              ' not available in MPAS domain'
                  call abor1_ftn(message)
               end if
            else
               write(message,*)'--> da_template_pool: ',trim(fieldname), &
                           ' not available in MPAS domain or geom % templated_fields'
               call abor1_ftn(message)
            end if

            call mpas_pool_template_field(template, allFields, fieldname, templatePool)

         end if
      end do

      call da_constant(templatePool, MPAS_JEDI_ZERO_kr)

   end subroutine da_template_pool


   !***********************************************************************
   !
   !  subroutine mpas_pool_template_field
   !
   !> \brief   Add a field to dstPool that is templated on a srcPool field
   !> \details
   !>  Duplicate srcFieldName from srcPool with the name dstFieldName
   !>  in dstPool
   !
   !-----------------------------------------------------------------------
   subroutine mpas_pool_template_field(srcFieldName, srcPool, dstFieldName, dstPool)

      ! Arguments
      type (mpas_pool_type), pointer, intent(in) :: srcPool
      character (len=*), intent(in) :: srcFieldName, dstFieldName
      type (mpas_pool_type), pointer, intent(inout) :: dstPool

      ! Local variables
      type (mpas_pool_iterator_type) :: poolItr
      type (field0DReal), pointer :: field0d_src, field0d_dst
      type (field1DReal), pointer :: field1d_src, field1d_dst
      type (field2DReal), pointer :: field2d_src, field2d_dst
      type (field3DReal), pointer :: field3d_src, field3d_dst, field3d
      type (field0DInteger), pointer :: ifield0d_src, ifield0d_dst
      type (field1DInteger), pointer :: ifield1d_src, ifield1d_dst
      type (field2DInteger), pointer :: ifield2d_src, ifield2d_dst
      type (field3DInteger), pointer :: ifield3d_src, ifield3d_dst

      !
      ! Iterate over srcPool, add the one that matches srcFieldName into dstPool
      !
      call mpas_pool_begin_iteration(srcPool)

      do while ( mpas_pool_get_next_member(srcPool, poolItr) )

         ! Only handle fields, ignore dimensions, namelist options, and other pools
         if ( poolItr % memberType == MPAS_POOL_FIELD .and. &
              trim(srcFieldName) == trim(poolItr % memberName) ) then

            ! Handle real and integer, ignore logical and char
            if (poolItr % dataType == MPAS_POOL_REAL) then

               ! Use correctly dimensioned field pointers
               if (poolItr % nDims == 0) then
                  call mpas_pool_get_field(srcPool, trim(srcFieldName), field0d_src)
                  call mpas_duplicate_field(field0d_src, field0d_dst)
                  field0d_dst % fieldName = dstFieldName
                  call mpas_pool_add_field(dstPool, trim(dstFieldName), field0d_dst)
               else if (poolItr % nDims == 1) then
                  call mpas_pool_get_field(srcPool, trim(srcFieldName), field1d_src)
                  call mpas_duplicate_field(field1d_src, field1d_dst)
                  field1d_dst % fieldName = dstFieldName
                  call mpas_pool_add_field(dstPool, trim(dstFieldName), field1d_dst)
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_field(srcPool, trim(srcFieldName), field2d_src)
                  call mpas_duplicate_field(field2d_src, field2d_dst)
                  field2d_dst % fieldName = dstFieldName
                  call mpas_pool_add_field(dstPool, trim(dstFieldName), field2d_dst)
               else if (poolItr % nDims == 3) then
                  call mpas_pool_get_field(srcPool, trim(srcFieldName), field3d_src)
                  call mpas_duplicate_field(field3d_src, field3d_dst)
                  field3d_dst % fieldName = dstFieldName
                  call mpas_pool_add_field(dstPool, trim(dstFieldName), field3d_dst)
               end if

            else if (poolItr % dataType == MPAS_POOL_INTEGER) then

               ! Use correctly dimensioned field pointers
               if (poolItr % nDims == 0) then
                  call mpas_pool_get_field(srcPool, trim(srcFieldName), ifield0d_src)
                  call mpas_duplicate_field(ifield0d_src, ifield0d_dst)
                  ifield0d_dst % fieldName = dstFieldName
                  call mpas_pool_add_field(dstPool, trim(dstFieldName), ifield0d_dst)
               else if (poolItr % nDims == 1) then
                  call mpas_pool_get_field(srcPool, trim(srcFieldName), ifield1d_src)
                  call mpas_duplicate_field(ifield1d_src, ifield1d_dst)
                  ifield1d_dst % fieldName = dstFieldName
                  call mpas_pool_add_field(dstPool, trim(dstFieldName), ifield1d_dst)
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_field(srcPool, trim(srcFieldName), ifield2d_src)
                  call mpas_duplicate_field(ifield2d_src, ifield2d_dst)
                  ifield2d_dst % fieldName = dstFieldName
                  call mpas_pool_add_field(dstPool, trim(dstFieldName), ifield2d_dst)
               else if (poolItr % nDims == 3) then
                  call mpas_pool_get_field(srcPool, trim(srcFieldName), ifield3d_src)
                  call mpas_duplicate_field(ifield3d_src, ifield3d_dst)
                  ifield3d_dst % fieldName = dstFieldName
                  call mpas_pool_add_field(dstPool, trim(dstFieldName), ifield3d_dst)
               end if
            end if

            return

         end if
      end do

      write(message,*)'--> mpas_pool_template_field: ',trim(srcFieldName), &
                  ' not available in srcPool'
      call abor1_ftn(message)

   end subroutine mpas_pool_template_field


   !***********************************************************************
   !
   !  function da_common_vars
   !
   !> \author  Gael Descombes
   !> \date    26 December 2017
   !> \details
   !>  Count the number of fields in a Pool related to a list of fields
   !
   !-----------------------------------------------------------------------
   function da_common_vars(pool_a, fieldname) result(nsize0)

      implicit none
      type (mpas_pool_type), pointer :: pool_a, pool_b
      type (mpas_pool_iterator_type) :: poolItr
      character (len=*) :: fieldname(:)
      integer :: ii, jj, nsize, nsize0

      nsize0 = 0
      nsize  = size(fieldname)
      call mpas_pool_begin_iteration(pool_a)

         do while ( mpas_pool_get_next_member(pool_a, poolItr) )
            ! Pools may in general contain dimensions, namelist options, fields, or other pools,
            ! so we select only those members of the pool that are fields
            if (poolItr % memberType == MPAS_POOL_FIELD) then
               ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
               if (poolItr % dataType == MPAS_POOL_REAL) then
                  do ii=1, nsize
                     if ( trim(fieldname(ii)).eq.(trim(poolItr % memberName)) ) then
                        nsize0 = nsize0 + 1
                     else if (match_scalar(trim(poolItr % memberName), trim(fieldname(ii)))) then
                        nsize0 = nsize0 + 1
                     end if
                  end do
               end if
            end if
         end do

   end function da_common_vars



   !***********************************************************************
   !
   !  subroutine da_random
   !
   !> \brief   Performs random for pool A
   !> \author  Gael Descombes
   !> \date    January 2018
   !> \details
   !
   !-----------------------------------------------------------------------
   subroutine da_random(pool_a, fld_select)

      implicit none

      type (mpas_pool_type), pointer, intent(inout) :: pool_a
      character (len=*), optional,    intent(in)    :: fld_select(:)

      integer, parameter :: rseed = 7
      type (mpas_pool_iterator_type) :: poolItr
      real (kind=RKIND), pointer :: r0d_ptr_a
      real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a
      real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a
      real (kind=RKIND), dimension(:,:,:), pointer :: r3d_ptr_a
      !loc
      real (kind=RKIND) :: zero, one

      !
      ! Iterate over all fields in pool_b, adding them to fields of the same
      ! name in pool_a
      !
      call mpas_pool_begin_iteration(pool_a)
      one=real(1.0, RKIND); zero=real(0.0, RKIND)

      do while ( mpas_pool_get_next_member(pool_a, poolItr) )

         if (present(fld_select)) then
            if (ufo_vars_getindex(fld_select,trim(poolItr % memberName)) < 0) cycle
         end if

         ! Pools may in general contain dimensions, namelist options, fields, or other pools,
         ! so we select only those members of the pool that are fields
         if (poolItr % memberType == MPAS_POOL_FIELD) then

            ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
            if (poolItr % dataType == MPAS_POOL_REAL) then

               ! Depending on the dimensionality of the field, we need to set pointers of
               ! the correct type
               if (poolItr % nDims == 0) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r0d_ptr_a)
               else if (poolItr % nDims == 1) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
                  call normal_distribution(r1d_ptr_a, zero, one, rseed)
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_a)
                  call normal_distribution(r2d_ptr_a, zero, one, rseed)
               else if (poolItr % nDims == 3) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r3d_ptr_a)
                  call normal_distribution(r3d_ptr_a, zero, one, rseed)
               end if

            end if
         end if
      end do

   end subroutine da_random

   !-----------------------------------------------------------------------
   !  subroutine da_operator
   !
   !> \brief   Performs A = A 'kind_op' B for pools A and B
   !> \author  Michael Duda
   !> \date    20 December 2017
   !> \details
   !>  Given two pools, A and B, where the fields in B are a subset of
   !>  the fields in A, this routine adds the fields in B to fields in A
   !>  with the same name. When A and B contain identical fields, this
   !>  is equivalent to A = A 'kind_op' B.
   !>  \modified by Gael DESCOMBES to apply diffferent operator
   !
   !-----------------------------------------------------------------------
   subroutine da_operator(kind_op, pool_a, pool_b, pool_c, fld_select)

      implicit none

      type (mpas_pool_type), pointer :: pool_a, pool_b
      type (mpas_pool_type), pointer, optional :: pool_c
      character (len=*) :: kind_op
      character (len=*), optional :: fld_select(:)

      type (mpas_pool_iterator_type) :: poolItr
      real (kind=RKIND), pointer :: r0d_ptr_a, r0d_ptr_b, r0d_ptr_c
      real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b, r1d_ptr_c
      real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b, r2d_ptr_c
      real (kind=RKIND), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b, r3d_ptr_c
      !
      ! Iterate over all fields in pool_b, adding them to fields of the same
      ! name in pool_a
      !
      call mpas_pool_begin_iteration(pool_b)

      write(message,*) 'Operator ',trim(kind_op)
      call fckit_log%debug(message)

      do while ( mpas_pool_get_next_member(pool_b, poolItr) )

         if (present(fld_select)) then
            if (ufo_vars_getindex(fld_select,trim(poolItr % memberName)) < 0) cycle
         end if

         ! Pools may in general contain dimensions, namelist options, fields, or other pools,
         ! so we select only those members of the pool that are fields
         if (poolItr % memberType == MPAS_POOL_FIELD) then

            ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
            if (poolItr % dataType == MPAS_POOL_REAL) then

               ! Depending on the dimensionality of the field, we need to set pointers of
               ! the correct type
               if (poolItr % nDims == 0) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r0d_ptr_a)
                  call mpas_pool_get_array(pool_b, trim(poolItr % memberName), r0d_ptr_b)
                  if (present(pool_c)) then
                     call mpas_pool_get_array(pool_c, trim(poolItr % memberName), r0d_ptr_c)
                     r0d_ptr_a = MPAS_JEDI_ZERO_kr
                  end if
                  if ( trim(kind_op).eq.'add' ) then
                     r0d_ptr_a = r0d_ptr_a + r0d_ptr_b
                     if (present(pool_c)) then
                       r0d_ptr_a = r0d_ptr_b + r0d_ptr_c
                     else
                        r0d_ptr_a = r0d_ptr_a + r0d_ptr_b
                     end if
                  else if ( trim(kind_op).eq.'schur' ) then
                     if (present(pool_c)) then
                        r0d_ptr_a = r0d_ptr_b * r0d_ptr_c
                     else
                        r0d_ptr_a = r0d_ptr_a * r0d_ptr_b
                     end if
                  else if ( trim(kind_op).eq.'sub' ) then
                     if (present(pool_c)) then
                        r0d_ptr_a = r0d_ptr_b - r0d_ptr_c
                     else
                        r0d_ptr_a = r0d_ptr_a - r0d_ptr_b
                     end if
                  end if

               else if (poolItr % nDims == 1) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
                  call mpas_pool_get_array(pool_b, trim(poolItr % memberName), r1d_ptr_b)
                  if (present(pool_c)) then
                     call mpas_pool_get_array(pool_c, trim(poolItr % memberName), r1d_ptr_c)
                     r1d_ptr_a = MPAS_JEDI_ZERO_kr
                  end if
                  if ( trim(kind_op).eq.'add' ) then
                     if (present(pool_c)) then
                        r1d_ptr_a = r1d_ptr_b + r1d_ptr_c
                     else
                        r1d_ptr_a = r1d_ptr_a + r1d_ptr_b
                     end if
                  else if ( trim(kind_op).eq.'schur' ) then
                     if (present(pool_c)) then
                        r1d_ptr_a = r1d_ptr_b * r1d_ptr_c
                     else
                        r1d_ptr_a = r1d_ptr_a * r1d_ptr_b
                     end if
                  else if ( trim(kind_op).eq.'sub' ) then
                     if (present(pool_c)) then
                        r1d_ptr_a = r1d_ptr_b - r1d_ptr_c
                     else
                        r1d_ptr_a = r1d_ptr_a - r1d_ptr_b
                     end if
                  end if

               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_a)
                  call mpas_pool_get_array(pool_b, trim(poolItr % memberName), r2d_ptr_b)
                  if (present(pool_c)) then
                     call mpas_pool_get_array(pool_c, trim(poolItr % memberName), r2d_ptr_c)
                     r2d_ptr_a = MPAS_JEDI_ZERO_kr
                  end if
                  if ( trim(kind_op).eq.'add' ) then
                     write(message,*) 'Operator_a add MIN/MAX: ',minval(r2d_ptr_a),maxval(r2d_ptr_a)
                     call fckit_log%debug(message)
                     write(message,*) 'Operator_b add MIN/MAX: ',minval(r2d_ptr_b),maxval(r2d_ptr_b)
                     call fckit_log%debug(message)
                     if (present(pool_c)) then
                        r2d_ptr_a = r2d_ptr_b + r2d_ptr_c
                     else
                        r2d_ptr_a = r2d_ptr_a + r2d_ptr_b
                     end if
                     write(message,*) 'Operator2 add MIN/MAX: ',minval(r2d_ptr_a),maxval(r2d_ptr_a)
                     call fckit_log%debug(message)
                  else if ( trim(kind_op).eq.'schur' ) then
                     if (present(pool_c)) then
                        r2d_ptr_a = r2d_ptr_b * r2d_ptr_c
                     else
                        r2d_ptr_a = r2d_ptr_a * r2d_ptr_b
                     end if
                  else if ( trim(kind_op).eq.'sub' ) then
                     if (present(pool_c)) then
                        r2d_ptr_a = r2d_ptr_b - r2d_ptr_c
                     else
                        r2d_ptr_a = r2d_ptr_a - r2d_ptr_b
                     end if
                  end if

               else if (poolItr % nDims == 3) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r3d_ptr_a)
                  call mpas_pool_get_array(pool_b, trim(poolItr % memberName), r3d_ptr_b)
                  if (present(pool_c)) then
                     call mpas_pool_get_array(pool_c, trim(poolItr % memberName), r3d_ptr_c)
                     r3d_ptr_a = MPAS_JEDI_ZERO_kr
                  end if
                  if ( trim(kind_op).eq.'add' ) then
                     if (present(pool_c)) then
                        r3d_ptr_a = r3d_ptr_b + r3d_ptr_c
                     else
                        r3d_ptr_a = r3d_ptr_a + r3d_ptr_b
                     end if
                  else if ( trim(kind_op).eq.'schur' ) then
                     if (present(pool_c)) then
                        r3d_ptr_a = r3d_ptr_b * r3d_ptr_c
                     else
                        r3d_ptr_a = r3d_ptr_a * r3d_ptr_b
                     end if
                  else if ( trim(kind_op).eq.'sub' ) then
                     if (present(pool_c)) then
                        r3d_ptr_a = r3d_ptr_b - r3d_ptr_c
                     else
                        r3d_ptr_a = r3d_ptr_a - r3d_ptr_b
                     end if
                  end if
               end if

            end if
         end if
      end do

   end subroutine da_operator

   !***********************************************************************
   !
   !  subroutine da_self_mult
   !
   !> \brief   Performs A = A * zz for pool A, zz a real number
   !> \author  Gael Descombes
   !> \date    22 December 2017
   !> \details
   !
   !-----------------------------------------------------------------------
   subroutine da_self_mult(pool_a, zz)

      implicit none

      type (mpas_pool_type), pointer :: pool_a
      real (kind=RKIND) :: zz

      type (mpas_pool_iterator_type) :: poolItr
      real (kind=RKIND), pointer :: r0d_ptr_a
      real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a
      real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a
      real (kind=RKIND), dimension(:,:,:), pointer :: r3d_ptr_a

      !
      ! Iterate over all fields in pool_b, adding them to fields of the same
      ! name in pool_a
      !
      call mpas_pool_begin_iteration(pool_a)

      do while ( mpas_pool_get_next_member(pool_a, poolItr) )

         ! Pools may in general contain dimensions, namelist options, fields, or other pools,
         ! so we select only those members of the pool that are fields
         if (poolItr % memberType == MPAS_POOL_FIELD) then

            ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
            if (poolItr % dataType == MPAS_POOL_REAL) then

               ! Depending on the dimensionality of the field, we need to set pointers of
               ! the correct type
               if (poolItr % nDims == 0) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r0d_ptr_a)
                  r0d_ptr_a = r0d_ptr_a * zz
               else if (poolItr % nDims == 1) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
                  r1d_ptr_a = r1d_ptr_a * zz
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_a)
                  r2d_ptr_a = r2d_ptr_a * zz
               else if (poolItr % nDims == 3) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r3d_ptr_a)
                  r3d_ptr_a = r3d_ptr_a * zz
               end if

            end if
         end if
      end do

   end subroutine da_self_mult


   !***********************************************************************
   !
   !  subroutine da_constant
   !
   !> \brief   Performs A = constant. for pool A
   !> \author  Gael Descombes
   !> \date    22 December 2017
   !> \details
   !
   !-----------------------------------------------------------------------
   subroutine da_constant(pool_a, realvalue, fld_select)

      implicit none

      type (mpas_pool_type), pointer, intent(inout) :: pool_a
      real (kind=kind_real),          intent(in)    :: realvalue
      character (len=*), optional,    intent(in)    :: fld_select(:)

      type (mpas_pool_iterator_type) :: poolItr
      real (kind=RKIND), pointer :: r0d_ptr_a
      real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a
      real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a
      real (kind=RKIND), dimension(:,:,:), pointer :: r3d_ptr_a

      call mpas_pool_begin_iteration(pool_a)

      do while ( mpas_pool_get_next_member(pool_a, poolItr) )

         if (present(fld_select)) then
            if (ufo_vars_getindex(fld_select,trim(poolItr % memberName)) < 0) cycle
         end if

         ! Pools may in general contain dimensions, namelist options, fields, or other pools,
         ! so we select only those members of the pool that are fields
         if (poolItr % memberType == MPAS_POOL_FIELD) then

            ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
            if (poolItr % dataType == MPAS_POOL_REAL) then

               ! Depending on the dimensionality of the field, we need to set pointers of
               ! the correct type
               if (poolItr % nDims == 0) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r0d_ptr_a)
                  r0d_ptr_a = realvalue
               else if (poolItr % nDims == 1) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
                  r1d_ptr_a = realvalue
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_a)
                  r2d_ptr_a = realvalue
               else if (poolItr % nDims == 3) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r3d_ptr_a)
                  r3d_ptr_a = realvalue
               end if

            end if
         end if
      end do

   end subroutine da_constant

   !***********************************************************************
   !
   !  subroutine da_posdef
   !
   !> \brief   Performs A = max(0.,A) for pool A
   !> \author  JJ Guerrette
   !> \date    12 July 2019
   !> \details
   !
   !-----------------------------------------------------------------------
   subroutine da_posdef(pool_a, fld_select)

      implicit none

      type (mpas_pool_type), pointer, intent(inout) :: pool_a
      character (len=*), optional,    intent(in)    :: fld_select(:)

      type (mpas_pool_iterator_type) :: poolItr
      real (kind=RKIND), pointer :: r0d_ptr_a
      real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a
      real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a
      real (kind=RKIND), dimension(:,:,:), pointer :: r3d_ptr_a

      !
      ! Iterate over all fields in pool_b, adding them to fields of the same
      ! name in pool_a
      !
      call mpas_pool_begin_iteration(pool_a)

      do while ( mpas_pool_get_next_member(pool_a, poolItr) )

         if (present(fld_select)) then
            if (ufo_vars_getindex(fld_select,trim(poolItr % memberName)) < 0) cycle
         end if

         ! Pools may in general contain dimensions, namelist options, fields, or other pools,
         ! so we select only those members of the pool that are fields
         if (poolItr % memberType == MPAS_POOL_FIELD) then

            ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
            if (poolItr % dataType == MPAS_POOL_REAL) then

               ! Depending on the dimensionality of the field, we need to set pointers of
               ! the correct type
               if (poolItr % nDims == 0) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r0d_ptr_a)
                  r0d_ptr_a = max(MPAS_JEDI_ZERO_kr, r0d_ptr_a)
               else if (poolItr % nDims == 1) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
                  r1d_ptr_a = max(MPAS_JEDI_ZERO_kr, r1d_ptr_a)
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_a)
                  r2d_ptr_a = max(MPAS_JEDI_ZERO_kr, r2d_ptr_a)
               else if (poolItr % nDims == 3) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r3d_ptr_a)
                  r3d_ptr_a = max(MPAS_JEDI_ZERO_kr, r3d_ptr_a)
               end if

            end if
         end if
      end do

   end subroutine da_posdef

   !***********************************************************************
   !
   !  subroutine da_setval
   !
   !> \brief   Performs A = Val_R. for pool A
   !> \author  Gael Descombes
   !> \date    22 December 2017
   !> \details
   !
   !-----------------------------------------------------------------------
   subroutine da_setval(pool_a,zz)

      implicit none

      type (mpas_pool_type), pointer :: pool_a
      real (kind=kind_real) :: zz

      type (mpas_pool_iterator_type) :: poolItr
      real (kind=RKIND), pointer :: r0d_ptr_a
      real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a
      real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a
      real (kind=RKIND), dimension(:,:,:), pointer :: r3d_ptr_a

      !
      ! Iterate over all fields in pool_b, adding them to fields of the same
      ! name in pool_a
      !
      call mpas_pool_begin_iteration(pool_a)

      do while ( mpas_pool_get_next_member(pool_a, poolItr) )

         ! Pools may in general contain dimensions, namelist options, fields, or other pools,
         ! so we select only those members of the pool that are fields
         if (poolItr % memberType == MPAS_POOL_FIELD) then

            ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
            if (poolItr % dataType == MPAS_POOL_REAL) then

               ! Depending on the dimensionality of the field, we need to set pointers of
               ! the correct type
               if (poolItr % nDims == 0) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r0d_ptr_a)
                  r0d_ptr_a = zz
               else if (poolItr % nDims == 1) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
                  r1d_ptr_a = zz
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_a)
                  r2d_ptr_a = zz
               else if (poolItr % nDims == 3) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r3d_ptr_a)
                  r3d_ptr_a = zz
               end if

            end if
         end if
      end do

   end subroutine da_setval





   !***********************************************************************
   !
   !  subroutine da_axpy
   !
   !> \brief   Performs A = A + B * zz for pools A and B
   !> \author  Gael Descombes
   !> \date    20 December 2017
   !> \details
   !>  Given two pools, A and B, where the fields in B are a subset of
   !>  the fields in A, this routine adds the fields in B to fields in A
   !>  with the same name. When A and B contain identical fields, this
   !>  is equivalent to A = A + B.
   !
   !-----------------------------------------------------------------------
   subroutine da_axpy(pool_a, pool_b, zz, fld_select)

      implicit none
      type (mpas_pool_type), pointer, intent(inout) :: pool_a
      type (mpas_pool_type), pointer, intent(in)    :: pool_b
      real (kind=RKIND), intent(in) :: zz
      character (len=*), optional, intent(in) :: fld_select(:)



      type (mpas_pool_iterator_type) :: poolItr
      real (kind=RKIND), pointer :: r0d_ptr_a, r0d_ptr_b
      real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
      real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
      real (kind=RKIND), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b

      !
      ! Iterate over all fields in pool_b, adding them to fields of the same
      ! name in pool_a
      !
      call mpas_pool_begin_iteration(pool_b)

      do while ( mpas_pool_get_next_member(pool_b, poolItr) )

         if (present(fld_select)) then
            if (ufo_vars_getindex(fld_select,trim(poolItr % memberName)) < 0) cycle
         end if

         ! Pools may in general contain dimensions, namelist options, fields, or other pools,
         ! so we select only those members of the pool that are fields
         if (poolItr % memberType == MPAS_POOL_FIELD) then

            ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
            if (poolItr % dataType == MPAS_POOL_REAL) then

               ! Depending on the dimensionality of the field, we need to set pointers of
               ! the correct type
               if (poolItr % nDims == 0) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r0d_ptr_a)
                  call mpas_pool_get_array(pool_b, trim(poolItr % memberName), r0d_ptr_b)
                  r0d_ptr_a = r0d_ptr_a + r0d_ptr_b * zz
               else if (poolItr % nDims == 1) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
                  call mpas_pool_get_array(pool_b, trim(poolItr % memberName), r1d_ptr_b)
                  r1d_ptr_a = r1d_ptr_a + r1d_ptr_b * zz
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_a)
                  call mpas_pool_get_array(pool_b, trim(poolItr % memberName), r2d_ptr_b)
                  r2d_ptr_a = r2d_ptr_a + r2d_ptr_b * zz
               else if (poolItr % nDims == 3) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r3d_ptr_a)
                  call mpas_pool_get_array(pool_b, trim(poolItr % memberName), r3d_ptr_b)
                  r3d_ptr_a = r3d_ptr_a + r3d_ptr_b * zz
               end if

            end if
         end if
      end do

   end subroutine da_axpy


   !***********************************************************************
   !
   !  subroutine da_gpnorm
   !
   !> \brief   Performs basic statistics min/max/norm given a pool
   !> \author  Gael Descombes
   !> \date    February 2018
   !> \details
   !>  Given a pool of fields, return min/max/norm array
   !
   !-----------------------------------------------------------------------

   subroutine da_gpnorm(pool_a, dminfo, nf, pstat, fld_select)

   implicit none
   type (mpas_pool_type), pointer, intent(in)  :: pool_a
   type (dm_info), pointer,        intent(in)  :: dminfo
   integer,                        intent(in)  :: nf
   character (len=*),              intent(in)  :: fld_select(nf)
   real(kind=RKIND),               intent(out) :: pstat(3, nf)

   type (mpas_pool_iterator_type) :: poolItr
   type (field1DInteger), pointer :: ifield1d
   type (field1DReal), pointer :: field1d
   type (field2DReal), pointer :: field2d
   type (field3DReal), pointer :: field3d
   real(kind=RKIND) :: globalSum, globalMin, globalMax, dimtot, dimtot_global, prodtot

   integer :: jj, ndims
   integer :: dim1, dim2, dim3
   integer, allocatable :: dimSizes(:)

   pstat = MPAS_JEDI_ZERO_kr

   !
   ! Iterate over all fields in pool_a
   ! name in pool_a
   !
   call mpas_pool_begin_iteration(pool_a)

   do while ( mpas_pool_get_next_member(pool_a, poolItr) )
      jj = ufo_vars_getindex(fld_select,trim(poolItr % memberName))
      if ( jj < 0 .or. jj > nf ) cycle

      ! Pools may in general contain dimensions, namelist options, fields, or other pools,
      ! so we select only those members of the pool that are fields
      if (poolItr % memberType == MPAS_POOL_FIELD) then

         ndims = poolItr % nDims
         dimSizes = getSolveDimSizes(pool_a, poolItr%memberName)

         ! Fields can be integer, logical, or real. Here, we operate only on real and integer fields.
         if (poolItr % dataType == MPAS_POOL_REAL) then

            ! Depending on the dimensionality of the field, we need to set pointers of
            ! the correct type
            if (ndims == 1) then
               dim1 = dimSizes(1)
               call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field1d)
               dimtot = real(dim1,kind_real)
               prodtot = sum(field1d % array(1:dim1)**2 )
               call mpas_dmpar_sum_real(dminfo, dimtot, dimtot_global)
               call mpas_dmpar_sum_real(dminfo, prodtot, globalSum)
               call mpas_dmpar_min_real(dminfo, minval(field1d % array(1:dim1)), globalMin)
               call mpas_dmpar_max_real(dminfo, maxval(field1d % array(1:dim1)), globalMax)
               pstat(1,jj) = globalMin
               pstat(2,jj) = globalMax
               pstat(3,jj) = sqrt( globalSum / dimtot_global )
            else if (ndims == 2) then
               dim1 = dimSizes(1)
               dim2 = dimSizes(2)
               call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field2d)
               dimtot  = real(dim1*dim2,kind_real)
               prodtot = sum(field2d % array(1:dim1,1:dim2)**2 )
               call mpas_dmpar_sum_real(dminfo, dimtot, dimtot_global)
               call mpas_dmpar_sum_real(dminfo, prodtot, globalSum)
               call mpas_dmpar_min_real(dminfo, minval(field2d % array(1:dim1,1:dim2)), globalMin)
               call mpas_dmpar_max_real(dminfo, maxval(field2d % array(1:dim1,1:dim2)), globalMax)
               pstat(1,jj) = globalMin
               pstat(2,jj) = globalMax
               pstat(3,jj) = sqrt( globalSum / dimtot_global )
            else if (ndims == 3) then
               dim1 = dimSizes(1)
               dim2 = dimSizes(2)
               dim3 = dimSizes(3)
               call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field3d)
               dimtot  = real(dim1*dim2*dim3,kind_real)
               prodtot = sum(field3d % array(1:dim1,1:dim2,1:dim3)**2 )
               call mpas_dmpar_sum_real(dminfo, dimtot, dimtot_global)
               call mpas_dmpar_sum_real(dminfo, prodtot, globalSum)
               call mpas_dmpar_min_real(dminfo, minval(field3d % array(1:dim1,1:dim2,1:dim3)), globalMin)
               call mpas_dmpar_max_real(dminfo, maxval(field3d % array(1:dim1,1:dim2,1:dim3)), globalMax)
               pstat(1,jj) = globalMin
               pstat(2,jj) = globalMax
               pstat(3,jj) = sqrt( globalSum / dimtot_global )
            end if
            deallocate(dimSizes)

         else if (poolItr % dataType == MPAS_POOL_INTEGER) then

            ! For now, we only handle the 1-dimensional integer variables.
            if (ndims == 1) then
               dim1 = dimSizes(1)
               call mpas_pool_get_field(pool_a, trim(poolItr % memberName), ifield1d)
               dimtot = real(dim1,kind_real)
               prodtot = sum( real(ifield1d % array(1:dim1))**2 )
               call mpas_dmpar_sum_real(dminfo, dimtot, dimtot_global)
               call mpas_dmpar_sum_real(dminfo, prodtot, globalSum)
               call mpas_dmpar_min_real(dminfo, real(minval(ifield1d % array(1:dim1)),RKIND), globalMin)
               call mpas_dmpar_max_real(dminfo, real(maxval(ifield1d % array(1:dim1)),RKIND), globalMax)
               pstat(1,jj) = globalMin
               pstat(2,jj) = globalMax
               pstat(3,jj) = sqrt( globalSum / dimtot_global )
            end if

         end if
      end if
   end do

   end subroutine da_gpnorm


   !***********************************************************************
   !
   !  subroutine da_fldrms
   !
   !> \brief   Performs basic statistics min/max/norm given a pool
   !> \author  Gael Descombes
   !> \date    February 2018
   !> \details
   !>  Given a pool of fields, return min/max/norm array
   !
   !-----------------------------------------------------------------------

   subroutine da_fldrms(pool_a, dminfo, fldrms, fld_select)

   implicit none
   type (mpas_pool_type), pointer, intent(in)  :: pool_a
   type (dm_info), pointer,        intent(in)  :: dminfo
   real(kind=RKIND),           intent(out) :: fldrms
   character (len=*), optional,    intent(in)  :: fld_select(:)

   type (mpas_pool_iterator_type) :: poolItr
   type (field1DReal), pointer :: field1d
   type (field2DReal), pointer :: field2d
   type (field3DReal), pointer :: field3d
   real(kind=RKIND) :: dimtot, dimtot_global, prodtot, prodtot_global

   integer :: ndims
   integer :: dim1, dim2, dim3
   integer, allocatable :: dimSizes(:)

   prodtot = MPAS_JEDI_ZERO_kr
   dimtot  = MPAS_JEDI_ZERO_kr

   !
   ! Iterate over all fields in pool_a
   ! named in pool_a
   !
   call mpas_pool_begin_iteration(pool_a)

   do while ( mpas_pool_get_next_member(pool_a, poolItr) )
      if (present(fld_select)) then
         if (ufo_vars_getindex(fld_select,trim(poolItr % memberName)) < 0) cycle
      end if
      if (poolItr % dataType == MPAS_POOL_REAL) then
         if (poolItr % memberType == MPAS_POOL_FIELD) then
            ndims = poolItr % nDims
            dimSizes = getSolveDimSizes(pool_a, poolItr%memberName)
            if (ndims == 1) then
               dim1 = dimSizes(1)
               call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field1d)
               dimtot  = dimtot + real(dim1,kind_real)
               prodtot = prodtot + sum( field1d % array(1:dim1)**2 )
            else if (ndims == 2) then
               dim1 = dimSizes(1)
               dim2 = dimSizes(2)
               call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field2d)
               dimtot  = dimtot + real(dim1*dim2,kind_real)
               prodtot = prodtot + sum( field2d % array(1:dim1,1:dim2)**2 )
            else if (ndims == 3) then
               dim1 = dimSizes(1)
               dim2 = dimSizes(2)
               dim3 = dimSizes(3)
               call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field3d)
               dimtot  = dimtot + real(dim1*dim2*dim3,kind_real)
               prodtot = prodtot + sum( field3d % array(1:dim1,1:dim2,1:dim3)**2 )
            end if
            deallocate(dimSizes)
         end if
      end if
   end do

   call mpas_dmpar_sum_real(dminfo, dimtot, dimtot_global)
   call mpas_dmpar_sum_real(dminfo, prodtot, prodtot_global)
   fldrms = sqrt(prodtot_global / dimtot_global)

   end subroutine da_fldrms


   !***********************************************************************
   !
   !  subroutine da_dot_product
   !
   !> \brief   Performs the dot_product given two pools of fields
   !> \author  Gael Descombes
   !> \date    February 2018
   !> \details
   !>  Given two pools of fields, compute the dot_product
   !
   !-----------------------------------------------------------------------

   subroutine da_dot_product(pool_a, pool_b, dminfo, zprod)

   implicit none
   type (mpas_pool_type), pointer, intent(in)  :: pool_a, pool_b
   type (dm_info), pointer,        intent(in)  :: dminfo
   real(kind=RKIND),           intent(out) :: zprod

   type (mpas_pool_iterator_type) :: poolItr
   type (field1DReal), pointer :: field1d_a, field1d_b
   type (field2DReal), pointer :: field2d_a, field2d_b
   type (field3DReal), pointer :: field3d_a, field3d_b
   real(kind=RKIND) :: fieldSum_local, zprod_local

   integer :: ndims
   integer :: dim1, dim2, dim3
   integer, allocatable :: dimSizes(:)

   !
   ! Iterate over all fields in pool_a
   ! named in pool_a
   !
   call mpas_pool_begin_iteration(pool_a)

   zprod_local = MPAS_JEDI_ZERO_kr

   do while ( mpas_pool_get_next_member(pool_a, poolItr) )
      if (poolItr % dataType == MPAS_POOL_REAL) then
         if (poolItr % memberType == MPAS_POOL_FIELD) then
            ndims = poolItr % nDims
            dimSizes = getSolveDimSizes(pool_a, poolItr%memberName)
            !TODO: add check that dimSizes are the same between pool_a and pool_b for poolItr
            if (ndims == 1) then
               dim1 = dimSizes(1)
               call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field1d_a)
               call mpas_pool_get_field(pool_b, trim(poolItr % memberName), field1d_b)
               fieldSum_local = sum(field1d_a % array(1:dim1) * field1d_b % array(1:dim1))
               zprod_local = zprod_local + fieldSum_local
            else if (ndims == 2) then
               dim1 = dimSizes(1)
               dim2 = dimSizes(2)
               call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field2d_a)
               call mpas_pool_get_field(pool_b, trim(poolItr % memberName), field2d_b)
               fieldSum_local = sum(field2d_a % array(1:dim1,1:dim2) &
                                  * field2d_b % array(1:dim1,1:dim2))
               zprod_local = zprod_local + fieldSum_local
            else if (ndims == 3) then
               dim1 = dimSizes(1)
               dim2 = dimSizes(2)
               dim3 = dimSizes(3)
               call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field3d_a)
               call mpas_pool_get_field(pool_b, trim(poolItr % memberName), field3d_b)
               fieldSum_local = sum(field3d_a % array(1:dim1,1:dim2,1:dim3) &
                                  * field3d_b % array(1:dim1,1:dim2,1:dim3))
               zprod_local = zprod_local + fieldSum_local
            end if
            deallocate(dimSizes)
         end if
      end if
   end do

   call mpas_dmpar_sum_real(dminfo, zprod_local, zprod)
   end subroutine da_dot_product


  subroutine cvt_oopsmpas_date(inString2,outString2,iconv)

     implicit none

     character (len=*), intent(in) :: inString2
     character (len=*), intent(inout) :: outString2
     integer, intent(in) :: iconv
     integer :: i, curLen
     integer :: year, month, day, hour, minute, second

     character (len=ShortStrKIND) :: timePart
     character (len=ShortStrKIND) :: yearFormat
     logical :: charExpand
     character (len=4) :: YYYY
     character (len=2) :: MM, DD, h, m, s
     character (len=21) :: outString, inString

     ! 2017-08-08T00:00:00Z OOPS/YAML format
     ! 2010-10-24_02.00.00  MPAS format
     ! iconv=1: MPAS --> OOPS/YAML
     ! iconv=-1: OOPS/YAML --> MPAS

     if (iconv.eq.-1) then
        YYYY = inString2(1:4)
        MM   = inString2(6:7)
        DD   = inString2(9:10)
        h    = inString2(12:13)
        m    = inString2(15:16)
        s    = inString2(18:19)
     else
        YYYY = inString2(1:4)
        MM   = inString2(6:7)
        DD   = inString2(9:10)
        h    = inString2(12:13)
        m    = inString2(15:16)
        s    = inString2(18:19)
     end if

     write(message,*) 'cvt_oopsmpas_date instring: ',trim(YYYY),trim(MM),trim(DD),trim(h),trim(m),trim(s)
     call fckit_log%debug(message)
     write(message,*) 'cvt_oopsmpas_date input ',trim(instring2)
     call fckit_log%debug(message)

     write(outString,*) ''
     instring = trim(outstring2)

     curLen = 0
     charExpand = .false.
     do i = 1, len_trim(inString)
           if (inString(i:i) == '$' ) then
               charExpand = .true.
           else if (inString(i:i) /= '$') then
               write(message,*) 'inString: ',trim(inString(i:i)),charExpand
               call fckit_log%debug(message)
               if (charExpand) then
                  select case (inString(i:i))
                     case ('Y')
                         outString = trim(outString) // trim(YYYY)
                     case ('M')
                         outString = trim(outString) // trim(MM)
                     case ('D')
                         outString = trim(outString) // trim(DD)
                     case ('h')
                         outString = trim(outString) // trim(h)
                     case ('m')
                         outString = trim(outString) // trim(m)
                     case ('s')
                         outString = trim(outString) // trim(s)
                     case default
                        call mpas_dmpar_global_abort('ERROR: mpas_timekeeping')
                  end select
                  curLen = len_trim(outString)
                  charExpand = .false.
                  write(message,*) 'outString: ',trim(outString)
                  call fckit_log%debug(message)
               else
                  outString(curLen+1:curLen+1) = outString2(i:i)
                  curLen = curLen+1
               end if
           end if
     end do

     outString2 = trim(outString)
     write(message,*) 'cvt_oopsmpas_date output ',trim(outstring2)
     call fckit_log%debug(message)

  end subroutine cvt_oopsmpas_date



! ------------------------------------------------------------------------
!  chunk of code from DART
! ------------------------------------------------------------------------

subroutine uv_cell_to_edges(domain, u_field, v_field, du, lonCell, latCell, &
                            &  nCells, edgeNormalVectors, nEdgesOnCell, edgesOnCell, nVertLevels)

   ! Project u, v wind increments at cell centers onto the edges.
   ! FIXME:
   !        we can hard-code R3 here since it comes from the (3d) x/y/z cartesian coordinate.
   !        We define nEdgesOnCell in get_grid_dims, and read edgesOnCell in get_grid.
   !        We read edgeNormalVectors in get_grid to use this subroutine.
   !        Here "U" is the prognostic variable in MPAS, and we update it with the wind
   !        increments at cell centers.

   implicit none

   type (domain_type), pointer, intent(inout) :: domain
   type (field2DReal), pointer, intent(in) :: u_field    ! u wind updated from filter
   type (field2DReal), pointer, intent(in) :: v_field    ! v wind updated from filter
   type (field2DReal), pointer, intent(inout) :: du       ! normal velocity increment on the edges
   real(RKIND), intent(in) :: lonCell(1:nCells), latCell(1:nCells) ! lon, lat at cell centers in radians
   real(RKIND), intent(in) :: edgeNormalVectors(:,:)
   integer, intent(in) :: nEdgesOnCell(:), edgesOnCell(:,:)
   integer, intent(in) :: nCells, nVertLevels

   ! Local variables
   integer, parameter :: R3 = 3
   real(RKIND), dimension(:,:), allocatable :: east, north
   integer  :: iCell, iEdge, jEdge, k

   ! allocation
   allocate(east(R3,nCells))
   allocate(north(R3,nCells))

   ! Initialization
   du%array(:,:) = MPAS_JEDI_ZERO_kr

   ! Compute unit vectors in east and north directions for each cell:
   do iCell = 1, nCells
       east(1,iCell) = -sin(lonCell(iCell))
       east(2,iCell) =  cos(lonCell(iCell))
       east(3,iCell) =  MPAS_JEDI_ZERO_kr
       call r3_normalize(east(1,iCell), east(2,iCell), east(3,iCell))
       north(1,iCell) = -cos(lonCell(iCell))*sin(latCell(iCell))
       north(2,iCell) = -sin(lonCell(iCell))*sin(latCell(iCell))
       north(3,iCell) =  cos(latCell(iCell))
       call r3_normalize(north(1,iCell), north(2,iCell), north(3,iCell))
   end do

   ! Project analysis increments from the cell centers to the edges

   do iCell = 1, nCells
      do jEdge = 1, nEdgesOnCell(iCell)
         iEdge = edgesOnCell(jEdge, iCell)
            do k = 1, nVertLevels
               du%array(k,iEdge) = du%array(k,iEdge) + MPAS_JEDI_HALF_kr * u_field%array(k,iCell)   &
                     * (edgeNormalVectors(1,iEdge) * east(1,iCell)  &
                     +  edgeNormalVectors(2,iEdge) * east(2,iCell)  &
                     +  edgeNormalVectors(3,iEdge) * east(3,iCell)) &
                     + MPAS_JEDI_HALF_kr * v_field%array(k,iCell)            &
                     * (edgeNormalVectors(1,iEdge) * north(1,iCell) &
                     +  edgeNormalVectors(2,iEdge) * north(2,iCell) &
                     +  edgeNormalVectors(3,iEdge) * north(3,iCell))
            end do
      end do
   end do

   ! deallocation
   deallocate(east)
   deallocate(north)

end subroutine uv_cell_to_edges

!----------------------------------------------------------------------

subroutine r3_normalize(ax, ay, az)

   implicit none

   real(RKIND), intent(inout) :: ax, ay, az
   real(RKIND) :: mi

   mi = MPAS_JEDI_ONE_kr / sqrt(ax**2 + ay**2 + az**2)

   ax = ax * mi
   ay = ay * mi
   az = az * mi

end subroutine r3_normalize

!===============================================================================================================

end module mpas4da_mod

