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

!oops
use kinds, only : kind_real

!saber?
use random_mod, only: normal_distribution

!ufo
use ufo_vars_mod

!MPAS-Model
use mpas_abort, only : mpas_dmpar_global_abort
use mpas_constants
use mpas_derived_types
use mpas_dmpar
use mpas_field_routines
use mpas_kind_types, only: ShortStrKIND
use mpas_pool_routines

!mpas-jedi
use mpas_constants_mod

   contains

   !***********************************************************************
   !
   !  function match_scalar_and_index_q
   !
   !> \brief   Test for a form of water vapor or hydrometeor of interest
   !> \author  Steven Vahl
   !> \date    11 July 2019
   !> \details
   !>  At various places in this module we wish to test for the case where
   !>  one string is 'scalars' and another string is one of several
   !>  'index_q?' values. Rather than repeat that logic many times, it
   !>  is encapsulated here.
   !
   !-----------------------------------------------------------------------
   pure function match_scalar_and_index_q(scalarName, indexName)

      implicit none

      character (len=*), intent(in) :: scalarName
      character (len=*), intent(in) :: indexName
      logical :: match_scalar_and_index_q

      match_scalar_and_index_q = scalarName.eq.'scalars' .and. &
                     (indexName.eq.'index_qv' .or. & ! water vapor
                      indexName.eq.'index_qc' .or. & ! cloud
                      indexName.eq.'index_qi' .or. & ! ice
                      indexName.eq.'index_qr' .or. & ! rain
                      indexName.eq.'index_qs' .or. & ! snow
                      indexName.eq.'index_qg' .or. & ! graupel
                      indexName.eq.'index_qh' .or. & ! hail
                      indexName.eq.'index_nr'      & ! number concentration of rain
                     )

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
      field % array(:,:) = 1.0
      write(0,*)'Dimensions Field: ',field % dimSizes(:)

      call mpas_pool_get_field(allFields, 'rho', field)
      call mpas_pool_add_field(da_state, 'rho', field)
      write(0,*) 'Now, max value of rho is ', maxval(field % array),minval(field % array)
      field % array(:,:) = 1.0

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
      real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
      real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
      real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
      real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b

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
!                  write(0,*)'Operator add MIN/MAX: ',minval(r2d_ptr_a),maxval(r2d_ptr_a)
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
   subroutine da_copy_all2sub_fields(domain, pool_a)

      implicit none

      type (mpas_pool_type), pointer, intent(inout) :: pool_a
      type (mpas_pool_type), pointer :: pool_b, state
      type (domain_type), pointer, intent(in) :: domain

      type (mpas_pool_iterator_type) :: poolItr_a, poolItr_b
      real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
      real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
      real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
      integer, pointer :: index_scalar

      type (field2DReal), pointer :: field2d
      type (field3DReal), pointer :: field3d


      pool_b => domain % blocklist % allFields
      call mpas_pool_get_subpool(domain % blocklist % structs,'state',state)
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

               if (( trim(poolItr_b % memberName)).eq.(trim(poolItr_a % memberName)) ) then

                  ! Depending on the dimensionality of the field, we need to set pointers of
                  ! the correct type
                  if (poolItr_b % nDims == 0) then
                     call mpas_pool_get_array(pool_a, trim(poolItr_a % memberName), r0d_ptr_a)
                     call mpas_pool_get_array(pool_b, trim(poolItr_b % memberName), r0d_ptr_b)
                     r0d_ptr_a = r0d_ptr_b
                  else if (poolItr_b % nDims == 1) then
                     call mpas_pool_get_array(pool_a, trim(poolItr_a % memberName), r1d_ptr_a)
                     call mpas_pool_get_array(pool_b, trim(poolItr_b % memberName), r1d_ptr_b)
                     r1d_ptr_a = r1d_ptr_b
                  else if (poolItr_b % nDims == 2) then
!                     write(0,*) 'tmp poolItr_b % memberName=',trim(poolItr_b % memberName)
                     call mpas_pool_get_array(pool_a, trim(poolItr_a % memberName), r2d_ptr_a)
                     call mpas_pool_get_array(pool_b, trim(poolItr_b % memberName), r2d_ptr_b)
                     r2d_ptr_a = r2d_ptr_b
!                     write(0,*)'Copy all2sub field ',trim(poolItr_b % memberName),' MIN/MAX: ',minval(r2d_ptr_a),maxval(r2d_ptr_a)
                  end if

               else if ( match_scalar_and_index_q(trim(poolItr_b % memberName), trim(poolItr_a % memberName)) ) then
!                 write(0,*)'Copy all2sub field: Looking at SCALARS now'
                  call mpas_pool_get_dimension(state, trim(poolItr_a % memberName), index_scalar)
                  if (index_scalar .gt. 0) then
                     call mpas_pool_get_field(pool_a, trim(poolItr_a % memberName), field2d)
                     call mpas_pool_get_field(pool_b, trim(poolItr_b % memberName), field3d)
                     field2d % array(:,:) = field3d % array(index_scalar,:,:)
!                        write(0,*)'Copy all2sub field ',trim(poolItr_a % memberName), &
!                                  minval(field2d % array), maxval(field2d % array)
                  else
                     write(0,*)'WARNING in Copy all2sub field; ',trim(poolItr_a % memberName), &
                                 'not available from MPAS'
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
   subroutine da_copy_sub2all_fields(domain, pool_a)

      implicit none

      type (mpas_pool_type), pointer, intent(in) :: pool_a
      type (mpas_pool_type), pointer :: pool_b, state
      type (domain_type), pointer, intent(inout) :: domain

      type (mpas_pool_iterator_type) :: poolItr_a, poolItr_b
      real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
      real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
      real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
      integer, pointer :: index_scalar

      type (field2DReal), pointer :: field2d
      type (field3DReal), pointer :: field3d


      pool_b => domain % blocklist % allFields
      call mpas_pool_get_subpool(domain % blocklist % structs,'state',state)
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

               if (( trim(poolItr_b % memberName)).eq.(trim(poolItr_a % memberName)) ) then

                  ! Depending on the dimensionality of the field, we need to set pointers of
                  ! the correct type
                  if (poolItr_b % nDims == 0) then
                     call mpas_pool_get_array(pool_a, trim(poolItr_a % memberName), r0d_ptr_a)
                     call mpas_pool_get_array(pool_b, trim(poolItr_b % memberName), r0d_ptr_b)
                     r0d_ptr_b = r0d_ptr_a
                  else if (poolItr_b % nDims == 1) then
                     call mpas_pool_get_array(pool_a, trim(poolItr_a % memberName), r1d_ptr_a)
                     call mpas_pool_get_array(pool_b, trim(poolItr_b % memberName), r1d_ptr_b)
                     r1d_ptr_b = r1d_ptr_a
!                     write(0,*)'Copy sub2all field MIN/MAX: ',trim(poolItr_b % memberName),minval(r1d_ptr_a),maxval(r1d_ptr_a)
                  else if (poolItr_b % nDims == 2) then
                     call mpas_pool_get_array(pool_a, trim(poolItr_a % memberName), r2d_ptr_a)
                     call mpas_pool_get_array(pool_b, trim(poolItr_b % memberName), r2d_ptr_b)
                     r2d_ptr_b = r2d_ptr_a
!                     write(0,*)'Copy sub2all field MIN/MAX: ',trim(poolItr_b % memberName),minval(r2d_ptr_a),maxval(r2d_ptr_a)
                  end if

               else if ( match_scalar_and_index_q(trim(poolItr_b % memberName), trim(poolItr_a % memberName)) ) then
                  !write(0,*)'Copy sub2all field: Looking at SCALARS now',trim(poolItr_a % memberName)
                  call mpas_pool_get_dimension(state, trim(poolItr_a % memberName), index_scalar)
                  if (index_scalar .gt. 0) then
                     call mpas_pool_get_field(pool_a, trim(poolItr_a % memberName), field2d)
                     call mpas_pool_get_field(pool_b, trim(poolItr_b % memberName), field3d)
                     field3d % array(index_scalar,:,:) = field2d % array(:,:)
!                             write(0,*)'Copy sub2all field ',trim(poolItr_a % memberName), &
!                                       minval(field2d % array), maxval(field2d % array)
                  else
                     write(0,*)'WARNING in Copy sub2all field; ',trim(poolItr_a % memberName), &
                              'not available from MPAS'
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
   !  subroutine da_make_subpool
   !
   !> \brief   Subset a pool from pools A to B
   !> \author  Gael Descombes
   !> \date    26 December 2017
   !> \details
   !>  Given pool A, create pool B as a subset of the fields in A
   !
   !-----------------------------------------------------------------------
   subroutine da_make_subpool(domain, subFields, nsize, fieldname, nfields)

      implicit none

      type (domain_type), pointer, intent(in) :: domain
      type (mpas_pool_type), pointer, intent(out) :: subFields
      character (len=*), intent(in) :: fieldname(:)
      integer, intent(out) :: nfields
      integer, intent(in) :: nsize
      type (mpas_pool_type), pointer :: allFields, state

      type (mpas_pool_iterator_type) :: poolItr
      type (field0DReal), pointer :: field0d_src, field0d_dst
      type (field1DReal), pointer :: field1d_src, field1d_dst
      type (field2DReal), pointer :: field2d_src, field2d_dst
      type (field3DReal), pointer :: field3d_src, field3d_dst, field3d
      type (field1DInteger), pointer :: ifield1d_src, ifield1d_dst

      integer, pointer :: index_scalar, dim0d
      integer :: ii
      integer, parameter :: ndims=10
      character(len=ShortStrKIND) :: dimnames(ndims)

      nfields = 0
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state',state) 
      allFields => domain % blocklist % allFields
 
      call mpas_pool_create_pool(subFields)

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
         call mpas_pool_get_dimension(domain % blocklist % dimensions, trim(dimnames(ii)), dim0d)
!         write(0,*)'Adding dimension ',trim(dimnames(ii)),': ',dim0d
         call mpas_pool_add_dimension(subFields, trim(dimnames(ii)), dim0d)
      end do

      !write(0,*)'Fieldname: ',nsize,fieldname(:)
      !
      ! Iterate over allFields and copy those with names matching fieldname into subFields
      !
      call mpas_pool_begin_iteration(allFields)
!      write(0,*)'Before iterating'

      do while ( mpas_pool_get_next_member(allFields, poolItr) )
         ! Pools may in general contain dimensions, namelist options, fields, or other pools,
         ! so we select only those members of the pool that are fields
         if (poolItr % memberType == MPAS_POOL_FIELD) then
            ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
            if (poolItr % dataType == MPAS_POOL_REAL) then
               do ii=1, nsize
                  if ( trim(fieldname(ii)).eq.(trim(poolItr % memberName)) ) then
!                     write(0,*)'Adding field in the pool da_make_subpool: '//trim(fieldname(ii))
                     ! Depending on the dimensionality of the field, we need to set pointers of
                     ! the correct type
                     if (poolItr % nDims == 0) then
                        call mpas_pool_get_field(allFields, trim(poolItr % memberName), field0d_src)
                        call mpas_duplicate_field(field0d_src, field0d_dst)
                        call mpas_pool_add_field(subFields, trim(poolItr % memberName), field0d_dst)
                     else if (poolItr % nDims == 1) then
                        call mpas_pool_get_field(allFields, trim(poolItr % memberName), field1d_src)
                        call mpas_duplicate_field(field1d_src, field1d_dst)
                        call mpas_pool_add_field(subFields, trim(poolItr % memberName), field1d_dst)
                     else if (poolItr % nDims == 2) then
                        call mpas_pool_get_field(allFields, trim(poolItr % memberName), field2d_src)
                        call mpas_duplicate_field(field2d_src, field2d_dst)
                        call mpas_pool_add_field(subFields, trim(poolItr % memberName), field2d_dst)
                     else if (poolItr % nDims == 3) then
                        call mpas_pool_get_field(allFields, trim(poolItr % memberName), field3d_src)
                        call mpas_duplicate_field(field3d_src, field3d_dst)
                        call mpas_pool_add_field(subFields, trim(poolItr % memberName), field3d_dst)
                     end if
                     nfields = nfields + 1
                  
                  else if (match_scalar_and_index_q(trim(poolItr % memberName), trim(fieldname(ii)))) then
                     call mpas_pool_get_dimension(state, trim(fieldname(ii)), index_scalar)
                     if (index_scalar .gt. 0) then
                        call mpas_pool_get_field(allFields, trim(poolItr % memberName), field3d)
                        call mpas_pool_get_field(allFields, 'theta_m', field2d_src)
                        call mpas_duplicate_field(field2d_src, field2d_dst)
                        field2d_dst % fieldName = trim(fieldname(ii))
                        field2d_dst % array(:,:) = field3d % array(index_scalar,:,:)
                        call mpas_pool_add_field(subFields, trim(fieldname(ii)), field2d_dst)

                        nfields = nfields + 1
                     else
                        write(0,*)'WARNING in da_make_subpool; ',trim(fieldname(ii)), &
                                    'not available from MPAS'
                     end if
                  end if
              end do
            end if
            if (poolItr % dataType == MPAS_POOL_INTEGER) then
              do ii=1, nsize
                  if ( trim(fieldname(ii)).eq.(trim(poolItr % memberName)) ) then
!                     write(0,*)'Adding field in the pool da_make_subpool: '//trim(fieldname(ii))
                     ! Depending on the dimensionality of the field, we need to set pointers of
                     ! the correct type
                     if (poolItr % nDims == 1) then
                        call mpas_pool_get_field(allFields, trim(poolItr % memberName), ifield1d_src)
                        call mpas_duplicate_field(ifield1d_src, ifield1d_dst)
                        call mpas_pool_add_field(subFields, trim(poolItr % memberName), ifield1d_dst)

!                        write(0,*) '1D MIN/MAX value: ', minval(ifield1d_dst % array),maxval(ifield1d_dst % array)
                        nfields = nfields + 1
                     end if
                  end if
              end do
            end if
          end if
      end do

      if ( nsize.ne.nfields ) then
        write(0,*)'Missing field in the pool da_make_subpool nsize different: ',nsize,nfields
      end if

!      write(0,*)'da_make_subpool new Pool size ',subFields % size

   end subroutine da_make_subpool

   
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
!      write(0,*)'da_common_vars size: ',nsize
      call mpas_pool_begin_iteration(pool_a)

         do while ( mpas_pool_get_next_member(pool_a, poolItr) )
            ! Pools may in general contain dimensions, namelist options, fields, or other pools,
            ! so we select only those members of the pool that are fields
            if (poolItr % memberType == MPAS_POOL_FIELD) then
               ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
               if (poolItr % dataType == MPAS_POOL_REAL) then
                  do ii=1, nsize
                     if ( trim(fieldname(ii)).eq.(trim(poolItr % memberName)) ) then
!                        write(0,*)'Common field: '//trim(fieldname(ii))
                        nsize0 = nsize0 + 1
                     else if (match_scalar_and_index_q(trim(poolItr % memberName), trim(fieldname(ii)))) then
!                        write(0,*)'Common field: '//trim(fieldname(ii))
                        nsize0 = nsize0 + 1 
                     end if
                  end do
               end if
            end if
         end do

!      write(0,*)'da_common_vars = ',nsize0

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
      real (kind=kind_real), pointer :: r0d_ptr_a
      real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a
      real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a
      real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a

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
                  !call normal_distribution(r0d_ptr_a, 0.0_kind_real, 1.0_kind_real, rseed)
               else if (poolItr % nDims == 1) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
                  call normal_distribution(r1d_ptr_a, 0.0_kind_real, 1.0_kind_real, rseed)
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_a)
                  call normal_distribution(r2d_ptr_a, 0.0_kind_real, 1.0_kind_real, rseed)
               else if (poolItr % nDims == 3) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r3d_ptr_a)
                  call normal_distribution(r3d_ptr_a, 0.0_kind_real, 1.0_kind_real, rseed)
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
      real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b, r0d_ptr_c
      real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b, r1d_ptr_c
      real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b, r2d_ptr_c
      real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b, r3d_ptr_c
      !
      ! Iterate over all fields in pool_b, adding them to fields of the same
      ! name in pool_a
      !
      call mpas_pool_begin_iteration(pool_b)

      !write(0,*)'-------------------------------------------------'
      !write(0,*)' Operator ',trim(kind_op)
      !write(0,*)'-------------------------------------------------'

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
                     r0d_ptr_a = 0.0_kind_real
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
                     r1d_ptr_a = 0.0_kind_real
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
                     r2d_ptr_a = 0.0_kind_real
                  end if
                  if ( trim(kind_op).eq.'add' ) then
!                     write(0,*)'Operator_a add MIN/MAX: ',minval(r2d_ptr_a),maxval(r2d_ptr_a) 
!                     write(0,*)'Operator_b add MIN/MAX: ',minval(r2d_ptr_b),maxval(r2d_ptr_b) 
                     if (present(pool_c)) then
                        r2d_ptr_a = r2d_ptr_b + r2d_ptr_c
                     else
!                        write(*,*)'regular addition'
                        r2d_ptr_a = r2d_ptr_a + r2d_ptr_b
                     end if
!                     write(0,*)'Operator2 add MIN/MAX: ',minval(r2d_ptr_a),maxval(r2d_ptr_a) 
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
                     r3d_ptr_a = 0.0_kind_real
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
      real (kind=kind_real) :: zz

      type (mpas_pool_iterator_type) :: poolItr
      real (kind=kind_real), pointer :: r0d_ptr_a
      real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a
      real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a
      real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a

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
   !  subroutine da_zeros
   !
   !> \brief   Performs A = 0. for pool A
   !> \author  Gael Descombes
   !> \date    22 December 2017
   !> \details
   !
   !-----------------------------------------------------------------------
   subroutine da_zeros(pool_a, fld_select)

      implicit none

      type (mpas_pool_type), pointer, intent(inout) :: pool_a
      character (len=*), optional,    intent(in)    :: fld_select(:)

      type (mpas_pool_iterator_type) :: poolItr
      real (kind=kind_real), pointer :: r0d_ptr_a
      real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a
      real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a
      real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a

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
                  r0d_ptr_a = 0.0_kind_real
               else if (poolItr % nDims == 1) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
                  r1d_ptr_a = 0.0_kind_real
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_a)
                  r2d_ptr_a = 0.0_kind_real
               else if (poolItr % nDims == 3) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r3d_ptr_a)
                  r3d_ptr_a = 0.0_kind_real
               end if

            end if
         end if
      end do

   end subroutine da_zeros

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
      real (kind=kind_real), pointer :: r0d_ptr_a
      real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a
      real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a
      real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a

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
                  r0d_ptr_a = max(0.0_kind_real, r0d_ptr_a)
               else if (poolItr % nDims == 1) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
                  r1d_ptr_a = max(0.0_kind_real, r1d_ptr_a)
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_a)
                  r2d_ptr_a = max(0.0_kind_real, r2d_ptr_a)
               else if (poolItr % nDims == 3) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r3d_ptr_a)
                  r3d_ptr_a = max(0.0_kind_real, r3d_ptr_a)
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
      real (kind=kind_real), pointer :: r0d_ptr_a
      real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a
      real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a
      real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a

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
      real (kind=kind_real), intent(in) :: zz
      character (len=*), optional, intent(in) :: fld_select(:)



      type (mpas_pool_iterator_type) :: poolItr
      real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
      real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
      real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
      real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b

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


   !*******************************************************************************
   !
   ! We can have information about which dims are decomposed in a single function
   !
   !*******************************************************************************

   logical function isDecomposed(dimName)

    implicit none

    character(len=*), intent(in) :: dimName

    if (trim(dimName) == 'nCells') then
       isDecomposed = .true.
       return
    else if (trim(dimName) == 'nEdges') then
       isDecomposed = .true.
       return
    else if (trim(dimName) == 'nVertices') then
       isDecomposed = .true.
       return
    else
       isDecomposed = .false.
       return
    end if
!    if ( isDecomposed ) then 
!       write(*,*)'is Decomposed ',trim(dimName)
!    else
!       write(*,*)'is not Decomposed ',trim(dimName)
!    end if

    end function isDecomposed
   

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
   real(kind=kind_real),           intent(out) :: pstat(3, nf)

   type (mpas_pool_iterator_type) :: poolItr
   type (field1DReal), pointer :: field1d
   type (field2DReal), pointer :: field2d
   type (field3DReal), pointer :: field3d
   real(kind=kind_real) :: globalSum, globalMin, globalMax, dimtot, dimtot_global, prodtot

   integer :: jj, ndims
   integer, pointer :: solveDim1, solveDim2, solveDim3
   !integer, pointer :: solveDim(:)

   pstat = 0.0_kind_real

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

            ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
            if (poolItr % dataType == MPAS_POOL_REAL) then

               ! Depending on the dimensionality of the field, we need to set pointers of
               ! the correct type
               ndims = poolItr % nDims

               !write(*,*)'gpnorm variable: ',trim(poolItr % memberName), ndims

               if (ndims == 1) then
                  call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field1d)
                  if ( isDecomposed(field1d % dimNames(ndims)) ) then
                     call mpas_pool_get_dimension(pool_a, trim(field1d % dimNames(ndims))//'Solve', solveDim1)
                  else
                     call mpas_pool_get_dimension(pool_a, trim(field1d % dimNames(ndims)), solveDim1)
                  end if
                  dimtot = real(solveDim1,kind_real)
                  prodtot = sum(field1d % array(1:solveDim1)**2 )
                  call mpas_dmpar_sum_real(dminfo, dimtot, dimtot_global)
                  call mpas_dmpar_sum_real(dminfo, prodtot, globalSum)
                  call mpas_dmpar_min_real(dminfo, minval(field1d % array(1:solveDim1)), globalMin)
                  call mpas_dmpar_max_real(dminfo, maxval(field1d % array(1:solveDim1)), globalMax)
                  pstat(1,jj) = globalMin
                  pstat(2,jj) = globalMax
                  pstat(3,jj) = sqrt( globalSum / dimtot_global )

               else if (ndims == 2) then
                  call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field2d)
                  if (isDecomposed(field2d % dimNames(ndims))) then
                      call mpas_pool_get_dimension(pool_a, trim(field2d % dimNames(ndims))//'Solve', solveDim1)
                  else
                      call mpas_pool_get_dimension(pool_a, trim(field2d % dimNames(ndims)), solveDim1)
                  end if
                  if (isDecomposed(field2d % dimNames(ndims-1))) then
                      call mpas_pool_get_dimension(pool_a, trim(field2d % dimNames(ndims-1))//'Solve', solveDim2)
                  else
                      call mpas_pool_get_dimension(pool_a, trim(field2d % dimNames(ndims-1)), solveDim2)
                  end if
                  dimtot  = real(solveDim1*solveDim2,kind_real)
                  prodtot = sum(field2d % array(1:solveDim2,1:solveDim1)**2 )
                  call mpas_dmpar_sum_real(dminfo, dimtot, dimtot_global)
                  !BJJ ?? call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field2d)
                  call mpas_dmpar_sum_real(dminfo, prodtot, globalSum)
                  call mpas_dmpar_min_real(dminfo, minval(field2d % array(1:solveDim2,1:solveDim1)), globalMin)
                  call mpas_dmpar_max_real(dminfo, maxval(field2d % array(1:solveDim2,1:solveDim1)), globalMax)
!                  write(*,*)'gpnorm prodtot: ',prodtot, dimtot
                  pstat(1,jj) = globalMin
                  pstat(2,jj) = globalMax
                  pstat(3,jj) = sqrt( globalSum / dimtot_global )

               else if (ndims == 3) then
                  call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field3d)
                  if (isDecomposed(field3d % dimNames(ndims))) then
                     call mpas_pool_get_dimension(pool_a, trim(field3d % dimNames(ndims))//'Solve', solveDim1)
                  else
                     call mpas_pool_get_dimension(pool_a, trim(field3d % dimNames(ndims)), solveDim1)
                  end if
                  if (isDecomposed(field3d % dimNames(ndims-1))) then
                     call mpas_pool_get_dimension(pool_a, trim(field3d % dimNames(ndims-1))//'Solve', solveDim2)
                  else
                     call mpas_pool_get_dimension(pool_a, trim(field3d % dimNames(ndims-1)), solveDim2)
                  end if
                  if (isDecomposed(field3d % dimNames(ndims-2))) then
                     call mpas_pool_get_dimension(pool_a, trim(field3d % dimNames(ndims-2))//'Solve', solveDim3)
                  else
                     call mpas_pool_get_dimension(pool_a, trim(field3d % dimNames(ndims-2)), solveDim3)
                  end if
                  dimtot  = real(solveDim1*solveDim2*solveDim3,kind_real)
                  prodtot = sum(field3d % array(1:solveDim3,1:solveDim2,1:solveDim1)**2 )
                  call mpas_dmpar_sum_real(dminfo, dimtot, dimtot_global)
                  !BJJ ?? call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field3d)
                  call mpas_dmpar_sum_real(dminfo, prodtot, globalSum)
                  call mpas_dmpar_min_real(dminfo, minval(field3d % array(1:solveDim3,1:solveDim2,1:solveDim1)), globalMin)
                  call mpas_dmpar_max_real(dminfo, maxval(field3d % array(1:solveDim3,1:solveDim2,1:solveDim1)), globalMax)
                  pstat(1,jj) = globalMin
                  pstat(2,jj) = globalMax
                  pstat(3,jj) = sqrt( globalSum / dimtot_global )
               end if

!               write(0,*)'Variable: ',trim(poolItr % memberName),jj
!               write(0,*)'Min/Max stat: ',pstat(1,jj),pstat(2,jj),pstat(3,jj)
!               write(0,*)'' 

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
   real(kind=kind_real),           intent(out) :: fldrms
   character (len=*), optional,    intent(in)  :: fld_select(:)

   type (mpas_pool_iterator_type) :: poolItr
   type (field1DReal), pointer :: field1d
   type (field2DReal), pointer :: field2d
   type (field3DReal), pointer :: field3d
   real(kind=kind_real) :: dimtot, dimtot_global, prodtot, prodtot_global 

   integer :: jj, ndims
   integer, pointer :: solveDim1, solveDim2, solveDim3

   prodtot = 0.0_kind_real
   dimtot  = 0.0_kind_real

   !
   ! Iterate over all fields in pool_a
   ! named in pool_a
   !
   call mpas_pool_begin_iteration(pool_a)
   jj = 1

   do while ( mpas_pool_get_next_member(pool_a, poolItr) )
         if (present(fld_select)) then
            if (ufo_vars_getindex(fld_select,trim(poolItr % memberName)) < 0) cycle
         end if

         if (poolItr % dataType == MPAS_POOL_REAL) then
            if (poolItr % memberType == MPAS_POOL_FIELD) then
            
               ndims = poolItr % nDims
!               write(*,*)'fldrms variable: ',trim(poolItr % memberName),ndims

               if (ndims == 1) then
                  call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field1d)
                  if (isDecomposed(field1d % dimNames(ndims))) then
                     call mpas_pool_get_dimension(pool_a, trim(field1d % dimNames(ndims))//'Solve', solveDim1)
                  else
                     call mpas_pool_get_dimension(pool_a, trim(field1d % dimNames(ndims)), solveDim1)
                  end if
                  dimtot  = dimtot + real(solveDim1,kind_real)
                  prodtot = prodtot + sum( field1d % array(1:solveDim1)**2 )

               else if (ndims == 2) then
                  call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field2d)
                  if (isDecomposed(field2d % dimNames(ndims))) then
                      call mpas_pool_get_dimension(pool_a, trim(field2d % dimNames(ndims))//'Solve', solveDim1)
                  else
                      call mpas_pool_get_dimension(pool_a, trim(field2d % dimNames(ndims)), solveDim1)
                  end if
                  if (isDecomposed(field2d % dimNames(ndims-1))) then
                      call mpas_pool_get_dimension(pool_a, trim(field2d % dimNames(ndims-1))//'Solve', solveDim2)
                  else
                      call mpas_pool_get_dimension(pool_a, trim(field2d % dimNames(ndims-1)), solveDim2)
                  end if
                  dimtot  = dimtot + real(solveDim1*solveDim2,kind_real)
                  prodtot = prodtot + sum( field2d % array(1:solveDim2,1:solveDim1)**2 )

!                  write(*,*)'fldrms dims: ',solveDim1, solveDim2
!                  write(*,*)'fldrms, dimtot, prodtot: ', dimtot, prodtot

                else if (ndims == 3) then
                  call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field3d)
!                  write(*,*)'dimNames ',trim(field3d % dimNames(ndims))
                  if (isDecomposed(field3d % dimNames(ndims))) then
                     call mpas_pool_get_dimension(pool_a, trim(field3d % dimNames(ndims))//'Solve', solveDim1)
                  else
                     call mpas_pool_get_dimension(pool_a, trim(field3d % dimNames(ndims)), solveDim1)
                  end if
                  if (isDecomposed(field3d % dimNames(ndims-1))) then
                     call mpas_pool_get_dimension(pool_a, trim(field3d % dimNames(ndims-1))//'Solve', solveDim2)
                  else
                     call mpas_pool_get_dimension(pool_a, trim(field3d % dimNames(ndims-1)), solveDim2)
                  end if
                  if (isDecomposed(field3d % dimNames(ndims-2))) then
                     call mpas_pool_get_dimension(pool_a, trim(field3d % dimNames(ndims-2))//'Solve', solveDim3)
                  else
                     call mpas_pool_get_dimension(pool_a, trim(field3d % dimNames(ndims-2)), solveDim3)
                  end if
                  dimtot  = dimtot + real(solveDim1*solveDim2*solveDim3,kind_real)
                  prodtot = prodtot + sum( field3d % array(1:solveDim3,1:solveDim2,1:solveDim1)**2 )

               end if
            end if
         end if
         jj = jj + 1

      end do

      call mpas_dmpar_sum_real(dminfo, dimtot, dimtot_global)
      call mpas_dmpar_sum_real(dminfo, prodtot, prodtot_global)
      fldrms = sqrt(prodtot_global / dimtot_global)
!      write(*,*)'fldrms = sqrt( prodtot_global / dimtot_global) : ', fldrms, prodtot_global, dimtot_global
      

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
   real(kind=kind_real),           intent(out) :: zprod

   type (mpas_pool_iterator_type) :: poolItr
   type (field1DReal), pointer :: field1d_a, field1d_b
   type (field2DReal), pointer :: field2d_a, field2d_b
   type (field3DReal), pointer :: field3d_a, field3d_b
   real(kind=kind_real) :: fieldSum_local, zprod_local

   integer :: jj, ndims
   integer, pointer :: solveDim1, solveDim2, solveDim3

   !
   ! Iterate over all fields in pool_a
   ! named in pool_a
   !
   call mpas_pool_begin_iteration(pool_a)

   zprod_local = 0.0_kind_real

   do while ( mpas_pool_get_next_member(pool_a, poolItr) )

         if (poolItr % dataType == MPAS_POOL_REAL) then
            if (poolItr % memberType == MPAS_POOL_FIELD) then

!               write(*,*)'variable: ',trim(poolItr % memberName)
               ndims = poolItr % nDims 

               if (ndims == 1) then
                  call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field1d_a)
                  if ( isDecomposed(field1d_a % dimNames(ndims)) ) then
                     call mpas_pool_get_dimension(pool_a, trim(field1d_a % dimNames(ndims))//'Solve', solveDim1)
                  else
                     call mpas_pool_get_dimension(pool_a, trim(field1d_a % dimNames(ndims)), solveDim1)
                  end if
                  call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field1d_a)
                  call mpas_pool_get_field(pool_b, trim(poolItr % memberName), field1d_b)
                  fieldSum_local = sum(field1d_a % array(1:solveDim1) * field1d_b % array(1:solveDim1))
                  zprod_local = zprod_local + fieldSum_local
!                  write(*,*)'dotprod: ',trim(field1d_a % dimNames(ndims)), zprod_local

               else if (ndims == 2) then
                  call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field2d_a)
                  if (isDecomposed(field2d_a % dimNames(ndims))) then
                      call mpas_pool_get_dimension(pool_a, trim(field2d_a % dimNames(ndims))//'Solve', solveDim1)
                  else
                      call mpas_pool_get_dimension(pool_a, trim(field2d_a % dimNames(ndims)), solveDim1)
                  end if
                  if (isDecomposed(field2d_a % dimNames(ndims-1))) then
                      call mpas_pool_get_dimension(pool_a, trim(field2d_a % dimNames(ndims-1))//'Solve', solveDim2)
                  else
                      call mpas_pool_get_dimension(pool_a, trim(field2d_a % dimNames(ndims-1)), solveDim2)
                  end if
                  call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field2d_a)
                  call mpas_pool_get_field(pool_b, trim(poolItr % memberName), field2d_b)
                  fieldSum_local = sum(field2d_a % array(1:solveDim2,1:solveDim1) * field2d_b % array(1:solveDim2,1:solveDim1))
                  zprod_local = zprod_local + fieldSum_local
!                  write(*,*)'dotprod: ',trim(field2d_a % dimNames(ndims)), zprod_local

               else if (ndims == 3) then
                  call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field3d_a)
!                  write(*,*)'dimNames ',trim(field3d_a % dimNames(ndims))
                  if (isDecomposed(field3d_a % dimNames(ndims))) then
                     call mpas_pool_get_dimension(pool_a, trim(field3d_a % dimNames(ndims))//'Solve', solveDim1)
                  else
                     call mpas_pool_get_dimension(pool_a, trim(field3d_a % dimNames(ndims)), solveDim1)
                  end if
                  if (isDecomposed(field3d_a % dimNames(ndims-1))) then
                     call mpas_pool_get_dimension(pool_a, trim(field3d_a % dimNames(ndims-1))//'Solve', solveDim2)
                  else
                     call mpas_pool_get_dimension(pool_a, trim(field3d_a % dimNames(ndims-1)), solveDim2)
                  end if
                  if (isDecomposed(field3d_a % dimNames(ndims-2))) then
                     call mpas_pool_get_dimension(pool_a, trim(field3d_a % dimNames(ndims-2))//'Solve', solveDim3)
                  else
                     call mpas_pool_get_dimension(pool_a, trim(field3d_a % dimNames(ndims-2)), solveDim3)
                  end if
                  call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field3d_a)
                  call mpas_pool_get_field(pool_b, trim(poolItr % memberName), field3d_b)
                  fieldSum_local = sum(field3d_a % array(1:solveDim3,1:solveDim2,1:solveDim1) &
                                    * field3d_b % array(1:solveDim3,1:solveDim2,1:solveDim1))
                  zprod_local = zprod_local + fieldSum_local
!                  write(*,*)'dotprod: ',trim(field2d_a % dimNames(ndims)), zprod_local

               end if

            end if
         end if

      end do

      call mpas_dmpar_sum_real(dminfo, zprod_local, zprod)
!      write(*,*)'dotprod: Final result = ',zprod

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

!     write(*,*)'cvt_oopsmpas_date instring: ',trim(YYYY),trim(MM),trim(DD),trim(h),trim(m),trim(s)
!     write(*,*)'cvt_oopsmpas_date input ',trim(instring2)

     write(outString,*) ''
     instring = trim(outstring2)     
     
     curLen = 0
     charExpand = .false.
     do i = 1, len_trim(inString)
           if (inString(i:i) == '$' ) then
               charExpand = .true.
           else if (inString(i:i) /= '$') then
               !write(*,*)'inString: ',trim(inString(i:i)),charExpand
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
                        write(*, *) 'ERROR: mpas_expand_string option $', inString(i:i), ' is not a valid expansion character.'
                        call mpas_dmpar_global_abort('ERROR: mpas_timekeeping')
                  end select
                  curLen = len_trim(outString)
                  charExpand = .false.
                  !write(*,*)'outString: ',trim(outString)
               else
                  outString(curLen+1:curLen+1) = outString2(i:i)
                  curLen = curLen+1
               end if
           end if
     end do

     outString2 = trim(outString)
!     write(*,*)'cvt_oopsmpas_date output ',trim(outstring2)

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
   real(kind_real), intent(in) :: lonCell(1:nCells), latCell(1:nCells)    ! normal velocity increment on the edges
   real(kind_real), intent(in) :: edgeNormalVectors(:,:)
   integer, intent(in) :: nEdgesOnCell(:), edgesOnCell(:,:)
   integer, intent(in) :: nCells, nVertLevels

   ! Local variables
   integer, parameter :: R3 = 3
   real(kind_real), dimension(:,:), allocatable :: east, north
   real(kind_real), dimension(:), allocatable :: lonCell_rad, latCell_rad
   integer  :: iCell, iEdge, jEdge, k

   ! allocation
   allocate(east(R3,nCells))
   allocate(north(R3,nCells))
   allocate(lonCell_rad(nCells))
   allocate(latCell_rad(nCells))

   ! Initialization
   du%array(:,:) = 0.0_kind_real

   ! Back to radians (locally)
   lonCell_rad = lonCell*deg2rad
   latCell_rad = latCell*deg2rad

   ! Compute unit vectors in east and north directions for each cell:
   do iCell = 1, nCells
       east(1,iCell) = -sin(lonCell_rad(iCell))
       east(2,iCell) =  cos(lonCell_rad(iCell))
       east(3,iCell) =  0.0_kind_real
       call r3_normalize(east(1,iCell), east(2,iCell), east(3,iCell))
       north(1,iCell) = -cos(lonCell_rad(iCell))*sin(latCell_rad(iCell))
       north(2,iCell) = -sin(lonCell_rad(iCell))*sin(latCell_rad(iCell))
       north(3,iCell) =  cos(latCell_rad(iCell))
       call r3_normalize(north(1,iCell), north(2,iCell), north(3,iCell))
   end do

   ! Project analysis increments from the cell centers to the edges

   do iCell = 1, nCells
      do jEdge = 1, nEdgesOnCell(iCell)
         iEdge = edgesOnCell(jEdge, iCell)
            do k = 1, nVertLevels
               du%array(k,iEdge) = du%array(k,iEdge) + 0.5_kind_real * u_field%array(k,iCell)   &
                     * (edgeNormalVectors(1,iEdge) * east(1,iCell)  &
                     +  edgeNormalVectors(2,iEdge) * east(2,iCell)  &
                     +  edgeNormalVectors(3,iEdge) * east(3,iCell)) &
                     + 0.5_kind_real * v_field%array(k,iCell)            &
                     * (edgeNormalVectors(1,iEdge) * north(1,iCell) &
                     +  edgeNormalVectors(2,iEdge) * north(2,iCell) &
                     +  edgeNormalVectors(3,iEdge) * north(3,iCell))
            end do
      end do
   end do

   ! deallocation
   deallocate(east)
   deallocate(north)
   deallocate(lonCell_rad)
   deallocate(latCell_rad)

end subroutine uv_cell_to_edges

!----------------------------------------------------------------------

subroutine r3_normalize(ax, ay, az)

   implicit none

   real(kind_real), intent(inout) :: ax, ay, az
   real(kind_real) :: mi

   mi = 1.0_kind_real / sqrt(ax**2 + ay**2 + az**2)

   ax = ax * mi
   ay = ay * mi
   az = az * mi

end subroutine r3_normalize

!===============================================================================================================

end module mpas4da_mod 

