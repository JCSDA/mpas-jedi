! Copyright (c) 2018, niversity Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html

module mpas2da_mod

   !***********************************************************************
   !
   !  Module mpas2da_mod to encapsulate operations needed for
   !  Data assimilation purpose.
   !  It can be used from /somewhere/MPAS/src/operators 
   !  or from /somewhere/mpas-bundle/mpas/model (OOPS) 
   !> \author  Gael Descombes/Mickael Duda NCAR/MMM
   !> \date    January 2018
   !
   !-----------------------------------------------------------------------

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_abort, only : mpas_dmpar_global_abort
   !use random_vectors_mod
 

   contains

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
                  write(0,*)'Operator add MIN/MAX: ',minval(r2d_ptr_a),maxval(r2d_ptr_a)
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
   !  subroutine da_make_subpool
   !
   !> \brief   Subset a pool from pools A to B
   !> \author  Gael Descombes
   !> \date    26 December 2017
   !> \details
   !>  Given pool A, create pool B as a subset of the fields in A
   !
   !-----------------------------------------------------------------------
   subroutine da_make_subpool(pool_a, pool_b, fieldname)

      implicit none

      type (mpas_pool_type), pointer, intent(in) :: pool_a
      type (mpas_pool_type), pointer, intent(out) :: pool_b
      character (len=*), intent(in) :: fieldname(:)

      type (mpas_pool_iterator_type) :: poolItr
      type (field0DReal), pointer :: field0d
      type (field1DReal), pointer :: field1d
      type (field2DReal), pointer :: field2d
      type (field3DReal), pointer :: field3d
      integer :: nsize, ii, jj

      jj = 0
      nsize = da_common_vars(pool_a, fieldname)
      !write(0,*)'--Create a sub Pool from list of variable'
      !call mpas_pool_create_pool(pool_b)
      call mpas_pool_create_pool(pool_b, nsize)
      write(0,*)'Fieldname: ',nsize,fieldname(:)
      !
      ! Iterate over all fields in pool_b, adding them to fields of the same
      ! name in pool_a
      call mpas_pool_begin_iteration(pool_a)

      do ii=1, nsize
         do while ( mpas_pool_get_next_member(pool_a, poolItr) )
            ! Pools may in general contain dimensions, namelist options, fields, or other pools,
            ! so we select only those members of the pool that are fields
            if (poolItr % memberType == MPAS_POOL_FIELD) then
               ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
               if (poolItr % dataType == MPAS_POOL_REAL) then
                  if ( trim(fieldname(ii)).eq.(trim(poolItr % memberName)) ) then
                     write(0,*)'Adding field in the pool da_make_subpool: '//trim(fieldname(ii))
                     ! Depending on the dimensionality of the field, we need to set pointers of
                     ! the correct type
                     if (poolItr % nDims == 0) then
                        call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field0d)
                        call mpas_pool_add_field(pool_b, trim(poolItr % memberName), field0d)
                     else if (poolItr % nDims == 1) then
                        call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field1d)
                        call mpas_pool_add_field(pool_b, trim(poolItr % memberName), field1d)
                        write(0,*) '1D MIN/MAX value: ', maxval(field1d % array),minval(field1d % array)
                     else if (poolItr % nDims == 2) then
                        call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field2d)
                        call mpas_pool_add_field(pool_b, trim(poolItr % memberName), field2d)
                        write(0,*) '2D MIN/MAX value: ', maxval(field2d % array),minval(field2d % array)
                     else if (poolItr % nDims == 3) then
                        call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field3d)
                        call mpas_pool_add_field(pool_b, trim(poolItr % memberName), field3d)
                        write(0,*) '3D MIN/MAX value: ', maxval(field3d % array),minval(field3d % array)
                     end if
                     jj = jj + 1
                  end if
               end if
            end if
         end do
         if ( ii.ne.jj ) then
            write(0,*)'Missing field in the pool da_make_subpool: '//trim(fieldname(ii))
            ! call mpas_dmpar_global_abort('Missing field in the pool da_make_subpool: '//trim(fieldname(ii))) 
         end if
      end do

      write(0,*)'da_make_subpool new Pool size ',Pool_b % size

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
      nsize = size(fieldname)
      call mpas_pool_begin_iteration(pool_a)

      do ii=1, nsize
         do while ( mpas_pool_get_next_member(pool_a, poolItr) )
            ! Pools may in general contain dimensions, namelist options, fields, or other pools,
            ! so we select only those members of the pool that are fields
            if (poolItr % memberType == MPAS_POOL_FIELD) then
               ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
               if (poolItr % dataType == MPAS_POOL_REAL) then
                  if ( trim(fieldname(ii)).eq.(trim(poolItr % memberName)) ) then
                     write(0,*)'Common field: '//trim(fieldname(ii))
                     nsize0 = nsize0 + 1
                  end if
               end if
            end if
         end do
      end do

      write(0,*)'common_vars = ',nsize0

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
   subroutine da_random(pool_a)

      implicit none

      type (mpas_pool_type), intent(inout),pointer :: pool_a

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
                  !call random_vector(r0d_ptr_a)
               else if (poolItr % nDims == 1) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
                  !call random_vector(r1d_ptr_a)
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_a)
                  !call random_vector(r2d_ptr_a)
               else if (poolItr % nDims == 3) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r3d_ptr_a)
                  !call random_vector(r3d_ptr_a)
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
   subroutine da_operator(kind_op, pool_a, pool_b, pool_c)

      implicit none

      type (mpas_pool_type), pointer :: pool_a, pool_b
      type (mpas_pool_type), pointer, optional :: pool_c
      character (len=*) :: kind_op

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

      !write(0,*)'-------------------------------------------------'
      !write(0,*)' Operator ',trim(kind_op)
      !write(0,*)'-------------------------------------------------'

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
                  if (present(pool_c)) then
                     call mpas_pool_get_array(pool_c, trim(poolItr % memberName), r0d_ptr_c)
                     r0d_ptr_a = 0.
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
                     r1d_ptr_a = 0.
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
                  write(0,*)'Operator_a add MIN/MAX: ',minval(r2d_ptr_a),maxval(r2d_ptr_a) 
                  write(0,*)'Operator_b add MIN/MAX: ',minval(r2d_ptr_b),maxval(r2d_ptr_b) 
                  if (present(pool_c)) then
                     call mpas_pool_get_array(pool_c, trim(poolItr % memberName), r2d_ptr_c)
                     r2d_ptr_a = 0.
                  end if
                  if ( trim(kind_op).eq.'add' ) then
                     if (present(pool_c)) then
                        r2d_ptr_a = r2d_ptr_b + r2d_ptr_c
                     else
                        write(*,*)'regular addition'
                        r2d_ptr_a = r2d_ptr_a + r2d_ptr_b
                     end if
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
                  write(0,*)'Operator2 add MIN/MAX: ',minval(r2d_ptr_a),maxval(r2d_ptr_a) 

               else if (poolItr % nDims == 3) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r3d_ptr_a)
                  call mpas_pool_get_array(pool_b, trim(poolItr % memberName), r3d_ptr_b)
                  if (present(pool_c)) then
                     call mpas_pool_get_array(pool_c, trim(poolItr % memberName), r3d_ptr_c)
                     r3d_ptr_a = 0.
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
   !  subroutine da_zeros
   !
   !> \brief   Performs A = 0. for pool A
   !> \author  Gael Descombes
   !> \date    22 December 2017
   !> \details
   !
   !-----------------------------------------------------------------------
   subroutine da_zeros(pool_a)

      implicit none

      type (mpas_pool_type), pointer :: pool_a

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
                  r0d_ptr_a = 0.
               else if (poolItr % nDims == 1) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
                  r1d_ptr_a = 0.
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_a)
                  r2d_ptr_a = 0.
               else if (poolItr % nDims == 3) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r3d_ptr_a)
                  r3d_ptr_a = 0.
               end if

            end if
         end if
      end do

   end subroutine da_zeros

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
   subroutine da_axpy(pool_a, pool_b, zz)

      implicit none

      type (mpas_pool_type), pointer :: pool_a, pool_b
      real (kind=RKIND) :: zz

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
   !> \brief   Performs basic statistics min/max/avg given a pool
   !> \author  Gael Descombes
   !> \date    January 2018
   !> \details
   !>  Given a pool of fields, return min/max/avg array
   !
   !-----------------------------------------------------------------------
   
   subroutine da_gpnorm(pool_a, nf, pstat)

   implicit none
   type (mpas_pool_type), intent(in),  pointer :: pool_a
   integer,              intent(in) :: nf
   real(kind=RKIND), intent(inout)  :: pstat(3, nf)

   type (mpas_pool_iterator_type) :: poolItr
   type (field1DReal), pointer :: r1d_ptr_a
   type (field2DReal), pointer :: r2d_ptr_a
   type (field3DReal), pointer :: r3d_ptr_a
   !real (field0DReal), pointer :: r0d_ptr_a
   !real (kind=RKIND), pointer :: r0d_ptr_a
   !real (kind=RKIND), dimension(:), pointer :: r1d_ptr_a
   !real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a
   !real (kind=RKIND), dimension(:,:,:), pointer :: r3d_ptr_a

   integer :: jj

   pstat = 0.

   !
   ! Iterate over all fields in pool_a
   ! name in pool_a
   !
   call mpas_pool_begin_iteration(pool_a)
   jj = 1

      do while ( mpas_pool_get_next_member(pool_a, poolItr) )

         ! Pools may in general contain dimensions, namelist options, fields, or other pools,
         ! so we select only those members of the pool that are fields
         if (poolItr % memberType == MPAS_POOL_FIELD) then

            ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
            if (poolItr % dataType == MPAS_POOL_REAL) then

               ! Depending on the dimensionality of the field, we need to set pointers of
               ! the correct type
               !if (poolItr % nDims == 0) then
               !   call mpas_pool_get_field(pool_a, trim(poolItr % memberName), r0d_ptr_a)
               !else if (poolItr % nDims == 1) then
               if (poolItr % nDims == 1) then
                  call mpas_pool_get_field(pool_a, trim(poolItr % memberName), r1d_ptr_a)
                  !pstat(1,jj)=minval(r1d_ptr_a%array)
                  !pstat(2,jj)=maxval(r1d_ptr_a%array)
                  !pstat(3,jj)=sqrt(sum(fld%fld(:,:,jj)**2) && /real(nl*nC,kind_real))
                  !deallocate(r1d_ptr_a)
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_field(pool_a, trim(poolItr % memberName), r2d_ptr_a)
                  pstat(1,jj)= minval(r2d_ptr_a%array)
                  pstat(2,jj)= maxval(r2d_ptr_a%array)
                  !!pstat(3,jj)=sqrt(sum(fld%fld(:,:,jj)**2) && /real(nl*nC,kind_real))
                  !deallocate(r2d_ptr_a) 
               else if (poolItr % nDims == 3) then
                  call mpas_pool_get_field(pool_a, trim(poolItr % memberName), r3d_ptr_a)
                  !pstat(1,jj)=minval(r3d_ptr_a)
                  !pstat(2,jj)=maxval(r3d_ptr_a)
                  !deallocate(r3d_ptr_a) 
                  !pstat(3,jj)=sqrt(sum(fld%fld(:,:,jj)**2) && /real(nl*nC,kind_real))
               end if
               write(0,*)'Variable: ',trim(poolItr % memberName),jj
               write(0,*)'Min/Max stat: ',pstat(1,jj),pstat(2,jj)
               write(0,*)'' 
            end if
         end if
         jj = jj + 1

      end do

   end subroutine da_gpnorm

end module 

