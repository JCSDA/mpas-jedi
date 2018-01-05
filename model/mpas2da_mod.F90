module mpas2da_mod

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_abort, only : mpas_dmpar_global_abort
   !use random_vectors_mod
 
   type (MPAS_Clock_type), pointer :: clock

   contains

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
      type (mpas_pool_type), pointer, intent(inout) :: pool_b
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
                     else if (poolItr % nDims == 2) then
                        call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field2d)
                        call mpas_pool_add_field(pool_b, trim(poolItr % memberName), field2d)
                     else if (poolItr % nDims == 3) then
                        call mpas_pool_get_field(pool_a, trim(poolItr % memberName), field3d)
                        call mpas_pool_add_field(pool_b, trim(poolItr % memberName), field3d)
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
      character (len=StrKIND) :: kind_op

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
                  if (present(pool_c)) then
                     call mpas_pool_get_array(pool_c, trim(poolItr % memberName), r2d_ptr_c)
                     r2d_ptr_a = 0.
                  end if
                  if ( trim(kind_op).eq.'add' ) then
                     if (present(pool_c)) then
                        r2d_ptr_a = r2d_ptr_b + r2d_ptr_c
                     else
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



end module 

