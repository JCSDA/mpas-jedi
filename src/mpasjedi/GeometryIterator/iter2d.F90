! (C) Copyright 2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Geometry iterator
module iter2d_mod

use iso_c_binding
use kinds, only : kind_real

use mpas_derived_types
use mpas_pool_routines
use mpas_kind_types, only: RKIND

use iter_mod
use mpas_constants_mod
use mpas_fields_mod, only: mpas_fields
use mpas_geom_mod, only: getSolveDimSizes

implicit none
private

public :: iter2d

type, extends(iter) :: iter2d
  contains
  procedure :: equals => equals_2d
  procedure :: current => current_2d
  procedure :: next => next_2d
  procedure :: getpoint => getpoint_2d
  procedure :: setpoint => setpoint_2d
end type iter2d

integer, parameter    :: max_string=8000
character(max_string) :: message

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Check for the geometry iterator equality (pointing to same cellIndex)
!!
!! \relates iter2d_mod::iter2d
logical(kind=c_bool) function equals_2d(self, other)
  class(iter2d), intent(in) :: self
  class(iter), pointer, intent(in) :: other
  equals_2d = (self%cellIndex==other%cellIndex)
end function equals_2d

! ------------------------------------------------------------------------------
!> Get geometry iterator current lat/lon/height
!!
!! \throws abor1_ftn aborts if iterator is out of bounds
!! \relates iter2d_mod::iter2d
subroutine current_2d(self, lon, lat, height)
  ! Passed variables
  class(iter2d), intent( in) :: self !< Geometry iterator
  real(kind_real), intent(out) :: lon    !< Longitude
  real(kind_real), intent(out) :: lat    !< Latitude
  real(kind_real), intent(out) :: height !< Height

  integer :: ii

  ! Check cellIndex
  if (self%cellIndex == -1) then
    ! special case of {-1,-1} means end of the grid
    ii = self%geom%nCellsSolve
  elseif (self%cellIndex < 1 .OR. self%cellIndex > self%geom%nCellsSolve) then
    ! outside of the grid
    call abor1_ftn('current_2d: lat/lon iterator out of bounds')
  else
    ! inside of the grid
    ii = self%cellIndex
  endif
  lat = self%geom%latCell(ii) * MPAS_JEDI_RAD2DEG_kr
  lon = self%geom%lonCell(ii) * MPAS_JEDI_RAD2DEG_kr
  height = -99999.
end subroutine current_2d

! ------------------------------------------------------------------------------
!> Update geometry iterator to next point
!!
!! \relates iter2d_mod::iter2d
subroutine next_2d(self)
  class(iter2d), intent(inout) :: self
  if (self%cellIndex.lt.self%geom%nCellsSolve) then
    self%cellIndex = self%cellIndex + 1
  elseif (self%cellIndex.eq.self%geom%nCellsSolve) then
    self%cellIndex=-1
  end if
end subroutine next_2d

! ------------------------------------------------------------------------------

subroutine getpoint_2d(self, fields, values, nval)

   implicit none

   ! Passed variables
   class(iter2d), intent(in) :: self
   class(mpas_fields), intent(in) :: fields
   integer(c_int) :: nval
   real(kind_real), intent(inout) :: values(nval)

   ! Local variables
   integer :: ii, nvert
   type (mpas_pool_iterator_type) :: poolItr
   type(mpas_pool_data_type), pointer :: fdata
   integer, allocatable :: dimSizes(:)

   ! Initialize
   ii = 0

   call mpas_pool_begin_iteration(fields%subFields)
   do while ( mpas_pool_get_next_member(fields%subFields, poolItr) )
      if (poolItr % memberType == MPAS_POOL_FIELD) then
         call fields%get(poolItr % memberName, fdata)
         dimSizes = getSolveDimSizes(fields%subFields, poolItr%memberName)
         if (poolItr % nDims == 1) then
            nvert = 1
            if (poolItr % dataType == MPAS_POOL_INTEGER) then
               values(ii+1) = real(fdata%i1%array(self%cellIndex), kind_real)
            else if (poolItr % dataType == MPAS_POOL_REAL) then
               values(ii+1) = real(fdata%r1%array(self%cellIndex), kind_real)
            endif
         else if (poolItr % nDims == 2) then
            nvert = dimSizes(1)
            if (poolItr % dataType == MPAS_POOL_INTEGER) then
               values(ii+1:ii+nvert) = real(fdata%i2%array(1:nvert, self%cellIndex), kind_real)
            else if (poolItr % dataType == MPAS_POOL_REAL) then
               values(ii+1:ii+nvert) = real(fdata%r2%array(1:nvert, self%cellIndex), kind_real)
            endif
         else
            write(message,*) '--> getpoint_2d: poolItr % nDims == ',poolItr % nDims,' not handled'
            call abor1_ftn(message)
         endif
         ii = ii + nvert
         deallocate(dimSizes)
      endif
   enddo

end subroutine getpoint_2d

! ------------------------------------------------------------------------------

subroutine setpoint_2d(self, fields, values, nval)

   implicit none

   ! Passed variables
   class(iter2d), intent(in) :: self
   class(mpas_fields), intent(inout) :: fields
   integer(c_int) :: nval
   real(kind_real), intent(in) :: values(nval)

   ! Local variables
   integer :: ii, nvert
   type (mpas_pool_iterator_type) :: poolItr
   type(mpas_pool_data_type), pointer :: fdata
   integer, allocatable :: dimSizes(:)

   ! Initialize
   ii = 0

   call mpas_pool_begin_iteration(fields%subFields)
   do while ( mpas_pool_get_next_member(fields%subFields, poolItr) )
      if (poolItr % memberType == MPAS_POOL_FIELD) then
         call fields%get(poolItr % memberName, fdata)
         dimSizes = getSolveDimSizes(fields%subFields, poolItr%memberName)
         if (poolItr % nDims == 1) then
            nvert = 1
            if (poolItr % dataType == MPAS_POOL_INTEGER) then
               fdata%i1%array(self%cellIndex) = nint(values(ii+1))
            else if (poolItr % dataType == MPAS_POOL_REAL) then
               fdata%r1%array(self%cellIndex) = real(values(ii+1), RKIND)
            endif
         else if (poolItr % nDims == 2) then
            nvert = dimSizes(1)
            if (poolItr % dataType == MPAS_POOL_INTEGER) then
               fdata%i2%array(1:nvert, self%cellIndex) = nint(values(ii+1:ii+nvert))
            else if (poolItr % dataType == MPAS_POOL_REAL) then
               fdata%r2%array(1:nvert, self%cellIndex) = real(values(ii+1:ii+nvert), RKIND)
            endif
         else
            write(message,*) '--> setpoint_2d: poolItr % nDims == ',poolItr % nDims,' not handled'
            call abor1_ftn(message)
         endif
         ii = ii + nvert
         deallocate(dimSizes)
      endif
   enddo

end subroutine setpoint_2d

! ------------------------------------------------------------------------------


end module iter2d_mod
