! (C) Copyright 2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Geometry iterator
module iter3d_mod

use iso_c_binding
use kinds, only : kind_real
use missing_values_mod

use mpas_derived_types
use mpas_pool_routines
use mpas_kind_types, only: RKIND

use iter_mod
use mpas_constants_mod
use mpas_fields_mod, only: mpas_fields
use mpas_geom_mod, only: getSolveDimSizes

implicit none
private

public :: iter3d
type, extends(iter) :: iter3d
  contains
  procedure :: equals => equals_3d
  procedure :: current => current_3d
  procedure :: next => next_3d
  procedure :: getpoint => getpoint_3d
  procedure :: setpoint => setpoint_3d
end type iter3d

integer, parameter    :: max_string=8000
character(max_string) :: message

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Check for the geometry iterator equality (pointing to same cellIndex/levIndex)
!!
!! \relates iter3d_mod::iter3d
logical(kind=c_bool) function equals_3d(self, other)
  class(iter3d), intent(in) :: self
  class(iter), pointer, intent(in) :: other
  equals_3d = ( (self%cellIndex==other%cellIndex) .and. (self%levIndex==other%levIndex) )
end function equals_3d

! ------------------------------------------------------------------------------
!> Get geometry iterator current lat/lon/height
!!
!! \throws abor1_ftn aborts if iterator is out of bounds
!! \relates iter3d_mod::iter3d
subroutine current_3d(self, lon, lat, height)
  ! Passed variables
  class(iter3d), intent( in) :: self !< Geometry iterator
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
    call abor1_ftn('current_3d: lat/lon iterator out of bounds')
  else
    ! inside of the grid
    ii = self%cellIndex
  endif
  lat = self%geom%latCell(ii) * MPAS_JEDI_RAD2DEG_kr
  lon = self%geom%lonCell(ii) * MPAS_JEDI_RAD2DEG_kr

  ! check levIndex
  if (self%levIndex == -1) then
    ! special case of {-1} means end of the grid
    height = self%geom%height(self%geom%nVertLevels, ii)
  elseif (self%levIndex == 0) then
    ! special case for surface fields; not currently used
    height = self%geom%zgrid(1, ii);
  elseif (self%levIndex < 1 .OR. self%levIndex > self%geom%nVertLevels) then
    ! out of range
    call abor1_ftn('current_3d: depth iterator out of bounds')
  else
    ! inside of the 3D grid
    height = self%geom%height(self%levIndex, ii)
  endif
end subroutine current_3d

! ------------------------------------------------------------------------------
!> Update geometry iterator to next point
!!
!! \relates iter3d_mod::iter3d
subroutine next_3d(self)
  class(iter3d), intent(inout) :: self
  if (self%cellIndex.lt.self%geom%nCellsSolve) then
    self%cellIndex = self%cellIndex + 1
  elseif (self%cellIndex.eq.self%geom%nCellsSolve) then
    self%cellIndex = 1
    self%levIndex = self%levIndex + 1
  end if !iloop

  if (self%levIndex > self%geom%nVertLevels) then
    self%cellIndex=-1
    self%levIndex=-1
  end if !kloop
end subroutine next_3d

! ------------------------------------------------------------------------------

subroutine getpoint_3d(self, fields, values, nval)

   implicit none

   ! Passed variables
   class(iter3d), intent(in) :: self
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
            if (self%levIndex == 1) then
              if (poolItr % dataType == MPAS_POOL_INTEGER) then
                 values(ii+1) = real(fdata%i1%array(self%cellIndex), kind_real)
              else if (poolItr % dataType == MPAS_POOL_REAL) then
                 values(ii+1) = real(fdata%r1%array(self%cellIndex), kind_real)
              endif
            else
              values(ii+1) = missing_value(values(1))
            end if
         else if (poolItr % nDims == 2) then
            nvert = dimSizes(1)
            values(ii+1:ii+nvert) = missing_value(values(1))
            if (poolItr % dataType == MPAS_POOL_INTEGER) then
               values(ii+self%levIndex) = real(fdata%i2%array(self%levIndex, self%cellIndex), kind_real)
            else if (poolItr % dataType == MPAS_POOL_REAL) then
               values(ii+self%levIndex) = real(fdata%r2%array(self%levIndex, self%cellIndex), kind_real)
            endif
         else
            write(message,*) '--> getpoint_3d: poolItr % nDims == ',poolItr % nDims,' not handled'
            call abor1_ftn(message)
         endif
         ii = ii + nvert
         deallocate(dimSizes)
      endif
   enddo

end subroutine getpoint_3d

! ------------------------------------------------------------------------------

subroutine setpoint_3d(self, fields, values, nval)

   implicit none

   ! Passed variables
   class(iter3d), intent(in) :: self
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
            if (self%levIndex == 1) then
              if (poolItr % dataType == MPAS_POOL_INTEGER) then
                 fdata%i1%array(self%cellIndex) = nint(values(ii+1))
              else if (poolItr % dataType == MPAS_POOL_REAL) then
                 fdata%r1%array(self%cellIndex) = real(values(ii+1), RKIND)
              endif
            endif
         else if (poolItr % nDims == 2) then
            nvert = dimSizes(1)
            if (poolItr % dataType == MPAS_POOL_INTEGER) then
               fdata%i2%array(self%levIndex, self%cellIndex) = nint(values(ii+self%levIndex))
            else if (poolItr % dataType == MPAS_POOL_REAL) then
               fdata%r2%array(self%levIndex, self%cellIndex) = real(values(ii+self%levIndex), RKIND)
            endif
         else
            write(message,*) '--> setpoint_3d: poolItr % nDims == ',poolItr % nDims,' not handled'
            call abor1_ftn(message)
         endif
         ii = ii + nvert
         deallocate(dimSizes)
      endif
   enddo

end subroutine setpoint_3d

! ------------------------------------------------------------------------------

end module iter3d_mod
