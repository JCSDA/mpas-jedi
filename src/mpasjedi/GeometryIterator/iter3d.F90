! (C) Copyright 2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Geometry iterator
module iter3d_mod

use iso_c_binding

!oops
use kinds, only : kind_real
use missing_values_mod

!MPAS-Model
use mpas_derived_types
use mpas_kind_types, only: RKIND

!mpas-jedi
use iter_mod
use mpas_constants_mod
use mpas_fields_mod, only: mpas_fields
use mpas_geom_mod, only: getVertLevels

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
  end if !cell loop

  if (self%levIndex > self%geom%nVertLevels) then
    self%cellIndex=-1
    self%levIndex=-1
  end if !lev loop
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
  integer :: ii, iVar
  type(mpas_pool_data_type), pointer :: fdata

  values(:) = missing_value(values(1))

  ii = 0

  do iVar = 1, fields % nf

    call fields%get(fields%fldnames(iVar), fdata)

    if (associated(fdata%r2)) then
      values(ii+1) = real(fdata%r2%array(self%levIndex, self%cellIndex), kind_real)

    else if (associated(fdata%r1)) then
      if (self%levIndex == 1) then
        values(ii+1) = real(fdata%r1%array(self%cellIndex), kind_real)
      end if

    else
      call abor1_ftn('getpoint_3d: only r1 and r2 data types are supported')
    end if

    ii = ii + 1

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
  integer :: ii, iVar
  type(mpas_pool_data_type), pointer :: fdata

  ii = 0

  do iVar = 1, fields % nf

    call fields%get(fields%fldnames(iVar), fdata)

    if (associated(fdata%r2)) then
      fdata%r2%array(self%levIndex, self%cellIndex) = real(values(ii+1), RKIND)

    else if (associated(fdata%r1)) then
      if (self%levIndex == 1) then
        fdata%r1%array(self%cellIndex) = real(values(ii+1), RKIND)
      end if

    else
      call abor1_ftn('setpoint_3d: only r1 and r2 data types are supported')
    end if

    ii = ii + 1

  enddo

end subroutine setpoint_3d

! ------------------------------------------------------------------------------

end module iter3d_mod
