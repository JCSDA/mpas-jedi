! (C) Copyright 2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Geometry iterator
module iter_mod

use mpas_geom_mod, only: mpas_geom

implicit none
private

public :: iter
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Geometry iterator
!!
!! When initialized, the iterator points to the first valid local grid cell.
!! Calls to iter::next() moves the iterator forward, and calls to
!! iter::current() retrieves the lat/lon/height of the current grid cell.
!! The iterator is mainly used by Increment::getLocal and Increment::setLocal.
type, abstract :: iter
  type(mpas_geom), pointer :: geom => null() !< Geometry

  integer :: cellIndex = 1 !< cell index of current grid point
  integer :: levIndex = 1  !< level index of current grid point

  contains

  !> \copybrief iter_setup \see iter_setup
  procedure :: setup => iter_setup

  ! methods below are deferred to extensions of iter

  !> \copybrief equals_iter
  procedure(equals_iter), deferred :: equals

  !> \copybrief current_iter
  procedure(current_iter), deferred :: current

  !> \copybrief next_iter
  procedure(next_iter), deferred :: next

  !> \copybrief getpoint_iter
  procedure(getpoint_iter), deferred :: getpoint

  !> \copybrief setpoint_iter
  procedure(setpoint_iter), deferred :: setpoint

end type iter

interface
  logical(kind=c_bool) function equals_iter(self, other)
    use iso_c_binding
    import :: iter
    class(iter), intent(in) :: self  !< Geometry iterator
    class(iter), pointer, intent(in) :: other !< Geometry iterator
  end function
  subroutine current_iter(self, lon, lat, height)
    use kinds, only : kind_real
    import :: iter
    class(iter), intent(in) :: self !< Geometry iterator
    real(kind_real), intent(out) :: lon    !< Longitude
    real(kind_real), intent(out) :: lat    !< Latitude
    real(kind_real), intent(out) :: height !< Height
  end subroutine current_iter
  subroutine next_iter(self)
    import :: iter
    class(iter), intent(inout) :: self !< Geometry iterator
  end subroutine next_iter
  subroutine getpoint_iter(self, fields, values, nval)
    use iso_c_binding
    use kinds, only : kind_real
    use mpas_fields_mod, only: mpas_fields
    import :: iter
    class(iter), intent(in) :: self
    class(mpas_fields), intent(in) :: fields
    integer(c_int) :: nval
    real(kind_real), intent(inout) :: values(nval)
  end subroutine getpoint_iter
  subroutine setpoint_iter(self, fields, values, nval)
    use iso_c_binding
    use kinds, only : kind_real
    use mpas_fields_mod, only: mpas_fields
    import :: iter
    class(iter), intent(in) :: self
    class(mpas_fields), intent(inout) :: fields
    integer(c_int) :: nval
    real(kind_real), intent(in) :: values(nval)
  end subroutine setpoint_iter
end interface

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Setup the geometry iterator
!!
!! \relates iter_mod::iter
subroutine iter_setup(self, geom, cellIndex, levIndex)
  class(iter),              intent(inout) :: self
  type(mpas_geom), pointer, intent(in)    :: geom !< Pointer to geometry
  integer,                  intent(in)    :: cellIndex, levIndex  !< starting index

  ! Associate geometry
  self%geom => geom

  ! Define cellIndex/levIndex for local tile
  self%cellIndex = cellIndex
  self%levIndex = levIndex
end subroutine iter_setup

! ------------------------------------------------------------------------------

end module iter_mod
