! (C) Copyright 2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Geometry iterator
module iter_mod

use iso_c_binding

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

  !> \copybrief iter_equals
  procedure(iter_equals), deferred :: equals

  !> \copybrief iter_current
  procedure(iter_current), deferred :: current

  !> \copybrief iter_next
  procedure(iter_next), deferred :: next

  !> \copybrief iter_getpoint
  procedure(iter_getpoint), deferred :: getpoint

  !> \copybrief iter_setpoint
  procedure(iter_setpoint), deferred :: setpoint

end type iter

interface
  logical(kind=c_bool) function iter_equals(self, other)
    use iso_c_binding
    import :: iter
    class(iter), intent(in) :: self  !< Geometry iterator
    class(iter), pointer, intent(in) :: other !< Geometry iterator
  end function
  subroutine iter_current(self, lon, lat, height)
    use kinds, only : kind_real
    import :: iter
    class(iter), intent(in) :: self !< Geometry iterator
    real(kind_real), intent(out) :: lon    !< Longitude
    real(kind_real), intent(out) :: lat    !< Latitude
    real(kind_real), intent(out) :: height !< Height
  end subroutine iter_current
  subroutine iter_next(self)
    import :: iter
    class(iter), intent(inout) :: self !< Geometry iterator
  end subroutine iter_next
  subroutine iter_getpoint(self, fields, values, nval)
    use iso_c_binding
    use kinds, only : kind_real
    use mpas_fields_mod, only: mpas_fields
    import :: iter
    class(iter), intent(in) :: self
    class(mpas_fields), intent(in) :: fields
    integer(c_int) :: nval
    real(kind_real), intent(inout) :: values(nval)
  end subroutine iter_getpoint
  subroutine iter_setpoint(self, fields, values, nval)
    use iso_c_binding
    use kinds, only : kind_real
    use mpas_fields_mod, only: mpas_fields
    import :: iter
    class(iter), intent(in) :: self
    class(mpas_fields), intent(inout) :: fields
    integer(c_int) :: nval
    real(kind_real), intent(in) :: values(nval)
  end subroutine iter_setpoint
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
  integer(c_int),           intent(in)    :: cellIndex, levIndex  !< starting index

  ! Associate geometry
  self%geom => geom

  ! Define cellIndex/levIndex for local tile
  self%cellIndex = cellIndex
  self%levIndex = levIndex
end subroutine iter_setup

! ------------------------------------------------------------------------------

end module iter_mod
