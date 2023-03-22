! (C) Copyright 2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Geometry iterator
module mpas_geom_iter_mod

use iso_c_binding

use kinds, only : kind_real

use mpas_geom_mod, only: mpas_geom
use iter_mod
use iter2d_mod
use iter3d_mod

implicit none
private

public :: mpas_geom_iter
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
!> Geometry iterator container
!!
!! Generic container for 2d and 3d iterators
type :: mpas_geom_iter
  class(iter), pointer :: iterator => null()
  class(iter2d), pointer :: i2d => null()
  class(iter3d), pointer :: i3d => null()

  contains

  !> \copybrief mpas_geom_iter_setup \see mpas_geom_iter_setup
  procedure :: setup => mpas_geom_iter_setup

  !> \copybrief mpas_geom_iter_clone \see mpas_geom_iter_clone
  procedure :: clone => mpas_geom_iter_clone

  !> \copybrief mpas_geom_iter_equals \see mpas_geom_iter_equals
  procedure :: equals => mpas_geom_iter_equals

  !> \copybrief mpas_geom_iter_current \see mpas_geom_iter_current
  procedure :: current => mpas_geom_iter_current

  !> \copybrief mpas_geom_iter_next \see mpas_geom_iter_next
  procedure :: next => mpas_geom_iter_next

  !> \copybrief mpas_geom_iter_delete \see mpas_geom_iter_delete
  procedure :: delete => mpas_geom_iter_delete

end type mpas_geom_iter

#define LISTED_TYPE mpas_geom_iter

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry for mpas_geom_iter instances
type(registry_t), public :: mpas_geom_iter_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
!> Setup the geometry iterator
!!
!! \relates mpas_geom_iter_mod::mpas_geom_iter
subroutine mpas_geom_iter_setup(self, geom, cellIndex, levIndex)
  class(mpas_geom_iter),    intent(inout) :: self
  type(mpas_geom), pointer, intent(in)    :: geom !< Pointer to geometry
  integer(c_int),           intent(in)    :: cellIndex, levIndex  !< starting index

  select case (geom%iterator_dimension)
    case (2)
      allocate(self%i2d)
      self % iterator => self%i2d
    case (3)
      allocate(self%i3d)
      self % iterator => self%i3d
    case default
      call abor1_ftn('mpas_geom_iter_setup: unknown geom%iterator_dimension')
  end select

  call self%iterator%setup(geom, cellIndex, levIndex)

end subroutine mpas_geom_iter_setup

! ------------------------------------------------------------------------------
!> Delete the geometry iterator
!!
!! \relates mpas_geom_iter_mod::mpas_geom_iter
subroutine mpas_geom_iter_delete(self)
  class(mpas_geom_iter), intent(inout) :: self
  if (associated(self%i2d)) deallocate(self%i2d)
  if (associated(self%i3d)) deallocate(self%i3d)
end subroutine mpas_geom_iter_delete

! ------------------------------------------------------------------------------
!> Clone the geometry iterator from \p other to \p self
!!
!! \relates mpas_geom_iter_mod::mpas_geom_iter
subroutine mpas_geom_iter_clone(self, other)
  class(mpas_geom_iter), intent(inout) :: self
  type(mpas_geom_iter), intent(in)     :: other !< Other geometry iterator to clone from
  call self%setup( &
    other%iterator%geom, &
    other%iterator%cellIndex, &
    other%iterator%levIndex)
end subroutine mpas_geom_iter_clone

! ------------------------------------------------------------------------------
!> Check for the geometry iterator equality (pointing to same cellIndex)
!!
!! \relates mpas_geom_iter_mod::mpas_geom_iter
logical(kind=c_bool) function mpas_geom_iter_equals(self, other)
  class(mpas_geom_iter), intent(inout) :: self
  type(mpas_geom_iter), intent(in)     :: other
  mpas_geom_iter_equals = self%iterator%equals(other%iterator)
end function mpas_geom_iter_equals

! ------------------------------------------------------------------------------
!> Get geometry iterator current lat/lon/height
!!
!! \throws abor1_ftn aborts if iterator is out of bounds
!! \relates mpas_geom_iter_mod::mpas_geom_iter
subroutine mpas_geom_iter_current(self, lon, lat, height)
  ! Passed variables
  class(mpas_geom_iter), intent( in) :: self !< Geometry iterator
  real(kind_real), intent(out) :: lon    !< Longitude
  real(kind_real), intent(out) :: lat    !< Latitude
  real(kind_real), intent(out) :: height !< Height
  call self%iterator%current(lon, lat, height)
end subroutine mpas_geom_iter_current

! ------------------------------------------------------------------------------
!> Update geometry iterator to next point (2d)
!!
!! \relates mpas_geom_iter_mod::mpas_geom_iter
subroutine mpas_geom_iter_next(self)
  class(mpas_geom_iter), intent(inout) :: self
  call self%iterator%next()
end subroutine mpas_geom_iter_next

! ------------------------------------------------------------------------------

end module mpas_geom_iter_mod
