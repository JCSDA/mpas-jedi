! (C) Copyright 2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> C++ interfaces for mpas_geom_iter_mod::mpas_geom_iter
module mpas_geom_iter_mod_c

use iso_c_binding
use kinds
use mpas_geom_iter_mod
use mpas_geom_mod, only: mpas_geom, mpas_geom_registry

implicit none
private

contains

! ------------------------------------------------------------------------------
!> C++ interface for mpas_geom_iter_mod::mpas_geom_iter::setup()
subroutine mpas_geom_iter_setup_c(c_key_self, c_key_geom, c_cellIndex, c_levIndex) bind(c, name='mpas_geom_iter_setup_f90')
  integer(c_int), intent(inout) :: c_key_self !< Geometry iterator
  integer(c_int), intent(   in) :: c_key_geom !< Geometry
  integer(c_int), intent(   in) :: c_cellIndex!< Index
  integer(c_int), intent(   in) :: c_levIndex !< Index

  ! Local variables
  type(mpas_geom_iter),     pointer :: self
  type(mpas_geom),          pointer :: geom

  ! Interface
  call mpas_geom_iter_registry%init()
  call mpas_geom_iter_registry%add(c_key_self)
  call mpas_geom_iter_registry%get(c_key_self, self)
  call mpas_geom_registry%get(c_key_geom, geom)

  ! Call Fortran
  call self%setup(geom, c_cellIndex, c_levIndex)

end subroutine mpas_geom_iter_setup_c


! ------------------------------------------------------------------------------
!> C++ interface for mpas_geom_iter_mod::mpas_geom_iter::clone()
subroutine mpas_geom_iter_clone_c(c_key_self, c_key_other) bind(c, name='mpas_geom_iter_clone_f90')
  integer(c_int), intent(inout) :: c_key_self  !< Geometry iterator
  integer(c_int), intent(   in) :: c_key_other !< Other geometry iterator

  ! Local variables
  type(mpas_geom_iter), pointer :: self, other

  ! Interface
  call mpas_geom_iter_registry%get(c_key_other, other)
  call mpas_geom_iter_registry%init()
  call mpas_geom_iter_registry%add(c_key_self)
  call mpas_geom_iter_registry%get(c_key_self, self)

  ! Call Fortran
  call self%clone(other)

end subroutine mpas_geom_iter_clone_c


! ------------------------------------------------------------------------------
!> !> C++ interface for deleting mpas_geom_iter_mod::mpas_geom_iter
subroutine mpas_geom_iter_delete_c(c_key_self) bind(c, name='mpas_geom_iter_delete_f90')
  integer(c_int), intent(inout) :: c_key_self !< Geometry iterator

  ! Local variables
  type(mpas_geom_iter), pointer :: self

  ! Clear interface
  call mpas_geom_iter_registry%get(c_key_self, self)
  call self%delete()
  call mpas_geom_iter_registry%remove(c_key_self)

end subroutine mpas_geom_iter_delete_c


! ------------------------------------------------------------------------------
!> C++ interface for mpas_geom_iter_mod::mpas_geom_iter::equals()
subroutine mpas_geom_iter_equals_c(c_key_self, c_key_other, c_equals) bind(c, name='mpas_geom_iter_equals_f90')
  integer(c_int), intent(inout)  :: c_key_self  !< Geometry iterator
  integer(c_int), intent(   in)  :: c_key_other !< Other geometry iterator
  logical(c_bool), intent(inout) :: c_equals    !< Equality boolean

  ! Local variables
  type(mpas_geom_iter), pointer :: self, other

  ! Interface
  call mpas_geom_iter_registry%get(c_key_self, self)
  call mpas_geom_iter_registry%get(c_key_other, other)

  ! Call Fortran
  c_equals = self%equals(other)

end subroutine mpas_geom_iter_equals_c


! ------------------------------------------------------------------------------
!> C++ interface for mpas_geom_iter_mod::mpas_geom_iter::current()
subroutine mpas_geom_iter_current_c(c_key_self, c_lon, c_lat, c_height) bind(c, name='mpas_geom_iter_current_f90')
  integer(c_int), intent(   in) :: c_key_self !< Geometry iterator
  real(c_double), intent(inout) :: c_lon      !< Longitude
  real(c_double), intent(inout) :: c_lat      !< Latitude
  real(c_double), intent(inout) :: c_height   !< Height

  ! Local variables
  type(mpas_geom_iter), pointer :: self

  ! Interface
  call mpas_geom_iter_registry%get(c_key_self, self)

  ! Call Fortran
  call self%current(c_lon, c_lat, c_height)

end subroutine mpas_geom_iter_current_c


! ------------------------------------------------------------------------------
!> C++ interface for mpas_geom_iter_mod::mpas_geom_iter::next()
subroutine mpas_geom_iter_next_c(c_key_self) bind(c, name='mpas_geom_iter_next_f90')
  integer(c_int), intent(in) :: c_key_self !< Geometry iterator

  ! Local variables
  type(mpas_geom_iter), pointer :: self

  ! Interface
  call mpas_geom_iter_registry%get(c_key_self, self)

  ! Call Fortran
  call self%next()
end subroutine mpas_geom_iter_next_c


! ------------------------------------------------------------------------------
!> C++ interface to get dimension of the iterator
subroutine mpas_geom_iter_dimension_c(c_key_self, c_val) bind(c, name='mpas_geom_iter_dimension_f90')
  integer(c_int), intent( in) :: c_key_self !< Geometry iterator
  integer(c_int), intent(out) :: c_val

  type(mpas_geom_iter), pointer :: self
  call mpas_geom_iter_registry%get(c_key_self, self)

  c_val = self%iterator%geom%iterator_dimension
end subroutine

end module mpas_geom_iter_mod_c
