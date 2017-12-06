! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_geom_mod

use iso_c_binding
use config_mod

implicit none
private
public :: mpas_geom, &
        & geo_setup, geo_clone, geo_delete, geo_info
public :: mpas_geom_registry


! ------------------------------------------------------------------------------

!> Fortran derived type to hold geometry definition
type :: mpas_geom
  integer :: nothing
end type mpas_geom

#define LISTED_TYPE mpas_geom

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_geom_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------
subroutine geo_setup(self, c_conf)
implicit none
type(mpas_geom)        :: self
type(c_ptr), intent(in) :: c_conf

end subroutine geo_setup

! ------------------------------------------------------------------------------

subroutine geo_clone(self, other)
implicit none
type(mpas_geom), pointer :: self, other

end subroutine geo_clone

! ------------------------------------------------------------------------------

subroutine geo_delete(self)

implicit none
type(mpas_geom), pointer :: self

end subroutine geo_delete

! ------------------------------------------------------------------------------

subroutine geo_info(self)
implicit none
type(mpas_geom), pointer :: self

end subroutine geo_info

! ------------------------------------------------------------------------------

end module mpas_geom_mod
