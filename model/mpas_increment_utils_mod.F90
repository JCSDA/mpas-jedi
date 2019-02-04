! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_increment_utils_mod

!MPAS-JEDI
use mpas_field_utils_mod

implicit none

private

public :: mpas_increment, mpas_increment_registry

! ------------------------------------------------------------------------------

   !> Fortran derived type to hold MPAS increment
   type, extends(mpas_field) :: mpas_increment
      private
      contains
   end type mpas_increment
 
#define LISTED_TYPE mpas_increment

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_increment_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

end module mpas_increment_utils_mod
