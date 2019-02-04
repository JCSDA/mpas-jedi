! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_state_utils_mod

use iso_c_binding

!UFO
use ufo_vars_mod

!MPAS-JEDI
use mpas_field_utils_mod
use mpas_geom_mod

implicit none

private

public :: mpas_state, mpas_state_registry, nf_aux

! ------------------------------------------------------------------------------

   !> Fortran derived type to hold MPAS state
   type, extends(mpas_field) :: mpas_state
      private
      contains
      procedure :: create => create_state
   end type mpas_state
 
#define LISTED_TYPE mpas_state

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_state_registry

integer, parameter :: nf_aux = 21
character(len=MAXVARLEN) :: fldnames_aux(nf_aux)

! ------------------------------------------------------------------------------
contains

! ------------------------------------------------------------------------------

!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine create_state(self, geom, vars)

    implicit none

    class(mpas_state), intent(inout) :: self
    type(mpas_geom),   intent(in), pointer :: geom
    type(ufo_vars),    intent(in)          :: vars

    type(ufo_vars) :: state_vars

    integer :: nfields
    character(len=1024) :: buf
    integer, save :: instances = 0

    if (instances == 0) then
       ! Temporary Variable Definitions
       !--- TODO: theta, rho, u should come from "State" variables in YAML
       !          + not enabled yet for variational runs in OOPS; 
       !          + var list for oops::State4D currently comes from
       !            cost_function/variables, but needs to come from
       !            cost_function/Jb/Background/state/variables
       !          + vars_ in ModelMPAS, TlmMPAS, and TlmIdMPAS will need to be 
       !            updated to reflect any changes
       !          + Note: it may be possible to merge State and Increment back
       !            into Fields class once variables are defined properly in OOPS
       !--- TODO: all other aux variables should come from "VariableChange" 
       !          variables in YAML (not enabled yet in OOPS)
       fldnames_aux = [ character(len=MAXVARLEN) :: &
          "theta", "rho", "u", &
          "landmask", "xice", "snowc", "skintemp", "ivgtyp", "isltyp", &
          "snowh", "vegfra", "u10", "v10", "lai", "smois", "tslb", "w", &
          "index_qc", "index_qi", "re_cloud", "re_ice" ]
    end if

    instances = instances + 1

    state_vars % nv = vars % nv + nf_aux
    allocate(state_vars % fldnames(state_vars % nv))
    state_vars % fldnames(1:vars%nv) = vars % fldnames(:)
    state_vars % fldnames(vars%nv+1:state_vars%nv) = fldnames_aux(:)
 
    call create_field(self, geom, state_vars)

    deallocate(self % fldnames_ci)
    self % nf_ci = vars % nv
    allocate(self % fldnames_ci(self % nf_ci))
    self % fldnames_ci(:) = vars % fldnames(:)

    call release_label(self % field_number)
    call get_new_label(self % field_number, .true.)
    write(*,*)'--> create_state: state number ',self % field_number

end subroutine create_state

! ------------------------------------------------------------------------------

end module mpas_state_utils_mod
