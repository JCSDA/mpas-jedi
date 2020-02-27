! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_covariance_mod

!oops
use kinds, only: kind_real

!ufo
use ufo_vars_mod, only: MAXVARLEN

implicit none

!> Fortran derived type to hold configuration data for the background/model covariance
type :: mpas_covar
  integer :: nothing_yet
  character(len=MAXVARLEN), allocatable :: var_scaling_variables(:)
  real (kind=kind_real), allocatable :: var_scaling_magnitudes(:)
end type mpas_covar

#define LISTED_TYPE mpas_covar

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_covar_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

!> Setup for the model's 3d error covariance matrices (B and Q_i)

!> This routine queries the configuration for the parameters that define the
!! covariance matrix, and stores the relevant values in the
!! error covariance structure.

subroutine mpas_covar_setup(self, geom, f_conf)
use fckit_configuration_module, only: fckit_configuration
use mpas_geom_mod
use iso_c_binding

implicit none
type(mpas_covar),          intent(inout) :: self    !< Covariance structure
type(fckit_configuration), intent(in)    :: f_conf  !< Configuration
type(mpas_geom),           intent(in)    :: geom    !< Geometry

character(kind=c_char,len=MAXVARLEN), allocatable :: char_array(:)
integer(c_size_t),parameter :: csize = MAXVARLEN
real(kind=c_float), allocatable :: real_array(:)

if (f_conf%has("var_scaling_variables") .and. f_conf%has("var_scaling_magnitudes")) then
   call f_conf%get_or_die("var_scaling_variables",csize,char_array)
   call f_conf%get_or_die("var_scaling_magnitudes",real_array)
   if(size(real_array) /= size(char_array)) then
      call abor1_ftn("--> mpas_b_setup_f90: var_scaling_variables and var_scaling_magnitudes have different sizes")
   end if

   allocate(self % var_scaling_variables(size(char_array)))
   self % var_scaling_variables = char_array
   allocate(self % var_scaling_magnitudes(size(real_array)))
   self % var_scaling_magnitudes = real(real_array,kind=kind_real)
else
   allocate(self % var_scaling_variables(0))
   allocate(self % var_scaling_magnitudes(0))
end if

end subroutine mpas_covar_setup

! ------------------------------------------------------------------------------

subroutine mpas_covar_delete(self)
implicit none
type(mpas_covar), intent(inout) :: self  !< Covariance structure

deallocate(self % var_scaling_variables)
deallocate(self % var_scaling_magnitudes)

end subroutine mpas_covar_delete

! ------------------------------------------------------------------------------

!> Multiply streamfunction by inverse(sqrt(C)), where C is 3d covariance matrix

subroutine mpas_covar_sqrt_inv_mult(self, xctl, xincr)
use mpas_field_utils_mod, only: mpas_field

implicit none
type(mpas_covar), intent(in)    :: self
real, intent(inout) :: xctl
type(mpas_field), intent(in)    :: xincr

end subroutine mpas_covar_sqrt_inv_mult

! ------------------------------------------------------------------------------

!> Multiply streamfunction by inverse(sqrt(C)) - Adjoint

subroutine mpas_covar_sqrt_inv_mult_ad(self, xctl, xincr)
use mpas_field_utils_mod, only: mpas_field

implicit none
type(mpas_covar), intent(in)    :: self
type(mpas_field), intent(inout) :: xincr
real, intent(in) :: xctl

end subroutine mpas_covar_sqrt_inv_mult_ad

! ------------------------------------------------------------------------------

!> Multiply streamfunction by sqrt(C), where C is a 3d covariance matrix

subroutine mpas_covar_sqrt_mult(self, xincr, xctl)
use mpas_field_utils_mod, only: mpas_field

implicit none
type(mpas_covar), intent(in)    :: self
type(mpas_field), intent(inout) :: xincr
real, intent(in) :: xctl

end subroutine mpas_covar_sqrt_mult

! ------------------------------------------------------------------------------

!> Multiply streamfunction by sqrt(C) - Adjoint

subroutine mpas_covar_sqrt_mult_ad(self, xincr, xctl)
use mpas_field_utils_mod, only: mpas_field

implicit none
type(mpas_covar), intent(in)    :: self
real, intent(inout) :: xctl
type(mpas_field), intent(in)    :: xincr

end subroutine mpas_covar_sqrt_mult_ad

! ------------------------------------------------------------------------------

end module mpas_covariance_mod
