! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_covariance_mod

implicit none

!> Fortran derived type to hold configuration data for the background/model covariance
type :: mpas_covar
  integer :: nothing_yet
end type mpas_covar

#define LISTED_TYPE mpas_covar

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_covar_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

!> Setup for the model's 3d error covariance matrices (B and Q_i)

!> This routine queries the configuration for the parameters that define the
!! covariance matrix, and stores the relevant values in the
!! error covariance structure.

subroutine mpas_covar_setup(self, geom, c_conf)
use iso_c_binding
use mpas_geom_mod
use config_mod

implicit none
type(mpas_covar), intent(inout) :: self    !< Covariance structure
type(c_ptr), intent(in)       :: c_conf  !< Configuration
type(mpas_geom), intent(in)     :: geom    !< Geometry

end subroutine mpas_covar_setup

! ------------------------------------------------------------------------------

subroutine mpas_covar_delete(self)
implicit none
type(mpas_covar), intent(inout) :: self  !< Covariance structure

end subroutine mpas_covar_delete

! ------------------------------------------------------------------------------

!> Multiply streamfunction by inverse(sqrt(C)), where C is 3d covariance matrix

subroutine mpas_covar_sqrt_inv_mult(self, xctl, xincr)
use mpas_increment_utils_mod, only: mpas_increment

implicit none
type(mpas_covar), intent(in)    :: self
real, intent(inout) :: xctl
type(mpas_increment), intent(in)    :: xincr

end subroutine mpas_covar_sqrt_inv_mult

! ------------------------------------------------------------------------------

!> Multiply streamfunction by inverse(sqrt(C)) - Adjoint

subroutine mpas_covar_sqrt_inv_mult_ad(self, xctl, xincr)
use mpas_increment_utils_mod, only: mpas_increment

implicit none
type(mpas_covar), intent(in)    :: self
type(mpas_increment), intent(inout) :: xincr
real, intent(in) :: xctl

end subroutine mpas_covar_sqrt_inv_mult_ad

! ------------------------------------------------------------------------------

!> Multiply streamfunction by sqrt(C), where C is a 3d covariance matrix

subroutine mpas_covar_sqrt_mult(self, xincr, xctl)
use mpas_increment_utils_mod, only: mpas_increment

implicit none
type(mpas_covar), intent(in)    :: self
type(mpas_increment), intent(inout) :: xincr
real, intent(in) :: xctl

end subroutine mpas_covar_sqrt_mult

! ------------------------------------------------------------------------------

!> Multiply streamfunction by sqrt(C) - Adjoint

subroutine mpas_covar_sqrt_mult_ad(self, xincr, xctl)
use mpas_increment_utils_mod, only: mpas_increment

implicit none
type(mpas_covar), intent(in)    :: self
real, intent(inout) :: xctl
type(mpas_increment), intent(in)    :: xincr

end subroutine mpas_covar_sqrt_mult_ad

! ------------------------------------------------------------------------------

end module mpas_covariance_mod
