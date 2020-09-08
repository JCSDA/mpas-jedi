! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_linvarcha_c2a_mod

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only : fckit_log
use iso_c_binding
use kinds

!mpas-jedi
use mpas_constants_mod
use mpas_field_utils_mod, only: mpas_field
use mpas_geom_mod, only: mpas_geom

!MPAS-Model
use mpas_derived_types
use mpas_field_routines
use mpas_pool_routines
use mpas_dmpar

implicit none
private
public :: mpas_linvarcha_c2a, mpas_linvarcha_c2a_registry, &
          mpas_linvarcha_c2a_setup, mpas_linvarcha_c2a_delete, &
          mpas_linvarcha_c2a_multiply, mpas_linvarcha_c2a_multiplyadjoint, &
          mpas_linvarcha_c2a_multiplyinverse, mpas_linvarcha_c2a_multiplyinverseadjoint

!> Fortran derived type to hold configuration data for the B mat variable change
type :: mpas_linvarcha_c2a
  integer :: wind_cvt_method  !< method for wind conversion
end type mpas_linvarcha_c2a

integer, parameter    :: max_string=800
character(max_string) :: err_msg

#define LISTED_TYPE mpas_linvarcha_c2a

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_linvarcha_c2a_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------

subroutine mpas_linvarcha_c2a_setup(self, bg, fg, geom, f_conf)

   implicit none
   type(mpas_linvarcha_c2a),  intent(inout) :: self    !< Change variable structure
   type(mpas_field), target, intent(in) :: bg
   type(mpas_field), target, intent(in) :: fg
   type(mpas_geom),          intent(in) :: geom
   type(fckit_configuration), intent(in)    :: f_conf  !< Configuration

   integer :: wind_cvt_method

   if (f_conf%has("wind_cvt_method")) then
      call f_conf%get_or_die("wind_cvt_method",wind_cvt_method)
      if( wind_cvt_method .ne. 1 .and. wind_cvt_method .ne. 2 ) &
         call abor1_ftn("LinVarChaC2AMPAS::LinVarChaC2AMPAS, mpas_linvarcha_c2a_setup: wind_cvt_method should be 1 or 2")
   else
      wind_cvt_method = 1 ! default method
   end if

   self % wind_cvt_method = wind_cvt_method
   write(err_msg,*) 'DEBUG:: mpas_linvarcha_c2a_setup, self % wind_cvt_method=',self % wind_cvt_method
   call fckit_log%debug(err_msg)

end subroutine mpas_linvarcha_c2a_setup

! ------------------------------------------------------------------------------

subroutine mpas_linvarcha_c2a_delete(self)

   implicit none
   type(mpas_linvarcha_c2a), intent(inout) :: self

   self % wind_cvt_method = 0

end subroutine mpas_linvarcha_c2a_delete

! ------------------------------------------------------------------------------

subroutine mpas_linvarcha_c2a_multiply(self,xctl,xana)

   implicit none
   type(mpas_linvarcha_c2a), intent(inout) :: self
   type(mpas_field), intent(inout) :: xctl
   type(mpas_field), intent(inout) :: xana

   type (field2DReal), pointer :: field2d_sf, field2d_vp, field2d_uRz, field2d_uRm
   type (field2DReal), pointer :: field2d_in, field2d_out
   type (field1DReal), pointer :: field1d_in, field1d_out
   integer :: k, ngrid

   ngrid = xctl % geom % nCells ! local + halo

   write(err_msg,*) "DEBUG: mpas_linvarcha_c2a_multiply: xana % fldnames(:) =",xana % fldnames(:)
   call fckit_log%debug(err_msg)
   write(err_msg,*) "DEBUG: mpas_linvarcha_c2a_multiply: xctl % fldnames(:) =",xctl % fldnames(:)
   call fckit_log%debug(err_msg)

   call mpas_pool_get_field(xctl % subFields,            'temperature', field2d_in)
   call mpas_pool_get_field(xana % subFields,            'temperature', field2d_out)
   field2d_out%array(:,:)=field2d_in%array(:,:)
   call mpas_pool_get_field(xctl % subFields,            'spechum', field2d_in)
   call mpas_pool_get_field(xana % subFields,            'spechum', field2d_out)
   field2d_out%array(:,:)=field2d_in%array(:,:)
   call mpas_pool_get_field(xctl % subFields,   'surface_pressure', field1d_in)
   call mpas_pool_get_field(xana % subFields,   'surface_pressure', field1d_out)
   field1d_out%array(:)=field1d_in%array(:)

   call mpas_pool_get_field(xctl % subFields,        'stream_function', field2d_sf)
   call mpas_pool_get_field(xctl % subFields,     'velocity_potential', field2d_vp)
   call mpas_pool_get_field(xana % subFields,      'uReconstructZonal', field2d_uRz)
   call mpas_pool_get_field(xana % subFields, 'uReconstructMeridional', field2d_uRm)

   call mpas_dmpar_exch_halo_field(field2d_sf)
   call mpas_dmpar_exch_halo_field(field2d_vp)

   if( self % wind_cvt_method .eq. 1) then
      do k=1,xctl%geom%nVertLevels
         call psichi_to_uv( xctl % geom, ngrid, field2d_sf%array(k,1:ngrid), field2d_vp%array(k,1:ngrid), &
                            field2d_uRz%array(k,1:ngrid), field2d_uRm%array(k,1:ngrid) )
      end do
   else if( self % wind_cvt_method .eq. 2) then
      do k=1,xctl%geom%nVertLevels
         call psichi_to_uv2( xctl % geom, field2d_sf%array(k,1:ngrid), field2d_vp%array(k,1:ngrid), &
                             field2d_uRz%array(k,1:ngrid), field2d_uRm%array(k,1:ngrid) )
      end do
   else
      call abor1_ftn("LinVarChaC2AMPAS:: This method is not implemented")
   end if

end subroutine mpas_linvarcha_c2a_multiply

! ------------------------------------------------------------------------------

subroutine mpas_linvarcha_c2a_multiplyadjoint(self,xana,xctl)

   implicit none
   type(mpas_linvarcha_c2a), intent(inout) :: self
   type(mpas_field), intent(inout) :: xana
   type(mpas_field), intent(inout) :: xctl

   type (field2DReal), pointer :: field2d_sf, field2d_vp, field2d_uRz, field2d_uRm
   type (field2DReal), pointer :: field2d_in, field2d_out
   type (field1DReal), pointer :: field1d_in, field1d_out
   integer :: k, ngrid

   ngrid = xctl % geom % nCells

   write(err_msg,*) "DEBUG: mpas_linvarcha_c2a_multiplyadjoint: xana % fldnames(:) =",xana % fldnames(:)
   call fckit_log%debug(err_msg)
   write(err_msg,*) "DEBUG: mpas_linvarcha_c2a_multiplyadjoint: xctl % fldnames(:) =",xctl % fldnames(:)
   call fckit_log%debug(err_msg)

   call mpas_pool_get_field(xana % subFields,            'temperature', field2d_in)
   call mpas_pool_get_field(xctl % subFields,            'temperature', field2d_out)
   field2d_out%array(:,:)=field2d_out%array(:,:)+field2d_in%array(:,:)
   field2d_in%array(:,:)=MPAS_JEDI_ZERO_kr
   call mpas_pool_get_field(xana % subFields,            'spechum', field2d_in)
   call mpas_pool_get_field(xctl % subFields,            'spechum', field2d_out)
   field2d_out%array(:,:)=field2d_out%array(:,:)+field2d_in%array(:,:)
   field2d_in%array(:,:)=MPAS_JEDI_ZERO_kr
   call mpas_pool_get_field(xana % subFields,   'surface_pressure', field1d_in)
   call mpas_pool_get_field(xctl % subFields,   'surface_pressure', field1d_out)
   field1d_out%array(:)=field1d_out%array(:)+field1d_in%array(:)
   field1d_in%array(:)=MPAS_JEDI_ZERO_kr

   call mpas_pool_get_field(xctl % subFields,        'stream_function', field2d_sf)
   call mpas_pool_get_field(xctl % subFields,     'velocity_potential', field2d_vp)
   call mpas_pool_get_field(xana % subFields,      'uReconstructZonal', field2d_uRz)
   call mpas_pool_get_field(xana % subFields, 'uReconstructMeridional', field2d_uRm)

   if( self % wind_cvt_method .eq. 1) then
      do k=1,xctl%geom%nVertLevels
         call psichi_to_uvAD( xctl % geom, ngrid, field2d_sf%array(k,1:ngrid), field2d_vp%array(k,1:ngrid), &
                              field2d_uRz%array(k,1:ngrid), field2d_uRm%array(k,1:ngrid) )
      end do
   else if( self % wind_cvt_method .eq. 2) then
      do k=1,xctl%geom%nVertLevels
         call psichi_to_uv2AD( xctl % geom, field2d_sf%array(k,1:ngrid), field2d_vp%array(k,1:ngrid), &
                               field2d_uRz%array(k,1:ngrid), field2d_uRm%array(k,1:ngrid) )
      end do
   else
      call abor1_ftn("LinVarChaC2AMPAS:: This method is not implemented")
   end if
   field2d_uRz%array(:,:)=MPAS_JEDI_ZERO_kr
   field2d_uRm%array(:,:)=MPAS_JEDI_ZERO_kr

   call mpas_dmpar_exch_halo_adj_field(field2d_sf)
   call mpas_dmpar_exch_halo_adj_field(field2d_vp)

end subroutine mpas_linvarcha_c2a_multiplyadjoint

! ------------------------------------------------------------------------------

subroutine mpas_linvarcha_c2a_multiplyinverse(self,xana,xctl)

   implicit none
   type(mpas_linvarcha_c2a), intent(inout) :: self
   type(mpas_field), intent(inout) :: xana
   type(mpas_field), intent(inout) :: xctl

   type (field2DReal), pointer :: field2d_in, field2d_out
   type (field1DReal), pointer :: field1d_in, field1d_out

   write(err_msg,*) "DEBUG: mpas_linvarcha_c2a_multiplyinverse: xana % fldnames(:) =",xana % fldnames(:)
   call fckit_log%debug(err_msg)
   write(err_msg,*) "DEBUG: mpas_linvarcha_c2a_multiplyinverse: xctl % fldnames(:) =",xctl % fldnames(:)
   call fckit_log%debug(err_msg)

   call mpas_pool_get_field(xana % subFields,            'temperature', field2d_in)
   call mpas_pool_get_field(xctl % subFields,            'temperature', field2d_out)
   field2d_out%array(:,:)=field2d_in%array(:,:)
   call mpas_pool_get_field(xana % subFields,            'temperature', field2d_in)
   call mpas_pool_get_field(xctl % subFields,            'temperature', field2d_out)
   field2d_out%array(:,:)=field2d_in%array(:,:)
   call mpas_pool_get_field(xana % subFields,            'spechum', field2d_in)
   call mpas_pool_get_field(xctl % subFields,            'spechum', field2d_out)
   field2d_out%array(:,:)=field2d_in%array(:,:)
   call mpas_pool_get_field(xana % subFields,      'uReconstructZonal', field2d_in)
   call mpas_pool_get_field(xctl % subFields,      'stream_function'  , field2d_out)
   field2d_out%array(:,:)=field2d_in%array(:,:)
   call mpas_pool_get_field(xana % subFields, 'uReconstructMeridional', field2d_in)
   call mpas_pool_get_field(xctl % subFields, 'velocity_potential'    , field2d_out)
   field2d_out%array(:,:)=field2d_in%array(:,:)
   call mpas_pool_get_field(xana % subFields,   'surface_pressure', field1d_in)
   call mpas_pool_get_field(xctl % subFields,   'surface_pressure', field1d_out)
   field1d_out%array(:)=field1d_in%array(:)

end subroutine mpas_linvarcha_c2a_multiplyinverse

! ------------------------------------------------------------------------------

subroutine mpas_linvarcha_c2a_multiplyinverseadjoint(self,xctl,xana)

   implicit none
   type(mpas_linvarcha_c2a), intent(inout) :: self
   type(mpas_field), intent(inout) :: xctl
   type(mpas_field), intent(inout) :: xana

   type (field2DReal), pointer :: field2d_in, field2d_out
   type (field1DReal), pointer :: field1d_in, field1d_out

   write(err_msg,*) "DEBUG: mpas_linvarcha_c2a_multiplyinverseadjoint: xana % fldnames(:) =",xana % fldnames(:)
   call fckit_log%debug(err_msg)
   write(err_msg,*) "DEBUG: mpas_linvarcha_c2a_multiplyinverseadjoint: xctl % fldnames(:) =",xctl % fldnames(:)
   call fckit_log%debug(err_msg)

   call mpas_pool_get_field(xctl % subFields,            'temperature', field2d_in)
   call mpas_pool_get_field(xana % subFields,            'temperature', field2d_out)
   field2d_out%array(:,:)=field2d_in%array(:,:)
   call mpas_pool_get_field(xctl % subFields,            'temperature', field2d_in)
   call mpas_pool_get_field(xana % subFields,            'temperature', field2d_out)
   field2d_out%array(:,:)=field2d_in%array(:,:)
   call mpas_pool_get_field(xctl % subFields,            'spechum', field2d_in)
   call mpas_pool_get_field(xana % subFields,            'spechum', field2d_out)
   field2d_out%array(:,:)=field2d_in%array(:,:)
   call mpas_pool_get_field(xctl % subFields,      'stream_function'  , field2d_in)
   call mpas_pool_get_field(xana % subFields,      'uReconstructZonal', field2d_out)
   field2d_out%array(:,:)=field2d_in%array(:,:)
   call mpas_pool_get_field(xctl % subFields, 'velocity_potential'    , field2d_in)
   call mpas_pool_get_field(xana % subFields, 'uReconstructMeridional', field2d_out)
   field2d_out%array(:,:)=field2d_in%array(:,:)
   call mpas_pool_get_field(xctl % subFields,   'surface_pressure', field1d_in)
   call mpas_pool_get_field(xana % subFields,   'surface_pressure', field1d_out)
   field1d_out%array(:)=field1d_in%array(:)

end subroutine mpas_linvarcha_c2a_multiplyinverseadjoint

!-------------------------------------------------------------------------------
! Adjoint code of subroutine mpas_reconstruct()
!                 in MPAS-Model/src/operators/mpas_vector_reconstruction.F
!
subroutine mpas_reconstruct_2dAD(meshPool, u, uReconstructX, uReconstructY, uReconstructZ,&
                                 uReconstructZonal, uReconstructMeridional, includeHalos)!{{{

    implicit none

    type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
    real (kind=RKIND), dimension(:,:), intent(inout) :: u !< Input: Velocity field on edges
    real (kind=RKIND), dimension(:,:), intent(inout) :: uReconstructX !< Output: X Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:,:), intent(inout) :: uReconstructY !< Output: Y Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:,:), intent(inout) :: uReconstructZ !< Output: Z Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:,:), intent(inout) :: uReconstructZonal !< Output: Zonal Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:,:), intent(inout) :: uReconstructMeridional !< Output: Meridional Component of velocity reconstructed to cell centers
    logical, optional, intent(in) :: includeHalos !< Input: Optional logical that allows reconstruction over halo regions

    !   temporary arrays needed in the compute procedure
    logical :: includeHalosLocal
    integer, pointer :: nCells
    integer, dimension(:,:), pointer :: edgesOnCell
    integer, dimension(:), pointer :: nEdgesOnCell
    integer :: iCell,iEdge, i
    real(kind=RKIND), dimension(:), pointer :: latCell, lonCell

    real (kind=RKIND), dimension(:,:,:), pointer :: coeffs_reconstruct

    logical, pointer :: on_a_sphere

    real (kind=RKIND) :: clat, slat, clon, slon

    if ( present(includeHalos) ) then
       includeHalosLocal = includeHalos
    else
       includeHalosLocal = .false.
    end if

    ! stored arrays used during compute procedure
    call mpas_pool_get_array(meshPool, 'coeffs_reconstruct', coeffs_reconstruct)

    ! temporary variables
    call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
    call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)

    if ( includeHalosLocal ) then
       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
    else
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCells)
    end if

    call mpas_pool_get_array(meshPool, 'latCell', latCell)
    call mpas_pool_get_array(meshPool, 'lonCell', lonCell)

    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

    if (on_a_sphere) then
      !$omp do schedule(runtime)
      do iCell = 1, nCells
        clat = cos(latCell(iCell))
        slat = sin(latCell(iCell))
        clon = cos(lonCell(iCell))
        slon = sin(lonCell(iCell))

        ! initialize the reconstructed vectors  !BJJ need to do something TODO
        uReconstructX(:,iCell) = 0.0
        uReconstructY(:,iCell) = 0.0
        uReconstructZ(:,iCell) = 0.0

        uReconstructX(:,iCell) = uReconstructX(:,iCell) - clon*slat * uReconstructMeridional(:,iCell)
        uReconstructY(:,iCell) = uReconstructY(:,iCell) - slon*slat * uReconstructMeridional(:,iCell)
        uReconstructZ(:,iCell) = uReconstructZ(:,iCell) + clat * uReconstructMeridional(:,iCell)

        uReconstructX(:,iCell) = uReconstructX(:,iCell) - slon * uReconstructZonal(:,iCell)
        uReconstructY(:,iCell) = uReconstructY(:,iCell) + clon * uReconstructZonal(:,iCell)
      end do
      !$omp end do
    else
      !$omp do schedule(runtime)
      do iCell = 1, nCells
        uReconstructY(:,iCell) = uReconstructY(:,iCell) + uReconstructMeridional(:,iCell)
        uReconstructX(:,iCell) = uReconstructX(:,iCell) + uReconstructZonal     (:,iCell)
      end do
      !$omp end do
    end if

    call mpas_threading_barrier()

    ! loop over cell centers
    !$omp do schedule(runtime)
    do iCell = 1, nCells
      ! a more efficient reconstruction where rbf_values*matrix_reconstruct
      ! has been precomputed in coeffs_reconstruct
      do i=1,nEdgesOnCell(iCell)
        iEdge = edgesOnCell(i,iCell)
        u(:,iEdge) = u(:,iEdge) + coeffs_reconstruct(1,i,iCell) * uReconstructX(:,iCell) &
                                + coeffs_reconstruct(2,i,iCell) * uReconstructY(:,iCell) &
                                + coeffs_reconstruct(3,i,iCell) * uReconstructZ(:,iCell)
      enddo
    enddo   ! iCell
    !$omp end do

end subroutine mpas_reconstruct_2dAD!}}}

!-------------------------------------------------------------------------------

subroutine mpas_reconstruct_1dAD(meshPool, u, uReconstructX, uReconstructY, uReconstructZ,&
                                 uReconstructZonal, uReconstructMeridional, includeHalos)!{{{

    implicit none

    type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
    real (kind=RKIND), dimension(:), intent(inout) :: u !< Input: Velocity field on edges
    real (kind=RKIND), dimension(:), intent(inout) :: uReconstructX !< Output: X Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:), intent(inout) :: uReconstructY !< Output: Y Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:), intent(inout) :: uReconstructZ !< Output: Z Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:), intent(inout) :: uReconstructZonal !< Output: Zonal Component of velocity reconstructed to cell centers
    real (kind=RKIND), dimension(:), intent(inout) :: uReconstructMeridional !< Output: Meridional Component of velocity reconstructed to cell centers
    logical, optional, intent(in) :: includeHalos !< Input: Optional logical that allows reconstruction over halo regions

    !   temporary arrays needed in the compute procedure
    logical :: includeHalosLocal
    integer, pointer :: nCells
    integer, dimension(:,:), pointer :: edgesOnCell
    integer, dimension(:), pointer :: nEdgesOnCell
    integer :: iCell,iEdge, i
    real(kind=RKIND), dimension(:), pointer :: latCell, lonCell

    real (kind=RKIND), dimension(:,:,:), pointer :: coeffs_reconstruct

    logical, pointer :: on_a_sphere

    real (kind=RKIND) :: clat, slat, clon, slon

    if ( present(includeHalos) ) then
       includeHalosLocal = includeHalos
    else
       includeHalosLocal = .false.
    end if

    ! stored arrays used during compute procedure
    call mpas_pool_get_array(meshPool, 'coeffs_reconstruct', coeffs_reconstruct)

    ! temporary variables
    call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
    call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)

    if ( includeHalosLocal ) then
       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
    else
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCells)
    end if

    call mpas_pool_get_array(meshPool, 'latCell', latCell)
    call mpas_pool_get_array(meshPool, 'lonCell', lonCell)

    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

    if (on_a_sphere) then
      !$omp do schedule(runtime)
      do iCell = 1, nCells
        clat = cos(latCell(iCell))
        slat = sin(latCell(iCell))
        clon = cos(lonCell(iCell))
        slon = sin(lonCell(iCell))

        ! initialize the reconstructed vectors  !BJJ need to do something TODO
        uReconstructX(iCell) = 0.0
        uReconstructY(iCell) = 0.0
        uReconstructZ(iCell) = 0.0

        uReconstructX(iCell) = uReconstructX(iCell) - clon*slat * uReconstructMeridional(iCell)
        uReconstructY(iCell) = uReconstructY(iCell) - slon*slat * uReconstructMeridional(iCell)
        uReconstructZ(iCell) = uReconstructZ(iCell) + clat * uReconstructMeridional(iCell)

        uReconstructX(iCell) = uReconstructX(iCell) - slon * uReconstructZonal(iCell)
        uReconstructY(iCell) = uReconstructY(iCell) + clon * uReconstructZonal(iCell)
      end do
      !$omp end do
    else
      !$omp do schedule(runtime)
      do iCell = 1, nCells
        uReconstructY(iCell) = uReconstructY(iCell) + uReconstructMeridional(iCell)
        uReconstructX(iCell) = uReconstructX(iCell) + uReconstructZonal     (iCell)
      end do
      !$omp end do
    end if

    call mpas_threading_barrier()

    ! loop over cell centers
    !$omp do schedule(runtime)
    do iCell = 1, nCells
      ! a more efficient reconstruction where rbf_values*matrix_reconstruct
      ! has been precomputed in coeffs_reconstruct
      do i=1,nEdgesOnCell(iCell)
        iEdge = edgesOnCell(i,iCell)
        u(iEdge) = u(iEdge) + coeffs_reconstruct(1,i,iCell) * uReconstructX(iCell) &
                                + coeffs_reconstruct(2,i,iCell) * uReconstructY(iCell) &
                                + coeffs_reconstruct(3,i,iCell) * uReconstructZ(iCell)
      enddo
    enddo   ! iCell
    !$omp end do

end subroutine mpas_reconstruct_1dAD!}}}

!-------------------------------------------------------------------------------
! input  : psi & chi @ cell center
! output : u & v @ cell center
subroutine psichi_to_uv(geom, ncells, psi, chi, u, v)

   implicit none
   type (mpas_geom),                         intent(in)  :: geom         !< geometry
   integer,                                  intent(in)  :: ncells
   real (kind=kind_real), dimension(ncells), intent(in)  :: psi, chi
   real (kind=kind_real), dimension(ncells), intent(out) :: u, v

   integer :: iC, iE, j
   real (kind=kind_real), dimension(geom%nCellsSolve) :: &
              psi_line_intg_dx, psi_line_intg_dy, chi_line_intg_dx, chi_line_intg_dy

   psi_line_intg_dx=MPAS_JEDI_ZERO_kr
   psi_line_intg_dy=MPAS_JEDI_ZERO_kr
   chi_line_intg_dx=MPAS_JEDI_ZERO_kr
   chi_line_intg_dy=MPAS_JEDI_ZERO_kr

   do iC = 1, geom%nCellsSolve
     do j = 1, geom%nEdgesOnCell(iC) ! or geom%maxEdges
       iE = geom%edgesOnCell(j,iC)
       psi_line_intg_dx(iC) = psi_line_intg_dx(iC) &
          + MPAS_JEDI_HALF_kr * ( psi(geom%cellsOnEdge(1,iE)) + psi(geom%cellsOnEdge(2,iE)) ) &
                              * geom%dvEdge(iE) * sin(MPAS_JEDI_ZERO_kr-geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC)
       psi_line_intg_dy(iC) = psi_line_intg_dy(iC) &
          + MPAS_JEDI_HALF_kr * ( psi(geom%cellsOnEdge(1,iE)) + psi(geom%cellsOnEdge(2,iE)) ) &
                              * geom%dvEdge(iE) * cos(geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC)
       chi_line_intg_dx(iC) = chi_line_intg_dx(iC) &
          + MPAS_JEDI_HALF_kr * ( chi(geom%cellsOnEdge(1,iE)) + chi(geom%cellsOnEdge(2,iE)) ) &
                              * geom%dvEdge(iE) * sin(MPAS_JEDI_ZERO_kr-geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC)
       chi_line_intg_dy(iC) = chi_line_intg_dy(iC) &
          + MPAS_JEDI_HALF_kr * ( chi(geom%cellsOnEdge(1,iE)) + chi(geom%cellsOnEdge(2,iE)) ) &
                              * geom%dvEdge(iE) * cos(geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC)
     enddo !- j
   enddo !- iC

   do iC=1, geom%nCellsSolve
     u(iC) = ( psi_line_intg_dx(iC) - chi_line_intg_dy(iC) ) / geom%areaCell(iC)
     v(iC) = ( psi_line_intg_dy(iC) + chi_line_intg_dx(iC) ) / geom%areaCell(iC)
   enddo

end subroutine psichi_to_uv

!-------------------------------------------------------------------------------

subroutine psichi_to_uvAD(geom, ncells, psi, chi, u, v)

   implicit none
   type (mpas_geom),                         intent(in)    :: geom         !< geometry
   integer,                                  intent(in)    :: ncells
   real (kind=kind_real), dimension(ncells), intent(inout) :: psi, chi
   real (kind=kind_real), dimension(ncells), intent(inout) :: u, v

   integer :: iC, iE, j
   real (kind=kind_real), dimension(geom%nCellsSolve) :: &
              psi_line_intg_dx, psi_line_intg_dy, chi_line_intg_dx, chi_line_intg_dy

   psi_line_intg_dx=MPAS_JEDI_ZERO_kr
   psi_line_intg_dy=MPAS_JEDI_ZERO_kr
   chi_line_intg_dx=MPAS_JEDI_ZERO_kr
   chi_line_intg_dy=MPAS_JEDI_ZERO_kr

   do iC=1, geom%nCellsSolve
     psi_line_intg_dx(iC) = psi_line_intg_dx(iC) + u(iC) / geom%areaCell(iC)
     chi_line_intg_dy(iC) = chi_line_intg_dy(iC) - u(iC) / geom%areaCell(iC)
     psi_line_intg_dy(iC) = psi_line_intg_dy(iC) + v(iC) / geom%areaCell(iC)
     chi_line_intg_dx(iC) = chi_line_intg_dx(iC) + v(iC) / geom%areaCell(iC)
   enddo

   do iC = 1, geom%nCellsSolve
     do j = 1, geom%nEdgesOnCell(iC) ! or geom%maxEdges
       iE = geom%edgesOnCell(j,iC)

       psi(geom%cellsOnEdge(1,iE)) = psi(geom%cellsOnEdge(1,iE)) + MPAS_JEDI_HALF_kr * geom%dvEdge(iE)&
 * sin(MPAS_JEDI_ZERO_kr-geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC) * psi_line_intg_dx(iC)
       psi(geom%cellsOnEdge(2,iE)) = psi(geom%cellsOnEdge(2,iE)) + MPAS_JEDI_HALF_kr * geom%dvEdge(iE)&
 * sin(MPAS_JEDI_ZERO_kr-geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC) * psi_line_intg_dx(iC)
       
       psi(geom%cellsOnEdge(1,iE)) = psi(geom%cellsOnEdge(1,iE)) + MPAS_JEDI_HALF_kr * geom%dvEdge(iE)&
 * cos(geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC) * psi_line_intg_dy(iC)
       psi(geom%cellsOnEdge(2,iE)) = psi(geom%cellsOnEdge(2,iE)) + MPAS_JEDI_HALF_kr * geom%dvEdge(iE)&
 * cos(geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC) * psi_line_intg_dy(iC)

       chi(geom%cellsOnEdge(1,iE)) = chi(geom%cellsOnEdge(1,iE)) + MPAS_JEDI_HALF_kr * geom%dvEdge(iE)&
 * sin(MPAS_JEDI_ZERO_kr-geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC) * chi_line_intg_dx(iC)
       chi(geom%cellsOnEdge(2,iE)) = chi(geom%cellsOnEdge(2,iE)) + MPAS_JEDI_HALF_kr * geom%dvEdge(iE)&
 * sin(MPAS_JEDI_ZERO_kr-geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC) * chi_line_intg_dx(iC)

       chi(geom%cellsOnEdge(1,iE)) = chi(geom%cellsOnEdge(1,iE)) + MPAS_JEDI_HALF_kr * geom%dvEdge(iE)&
 * cos(geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC) * chi_line_intg_dy(iC)
       chi(geom%cellsOnEdge(2,iE)) = chi(geom%cellsOnEdge(2,iE)) + MPAS_JEDI_HALF_kr * geom%dvEdge(iE)&
 * cos(geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC) * chi_line_intg_dy(iC)
     enddo !- j
   enddo !- iC

end subroutine psichi_to_uvAD

!-------------------------------------------------------------------------------
! input  : psi & chi @ cell center
! output : u & v @ cell center
subroutine psichi_to_uv2(geom, psi, chi, u, v)

   use mpas_vector_reconstruction

   implicit none
   type (mpas_geom),                              intent(in)  :: geom         !< geometry
   real (kind=kind_real), dimension(geom%nCells), intent(in)  :: psi, chi
   real (kind=kind_real), dimension(geom%nCells), intent(out) :: u, v

   integer :: iC, iE, iV, j
   real (kind=kind_real), dimension(:), allocatable :: psi_v, edge_wind, &
                                        uReconstructX, uReconstructY, uReconstructZ
   type (mpas_pool_type), pointer :: mesh

   allocate(psi_v(geom%nVertices))
   allocate(edge_wind(geom%nEdges))
   allocate(uReconstructX(geom%nCells))
   allocate(uReconstructY(geom%nCells))
   allocate(uReconstructZ(geom%nCells))

   psi_v=MPAS_JEDI_ZERO_kr
   edge_wind=MPAS_JEDI_ZERO_kr

   ! Interpolate psi/chi in cell center to vertice
   do iV = 1, geom%nVerticesSolve !nVertices
     do j = 1, geom%vertexDegree
       iC = geom%cellsOnVertex(j,iV)
       psi_v(iV) = psi_v(iV) + geom%kiteAreasOnVertex(j,iV) * psi(iC)
     enddo
     psi_v(iV) = psi_v(iV) / geom%areaTriangle(iV)
   enddo

   !Halo exchange @ vertices ???

   !get edge_wind
   do iE = 1, geom%nEdgesSolve !nEdges
      edge_wind(iE) = edge_wind(iE) - &
       ( chi(geom%cellsOnEdge(2,iE)) - chi(geom%cellsOnEdge(1,iE)) ) / geom%dcEdge(iE) - &
       ( psi_v(geom%verticesOnEdge(2,iE)) - psi_v(geom%verticesOnEdge(1,iE)) ) / geom%dvEdge(iE)
   enddo

   !Halo exchange @ edge ???

   !edge_wind -> u/v @ cell center
   call mpas_pool_get_subpool( geom % domain % blocklist % structs, 'mesh', mesh)

   call mpas_reconstruct(mesh, edge_wind,           &
                         uReconstructX,             &
                         uReconstructY,             &
                         uReconstructZ,             &
                         u,                         &
                         v, .True.                  &
                        )

   deallocate(psi_v, edge_wind, uReconstructX, uReconstructY, uReconstructZ)

end subroutine psichi_to_uv2

!-------------------------------------------------------------------------------

subroutine psichi_to_uv2AD(geom, psi ,chi, u, v)

   implicit none
   type (mpas_geom),                              intent(in)    :: geom         !< geometry
   real (kind=kind_real), dimension(geom%nCells), intent(inout) :: psi, chi
   real (kind=kind_real), dimension(geom%nCells), intent(inout) :: u, v

   integer :: iC, iE, iV, j
   real (kind=kind_real), dimension(:), allocatable :: psi_v, edge_wind, &
                                        uReconstructX, uReconstructY, uReconstructZ
   type (mpas_pool_type), pointer :: mesh

   allocate(psi_v(geom%nVertices))
   allocate(edge_wind(geom%nEdges))
   allocate(uReconstructX(geom%nCells))
   allocate(uReconstructY(geom%nCells))
   allocate(uReconstructZ(geom%nCells))

   psi_v=MPAS_JEDI_ZERO_kr
   edge_wind=MPAS_JEDI_ZERO_kr

   !-
   call mpas_pool_get_subpool( geom % domain % blocklist % structs, 'mesh', mesh)
   call mpas_reconstruct_1dAD(mesh, edge_wind,      &
                         uReconstructX,             &
                         uReconstructY,             &
                         uReconstructZ,             &
                         u,                         &
                         v, .True.                  &
                        )

   !-
   do iE = 1, geom%nEdgesSolve !nEdges
      chi(geom%cellsOnEdge(2,iE))      = chi(geom%cellsOnEdge(2,iE))      - edge_wind(iE) / geom%dcEdge(iE)
      chi(geom%cellsOnEdge(1,iE))      = chi(geom%cellsOnEdge(1,iE))      + edge_wind(iE) / geom%dcEdge(iE)
      psi_v(geom%verticesOnEdge(2,iE)) = psi_v(geom%verticesOnEdge(2,iE)) - edge_wind(iE) / geom%dvEdge(iE)
      psi_v(geom%verticesOnEdge(1,iE)) = psi_v(geom%verticesOnEdge(1,iE)) + edge_wind(iE) / geom%dvEdge(iE)
   enddo

   !-
   do iV = 1, geom%nVerticesSolve !nVertices
     psi_v(iV) = psi_v(iV) / geom%areaTriangle(iV)
     do j = 1, geom%vertexDegree
       iC = geom%cellsOnVertex(j,iV)
       psi(iC) = psi(iC) + geom%kiteAreasOnVertex(j,iV) * psi_v(iV)
     enddo
   enddo

   deallocate(psi_v, edge_wind, uReconstructX, uReconstructY, uReconstructZ)

end subroutine psichi_to_uv2AD

! ------------------------------------------------------------------------------

end module mpas_linvarcha_c2a_mod
