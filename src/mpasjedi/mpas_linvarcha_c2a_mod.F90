! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_linvarcha_c2a_mod

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only : fckit_log
use iso_c_binding
use kinds

!oops
use oops_variables_mod, only: oops_variables

!mpas-jedi
use mpas_constants_mod
use mpas_fields_mod, only: mpas_fields
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
  type (mpas_pool_type), pointer :: trajectories => null()
end type mpas_linvarcha_c2a

integer, parameter    :: max_string=800
character(max_string) :: message
character(len=1024) :: buf

#define LISTED_TYPE mpas_linvarcha_c2a

!> Linked list interface - defines registry_t type
#include <oops/util/linkedList_i.f>

!> Global registry
type(registry_t) :: mpas_linvarcha_c2a_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include <oops/util/linkedList_c.f>
! ------------------------------------------------------------------------------

subroutine mpas_linvarcha_c2a_setup(self, bg, fg, geom, f_conf, vars)

   implicit none
   type(mpas_linvarcha_c2a),  intent(inout) :: self    !< Change variable structure
   type(mpas_fields), target,  intent(in)    :: bg
   type(mpas_fields), target,  intent(in)    :: fg
   type(mpas_geom),           intent(in)    :: geom
   type(fckit_configuration), intent(in)    :: f_conf  !< Configuration
   type(oops_variables),      intent(in)    :: vars

   type (field2DReal), pointer :: fld2d_t, fld2d_p, fld2d_qs

   integer :: ngrid

   if ( vars % has ('relhum') ) then
      !-- set trajectories for linear variable change
      ngrid = geom % nCells ! local + halo

      call mpas_pool_create_pool(self % trajectories)

      if ( .not. fg % has ('temperature') .or. .not. fg % has ('pressure') ) &
         call abor1_ftn("LinVarChaC2AMPAS::LinVarChaC2AMPAS, mpas_linvarcha_c2a_setup: Trajectory failed")

      call mpas_pool_get_field(fg % subFields, 'temperature', fld2d_t)
      call mpas_pool_get_field(fg % subFields, 'pressure'   , fld2d_p)

      call mpas_duplicate_field(fld2d_t, fld2d_qs) ! for saturation specific humidity

      call da_tp_to_qs( fld2d_t % array(:,1:ngrid), fld2d_p % array(:,1:ngrid), fld2d_qs% array(:,1:ngrid))

      call mpas_pool_add_field(self % trajectories, 'spechum_sat', fld2d_qs)
   end if

end subroutine mpas_linvarcha_c2a_setup

! ------------------------------------------------------------------------------

subroutine mpas_linvarcha_c2a_delete(self)

   implicit none
   type(mpas_linvarcha_c2a), intent(inout) :: self

   if (associated(self % trajectories)) then
      call mpas_pool_destroy_pool(self % trajectories)
   end if

end subroutine mpas_linvarcha_c2a_delete

! ------------------------------------------------------------------------------

subroutine mpas_linvarcha_c2a_multiply(self,xctl,xana)

   implicit none
   type(mpas_linvarcha_c2a), intent(inout) :: self
   type(mpas_fields), intent(inout) :: xctl
   type(mpas_fields), intent(inout) :: xana

   type (mpas_pool_iterator_type) :: poolItr
   type (field2DReal), pointer :: field2d_sf, field2d_vp, field2d_uRz, field2d_uRm
   type (field2DReal), pointer :: field2d_sh, field2d_rh, field2d_traj_qs
   type (field2DReal), pointer :: field2d_v_src, field2d_sf_v, field2d_e_src, field2d_u

   type (field2DReal), pointer :: field2d_ctl, field2d_ana
   type (field1DReal), pointer :: field1d_ctl, field1d_ana
   integer :: k, ngrid, nCells, nVertices, nEdges

   write(message,*) "DEBUG: mpas_linvarcha_c2a_multiply: xana % fldnames(:) =",xana % fldnames(:)
   call fckit_log%debug(message)
   write(message,*) "DEBUG: mpas_linvarcha_c2a_multiply: xctl % fldnames(:) =",xctl % fldnames(:)
   call fckit_log%debug(message)

   ngrid = xctl % geom % nCellsSolve ! local 

   ! no variable change
   call mpas_pool_begin_iteration(xctl%subFields)
   do while ( mpas_pool_get_next_member(xctl%subFields, poolItr) )
      ! NOTE: only 1d or 2d real fields are handled.
      if (poolItr % memberType == MPAS_POOL_FIELD .and. poolItr % dataType == MPAS_POOL_REAL) then 
         if( xana%has( trim(poolItr % memberName) ) ) then ! check if the variable names are the same.
            if (poolItr % nDims == 1) then
               call mpas_pool_get_field(xctl%subFields, trim(poolItr % memberName), field1d_ctl)
               call mpas_pool_get_field(xana%subFields, trim(poolItr % memberName), field1d_ana)
               field1d_ana%array(1:ngrid)=field1d_ctl%array(1:ngrid)
            else if (poolItr % nDims == 2) then
               call mpas_pool_get_field(xctl%subFields, trim(poolItr % memberName), field2d_ctl)
               call mpas_pool_get_field(xana%subFields, trim(poolItr % memberName), field2d_ana)
               field2d_ana%array(:,1:ngrid)=field2d_ctl%array(:,1:ngrid)
            else
               write(buf,*) '--> mpas_linvarcha_c2a_multiply: poolItr % nDims == ',poolItr % nDims,' not handled'
               call abor1_ftn(buf)
            end if
         end if
      end if
   end do !- end of pool iteration

   if( xctl%has('relhum') .and. xana%has('spechum') ) then
      call mpas_pool_get_field(xctl % subFields,             'relhum', field2d_rh)
      call mpas_pool_get_field(xana % subFields,            'spechum', field2d_sh)
      call mpas_pool_get_field(self % trajectories,     'spechum_sat', field2d_traj_qs)
      call pseudorh_to_spechum( field2d_rh % array(:,1:ngrid), field2d_sh % array(:,1:ngrid), &
                                field2d_traj_qs % array(:,1:ngrid) )
   end if

   if( xctl%has('stream_function') .and. xctl%has('velocity_potential') .and. &
       xana%has('uReconstructZonal') .and. xana%has('uReconstructMeridional') ) then
      call mpas_pool_get_field(xctl % subFields,        'stream_function', field2d_sf)
      call mpas_pool_get_field(xctl % subFields,     'velocity_potential', field2d_vp)
      call mpas_pool_get_field(xana % subFields,      'uReconstructZonal', field2d_uRz)
      call mpas_pool_get_field(xana % subFields, 'uReconstructMeridional', field2d_uRm)

      call mpas_dmpar_exch_halo_field(field2d_sf)
      call mpas_dmpar_exch_halo_field(field2d_vp)

      nCells = xctl % geom % nCells    ! local + halo
      nVertices = xctl % geom % nVertices ! local + halo
      nEdges = xctl % geom % nEdges    ! local + halo

      ! duplicate two temporary working spaces
      call mpas_pool_get_field( xctl % geom % domain % blocklist % allFields, 'vorticity', field2d_v_src)
      call mpas_duplicate_field(field2d_v_src, field2d_sf_v)
      field2d_sf_v % fieldName = 'stream_function at vertices'
      call mpas_pool_get_field( xctl % geom % domain % blocklist % allFields, 'u', field2d_e_src)
      call mpas_duplicate_field(field2d_e_src, field2d_u)
      field2d_u % fieldName = 'Horizontal normal velocity at edges'

      call psichi_to_uv_edge_step1( xctl % geom, field2d_sf%array(:,1:nCells), field2d_sf_v%array(:,1:nVertices) )
      call mpas_dmpar_exch_halo_field(field2d_sf_v)
      call psichi_to_uv_edge_step2( xctl % geom, field2d_sf_v%array(:,1:nVertices), field2d_vp%array(:,1:nCells), &
                                   field2d_u%array(:,1:nEdges) )
      call mpas_dmpar_exch_halo_field(field2d_u)
      call psichi_to_uv_edge_step3( xctl % geom, field2d_u%array(:,1:nEdges), &
                                   field2d_uRz%array(:,1:nCells), field2d_uRm%array(:,1:nCells) )

      call mpas_deallocate_field(field2d_sf_v)
      call mpas_deallocate_field(field2d_u)
   end if

end subroutine mpas_linvarcha_c2a_multiply

! ------------------------------------------------------------------------------

subroutine mpas_linvarcha_c2a_multiplyadjoint(self,xana,xctl)

   implicit none
   type(mpas_linvarcha_c2a), intent(inout) :: self
   type(mpas_fields), intent(inout) :: xana
   type(mpas_fields), intent(inout) :: xctl

   type (mpas_pool_iterator_type) :: poolItr
   type (field2DReal), pointer :: field2d_sf, field2d_vp, field2d_uRz, field2d_uRm
   type (field2DReal), pointer :: field2d_rh, field2d_sh, field2d_traj_qs
   type (field2DReal), pointer :: field2d_v_src, field2d_sf_v, field2d_e_src, field2d_u

   type (field2DReal), pointer :: field2d_ctl, field2d_ana
   type (field1DReal), pointer :: field1d_ctl, field1d_ana
   integer :: k, ngrid, nCells, nVertices, nEdges

   write(message,*) "DEBUG: mpas_linvarcha_c2a_multiplyadjoint: xana % fldnames(:) =",xana % fldnames(:)
   call fckit_log%debug(message)
   write(message,*) "DEBUG: mpas_linvarcha_c2a_multiplyadjoint: xctl % fldnames(:) =",xctl % fldnames(:)
   call fckit_log%debug(message)

   ngrid = xctl % geom % nCellsSolve

   ! no variable change
   call mpas_pool_begin_iteration(xctl%subFields)
   do while ( mpas_pool_get_next_member(xctl%subFields, poolItr) )
      ! NOTE: only 1d or 2d real fields are handled.
      if (poolItr % memberType == MPAS_POOL_FIELD .and. poolItr % dataType == MPAS_POOL_REAL) then
         if( xana%has( trim(poolItr % memberName) ) ) then ! check if the variable names are the same.
            if (poolItr % nDims == 1) then
               call mpas_pool_get_field(xctl%subFields, trim(poolItr % memberName), field1d_ctl)
               call mpas_pool_get_field(xana%subFields, trim(poolItr % memberName), field1d_ana)
               field1d_ctl%array(1:ngrid)=field1d_ctl%array(1:ngrid)+field1d_ana%array(1:ngrid)
               field1d_ana%array(1:ngrid)=MPAS_JEDI_ZERO_kr
            else if (poolItr % nDims == 2) then
               call mpas_pool_get_field(xctl%subFields, trim(poolItr % memberName), field2d_ctl)
               call mpas_pool_get_field(xana%subFields, trim(poolItr % memberName), field2d_ana)
               field2d_ctl%array(:,1:ngrid)=field2d_ctl%array(:,1:ngrid)+field2d_ana%array(:,1:ngrid)
               field2d_ana%array(:,1:ngrid)=MPAS_JEDI_ZERO_kr
            else
               write(buf,*) '--> mpas_linvarcha_c2a_multiplyadjoint: poolItr % nDims == ',poolItr % nDims,' not handled'
               call abor1_ftn(buf)
            end if
         end if
      end if
   end do !- end of pool iteration

   if( xctl%has('relhum') .and. xana%has('spechum') ) then
      call mpas_pool_get_field(xctl % subFields,             'relhum', field2d_rh)
      call mpas_pool_get_field(xana % subFields,            'spechum', field2d_sh)
      call mpas_pool_get_field(self % trajectories,     'spechum_sat', field2d_traj_qs)
      call pseudorh_to_spechumAD( field2d_rh % array(:,1:ngrid), field2d_sh % array(:,1:ngrid), &
                                  field2d_traj_qs % array(:,1:ngrid) )
   end if

   if( xctl%has('stream_function') .and. xctl%has('velocity_potential') .and. &
       xana%has('uReconstructZonal') .and. xana%has('uReconstructMeridional') ) then
      call mpas_pool_get_field(xctl % subFields,        'stream_function', field2d_sf)
      call mpas_pool_get_field(xctl % subFields,     'velocity_potential', field2d_vp)
      call mpas_pool_get_field(xana % subFields,      'uReconstructZonal', field2d_uRz)
      call mpas_pool_get_field(xana % subFields, 'uReconstructMeridional', field2d_uRm)

      nCells = xctl % geom % nCells    ! local + halo
      nVertices = xctl % geom % nVertices ! local + halo
      nEdges = xctl % geom % nEdges    ! local + halo

      ! duplicate two temporary working spaces
      call mpas_pool_get_field( xctl % geom % domain % blocklist % allFields, 'vorticity', field2d_v_src)
      call mpas_duplicate_field(field2d_v_src, field2d_sf_v)
      field2d_sf_v % fieldName = 'stream_function at vertices'
      call mpas_pool_get_field( xctl % geom % domain % blocklist % allFields, 'u', field2d_e_src)
      call mpas_duplicate_field(field2d_e_src, field2d_u)
      field2d_u % fieldName = 'Horizontal normal velocity at edges'

      call psichi_to_uv_edge_step3AD( xctl % geom, field2d_u%array(:,1:nEdges), &
                                     field2d_uRz%array(:,1:nCells), field2d_uRm%array(:,1:nCells) )
      call mpas_dmpar_exch_halo_adj_field(field2d_u)
      call psichi_to_uv_edge_step2AD( xctl % geom, field2d_sf_v%array(:,1:nVertices), field2d_vp%array(:,1:nCells), &
                                     field2d_u%array(:,1:nEdges) )
      call mpas_dmpar_exch_halo_adj_field(field2d_sf_v)
      call psichi_to_uv_edge_step1AD( xctl % geom, field2d_sf%array(:,1:nCells), field2d_sf_v%array(:,1:nVertices) )

      call mpas_deallocate_field(field2d_sf_v)
      call mpas_deallocate_field(field2d_u)

      field2d_uRz%array(:,1:xana%geom%nCellsSolve)=MPAS_JEDI_ZERO_kr
      field2d_uRm%array(:,1:xana%geom%nCellsSolve)=MPAS_JEDI_ZERO_kr

      call mpas_dmpar_exch_halo_adj_field(field2d_sf)
      call mpas_dmpar_exch_halo_adj_field(field2d_vp)
   end if

end subroutine mpas_linvarcha_c2a_multiplyadjoint

! ------------------------------------------------------------------------------

subroutine mpas_linvarcha_c2a_multiplyinverse(self,xana,xctl)

   implicit none
   type(mpas_linvarcha_c2a), intent(inout) :: self
   type(mpas_fields), intent(inout) :: xana
   type(mpas_fields), intent(inout) :: xctl

   type (mpas_pool_iterator_type) :: poolItr
   type (field2DReal), pointer :: field2d_rh, field2d_sh, field2d_traj_qs

   type (field2DReal), pointer :: field2d_ctl, field2d_ana
   type (field1DReal), pointer :: field1d_ctl, field1d_ana
   integer :: ngrid

   write(message,*) "DEBUG: mpas_linvarcha_c2a_multiplyinverse: xana % fldnames(:) =",xana % fldnames(:)
   call fckit_log%debug(message)
   write(message,*) "DEBUG: mpas_linvarcha_c2a_multiplyinverse: xctl % fldnames(:) =",xctl % fldnames(:)
   call fckit_log%debug(message)

   ngrid = xctl % geom % nCellsSolve ! local

   ! no variable change
   call mpas_pool_begin_iteration(xctl%subFields)
   do while ( mpas_pool_get_next_member(xctl%subFields, poolItr) )
      ! only 1d or 2d real fields are handled.
      if (poolItr % memberType == MPAS_POOL_FIELD .and. poolItr % dataType == MPAS_POOL_REAL) then
         if( xana%has( trim(poolItr % memberName) ) ) then ! check if the variable names are the same.
            if (poolItr % nDims == 1) then
               call mpas_pool_get_field(xctl%subFields, trim(poolItr % memberName), field1d_ctl)
               call mpas_pool_get_field(xana%subFields, trim(poolItr % memberName), field1d_ana)
               field1d_ctl%array(1:ngrid)=field1d_ana%array(1:ngrid)
            else if (poolItr % nDims == 2) then
               call mpas_pool_get_field(xctl%subFields, trim(poolItr % memberName), field2d_ctl)
               call mpas_pool_get_field(xana%subFields, trim(poolItr % memberName), field2d_ana)
               field2d_ctl%array(:,1:ngrid)=field2d_ana%array(:,1:ngrid)
            else
               write(buf,*) '--> mpas_linvarcha_c2a_multiplyinverse: poolItr % nDims == ',poolItr % nDims,' not handled'
               call abor1_ftn(buf)
            end if
         end if
      end if
   end do !- end of pool iteration

   if( xctl%has('relhum') .and. xana%has('spechum') ) then
      call mpas_pool_get_field(xana % subFields,            'spechum', field2d_sh)
      call mpas_pool_get_field(xctl % subFields,             'relhum', field2d_rh)
      call mpas_pool_get_field(self % trajectories,     'spechum_sat', field2d_traj_qs)
      call pseudorh_to_spechum_inverse( field2d_rh % array(:,1:ngrid), field2d_sh % array(:,1:ngrid), &
                                        field2d_traj_qs % array(:,1:ngrid) )
   end if

   !-- dummy inverse operator
   if( xctl%has('stream_function') .and. xctl%has('velocity_potential') .and. &
       xana%has('uReconstructZonal') .and. xana%has('uReconstructMeridional') ) then
      call mpas_pool_get_field(xana % subFields,      'uReconstructZonal', field2d_ana)
      call mpas_pool_get_field(xctl % subFields,      'stream_function'  , field2d_ctl)
      field2d_ctl%array(:,1:ngrid)=field2d_ana%array(:,1:ngrid)

      call mpas_pool_get_field(xana % subFields, 'uReconstructMeridional', field2d_ana)
      call mpas_pool_get_field(xctl % subFields, 'velocity_potential'    , field2d_ctl)
      field2d_ctl%array(:,1:ngrid)=field2d_ana%array(:,1:ngrid)
   end if

end subroutine mpas_linvarcha_c2a_multiplyinverse

! ------------------------------------------------------------------------------

subroutine mpas_linvarcha_c2a_multiplyinverseadjoint(self,xctl,xana)

   implicit none
   type(mpas_linvarcha_c2a), intent(inout) :: self
   type(mpas_fields), intent(inout) :: xctl
   type(mpas_fields), intent(inout) :: xana

   type (mpas_pool_iterator_type) :: poolItr
   type (field2DReal), pointer :: field2d_rh, field2d_sh, field2d_traj_qs

   type (field2DReal), pointer :: field2d_ctl, field2d_ana
   type (field1DReal), pointer :: field1d_ctl, field1d_ana
   integer :: ngrid

   write(message,*) "DEBUG: mpas_linvarcha_c2a_multiplyinverseadjoint: xana % fldnames(:) =",xana % fldnames(:)
   call fckit_log%debug(message)
   write(message,*) "DEBUG: mpas_linvarcha_c2a_multiplyinverseadjoint: xctl % fldnames(:) =",xctl % fldnames(:)
   call fckit_log%debug(message)

   ngrid = xctl % geom % nCellsSolve ! local

   ! no variable change
   call mpas_pool_begin_iteration(xctl%subFields)
   do while ( mpas_pool_get_next_member(xctl%subFields, poolItr) )
      ! only 1d or 2d real fields are handled.
      if (poolItr % memberType == MPAS_POOL_FIELD .and. poolItr % dataType == MPAS_POOL_REAL) then
         if( xana%has( trim(poolItr % memberName) ) ) then ! check if the variable names are the same.
            if (poolItr % nDims == 1) then
               call mpas_pool_get_field(xctl%subFields, trim(poolItr % memberName), field1d_ctl)
               call mpas_pool_get_field(xana%subFields, trim(poolItr % memberName), field1d_ana)
               field1d_ana%array(1:ngrid)=field1d_ana%array(1:ngrid)+field1d_ctl%array(1:ngrid)
               field1d_ctl%array(1:ngrid)=MPAS_JEDI_ZERO_kr
            else if (poolItr % nDims == 2) then
               call mpas_pool_get_field(xctl%subFields, trim(poolItr % memberName), field2d_ctl)
               call mpas_pool_get_field(xana%subFields, trim(poolItr % memberName), field2d_ana)
               field2d_ana%array(:,1:ngrid)=field2d_ana%array(:,1:ngrid)+field2d_ctl%array(:,1:ngrid)
               field2d_ctl%array(:,1:ngrid)=MPAS_JEDI_ZERO_kr
            else
               write(buf,*) '--> mpas_linvarcha_c2a_multiplyinverse: poolItr % nDims == ',poolItr % nDims,' not handled'
               call abor1_ftn(buf)
            end if
         end if
      end if
   end do !- end of pool iteration

   if( xctl%has('relhum') .and. xana%has('spechum') ) then
      call mpas_pool_get_field(xana % subFields,            'spechum', field2d_sh)
      call mpas_pool_get_field(xctl % subFields,             'relhum', field2d_rh)
      call mpas_pool_get_field(self % trajectories,     'spechum_sat', field2d_traj_qs)
      call pseudorh_to_spechum_inverseAD( field2d_rh % array(:,1:ngrid), field2d_sh % array(:,1:ngrid), &
                                          field2d_traj_qs % array(:,1:ngrid) )
   end if

   !-- dummy inverseAD operator
   if( xctl%has('stream_function') .and. xctl%has('velocity_potential') .and. &
       xana%has('uReconstructZonal') .and. xana%has('uReconstructMeridional') ) then
      call mpas_pool_get_field(xana % subFields,      'uReconstructZonal', field2d_ana)
      call mpas_pool_get_field(xctl % subFields,      'stream_function'  , field2d_ctl)
      field2d_ana%array(:,1:ngrid)=field2d_ana%array(:,1:ngrid)+field2d_ctl%array(:,1:ngrid)
      field2d_ctl%array(:,1:ngrid)=MPAS_JEDI_ZERO_kr

      call mpas_pool_get_field(xana % subFields, 'uReconstructMeridional', field2d_ana)
      call mpas_pool_get_field(xctl % subFields, 'velocity_potential'    , field2d_ctl)
      field2d_ana%array(:,1:ngrid)=field2d_ana%array(:,1:ngrid)+field2d_ctl%array(:,1:ngrid)
      field2d_ctl%array(:,1:ngrid)=MPAS_JEDI_ZERO_kr
   end if

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

    ! initialize the reconstructed vectors
    uReconstructX = MPAS_JEDI_ZERO_kr
    uReconstructY = MPAS_JEDI_ZERO_kr
    uReconstructZ = MPAS_JEDI_ZERO_kr
    if (on_a_sphere) then
      do iCell = 1, nCells
        clat = cos(latCell(iCell))
        slat = sin(latCell(iCell))
        clon = cos(lonCell(iCell))
        slon = sin(lonCell(iCell))

        uReconstructX(:,iCell) = uReconstructX(:,iCell) - clon*slat * uReconstructMeridional(:,iCell)
        uReconstructY(:,iCell) = uReconstructY(:,iCell) - slon*slat * uReconstructMeridional(:,iCell)
        uReconstructZ(:,iCell) = uReconstructZ(:,iCell) + clat * uReconstructMeridional(:,iCell)

        uReconstructX(:,iCell) = uReconstructX(:,iCell) - slon * uReconstructZonal(:,iCell)
        uReconstructY(:,iCell) = uReconstructY(:,iCell) + clon * uReconstructZonal(:,iCell)
      end do
    else
      do iCell = 1, nCells
        uReconstructY(:,iCell) = uReconstructY(:,iCell) + uReconstructMeridional(:,iCell)
        uReconstructX(:,iCell) = uReconstructX(:,iCell) + uReconstructZonal     (:,iCell)
      end do
    end if

    ! loop over cell centers
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

    ! initialize the reconstructed vectors
    uReconstructX = MPAS_JEDI_ZERO_kr
    uReconstructY = MPAS_JEDI_ZERO_kr
    uReconstructZ = MPAS_JEDI_ZERO_kr
    if (on_a_sphere) then
      do iCell = 1, nCells
        clat = cos(latCell(iCell))
        slat = sin(latCell(iCell))
        clon = cos(lonCell(iCell))
        slon = sin(lonCell(iCell))

        uReconstructX(iCell) = uReconstructX(iCell) - clon*slat * uReconstructMeridional(iCell)
        uReconstructY(iCell) = uReconstructY(iCell) - slon*slat * uReconstructMeridional(iCell)
        uReconstructZ(iCell) = uReconstructZ(iCell) + clat * uReconstructMeridional(iCell)

        uReconstructX(iCell) = uReconstructX(iCell) - slon * uReconstructZonal(iCell)
        uReconstructY(iCell) = uReconstructY(iCell) + clon * uReconstructZonal(iCell)
      end do
    else
      do iCell = 1, nCells
        uReconstructY(iCell) = uReconstructY(iCell) + uReconstructMeridional(iCell)
        uReconstructX(iCell) = uReconstructX(iCell) + uReconstructZonal     (iCell)
      end do
    end if

    ! loop over cell centers
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

end subroutine mpas_reconstruct_1dAD!}}}

!-------------------------------------------------------------------------------
! input  : psi & chi @ cell center
! output : u & v @ cell center
subroutine psichi_to_uv_center(geom, psi, chi, u, v)

   implicit none
   type (mpas_geom),                                               intent(in)  :: geom         !< geometry
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nCells), intent(in)  :: psi, chi
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nCells), intent(out) :: u, v

   integer :: iC, iE, j
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nCellsSolve) :: &
              psi_line_intg_dx, psi_line_intg_dy, chi_line_intg_dx, chi_line_intg_dy

   psi_line_intg_dx=MPAS_JEDI_ZERO_kr
   psi_line_intg_dy=MPAS_JEDI_ZERO_kr
   chi_line_intg_dx=MPAS_JEDI_ZERO_kr
   chi_line_intg_dy=MPAS_JEDI_ZERO_kr

   do iC = 1, geom%nCellsSolve
     do j = 1, geom%nEdgesOnCell(iC) ! or geom%maxEdges
       iE = geom%edgesOnCell(j,iC)
       psi_line_intg_dx(:,iC) = psi_line_intg_dx(:,iC) &
          + MPAS_JEDI_HALF_kr * ( psi(:,geom%cellsOnEdge(1,iE)) + psi(:,geom%cellsOnEdge(2,iE)) ) &
                              * geom%dvEdge(iE) * sin(MPAS_JEDI_ZERO_kr-geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC)
       psi_line_intg_dy(:,iC) = psi_line_intg_dy(:,iC) &
          + MPAS_JEDI_HALF_kr * ( psi(:,geom%cellsOnEdge(1,iE)) + psi(:,geom%cellsOnEdge(2,iE)) ) &
                              * geom%dvEdge(iE) * cos(geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC)
       chi_line_intg_dx(:,iC) = chi_line_intg_dx(:,iC) &
          + MPAS_JEDI_HALF_kr * ( chi(:,geom%cellsOnEdge(1,iE)) + chi(:,geom%cellsOnEdge(2,iE)) ) &
                              * geom%dvEdge(iE) * sin(MPAS_JEDI_ZERO_kr-geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC)
       chi_line_intg_dy(:,iC) = chi_line_intg_dy(:,iC) &
          + MPAS_JEDI_HALF_kr * ( chi(:,geom%cellsOnEdge(1,iE)) + chi(:,geom%cellsOnEdge(2,iE)) ) &
                              * geom%dvEdge(iE) * cos(geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC)
     enddo !- j
   enddo !- iC

   do iC=1, geom%nCellsSolve
     u(:,iC) = ( psi_line_intg_dx(:,iC) - chi_line_intg_dy(:,iC) ) / geom%areaCell(iC)
     v(:,iC) = ( psi_line_intg_dy(:,iC) + chi_line_intg_dx(:,iC) ) / geom%areaCell(iC)
   enddo

end subroutine psichi_to_uv_center

!-------------------------------------------------------------------------------

subroutine psichi_to_uv_centerAD(geom, psi, chi, u, v)

   implicit none
   type (mpas_geom),                                               intent(in)    :: geom         !< geometry
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nCells), intent(inout) :: psi, chi
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nCells), intent(inout) :: u, v

   integer :: iC, iE, j
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nCellsSolve) :: &
              psi_line_intg_dx, psi_line_intg_dy, chi_line_intg_dx, chi_line_intg_dy

   psi_line_intg_dx=MPAS_JEDI_ZERO_kr
   psi_line_intg_dy=MPAS_JEDI_ZERO_kr
   chi_line_intg_dx=MPAS_JEDI_ZERO_kr
   chi_line_intg_dy=MPAS_JEDI_ZERO_kr

   do iC=1, geom%nCellsSolve
     psi_line_intg_dx(:,iC) = psi_line_intg_dx(:,iC) + u(:,iC) / geom%areaCell(iC)
     chi_line_intg_dy(:,iC) = chi_line_intg_dy(:,iC) - u(:,iC) / geom%areaCell(iC)
     psi_line_intg_dy(:,iC) = psi_line_intg_dy(:,iC) + v(:,iC) / geom%areaCell(iC)
     chi_line_intg_dx(:,iC) = chi_line_intg_dx(:,iC) + v(:,iC) / geom%areaCell(iC)
   enddo

   do iC = 1, geom%nCellsSolve
     do j = 1, geom%nEdgesOnCell(iC) ! or geom%maxEdges
       iE = geom%edgesOnCell(j,iC)

       psi(:,geom%cellsOnEdge(1,iE)) = psi(:,geom%cellsOnEdge(1,iE)) + MPAS_JEDI_HALF_kr * geom%dvEdge(iE)&
 * sin(MPAS_JEDI_ZERO_kr-geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC) * psi_line_intg_dx(:,iC)
       psi(:,geom%cellsOnEdge(2,iE)) = psi(:,geom%cellsOnEdge(2,iE)) + MPAS_JEDI_HALF_kr * geom%dvEdge(iE)&
 * sin(MPAS_JEDI_ZERO_kr-geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC) * psi_line_intg_dx(:,iC)
       
       psi(:,geom%cellsOnEdge(1,iE)) = psi(:,geom%cellsOnEdge(1,iE)) + MPAS_JEDI_HALF_kr * geom%dvEdge(iE)&
 * cos(geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC) * psi_line_intg_dy(:,iC)
       psi(:,geom%cellsOnEdge(2,iE)) = psi(:,geom%cellsOnEdge(2,iE)) + MPAS_JEDI_HALF_kr * geom%dvEdge(iE)&
 * cos(geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC) * psi_line_intg_dy(:,iC)

       chi(:,geom%cellsOnEdge(1,iE)) = chi(:,geom%cellsOnEdge(1,iE)) + MPAS_JEDI_HALF_kr * geom%dvEdge(iE)&
 * sin(MPAS_JEDI_ZERO_kr-geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC) * chi_line_intg_dx(:,iC)
       chi(:,geom%cellsOnEdge(2,iE)) = chi(:,geom%cellsOnEdge(2,iE)) + MPAS_JEDI_HALF_kr * geom%dvEdge(iE)&
 * sin(MPAS_JEDI_ZERO_kr-geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC) * chi_line_intg_dx(:,iC)

       chi(:,geom%cellsOnEdge(1,iE)) = chi(:,geom%cellsOnEdge(1,iE)) + MPAS_JEDI_HALF_kr * geom%dvEdge(iE)&
 * cos(geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC) * chi_line_intg_dy(:,iC)
       chi(:,geom%cellsOnEdge(2,iE)) = chi(:,geom%cellsOnEdge(2,iE)) + MPAS_JEDI_HALF_kr * geom%dvEdge(iE)&
 * cos(geom%angleEdge(iE)) * geom%edgesOnCell_sign(j,iC) * chi_line_intg_dy(:,iC)
     enddo !- j
   enddo !- iC

end subroutine psichi_to_uv_centerAD

!-------------------------------------------------------------------------------
! input  : psi @ cell center
! output : psi @ vertices
subroutine psichi_to_uv_edge_step1(geom, psi, psi_v)

   implicit none
   type (mpas_geom),                                                  intent(in)    :: geom         !< geometry
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nCells),    intent(in)    :: psi
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nVertices), intent(inout) :: psi_v

   integer :: iC, iV, j

   psi_v=MPAS_JEDI_ZERO_kr

   ! Interpolate psi in cell center to vertice
   do iV = 1, geom%nVerticesSolve ! local
     do j = 1, geom%vertexDegree
       iC = geom%cellsOnVertex(j,iV)
       psi_v(:,iV) = psi_v(:,iV) + geom%kiteAreasOnVertex(j,iV) * psi(:,iC)
     enddo
     psi_v(:,iV) = psi_v(:,iV) / geom%areaTriangle(iV)
   enddo

end subroutine psichi_to_uv_edge_step1

!-------------------------------------------------------------------------------
! input  : psi @ vertices, chi @ cell center
! output : edge_normal_wind @ edges
subroutine psichi_to_uv_edge_step2(geom, psi_v, chi, edge_normal_wind)

   implicit none
   type (mpas_geom),                                                  intent(in)    :: geom         !< geometry
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nVertices), intent(in)    :: psi_v
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nCells),    intent(in)    :: chi
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nEdges),    intent(inout) :: edge_normal_wind

   integer :: iE

   edge_normal_wind=MPAS_JEDI_ZERO_kr

   !get edge_normal_wind
   do iE = 1, geom%nEdgesSolve ! local
     edge_normal_wind(:,iE) = edge_normal_wind(:,iE) - &
                ( chi(:,geom%cellsOnEdge(2,iE)) - chi(:,geom%cellsOnEdge(1,iE)) ) / geom%dcEdge(iE) - &
                ( psi_v(:,geom%verticesOnEdge(2,iE)) - psi_v(:,geom%verticesOnEdge(1,iE)) ) / geom%dvEdge(iE)
   enddo

end subroutine psichi_to_uv_edge_step2

!-------------------------------------------------------------------------------
! input  : edge_normal_wind @ edges
! output : zonal- and meridional- wind @ cell center
subroutine psichi_to_uv_edge_step3(geom, edge_normal_wind, u, v)

   use mpas_vector_reconstruction

   implicit none
   type (mpas_geom),                                               intent(in)    :: geom         !< geometry
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nEdges), intent(in)    :: edge_normal_wind
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nCells), intent(inout) :: u, v

   real (kind=kind_real), dimension(:,:), allocatable :: &
                                        uReconstructX, uReconstructY, uReconstructZ
   type (mpas_pool_type), pointer :: mesh

   allocate(uReconstructX(geom%nVertLevels,geom%nCells))
   allocate(uReconstructY(geom%nVertLevels,geom%nCells))
   allocate(uReconstructZ(geom%nVertLevels,geom%nCells))

   !edge_wind -> u/v @ cell center
   call mpas_pool_get_subpool( geom % domain % blocklist % structs, 'mesh', mesh)
   call mpas_reconstruct(mesh, edge_normal_wind,    &
                         uReconstructX,             &
                         uReconstructY,             &
                         uReconstructZ,             &
                         u,                         &
                         v, .False. ) ! local only, no halo calculation

   deallocate(uReconstructX, uReconstructY, uReconstructZ)

end subroutine psichi_to_uv_edge_step3

!-------------------------------------------------------------------------------
subroutine psichi_to_uv_edge_step1AD(geom, psi, psi_v)

   implicit none
   type (mpas_geom),                                                  intent(in)    :: geom         !< geometry
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nCells),    intent(inout) :: psi
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nVertices), intent(inout) :: psi_v

   integer :: iC, iV, j

   psi=MPAS_JEDI_ZERO_kr

   ! Interpolate psi in cell center to vertice
   do iV = 1, geom%nVerticesSolve ! local
     psi_v(:,iV) = psi_v(:,iV) / geom%areaTriangle(iV)
     do j = 1, geom%vertexDegree
       iC = geom%cellsOnVertex(j,iV)
       psi(:,iC) = psi(:,iC) + geom%kiteAreasOnVertex(j,iV) * psi_v(:,iV)
     enddo
   enddo

end subroutine psichi_to_uv_edge_step1AD

!-------------------------------------------------------------------------------
subroutine psichi_to_uv_edge_step2AD(geom, psi_v, chi, edge_normal_wind)

   implicit none
   type (mpas_geom),                                                  intent(in)    :: geom         !< geometry
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nVertices), intent(inout) :: psi_v
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nCells),    intent(inout) :: chi
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nEdges),    intent(inout) :: edge_normal_wind

   integer :: iE

   chi=MPAS_JEDI_ZERO_kr
   psi_v=MPAS_JEDI_ZERO_kr

   !get edge_normal_wind
   do iE = 1, geom%nEdgesSolve ! local
     chi(:,geom%cellsOnEdge(2,iE))      = chi(:,geom%cellsOnEdge(2,iE))      - edge_normal_wind(:,iE) / geom%dcEdge(iE)
     chi(:,geom%cellsOnEdge(1,iE))      = chi(:,geom%cellsOnEdge(1,iE))      + edge_normal_wind(:,iE) / geom%dcEdge(iE)
     psi_v(:,geom%verticesOnEdge(2,iE)) = psi_v(:,geom%verticesOnEdge(2,iE)) - edge_normal_wind(:,iE) / geom%dvEdge(iE)
     psi_v(:,geom%verticesOnEdge(1,iE)) = psi_v(:,geom%verticesOnEdge(1,iE)) + edge_normal_wind(:,iE) / geom%dvEdge(iE)
   enddo

end subroutine psichi_to_uv_edge_step2AD

!-------------------------------------------------------------------------------
subroutine psichi_to_uv_edge_step3AD(geom, edge_normal_wind, u, v)

   use mpas_vector_reconstruction

   implicit none
   type (mpas_geom),                                               intent(in)    :: geom         !< geometry
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nEdges), intent(inout) :: edge_normal_wind
   real (kind=kind_real), dimension(geom%nVertLevels,geom%nCells), intent(inout) :: u, v

   real (kind=kind_real), dimension(:,:), allocatable :: &
                                        uReconstructX, uReconstructY, uReconstructZ
   type (mpas_pool_type), pointer :: mesh

   edge_normal_wind=MPAS_JEDI_ZERO_kr

   allocate(uReconstructX(geom%nVertLevels,geom%nCells))
   allocate(uReconstructY(geom%nVertLevels,geom%nCells))
   allocate(uReconstructZ(geom%nVertLevels,geom%nCells))

   call mpas_pool_get_subpool( geom % domain % blocklist % structs, 'mesh', mesh)
   call mpas_reconstruct_2dAD(mesh, edge_normal_wind, &
                         uReconstructX,             &
                         uReconstructY,             &
                         uReconstructZ,             &
                         u,                         &
                         v, .False. ) ! local only, no halo calculation

   deallocate(uReconstructX, uReconstructY, uReconstructZ)

end subroutine psichi_to_uv_edge_step3AD

! ------------------------------------------------------------------------------

elemental subroutine pseudorh_to_spechum(pseudorh,spechum,saturation_spechum)

   implicit none
   real (kind=kind_real), intent(in)  :: pseudorh
   real (kind=kind_real), intent(out) :: spechum
   real (kind=kind_real), intent(in)  :: saturation_spechum

   spechum = pseudorh * saturation_spechum / MPAS_JEDI_HUNDRED_kr

end subroutine pseudorh_to_spechum

elemental subroutine pseudorh_to_spechumAD(pseudorh,spechum,saturation_spechum)

   implicit none
   real (kind=kind_real), intent(inout) :: pseudorh
   real (kind=kind_real), intent(inout) :: spechum
   real (kind=kind_real), intent(in)    :: saturation_spechum

   pseudorh = pseudorh + spechum * saturation_spechum / MPAS_JEDI_HUNDRED_kr
   spechum = MPAS_JEDI_ZERO_kr

end subroutine pseudorh_to_spechumAD

elemental subroutine pseudorh_to_spechum_inverse(pseudorh,spechum,saturation_spechum)

   implicit none
   real (kind=kind_real), intent(out) :: pseudorh
   real (kind=kind_real), intent(in)  :: spechum
   real (kind=kind_real), intent(in)  :: saturation_spechum

   pseudorh = spechum / saturation_spechum * MPAS_JEDI_HUNDRED_kr

end subroutine pseudorh_to_spechum_inverse

elemental subroutine pseudorh_to_spechum_inverseAD(pseudorh,spechum,saturation_spechum)

   implicit none
   real (kind=kind_real), intent(inout) :: pseudorh
   real (kind=kind_real), intent(inout) :: spechum
   real (kind=kind_real), intent(in)    :: saturation_spechum

   spechum = spechum + pseudorh / saturation_spechum * MPAS_JEDI_HUNDRED_kr
   pseudorh = MPAS_JEDI_ZERO_kr

end subroutine pseudorh_to_spechum_inverseAD

elemental subroutine da_tp_to_qs( t, p, qs)

   !---------------------------------------------------------------------------
   ! Purpose: Convert T/p to saturation specific humidity.
   !
   !  Method: qs = es_alpha * es / ( p - ( 1 - rd_over_rv ) * es ).
   !          use Rogers & Yau (1989) formula: es = a exp( bTc / (T_c + c) )
   !--------------------------------------------------------------------------

   implicit none

   real (kind=kind_real), intent(in)  :: t, p
   real (kind=kind_real), intent(out) :: qs

   real (kind=kind_real) :: es
   real (kind=kind_real) :: t_c              ! T in degreesC.

   !---------------------------------------------------------------------------
   ! [1.0] initialise:
   !---------------------------------------------------------------------------
   t_c = t - t_kelvin

   !---------------------------------------------------------------------------
   ! [2.0] Calculate saturation vapour pressure:
   !---------------------------------------------------------------------------
   es = es_alpha * exp( es_beta * t_c / ( t_c + es_gamma ) )

   !---------------------------------------------------------------------------
   ! [3.0] Calculate saturation specific humidity:
   !---------------------------------------------------------------------------
   qs = rd_over_rv * es / ( p - rd_over_rv1 * es )

end subroutine da_tp_to_qs

end module mpas_linvarcha_c2a_mod
