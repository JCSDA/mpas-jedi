! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module mpasjedi_lineargetvalues_mod

! oops
use datetime_mod,                   only: datetime
use type_bump,                      only: bump_type
use kinds,                          only: kind_real

! ufo
use ufo_locs_mod,                   only: ufo_locs, ufo_locs_time_mask
use ufo_geovals_mod,                only: ufo_geovals
use ufo_vars_mod

!MPAS-Model
use mpas_constants
use mpas_derived_types
use mpas_field_routines
use mpas_kind_types, only: StrKIND
use mpas_pool_routines

!mpas-jedi
use mpas_constants_mod
use mpas_geom_mod
use mpas_field_utils_mod
use mpas2ufo_vars_mod
use mpas4da_mod
use mpasjedi_getvalues_mod

! --------------------------------------------------------------------------------------------------

implicit none
private
public :: mpasjedi_lineargetvalues
public :: mpas_lineargetvalues_registry

type, extends(mpasjedi_getvalues_base) :: mpasjedi_lineargetvalues
   type (mpas_pool_type), pointer :: pool_traj
   contains
      procedure, public :: create
      procedure, public :: delete
      procedure, public :: set_trajectory
      procedure, public :: fill_geovals_tl
      procedure, public :: fill_geovals_ad
end type mpasjedi_lineargetvalues

#define LISTED_TYPE mpasjedi_lineargetvalues

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_lineargetvalues_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

! TODO: Slightly more appropriate to move create & delete to base class for
! bump init/dealloc, and just extend them here to also handle the trajectory pool.
subroutine create(self, geom, locs)
   implicit none
   class(mpasjedi_lineargetvalues), intent(inout) :: self
   type(mpas_geom),                 intent(in)    :: geom
   type(ufo_locs),                  intent(in)    :: locs

   call initialize_bump(self, geom, locs)
   call mpas_pool_create_pool(self%pool_traj)

 end subroutine create


! ------------------------------------------------------------------------------

 subroutine delete(self)

   class(mpasjedi_lineargetvalues), intent(inout) :: self

   call mpas_pool_destroy_pool(self%pool_traj)
   call self%bump%dealloc()

 end subroutine delete

! ------------------------------------------------------------------------------


subroutine set_trajectory (self, geom, fields, t1, t2, locs, gom)

   class(mpasjedi_lineargetvalues), intent(inout) :: self
   type(mpas_geom),                intent(in)    :: geom
   type(mpas_field),               intent(in)    :: fields
   type(datetime),                 intent(in)    :: t1
   type(datetime),                 intent(in)    :: t2
   type(ufo_locs),                 intent(in)    :: locs
   type(ufo_geovals),              intent(inout) :: gom

   type (mpas_pool_type), pointer :: pool_tmp  !< temporary pool for setting trajectory
   type (field2DReal), pointer :: field2d_src => null() !< for setting trajectory
   type (field1DReal), pointer :: field1d_src => null() !< for setting trajectory

   ! Re-set trajectory pool (maybe there's a better way?)
   ! ---------------------------------------
   call mpas_pool_destroy_pool(self%pool_traj)
   call mpas_pool_create_pool(self%pool_traj)

   ! Initialize the interpolation trajectory
   ! ---------------------------------------
   call mpas_pool_create_pool( pool_tmp )

   call mpas_pool_get_field(fields % subFields, 'temperature', field2d_src)
   call mpas_pool_add_field(pool_tmp, 'temperature', field2d_src)
   call mpas_pool_get_field(fields % subFields, 'spechum', field2d_src)
   call mpas_pool_add_field(pool_tmp, 'spechum', field2d_src)
   call mpas_pool_get_field(fields % subFields, 'pressure', field2d_src)
   call mpas_pool_add_field(pool_tmp, 'pressure', field2d_src)
   call mpas_pool_get_field(fields % subFields, 'surface_pressure', field1d_src)
   call mpas_pool_add_field(pool_tmp, 'surface_pressure', field1d_src)

   call mpas_pool_clone_pool(pool_tmp, self%pool_traj)

   call mpas_pool_empty_pool(pool_tmp)
   call mpas_pool_destroy_pool(pool_tmp)

   ! NOTE: this call may get removed soon. See OOPS PR #577.
   call fill_geovals(self, geom, fields, t1, t2, locs, gom)

end subroutine set_trajectory

subroutine fill_geovals_tl(self, geom, fields, t1, t2, locs, gom)

    class(mpasjedi_lineargetvalues), intent(inout) :: self
    type(mpas_geom),                intent(in)    :: geom
    type(mpas_field),               intent(inout) :: fields
    type(datetime),                 intent(in)    :: t1
    type(datetime),                 intent(in)    :: t2
    type(ufo_locs),                 intent(in)    :: locs
    type(ufo_geovals),              intent(inout) :: gom

    character(len=*), parameter :: myname = 'fill_geovals_tl'

    logical, allocatable :: time_mask(:)
    integer :: jvar, jlev, ilev, jloc, ngrid, maxlevels, nlevels, nlocs, ivar
    real(kind=kind_real), allocatable :: mod_field(:,:)
    real(kind=kind_real), allocatable :: obs_field(:,:)

    type (mpas_pool_type), pointer :: pool_ufo
    type (mpas_pool_iterator_type) :: poolItr
    real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
    real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
    real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
    real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
    character(len=MAXVARLEN), allocatable :: ufo_vars(:)

   ! Get grid dimensions and checks
   ! ------------------------------
    ngrid = geom%nCellsSolve !or traj%ngrid
    nlocs = locs%nlocs           !or traj%nobs
 !   write(*,*)'fill_geovals_tl: ngrid, nlocs = : ',ngrid, nlocs
    !call interp_checks("tl", inc, locs, vars, gom)

    !Make sure the return values are allocated and set
    !-------------------------------------------------
    do jvar=1,gom%nvar
       if( .not. allocated(gom%geovals(jvar)%vals) )then
          ! air_pressure is required for GNSSRO, but mpas-jedi does not have corresponding TL/AD
          ! code in mpas2ufo_vars_mod.  Temporarily perform allocation here for all "vars".
          ! TODO: move this allocation for 2D variables into loop over mpas_pool members after
          !       either
          !   (1) adding TL/AD for air_pressure when added as analysis variable and move
          !       this allocation for 2D variables into loop over mpas_pool members
          !    OR
          !   (2) implementing increment variable change for mpas-jedi, including
          !       dependencies on non-analyzed variables
          !if ( geom%nVertLevels > 1) then
          gom%geovals(jvar)%nval = geom%nVertLevels
          !else
          !   gom%geovals(jvar)%nval = 1
          !end if
          allocate( gom%geovals(jvar)%vals(gom%geovals(jvar)%nval,gom%geovals(jvar)%nlocs) )
          gom%geovals(jvar)%vals = MPAS_JEDI_ZERO_kr
 !         write(*,*) ' gom%geovals(n)%vals allocated'
          gom%linit = .true.
       endif
    enddo

  ! Get mask for locations in this time window
  ! ------------------------------------------
    call ufo_locs_time_mask(locs, t1, t2, time_mask)

    !Interpolate fields to obs locations using pre-calculated weights
    !----------------------------------------------------------------
    allocate(ufo_vars(gom%nvar))
    do ivar = 1, gom%nvar
       ufo_vars(ivar) = trim(gom%variables(ivar))
    end do
 !   write(0,*)'fill_geovals_tl: gom%nvar   : ',gom%nvar
 !   write(0,*)'fill_geovals_tl: vars%varlist : ',ufo_vars
 !   write(0,*)'fill_geovals_tl: nlocs        : ',nlocs

    !------- need some table matching UFO_Vars & related MPAS_Vars
    !------- for example, Tv @ UFO may require Theta, Pressure, Qv.
    !-------                               or  Theta_m, exner_base, Pressure_base, Scalar(&index_qv)
    call convert_mpas_field2ufoTL(geom, self%pool_traj, fields % subFields, pool_ufo, ufo_vars, gom%nvar, ngrid) !--pool_ufo is new pool with ufo_vars

    maxlevels = geom%nVertLevelsP1
    allocate(mod_field(ngrid,maxlevels))
    allocate(obs_field(nlocs,maxlevels))
    mod_field = MPAS_JEDI_ZERO_kr
    obs_field = MPAS_JEDI_ZERO_kr

    call mpas_pool_begin_iteration(pool_ufo)
    do while ( mpas_pool_get_next_member(pool_ufo, poolItr) )
       if (poolItr % memberType == MPAS_POOL_FIELD) then
          ivar = ufo_vars_getindex( ufo_vars,trim(poolItr % memberName))
          if ( ivar == -1 ) cycle

 !         write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)

          nlevels = gom%geovals(ivar)%nval
          if (poolItr % nDims == 1) then
             ! Temporarily treated for 2D surface_pressure: we use ps as a control variable now, hofx is identity. In the future, if we use 3D P as the
             ! control variable, this part will be modified, together with 'air_pressure' in convert_mpas_field2ufoTL.
            if ( .not.allocated(gom%geovals(ivar)%vals) .OR. &
               size(gom%geovals(ivar)%vals,1) /= 1 .OR. &
               gom%geovals(ivar)%nval /= 1 ) then
                  deallocate(gom%geovals(ivar)%vals)
                  gom%geovals(ivar)%nval = 1
                  allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
                  gom%geovals(ivar)%vals = MPAS_JEDI_ZERO_kr
                  gom%linit = .true.
             endif
             nlevels = gom%geovals(ivar)%nval
             if (poolItr % dataType == MPAS_POOL_REAL) then
                call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r1d_ptr_a)
                mod_field(:,1) = r1d_ptr_a(1:ngrid)
             endif
          else if (poolItr % nDims == 2) then
 
             if (poolItr % dataType == MPAS_POOL_REAL) then
                call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r2d_ptr_a)
 !              write(*,*) "fill_geovals_tl: var=",trim(poolItr % memberName)
 !              write(*,*) "ufo_vars_getindex: ivar=",ivar

                do jlev = 1, nlevels
                   mod_field(:,jlev) = r2d_ptr_a(jlev,1:ngrid)
                end do
             end if
 !         else if (poolItr % nDims == 3) then

          end if

          !TODO - JJG: Reduce wall-time of getvalues/_tl/_ad
          ! + apply_obsop takes ~99.9% of wall-time of getvalues_tl on cheyenne and
          !   scales with node count. Seems to have MPI-related issue.
          !
          self%bump%geom%nl0 = nlevels
          call self%bump%apply_obsop(mod_field(:,1:nlevels),obs_field(:,1:nlevels))
 
          do jlev = 1, nlevels
             ilev = nlevels - jlev + 1
             do jloc = 1, nlocs
                !BJJ-tmp vertical flip, top-to-bottom for CRTM geoval
                ! only selected obs (using locs%indx()) are filling "geovals"
                if (time_mask(jloc)) gom%geovals(ivar)%vals(ilev,locs%indx(jloc)) = obs_field(jloc,jlev)
             end do
          end do

       end if
    end do
    deallocate(mod_field)
    deallocate(obs_field)

    self%bump%geom%nl0 = 1
 !   allocate(mod_field(ngrid,1))
 !   allocate(obs_field(nlocs,1))

 !    ! Special cases go here

 !   deallocate(mod_field)
 !   deallocate(obs_field)

    call mpas_pool_destroy_pool(pool_ufo)

    if (allocated(ufo_vars)) deallocate(ufo_vars)

 !   write(*,*) '---- Leaving fill_geovals_tl ---'

end subroutine fill_geovals_tl

subroutine fill_geovals_ad(self, geom, fields, t1, t2, locs, gom)

    class(mpasjedi_lineargetvalues), intent(inout) :: self
    type(mpas_geom),                intent(in)     :: geom
    type(mpas_field),               intent(inout)  :: fields
    type(datetime),                 intent(in)     :: t1
    type(datetime),                 intent(in)     :: t2
    type(ufo_locs),                 intent(in)     :: locs
    type(ufo_geovals),              intent(inout)  :: gom

    character(len=*), parameter :: myname = 'fill_geovals_ad'

    logical, allocatable :: time_mask(:)
    integer :: jvar, jlev, ilev, jloc, ngrid, maxlevels, nlevels, nlocs, ivar
    real(kind=kind_real), allocatable :: mod_field(:,:)
    real(kind=kind_real), allocatable :: obs_field(:,:)

    type (mpas_pool_type), pointer :: pool_ufo
    type (mpas_pool_iterator_type) :: poolItr
    real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
    real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
    real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
    real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
    character(len=MAXVARLEN), allocatable :: ufo_vars(:)

   ! If no observations can early exit
   ! ---------------------------------
!    if (self%traj%noobs)  return

  ! Get mask for locations in this time window
  ! ------------------------------------------
    call ufo_locs_time_mask(locs, t1, t2, time_mask)

   ! Get grid dimensions and checks
   ! ------------------------------
    ngrid = geom%nCellsSolve !or traj%ngrid
    nlocs = locs%nlocs       !or traj%nlocs

   !Interpolate fields to obs locations using pre-calculated weights
   !----------------------------------------------------------------
    allocate(ufo_vars(gom%nvar))
    do ivar = 1, gom%nvar
       ufo_vars(ivar) = trim(gom%variables(ivar))
    end do
 !   write(0,*)'getvalues_ad: nlocs        : ',nlocs

    !NOTE: This TL routine is called JUST to create "pool_ufo". Their values from TL routine doesn't matter.
    !    : Actually their values are initialized as "zero" in following "do while" loop.
    call convert_mpas_field2ufoTL(geom, self%pool_traj, fields % subFields, pool_ufo, ufo_vars, gom%nvar, ngrid) !--pool_ufo is new pool with ufo_vars

    maxlevels = geom%nVertLevelsP1
    allocate(mod_field(ngrid,maxlevels))
    allocate(obs_field(nlocs,maxlevels))
    mod_field = MPAS_JEDI_ZERO_kr
    obs_field = MPAS_JEDI_ZERO_kr

    call mpas_pool_begin_iteration(pool_ufo)
    do while ( mpas_pool_get_next_member(pool_ufo, poolItr) )
       if (poolItr % memberType == MPAS_POOL_FIELD) then
          ivar = ufo_vars_getindex( ufo_vars,trim(poolItr % memberName))
          if ( ivar == -1 ) cycle

 !         write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)
          nlevels = gom%geovals(ivar)%nval
          do jlev = 1, nlevels
             !ORG- obs_field(:,jlev) = gom%geovals(ivar)%vals(jlev,:)
             !BJJ-tmp vertical flip, top-to-bottom for CRTM geoval
             ! only selected obs (using locs%indx()) are filling "geovals"
             ilev = nlevels - jlev + 1
             do jloc = 1, nlocs
                if (time_mask(jloc)) then
                  obs_field(jloc,jlev) = gom%geovals(ivar)%vals(ilev, locs%indx(jloc))
                  gom%geovals(ivar)%vals(ilev, locs%indx(jloc)) = MPAS_JEDI_ZERO_kr
                end if
             end do
          end do
 
          !TODO - JJG: Reduce wall-time of getvalues/_tl/_ad
          ! + apply_obsop takes ~99.9% of wall-time of getvalues_ad on cheyenne and
          !   scales with node count. Seems to have MPI-related issue.
          !
          self%bump%geom%nl0 = nlevels
          call self%bump%apply_obsop_ad(obs_field(:,1:nlevels),mod_field(:,1:nlevels))
 
          if (poolItr % nDims == 1) then
             call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r1d_ptr_a)
             r1d_ptr_a=MPAS_JEDI_ZERO_kr
             r1d_ptr_a(1:ngrid) = r1d_ptr_a(1:ngrid) + mod_field(:,1)
          else if (poolItr % nDims == 2) then
 
             if (poolItr % dataType == MPAS_POOL_REAL) then
                call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r2d_ptr_a)
                r2d_ptr_a=MPAS_JEDI_ZERO_kr
 !               write(*,*) "Interp. var=",trim(poolItr % memberName)
 !               write(*,*) "ufo_vars_getindex, ivar=",ivar
                do jlev = 1, nlevels
                   r2d_ptr_a(jlev,1:ngrid) = r2d_ptr_a(jlev,1:ngrid) + mod_field(:,jlev)
                end do
             end if
 !         else if (poolItr % nDims == 3) then
          end if
       end if
    end do
    deallocate(mod_field)
    deallocate(obs_field)

    call convert_mpas_field2ufoAD(geom, self%pool_traj, fields % subFields, pool_ufo, ufo_vars, gom%nvar, ngrid) !--pool_ufo is new pool with ufo_vars

    self%bump%geom%nl0 = 1
 !   allocate(mod_field(ngrid,1))
 !   allocate(obs_field(nlocs,1))

 !    ! Special cases go here

 !   deallocate(mod_field)
 !   deallocate(obs_field)

    call mpas_pool_destroy_pool(pool_ufo)

    if (allocated(ufo_vars)) deallocate(ufo_vars)

 !   write(*,*) '---- Leaving fill_geovals_ad ---'

end subroutine fill_geovals_ad

end module mpasjedi_lineargetvalues_mod
