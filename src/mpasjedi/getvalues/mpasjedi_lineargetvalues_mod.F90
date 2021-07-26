! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module mpasjedi_lineargetvalues_mod

use iso_c_binding

! fckit
use fckit_mpi_module, only: fckit_mpi_sum
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log

! oops
use datetime_mod, only: datetime, datetime_to_string
use kinds, only: kind_real

! saber
use interpolatorbump_mod, only: bump_interpolator

! ufo
use ufo_locations_mod
use ufo_geovals_mod, only: ufo_geovals
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
use mpas_fields_mod
use mpas2ufo_vars_mod
use mpas4da_mod
use mpasjedi_getvalues_mod

! --------------------------------------------------------------------------------------------------

implicit none
private
public :: mpasjedi_lineargetvalues
public :: mpas_lineargetvalues_registry

type, extends(mpasjedi_getvalues_base) :: mpasjedi_lineargetvalues
  private
  real(kind=kind_real), allocatable :: obs_field(:,:), mod_field(:,:)
  logical(c_bool), allocatable :: time_mask(:)
  contains
  procedure, public :: create
  procedure, public :: delete
  procedure, public :: fill_geovals_tl
  procedure, public :: fill_geovals_ad
end type mpasjedi_lineargetvalues

integer, parameter    :: max_string=8000
character(max_string) :: message

#define LISTED_TYPE mpasjedi_lineargetvalues

!> Linked list interface - defines registry_t type
#include <oops/util/linkedList_i.f>

!> Global registry
type(registry_t) :: mpas_lineargetvalues_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------
!> Linked list implementation
#include <oops/util/linkedList_c.f>

! --------------------------------------------------------------------------------------------------

!> \brief LinearGetValues class 'create' logic
!!
!! \details **create** This subroutine constructs an mpasjedi_lineargetvalues object
!! class instance.
subroutine create(self, geom, locs, f_conf)
  implicit none
  class(mpasjedi_lineargetvalues), intent(inout) :: self   !< lineargetvalues self
  type(mpas_geom),                 intent(in)    :: geom   !< geometry (mpas mesh)
  type(ufo_locations),             intent(in)    :: locs   !< ufo geovals (obs) locations
  type(fckit_configuration),       intent(in)    :: f_conf !< configuration

  integer :: nlocs, maxlevels

  call getvalues_base_create(self, geom, locs, f_conf)

  ! Get grid dimensions
  ! -------------------
  nlocs = locs % nlocs() ! # of location for entire window
  maxlevels = geom%nVertLevelsP1

  if (allocated(self%time_mask) .or. allocated(self%obs_field) .or. allocated(self%mod_field)) then
    ! If we're in here this subroutine has not been called as it is intended and someone should know.
    call abor1_ftn('--> mpasjedi_lineargetvalues::create subroutine called when member arrays already allocated.')
  end if

  ! These working arrays, used in fill_geovals_tl and fill_geovals_ad, are not repeatedly allocated
  ! and deallocated during the inner minimization loop by instead allocating them here.
  ! They are deallocated in the delete subroutine.
  allocate(self%time_mask(nlocs))
  allocate(self%obs_field(nlocs,maxlevels))
  if (.not. self%use_bump_interp) then
    allocate(self%mod_field(geom%nCellsSolve, maxlevels))
  end if

end subroutine create


! --------------------------------------------------------------------------------------------------

!> \brief LinearGetValues class 'delete' logic
!!
!! \details **delete** This subroutine deletes (frees memory) for an mpas_jedilineargetvalues object
!! class instance.
subroutine delete(self)
  implicit none
  class(mpasjedi_lineargetvalues), intent(inout) :: self !< lineargetvalues self

  if (allocated(self%time_mask)) deallocate(self%time_mask)
  if (allocated(self%obs_field)) deallocate(self%obs_field)
  if (allocated(self%mod_field)) deallocate(self%mod_field)
  call getvalues_base_delete(self)

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine fill_geovals_tl(self, geom, inc, t1, t2, locs, gom)
  implicit none
  class(mpasjedi_lineargetvalues), intent(inout) :: self !< lineargetvalues self
  type(mpas_geom),                 intent(in)    :: geom !< geometry (mpas mesh)
  class(mpas_fields),              intent(in)    :: inc  !< increment containing geovars
  type(datetime),                  intent(in)    :: t1   !< time window begin
  type(datetime),                  intent(in)    :: t2   !< time window end
  type(ufo_locations),             intent(in)    :: locs !< observation locations
  type(ufo_geovals),               intent(inout) :: gom  !< geovals

  character(len=*), parameter :: myname = 'fill_geovals_tl'

  integer :: jvar, jlev, ilev, jloc, nDims
  integer :: nCells, nlevels, nlocs, nlocsg

  type (mpas_pool_data_type), pointer :: gdata
  type (mpas_pool_iterator_type) :: poolItr
  character (len=MAXVARLEN) :: geovar

  logical :: allocateGeo

  ! Get grid dimensions and checks
  ! ------------------------------
  nCells = geom % nCellsSolve
  nlocs = locs % nlocs() ! # of location for entire window

  ! If no observations can early exit
  ! ---------------------------------
  call geom%f_comm%allreduce(nlocs,nlocsg,fckit_mpi_sum())
  if (nlocsg == 0) then
    return
  endif

  ! Get mask for locations in this time window
  ! ------------------------------------------
  call locs%get_timemask(t1, t2, self%time_mask)

  ! TL of interpolate fields to obs locations using pre-calculated weights
  ! ----------------------------------------------------------------------
  call mpas_pool_begin_iteration(inc%subFields)
  do while ( mpas_pool_get_next_member(inc%subFields, poolItr) )
    if (poolItr % memberType == MPAS_POOL_FIELD) then

      geovar = trim(poolItr % memberName)
      nDims = poolItr % nDims
      call inc%get(geovar, gdata)

      jvar = ufo_vars_getindex(gom%variables, geovar)
      if ( jvar < 1 ) cycle

      self%obs_field = MPAS_JEDI_ZERO_kr

      nlevels = gom%geovals(jvar)%nval

      if (poolItr % dataType == MPAS_POOL_REAL) then
        if (nDims == 1) then
          if (self%use_bump_interp) then
            call self%bumpinterp%apply(gdata%r1%array(1:nCells), &
                                       self%obs_field(:,1))
          else
            call self%unsinterp%apply(gdata%r1%array(1:nCells), &
                                      self%obs_field(:,1))
          endif
        else if (nDims == 2) then
          if (self%use_bump_interp) then
            call self%bumpinterp%apply(gdata%r2%array(1:nlevels,1:nCells), &
                                       self%obs_field(:,1:nlevels), &
                                       trans_in=.true.)
          else
            ! Transpose to get the slices we pass to 'apply' contiguous in memory
            self%mod_field(:,1:nlevels) = transpose(gdata%r2%array(1:nlevels,1:nCells))
            do jlev = 1, nlevels
              call self%unsinterp%apply(self%mod_field(:,jlev), &
                                        self%obs_field(:,jlev))
            enddo
          endif
        end if
      else
        call abor1_ftn('--> fill_geovals_tl: not a real field')
      end if

      do jlev = 1, nlevels
        ilev = nlevels - jlev + 1
        do jloc = 1, nlocs
          !BJJ-tmp vertical flip, top-to-bottom for CRTM geoval
          if (self%time_mask(jloc)) gom%geovals(jvar)%vals(ilev,jloc) = self%obs_field(jloc,jlev)
        end do
      end do
    endif
  end do

end subroutine fill_geovals_tl

subroutine fill_geovals_ad(self, geom, inc, t1, t2, locs, gom)
  implicit none
  class(mpasjedi_lineargetvalues), intent(inout) :: self !< lineargetvalues self
  type(mpas_geom),                 intent(in)    :: geom !< geometry (mpas mesh)
  class(mpas_fields),              intent(inout) :: inc  !< increment containing geovars
  type(datetime),                  intent(in)    :: t1   !< time window begin
  type(datetime),                  intent(in)    :: t2   !< time window end
  type(ufo_locations),             intent(in)    :: locs !< observation locations
  type(ufo_geovals),               intent(in)    :: gom  !< geovals

  character(len=*), parameter :: myname = 'fill_geovals_ad'

  integer :: jvar, jlev, ilev, jloc, nDims
  integer :: nCells, nlevels, nlocs, nlocsg

  type (mpas_pool_data_type), pointer :: gdata
  type (mpas_pool_iterator_type) :: poolItr
  character (len=MAXVARLEN) :: geovar

  ! Get grid dimensions and checks
  ! ------------------------------
  nCells = geom % nCellsSolve
  nlocs = locs % nlocs() ! # of location for entire window
  write(message,*) 'fill_geovals_ad: nlocs        : ',nlocs
  call fckit_log%debug(message)

  ! Get mask for locations in this time window
  ! ------------------------------------------
  call locs%get_timemask(t1, t2, self%time_mask)

  ! If no observations can early exit
  ! ---------------------------------
  call geom%f_comm%allreduce(nlocs,nlocsg,fckit_mpi_sum())
  if (nlocsg == 0) then
    return
  endif

  ! zero out adjoint geovar fields
  call inc%zeros()

  ! Adjoint of interpolate fields to obs locations using pre-calculated weights
  ! ---------------------------------------------------------------------------
  call mpas_pool_begin_iteration(inc%subFields)
  do while ( mpas_pool_get_next_member(inc%subFields, poolItr) )
    if (poolItr % memberType == MPAS_POOL_FIELD) then

      geovar = trim(poolItr % memberName)
      nDims = poolItr % nDims
      call inc%get(geovar, gdata)

      jvar = ufo_vars_getindex(gom%variables, geovar)
      if ( jvar < 1 ) cycle

      self%obs_field = MPAS_JEDI_ZERO_kr

      write(message,*) 'fill_geovals_ad: nDims, geovar =', nDims , geovar
      call fckit_log%debug(message)
      nlevels = gom%geovals(jvar)%nval
      do jlev = 1, nlevels
        !ORG- obs_field(:,jlev) = gom%geovals(jvar)%vals(jlev,:)
        !BJJ-tmp vertical flip, top-to-bottom for CRTM geoval
        ilev = nlevels - jlev + 1
        do jloc = 1, nlocs
          if (self%time_mask(jloc)) then
            self%obs_field(jloc,jlev) = gom%geovals(jvar)%vals(ilev, jloc)
          endif
        end do
      end do
      if (poolItr % dataType == MPAS_POOL_REAL) then
        if (self%use_bump_interp) then
          if (nDims == 1) then
            call self%bumpinterp%apply_ad(self%obs_field(:,1), &
                                          gdata%r1%array(1:nCells))
          else if (nDims == 2) then
            call self%bumpinterp%apply_ad(self%obs_field(:,1:nlevels), &
                                          gdata%r2%array(1:nlevels,1:nCells), &
                                          trans_in=.true.)
          end if
        else ! Unstructured interpolation
          do jlev = 1, nlevels
            call self%unsinterp%apply_ad(self%mod_field(:,jlev), &
                                         self%obs_field(:,jlev))
          end do
          if (nDims == 1) then
            gdata%r1%array(1:nCells) = self%mod_field(:,1)
          else if (nDims == 2) then
            gdata%r2%array(1:nlevels,1:nCells) = transpose(self%mod_field(:,1:nlevels))
          end if
        endif ! self%use_bump_interp
      else
        call abor1_ftn('--> fill_geovals_ad: not a real field')
      end if

    endif ! poolItr % memberType == MPAS_POOL_FIELD
  end do

end subroutine fill_geovals_ad

end module mpasjedi_lineargetvalues_mod
