! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module mpasjedi_getvalues_mod

use iso_c_binding

! fckit
use fckit_mpi_module,               only: fckit_mpi_sum
use fckit_configuration_module,     only: fckit_configuration

! oops
use datetime_mod,                   only: datetime
use kinds,                          only: kind_real
use unstructured_interpolation_mod

! saber
use interpolatorbump_mod,           only: bump_interpolator

! ufo
use ufo_locations_mod
use ufo_geovals_mod,                only: ufo_geovals
use ufo_vars_mod

!MPAS-Model
use mpas_constants
use mpas_derived_types
use mpas_field_routines
use mpas_kind_types, only: StrKIND
use mpas_pool_routines
use mpas_dmpar, only: mpas_dmpar_exch_halo_field
use mpasjedi_unstructured_interp_mod


!mpas-jedi
use mpas_constants_mod
use mpas_geom_mod
!use mpas_fields_mod
use mpas_fields_mod, only: mpas_fields
use mpas2ufo_vars_mod
use mpas4da_mod

! --------------------------------------------------------------------------------------------------

implicit none
private
public :: mpasjedi_getvalues, mpasjedi_getvalues_base
public :: mpas_getvalues_registry
public :: fill_geovals, getvalues_base_create, getvalues_base_delete

type, abstract :: mpasjedi_getvalues_base
  private
  logical, public :: use_bump_interp
  type(bump_interpolator), public :: bumpinterp
  type(unstrc_interp), public     :: unsinterp
  contains
  procedure :: initialize_uns_interp
  procedure, public :: fill_geovals
  generic, public :: set_trajectory => fill_geovals
  procedure :: integer_interpolation_bump
  procedure :: integer_interpolation_unstructured
end type mpasjedi_getvalues_base

type, extends(mpasjedi_getvalues_base) :: mpasjedi_getvalues
  private
  contains
  procedure, public :: create
  procedure, public :: delete
end type mpasjedi_getvalues

#define LISTED_TYPE mpasjedi_getvalues

!> Linked list interface - defines registry_t type
#include <oops/util/linkedList_i.f>

!> Global registry
type(registry_t) :: mpas_getvalues_registry

! --------------------------------------------------------------------------------------------------

character (len=1024) :: message

contains

! --------------------------------------------------------------------------------------------------
!> Linked list implementation
#include <oops/util/linkedList_c.f>

! ------------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> \brief GetValues base class 'create' logic
!!
!! \details **getvalues_base_create** This subroutine populates the getvalues_base
!! class members. This subroutine is called from the 'create' subroutines of all
!! derived classes. (i.e. getvalues and lineargetvalues)
subroutine getvalues_base_create(self, geom, locs, f_conf)
  implicit none
  class(mpasjedi_getvalues_base), intent(inout) :: self   !< getvalues_base self
  type(mpas_geom),                intent(in)    :: geom   !< geometry (mpas mesh)
  type(ufo_locations),            intent(in)    :: locs   !< ufo geovals (obs) locations
  type(fckit_configuration),      intent(in)    :: f_conf !< configuration

  real(kind=kind_real), allocatable :: lons(:), lats(:)
  integer :: nlocs
  character (len=:), allocatable    :: interp_type

  nlocs = locs%nlocs()
  allocate(lons(nlocs), lats(nlocs))
  call locs%get_lons(lons)
  call locs%get_lats(lats)

  if (f_conf%get("interpolation type", interp_type)) then
    select case (interp_type)
      case ('bump')
        self%use_bump_interp = .True.
      case ('unstructured')
        self%use_bump_interp = .False.
      case default
        write(message,*) '--> getvalues_base_create: interpolation type: ',interp_type,' not implemented'
        call abor1_ftn(message)
    end select
  else
    self%use_bump_interp = .True. ! BUMP is default interpolation
  end if

  if (self%use_bump_interp) then
    call self%bumpinterp%init(geom%f_comm, afunctionspace_in=geom%afunctionspace, lon_out=lons, lat_out=lats, &
      & nl0=geom%nVertLevels)
  else
    call initialize_uns_interp(self, geom, lats, lons)
  endif

  if (allocated(interp_type)) deallocate(interp_type)
  deallocate(lons, lats)

end subroutine getvalues_base_create

! --------------------------------------------------------------------------------------------------

!> \brief GetValues class 'create' logic
!!
!! \details **create** This subroutine contstructs an mpasjedi_getvalues object
!! class instance.
subroutine create(self, geom, locs, f_conf)
  implicit none
  class(mpasjedi_getvalues),      intent(inout) :: self   !< getvalues self
  type(mpas_geom),                intent(in)    :: geom   !< geometry (mpas mesh)
  type(ufo_locations),            intent(in)    :: locs   !< ufo geovals (obs) locations
  type(fckit_configuration),      intent(in)    :: f_conf !< configuration
  call getvalues_base_create(self, geom, locs, f_conf)
end subroutine create

! --------------------------------------------------------------------------------------------------

!> \brief GetValues base class 'delete' logic
!!
!! \details **getvalues_base_delete** This subroutine deletes (frees memory) for
!! the getvalues_base class members. This subroutine is called from the 'delete'
!! subroutines of all derived classes. (i.e. getvalues and lineargetvalues)
subroutine getvalues_base_delete(self)
  class(mpasjedi_getvalues_base), intent(inout) :: self  !< getvalues_base self
  if (self%use_bump_interp) then
    call self%bumpinterp%delete()
  else
    call self%unsinterp%delete()
  endif
end subroutine getvalues_base_delete

! --------------------------------------------------------------------------------------------------

!> \brief GetValues class 'delete' logic
!!
!! \details **delete** This subroutine deletes (frees memory) for an mpasjedi_getvalues object
!! class instance.
subroutine delete(self)
  class(mpasjedi_getvalues), intent(inout) :: self !< getvalues self
  call getvalues_base_delete(self)
end subroutine delete

! --------------------------------------------------------------------------------------------------

!> \brief Interpolates from geovar mpas_fields to populate ufo_geovals
!!
!! \details **fill_geovals** This subroutine populates the variables in a
!! ufo_geovals object by interpolating the state variables in an mpas_fields object.
!! This is the non-linear subroutine used in both GetValues and LinearGetValues classes
subroutine fill_geovals(self, geom, state, t1, t2, locs, gom)
  implicit none
  class(mpasjedi_getvalues_base), intent(inout) :: self    !< getvalues_base self
  type(mpas_geom),                intent(in)    :: geom    !< geometry (mpas mesh)
  type(mpas_fields),              intent(in)    :: state   !< state containing geovars
  type(datetime),                 intent(in)    :: t1      !< time window begin
  type(datetime),                 intent(in)    :: t2      !< time window end
  type(ufo_locations),            intent(in)    :: locs    !< observation locations
  type(ufo_geovals),              intent(inout) :: gom     !< geovals

  logical(c_bool), allocatable :: time_mask(:)
  integer :: jvar, jlev, ilev, jloc, nDims
  integer :: nCells, maxlevels, nlevels, nlocs, nlocsg
  integer, allocatable ::obs_field_int(:,:)
  real(kind=kind_real), allocatable :: mod_field(:,:), obs_field(:,:)

  character(len=MAXVARLEN) :: geovar

  type(mpas_pool_iterator_type) :: poolItr
  real(kind=kind_real), pointer :: ptrr1(:)
  real(kind=kind_real), pointer :: ptrr2(:,:)
  integer, pointer :: ptri1(:)
  integer, pointer :: ptri2(:,:)

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
  allocate(time_mask(nlocs))
  call locs%get_timemask(t1, t2, time_mask)

  ! Interpolate state to obs locations using pre-calculated weights
  ! ----------------------------------------------------------------
  maxlevels = geom%nVertLevelsP1
  allocate(mod_field(nCells,maxlevels))
  allocate(obs_field(nlocs,maxlevels))
  allocate(obs_field_int(nlocs,1))

  call mpas_pool_begin_iteration(state%subFields)
  do while ( mpas_pool_get_next_member(state%subFields, poolItr) )
    if (poolItr % memberType == MPAS_POOL_FIELD) then

      geovar = trim(poolItr % memberName)
      nDims = poolItr % nDims

      jvar = ufo_vars_getindex(gom%variables, geovar)
      if ( jvar < 1 ) cycle

      if (poolItr % dataType == MPAS_POOL_REAL) then
        nlevels = gom%geovals(jvar)%nval

        if (nDims == 1) then
          call state%get(geovar, ptrr1)
          mod_field(:,1) = ptrr1(1:nCells)
        else if (nDims == 2) then
          call state%get(geovar, ptrr2)
          mod_field(:,1:nlevels) = transpose(ptrr2(1:nlevels,1:nCells))
        else
          write(message,*) '--> fill_geovals: nDims == ',nDims,' not handled for reals'
          call abor1_ftn(message)
        endif

        !TODO: make a generic self%apply that encapsulates all this logic and looping:
        ! call self%apply(mod_field(:,1:nlevels), obs_field(:,1:nlevels), nlevels, nCells, nlocs)

        if (self%use_bump_interp) then
          call self%bumpinterp%apply(mod_field(1:nCells,1:nlevels), &
                                     obs_field(:,1:nlevels), &
                                     trans_in=.false.)
        else
          do jlev = 1, nlevels
            call self%unsinterp%apply(mod_field(:,jlev), &
                                      obs_field(:,jlev))
          end do
        endif
        do jloc = 1, nlocs
          if (time_mask(jloc)) then
            do jlev = 1, nlevels
              !BJJ-tmp vertical flip, top-to-bottom for CRTM geoval
              ilev = nlevels - jlev + 1
              gom%geovals(jvar)%vals(ilev,jloc) = obs_field(jloc,jlev)
            end do
          endif
        end do
      else if (poolItr % dataType == MPAS_POOL_INTEGER) then
        if (nDims == 1) then
          call state%get(geovar, ptri1)
          mod_field(:,1) = real(ptri1(1:nCells), kind_real)
        else
          write(message,*) '--> fill_geovals: nDims == ',nDims,' not handled for integers'
          call abor1_ftn(message)
        endif

        jvar = ufo_vars_getindex(gom%variables, poolItr % memberName)
        if (self%use_bump_interp) then
          call self%integer_interpolation_bump(nCells, nlocs, &
            ptri1, obs_field_int, gom, jvar, time_mask)
        else
          call self%integer_interpolation_unstructured(nCells, nlocs, &
            ptri1, mod_field(:,1), gom, jvar, time_mask)
        endif
      end if

    endif
  end do !- end of pool iteration
  deallocate(mod_field)
  deallocate(obs_field)
  deallocate(obs_field_int)
  deallocate(time_mask)

end subroutine fill_geovals

! --------------------------------------------------------------------------------------------------

!> \brief Initializes an unstructured interpolation type
!!
!! \details **initialize_uns_interp** This subroutine calls unsinterp%create,
!! which calculates the barycentric weights used to interpolate data between the
!! mpas mesh locations (grid) and the observation locations.
subroutine initialize_uns_interp(self, grid, lats_obs, lons_obs)

  implicit none
  class(mpasjedi_getvalues_base), intent(inout) :: self        !< self
  type(mpas_geom),          intent(in)          :: grid        !< mpas mesh data
  real(kind=kind_real), allocatable, intent(in) :: lats_obs(:) !< latitudes of obs
  real(kind=kind_real), allocatable, intent(in) :: lons_obs(:) !< longitudes of obs

  integer :: nn, ngrid_in
  character(len=8) :: wtype = 'barycent'
  real(kind=kind_real), allocatable :: lats_in(:), lons_in(:)

  ! Get the Solution dimensions
  ! ---------------------------
  ngrid_in = grid%nCellsSolve

  !Calculate interpolation weight
  !------------------------------------------
  allocate( lats_in(ngrid_in) )
  allocate( lons_in(ngrid_in) )
  lats_in(:) = grid%latCell( 1:ngrid_in ) * MPAS_JEDI_RAD2DEG_kr !- to Degrees
  lons_in(:) = grid%lonCell( 1:ngrid_in ) * MPAS_JEDI_RAD2DEG_kr !- to Degrees

  ! Initialize unsinterp
  ! ---------------
  nn = 3 ! number of nearest neigbors
  call self%unsinterp%create(grid%f_comm, nn, wtype, &
                            ngrid_in, lats_in, lons_in, &
                            size(lats_obs), lats_obs, lons_obs)

  ! Release memory
  ! --------------
  deallocate(lats_in)
  deallocate(lons_in)

end subroutine initialize_uns_interp

! ------------------------------------------------------------------------------

!> \brief Performs interpolation of integer fields using BUMP
!!
!! \details **integer_interpolation_bump** This subroutine performs the interpolation
!! of integer-valued fields (i.e. types) using BUMP
subroutine integer_interpolation_bump(self, ngrid, nlocs, &
                                      data_in, obs_field_int, gom, jvar, time_mask)
  implicit none

  class(mpasjedi_getvalues_base), intent(inout) :: self     !< self
  integer, intent(in)                 :: ngrid              !< number of cells in model mesh
  integer, intent(in)                 :: nlocs              !< number of locations for obs
  integer, dimension(:), intent(in)   :: data_in            !< data to interpolate
  integer, allocatable, intent(inout) :: obs_field_int(:,:) !< output array of interpolated data
  type(ufo_geovals), intent(inout)    :: gom                !< output geoVaLs
  integer, intent(in)                 :: jvar               !< gom % geovals index
  logical(c_bool), intent(in)         :: time_mask(nlocs)   !< mask for time window

  integer :: jloc

  ! This code assumes time_mask already allocated to size nlocs and poplulated earlier.

  ! by default, bumpinterp%apply uses nearest neighbor interpolation for integer types
  call self%bumpinterp%apply(data_in(1:ngrid), obs_field_int)
  do jloc = 1, nlocs
    if (time_mask(jloc)) then
! JJG: apply was here, but moved outside loop as it should only get called once
!      call self%bumpinterp%apply(data_in(1:ngrid), obs_field_int)
      gom%geovals(jvar)%vals(1,jloc) = real(obs_field_int(jloc,1), kind_real)
    endif
  enddo

end subroutine integer_interpolation_bump

! ------------------------------------------------------------------------------

!> \brief Performs interpolation of integer fields using unstructured interpolation
!!
!! \details **integer_interpolation_unstructured** This subroutine performs the interpolation
!! of integer-valued fields (i.e. types) using unstructured interpolation
subroutine integer_interpolation_unstructured(self, ngrid, nlocs, &
                                           data_in, work_field, gom, jvar, time_mask)
  implicit none

  class(mpasjedi_getvalues_base), intent(inout) :: self    !< self
  integer, intent(in)                 :: ngrid             !< number of cells in model mesh
  integer, intent(in)                 :: nlocs             !< number of locations for obs
  integer, dimension(:), intent(in)   :: data_in           !< data to interpolate
  real(kind=kind_real), intent(inout) :: work_field(ngrid) !< (ngrid,1) array
  type(ufo_geovals), intent(inout)    :: gom               !< output geoVaLs
  integer, intent(in)                 :: jvar              !< gom % geovals index
  logical(c_bool),  intent(in)        :: time_mask(nlocs)  !< mask for time window

  integer :: jloc
  real(kind=kind_real), dimension(nlocs) :: interpolated_data

  ! This code assumes time_mask already allocated to size nlocs and poplulated earlier.
  ! Also, that work_field already allocated to size (ngrid, 1)

  work_field = real(data_in(1:ngrid), kind_real)
  call unsinterp_integer_apply(self%unsinterp, work_field, interpolated_data)
  do jloc = 1, nlocs
    if (time_mask(jloc)) then
      gom%geovals(jvar)%vals(1,jloc) = interpolated_data(jloc)
    endif
  enddo

end subroutine integer_interpolation_unstructured

end module mpasjedi_getvalues_mod
