! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module mpasjedi_lineargetvalues_mod

! fckit
use fckit_mpi_module,               only: fckit_mpi_sum

! oops
use datetime_mod,                   only: datetime, datetime_to_string
use kinds,                          only: kind_real

! saber
use interpolatorbump_mod,         only: bump_interpolator

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
  type (mpas_pool_type), pointer :: trajectories
  contains
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: set_trajectory
    procedure, public :: fill_geovals_tl
    procedure, public :: fill_geovals_ad
end type mpasjedi_lineargetvalues

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

! TODO: Slightly more appropriate to move create & delete to base class for
! bump init/dealloc, and just extend them here to also handle the trajectory pool.
subroutine create(self, geom, locs)
  implicit none
  class(mpasjedi_lineargetvalues), intent(inout) :: self
  type(mpas_geom),                 intent(in)    :: geom
  type(ufo_locs),                  intent(in)    :: locs

  call self%bumpinterp%init(geom%f_comm, afunctionspace_in=geom%afunctionspace, lon_out=locs%lon, lat_out=locs%lat, &
     & nl=geom%nVertLevels)
  call mpas_pool_create_pool(self % trajectories)

end subroutine create


! ------------------------------------------------------------------------------

subroutine delete(self)
  implicit none
  class(mpasjedi_lineargetvalues), intent(inout) :: self

  call mpas_pool_destroy_pool(self % trajectories)

  call self%bumpinterp%delete()

end subroutine delete

! ------------------------------------------------------------------------------


subroutine set_trajectory (self, geom, fields, t1, t2, locs, gom)
  implicit none
  class(mpasjedi_lineargetvalues), intent(inout) :: self
  type(mpas_geom),                intent(in)    :: geom
  type(mpas_field),               intent(in)    :: fields
  type(datetime),                 intent(in)    :: t1
  type(datetime),                 intent(in)    :: t2
  type(ufo_locs),                 intent(in)    :: locs
  type(ufo_geovals),              intent(inout) :: gom

  type (mpas_pool_type), pointer :: windowtraj  !< trajectory pool for t1 to t2
  character(len=41) :: windowkey

  type (field2DReal), pointer :: field2d_src, field2d_dst
  type (field1DReal), pointer :: field1d_src, field1d_dst


  ! Initialize this window interpolation trajectory
  ! -----------------------------------------------
  call mpas_pool_create_pool( windowtraj )


  ! Add trajectory fields
  ! ---------------------
  call mpas_pool_get_field(fields % subFields, 'temperature', field2d_src)
  call mpas_duplicate_field(field2d_src, field2d_dst)
  call mpas_pool_add_field(windowtraj, 'temperature', field2d_dst)

  call mpas_pool_get_field(fields % subFields, 'spechum', field2d_src)
  call mpas_duplicate_field(field2d_src, field2d_dst)
  call mpas_pool_add_field(windowtraj, 'spechum', field2d_dst)

  call mpas_pool_get_field(fields % subFields, 'pressure', field2d_src)
  call mpas_duplicate_field(field2d_src, field2d_dst)
  call mpas_pool_add_field(windowtraj, 'pressure', field2d_dst)

  call mpas_pool_get_field(fields % subFields, 'surface_pressure', field1d_src)
  call mpas_duplicate_field(field1d_src, field1d_dst)
  call mpas_pool_add_field(windowtraj, 'surface_pressure', field1d_dst)


  ! Push to self%trajectories
  ! -------------------------
  call get_windowkey(t1, t2, windowkey)
  call mpas_pool_add_subpool(self%trajectories, windowkey, windowtraj)


  ! fill geovals at this state
  ! --------------------------
  ! NOTE: this call may get removed soon. See OOPS PR #577.
  call fill_geovals(self, geom, fields, t1, t2, locs, gom)

end subroutine set_trajectory

subroutine fill_geovals_tl(self, geom, fields, t1, t2, locs, gom)
  implicit none
  class(mpasjedi_lineargetvalues), intent(inout) :: self
  type(mpas_geom),                intent(in)    :: geom
  type(mpas_field),               intent(inout) :: fields
  type(datetime),                 intent(in)    :: t1
  type(datetime),                 intent(in)    :: t2
  type(ufo_locs),                 intent(in)    :: locs
  type(ufo_geovals),              intent(inout) :: gom

  character(len=*), parameter :: myname = 'fill_geovals_tl'

  logical, allocatable :: time_mask(:)
  integer :: jvar, jlev, ilev, jloc, ngrid, maxlevels, nlevels, nlocs, nlocsg
  real(kind=kind_real), allocatable :: obs_field(:,:)

  type (mpas_pool_type), pointer :: windowtraj  !< trajectory pool for t1 to t2
  character(len=41) :: windowkey

  type (mpas_pool_type), pointer :: pool_ufo
  type (mpas_pool_iterator_type) :: poolItr
  real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
  real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
  real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
  real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b

  logical :: allocateGeo
  character(len=1024) :: buf

  ! Get grid dimensions and checks
  ! ------------------------------
  ngrid = geom % nCellsSolve
  nlocs = locs % nlocs ! # of location for entire window

  !write(*,*)'fill_geovals_tl: ngrid, nlocs = : ',ngrid, nlocs
  !call interp_checks("tl", inc, locs, vars, gom)

  ! If no observations can early exit
  ! ---------------------------------
  call geom%f_comm%allreduce(nlocs,nlocsg,fckit_mpi_sum())
  if (nlocsg == 0) then
    return
  endif

  ! Get mask for locations in this time window
  ! ------------------------------------------
  call ufo_locs_time_mask(locs, t1, t2, time_mask)

  !write(0,*)'fill_geovals_tl: gom%nvar        : ',gom%nvar
  !write(0,*)'fill_geovals_tl: gom%variables   : ',gom%variables
  !write(0,*)'fill_geovals_tl: nlocs           : ',nlocs
  !write(0,*)'fill_geovals_tl: size(time_mask) : ',size(time_mask)

  ! Allocate intermediate pool of fields w/ geovals vars
  ! and TL of variable conversion
  ! ----------------------------------------------------
  call get_windowkey(t1, t2, windowkey)
  call mpas_pool_get_subpool(self%trajectories, windowkey, windowtraj)

  call convert_mpas_field2ufoTL(geom, windowtraj, fields % subFields, pool_ufo, gom%variables, gom%nvar, ngrid)


  ! Allocate and initialize output geovals
  ! --------------------------------------
  call mpas_pool_begin_iteration(pool_ufo)
  do while ( mpas_pool_get_next_member(pool_ufo, poolItr) )
    if (poolItr % memberType == MPAS_POOL_FIELD) then
      jvar = ufo_vars_getindex( gom%variables,trim(poolItr % memberName))
      if ( jvar < 1 ) cycle

      if (poolItr % nDims == 1) then
        gom%geovals(jvar)%nval = 1
      else if (poolItr % nDims == 2) then
        call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r2d_ptr_a)
        gom%geovals(jvar)%nval = size(r2d_ptr_a,1)
      endif

      allocateGeo = .true.
      if ( allocated(gom%geovals(jvar)%vals) ) then
        if ( size(gom%geovals(jvar)%vals,1) /= gom%geovals(jvar)%nval ) then
          deallocate(gom%geovals(jvar)%vals)
        else
          allocateGeo = .false.
        end if
      endif

      if ( allocateGeo ) then
        allocate( gom%geovals(jvar)%vals(gom%geovals(jvar)%nval,gom%geovals(jvar)%nlocs) )
        gom%geovals(jvar)%vals = MPAS_JEDI_ZERO_kr
      endif
    endif
  end do
  do jvar=1,gom%nvar
    if( .not. allocated(gom%geovals(jvar)%vals) )then
      ! air_pressure is required for GNSSRO, but mpas-jedi does not have corresponding TL/AD
      ! code in mpas2ufo_vars_mod.  Temporarily perform allocation here for excluded gom%variables.
      ! TODO: move this allocation for 2D variables into loop over mpas_pool members (above)
      !   Either
      !   (1) add TL/AD for air_pressure if added as analysis variable and move
      !       this allocation for 2D variables
      !    OR
      !   (2) implement increment variable change for mpas-jedi, including
      !       dependencies of UFO TL/AD on non-analyzed variables
      gom%geovals(jvar)%nval = geom%nVertLevels
      allocate( gom%geovals(jvar)%vals(gom%geovals(jvar)%nval,gom%geovals(jvar)%nlocs) )
      gom%geovals(jvar)%vals = MPAS_JEDI_ZERO_kr
      !write(buf,*) "--> fill_geovals_tl: geoval is not allocated, ", gom%variables(jvar)
      !call abor1_ftn(buf)
    endif
  enddo
  gom%linit = .true.


  ! TL of interpolate fields to obs locations using pre-calculated weights
  ! ----------------------------------------------------------------------
  maxlevels = geom%nVertLevelsP1
  allocate(obs_field(nlocs,maxlevels))

  call mpas_pool_begin_iteration(pool_ufo)
  do while ( mpas_pool_get_next_member(pool_ufo, poolItr) )
    if (poolItr % memberType == MPAS_POOL_FIELD) then
      jvar = ufo_vars_getindex( gom%variables,trim(poolItr % memberName))
      if ( jvar < 1 ) cycle

      !write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)

      obs_field = MPAS_JEDI_ZERO_kr

      nlevels = gom%geovals(jvar)%nval

      !TODO - JJG: Reduce wall-time of getvalues/_tl/_ad
      ! + apply_obsop takes ~99.9% of wall-time of getvalues_tl on cheyenne and
      !   scales with node count. Seems to have MPI-related issue.
      !
      if (poolItr % dataType == MPAS_POOL_REAL) then
        if (poolItr%nDims==1) then
          call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r1d_ptr_a)
          call self%bumpinterp%apply(r1d_ptr_a(1:ngrid), obs_field(:,1))
        else if (poolItr%nDims==2) then
          call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r2d_ptr_a)
          call self%bumpinterp%apply(r2d_ptr_a(1:nlevels,1:ngrid), obs_field(:,1:nlevels),trans_in=.true.)
        end if
      else
        call abor1_ftn('--> fill_geovals_tl: not a real field')
      end if

      do jlev = 1, nlevels
        ilev = nlevels - jlev + 1
        do jloc = 1, nlocs
          !BJJ-tmp vertical flip, top-to-bottom for CRTM geoval
          if (time_mask(jloc)) gom%geovals(jvar)%vals(ilev,jloc) = obs_field(jloc,jlev)
        end do
      end do
    endif
  end do
  deallocate(obs_field)

  !allocate(mod_field(ngrid,1))
  !allocate(obs_field(nlocs,1))

  !! Special cases go here

  !deallocate(mod_field)
  !deallocate(obs_field)

  call mpas_pool_destroy_pool(pool_ufo)

  !write(*,*) '---- Leaving fill_geovals_tl ---'

end subroutine fill_geovals_tl

subroutine fill_geovals_ad(self, geom, fields, t1, t2, locs, gom)
  implicit none
  class(mpasjedi_lineargetvalues), intent(inout) :: self
  type(mpas_geom),                intent(in)     :: geom
  type(mpas_field),               intent(inout)  :: fields
  type(datetime),                 intent(in)     :: t1
  type(datetime),                 intent(in)     :: t2
  type(ufo_locs),                 intent(in)     :: locs
  type(ufo_geovals),              intent(in)     :: gom

  character(len=*), parameter :: myname = 'fill_geovals_ad'

  logical, allocatable :: time_mask(:)
  integer :: jvar, jlev, ilev, jloc, ngrid, maxlevels, nlevels, nlocs, nlocsg
  real(kind=kind_real), allocatable :: obs_field(:,:)

  type (mpas_pool_type), pointer :: windowtraj  !< trajectory pool for t1 to t2
  character(len=41) :: windowkey

  type (mpas_pool_type), pointer :: pool_ufo
  type (mpas_pool_iterator_type) :: poolItr
  real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
  real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
  real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
  real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b

  ! If no observations can early exit
  ! ---------------------------------
  !if (self%traj%noobs)  return


  ! Get mask for locations in this time window
  ! ------------------------------------------
  call ufo_locs_time_mask(locs, t1, t2, time_mask)


  ! Get grid dimensions and checks
  ! ------------------------------
  ngrid = geom % nCellsSolve
  nlocs = locs % nlocs ! # of location for entire window
  !write(0,*)'getvalues_ad: nlocs        : ',nlocs

  ! If no observations can early exit
  ! ---------------------------------
  call geom%f_comm%allreduce(nlocs,nlocsg,fckit_mpi_sum())
  if (nlocsg == 0) then
    return
  endif

  ! Allocate intermediate pool of fields w/ geovals vars
  ! ----------------------------------------------------
  call get_windowkey(t1, t2, windowkey)
  call mpas_pool_get_subpool(self%trajectories, windowkey, windowtraj)

  !NOTE: this TL routine creates "pool_ufo". The field values from TL routine are unused.
  !    : The fields are initialized as "zero" in following "do while" loop.
  call convert_mpas_field2ufoTL(geom, windowtraj, fields % subFields, pool_ufo, gom%variables, gom%nvar, ngrid)


  ! Adjoint of interpolate fields to obs locations using pre-calculated weights
  ! ---------------------------------------------------------------------------
  maxlevels = geom%nVertLevelsP1
  allocate(obs_field(nlocs,maxlevels))

  call mpas_pool_begin_iteration(pool_ufo)
  do while ( mpas_pool_get_next_member(pool_ufo, poolItr) )
    if (poolItr % memberType == MPAS_POOL_FIELD) then
      jvar = ufo_vars_getindex( gom%variables,trim(poolItr % memberName))
      if ( jvar < 1 ) cycle

      obs_field = MPAS_JEDI_ZERO_kr

      !write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)
      nlevels = gom%geovals(jvar)%nval
      do jlev = 1, nlevels
        !ORG- obs_field(:,jlev) = gom%geovals(jvar)%vals(jlev,:)
        !BJJ-tmp vertical flip, top-to-bottom for CRTM geoval
        ilev = nlevels - jlev + 1
        do jloc = 1, nlocs
          if (time_mask(jloc)) then
            obs_field(jloc,jlev) = gom%geovals(jvar)%vals(ilev, jloc)
          endif
        end do
      end do

      !TODO - JJG: Reduce wall-time of getvalues/_tl/_ad
      ! + apply_obsop takes ~99.9% of wall-time of getvalues_ad on cheyenne and
      !   scales with node count. Seems to have MPI-related issue.
      !
      ! Copy data
      if (poolItr % dataType == MPAS_POOL_REAL) then
        if (poolItr%nDims==1) then
          call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r1d_ptr_a)
          r1d_ptr_a = MPAS_JEDI_ZERO_kr
          call self%bumpinterp%apply_ad(obs_field(:,1),r1d_ptr_a(1:ngrid))
        else if (poolItr%nDims==2) then
          call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r2d_ptr_a)
          r2d_ptr_a = MPAS_JEDI_ZERO_kr
          call self%bumpinterp%apply_ad(obs_field,r2d_ptr_a(1:nlevels,1:ngrid),trans_in=.true.)
        end if
      else
        call abor1_ftn('--> fill_geovals_ad: not a real field')
      end if
    endif
  end do
  deallocate(obs_field)

  ! AD of variable conversion
  ! -------------------------
  call convert_mpas_field2ufoAD(geom, windowtraj, fields % subFields, pool_ufo, gom%variables, gom%nvar, ngrid)

  !allocate(mod_field(ngrid,1))
  !allocate(obs_field(nlocs,1))

  !! Special cases go here

  !deallocate(mod_field)
  !deallocate(obs_field)

  call mpas_pool_destroy_pool(pool_ufo)

  !write(*,*) '---- Leaving fill_geovals_ad ---'

end subroutine fill_geovals_ad

subroutine get_windowkey(t1, t2, key)
  implicit none
  type(datetime),    intent(in)  :: t1, t2
  character(len=41), intent(out) :: key

  character(len=20) :: t1str, t2str

  call datetime_to_string(t1, t1str)
  call datetime_to_string(t2, t2str)
  key = trim(t1str)//'_'//trim(t2str)
end subroutine get_windowkey

end module mpasjedi_lineargetvalues_mod
