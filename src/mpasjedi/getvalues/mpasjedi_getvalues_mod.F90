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
use mpas_fields_mod
use mpas2ufo_vars_mod
use mpas4da_mod

! --------------------------------------------------------------------------------------------------

implicit none
private
public :: mpasjedi_getvalues, mpasjedi_getvalues_base
public :: mpas_getvalues_registry
public :: fill_geovals, getvalues_base_create, getvalues_base_delete

type, abstract :: mpasjedi_getvalues_base
logical        :: use_bump_interp
type(bump_interpolator) :: bumpinterp
type(unstrc_interp)     :: unsinterp
contains
    procedure, private :: initialize_uns_interp
    procedure, public  :: fill_geovals
    procedure, private :: integer_interpolation_bump
    procedure, public  :: getvalues_base_create
    procedure, public  :: getvalues_base_delete
end type mpasjedi_getvalues_base

type, extends(mpasjedi_getvalues_base) :: mpasjedi_getvalues
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
  character (len=1024) :: buf

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
        write(buf,*) '--> getvalues_base_create: interpolation type: ',interp_type,' not implemented'
        call abor1_ftn(buf)
    end select
  else
    self%use_bump_interp = .True. ! BUMP is default interpolation
  end if

  if (self%use_bump_interp) then
    call self%bumpinterp%init(geom%f_comm, afunctionspace_in=geom%afunctionspace, lon_out=lons, lat_out=lats, &
      & nl=geom%nVertLevels)
  else
    call initialize_uns_interp(self, geom, lats, lons)
  endif

  if (allocated(interp_type)) deallocate(interp_type)
  deallocate(lons, lats)

end subroutine getvalues_base_create

! --------------------------------------------------------------------------------------------------

! -----------------------------------------------------------------------------
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

! -----------------------------------------------------------------------------
!> \brief GetValues class 'create' logic
!!
!! \details **create** This subroutine populates a functional getvalues
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

!> \brief GetValues class 'delete' logic
!!
!! \details **delete** This subroutine deletes (frees memory) for a getvalues
!! class instance.
subroutine delete(self)

  class(mpasjedi_getvalues), intent(inout) :: self !< getvalues self
  call getvalues_base_delete(self)

end subroutine delete

! --------------------------------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!> \brief Interpolates and converts from mpas_fields to populate ufo_geovals
!!
!! \details **fill_geovals** This subroutine populates the variables in a 
!! ufo_geovals object by interpolating and converting the state or increment
!! variables in an mpas_fields object. This is the non-linear subroutine used in
!! both GetValues and LinearGetValues classes
subroutine fill_geovals(self, geom, fields, t1, t2, locs, gom)
  implicit none
  class(mpasjedi_getvalues_base), intent(inout) :: self    !< getvalues_base self
  type(mpas_geom),                intent(in)    :: geom    !< geometry (mpas mesh)
  type(mpas_fields),              intent(in)    :: fields  !< state or increment
  type(datetime),                 intent(in)    :: t1      !< time window begin
  type(datetime),                 intent(in)    :: t2      !< time window end
  type(ufo_locations),            intent(in)    :: locs    !< observation locations
  type(ufo_geovals),              intent(inout) :: gom     !< geovals

  logical(c_bool), allocatable :: time_mask(:)
  integer :: jj, jvar, jlev, ilev, jloc, ngrid, maxlevels, nlevels, ivar, nlocs, nlocsg
  integer, allocatable ::obs_field_int(:,:)
  real(kind=kind_real), allocatable :: mod_field(:,:)
  real(kind=kind_real), allocatable :: obs_field(:,:)

  type (mpas_pool_type), pointer :: pool_ufo  !< pool with ufo variables
  type (mpas_pool_iterator_type) :: poolItr
  real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
  real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
  integer, dimension(:), pointer :: i1d_ptr_a, i1d_ptr_b
  integer, dimension(:,:), pointer :: i2d_ptr_a

  integer :: ivarw, ivarl, ivari, ivars  !< for sfc fraction indices

  logical :: allocateGeo
  character(len=1024) :: buf

  ! Get grid dimensions and checks
  ! ------------------------------
  ngrid = geom % nCellsSolve
  nlocs = locs % nlocs() ! # of location for entire window

  ! If no observations can early exit
  ! ---------------------------------
  call geom%f_comm%allreduce(nlocs,nlocsg,fckit_mpi_sum())
  if (nlocsg == 0) then
    return
  endif

  !TODO: implement interp_checks if still necessary
  !call interp_checks("nl", fields, locs, gom%variables, gom)

  ! Get mask for locations in this time window
  ! ------------------------------------------
  allocate(time_mask(nlocs))
  call locs%get_timemask(t1, t2, time_mask)

  !write(0,*)'fill_geovals: nlocs           : ',nlocs
  !write(0,*)'fill_geovals: size(time_mask) : ',size(time_mask)

  ! Allocate intermediate pool of fields w/ geovals vars
  ! and NL variable conversion
  ! ---------------------------------------------------------
  call convert_mpas_field2ufo(geom, fields % subFields, pool_ufo, gom%variables, gom%nvar, ngrid)


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
        if (poolItr % dataType == MPAS_POOL_INTEGER) then
          call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), i2d_ptr_a)
          gom%geovals(jvar)%nval = size(i2d_ptr_a,1)
        else if (poolItr % dataType == MPAS_POOL_REAL) then
          call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r2d_ptr_a)
          gom%geovals(jvar)%nval = size(r2d_ptr_a,1)
        endif
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
  gom%linit = .true.


  ! Interpolate fields to obs locations using pre-calculated weights
  ! ----------------------------------------------------------------
  maxlevels = geom%nVertLevelsP1
  allocate(mod_field(ngrid,maxlevels))
  allocate(obs_field(nlocs,maxlevels))

  call mpas_pool_begin_iteration(pool_ufo)
  do while ( mpas_pool_get_next_member(pool_ufo, poolItr) )
    if (poolItr % memberType == MPAS_POOL_FIELD) then
      ivar = ufo_vars_getindex( gom%variables,trim(poolItr % memberName))
      if ( ivar == -1 ) cycle

      !write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)
      nlevels = gom%geovals(ivar)%nval
      if (poolItr % nDims == 1) then
        if (poolItr % dataType == MPAS_POOL_INTEGER) then
          call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), i1d_ptr_a)
          mod_field(:,1) = real( i1d_ptr_a(1:ngrid), kind_real)
        else if (poolItr % dataType == MPAS_POOL_REAL) then
          call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r1d_ptr_a)
          mod_field(:,1) = r1d_ptr_a(1:ngrid)
        endif
      else if (poolItr % nDims == 2) then
        if (poolItr % dataType == MPAS_POOL_INTEGER) then
          call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), i2d_ptr_a)
          mod_field(:,1:nlevels) = real( transpose (i2d_ptr_a(1:nlevels,1:ngrid)), kind_real )
        else if (poolItr % dataType == MPAS_POOL_REAL) then
          call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r2d_ptr_a)
          mod_field(:,1:nlevels) = transpose(r2d_ptr_a(1:nlevels,1:ngrid))
        endif
      else
        write(buf,*) '--> fill_geovals: poolItr % nDims == ',poolItr % nDims,' not handled'
        call abor1_ftn(buf)
      endif
      !write(*,*) "interp: var, ufo_var_index = ",trim(poolItr % memberName), ivar
      !write(*,*) 'MIN/MAX of ',trim(poolItr % memberName), &
      !           minval(mod_field(:,1:nlevels)), &
      !           maxval(mod_field(:,1:nlevels))

      if (self%use_bump_interp) then
        call self%bumpinterp%apply(mod_field(1:ngrid,1:nlevels), obs_field(:,1:nlevels),trans_in=.false.)
      else
        do jlev = 1, nlevels
          call self%unsinterp%apply(mod_field(:,jlev),obs_field(:,jlev))
        end do
      endif

      do jlev = 1, nlevels
        !BJJ-tmp vertical flip, top-to-bottom for CRTM geoval
        ilev = nlevels - jlev + 1
        do jloc = 1, nlocs
          if (time_mask(jloc)) gom%geovals(ivar)%vals(ilev, jloc) = obs_field(jloc,jlev)
        end do
      end do
      !write(*,*) 'MIN/MAX of ',trim(poolItr % memberName),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)

    endif
  end do !- end of pool iteration
  deallocate(mod_field)
  deallocate(obs_field)

  allocate(mod_field(ngrid,1))
  allocate(obs_field(nlocs,1))
  allocate(obs_field_int(nlocs,1))

  ! Special cases --> nearest neighbor (surface integer values)
  ! -----------------------------------------------------------
  !---add special cases: var_sfc_landtyp, var_sfc_vegtyp, var_sfc_soiltyp
  if ( (ufo_vars_getindex(gom%variables,var_sfc_landtyp) > 0) .or. &
       (ufo_vars_getindex(gom%variables,var_sfc_vegtyp)  > 0) .or. &
       (ufo_vars_getindex(gom%variables,var_sfc_soiltyp) > 0) ) then

    !write(*,*) ' BJJ: special cases: var_sfc_landtyp, var_sfc_vegtyp, or var_sfc_soiltyp'

    call mpas_pool_get_array(fields % subFields, "ivgtyp", i1d_ptr_a)
    !write(*,*) 'MIN/MAX of ivgtyp=',minval(i1d_ptr_a),maxval(i1d_ptr_a)

    call mpas_pool_get_array(fields % subFields, "isltyp", i1d_ptr_b)
    !write(*,*) 'MIN/MAX of isltyp=',minval(i1d_ptr_b),maxval(i1d_ptr_b)

    if (self%use_bump_interp) then
      !- allocate geoval & put values for var_sfc_landtyp, var_sfc_vegtyp, var_sfc_soiltyp
      call integer_interpolation_bump(self, var_sfc_landtyp, ngrid, nlocs, time_mask, i1d_ptr_a, obs_field_int, gom)
      call integer_interpolation_bump(self, var_sfc_vegtyp,  ngrid, nlocs, time_mask, i1d_ptr_a, obs_field_int, gom)
      call integer_interpolation_bump(self, var_sfc_soiltyp, ngrid, nlocs, time_mask, i1d_ptr_b, obs_field_int, gom)
    else
      call integer_interpolation_unsinterp(self, var_sfc_landtyp, ngrid, nlocs, time_mask, i1d_ptr_a, mod_field, gom)
      call integer_interpolation_unsinterp(self, var_sfc_vegtyp,  ngrid, nlocs, time_mask, i1d_ptr_a, mod_field, gom)
      call integer_interpolation_unsinterp(self, var_sfc_soiltyp, ngrid, nlocs, time_mask, i1d_ptr_b, mod_field, gom)
    endif
  endif  !---end surface integer values

  ! Special cases --> interdependent geovals (surface fractions)
  ! ------------------------------------------------------------
  !---add special cases: var_sfc_wfrac, var_sfc_lfrac, var_sfc_ifrac, var_sfc_sfrac
  !---    simple interpolation now, but can be more complex: Consider FOV ??
  if ( (ufo_vars_getindex(gom%variables,var_sfc_wfrac) > 0) .or. &
       (ufo_vars_getindex(gom%variables,var_sfc_lfrac) > 0) .or. &
       (ufo_vars_getindex(gom%variables,var_sfc_ifrac) > 0) .or. &
       (ufo_vars_getindex(gom%variables,var_sfc_sfrac) > 0) ) then

    !write(*,*) ' BJJ: special cases: var_sfc_wfrac, var_sfc_lfrac, var_sfc_ifrac, or var_sfc_sfrac'

    call mpas_pool_get_array(fields % subFields, "landmask", i1d_ptr_a) !"land-ocean mask (1=land ; 0=ocean)"
    !write(*,*) 'MIN/MAX of landmask=',minval(i1d_ptr_a),maxval(i1d_ptr_a)
    call mpas_pool_get_array(fields % subFields, "xice", r1d_ptr_a)     !"fractional area coverage of sea-ice"
    !write(*,*) 'MIN/MAX of xice=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
    call mpas_pool_get_array(fields % subFields, "snowc", r1d_ptr_b)    !"flag for snow on ground (=0 no snow; =1,otherwise"
    !write(*,*) 'MIN/MAX of snowc=',minval(r1d_ptr_b),maxval(r1d_ptr_b)
    !write(*,*) 'MIN/MAX of lnad+xice+snowc=', &
  !                minval(real(i1d_ptr_a)+r1d_ptr_a+r1d_ptr_b),maxval(real(i1d_ptr_a)+r1d_ptr_a+r1d_ptr_b)

    ivarw = ufo_vars_getindex(gom%variables,var_sfc_wfrac)
    ivarl = ufo_vars_getindex(gom%variables,var_sfc_lfrac)
    ivari = ufo_vars_getindex(gom%variables,var_sfc_ifrac)
    ivars = ufo_vars_getindex(gom%variables,var_sfc_sfrac)

    !--- Land first. will be adjusted later
    ivar = ivarl
    if(ivar > 0) then
      if( .not. allocated(gom%geovals(ivar)%vals) )then
        gom%geovals(ivar)%nval = 1
        allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
      endif
      if (self%use_bump_interp) then
        call self%bumpinterp%apply(real(i1d_ptr_a(1:ngrid),kind_real),obs_field)
      else
        mod_field(:,1) = r1d_ptr_a(1:ngrid)
        call self%unsinterp%apply(mod_field,obs_field)
      endif
      do jloc = 1, nlocs
        if (time_mask(jloc)) gom%geovals(ivar)%vals(1,jloc) = obs_field(jloc,1)
      enddo
      !write(*,*) 'MIN/MAX of ',trim(var_sfc_lfrac),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
    endif
    !--- determine ICE
    ivar = ivari
    if(ivar > 0) then
      if( .not. allocated(gom%geovals(ivar)%vals) )then
        gom%geovals(ivar)%nval = 1
        allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
      endif
      if (self%use_bump_interp) then
        call self%bumpinterp%apply(r1d_ptr_a(1:ngrid), obs_field)
      else
        mod_field(:,1) = r1d_ptr_a(1:ngrid)
        call self%unsinterp%apply(mod_field,obs_field)
      endif
      do jloc = 1, nlocs
        if (time_mask(jloc)) gom%geovals(ivar)%vals(1,jloc) = obs_field(jloc,1)
      enddo
      !write(*,*) 'MIN/MAX of ',trim(var_sfc_ifrac),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
    endif
    !--- detemine/adjust SNOW & SEA
    ivar = ivars
    if(ivar > 0) then
      if( .not. allocated(gom%geovals(ivar)%vals) )then
        gom%geovals(ivar)%nval = 1
        allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
      endif
      if( .not. allocated(gom%geovals(ivarw)%vals) )then
        gom%geovals(ivarw)%nval = 1
        allocate( gom%geovals(ivarw)%vals(gom%geovals(ivarw)%nval,gom%geovals(ivarw)%nlocs) )
      endif
      if (self%use_bump_interp) then
        call self%bumpinterp%apply(r1d_ptr_b(1:ngrid), obs_field)
      else
        mod_field(:,1) = r1d_ptr_b(1:ngrid)
        call self%unsinterp%apply(mod_field,obs_field)
      endif
      do jloc = 1, nlocs
        if (time_mask(jloc)) gom%geovals(ivar)%vals(1,jloc) = obs_field(jloc,1)
      enddo
      do jloc = 1, nlocs
        if (time_mask(jloc)) then
          if (gom%geovals(ivari)%vals(1,jloc) > MPAS_JEDI_ZERO_kr) then
            gom%geovals(ivar)%vals(1,jloc)  = &
              max(min(gom%geovals(ivar)%vals(1,jloc), &
                      MPAS_JEDI_ONE_kr - &
                      gom%geovals(ivari)%vals(1,jloc)), &
                  MPAS_JEDI_ZERO_kr)
            gom%geovals(ivarw)%vals(1,jloc) = &
              max(MPAS_JEDI_ONE_kr - &
                  gom%geovals(ivari)%vals(1,jloc) - &
                  gom%geovals(ivar)%vals(1,jloc), &
                  MPAS_JEDI_ZERO_kr)
          else if (gom%geovals(ivar)%vals(1,jloc) > MPAS_JEDI_ZERO_kr) then
            gom%geovals(ivarw)%vals(1,jloc) = &
              max(MPAS_JEDI_ONE_kr - &
                  gom%geovals(ivar)%vals(1,jloc), &
                  MPAS_JEDI_ZERO_kr)
            gom%geovals(ivari)%vals(1,jloc) = MPAS_JEDI_ZERO_kr
          else
            gom%geovals(ivarw)%vals(1,jloc) = &
              max(MPAS_JEDI_ONE_kr - &
                  gom%geovals(ivarl)%vals(1,jloc), &
                  MPAS_JEDI_ZERO_kr)
            gom%geovals(ivari)%vals(1,jloc) = MPAS_JEDI_ZERO_kr
            gom%geovals(ivar)%vals(1,jloc)  = MPAS_JEDI_ZERO_kr
          endif
        endif
      enddo
      !write(*,*) 'MIN/MAX of ',trim(var_sfc_sfrac),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
      !write(*,*) 'MIN/MAX of ',trim(var_sfc_wfrac),minval(gom%geovals(ivarw)%vals),maxval(gom%geovals(ivarw)%vals)
    endif
    !--- Final adjust LAND
    ivar = ivarl
    if(ivar > 0) then
      do jloc = 1, nlocs
        if (time_mask(jloc)) then
          gom%geovals(ivar)%vals(1,jloc) = &
            max(MPAS_JEDI_ONE_kr - &
                gom%geovals(ivarw)%vals(1,jloc) - &
                gom%geovals(ivari)%vals(1,jloc) - &
                gom%geovals(ivars)%vals(1,jloc), &
                MPAS_JEDI_ZERO_kr)
        endif
      enddo
      !write(*,*) 'MIN/MAX of ',trim(var_sfc_lfrac),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
    endif
    !do ii=17,19  !1,nlocs
    !write(*,*) gom%geovals(ivarl)%vals(1,ii), gom%geovals(ivarw)%vals(1,ii), &
    !           gom%geovals(ivari)%vals(1,ii), gom%geovals(ivars)%vals(1,ii)
    !enddo

    !--- OMG: adjust between sfc coverage & sfc type
    !         see wrfda/da_get_innov_vector_crtm.inc#L521
    if ( (ufo_vars_getindex(gom%variables,var_sfc_landtyp) > 0) .or. &
         (ufo_vars_getindex(gom%variables,var_sfc_vegtyp)  > 0) .or. &
         (ufo_vars_getindex(gom%variables,var_sfc_soiltyp) > 0) ) then
      do jloc = 1, nlocs
        if (time_mask(jloc)) then
          if(gom%geovals(ivarl)%vals(1,jloc) > MPAS_JEDI_ZERO_kr) then
            if(nint(gom%geovals(ufo_vars_getindex(gom%variables,var_sfc_soiltyp))%vals(1,jloc)) .eq. 9 .or. &
              nint(gom%geovals(ufo_vars_getindex(gom%variables,var_sfc_vegtyp))%vals(1,jloc)) .eq. 13 ) then
              gom%geovals(ivari)%vals(1,jloc) = min( gom%geovals(ivari)%vals(1,jloc) + &
                                                   gom%geovals(ivarl)%vals(1,jloc), &
                                                   MPAS_JEDI_ONE_kr )
              gom%geovals(ivarl)%vals(1,jloc) = MPAS_JEDI_ZERO_kr
            endif
          endif
        endif
      enddo
    endif  !--- OMG: end
  endif  !---end special cases

  ! Deallocate local memory
  ! -----------------------
  deallocate(mod_field)
  deallocate(obs_field)
  deallocate(time_mask)

  call mpas_pool_destroy_pool(pool_ufo)

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

  integer :: nn, ngrid_in, ngrid_out
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
subroutine integer_interpolation_bump(self, varname, ngrid, nlocs, time_mask, &
                                      data_in, obs_field_int, gom)

  implicit none

  class(mpasjedi_getvalues_base), intent(inout) :: self     !< self
  character(len=*), intent(in)             :: varname            !< name of interp variable
  integer, intent(in)                      :: ngrid              !< number of cells in model mesh
  integer, intent(in)                      :: nlocs              !< number of locations for obs
  logical(c_bool), allocatable, intent(in) :: time_mask(:)       !< mask for time window
  integer, dimension(:), intent(in)        :: data_in            !< data to interpolate
  integer, allocatable, intent(inout)      :: obs_field_int(:,:) !< output array of interpolated data
  type(ufo_geovals), intent(inout)         :: gom                !< output geoVaLs

  integer :: jloc, ivar

  ! This code assumes time_mask already allocated to size nlocs and poplulated earlier.

  !- allocate geoval & put values for 'varname'
  ivar = ufo_vars_getindex(gom%variables, varname)
  if(ivar > 0) then
    if( .not. allocated(gom%geovals(ivar)%vals) )then
      gom%geovals(ivar)%nval = 1
      allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
    endif
    call self%bumpinterp%apply(data_in(1:ngrid),obs_field_int)
    do jloc = 1, nlocs
      if (time_mask(jloc)) then
        select case (varname)
          case (var_sfc_landtyp)
            gom%geovals(ivar)%vals(1,jloc) = real(obs_field_int(jloc,1), kind_real)

          case (var_sfc_vegtyp)
            gom%geovals(ivar)%vals(1,jloc) = &
              real(convert_type_veg(obs_field_int(jloc, 1)), kind_real)

          case (var_sfc_soiltyp)
            gom%geovals(ivar)%vals(1,jloc) = &
              real(convert_type_soil(obs_field_int(jloc, 1)), kind_real)

          case default
            call abor1_ftn('integer_interpolation_bump: unimplemented varname')
        end select
      endif
    enddo
    !write(*,*) 'MIN/MAX of ',trim(varname),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
  endif

end subroutine integer_interpolation_bump

! ------------------------------------------------------------------------------

!> \brief Performs interpolation of integer fields using unstructured interpolation
!!
!! \details **integer_interpolation_unsinterp** This subroutine performs the interpolation
!! of integer-valued fields (i.e. types) using unstructured interpolation
subroutine integer_interpolation_unsinterp(self, varname, ngrid, nlocs, time_mask, &
                                          data_in, mod_field, gom)

  implicit none

  class(mpasjedi_getvalues_base), intent(inout)    :: self           !< self
  character(len=*), intent(in)                     :: varname        !< name of interp variable
  integer, intent(in)                              :: ngrid          !< number of cells in model mesh
  integer, intent(in)                              :: nlocs          !< number of locations for obs
  logical(c_bool), allocatable, intent(in)         :: time_mask(:)   !< mask for time window
  integer, dimension(:), intent(in)                :: data_in        !< data to interpolate
  real(kind=kind_real), allocatable, intent(inout) :: mod_field(:,:) !< pre-allocated 2-d (ngrid,1) array
  type(ufo_geovals), intent(inout)                 :: gom            !< output geoVaLs

  integer :: jloc, ivar
  real(kind=kind_real), dimension(nlocs) :: interpolated_data

  ! This code assumes time_mask already allocated to size nlocs and poplulated earlier.
  ! Also, that mod_field already allocated to size (ngrid, 1)

  !- allocate geoval & put values for 'varname'
  ivar = ufo_vars_getindex(gom%variables, varname)
  if(ivar > 0) then
    if( .not. allocated(gom%geovals(ivar)%vals) )then
      gom%geovals(ivar)%nval = 1
      allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
    endif
    mod_field(:,1) = real( data_in(1:ngrid), kind_real)
    call unsinterp_integer_apply(self%unsinterp, mod_field(:,1), interpolated_data)
    do jloc = 1, nlocs
      if (time_mask(jloc)) then
        select case (varname)
          case (var_sfc_landtyp)
            gom%geovals(ivar)%vals(1,jloc) = interpolated_data(jloc)

          case (var_sfc_vegtyp)
            gom%geovals(ivar)%vals(1,jloc) = &
              real(convert_type_veg( int(interpolated_data(jloc))), kind_real)

          case (var_sfc_soiltyp)
            gom%geovals(ivar)%vals(1,jloc) = &
              real(convert_type_soil( int(interpolated_data(jloc))), kind_real)

          case default
            call abor1_ftn('integer_interpolation_unsinterp: unimplemented varname')
        end select
      endif
    enddo
    !write(*,*) 'MIN/MAX of ',trim(varname),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
  endif

end subroutine integer_interpolation_unsinterp

end module mpasjedi_getvalues_mod
