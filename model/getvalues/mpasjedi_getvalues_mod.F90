! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module mpasjedi_getvalues_mod

! fckit
use fckit_mpi_module,              only: fckit_mpi_comm
use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum

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
use mpas_dmpar, only: mpas_dmpar_exch_halo_field


!mpas-jedi
use mpas_constants_mod
use mpas_geom_mod
use mpas_field_utils_mod
use mpas2ufo_vars_mod
use mpas4da_mod

! --------------------------------------------------------------------------------------------------

implicit none
private
public :: mpasjedi_getvalues, mpasjedi_getvalues_base
public :: mpas_getvalues_registry
public :: fill_geovals, initialize_bump

type, abstract :: mpasjedi_getvalues_base
  type(bump_type) :: bump
  contains
    procedure, public :: initialize_bump
    procedure, public :: fill_geovals
end type mpasjedi_getvalues_base

type, extends(mpasjedi_getvalues_base) :: mpasjedi_getvalues
  contains
    procedure, public :: create
    procedure, public :: delete
end type mpasjedi_getvalues

#define LISTED_TYPE mpasjedi_getvalues

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_getvalues_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine create(self, geom, locs)
  implicit none
  class(mpasjedi_getvalues),      intent(inout) :: self
  type(mpas_geom),                intent(in)    :: geom
  type(ufo_locs),                 intent(in)    :: locs

  call initialize_bump(self, geom, locs)

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

  class(mpasjedi_getvalues), intent(inout) :: self

  call self%bump%dealloc()

end subroutine delete

! --------------------------------------------------------------------------------------------------

! Non-linear fill_geovals used in both GetValues and LinearGetValues classes
subroutine fill_geovals(self, geom, fields, t1, t2, locs, gom)

  class(mpasjedi_getvalues_base), intent(inout) :: self
  type(mpas_geom),                intent(in)    :: geom
  type(mpas_field),               intent(in)    :: fields
  type(datetime),                 intent(in)    :: t1
  type(datetime),                 intent(in)    :: t2
  type(ufo_locs),                 intent(in)    :: locs
  type(ufo_geovals),              intent(inout) :: gom

  logical, allocatable :: time_mask(:)
  type(fckit_mpi_comm)     :: f_comm
  integer :: jj, jvar, jlev, ilev, jloc, ngrid, maxlevels, nlevels, ivar, nlocs, nlocsg
  real(kind=kind_real), allocatable :: mod_field(:,:), mod_field_ext(:,:)
  real(kind=kind_real), allocatable :: obs_field(:,:)
  real(kind=kind_real), allocatable :: tmp_field(:,:)  !< for wspeed/wdir

  type (mpas_pool_type), pointer :: pool_ufo  !< pool with ufo variables
  type (mpas_pool_iterator_type) :: poolItr
  real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
  real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
  integer, dimension(:), pointer :: i1d_ptr_a, i1d_ptr_b
  integer, dimension(:,:), pointer :: i2d_ptr_a
  integer, allocatable :: index_nn(:)
  real (kind=kind_real), allocatable :: weight_nn(:)

  real(kind=kind_real) :: wdir           !< for wind direction
  integer :: ivarw, ivarl, ivari, ivars  !< for sfc fraction indices

  logical :: allocateGeo
  character(len=1024) :: buf

  ! Get grid dimensions and checks
  ! ------------------------------
  ngrid = geom % nCellsSolve
  nlocs = locs % nlocs ! # of location for entire window

  ! If no observations can early exit
  ! ---------------------------------
  f_comm = geom%f_comm
  call f_comm%allreduce(nlocs,nlocsg,fckit_mpi_sum())
  if (nlocsg == 0) then
    return
  endif

  !TODO: implement interp_checks if still necessary
  !call interp_checks("nl", fields, locs, gom%variables, gom)

  ! Get mask for locations in this time window
  ! ------------------------------------------
  call ufo_locs_time_mask(locs, t1, t2, time_mask)

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
          mod_field(:,1:nlevels) = real( transpose(i2d_ptr_a(1:nlevels,1:ngrid)), kind_real )
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

      !TODO - JJG: Reduce wall-time of getvalues/_tl/_ad
      ! + apply_obsop takes ~50% of wall-time of getvalues on cheyenne and
      !   scales with node count. Seems to have MPI-related issue.
      !
      ! + initialize_bump takes other ~50% on cheyennne
      !
      self%bump%geom%nl0 = nlevels
      call self%bump%apply_obsop(mod_field(:,1:nlevels),obs_field(:,1:nlevels))

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

  self%bump%geom%nl0 = 1
  allocate(mod_field(ngrid,1))
  allocate(obs_field(nlocs,1))


  ! Special cases --> interpolation, then conversion (wspeed, wdir)
  ! -------------------------------------------------
  !---add special cases: var_sfc_wspeed and/or var_sfc_wdir
  if ( (ufo_vars_getindex(gom%variables,var_sfc_wspeed) > 0) .or. &
       (ufo_vars_getindex(gom%variables,var_sfc_wdir)   > 0) ) then

    !write(*,*) ' BJJ: special cases: var_sfc_wspeed and/or var_sfc_wdir'

    !- allocate
    allocate(tmp_field(nlocs,2))

    !- read/interp.
    call mpas_pool_get_array(fields % subFields, "u10", r1d_ptr_a)
    mod_field(:,1) = r1d_ptr_a(1:ngrid)
    !write(*,*) 'MIN/MAX of u10=',minval(mod_field(:,1)),maxval(mod_field(:,1))
    call self%bump%apply_obsop(mod_field,obs_field)
    tmp_field(:,1)=obs_field(:,1)
    call mpas_pool_get_array(fields % subFields, "v10", r1d_ptr_a)
    mod_field(:,1) = r1d_ptr_a(1:ngrid)
    !write(*,*) 'MIN/MAX of v10=',minval(mod_field(:,1)),maxval(mod_field(:,1))
    call self%bump%apply_obsop(mod_field,obs_field)
    tmp_field(:,2)=obs_field(:,1)

    !- allocate geoval & put values for var_sfc_wspeed
    ivar = ufo_vars_getindex(gom%variables,var_sfc_wspeed)
    if(ivar > 0) then
      if( .not. allocated(gom%geovals(ivar)%vals) )then
        gom%geovals(ivar)%nval = 1
        allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
      endif
      do jloc = 1, nlocs
        if (time_mask(jloc)) then
          gom%geovals(ivar)%vals(1,jloc) = sqrt( tmp_field(jloc,1)**2 + tmp_field(jloc,2)**2 ) ! ws = sqrt(u**2+v**2) [m/s]
        endif
      end do
      !write(*,*) 'MIN/MAX of ',trim(var_sfc_wspeed),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
    endif

    !- allocate geoval & put values for var_sfc_wdir
    ivar = ufo_vars_getindex(gom%variables,var_sfc_wdir)
    if(ivar > 0) then
      if( .not. allocated(gom%geovals(ivar)%vals) )then
        gom%geovals(ivar)%nval = 1
        allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
      endif
      do jloc = 1, nlocs
        call uv_to_wdir(tmp_field(jloc,1), tmp_field(jloc,2), wdir) ! uu, vv, wind10_direction in radian
        if (time_mask(jloc)) gom%geovals(ivar)%vals(1,jloc) = wdir * MPAS_JEDI_RAD2DEG_kr  ! radian -> degree
      enddo
      !write(*,*) 'MIN/MAX of ',trim(var_sfc_wdir),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
    endif

    !- deallocate
    deallocate(tmp_field)

  endif  !---end special cases


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


    !initialize vector of nearest neighbor indices
    allocate( index_nn(nlocs) )
    allocate( weight_nn(self%bump%obsop%h%n_s) )

    do jloc = 1, nlocs
      !Picks index of self%bump%obsop%h%S containing maxium weight for obs jloc
      !Generic method for any interpolation scheme
      weight_nn = MPAS_JEDI_ZERO_kr
      where ( self%bump%obsop%h%row .eq. jloc )
        weight_nn = self%bump%obsop%h%S
      end where
      jj = maxloc(weight_nn,1)

  !       !Cheaper method that works for BUMP unstructured "triangular mesh" ( 3 vertices per obs ) with Bilinear interp.
  !       jj=3*(jloc-1) + maxloc(self%bump%obsop%h%S( 3*(jloc-1)+1:3*(jloc-1)+3 ),1) !nearest-interp. / maximum-weight specified.

      !Store index of BUMP extended vector
      index_nn(jloc) = self%bump%obsop%h%col(jj)
    enddo

    deallocate(weight_nn)

    !- allocate geoval & put values for var_sfc_landtyp
    ivar = ufo_vars_getindex(gom%variables,var_sfc_landtyp)
    if(ivar > 0) then
      if( .not. allocated(gom%geovals(ivar)%vals) )then
        gom%geovals(ivar)%nval = 1
        allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
      endif
      mod_field(:,1) = real( i1d_ptr_a(1:ngrid), kind_real)
      allocate( mod_field_ext(self%bump%obsop%nc0b,1) )
      call self%bump%obsop%com%ext(self%bump%mpl,1,mod_field,mod_field_ext)
      do jloc = 1, nlocs
        if (time_mask(jloc)) gom%geovals(ivar)%vals(1,jloc) = mod_field_ext( index_nn(jloc), 1 )
      enddo
      deallocate( mod_field_ext )
      !write(*,*) 'MIN/MAX of ',trim(var_sfc_landtyp),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
    endif

    !- allocate geoval & put values for var_sfc_vegtyp
    ivar = ufo_vars_getindex(gom%variables,var_sfc_vegtyp)
    if(ivar > 0) then
      if( .not. allocated(gom%geovals(ivar)%vals) )then
        gom%geovals(ivar)%nval = 1
        allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
      endif
      mod_field(:,1) = real( i1d_ptr_a(1:ngrid), kind_real)
      allocate( mod_field_ext(self%bump%obsop%nc0b,1) )
      call self%bump%obsop%com%ext(self%bump%mpl,1,mod_field,mod_field_ext)
      do jloc = 1, nlocs
        if (time_mask(jloc)) then
          gom%geovals(ivar)%vals(1,jloc) = real(convert_type_veg( int(mod_field_ext( index_nn(jloc), 1 ))), kind_real)
        endif
      enddo
      deallocate( mod_field_ext )
      !write(*,*) 'MIN/MAX of ',trim(var_sfc_vegtyp),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
    endif

    !- allocate geoval & put values for var_sfc_soiltyp
    ivar = ufo_vars_getindex(gom%variables,var_sfc_soiltyp)
    if(ivar > 0) then
      if( .not. allocated(gom%geovals(ivar)%vals) )then
        gom%geovals(ivar)%nval = 1
        allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
      endif
      mod_field(:,1) = real( i1d_ptr_b(1:ngrid), kind_real)
      allocate( mod_field_ext(self%bump%obsop%nc0b,1) )
      call self%bump%obsop%com%ext(self%bump%mpl,1,mod_field,mod_field_ext)
      do jloc = 1, nlocs
        if (time_mask(jloc)) then
          gom%geovals(ivar)%vals(1,jloc) = real(convert_type_soil( int(mod_field_ext( index_nn(jloc), 1 ))), kind_real)
        endif
      enddo
      deallocate( mod_field_ext )
      !write(*,*) 'MIN/MAX of ',trim(var_sfc_soiltyp),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
    endif

    deallocate(index_nn)

  endif  !---end special cases


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
      mod_field(:,1) = real(i1d_ptr_a(1:ngrid))
      call self%bump%apply_obsop(mod_field,obs_field)
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
      mod_field(:,1) = r1d_ptr_a(1:ngrid)
      call self%bump%apply_obsop(mod_field,obs_field)
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
      mod_field(:,1) = r1d_ptr_b(1:ngrid)
      call self%bump%apply_obsop(mod_field,obs_field)
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

  call mpas_pool_destroy_pool(pool_ufo)

end subroutine fill_geovals

! --------------------------------------------------------------------------------------------------

subroutine initialize_bump(self, grid, locs)

  implicit none
  class(mpasjedi_getvalues_base), intent(inout) :: self
  type(mpas_geom),          intent(in)  :: grid
  type(ufo_locs),           intent(in)  :: locs

  integer, save                     :: bumpid = 1000
  integer                           :: mod_nz,mod_num
  real(kind=kind_real), allocatable :: mod_lat(:), mod_lon(:)
  real(kind=kind_real), allocatable :: area(:),vunit(:,:)
  logical, allocatable              :: lmask(:,:)

  character(len=5)   :: cbumpcount
  character(len=255) :: bump_nam_prefix

  type(fckit_mpi_comm) :: f_comm

  ! Each bump%nam%prefix must be distinct
  ! -------------------------------------
  bumpid = bumpid + 1
  write(cbumpcount,"(I0.5)") bumpid
  bump_nam_prefix = 'mpas_bump_data_'//cbumpcount

  ! Get the Solution dimensions
  ! ---------------------------
  mod_num = grid%nCellsSolve
  !write(*,*)'initialize_bump mod_num,obs_num = ', mod_num, locs%nlocs

  !Calculate interpolation weight using BUMP
  !------------------------------------------
  allocate( mod_lat(mod_num) )
  allocate( mod_lon(mod_num) )
  mod_lat(:) = grid%latCell( 1:mod_num ) * MPAS_JEDI_RAD2DEG_kr !- to Degrees
  mod_lon(:) = grid%lonCell( 1:mod_num ) * MPAS_JEDI_RAD2DEG_kr !- to Degrees

  ! Namelist options
  ! ----------------

  !Important namelist options
  f_comm = grid%f_comm
  call self%bump%nam%init(f_comm%size())

  !Less important namelist options (should not be changed)
  self%bump%nam%prefix       = trim(bump_nam_prefix)  ! Prefix for files output
  self%bump%nam%default_seed = .true.
  self%bump%nam%new_obsop    = .true.
  self%bump%nam%write_obsop  = .false.
  self%bump%nam%verbosity    = "none"

  ! Initialize geometry
  ! -------------------
  allocate(area(mod_num))
  allocate(vunit(mod_num,1))
  allocate(lmask(mod_num,1))

  area  = MPAS_JEDI_ONE_kr          ! Dummy area, unit [m^2]
  vunit = MPAS_JEDI_ONE_kr          ! Dummy vertical unit
  lmask = .true.       ! Mask

  ! Initialize BUMP
  ! ---------------
  call self%bump%setup_online(f_comm,mod_num,1,1,1,mod_lon,mod_lat,area,vunit,lmask, &
                              nobs=locs%nlocs,lonobs=locs%lon(:),latobs=locs%lat(:))

  ! Run BUMP drivers
  call self%bump%run_drivers

  ! Partial deallocate option
  call self%bump%partial_dealloc

  ! Release memory
  ! --------------
  deallocate(area)
  deallocate(vunit)
  deallocate(lmask)
  deallocate(mod_lat)
  deallocate(mod_lon)

end subroutine initialize_bump

end module mpasjedi_getvalues_mod
