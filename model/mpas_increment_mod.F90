! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_increment_mod

use iso_c_binding 

!oops
use config_mod
use datetime_mod
use kinds, only: kind_real
use variables_mod

!ufo
use ufo_locs_mod
use ufo_geovals_mod
use ufo_vars_mod

!MPAS-Model
use mpas_constants
use mpas_dmpar
use mpas_derived_types
use mpas_kind_types, only: StrKIND
use mpas_pool_routines

!mpas-jedi
use mpas_constants_mod
use mpas_geom_mod
use mpas_getvaltraj_mod, only: mpas_getvaltraj
use mpas_field_utils_mod
use mpas_increment_utils_mod
use mpas_state_utils_mod, only: mpas_state
use mpas2ufo_vars_mod
use mpas4da_mod

implicit none

private

public :: diff_incr, dirac, &
        & getvalues_tl, getvalues_ad, &
        & ug_coord, increment_to_ug, increment_from_ug

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine diff_incr(lhs,x1,x2)

   implicit none
   class(mpas_increment), intent(inout) :: lhs
   class(mpas_state),     intent(in)    :: x1
   class(mpas_state),     intent(in)    :: x2
   character(len=StrKIND) :: kind_op

   call lhs%zeros()
   if (x1%geom%nCells==x2%geom%nCells .and. x1%geom%nVertLevels==x2%geom%nVertLevels) then
     if (lhs%geom%nCells==x1%geom%nCells .and. lhs%geom%nVertLevels==x1%geom%nVertLevels) then
        kind_op = 'sub'
        call da_operator(trim(kind_op), lhs % subFields, x1 % subFields, x2 % subFields, fld_select = lhs % fldnames_ci)
     else
       call abor1_ftn("mpas_increment:diff_incr: dimension mismatch between the two variables.")
     endif
   else
     call abor1_ftn("mpas_increment:diff_incr: states not at same resolution")
   endif

   return

end subroutine diff_incr

! ------------------------------------------------------------------------------

subroutine dirac(self, c_conf)

   implicit none
   class(mpas_increment), intent(inout) :: self
   type(c_ptr),           intent(in)    :: c_conf   !< Configuration
   integer                :: ndir, idir, ildir, ndirlocal
   character(len=3)       :: idirchar
   character(len=StrKIND) :: dirvar
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a
   integer :: nearestCell
   integer, allocatable, dimension(:) :: dirOwned, dirOwnedGlobal
   real (kind=kind_real), allocatable, dimension(:) :: dirLats
   real (kind=kind_real), allocatable, dimension(:) :: dirLons
   integer, allocatable, dimension(:) :: dirCells

   ! Get number and positions of Diracs
   ndir = config_get_int(c_conf,"ndir")

   allocate( dirOwned(ndir) )
   allocate( dirLats(ndir) )
   allocate( dirLons(ndir) )
   allocate( dirCells(ndir) )

   do idir=1,ndir
      write(idirchar,'(i3)') idir
      dirLats(idir) = config_get_real(c_conf,"dirLats("//trim(adjustl(idirchar))//")")
      dirLons(idir) = config_get_real(c_conf,"dirLons("//trim(adjustl(idirchar))//")")
   end do
   ildir = config_get_int(c_conf,"ildir")
   dirvar = config_get_string(c_conf,len(dirvar),"dirvar")

   !Test if dir is owned and find the nearest local cell
   ! (repurposed from MPAS-Release/src/core_atmosphere/diagnostics/soundings.F)

   ndirlocal = 0
   do idir=1,ndir
      nearestCell = self % geom % nCellsSolve
      nearestCell = nearest_cell( (dirLats(idir) * deg2rad), &
                                  (dirLons(idir) * deg2rad), &
                                  nearestCell, self % geom % nCells, self % geom % maxEdges, &
                                  self % geom % nEdgesOnCell, self % geom % cellsOnCell, &
                                  self % geom % latCell, self % geom % lonCell )

      if (nearestCell <= self % geom % nCellsSolve) then
          dirOwned(idir) = 1
          dirCells(idir) = nearestCell
          ndirlocal = ndirlocal + 1
      else
          dirOwned(idir) = 0
          dirCells(idir) = self % geom % nCells + 1
      end if
   end do

   write(*,*) ' This processor owns ',ndirlocal, &
              ' dirac forcing locations'

   ! Check
   if (ndir<1) call abor1_ftn("mpas_increment:dirac non-positive ndir")

   allocate( dirOwnedGlobal(ndir) )
   call mpas_dmpar_max_int_array( self % geom % domain % dminfo, ndir, dirOwned, dirOwnedGlobal)
   if ( any(dirOwnedGlobal.lt.1) ) then
         call abor1_ftn("mpas_increment:dirac invalid Lat/Lon")
   end if

   call mpas_dmpar_sum_int_array( self % geom % domain % dminfo, ndir, dirOwned, dirOwnedGlobal)
   if ( any(dirOwnedGlobal.gt.1) ) then
         call abor1_ftn("mpas_increment:duplicated dirac on >1 processors")
   end if
   deallocate( dirOwnedGlobal )

  if ((ildir < 1) .or. (ildir > self % geom % nVertLevels)) then
      call abor1_ftn("mpas_increment:dirac invalid ildir")
   endif

   ! Setup Diracs
   call self%zeros()

   call mpas_pool_begin_iteration(self % subFields)

   do while ( mpas_pool_get_next_member(self % subFields, poolItr) )
      ! Pools may in general contain dimensions, namelist options, fields, or other pools,
      ! so we select only those members of the pool that are fields
      if (poolItr % memberType == MPAS_POOL_FIELD) then
         ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
         write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , poolItr % memberName

         if (poolItr % dataType == MPAS_POOL_REAL) then
            ! Depending on the dimensionality of the field, we need to set pointers of
            ! the correct type
            if (poolItr % nDims == 1) then
               write(*,*)'Not implemented yet'
            else if (poolItr % nDims == 2) then
               call mpas_pool_get_array(self % subFields, trim(poolItr % memberName), r2d_ptr_a)
               if( trim(dirvar) .eq. trim(poolItr % memberName) ) then
                 ndirlocal = 0
                 do idir=1, ndir
                    if ( dirOwned(idir).eq.1 ) then
                       r2d_ptr_a( ildir, dirCells(idir) ) = 1.0_kind_real
                       ndirlocal = ndirlocal + 1
                    end if
                 end do
                 write(*,*) ' Dirac is set in ',ndirlocal,'locations for',trim(poolItr % memberName)
               end if
            else if (poolItr % nDims == 3) then
               write(*,*)'Not implemented yet'
            end if
         end if
      end if
   end do

   deallocate( dirOwned )
   deallocate( dirLats )
   deallocate( dirLons )
   deallocate( dirCells )

end subroutine dirac

! ------------------------------------------------------------------------------

!!TODO: Alternatively could make the function nearest_cell public in MPAS
!!      and then define an interface to it within this module (cleaner)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finds the MPAS grid cell nearest to (target_lat, target_lon)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function nearest_cell(target_lat, target_lon, start_cell, nCells, maxEdges, &
                              nEdgesOnCell, cellsOnCell, latCell, lonCell)

   implicit none

   real (kind=kind_real), intent(in) :: target_lat, target_lon
   integer, intent(in) :: start_cell
   integer, intent(in) :: nCells, maxEdges
   integer, dimension(nCells), intent(in) :: nEdgesOnCell
   integer, dimension(maxEdges,nCells), intent(in) :: cellsOnCell
   real (kind=kind_real), dimension(nCells), intent(in) :: latCell, lonCell

   integer :: i
   integer :: iCell
   integer :: current_cell
   real (kind=kind_real) :: current_distance, d
   real (kind=kind_real) :: nearest_distance

   nearest_cell = start_cell
   current_cell = -1

   do while (nearest_cell /= current_cell)
      current_cell = nearest_cell
      current_distance = sphere_distance(latCell(current_cell), lonCell(current_cell), target_lat, &
                                         target_lon, 1.0_kind_real)
      nearest_cell = current_cell
      nearest_distance = current_distance
      do i = 1, nEdgesOnCell(current_cell)
         iCell = cellsOnCell(i,current_cell)
         if (iCell <= nCells) then
            d = sphere_distance(latCell(iCell), lonCell(iCell), target_lat, target_lon, 1.0_kind_real)
            if (d < nearest_distance) then
               nearest_cell = iCell
               nearest_distance = d
            end if
         end if
      end do
   end do
end function nearest_cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) 
!    on a sphere with given radius.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (kind=kind_real) function sphere_distance(lat1, lon1, lat2, lon2, radius)

   implicit none

   real (kind=kind_real), intent(in) :: lat1, lon1, lat2, lon2, radius
   real (kind=kind_real) :: arg1

   arg1 = sqrt( sin(0.5*(lat2-lat1))**2 + cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
   sphere_distance = 2.0 * radius * asin(arg1)

end function sphere_distance

! ------------------------------------------------------------------------------

subroutine getvalues_tl(inc, locs, vars, gom, traj)

   implicit none
   class(mpas_increment), intent(inout) :: inc
   type(ufo_locs),        intent(in)    :: locs
   type(oops_vars),        intent(in)    :: vars
   type(ufo_geovals),     intent(inout) :: gom
   type(mpas_getvaltraj), intent(inout) :: traj

   character(len=*), parameter :: myname = 'getvalues_tl'
   
   integer :: jvar, jlev, ilev, jloc, ngrid, maxlevels, nlevels, nlocs, ivar
   real(kind=kind_real), allocatable :: mod_field(:,:)
   real(kind=kind_real), allocatable :: obs_field(:,:)
   
   type (mpas_pool_type), pointer :: pool_ufo
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   character(len=MAXVARLEN) :: ufo_var_name

   ! Check traj is implemented
   ! -------------------------
   if (.not.traj%lalloc) &
   call abor1_ftn(trim(myname)//" trajectory for this obs op not found")

   ! If no observations can early exit
   ! ---------------------------------
   if (traj%noobs)  return

   ! Get grid dimensions and checks
   ! ------------------------------
   ngrid = inc%geom%nCellsSolve !or traj%ngrid
   nlocs = locs%nlocs           !or traj%nobs
!   write(*,*)'getvalues_tl: ngrid, nlocs = : ',ngrid, nlocs
   call interp_checks("tl", inc, locs, vars, gom)
   
   !Make sure the return values are allocated and set
   !-------------------------------------------------
   do jvar=1,vars%nv
      if( .not. allocated(gom%geovals(jvar)%vals) )then
         gom%geovals(jvar)%nval = inc%geom%nVertLevels
         allocate( gom%geovals(jvar)%vals(gom%geovals(jvar)%nval,gom%geovals(jvar)%nlocs) )
         gom%geovals(jvar)%vals = 0.0_kind_real
!         write(*,*) ' gom%geovals(n)%vals allocated'
         gom%linit = .true.
      endif
   enddo
     
   !Interpolate fields to obs locations using pre-calculated weights
   !----------------------------------------------------------------
   write(0,*)'getvalues_tl: vars%nv       : ',vars%nv
   write(0,*)'getvalues_tl: vars%fldnames : ',vars%fldnames
   write(0,*)'getvalues_tl: nlocs         : ',nlocs

   !------- need some table matching UFO_Vars & related MPAS_Vars
   !------- for example, Tv @ UFO may require Theta, Pressure, Qv.
   !-------                               or  Theta_m, exner_base, Pressure_base, Scalar(&index_qv)
   call convert_mpas_field2ufoTL(inc % geom, traj % pool_traj, inc % subFields, pool_ufo, vars % fldnames, vars % nv, ngrid) !--pool_ufo is new pool with ufo_vars

   maxlevels = inc%geom%nVertLevelsP1
   allocate(mod_field(ngrid,maxlevels))
   allocate(obs_field(nlocs,maxlevels))

   call mpas_pool_begin_iteration(pool_ufo)
   do while ( mpas_pool_get_next_member(pool_ufo, poolItr) )
      if (poolItr % memberType == MPAS_POOL_FIELD) then
         ufo_var_name = trim(poolItr % memberName)
         ivar = str_match(ufo_var_name, vars%fldnames)
         if ( ivar == -1 ) cycle

!         write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)

         nlevels = gom%geovals(ivar)%nval
         if (poolItr % nDims == 1) then

         else if (poolItr % nDims == 2) then

            if (poolItr % dataType == MPAS_POOL_REAL) then
               call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r2d_ptr_a)
!              write(*,*) "getvalues_tl: var=",trim(poolItr % memberName)
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
         traj%bump%geom%nl0 = nlevels
         call traj%bump%apply_obsop(mod_field(:,1:nlevels),obs_field(:,1:nlevels))

         do jlev = 1, nlevels
            ilev = nlevels - jlev + 1
            do jloc = 1, nlocs
               !BJJ-tmp vertical flip, top-to-bottom for CRTM geoval
               ! only selected obs (using locs%indx()) are filling "geovals"
               gom%geovals(ivar)%vals(ilev,locs%indx(jloc)) = obs_field(jloc,jlev)
            end do
         end do

      end if
   end do
   deallocate(mod_field)
   deallocate(obs_field)

   traj%bump%geom%nl0 = 1
!   allocate(mod_field(ngrid,1))
!   allocate(obs_field(nlocs,1))

!    ! Special cases go here

!   deallocate(mod_field)
!   deallocate(obs_field)

   call mpas_pool_destroy_pool(pool_ufo)

!   write(*,*) '---- Leaving getvalues_tl ---'
end subroutine getvalues_tl

! ------------------------------------------------------------------------------

subroutine getvalues_ad(inc, locs, vars, gom, traj)

   implicit none
   class(mpas_increment), intent(inout) :: inc
   type(ufo_locs),        intent(in)    :: locs
   type(oops_vars),        intent(in)    :: vars
   type(ufo_geovals),     intent(inout) :: gom
   type(mpas_getvaltraj), intent(inout) :: traj

   character(len=*), parameter :: myname = 'getvalues_ad'

   integer :: jvar, jlev, ilev, jloc, ngrid, maxlevels, nlevels, nlocs, ivar
   real(kind=kind_real), allocatable :: mod_field(:,:)
   real(kind=kind_real), allocatable :: obs_field(:,:)

   type (mpas_pool_type), pointer :: pool_ufo
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   character(len=MAXVARLEN) :: ufo_var_name

   ! Check traj is implemented
   ! -------------------------
   if (.not.traj%lalloc) &
   call abor1_ftn(trim(myname)//" trajectory for this obs op not found")

   ! If no observations can early exit
   ! ---------------------------------
   if (traj%noobs)  return

   ! Get grid dimensions and checks
   ! ------------------------------
   ngrid = inc%geom%nCellsSolve !or traj%ngrid
   nlocs = locs%nlocs           !or traj%nlocs

   !TODO: make this work for "ad"
   !call interp_checks("ad", inc, locs, vars, gom)

   !Interpolate fields to obs locations using pre-calculated weights
   !----------------------------------------------------------------
   write(0,*)'getvalues_ad: vars%nv       : ',vars%nv
   write(0,*)'getvalues_ad: vars%fldnames : ',vars%fldnames
   write(0,*)'getvalues_ad: nlocs         : ',nlocs

   !NOTE: This TL routine is called JUST to create "pool_ufo". Their values from TL routine doesn't matter.
   !    : Actually their values are initialized as "zero" in following "do while" loop.
   call convert_mpas_field2ufoTL(inc % geom, traj % pool_traj, inc % subFields, pool_ufo, vars % fldnames, vars % nv, ngrid) !--pool_ufo is new pool with ufo_vars

   maxlevels = inc%geom%nVertLevelsP1
   allocate(mod_field(ngrid,maxlevels))
   allocate(obs_field(nlocs,maxlevels))

   call mpas_pool_begin_iteration(pool_ufo)
   do while ( mpas_pool_get_next_member(pool_ufo, poolItr) )
      if (poolItr % memberType == MPAS_POOL_FIELD) then
         ufo_var_name = trim(poolItr % memberName)
         ivar = str_match(ufo_var_name, vars%fldnames)
         if ( ivar == -1 ) cycle

!         write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)
         nlevels = gom%geovals(ivar)%nval
         do jlev = 1, nlevels
            !ORG- obs_field(:,jlev) = gom%geovals(ivar)%vals(jlev,:)
            !BJJ-tmp vertical flip, top-to-bottom for CRTM geoval
            ! only selected obs (using locs%indx()) are filling "geovals"
            ilev = nlevels - jlev + 1
            do jloc = 1, nlocs
               obs_field(jloc,jlev) = gom%geovals(ivar)%vals(ilev, locs%indx(jloc))
               gom%geovals(ivar)%vals(ilev, locs%indx(jloc)) = 0.0_kind_real
            end do
         end do

         !TODO - JJG: Reduce wall-time of getvalues/_tl/_ad
         ! + apply_obsop takes ~99.9% of wall-time of getvalues_ad on cheyenne and
         !   scales with node count. Seems to have MPI-related issue.
         !
         traj%bump%geom%nl0 = nlevels
         call traj%bump%apply_obsop_ad(obs_field(:,1:nlevels),mod_field(:,1:nlevels))

         if (poolItr % nDims == 1) then

         else if (poolItr % nDims == 2) then

            if (poolItr % dataType == MPAS_POOL_REAL) then
               call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r2d_ptr_a)
               r2d_ptr_a=0.0_kind_real
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

   call convert_mpas_field2ufoAD(inc % geom, traj % pool_traj, inc % subFields, pool_ufo, vars % fldnames, vars % nv, ngrid) !--pool_ufo is new pool with ufo_vars

   traj%bump%geom%nl0 = 1
!   allocate(mod_field(ngrid,1))
!   allocate(obs_field(nlocs,1))

!    ! Special cases go here

!   deallocate(mod_field)
!   deallocate(obs_field)

   call mpas_pool_destroy_pool(pool_ufo)

!   write(*,*) '---- Leaving getvalues_ad ---' 
end subroutine getvalues_ad

! ------------------------------------------------------------------------------

subroutine ug_size(self, ug)

   use unstructured_grid_mod
   
   implicit none
   class(mpas_increment),   intent(in)    :: self
   type(unstructured_grid), intent(inout) :: ug
   integer :: igrid
   
   ! Set number of grids
   if (ug%colocated==1) then
      ! Colocatd
      ug%ngrid = 1
   else
      ! Not colocatedd
      ug%ngrid = 1
   end if

   ! Allocate grid instances
   if (.not.allocated(ug%grid)) allocate(ug%grid(ug%ngrid))

   if (ug%colocated==1) then ! colocated

      ! Set local number of points
      ug%grid(1)%nmga = self%geom%nCellsSolve

      ! Set number of levels
      ug%grid(1)%nl0 = self%geom%nVertLevels

      ! Set number of variables
      ug%grid(1)%nv = self%nf

      ! Set number of timeslots
      ug%grid(1)%nts = ug%nts

   else ! Not colocated

      do igrid=1,ug%ngrid
         ! Set local number of points
         ug%grid(igrid)%nmga = self%geom%nCellsSolve

         ! Set number of levels
         ug%grid(igrid)%nl0 = self%geom%nVertLevels

         ! Set number of variables
         ug%grid(igrid)%nv = self%nf

         ! Set number of timeslots
         ug%grid(igrid)%nts = ug%nts
      enddo
   end if
end subroutine ug_size

! ------------------------------------------------------------------------------

subroutine ug_coord(self, ug)

   use unstructured_grid_mod
   
   implicit none
   class(mpas_increment),   intent(in)    :: self
   type(unstructured_grid), intent(inout) :: ug
   
   integer :: jl, igrid
   
   ! Define size
   call ug_size(self, ug)

   ! Alocate unstructured grid coordinates
   call allocate_unstructured_grid_coord(ug)

   ! Copy coordinates
   if (ug%colocated==1) then ! colocated
     ug%grid(1)%lon = self%geom%lonCell(1:ug%grid(1)%nmga) / deg2rad !- to Degrees
     ug%grid(1)%lat = self%geom%latCell(1:ug%grid(1)%nmga) / deg2rad !- to Degrees
     ug%grid(1)%area = self%geom%areaCell(1:ug%grid(1)%nmga)
     do jl=1,self%geom%nVertLevels
       ug%grid(1)%vunit(:,jl) = real(jl,kind=kind_real)
       ug%grid(1)%lmask(:,jl) = .true.
     enddo
   else ! Not colocated
     do igrid=1,ug%ngrid
       ug%grid(igrid)%lon = self%geom%lonCell(1:ug%grid(igrid)%nmga) / deg2rad !- to Degrees
       ug%grid(igrid)%lat = self%geom%latCell(1:ug%grid(igrid)%nmga) / deg2rad !- to Degrees
       ug%grid(igrid)%area = self%geom%areaCell(1:ug%grid(igrid)%nmga)
       do jl=1,self%geom%nVertLevels
         ug%grid(igrid)%vunit(:,jl) = real(jl,kind=kind_real)
         ug%grid(igrid)%lmask(:,jl) = .true.
       enddo
     enddo
   endif

end subroutine ug_coord

! ------------------------------------------------------------------------------

subroutine increment_to_ug(self, ug, its)

   use mpas_pool_routines
   use unstructured_grid_mod
   
   implicit none
   class(mpas_increment),   intent(in)    :: self
   type(unstructured_grid), intent(inout) :: ug
   integer,                 intent(in)    :: its
   
   integer :: idx_var,jC,jl  
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   type(oops_vars) :: vars ! temporary to access variable "index" easily
   character(len=MAXVARLEN) :: ufo_var_name

   ! Set list of variables
   vars % nv = self % nf
   allocate(vars % fldnames(vars % nv))
   vars % fldnames(:) = self % fldnames(:)

   ! Define size
   call ug_size(self, ug)

   ! Allocate unstructured grid field
   call allocate_unstructured_grid_field(ug)

   ! Copy field
   call mpas_pool_begin_iteration(self % subFields)
   
   do while ( mpas_pool_get_next_member(self % subFields, poolItr) )
      ! Pools may in general contain dimensions, namelist options, fields, or other pools,
      ! so we select only those members of the pool that are fields
      if (poolItr % memberType == MPAS_POOL_FIELD) then
         idx_var = -999
         ufo_var_name = trim(poolItr % memberName)
         idx_var = str_match(ufo_var_name, vars%fldnames)
         write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)

         ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
         if (poolItr % dataType == MPAS_POOL_REAL) then
            ! Depending on the dimensionality of the field, we need to set pointers of
            ! the correct type
            if (poolItr % nDims == 1) then
               write(*,*)'Not implemented yet'
            else if (poolItr % nDims == 2) then
               call mpas_pool_get_array(self % subFields, trim(poolItr % memberName), r2d_ptr_a)
               if(idx_var.gt.0) then
!                  write(*,*) '  sub. increment_to_ug, poolItr % memberName=',trim(poolItr % memberName)
!                  write(*,*) '  sub. increment_to_ug, idx_var=',idx_var
                  do jC=1,self%geom%nCellsSolve
                    do jl=1,self%geom%nVertLevels
                      ug%grid(1)%fld(jC,jl,idx_var,its) = r2d_ptr_a(jl,jC)
                    enddo
                  enddo
               endif
            else if (poolItr % nDims == 3) then
               write(*,*)'Not implemented yet'
            end if
         end if
      end if
   end do

   ! Cleanup
   call oops_vars_delete(vars)

end subroutine increment_to_ug

! -----------------------------------------------------------------------------

subroutine increment_from_ug(self, ug, its)

   use mpas_pool_routines
   use unstructured_grid_mod

   implicit none
   class(mpas_increment),   intent(inout) :: self
   type(unstructured_grid), intent(in)    :: ug
   integer,                 intent(in)    :: its

   integer :: idx_var,jC,jl
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   type(oops_vars) :: vars ! temporary to access variable "index" easily
   character(len=MAXVARLEN) :: ufo_var_name

   ! Set list of variables
   vars % nv = self % nf
   allocate(vars % fldnames(vars % nv))
   vars % fldnames(:) = self % fldnames(:)

   ! Copy field
   call mpas_pool_begin_iteration(self % subFields)
   
   do while ( mpas_pool_get_next_member(self % subFields, poolItr) )
      ! Pools may in general contain dimensions, namelist options, fields, or other pools,
      ! so we select only those members of the pool that are fields
      if (poolItr % memberType == MPAS_POOL_FIELD) then
         ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
         write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)

         if (poolItr % dataType == MPAS_POOL_REAL) then
            ! Depending on the dimensionality of the field, we need to set pointers of
            ! the correct type
            idx_var = -999
            ufo_var_name = trim(poolItr % memberName)
            idx_var = str_match(ufo_var_name, vars%fldnames)

            if(idx_var.gt.0) then
               if (poolItr % nDims == 1) then
                  write(*,*)'Not implemented yet'
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_array(self % subFields, trim(poolItr % memberName), r2d_ptr_a)

!                  write(*,*) '  sub. convert_from_ug, poolItr % memberName=',trim(poolItr % memberName)
                  do jC=1,self%geom%nCellsSolve
                     do jl=1,self%geom%nVertLevels
                        r2d_ptr_a(jl,jC) = ug%grid(1)%fld(jC,jl,idx_var,its)
                     enddo
                  enddo
               else if (poolItr % nDims == 3) then
                  write(*,*)'Not implemented yet'
               end if
            end if
         end if
      end if
   end do

   ! Cleanup
   call oops_vars_delete(vars)

   ! TODO: Since only local locations are updated/transferred from ug, 
   !       need MPAS HALO comms before using these fields in MPAS

end subroutine increment_from_ug

! ------------------------------------------------------------------------------

end module mpas_increment_mod
