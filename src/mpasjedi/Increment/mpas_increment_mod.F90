! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_increment_mod

use atlas_module, only: atlas_geometry
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log

!oops
use datetime_mod
use kinds, only: kind_real
use oops_variables_mod

!ufo
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
use mpas_fields_mod
use mpas2ufo_vars_mod
use mpas4da_mod

implicit none

private

public :: diff_incr, dirac

integer, parameter    :: max_string=8000
character(max_string) :: message

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine diff_incr(lhs,x1,x2)

   implicit none
   class(mpas_fields), intent(inout) :: lhs
   class(mpas_fields), intent(in)    :: x1
   class(mpas_fields), intent(in)    :: x2
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

subroutine dirac(self, f_conf)

   implicit none
   class(mpas_fields),        intent(inout) :: self
   type(fckit_configuration), intent(in)    :: f_conf   !< Configuration
   character(len=:), allocatable :: str
   integer                :: ndir, idir, ildir, ndirlocal
   character(len=3)       :: idirchar
   character(len=StrKIND) :: dirvar
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=RKIND), dimension(:,:), pointer :: r2d_ptr_a
   real (kind=RKIND), dimension(:), pointer   :: r1d_ptr_a
   integer :: nearestCell
   integer, allocatable, dimension(:) :: dirOwned, dirOwnedGlobal
   real (kind=kind_real), allocatable, dimension(:) :: dirLats
   real (kind=kind_real), allocatable, dimension(:) :: dirLons
   integer, allocatable, dimension(:) :: dirCells
   real (kind=RKIND) :: x1, y1

   ! Get number and positions of Diracs
   call f_conf%get_or_die("ndir",ndir)

   allocate( dirOwned(ndir) )
   allocate( dirCells(ndir) )

   call f_conf%get_or_die("dirLats",dirLats)
   call f_conf%get_or_die("dirLons",dirLons)

   call f_conf%get_or_die("ildir",ildir)
   call f_conf%get_or_die("dirvar",str)
   dirvar = str

   !Test if dir is owned and find the nearest local cell
   ! (repurposed from MPAS-Release/src/core_atmosphere/diagnostics/soundings.F)

   ndirlocal = 0
   do idir=1,ndir
      nearestCell = self % geom % nCellsSolve
      x1=real(dirLats(idir) * MPAS_JEDI_DEG2RAD_kr, kind=RKIND)
      y1=real(dirLons(idir) * MPAS_JEDI_DEG2RAD_kr, kind=RKIND)
      nearestCell = nearest_cell( x1, y1, &
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

   write(message,*) 'This processor owns ',ndirlocal,' dirac forcing locations'
   call fckit_log%debug(message)

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
         write(message,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , poolItr % memberName
         call fckit_log%debug(message)

         if (poolItr % dataType == MPAS_POOL_REAL) then
            ! Depending on the dimensionality of the field, we need to set pointers of
            ! the correct type
            if( trim(dirvar) .eq. trim(poolItr % memberName) ) then
               ndirlocal = 0
               do idir=1, ndir
                  if ( dirOwned(idir).eq.1 ) then
                     if (poolItr % nDims == 1) then
                        call mpas_pool_get_array(self % subFields, trim(poolItr % memberName), r1d_ptr_a)
                        r1d_ptr_a( dirCells(idir) ) = MPAS_JEDI_ONE_kr
                     else if (poolItr % nDims == 2) then
                        call mpas_pool_get_array(self % subFields, trim(poolItr % memberName), r2d_ptr_a)
                        r2d_ptr_a( ildir, dirCells(idir) ) = MPAS_JEDI_ONE_kr
                     else if (poolItr % nDims == 3) then
                        call fckit_log%info ('Not implemented yet')
                     end if
                     ndirlocal = ndirlocal + 1
                  end if
               end do
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

   real (kind=RKIND), intent(in) :: target_lat, target_lon
   integer, intent(in) :: start_cell
   integer, intent(in) :: nCells, maxEdges
   integer, dimension(nCells), intent(in) :: nEdgesOnCell
   integer, dimension(maxEdges,nCells), intent(in) :: cellsOnCell
   real (kind=RKIND), dimension(nCells), intent(in) :: latCell, lonCell

   integer :: i
   integer :: iCell
   integer :: current_cell
   real (kind=kind_real) :: current_distance, d
   real (kind=kind_real) :: nearest_distance
   type (atlas_geometry) :: ageometry
   
   real (kind=kind_real) :: x1, y1, x2, y2

   nearest_cell = start_cell
   current_cell = -1

   ageometry = atlas_geometry("UnitSphere")

   do while (nearest_cell /= current_cell)
      current_cell = nearest_cell
!      current_distance = ageometry%distance(lonCell(current_cell), latCell(current_cell),  target_lon, target_lat)
      x1=lonCell(current_cell); y1=latCell(current_cell); x2=target_lon; y2=target_lat
      current_distance = ageometry%distance(x1,y1,x2,y2)
      nearest_cell = current_cell
      nearest_distance = current_distance
      do i = 1, nEdgesOnCell(current_cell)
         iCell = cellsOnCell(i,current_cell)
         if (iCell <= nCells) then
            x1=lonCell(iCell); y1=latCell(iCell); x2=target_lon; y2=target_lat
            !d = ageometry%distance(lonCell(iCell), latCell(iCell), target_lon, target_lat)
            d = ageometry%distance(x1,y1,x2,y2)
            if (d < nearest_distance) then
               nearest_cell = iCell
               nearest_distance = d
            end if
         end if
      end do
   end do
end function nearest_cell

! ------------------------------------------------------------------------------

end module mpas_increment_mod
