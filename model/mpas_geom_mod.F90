! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_geom_mod

use iso_c_binding
use config_mod
use netcdf
use tools_nc
use type_mpl, only: mpl

use mpas_derived_types
use mpas_kind_types
use mpas_constants
use mpas_kinds, only : kind_real

use mpas_subdriver
use atm_core
use mpas_pool_routines

use mpi ! only MPI_COMM_WORLD

implicit none
private
public :: mpas_geom, &
        & geo_setup, geo_clone, geo_delete, geo_info
public :: mpas_geom_registry


! ------------------------------------------------------------------------------

!> Fortran derived type to hold geometry definition
type :: mpas_geom
   integer :: nCellsGlobal !Global count
   integer :: nEdgesGlobal !Global count
   integer :: nVerticesGlobal !Global count
   integer :: nCells !Memory count
   integer :: nEdges !Memory count
   integer :: nVertices !Memory count
   integer :: nCellsSolve !Local count
   integer :: nEdgesSolve !Local count
   integer :: nVerticesSolve !Global count
   integer :: nVertLevels
   integer :: nVertLevelsP1
   integer :: nSoilLevels
   integer :: vertexDegree
   integer :: maxEdges
   character(len=StrKIND) :: gridfname
   real(kind=kind_real), DIMENSION(:),   ALLOCATABLE :: latCell, lonCell, xland
   real(kind=kind_real), DIMENSION(:),   ALLOCATABLE :: areaCell
   real(kind=kind_real), DIMENSION(:),   ALLOCATABLE :: latEdge, lonEdge
   real(kind=kind_real), DIMENSION(:,:), ALLOCATABLE :: edgeNormalVectors
   real(kind=kind_real), DIMENSION(:,:), ALLOCATABLE :: zgrid
!   integer, allocatable :: nEdgesOnCell(:)
!   integer, allocatable :: edgesOnCell(:,:)

   type (domain_type), pointer :: domain => null() 
   type (core_type), pointer :: corelist => null()
end type mpas_geom

#define LISTED_TYPE mpas_geom

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_geom_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------
subroutine geo_setup(self, c_conf)
implicit none
type(mpas_geom), intent(inout) :: self
type(c_ptr), intent(in) :: c_conf
character(len=StrKIND) :: string1
integer :: ncid, dimid, varid
real(kind=kind_real), parameter :: deg2rad = pii/180.0_kind_real

type (mpas_pool_type), pointer :: meshPool, fg
type (block_type), pointer :: block_ptr

real (kind=kind_real), pointer :: r1d_ptr(:), r2d_ptr(:,:)
integer, pointer :: i1d_ptr(:), i0d_ptr

write(*,*)' ==> create geom'
!> MPAS subdriver
call mpas_init(self % corelist, self % domain, mpi_comm=MPI_COMM_WORLD)
if (associated(self % domain)) then
    write(*,*)'inside geom: geom % domain associated'
end if
if (associated(self % corelist)) then
    write(*,*)'inside geom: geom % corelist associated'
else
    write(*,*)'inside geom: geom % corelist not associated'
end if

! The mesh info reading above could be replaced with accesses to the pool
!       block_ptr % structs (after mpas_init)
!  These pool accesses would refer to memory (local+halo) instead of the global sizes and arrays above,
!   requiring updates in mpas_fields_mod, mpas_model_mod, etc...  
!   That effort would streamline code for future parallelization.

!  The code below assumes only a single block is located on each processor.
!   (could remove "Mem" labels if these become the standard vars)

    block_ptr => self % domain % blocklist

    call mpas_pool_get_subpool ( block_ptr % structs, 'mesh', meshPool)

    call mpas_pool_get_dimension ( block_ptr % dimensions, 'nCells', i0d_ptr )         
    self % nCells = i0d_ptr
    call mpas_pool_get_dimension ( block_ptr % dimensions, 'nEdges', i0d_ptr )         
    self % nEdges = i0d_ptr
    call mpas_pool_get_dimension ( block_ptr % dimensions, 'nVertices', i0d_ptr )      
    self % nVertices = i0d_ptr
    call mpas_pool_get_dimension ( block_ptr % dimensions, 'nCellsSolve', i0d_ptr )    
    self % nCellsSolve = i0d_ptr
    call mpas_pool_get_dimension ( block_ptr % dimensions, 'nEdgesSolve', i0d_ptr )    
    self % nEdgesSolve = i0d_ptr
    call mpas_pool_get_dimension ( block_ptr % dimensions, 'nVerticesSolve', i0d_ptr ) 
    self % nVerticesSolve = i0d_ptr

    call mpas_pool_get_dimension ( block_ptr % dimensions, 'nVertLevels', i0d_ptr )    
    self % nVertLevels = i0d_ptr
    call mpas_pool_get_dimension ( block_ptr % dimensions, 'nVertLevelsP1', i0d_ptr )  
    self % nVertLevelsP1 = i0d_ptr
    call mpas_pool_get_dimension ( block_ptr % dimensions, 'nSoilLevels', i0d_ptr )    
    self % nSoilLevels = i0d_ptr
    call mpas_pool_get_dimension ( block_ptr % dimensions, 'vertexDegree', i0d_ptr )   
    self % vertexDegree = i0d_ptr
    call mpas_pool_get_dimension ( block_ptr % dimensions, 'maxEdges', i0d_ptr )       
    self % maxEdges = i0d_ptr


!  Could make this more flexible/clean by using pointers for array variables in mpas_geom 
!  + all later value modifications must be consistent with MPAS modeling
!  + should nullify instead of allocate/dellocate

    allocate ( self%latCell ( self%nCells ) )
    allocate ( self%lonCell ( self%nCells ) )
    allocate ( self%latEdge ( self%nEdges ) )
    allocate ( self%lonEdge ( self%nEdges ) )
    allocate ( self%xland ( self%nCells ) )
    allocate ( self%areaCell ( self%nCells ) )
    allocate ( self%edgeNormalVectors (3,  self%nEdges ) )
    allocate ( self%zgrid ( self%nVertLevelsP1, self%nCells ) )

    call mpas_pool_get_array ( meshPool, 'latCell', r1d_ptr )           
    self % latCell = r1d_ptr
    call mpas_pool_get_array ( meshPool, 'lonCell', r1d_ptr )           
    self % lonCell = r1d_ptr
    call mpas_pool_get_array ( meshPool, 'areaCell', r1d_ptr )          
    self % areaCell = r1d_ptr
    call mpas_pool_get_array ( meshPool, 'latEdge', r1d_ptr )           
    self % latEdge = r1d_ptr
    call mpas_pool_get_array ( meshPool, 'lonEdge', r1d_ptr )           
    self % lonEdge = r1d_ptr
    call mpas_pool_get_array ( meshPool, 'edgeNormalVectors', r2d_ptr ) 
    self % edgeNormalVectors = r2d_ptr ( 1:3,1:self%nEdges )
    call mpas_pool_get_array ( meshPool, 'zgrid', r2d_ptr )             
    self % zgrid = r2d_ptr ( 1:self%nVertLevelsP1,1:self%nCells )

!   !Not using these edges yet
!    allocate ( self%nEdgesOnCell ( self%nCells ) )
!    allocate ( self%edgesOnCell ( self%maxEdges,self%nCells ) )
!    call mpas_pool_get_array ( mesh, 'nEdgesOnCell', i1d_ptr )
!    self % nEdgesOnCell = i1d_ptr
!    call mpas_pool_get_array ( mesh, 'edgesOnCell', i1d_ptr )
!    self % edgesOnCell = i1d_ptr ( 1:self%maxEdges,1:self%nCells )

!XLAND MIGHT NEED TO BE DONE AFTER FG IS INITIALIZED...WHEN IS THAT?
!    call mpas_pool_get_subpool(block_ptr % structs, 'fg', fg)
!    call mpas_pool_get_array(fg, 'xland', r1d_ptr)
!    self % xland = r1d_ptr

!> radians to degrees
self%latCell = self%latCell !/ deg2rad
self%lonCell = self%lonCell !/ deg2rad
self%latEdge = self%latEdge !/ deg2rad
self%lonEdge = self%lonEdge !/ deg2rad


!> Open a grid mesh file
self%gridfname = config_get_string(c_conf, StrKIND, "gridfname")
string1 = self%gridfname
call ncerr(string1, nf90_open(trim(self%gridfname),nf90_nowrite,ncid))

!> Grid dimensions
call ncerr(string1, nf90_inq_dimid        (ncid,'nCells',dimid))
call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=self%nCellsGlobal))
call ncerr(string1, nf90_inq_dimid        (ncid,'nEdges',dimid))
call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=self%nEdgesGlobal))
call ncerr(string1, nf90_inq_dimid        (ncid,'nVertices',dimid))
call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=self%nVerticesGlobal))

!Still need to figure out where/when to get land from block_ptr%structs =>
allocate(r1d_ptr(self%nCellsGlobal))
call ncerr(string1, nf90_inq_varid(ncid,'xland',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,r1d_ptr))
call mpas_pool_get_array( block_ptr % allFields, 'indexToCellID', i1d_ptr)
self % xland = r1d_ptr(i1d_ptr(1:self%nCells))
deallocate(r1d_ptr)



!> close file
call ncerr(string1,nf90_close(ncid))


!call geo_to_local ( self )
 
write(*,*)'End of geo_setup'

end subroutine geo_setup

! ------------------------------------------------------------------------------

!subroutine geo_to_local ( self )
!   ! Description : This subroutine transfers creates index converters from
!   !               global and memory arrays to local arrays used by BUMP
!   !               It is required to loop through the local MPAS blocklist, each member
!   !               of which contains a subset of the local geometry (if nblock>1).
!   ! See the subroutine mpas_block_creator_finalize_block_phase1 in mpas_block_creator.F 
!   ! for more details.
!
!   implicit none
!
!   type (mpas_geom), intent(inout) :: self
!   type (block_type), pointer :: block_ptr
!   integer, pointer :: indexToCellID(:)
!   integer, allocatable :: CellsGlobalToHalo(:)
!   integer, pointer :: nCells_blk, nCellsSolve_blk
!   integer          :: nCellsHalo, CellsStart, CellsEnd, CellsSolveStart, CellsSolveEnd, ii, jj, nblock
!
!   if (.not. allocated(self % CellsGlobalToLocal)) then
!      block_ptr => self % domain % blocklist
!      nCellsHalo = 0
!      self % nCellsLocal = 0
!      nblock = 0
!      do while(associated(block_ptr))
!         nblock = nblock + 1
!         call mpas_pool_get_dimension(block_ptr % dimensions, 'nCellsSolve', nCellsSolve_blk)
!         self % nCellsLocal = self % nCellsLocal + nCellsSolve_blk
!
!         !Useful if local+halo cells are needed at a later time
!         call mpas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells_blk)
!         nCellsHalo = nCellsHalo + nCells_blk
!
!         block_ptr => block_ptr % next
!      end do
!
!      allocate(CellsGlobalToHalo(nCellsHalo))
!      allocate(self % CellsGlobalToLocal(self % nCellsLocal))
!
!      block_ptr => self % domain % blocklist
!      CellsEnd = 0
!      CellsSolveEnd = 0
!      do while(associated(block_ptr))
!         call mpas_pool_get_array(block_ptr % allFields, 'indexToCellID', indexToCellID)
!
!         call mpas_pool_get_dimension(block_ptr % dimensions, 'nCellsSolve', nCellsSolve_blk)
!         CellsSolveStart = CellsSolveEnd + 1
!         CellsSolveEnd = CellsSolveStart + nCellsSolve_blk - 1
!         self % CellsGlobalToLocal(CellsSolveStart:CellsSolveEnd) = indexToCellID(1:nCellsSolve_blk)
!
!         call mpas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells_blk)
!         CellsStart = CellsEnd + 1
!         CellsEnd = CellsStart + nCells_blk - 1
!         CellsGlobalToHalo(CellsStart:CellsEnd) = indexToCellID(1:nCells_blk)
!
!         block_ptr => block_ptr % next
!      end do
!
!      ! Generate index converter from Local+Halo array to Local array
!      allocate(self % CellsMemToLocal(self % nCellsLocal))
!      if (nblock .eq. 1) then
!         self % CellsMemToLocal = self % CellsGlobalToLocal         
!      else
!         ! More complicated when there is more than one block (Is that ever the case?)
!         do ii = 1, nCellsHalo
!            do jj = 1, self % nCellsLocal
!               if ( CellsGlobalToHalo(ii) .eq. self % CellsGlobalToLocal(jj) ) then
!                  self % CellsMemToLocal(jj) = ii
!                  exit
!               end if
!            end do
!         end do
!      end if
!      deallocate(CellsGlobalToHalo)
!
!!write(*,*) 'self % CellsGlobalToLocal = ', self % CellsGlobalToLocal
!!write(*,*) 'self % CellsMemToLocal = ', self % CellsMemToLocal
!
!   end if
!end subroutine geo_to_local   

! ------------------------------------------------------------------------------

subroutine geo_clone(self, other)
implicit none
type(mpas_geom), intent(in) :: self
type(mpas_geom), intent(inout) :: other

write(*,*)'====> copy of geom array'
if (allocated(other%latCell)) then 
   write(*,*)'Allocated array other%latCell'
else
   write(*,*)'Not Allocated array other%latCell'
end if   

other%nCellsGlobal  = self%nCellsGlobal
other%nCells        = self%nCells
other%nCellsSolve   = self%nCellsSolve

other%nEdgesGlobal  = self%nEdgesGlobal
other%nEdges        = self%nEdges
other%nEdgesSolve   = self%nEdgesSolve

other%nVerticesGlobal  = self%nVerticesGlobal
other%nVertices        = self%nVertices
other%nVerticesSolve   = self%nVerticesSolve

other%nVertLevels   = self%nVertLevels
other%nVertLevelsP1 = self%nVertLevelsP1
other%nSoilLevels   = self%nSoilLevels 
other%vertexDegree  = self%vertexDegree
other%maxEdges      = self%maxEdges

if (.not.allocated(other%latCell)) allocate(other%latCell(self%nCells))
if (.not.allocated(other%lonCell)) allocate(other%lonCell(self%nCells))
if (.not.allocated(other%latEdge)) allocate(other%latEdge(self%nEdges))
if (.not.allocated(other%lonEdge)) allocate(other%lonEdge(self%nEdges))
if (.not.allocated(other%xland)) allocate(other%xland(self%nCells))
if (.not.allocated(other%areaCell)) allocate(other%areaCell(self%nCells))
if (.not.allocated(other%edgeNormalVectors)) allocate(other%edgeNormalVectors(3, self%nEdges))
if (.not.allocated(other%zgrid)) allocate(other%zgrid(self%nVertLevelsP1, self%nCells))
!if (.not.allocated(other%edgesOnCell)) allocate(other%edgesOnCell(self%maxEdges,self%nCells))
!if (.not.allocated(other%nEdgesOnCell)) allocate(other%nEdgesOnCell(self%nCells))

other%latCell       = self%latCell
other%lonCell       = self%lonCell
other%areaCell      = self%areaCell
other%latEdge       = self%latEdge
other%lonEdge       = self%lonEdge
other%xland         = self%xland
other%edgeNormalVectors = self%edgeNormalVectors
other%zgrid         = self%zgrid
!other%edgesOnCell   = self%edgesOnCell
!other%nEdgesOnCell  = self%nEdgesOnCell

write(*,*)'====> copy of geom corelist and domain'

if ((associated(other % corelist)).and.(associated(other % domain))) then 
   write(*,*)'associated(other % corelist), associated(other % domain)'
else
   write(*,*)'not associated(other % corelist), associated(other % domain)'
   call mpas_init(other % corelist, other % domain, mpi_comm=MPI_COMM_WORLD)
end if

write(*,*)'====> copy of geom done'

end subroutine geo_clone

! ------------------------------------------------------------------------------

subroutine geo_delete(self)

implicit none
type(mpas_geom), intent(inout) :: self

write(*,*)'==> delete geom array'
if (allocated(self%latCell)) deallocate(self%latCell)
if (allocated(self%lonCell)) deallocate(self%lonCell)
if (allocated(self%latEdge)) deallocate(self%latEdge)
if (allocated(self%lonEdge)) deallocate(self%lonEdge)
if (allocated(self%xland)) deallocate(self%xland)
if (allocated(self%areaCell)) deallocate(self%areaCell)
if (allocated(self%edgeNormalVectors)) deallocate(self%edgeNormalVectors)
if (allocated(self%zgrid)) deallocate(self%zgrid)
!if (allocated(self%nEdgesOnCell)) deallocate(self%nEdgesOnCell)
!if (allocated(self%edgesOnCell)) deallocate(self%edgesOnCell)

if ((associated(self % corelist)).and.(associated(self % domain))) then
   write(*,*)'==> delete geom corelist and domain'
   call mpas_timer_set_context( self % domain )
   call mpas_finalize(self % corelist, self % domain)
end if
write(*,*)'==> delete geom done'
!self % corelist => null()
!self % domain => null()

end subroutine geo_delete

! ------------------------------------------------------------------------------

subroutine geo_info(self, nCellsGlobal, nCells, nCellsSolve, &
                          nEdgesGlobal, nEdges, nEdgesSolve, &
                          nVertLevels, nVertLevelsP1)

implicit none
type(mpas_geom), intent(in) :: self
integer, intent(inout) :: nCellsGlobal, nCells, nCellsSolve
integer, intent(inout) :: nEdgesGlobal, nEdges, nEdgesSolve
integer, intent(inout) :: nVertLevels
integer, intent(inout) :: nVertLevelsP1

nCellsGlobal  = self%nCellsGlobal
nCells        = self%nCells
nCellsSolve   = self%nCellsSolve

nEdgesGlobal  = self%nEdgesGlobal
nEdges        = self%nEdges
nEdgesSolve   = self%nEdgesSolve

nVertLevels   = self%nVertLevels
nVertLevelsP1 = self%nVertLevelsP1

end subroutine geo_info

! ------------------------------------------------------------------------------

end module mpas_geom_mod
