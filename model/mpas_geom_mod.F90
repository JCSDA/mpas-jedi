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
   integer :: nCells
   integer :: nEdges
   integer :: nVertices
   integer :: nVertLevels
   integer :: nVertLevelsP1
   integer :: nSoilLevels
   integer :: vertexDegree
   integer :: maxEdges
   integer :: nCellsLocal
   integer, allocatable :: CellsLocal(:)
!   integer :: nEdgesLocal
!   integer, allocatable :: EdgesLocal(:)
!   integer :: nVerticesLocal
!   integer, allocatable :: VerticesLocal(:)

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

write(*,*)' ==> create geom'

!> Open a grid mesh file
self%gridfname = config_get_string(c_conf, StrKIND, "gridfname")
string1 = self%gridfname
call ncerr(string1, nf90_open(trim(self%gridfname),nf90_nowrite,ncid))

!> Grid dimensions
call ncerr(string1, nf90_inq_dimid        (ncid,'nCells',dimid))
call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=self%nCells))
call ncerr(string1, nf90_inq_dimid        (ncid,'nEdges',dimid))
call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=self%nEdges))
call ncerr(string1, nf90_inq_dimid        (ncid,'nVertices',dimid))
call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=self%nVertices))
call ncerr(string1, nf90_inq_dimid        (ncid,'nVertLevels',dimid))
call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=self%nVertLevels))
call ncerr(string1, nf90_inq_dimid        (ncid,'nVertLevelsP1',dimid))
call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=self%nVertLevelsP1))
call ncerr(string1, nf90_inq_dimid        (ncid,'nSoilLevels',dimid))
call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=self%nSoilLevels))
call ncerr(string1, nf90_inq_dimid        (ncid,'vertexDegree',dimid))
call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=self%vertexDegree))
call ncerr(string1, nf90_inq_dimid        (ncid,'maxEdges',dimid))
call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=self%maxEdges))
!self%nCells       = config_get_int(c_conf, "nCells")

!> Allocate memory for mesh variables
allocate(self%latCell(self%nCells))
allocate(self%lonCell(self%nCells))
allocate(self%latEdge(self%nEdges))
allocate(self%lonEdge(self%nEdges))
allocate(self%xland(self%nCells))
allocate(self%areaCell(self%nCells))
allocate(self%edgeNormalVectors(3, self%nEdges))
allocate(self%zgrid(self%nVertLevelsP1, self%nCells))
!allocate(self%nEdgesOnCell(self%nCells))
!allocate(self%edgesOnCell(self% maxEdges,self % nCells))

!> Read mesh variables
call ncerr(string1, nf90_inq_varid(ncid,'latCell',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%latCell)) !,1,self%nCells))
call ncerr(string1, nf90_inq_varid(ncid,'lonCell',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%lonCell)) !,1,self%nCells))
call ncerr(string1, nf90_inq_varid(ncid,'areaCell',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%areaCell)) !,1,self%nCells))
write(0,*)'BJJ areaCell MIN/MAX: ',minval(self%areaCell),maxval(self%areaCell)
call ncerr(string1, nf90_inq_varid(ncid,'latEdge',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%latEdge)) !,1,self%nEdges))
call ncerr(string1, nf90_inq_varid(ncid,'lonEdge',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%lonEdge)) !,1,self%nEdges))
call ncerr(string1, nf90_inq_varid(ncid,'xland',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%xland)) !,1,self%nCells))
call ncerr(string1, nf90_inq_varid(ncid,'edgeNormalVectors',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%edgeNormalVectors,(/1,1/),(/3,self%nEdges/)))
call ncerr(string1, nf90_inq_varid(ncid,'zgrid',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%zgrid,(/1,1/),(/self%nVertLevelsP1,self%nCells/)))
!call ncerr(string1, nf90_inq_varid(ncid,'nEdgesOnCell',varid))
!call ncerr(string1, nf90_get_var  (ncid,varid,self%nEdgesOnCell))
!call ncerr(string1, nf90_inq_varid(ncid,'edgesOnCell',varid))
!call ncerr(string1, nf90_get_var  (ncid,varid,self%edgesOnCell,(/1,1/),(/self%maxEdges,self%nCells/)))

!> radians to degrees
self%latCell = self%latCell !/ deg2rad
self%lonCell = self%lonCell !/ deg2rad
self%latEdge = self%latEdge !/ deg2rad
self%lonEdge = self%lonEdge !/ deg2rad

!> close file
call ncerr(string1,nf90_close(ncid))

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

call geo_get_local ( self )
 
write(*,*)'End of geo_setup'

end subroutine geo_setup

! ------------------------------------------------------------------------------

subroutine geo_get_local ( self )
   ! Description : This subroutine transfers the local geom info to an mpas_geom type.
   !               It is required to loop through the local MPAS blocklist, each member
   !               of which contains a subset of the local geometry.
   ! See the subroutine mpas_block_creator_finalize_block_phase1 in mpas_block_creator.F 
   ! for more details.

   implicit none

   type (mpas_geom), intent(inout) :: self
   type (block_type), pointer :: block_ptr
   integer, pointer :: nCells_blk, indexToCellIDPool(:)!, nCellsSolveArray_blk(:)
   integer          :: CellStart, CellEnd
!   integer, pointer :: nEdges_blk, indexToEdgeIDPool(:)
!   integer          :: EdgeStart, EdgeEnd
!   integer, pointer :: nVertices_blk, indexToVerticeIDPool(:)
!   integer          :: VerticeStart, VerticeEnd
   integer, pointer :: dummy_1d(:), dummy_2d(:,:)

   if (.not. allocated(self % CellsLocal)) then
      block_ptr => self % domain % blocklist
      self % nCellsLocal = 0
      do while(associated(block_ptr))
         call mpas_pool_get_dimension(block_ptr % dimensions, 'nCellsSolve', nCells_blk)
write(*,*) 'nCellsSolve_blk = ', nCells_blk
         call mpas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells_blk)
write(*,*) 'nCells = ', nCells_blk

         self % nCellsLocal = self % nCellsLocal + nCells_blk


         call mpas_pool_get_dimension(block_ptr % dimensions, 'nCellsArray', dummy_1d)
write(*,*) 'size(nCellsArray) = ', size(dummy_1d)
write(*,*) 'nCellsArray = ', dummy_1d
         call mpas_pool_get_array(block_ptr % dimensions, 'indexToCellID_blk', dummy_1d)
write(*,*) 'size(indexToCellID_blk) = ', size(dummy_1d)
write(*,*) 'min(indexToCellID_blk) = ', min(dummy_1d)
write(*,*) 'max(indexToCellID_blk) = ', max(dummy_1d)
write(*,*) 'count(indexToCellID_blk.gt.0) = ', count(dummy_1d.gt.0)
write(*,*) 'indexToCellID_blk = ', dummy_1d
         call mpas_pool_get_array(block_ptr % dimensions, 'indexToCellID', dummy_1d)
write(*,*) 'size(indexToCellID) = ', size(dummy_1d)
write(*,*) 'min(indexToCellID) = ', min(dummy_1d)
write(*,*) 'max(indexToCellID) = ', max(dummy_1d)
write(*,*) 'count(indexToCellID.gt.0) = ', count(dummy_1d.gt.0)
write(*,*) 'indexToCellID = ', dummy_1d
         call mpas_pool_get_array(block_ptr % dimensions, 'cellsOnCell', dummy_2d)
write(*,*) 'size(cellsOnCell,1) = ', size(dummy_2d,1)
write(*,*) 'size(cellsOnCell,2) = ', size(dummy_2d,2)
write(*,*) 'count(cellsOnCell.gt.0) = ', count(dummy_2d.gt.0)
!write(*,*) 'sum(cellsOnCell) = ', sum(dummy_2d)
!write(*,*) 'cellsOnCell = ', dummy_2d

!         call mpas_pool_get_dimension(block_ptr % dimensions, 'nEdgesSolve', nEdgesSolve_blk)
!         self % nEdgesLocal = self % nEdgesLocal + nEdgesSolve_blk

!         call mpas_pool_get_dimension(block_ptr % dimensions, 'nVerticesSolve', nVerticesSolve_blk)
!         self % nVerticesLocal = self % nVerticesLocal + nVerticesSolve_blk

         block_ptr => block_ptr % next
      end do

      allocate(self % CellsLocal(self % nCellsLocal))

      block_ptr => self % domain % blocklist
      CellEnd = 0
      do while(associated(block_ptr))
         call mpas_pool_get_dimension(block_ptr % dimensions, 'nCellsSolve', nCellsSolve_blk)
         call mpas_pool_get_array(block_ptr % allFields, 'indexToCellID_blk', indexToCellIDPool)
         CellStart = CellEnd + 1
         CellEnd = CellStart + nCellsSolve_blk - 1
         self % CellsLocal(CellStart:CellEnd) = indexToCellIDPool(1:nCellsSolve_blk)

!         call mpas_pool_get_dimension(block_ptr % dimensions, 'nEdgesSolve', nEdgesSolve_blk)
!         call mpas_pool_get_array(block_ptr % allFields, 'indexToEdgeID_blk', indexToEdgeIDPool)
!         EdgeStart = EdgeEnd + 1
!         EdgeEnd = EdgeStart + nEdgesSolve_blk - 1
!         self % EdgesLocal(EdgeStart:EdgeEnd) = indexToEdgeIDPool(1:nEdgesSolve_blk)

!         call mpas_pool_get_dimension(block_ptr % dimensions, 'nVerticesSolve', nVerticesSolve_blk)
!         call mpas_pool_get_array(block_ptr % allFields, 'indexToVerticeID_blk', indexToVerticeIDPool)
!         VerticeStart = VerticeEnd + 1
!         VerticeEnd = VerticeStart + nVerticesSolve_blk - 1
!         self % VerticesLocal(VerticeStart:VerticeEnd) = indexToVerticeIDPool(1:nVerticesSolve_blk)

         block_ptr => block_ptr % next
      end do
   end if
end subroutine geo_get_local   

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

other%nCells        = self%nCells
other%nEdges        = self%nEdges
other%nVertices     = self%nVertices
other%nVertLevels   = self%nVertLevels
other%nVertLevelsP1 = self%nVertLevelsP1
other%nSoilLevels   = self%nSoilLevels 
other%vertexDegree  = self%vertexDegree
other%maxEdges      = self%maxEdges
other%nCellsLocal   = self%nCellsLocal
if (.not.allocated(other%CellsLocal)) allocate(other%CellsLocal(self%nCellsLocal))
other%CellsLocal    = self%CellsLocal

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
if (allocated(self%CellsLocal)) deallocate(self%CellsLocal)
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

subroutine geo_info(self, nCells, nEdges, nVertLevels, nVertLevelsP1, nCellsLocal)

implicit none
type(mpas_geom), intent(in) :: self
integer, intent(inout) :: nCells
integer, intent(inout) :: nEdges
integer, intent(inout) :: nVertLevels
integer, intent(inout) :: nVertLevelsP1
integer, intent(inout) :: nCellsLocal

nCells        = self%nCells
nEdges        = self%nEdges
nVertLevels   = self%nVertLevels
nVertLevelsP1 = self%nVertLevelsP1
nCellsLocal   = self%nCellsLocal

end subroutine geo_info

! ------------------------------------------------------------------------------

end module mpas_geom_mod
