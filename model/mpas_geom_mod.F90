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
   integer :: nEdges !Global count
   integer :: nVertices !Global count
   integer :: nVertLevels
   integer :: nVertLevelsP1
   integer :: nSoilLevels
   integer :: vertexDegree
   integer :: maxEdges
   integer :: nCellsLocal !Local count
   integer, allocatable :: CellsGlobalToLocal(:) !Local cell global indices (for global MPAS arrays)
   integer, allocatable :: CellsMemToLocal(:) !Local + Halo cell global indices (for local MPAS arrays)
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
call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=self%nCellsGlobal))
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
!self%nCellsGlobal       = config_get_int(c_conf, "nCells")

!> Allocate memory for mesh variables
allocate(self%latCell(self%nCellsGlobal))
allocate(self%lonCell(self%nCellsGlobal))
allocate(self%latEdge(self%nEdges))
allocate(self%lonEdge(self%nEdges))
allocate(self%xland(self%nCellsGlobal))
allocate(self%areaCell(self%nCellsGlobal))
allocate(self%edgeNormalVectors(3, self%nEdges))
allocate(self%zgrid(self%nVertLevelsP1, self%nCellsGlobal))
!allocate(self%nEdgesOnCell(self%nCellsGlobal))
!allocate(self%edgesOnCell(self% maxEdges,self % nCellsGlobal))

!> Read mesh variables
call ncerr(string1, nf90_inq_varid(ncid,'latCell',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%latCell)) !,1,self%nCellsGlobal))
call ncerr(string1, nf90_inq_varid(ncid,'lonCell',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%lonCell)) !,1,self%nCellsGlobal))
call ncerr(string1, nf90_inq_varid(ncid,'areaCell',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%areaCell)) !,1,self%nCellsGlobal))
write(0,*)'BJJ areaCell MIN/MAX: ',minval(self%areaCell),maxval(self%areaCell)
call ncerr(string1, nf90_inq_varid(ncid,'latEdge',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%latEdge)) !,1,self%nEdges))
call ncerr(string1, nf90_inq_varid(ncid,'lonEdge',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%lonEdge)) !,1,self%nEdges))
call ncerr(string1, nf90_inq_varid(ncid,'xland',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%xland)) !,1,self%nCellsGlobal))
call ncerr(string1, nf90_inq_varid(ncid,'edgeNormalVectors',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%edgeNormalVectors,(/1,1/),(/3,self%nEdges/)))
call ncerr(string1, nf90_inq_varid(ncid,'zgrid',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%zgrid,(/1,1/),(/self%nVertLevelsP1,self%nCellsGlobal/)))
!call ncerr(string1, nf90_inq_varid(ncid,'nEdgesOnCell',varid))
!call ncerr(string1, nf90_get_var  (ncid,varid,self%nEdgesOnCell))
!call ncerr(string1, nf90_inq_varid(ncid,'edgesOnCell',varid))
!call ncerr(string1, nf90_get_var  (ncid,varid,self%edgesOnCell,(/1,1/),(/self%maxEdges,self%nCellsGlobal/)))

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

! The mesh info reading above could be replaced with accesses to the pool
!       self % domain % blocklist % structs (after mpas_init)
!  These pool accesses would refer to memory (local+halo) instead of the global sizes and arrays above,
!   requiring updates in mpas_fields_mod, mpas_model_mod, etc...  
!   That effort would streamline code for future parallelization.

!  The code below assumes only a single block is located on each processor.
!   (could remove "Mem" labels if these become the standard vars)

!    call mpas_pool_get_subpool(self % domain % blocklist % structs, 'mesh', meshPool)

!    call mpas_pool_get_dimension(meshPool, 'nCells', dim0d ); self % nCellsMem = dim0d
!    call mpas_pool_get_dimension(meshPool, 'nEdges', dim0d ); self % nEdgesMem = dim0d
!    call mpas_pool_get_dimension(meshPool, 'nVertices', dim0d ); self % nVerticesMem = dim0d

!    call mpas_pool_get_dimension(meshPool, 'nVertLevels', dim0d ); self % nVertLevels = dim0d
!    call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', dim0d ); self % nVertLevelsP1 = dim0d
!    call mpas_pool_get_dimension(meshPool, 'nSoilLevels', dim0d ); self % nSoilLevels = dim0d
!    call mpas_pool_get_dimension(meshPool, 'vertexDegree', dim0d ); self % vertexDegree = dim0d
!    call mpas_pool_get_dimension(meshPool, 'maxEdges', dim0d ); self % maxEdges = dim0d

!  Could make this more flexible/clean by using pointers for array variables in mpas_geom 
!  + value modifications must be consistent with MPAS modeling
!  + should nullify instead of allocate/dellocate

!    allocate(self%latCell(self%nCellsMem))
!    allocate(self%lonCell(self%nCellsMem))
!    allocate(self%latEdge(self%nEdgesMem))
!    allocate(self%lonEdge(self%nEdgesMem))
!    allocate(self%xland(self%nCellsMem))
!    allocate(self%areaCell(self%nCellsMem))
!    allocate(self%edgeNormalVectors(3, self%nEdgesMem))
!    allocate(self%zgrid(self%nVertLevelsP1, self%nCellsMem))

!    call mpas_pool_get_array(meshPool, 'latCell', dummy_1d ); self % latCellMem = dummy_1d
!    call mpas_pool_get_array(meshPool, 'lonCell', dummy_1d ); self % lonCellMem = dummy_1d
!    call mpas_pool_get_array(meshPool, 'areaCell', dummy_1d ); self % areaCellMem = dummy_1d
!    call mpas_pool_get_array(meshPool, 'latEdge', dummy_1d ); self % latEdgeMem = dummy_1d
!    call mpas_pool_get_array(meshPool, 'lonEdge', dummy_1d ); self % lonEdgeMem = dummy_1d
!    call mpas_pool_get_array(meshPool, 'edgeNormalVectors', dummy_2d ); self % edgeNormalVectorsMem = dummy_2d(1:3,1:self%nEdgesMem)
!    call mpas_pool_get_array(meshPool, 'zgrid', dummy_2d ); self % zgridMem = dummy_2d(1:self%nVertLevelsP1,1:self%nCellsMem)

!!    call mpas_pool_get_array(mesh, 'nEdgesOnCell', dim0d ); nEdgesOnCellMem = dim0d
!!    call mpas_pool_get_array(mesh, 'edgesOnCell', dummy_1d ); edgesOnCellMem = dummy_1d(1:self%maxEdges,1:self%nCellsMem)

!    call mpas_pool_get_subpool(self % domain % blocklist % structs, 'fg', fg)
!    call mpas_pool_get_array(fg, 'xland', dummy_1d); xlandMem = dummy_1d

call geo_to_local ( self )
 
write(*,*)'End of geo_setup'

end subroutine geo_setup

! ------------------------------------------------------------------------------

subroutine geo_to_local ( self )
   ! Description : This subroutine creates index converters from
   !               global and memory arrays to local arrays used by BUMP
   !               It is required to loop through the local MPAS blocklist, each member
   !               of which contains a subset of the local geometry (if nblock>1).
   ! See the subroutine mpas_block_creator_finalize_block_phase1 in mpas_block_creator.F 
   ! for more details.

   implicit none

   type (mpas_geom), intent(inout) :: self
   type (block_type), pointer :: block_ptr
   integer, pointer :: indexToCellID(:)
   integer, allocatable :: CellsGlobalToHalo(:)
   integer, pointer :: nCells_blk, nCellsSolve_blk
   integer          :: nCellsHalo, CellsStart, CellsEnd, CellsSolveStart, CellsSolveEnd, ii, jj, nblock

   if (.not. allocated(self % CellsGlobalToLocal)) then
      block_ptr => self % domain % blocklist
      nCellsHalo = 0
      self % nCellsLocal = 0
      nblock = 0
      do while(associated(block_ptr))
         nblock = nblock + 1
         call mpas_pool_get_dimension(block_ptr % dimensions, 'nCellsSolve', nCellsSolve_blk)
         self % nCellsLocal = self % nCellsLocal + nCellsSolve_blk

         !Useful if local+halo cells are needed at a later time
         call mpas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells_blk)
         nCellsHalo = nCellsHalo + nCells_blk

         block_ptr => block_ptr % next
      end do

      allocate(CellsGlobalToHalo(nCellsHalo))
      allocate(self % CellsGlobalToLocal(self % nCellsLocal))

      block_ptr => self % domain % blocklist
      CellsEnd = 0
      CellsSolveEnd = 0
      do while(associated(block_ptr))
         call mpas_pool_get_array(block_ptr % allFields, 'indexToCellID', indexToCellID)

         call mpas_pool_get_dimension(block_ptr % dimensions, 'nCellsSolve', nCellsSolve_blk)
         CellsSolveStart = CellsSolveEnd + 1
         CellsSolveEnd = CellsSolveStart + nCellsSolve_blk - 1
         self % CellsGlobalToLocal(CellsSolveStart:CellsSolveEnd) = indexToCellID(1:nCellsSolve_blk)

         call mpas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells_blk)
         CellsStart = CellsEnd + 1
         CellsEnd = CellsStart + nCells_blk - 1
         CellsGlobalToHalo(CellsStart:CellsEnd) = indexToCellID(1:nCells_blk)

         block_ptr => block_ptr % next
      end do

      ! Generate index converter from Local+Halo array to Local array
      allocate(self % CellsMemToLocal(self % nCellsLocal))
      if (nblock .eq. 1) then
         self % CellsMemToLocal = self % CellsGlobalToLocal         
      else
         ! More complicated when there is more than one block (Is that ever the case?)
         do ii = 1, nCellsHalo
            do jj = 1, self % nCellsLocal
               if ( CellsGlobalToHalo(ii) .eq. self % CellsGlobalToLocal(jj) ) then
                  self % CellsMemToLocal(jj) = ii
                  exit
               end if
            end do
         end do
      end if
      deallocate(CellsGlobalToHalo)

!write(*,*) 'self % CellsGlobalToLocal = ', self % CellsGlobalToLocal
!write(*,*) 'self % CellsMemToLocal = ', self % CellsMemToLocal

   end if
end subroutine geo_to_local   

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

other%nCellsGlobal        = self%nCellsGlobal
other%nEdges        = self%nEdges
other%nVertices     = self%nVertices
other%nVertLevels   = self%nVertLevels
other%nVertLevelsP1 = self%nVertLevelsP1
other%nSoilLevels   = self%nSoilLevels 
other%vertexDegree  = self%vertexDegree
other%maxEdges      = self%maxEdges
other%nCellsLocal   = self%nCellsLocal
if (.not.allocated(other%CellsGlobalToLocal)) allocate(other%CellsGlobalToLocal(self%nCellsLocal))
if (.not.allocated(other%CellsMemToLocal)) allocate(other%CellsMemToLocal(self%nCellsLocal))
other%CellsGlobalToLocal    = self%CellsGlobalToLocal
other%CellsMemToLocal    = self%CellsMemToLocal

if (.not.allocated(other%latCell)) allocate(other%latCell(self%nCellsGlobal))
if (.not.allocated(other%lonCell)) allocate(other%lonCell(self%nCellsGlobal))
if (.not.allocated(other%latEdge)) allocate(other%latEdge(self%nEdges))
if (.not.allocated(other%lonEdge)) allocate(other%lonEdge(self%nEdges))
if (.not.allocated(other%xland)) allocate(other%xland(self%nCellsGlobal))
if (.not.allocated(other%areaCell)) allocate(other%areaCell(self%nCellsGlobal))
if (.not.allocated(other%edgeNormalVectors)) allocate(other%edgeNormalVectors(3, self%nEdges))
if (.not.allocated(other%zgrid)) allocate(other%zgrid(self%nVertLevelsP1, self%nCellsGlobal))
!if (.not.allocated(other%edgesOnCell)) allocate(other%edgesOnCell(self%maxEdges,self%nCellsGlobal))
!if (.not.allocated(other%nEdgesOnCell)) allocate(other%nEdgesOnCell(self%nCellsGlobal))

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
if (allocated(self%CellsGlobalToLocal)) deallocate(self%CellsGlobalToLocal)
if (allocated(self%CellsMemToLocal)) deallocate(self%CellsMemToLocal)
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

subroutine geo_info(self, nCellsGlobal, nEdges, nVertLevels, nVertLevelsP1, nCellsLocal)

implicit none
type(mpas_geom), intent(in) :: self
integer, intent(inout) :: nCellsGlobal
integer, intent(inout) :: nEdges
integer, intent(inout) :: nVertLevels
integer, intent(inout) :: nVertLevelsP1
integer, intent(inout) :: nCellsLocal

nCellsGlobal        = self%nCellsGlobal
nEdges        = self%nEdges
nVertLevels   = self%nVertLevels
nVertLevelsP1 = self%nVertLevelsP1
nCellsLocal   = self%nCellsLocal

end subroutine geo_info

! ------------------------------------------------------------------------------

end module mpas_geom_mod
