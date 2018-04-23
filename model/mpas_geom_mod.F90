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

write(*,*)'End of geo_setup'

end subroutine geo_setup

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

subroutine geo_info(self, nCells, nEdges, nVertLevels, nVertLevelsP1)

implicit none
type(mpas_geom), intent(in) :: self
integer, intent(inout) :: nCells
integer, intent(inout) :: nEdges
integer, intent(inout) :: nVertLevels
integer, intent(inout) :: nVertLevelsP1

nCells        = self%nCells
nEdges        = self%nEdges
nVertLevels   = self%nVertLevels
nVertLevelsP1 = self%nVertLevelsP1

end subroutine geo_info

! ------------------------------------------------------------------------------

end module mpas_geom_mod
