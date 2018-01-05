! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_geom_mod

use iso_c_binding
use config_mod
use netcdf
use tools_nc

use mpas_derived_types
use mpas_kind_types
use mpas_constants

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
   real(kind=RKIND), DIMENSION(:),   ALLOCATABLE :: latCell, lonCell, xland
   real(kind=RKIND), DIMENSION(:),   ALLOCATABLE :: areaCell
   real(kind=RKIND), DIMENSION(:),   ALLOCATABLE :: latEdge, lonEdge
   real(kind=RKIND), DIMENSION(:,:), ALLOCATABLE :: edgeNormalVectors
   real(kind=RKIND), DIMENSION(:,:), ALLOCATABLE :: zgrid
   type (dm_info), pointer :: dminfo
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
type(mpas_geom)        :: self
type(c_ptr), intent(in) :: c_conf
character(len=StrKIND) :: string1
integer :: ncid, dimid, varid
real(kind=RKIND), parameter :: deg2rad = pii/180.      

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

!> Read mesh variables
call ncerr(string1, nf90_inq_varid(ncid,'latCell',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%latCell)) !,1,self%nCells))
call ncerr(string1, nf90_inq_varid(ncid,'lonCell',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%lonCell)) !,1,self%nCells))
call ncerr(string1, nf90_inq_varid(ncid,'areaCell',varid))
call ncerr(string1, nf90_get_var  (ncid,varid,self%areaCell)) !,1,self%nCells))
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

!> radians to degrees
self%latCell = self%latCell / deg2rad
self%lonCell = self%lonCell / deg2rad
self%latEdge = self%latEdge / deg2rad
self%lonEdge = self%lonEdge / deg2rad

!> close file
call ncerr(string1,nf90_close(ncid))

end subroutine geo_setup

! ------------------------------------------------------------------------------

subroutine geo_clone(self, other)
implicit none
type(mpas_geom), pointer :: self, other

other%nCells        = self%nCells
other%nEdges        = self%nEdges
other%nVertices     = self%nVertices
other%nVertLevels   = self%nVertLevels
other%nVertLevelsP1 = self%nVertLevelsP1
other%nSoilLevels   = self%nSoilLevels 
other%vertexDegree  = self%vertexDegree
other%maxEdges      = self%maxEdges
other%latCell       = self%latCell
other%lonCell       = self%lonCell
other%areaCell      = self%areaCell
other%latEdge       = self%latEdge
other%lonEdge       = self%lonEdge
other%xland         = self%xland
other%edgeNormalVectors = self%edgeNormalVectors
other%zgrid         = self%zgrid

end subroutine geo_clone

! ------------------------------------------------------------------------------

subroutine geo_delete(self)

implicit none
type(mpas_geom), pointer :: self

if (associated(self % dminfo)) nullify(self % dminfo) 

end subroutine geo_delete

! ------------------------------------------------------------------------------

subroutine geo_info(self, nCells, nEdges, nVertLevels, nVertLevelsP1)

implicit none
type(mpas_geom), pointer :: self
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
