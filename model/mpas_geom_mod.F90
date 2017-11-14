
!> Fortran module handling geometry for the MPAS-A model

module mpas_geom_mod

use iso_c_binding
use config_mod
use kinds
use tools_nc
use netcdf
use mpas_constants, only: deg2rad

implicit none
private
public :: mpas_geom
public :: mpas_geom_registry

! ------------------------------------------------------------------------------
 integer, parameter :: max_string_length = 128
! ------------------------------------------------------------------------------

!> Fortran derived type to hold geometry data for the MPAS-A model

!> Grid dimensions
type :: mpas_geom
 integer :: nCells
 integer :: nEdges
 integer :: nVertices
 integer :: nVertLevels
 integer :: nVertLevelsP1
 integer :: nSoilLevels
 integer :: vertexDegree
 integer :: maxEdges
 character(len=max_string_length) :: gridfname
 real(kind=kind_real), DIMENSION(:),   ALLOCATABLE :: latCell, lonCell, xland
 real(kind=kind_real), DIMENSION(:),   ALLOCATABLE :: areaCell
 real(kind=kind_real), DIMENSION(:),   ALLOCATABLE :: latEdge, lonEdge
 real(kind=kind_real), DIMENSION(:,:), ALLOCATABLE :: edgeNormalVectors
 real(kind=kind_real), DIMENSION(:,:), ALLOCATABLE :: zgrid
end type mpas_geom

#define LISTED_TYPE mpas_geom

!> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_geom_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_mpas_geo_setup(c_key_self, c_conf) bind(c,name='mpas_geo_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)    :: c_conf

type(mpas_geom), pointer :: self

character(len=max_string_length) :: string1
integer :: ncid, dimid, varid

call mpas_geom_registry%init()
call mpas_geom_registry%add(c_key_self)
call mpas_geom_registry%get(c_key_self,self)

!> Open a grid mesh file
self%gridfname = config_get_string(c_conf, max_string_length, "gridfname")
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

end subroutine c_mpas_geo_setup

! ------------------------------------------------------------------------------

subroutine c_mpas_geo_clone(c_key_self, c_key_other) bind(c,name='mpas_geo_clone_f90')
implicit none
integer(c_int), intent(in   ) :: c_key_self
integer(c_int), intent(inout) :: c_key_other

type(mpas_geom), pointer :: self, other

call mpas_geom_registry%add(c_key_other)
call mpas_geom_registry%get(c_key_other, other)
call mpas_geom_registry%get(c_key_self , self )

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

end subroutine c_mpas_geo_clone

! ------------------------------------------------------------------------------

subroutine c_mpas_geo_delete(c_key_self) bind(c,name='mpas_geo_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self     

call mpas_geom_registry%remove(c_key_self)

end subroutine c_mpas_geo_delete

! ------------------------------------------------------------------------------

subroutine c_mpas_geo_info(c_key_self, c_nCells, c_nEdges, c_nVertLevels, c_nVertLevelsP1 ) bind(c,name='mpas_geo_info_f90')
implicit none
integer(c_int), intent(in   ) :: c_key_self
integer(c_int), intent(inout) :: c_nCells
integer(c_int), intent(inout) :: c_nEdges
integer(c_int), intent(inout) :: c_nVertLevels
integer(c_int), intent(inout) :: c_nVertLevelsP1
type(mpas_geom), pointer :: self

call mpas_geom_registry%get(c_key_self , self )
c_nCells        = self%nCells
c_nEdges        = self%nEdges
c_nVertLevels   = self%nVertLevels
c_nVertLevelsP1 = self%nVertLevelsP1

end subroutine c_mpas_geo_info

! ------------------------------------------------------------------------------

end module mpas_geom_mod
