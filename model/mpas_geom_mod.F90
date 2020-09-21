! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module mpas_geom_mod

use atlas_module, only: atlas_functionspace, atlas_fieldset, atlas_field, atlas_real
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm
use iso_c_binding

!MPAS-Model
use mpas_derived_types
use mpas_kind_types
use mpas_constants
use kinds, only : kind_real
use mpas_dmpar, only: mpas_dmpar_sum_int
use mpas_subdriver
use atm_core
use mpas_pool_routines

!mpas_jedi
use mpas_constants_mod

implicit none
private
public :: mpas_geom, &
        & geo_setup, geo_set_atlas_lonlat, geo_fill_atlas_fieldset, geo_clone, geo_delete, geo_info
public :: mpas_geom_registry


! ------------------------------------------------------------------------------

!> Fortran derived type to hold geometry definition
type :: mpas_geom
   integer :: nCellsGlobal !Global count
   integer :: nEdgesGlobal !Global count
   integer :: nVerticesGlobal !Global count
   integer :: nCells !Memory count (Local + Halo)
   integer :: nEdges !Memory count (Local + Halo)
   integer :: nVertices !Memory count (Local + Halo)
   integer :: nCellsSolve !Local count
   integer :: nEdgesSolve !Local count
   integer :: nVerticesSolve !Local count
   integer :: nVertLevels
   integer :: nVertLevelsP1
   integer :: nSoilLevels
   integer :: vertexDegree
   integer :: maxEdges
   character(len=StrKIND) :: gridfname
   real(kind=kind_real), dimension(:),   allocatable :: latCell, lonCell
   real(kind=kind_real), dimension(:),   allocatable :: areaCell
   real(kind=kind_real), dimension(:),   allocatable :: latEdge, lonEdge
   real(kind=kind_real), dimension(:,:), allocatable :: edgeNormalVectors
   real(kind=kind_real), dimension(:,:), allocatable :: zgrid
   integer, allocatable :: nEdgesOnCell(:)
   integer, allocatable :: cellsOnCell(:,:)
   integer, allocatable :: edgesOnCell(:,:)
   integer, allocatable :: cellsOnVertex(:,:)
   integer, allocatable :: cellsOnEdge(:,:)
   integer, allocatable :: verticesOnEdge(:,:)
   real(kind=kind_real), DIMENSION(:), ALLOCATABLE :: dcEdge, dvEdge
   real(kind=kind_real), DIMENSION(:), ALLOCATABLE :: areaTriangle, angleEdge
   real(kind=kind_real), DIMENSION(:,:), ALLOCATABLE :: kiteAreasOnVertex, edgesOnCell_sign

   type (domain_type), pointer :: domain => null()
   type (core_type), pointer :: corelist => null()

   type(fckit_mpi_comm) :: f_comm

   type(atlas_functionspace) :: afunctionspace
end type mpas_geom

type :: idcounter
   integer :: id, counter
end type idcounter

type(idcounter), allocatable :: geom_count(:)

#define LISTED_TYPE mpas_geom

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_geom_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------
subroutine geo_setup(self, f_conf, f_comm)

   implicit none

   type(mpas_geom),           intent(inout) :: self
   type(fckit_configuration), intent(in)    :: f_conf
   type(fckit_mpi_comm),   intent(in)    :: f_comm

   character(len=StrKIND) :: string1
   type (mpas_pool_type), pointer :: meshPool, fg
   type (block_type), pointer :: block_ptr

   real (kind=kind_real), pointer :: r1d_ptr(:), r2d_ptr(:,:)
   integer, pointer :: i0d_ptr, i1d_ptr(:), i2d_ptr(:,:)

   type(idcounter), allocatable :: prev_count(:)
   integer :: ii, nprev
!   write(*,*)' ==> create geom'

   ! MPI communicator
   self%f_comm = f_comm

   ! Domain decomposition and templates for state/increment variables
   call mpas_init( self % corelist, self % domain, mpi_comm = self%f_comm%communicator() )

   if (allocated(geom_count)) then
      nprev = size(geom_count)
      allocate(prev_count(nprev))
      do ii = 1, nprev
         prev_count(ii) = geom_count(ii)
         if (prev_count(ii)%id == self%domain%domainID) then
            call abor1_ftn("domainID already used")
         end if
      end do
      deallocate(geom_count)
   else
      nprev = 0
   end if
   allocate(geom_count(nprev+1))
   do ii = 1, nprev
      geom_count(ii) = prev_count(ii)
   end do
   geom_count(nprev+1)%id = self%domain%domainID
   geom_count(nprev+1)%counter = 1

!   if (associated(self % domain)) then
!       write(*,*)'inside geom: geom % domain associated for domainID = ', self % domain % domainID
!   end if
!   if (associated(self % corelist)) then
!       write(*,*)'inside geom: geom % corelist associated'
!   else
!       write(*,*)'inside geom: geom % corelist not associated'
!   end if

   !  These pool accesses refer to memory (local+halo) for a single MPAS block (standard)
   block_ptr => self % domain % blocklist

   call mpas_pool_get_subpool ( block_ptr % structs, 'mesh', meshPool )

   call mpas_pool_get_dimension ( block_ptr % dimensions, 'nCells', i0d_ptr )
   self % nCells = i0d_ptr
   call mpas_pool_get_dimension ( block_ptr % dimensions, 'nCellsSolve', i0d_ptr )
   self % nCellsSolve = i0d_ptr
   call mpas_dmpar_sum_int ( self % domain % dminfo, &
                             self % nCellsSolve, self % nCellsGlobal )

   call mpas_pool_get_dimension ( block_ptr % dimensions, 'nEdges', i0d_ptr )
   self % nEdges = i0d_ptr
   call mpas_pool_get_dimension ( block_ptr % dimensions, 'nEdgesSolve', i0d_ptr )
   self % nEdgesSolve = i0d_ptr
   call mpas_dmpar_sum_int ( self % domain % dminfo, &
                             self % nEdgesSolve, self % nEdgesGlobal )

   call mpas_pool_get_dimension ( block_ptr % dimensions, 'nVertices', i0d_ptr )
   self % nVertices = i0d_ptr
   call mpas_pool_get_dimension ( block_ptr % dimensions, 'nVerticesSolve', i0d_ptr )
   self % nVerticesSolve = i0d_ptr
   call mpas_dmpar_sum_int ( self % domain % dminfo, &
                             self % nVerticesSolve, self % nVerticesGlobal )

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

   allocate ( self % latCell ( self % nCells ) )
   allocate ( self % lonCell ( self % nCells ) )
   allocate ( self % latEdge ( self % nEdges ) )
   allocate ( self % lonEdge ( self % nEdges ) )
   allocate ( self % areaCell ( self % nCells ) )
   allocate ( self % edgeNormalVectors (3,  self % nEdges ) )
   allocate ( self % zgrid ( self % nVertLevelsP1, self % nCells ) )
   allocate ( self % nEdgesOnCell ( self % nCells ) )
   allocate ( self % edgesOnCell ( self % maxEdges, self % nCells ) )
   allocate ( self % cellsOnCell ( self % maxEdges, self % nCells ) )
   allocate ( self % cellsOnVertex ( self % vertexDegree, self % nVertices ) )
   allocate ( self % cellsOnEdge ( 2, self % nEdges ) )
   allocate ( self % verticesOnEdge ( 2, self % nEdges ) )
   allocate ( self % dcEdge ( self % nEdges ) )
   allocate ( self % dvEdge ( self % nEdges ) )
   allocate ( self % kiteAreasOnVertex ( self % vertexDegree, self % nVertices ) )
   allocate ( self % edgesOnCell_sign ( self % maxEdges, self % nCells ) )
   allocate ( self % areaTriangle ( self % nVertices ) )
   allocate ( self % angleEdge ( self % nEdges ) )

   call mpas_pool_get_array ( meshPool, 'latCell', r1d_ptr )
   self % latCell = r1d_ptr(1:self % nCells)
   where (self % latCell > MPAS_JEDI_PIIo2_kr)
       self % latCell = MPAS_JEDI_PIIo2_kr
   end where
   where (self % latCell < - MPAS_JEDI_PIIo2_kr)
       self % latCell = - MPAS_JEDI_PIIo2_kr
   end where

   call mpas_pool_get_array ( meshPool, 'lonCell', r1d_ptr )
   self % lonCell = r1d_ptr(1:self % nCells)
   call mpas_pool_get_array ( meshPool, 'areaCell', r1d_ptr )
   self % areaCell = r1d_ptr(1:self % nCells)
   call mpas_pool_get_array ( meshPool, 'latEdge', r1d_ptr )
   self % latEdge = r1d_ptr(1:self % nEdges)
   call mpas_pool_get_array ( meshPool, 'lonEdge', r1d_ptr )
   self % lonEdge = r1d_ptr(1:self % nEdges)
   call mpas_pool_get_array ( meshPool, 'edgeNormalVectors', r2d_ptr )
   self % edgeNormalVectors = r2d_ptr ( 1:3, 1:self % nEdges )
   call mpas_pool_get_array ( meshPool, 'nEdgesOnCell', i1d_ptr )
   self % nEdgesOnCell = i1d_ptr(1:self % nCells)
   call mpas_pool_get_array ( meshPool, 'edgesOnCell', i2d_ptr )
   self % edgesOnCell = i2d_ptr ( 1:self % maxEdges, 1:self % nCells )
   call mpas_pool_get_array ( meshPool, 'cellsOnCell', i2d_ptr )
   self % cellsOnCell = i2d_ptr( 1:self % maxEdges, 1:self % nCells )
   call mpas_pool_get_array ( meshPool, 'cellsOnVertex', i2d_ptr )
   self % cellsOnVertex = i2d_ptr( 1:self % vertexDegree, 1:self % nVertices )
   call mpas_pool_get_array ( meshPool, 'cellsOnEdge', i2d_ptr )
   self % cellsOnEdge = i2d_ptr( 1:2, 1:self % nEdges )
   call mpas_pool_get_array ( meshPool, 'verticesOnEdge', i2d_ptr )
   self % verticesOnEdge = i2d_ptr( 1:2, 1:self % nEdges )
   call mpas_pool_get_array ( meshPool, 'dcEdge', r1d_ptr )           
   self % dcEdge = r1d_ptr(1:self % nEdges)
   call mpas_pool_get_array ( meshPool, 'dvEdge', r1d_ptr )           
   self % dvEdge = r1d_ptr(1:self % nEdges)
   call mpas_pool_get_array ( meshPool, 'kiteAreasOnVertex', r2d_ptr )             
   self % kiteAreasOnVertex = r2d_ptr ( 1:self % vertexDegree, 1:self % nVertices )
   call mpas_pool_get_array ( meshPool, 'edgesOnCell_sign', r2d_ptr )             
   self % edgesOnCell_sign = r2d_ptr ( 1:self % maxEdges, 1:self % nCells )
   call mpas_pool_get_array ( meshPool, 'areaTriangle', r1d_ptr )           
   self % areaTriangle = r1d_ptr(1:self % nVertices)
   call mpas_pool_get_array ( meshPool, 'angleEdge', r1d_ptr )           
   self % angleEdge = r1d_ptr(1:self % nEdges)

   call mpas_pool_get_array ( meshPool, 'zgrid', r2d_ptr )
   self % zgrid = r2d_ptr ( 1:self % nVertLevelsP1, 1:self % nCells )

!   write(*,*)'End of geo_setup'

end subroutine geo_setup

! --------------------------------------------------------------------------------------------------

subroutine geo_set_atlas_lonlat(self, afieldset)

   implicit none

   type(mpas_geom),  intent(inout) :: self
   type(atlas_fieldset), intent(inout) :: afieldset

   real(kind_real), pointer :: real_ptr(:,:)
   type(atlas_field) :: afield

   ! Create lon/lat field
   afield = atlas_field(name="lonlat", kind=atlas_real(kind_real), shape=(/2,self%nCellsSolve/))
   call afield%data(real_ptr)
   real_ptr(1,:) = self%lonCell(1:self%nCellsSolve) * MPAS_JEDI_RAD2DEG_kr
   real_ptr(2,:) = self%latCell(1:self%nCellsSolve) * MPAS_JEDI_RAD2DEG_kr
   call afieldset%add(afield)

end subroutine geo_set_atlas_lonlat

! --------------------------------------------------------------------------------------------------

subroutine geo_fill_atlas_fieldset(self, afieldset)

   implicit none

   type(mpas_geom),  intent(inout) :: self
   type(atlas_fieldset), intent(inout) :: afieldset

   integer :: i, jz
   real(kind=kind_real), pointer :: real_ptr_1(:), real_ptr_2(:,:)
   type(atlas_field) :: afield

   ! Add area
   afield = self%afunctionspace%create_field(name='area', kind=atlas_real(kind_real), levels=0)
   call afield%data(real_ptr_1)
   real_ptr_1 = self%areaCell(1:self%nCellsSolve)
   call afieldset%add(afield)
   call afield%final()

   ! Add vertical unit
   afield = self%afunctionspace%create_field(name='vunit', kind=atlas_real(kind_real), levels=self%nVertLevels)
   call afield%data(real_ptr_2)
   do jz=1,self%nVertLevels
      real_ptr_2(jz,:) = real(jz, kind_real)
      ! Try different vertical coordinate
      !--For height
      !real_ptr_2(jz,1:self%nCellsSolve) = MPAS_JEDI_HALF_kr * &
      !   ( self%zgrid(jz,1:self%nCellsSolve) + self%zgrid(jz+1,1:self%nCellsSolve) )
      !--Similarly for ln of pressure_base. I don't know if we can access to "total pressure" here.
      !real (kind=kind_real), dimension(:,:), pointer :: pressure_base
      !call mpas_pool_get_array(self%domain%blocklist%allFields,'pressure_base', pressure_base)
      !real_ptr_2(jz,1:self%nCellsSolve) = log( pressure_base(jz,1:self%nCellsSolve) )
   end do
   call afieldset%add(afield)
   call afield%final()

end subroutine geo_fill_atlas_fieldset

! ------------------------------------------------------------------------------

subroutine geo_clone(self, other)

   implicit none

   type(mpas_geom), intent(inout) :: self
   type(mpas_geom), intent(in) :: other

   integer :: ii

   ! Clone communicator
   self%f_comm = other%f_comm

!   write(*,*)'====> copy of geom array'

   self % nCellsGlobal  = other % nCellsGlobal
   self % nCells        = other % nCells
   self % nCellsSolve   = other % nCellsSolve

   self % nEdgesGlobal  = other % nEdgesGlobal
   self % nEdges        = other % nEdges
   self % nEdgesSolve   = other % nEdgesSolve

   self % nVerticesGlobal  = other % nVerticesGlobal
   self % nVertices        = other % nVertices
   self % nVerticesSolve   = other % nVerticesSolve

   self % nVertLevels   = other % nVertLevels
   self % nVertLevelsP1 = other % nVertLevelsP1
   self % nSoilLevels   = other % nSoilLevels
   self % vertexDegree  = other % vertexDegree
   self % maxEdges      = other % maxEdges

   if (.not.allocated(self % latCell)) allocate(self % latCell(other % nCells))
   if (.not.allocated(self % lonCell)) allocate(self % lonCell(other % nCells))
   if (.not.allocated(self % latEdge)) allocate(self % latEdge(other % nEdges))
   if (.not.allocated(self % lonEdge)) allocate(self % lonEdge(other % nEdges))
   if (.not.allocated(self % areaCell)) allocate(self % areaCell(other % nCells))
   if (.not.allocated(self % edgeNormalVectors)) allocate(self % edgeNormalVectors(3, other % nEdges))
   if (.not.allocated(self % zgrid)) allocate(self % zgrid(other % nVertLevelsP1, other % nCells))
   if (.not.allocated(self % nEdgesOnCell)) allocate(self % nEdgesOnCell(other % nCells))
   if (.not.allocated(self % edgesOnCell)) allocate(self % edgesOnCell(other % maxEdges, other % nCells))
   if (.not.allocated(self % cellsOnCell)) allocate(self % cellsOnCell (other % maxEdges, other % nCells))
   if (.not.allocated(self % cellsOnVertex)) allocate(self % cellsOnVertex(self % vertexDegree, self % nVertices))
   if (.not.allocated(self % cellsOnEdge)) allocate(self % cellsOnEdge (2, self % nEdges))
   if (.not.allocated(self % verticesOnEdge)) allocate(self % verticesOnEdge (2, self % nEdges))
   if (.not.allocated(self % dcEdge)) allocate(self % dcEdge (self % nEdges))
   if (.not.allocated(self % dvEdge)) allocate(self % dvEdge (self % nEdges))
   if (.not.allocated(self % kiteAreasOnVertex)) allocate(self % kiteAreasOnVertex(self % vertexDegree, self % nVertices))
   if (.not.allocated(self % edgesOnCell_sign)) allocate(self % edgesOnCell_sign(self % maxEdges, self % nCells))
   if (.not.allocated(self % areaTriangle)) allocate(self % areaTriangle(self % nVertices))
   if (.not.allocated(self % angleEdge)) allocate(self % angleEdge(self % nEdges))

   self % latCell           = other % latCell
   self % lonCell           = other % lonCell
   self % areaCell          = other % areaCell
   self % latEdge           = other % latEdge
   self % lonEdge           = other % lonEdge
   self % edgeNormalVectors = other % edgeNormalVectors
   self % zgrid             = other % zgrid
   self % nEdgesOnCell      = other % nEdgesOnCell
   self % edgesOnCell       = other % edgesOnCell
   self % cellsOnCell       = other % cellsOnCell
   self % cellsOnVertex     = other % cellsOnVertex
   self % cellsOnEdge       = other % cellsOnEdge
   self % verticesOnEdge    = other % verticesOnEdge
   self % dcEdge            = other % dcEdge
   self % dvEdge            = other % dvEdge
   self % kiteAreasOnVertex = other % kiteAreasOnVertex
   self % edgesOnCell_sign  = other % edgesOnCell_sign
   self % areaTriangle      = other % areaTriangle
   self % angleEdge         = other % angleEdge

!   write(*,*)'====> copy of geom corelist and domain'

   self % corelist => other % corelist
   self % domain   => other % domain
!   write(*,*)'inside geo_clone: other % domain % domainID = ', other % domain % domainID

   do ii = 1, size(geom_count)
      if (geom_count(ii)%id == self%domain%domainID) then
         geom_count(ii)%counter = geom_count(ii)%counter + 1
      end if
   end do

!   write(*,*)'====> copy of geom done'

end subroutine geo_clone

! ------------------------------------------------------------------------------

subroutine geo_delete(self)

   implicit none

   type(mpas_geom), intent(inout) :: self

   integer :: ii

   if (allocated(self%latCell)) deallocate(self%latCell)
   if (allocated(self%lonCell)) deallocate(self%lonCell)
   if (allocated(self%latEdge)) deallocate(self%latEdge)
   if (allocated(self%lonEdge)) deallocate(self%lonEdge)
   if (allocated(self%areaCell)) deallocate(self%areaCell)
   if (allocated(self%edgeNormalVectors)) deallocate(self%edgeNormalVectors)
   if (allocated(self%zgrid)) deallocate(self%zgrid)
   if (allocated(self%nEdgesOnCell)) deallocate(self%nEdgesOnCell)
   if (allocated(self%edgesOnCell)) deallocate(self%edgesOnCell)
   if (allocated(self%cellsOnCell)) deallocate(self%cellsOnCell)
   if (allocated(self%cellsOnVertex)) deallocate(self%cellsOnVertex)
   if (allocated(self%cellsOnEdge)) deallocate(self%cellsOnEdge)
   if (allocated(self%verticesOnEdge)) deallocate(self%verticesOnEdge)
   if (allocated(self%dcEdge)) deallocate(self%dcEdge)
   if (allocated(self%dvEdge)) deallocate(self%dvEdge)
   if (allocated(self%kiteAreasOnVertex)) deallocate(self%kiteAreasOnVertex)
   if (allocated(self%edgesOnCell_sign)) deallocate(self%edgesOnCell_sign)
   if (allocated(self%areaTriangle)) deallocate(self%areaTriangle)
   if (allocated(self%angleEdge)) deallocate(self%angleEdge)

   do ii = 1, size(geom_count)
      if (geom_count(ii)%id == self%domain%domainID) then
         geom_count(ii)%counter = geom_count(ii)%counter - 1
         if ((associated(self % corelist)).and.(associated(self % domain))) then
!            if (geom_count(ii)%counter == 0) then
!               ! Completely destroy corelist and domain
!               !  + can only do this if they are not used by another copy of self
!               !  + otherwise the persistance of the underlying domain
!               !    and corelist is a known memory leak
!!               write(*,*)'==> delete model corelist and domain'
!               call mpas_timer_set_context( self % domain )
!               TODO(JJG): divide mpas_finalize into smaller components to avoid
!               MPI finalize errors
!               call mpas_finalize(self % corelist, self % domain)
!            else
               nullify(self % corelist)
               nullify(self % domain)
!            end if
            exit
         end if
      end if
   end do

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
