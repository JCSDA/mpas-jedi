! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_geom_mod

use atlas_module, only: atlas_functionspace, atlas_fieldset, atlas_field, atlas_real
use fckit_configuration_module, only: fckit_configuration
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
use mpas_run_mod, only : run_corelist=>corelist, run_domain=>domain

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
   real(kind=kind_real), DIMENSION(:),   ALLOCATABLE :: latCell, lonCell
   real(kind=kind_real), DIMENSION(:),   ALLOCATABLE :: areaCell
   real(kind=kind_real), DIMENSION(:),   ALLOCATABLE :: latEdge, lonEdge
   real(kind=kind_real), DIMENSION(:,:), ALLOCATABLE :: edgeNormalVectors
   real(kind=kind_real), DIMENSION(:,:), ALLOCATABLE :: zgrid
   integer, allocatable :: nEdgesOnCell(:)
   integer, allocatable :: cellsOnCell(:,:)
   integer, allocatable :: edgesOnCell(:,:)

   type (domain_type), pointer :: domain => null() 
   type (core_type), pointer :: corelist => null()
   type(atlas_functionspace) :: afunctionspace
end type mpas_geom

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
subroutine geo_setup(self, f_conf)

   implicit none

   type(mpas_geom),           intent(inout) :: self
   type(fckit_configuration), intent(in)    :: f_conf

   character(len=StrKIND) :: string1
   type (mpas_pool_type), pointer :: meshPool, fg
   type (block_type), pointer :: block_ptr

   real (kind=kind_real), pointer :: r1d_ptr(:), r2d_ptr(:,:)
   integer, pointer :: i0d_ptr, i1d_ptr(:), i2d_ptr(:,:)

!   write(*,*)' ==> create geom'

   self % corelist => run_corelist
   self % domain   => run_domain

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


!  Could make this more flexible/clean by using pointers for array variables in mpas_geom 
!  + any later value modifications need to be consistent with MPAS model (e.g., no unit conversions)
!  + would need to nullify instead of allocate/dellocate
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
   end do
   call afieldset%add(afield)
   call afield%final()

end subroutine geo_fill_atlas_fieldset

! ------------------------------------------------------------------------------

subroutine geo_clone(self, other)

   implicit none

   type(mpas_geom), intent(in) :: self
   type(mpas_geom), intent(inout) :: other

!   write(*,*)'====> copy of geom array'
!   if (allocated(other%latCell)) then 
!      write(*,*)'Allocated array other%latCell'
!   else
!      write(*,*)'Not Allocated array other%latCell'
!   end if   

   other % nCellsGlobal  = self % nCellsGlobal
   other % nCells        = self % nCells
   other % nCellsSolve   = self % nCellsSolve

   other % nEdgesGlobal  = self % nEdgesGlobal
   other % nEdges        = self % nEdges
   other % nEdgesSolve   = self % nEdgesSolve

   other % nVerticesGlobal  = self % nVerticesGlobal
   other % nVertices        = self % nVertices
   other % nVerticesSolve   = self % nVerticesSolve

   other % nVertLevels   = self % nVertLevels
   other % nVertLevelsP1 = self % nVertLevelsP1
   other % nSoilLevels   = self % nSoilLevels 
   other % vertexDegree  = self % vertexDegree
   other % maxEdges      = self % maxEdges

   if (.not.allocated(other % latCell)) allocate(other % latCell(self % nCells))
   if (.not.allocated(other % lonCell)) allocate(other % lonCell(self % nCells))
   if (.not.allocated(other % latEdge)) allocate(other % latEdge(self % nEdges))
   if (.not.allocated(other % lonEdge)) allocate(other % lonEdge(self % nEdges))
   if (.not.allocated(other % areaCell)) allocate(other % areaCell(self % nCells))
   if (.not.allocated(other % edgeNormalVectors)) allocate(other % edgeNormalVectors(3, self % nEdges))
   if (.not.allocated(other % zgrid)) allocate(other % zgrid(self % nVertLevelsP1, self % nCells))
   if (.not.allocated(other % nEdgesOnCell)) allocate(other % nEdgesOnCell(self % nCells))
   if (.not.allocated(other % edgesOnCell)) allocate(other % edgesOnCell(self % maxEdges, self % nCells))
   if (.not.allocated(other % cellsOnCell)) allocate (other % cellsOnCell ( self % maxEdges, self % nCells ) )

   other % latCell           = self % latCell
   other % lonCell           = self % lonCell
   other % areaCell          = self % areaCell
   other % latEdge           = self % latEdge
   other % lonEdge           = self % lonEdge
   other % edgeNormalVectors = self % edgeNormalVectors
   other % zgrid             = self % zgrid
   other % nEdgesOnCell      = self % nEdgesOnCell
   other % edgesOnCell       = self % edgesOnCell
   other % cellsOnCell       = self % cellsOnCell

!   write(*,*)'====> copy of geom corelist and domain'

   if ((associated(other % corelist)).and.(associated(other % domain))) then 
!      write(*,*)'associated(other % corelist), associated(other % domain)'
   else
!      write(*,*)'not associated(other % corelist), associated(other % domain)'
      other % corelist => run_corelist
      other % domain   => run_domain
   end if
!   write(*,*)'inside geo_clone: other % domain % domainID = ', other % domain % domainID

!   write(*,*)'====> copy of geom done'

end subroutine geo_clone

! ------------------------------------------------------------------------------

subroutine geo_delete(self)

   implicit none

   type(mpas_geom), intent(inout) :: self

!   write(*,*)'==> delete geom array'
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

   !call mpas_timer_set_context( self % domain )
   if ((associated(self % corelist)).and.(associated(self % domain))) then
      nullify(self % corelist)
      nullify(self % domain)
!      write(*,*)'==> nullify geom corelist and domain'
   end if
!   write(*,*)'==> delete geom done'

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
