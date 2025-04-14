! (C) Copyright 2017-2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module mpas_geom_mod

use atlas_module, only: atlas_field, atlas_fieldset, atlas_integer, atlas_real, atlas_functionspace
use fckit_configuration_module, only: fckit_configuration, fckit_YAMLConfiguration
use fckit_pathname_module, only: fckit_pathname
use fckit_log_module, only: fckit_log
use fckit_mpi_module, only: fckit_mpi_comm
use iso_c_binding

!oops
use oops_variables_mod, only: oops_variables
use kinds, only : kind_int, kind_real

!ufo
use ufo_vars_mod, only: MAXVARLEN, ufo_vars_getindex

!MPAS-Model
use mpas_derived_types
use mpas_kind_types
use mpas_constants
use mpas_dmpar, only: mpas_dmpar_sum_int, mpas_dmpar_sum_real
use mpas_subdriver
use atm_core
use mpas_pool_routines

!mpas_jedi
use mpas_constants_mod
use module_mp_thompson_cldfra3_saca, only: saca_param_type

#ifdef MPAS_EXTERNAL_ESMF_LIB
!external ESMF
use ESMF
#endif

implicit none
private
public :: mpas_geom, &
          geo_setup, geo_clone, geo_delete, geo_info, geo_is_equal, &
          geo_set_lonlat, geo_fill_geometry_fields, pool_has_field, &
          getSolveDimSizes, getSolveDimNames, getVertLevels, geo_vert_coord

public :: mpas_geom_registry


! ------------------------------------------------------------------------------

type :: templated_field
   character(len=MAXVARLEN) :: name
   character(len=MAXVARLEN) :: template
   character(len=MAXVARLEN) :: identity
   character(len=MAXVARLEN) :: io_name
   real(kind=kind_real)     :: io_scaling_factor
end type templated_field

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
   logical :: deallocate_nonda_fields
   character(len=StrKIND) :: bump_vunit
   real(kind=RKIND), dimension(:),   allocatable :: latCell, lonCell
   real(kind=RKIND), dimension(:),   allocatable :: areaCell
   real(kind=RKIND), dimension(:),   allocatable :: latEdge, lonEdge
   real(kind=RKIND), dimension(:,:), allocatable :: edgeNormalVectors
   real(kind=RKIND), dimension(:,:), allocatable :: zgrid
   real(kind=RKIND), dimension(:,:), allocatable :: height
   real(kind=RKIND), dimension(:,:), allocatable :: scaleheight
   integer, allocatable :: nEdgesOnCell(:)
   integer, allocatable :: cellsOnCell(:,:)
   integer, allocatable :: edgesOnCell(:,:)
   integer, allocatable :: cellsOnVertex(:,:)
   integer, allocatable :: cellsOnEdge(:,:)
   integer, allocatable :: verticesOnEdge(:,:)
   real(kind=RKIND), DIMENSION(:), ALLOCATABLE :: dcEdge, dvEdge
   real(kind=RKIND), DIMENSION(:), ALLOCATABLE :: areaTriangle, angleEdge
   real(kind=RKIND), DIMENSION(:,:), ALLOCATABLE :: kiteAreasOnVertex, edgesOnCell_sign
   logical :: is_regional = .false. ! global mesh by default
   integer, dimension(:), allocatable :: bdyMaskCell, bdyMaskEdge, bdyMaskVertex

   type (domain_type), pointer :: domain => null()
   type (core_type), pointer :: corelist => null()

   type(fckit_mpi_comm) :: f_comm

   type(atlas_functionspace) :: afunctionspace
   
   type(templated_field), allocatable :: templated_fields(:)

   integer :: iterator_dimension

   type (saca_param_type), public :: saca_params

   contains

   procedure, public :: is_templated => field_is_templated
   procedure, public :: template => template_fieldname
   procedure, public :: has_identity => field_has_identity
   procedure, public :: identity => identity_fieldname
   procedure, public :: has_io_name => field_has_io_name
   procedure, public :: io_name => io_name_fieldname
   procedure, public :: has_io_scaling_factor => field_has_io_scaling_factor
   procedure, public :: io_scaling_factor => io_scaling_factor_value
   procedure, public :: full_to_half => full_to_half_levels
   procedure, public :: get_num_nodes_and_elements
   procedure, public :: get_coords_and_connectivities
   generic, public :: nlevels => variables_nlevels, var_nlevels
   procedure :: variables_nlevels
   procedure :: var_nlevels


end type mpas_geom

type :: idcounter
   integer :: id, counter
end type idcounter

type(idcounter), allocatable :: geom_count(:)

character(len=1024) :: message

character(len=MAXVARLEN), parameter :: &
   MPASVerticalCoordinates(6) = &
      [character(len=MAXVARLEN) :: 'nVertLevels', 'nVertLevelsP1', &
         'nVertLevelsP2', 'nSoilLevels', 'nOznLevels', 'nAerLevels']

#define LISTED_TYPE mpas_geom

!> Linked list interface - defines registry_t type
#include <oops/util/linkedList_i.f>

!> Global registry
type(registry_t) :: mpas_geom_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include <oops/util/linkedList_c.f>

! ------------------------------------------------------------------------------
subroutine geo_setup(self, f_conf, f_comm)

   implicit none

   type(mpas_geom),           intent(inout) :: self
   type(fckit_configuration), intent(in)    :: f_conf
   type(fckit_mpi_comm),   intent(in)    :: f_comm

   character(len=StrKIND) :: string1
   type (mpas_pool_type), pointer :: meshPool, fg
   type (block_type), pointer :: block_ptr
   !character(len=120) :: fn
   character(len=512) :: nml_file, streams_file, fields_file
   character(len=:), allocatable :: str
   real(kind=kind_real) :: scaling_factor
   type(fckit_configuration) :: template_conf
   type(fckit_configuration), allocatable :: fields_conf(:)

   real (kind=RKIND), pointer :: r1d_ptr(:), r2d_ptr(:,:)
   integer, pointer :: i0d_ptr, i1d_ptr(:), i2d_ptr(:,:)
   logical, pointer :: config_apply_lbcs

   type(idcounter), allocatable :: prev_count(:)
   integer :: ii, nprev

   logical :: deallocate_fields
   logical :: bump_interp

   call fckit_log%info('==> create geom')

   ! MPI communicator
   self%f_comm = f_comm

   ! "nml_file" and "streams_file" are mandatory for mpas_init.
   call f_conf%get_or_die("nml_file",str)
   nml_file = str
   call f_conf%get_or_die("streams_file",str)
   streams_file = str

#ifdef MPAS_EXTERNAL_ESMF_LIB
   !external ESMF - initialize on behalf of MPAS; ok to hardcode calendar?
   call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN)
#endif

   ! Domain decomposition and templates for state/increment variables
   call mpas_init( self % corelist, self % domain, external_comm = self%f_comm%communicator(), &
                 & namelistFileParam = trim(nml_file), streamsFileParam = trim(streams_file))

   !Deallocate not-used fields for memory reduction
   call f_conf%get_or_die("deallocate non-da fields",deallocate_fields)
   self % deallocate_nonda_fields = deallocate_fields
   if (self % deallocate_nonda_fields) call geo_deallocate_nonda_fields (f_conf, self % domain)

   ! Set up the vertical coordinate for bump
   call f_conf%get_or_die("bump vunit",str)
   self % bump_vunit = str

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

   ! first read array of templated field configurations
   call f_conf%get_or_die('template fields file',str)
   fields_file = str
   template_conf = fckit_YAMLConfiguration(fckit_pathname(fields_file))
   call template_conf%get_or_die('fields',fields_conf)

   ! create templated field descriptors
   allocate(self%templated_fields(size(fields_conf)))
   do ii = 1, size(fields_conf)
      call fields_conf(ii)%get_or_die('field name',str)
      self%templated_fields(ii)%name = trim(str)
      call fields_conf(ii)%get_or_die('mpas template field',str)
      self%templated_fields(ii)%template = trim(str)
      if (fields_conf(ii)%get('mpas identity field',str)) then
        self%templated_fields(ii)%identity = trim(str)
      else
        self%templated_fields(ii)%identity = 'none'
      end if
      if (fields_conf(ii)%get('mpas io name',str)) then
        self%templated_fields(ii)%io_name = trim(str)
      else
        self%templated_fields(ii)%io_name = 'none'
      end if
      if (fields_conf(ii)%get('mpas io scaling factor',scaling_factor)) then
        self%templated_fields(ii)%io_scaling_factor = scaling_factor
      else
        self%templated_fields(ii)%io_scaling_factor = -1.0_kind_real
      end if
   end do
   deallocate(fields_conf)

  ! retrieve iterator dimension from config
  if ( .not. f_conf%get("iterator dimension", self%iterator_dimension) ) &
      self%iterator_dimension = 2

   call f_conf%get_or_die('l_build_madwrf',self%saca_params%l_build_madwrf)
   call f_conf%get_or_die('l_build_gsdcloud',self%saca_params%l_build_gsdcloud)
   call f_conf%get_or_die('l_saturate_qv',self%saca_params%l_saturate_qv)
   call f_conf%get_or_die('l_conserve_thetaV',self%saca_params%l_conserve_thetaV)
   call f_conf%get_or_die('cldfra_def',self%saca_params%cldfra_def)
   call f_conf%get_or_die('cldfra_thresh',self%saca_params%cldfra_thresh)
   call f_conf%get_or_die('cldmask_thresh',self%saca_params%cldmask_thresh)
   call f_conf%get_or_die('cld_bld_hgt',self%saca_params%cld_bld_hgt)
   if ( (self%saca_params%l_build_madwrf .and. self%saca_params%l_build_gsdcloud) .or. &
        (.not.self%saca_params%l_build_madwrf .and. .not.self%saca_params%l_build_gsdcloud) ) then
      write(message,*) '--> only one of l_build_madwrf OR l_build_gsdcloud should be set as .true. in mpas_geom%params%'
      call abor1_ftn(message)
   end if
   if ( (self%saca_params%l_saturate_qv .and. self%saca_params%l_conserve_thetaV) .or. &
        (.not.self%saca_params%l_saturate_qv .and. .not.self%saca_params%l_conserve_thetaV) ) then
      write(message,*) '--> only one of l_saturate_qv OR l_conserve_thetaV should be set as .true. in mpas_geom%params%'
      call abor1_ftn(message)
   end if
!   if (associated(self % domain)) then
!       write(message,*) 'inside geom: geom % domain associated for domainID = ', self % domain % domainID
!       call fckit_log%debug(message)
!   end if
!   if (associated(self % corelist)) then
!       call fckit_log%debug('inside geom: geom % corelist associated')
!   else
!       call fckit_log%debug('inside geom: geom % corelist not associated')
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
   allocate ( self % height ( self % nVertLevels, self % nCells ) )
   allocate ( self % scaleheight ( self % nVertLevels, self % nCells ) )
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
   allocate ( self % bdyMaskCell ( self % nCells ) )
   allocate ( self % bdyMaskEdge ( self % nEdges ) )
   allocate ( self % bdyMaskVertex ( self % nVertices ) )

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
   call mpas_pool_get_config( self%domain%blocklist%configs, 'config_apply_lbcs', config_apply_lbcs )
   self % is_regional = config_apply_lbcs
   call mpas_pool_get_array ( meshPool, 'bdyMaskCell', i1d_ptr )
   self % bdyMaskCell = i1d_ptr(1:self % nCells)
   call mpas_pool_get_array ( meshPool, 'bdyMaskEdge', i1d_ptr )
   self % bdyMaskEdge = i1d_ptr(1:self % nEdges)
   call mpas_pool_get_array ( meshPool, 'bdyMaskVertex', i1d_ptr )
   self % bdyMaskVertex = i1d_ptr(1:self % nVertices)

   call mpas_pool_get_array ( meshPool, 'zgrid', r2d_ptr )
   self % zgrid = r2d_ptr ( 1:self % nVertLevelsP1, 1:self % nCells )

   call self%full_to_half(self % zgrid, self % height, self % nCells)


   call fckit_log%debug('End of geo_setup')
   if (allocated(prev_count)) deallocate(prev_count)
   if (allocated(str)) deallocate(str)

end subroutine geo_setup

! --------------------------------------------------------------------------------------------------

subroutine full_to_half_levels(self, full, half, nx)

   implicit none

   class(mpas_geom), intent(in) :: self
   integer, intent(in) :: nx
   real (kind=RKIND), dimension(self%nVertLevels+1,nx), intent(in) :: full
   real (kind=RKIND), dimension(self%nVertLevels,nx), intent(out) :: half
   integer :: i, k

   !  calculate half level values of a variable at full levels:
   do k = 1, self%nVertLevels
      half(k,:) = ( full(k,:) + full(k+1,:) ) * MPAS_JEDI_HALF_kr
   enddo
end subroutine full_to_half_levels

! --------------------------------------------------------------------------------------------------

subroutine geo_set_lonlat(self, afieldset, include_halo)

   implicit none

   type(mpas_geom),  intent(inout) :: self
   type(atlas_fieldset), intent(inout) :: afieldset
   logical,              intent(in) :: include_halo
   
   !Locals
   real(kind_real), pointer :: real_ptr(:,:)
   type(atlas_field) :: afield, afield_incl_halo
   integer :: nx
   
   ! Create lon/lat field
   nx=self%nCellsSolve
   afield = atlas_field(name="lonlat", kind=atlas_real(kind_real), shape=(/2,nx/))   
   call afield%data(real_ptr)
   real_ptr(1,:) = self%lonCell(1:nx) * MPAS_JEDI_RAD2DEG_kr
   real_ptr(2,:) = self%latCell(1:nx) * MPAS_JEDI_RAD2DEG_kr
   call afieldset%add(afield)

   if (include_halo) then
      nullify(real_ptr)
      nx=self%nCells      
      ! Create an additional lon/lat field containing owned points (as above) and also halo
      afield_incl_halo = atlas_field(name="lonlat_including_halo", kind=atlas_real(kind_real), &
           shape=(/2,nx/))
      call afield_incl_halo%data(real_ptr)
      real_ptr(1,:) = self%lonCell(1:nx) * MPAS_JEDI_RAD2DEG_kr
      real_ptr(2,:) = self%latCell(1:nx) * MPAS_JEDI_RAD2DEG_kr
      call afieldset%add(afield_incl_halo)
   endif
   
end subroutine geo_set_lonlat

! --------------------------------------------------------------------------------------------------

subroutine geo_fill_geometry_fields(self, afieldset)

   implicit none

   type(mpas_geom),  intent(inout) :: self
   type(atlas_fieldset), intent(inout) :: afieldset

   integer :: i, iz, jz
   real(kind=kind_real), pointer :: real_ptr(:,:)
   real(kind=RKIND), pointer :: pressure_base(:,:)
   real(kind=RKIND) :: var_local, var_global
   type(atlas_field) :: afield

   integer :: nx, nx2
   integer, pointer :: int_ptr(:,:)
   ! subroutine for saber, always include halo
   
   nx = self%nCells        ! rank with halo
   nx2= self%nCellsSolve   !      without

   ! Add owned vs halo/BC field
   afield = self%afunctionspace%create_field &
        (name='owned', kind=atlas_integer(kind_int), levels=1)
   call afield%data(int_ptr)
   int_ptr(1,:)=0
   int_ptr(1,1:nx2) = 1
   call afieldset%add(afield)
   call afield%final()
   
   ! Add area
   afield = self%afunctionspace%create_field &
        (name='area', kind=atlas_real(kind_real), levels=1)
   call afield%data(real_ptr)
   real_ptr(1,1:nx) = real(self%areaCell(1:nx), kind_real)
   call afieldset%add(afield)
   call afield%final()

   
   ! Add vertical unit
   afield = self%afunctionspace%create_field &
        (name='vert_coord', kind=atlas_real(kind_real), levels=self%nVertLevels)
   call afield%data(real_ptr)
   do jz=1,self%nVertLevels
      iz = self%nVertLevels - jz + 1
      if (trim(self % bump_vunit) .eq. 'modellevel') then
         real_ptr(jz,1:nx) = real(iz, kind_real)
      else if (trim(self % bump_vunit) .eq. 'height') then
         real_ptr(jz,1:nx) = self%height(iz,1:nx)
      else if (trim(self % bump_vunit) .eq. 'avgheight') then
         ! find avgheight by averaging across MPI processes, of course without halo
         var_local = MPAS_JEDI_HALF_kr * ( sum(self%zgrid(iz,  1:nx2)) &
                                         + sum(self%zgrid(iz+1,1:nx2)) )
         call mpas_dmpar_sum_real(self%domain%dminfo,var_local,var_global)
         real_ptr(jz,1:nx) = real(var_global / self%nCellsGlobal , kind_real)
      else if (trim(self % bump_vunit) .eq. 'scaleheight' .or. trim(self % bump_vunit) .eq. 'logp') then
         ! defined as a natural logarithm of pressure_base (base state pressure)
         call mpas_pool_get_array(self % domain % blocklist % allFields,'pressure_base', pressure_base)
         real_ptr(jz,1:nx) = log(real(pressure_base(iz,1:nx), kind_real))
      else
         write(message,*) '--> geo_fill_geometry_fields: bump_vunit "'&
                          //trim(self % bump_vunit)//'" is not implemented'
         call abor1_ftn(message)
      end if
   end do

   call afieldset%add(afield)
   call afield%final()

end subroutine geo_fill_geometry_fields

! ------------------------------------------------------------------------------
subroutine geo_deallocate_nonda_fields(f_conf, domain)

   implicit none
   type(fckit_configuration),      intent(in)    :: f_conf
   type (domain_type), pointer,    intent(inout) :: domain
   type (mpas_pool_type), pointer                :: pool_a, pool_b
   type (mpas_pool_data_type), pointer           :: mem
   type (mpas_pool_iterator_type)                :: poolItr_b

   type (field0DReal), pointer :: field0d
   type (field1DReal), pointer :: field1d
   type (field2DReal), pointer :: field2d

   integer            :: i
   character (len=22), allocatable :: poolname_a(:), poolname_b(:)
   character (len=22), allocatable :: da_fieldnames(:)

   type(fckit_configuration) :: template_conf
   character(len=:), allocatable :: str, str_arr(:)
   character(len=512) :: fields_file

   allocate(poolname_a(3))
   poolname_a(1)='tend'
   poolname_a(2)='tend_physics'
   poolname_a(3)='atm_input'

   allocate(poolname_b(3))
   poolname_b(1)='diag_physics'
   poolname_b(2)='diag'
   poolname_b(3)='sfc_input'

   ! first read array of templated field configurations
   call f_conf%get_or_die('kept fields file',str)
   fields_file = str
   template_conf = fckit_YAMLConfiguration(fckit_pathname(fields_file))
   call template_conf%get_or_die('fields',str_arr)

   ! create da_fieldnames array to be kept
   allocate(da_fieldnames(size(str_arr)))
   do i = 1, size(str_arr)
      da_fieldnames(i)=trim(str_arr(i))
   end do
   deallocate(str_arr)

   do i=1,size(poolname_a)
      mem => pool_get_member(domain % blocklist % structs, poolname_a(i), MPAS_POOL_SUBPOOL)
      if (associated(mem)) then
         call mpas_pool_get_subpool(domain % blocklist % structs, poolname_a(i), pool_a)
         call mpas_pool_destroy_pool(pool_a)
         call mpas_pool_remove_subpool(domain % blocklist % structs, poolname_a(i))
      end if
   end do

   do i=1,size(poolname_b)

      mem => pool_get_member(domain % blocklist % structs, poolname_b(i), MPAS_POOL_SUBPOOL)
      if (associated(mem)) then
         call mpas_pool_get_subpool(domain % blocklist % structs, poolname_b(i), pool_b)
         call mpas_pool_begin_iteration(pool_b)

         do while ( mpas_pool_get_next_member(pool_b, poolItr_b) )

            if (poolItr_b % memberType == MPAS_POOL_FIELD .and. poolItr_b % dataType == MPAS_POOL_REAL .and. &
               ufo_vars_getindex(da_fieldnames, poolItr_b % memberName) < 1) then

               if (poolItr_b % nDims == 0) then
                  call mpas_pool_get_field(pool_b,trim(poolItr_b % memberName),field0d)
                  if (associated(field0d)) then
                     call mpas_deallocate_field(field0d)
                     call mpas_pool_remove_field(pool_b, trim(poolItr_b % memberName))
                  end if
                  nullify(field0d)
               else if (poolItr_b % nDims == 1) then
                  call mpas_pool_get_field(pool_b,trim(poolItr_b % memberName),field1d)
                  if (associated(field1d)) then
                     call mpas_deallocate_field(field1d)
                     call mpas_pool_remove_field(pool_b, trim(poolItr_b % memberName))
                  end if
                  nullify(field1d)
               else if (poolItr_b % nDims == 2) then
                  call mpas_pool_get_field(pool_b,trim(poolItr_b % memberName),field2d)
                  if (associated(field2d)) then
                     call mpas_deallocate_field(field2d)
                     call mpas_pool_remove_field(pool_b, trim(poolItr_b % memberName))
                  end if
                  nullify(field2d)
               endif

            end if

         end do

      end if

   end do

   deallocate(poolname_a)
   deallocate(poolname_b)
   deallocate(da_fieldnames)

end subroutine geo_deallocate_nonda_fields

! ------------------------------------------------------------------------------

subroutine geo_clone(self, other)

   implicit none

   type(mpas_geom), intent(inout) :: self
   type(mpas_geom), intent(in) :: other

   integer :: ii

   ! Clone communicator
   self%f_comm = other%f_comm

   call fckit_log%debug('====> copy of geom array')

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

   self % iterator_dimension = other % iterator_dimension

   self % saca_params = other % saca_params

   if (.not.allocated(self % latCell)) allocate(self % latCell(other % nCells))
   if (.not.allocated(self % lonCell)) allocate(self % lonCell(other % nCells))
   if (.not.allocated(self % latEdge)) allocate(self % latEdge(other % nEdges))
   if (.not.allocated(self % lonEdge)) allocate(self % lonEdge(other % nEdges))
   if (.not.allocated(self % areaCell)) allocate(self % areaCell(other % nCells))
   if (.not.allocated(self % edgeNormalVectors)) allocate(self % edgeNormalVectors(3, other % nEdges))
   if (.not.allocated(self % zgrid)) allocate(self % zgrid(other % nVertLevelsP1, other % nCells))
   if (.not.allocated(self % height)) allocate(self % height(other % nVertLevels, other % nCells))
   if (.not.allocated(self % scaleheight)) allocate(self % scaleheight(other % nVertLevels, other % nCells))
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
   if (.not.allocated(self % bdyMaskCell)) allocate ( self % bdyMaskCell ( self % nCells ) )
   if (.not.allocated(self % bdyMaskEdge)) allocate ( self % bdyMaskEdge ( self % nEdges ) )
   if (.not.allocated(self % bdyMaskVertex)) allocate ( self % bdyMaskVertex ( self % nVertices ) )

   self % templated_fields  = other % templated_fields
   self % latCell           = other % latCell
   self % lonCell           = other % lonCell
   self % areaCell          = other % areaCell
   self % latEdge           = other % latEdge
   self % lonEdge           = other % lonEdge
   self % edgeNormalVectors = other % edgeNormalVectors
   self % zgrid             = other % zgrid
   self % height            = other % height
   self % scaleheight       = other % scaleheight
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
   self % is_regional       = other % is_regional
   self % bdyMaskCell       = other % bdyMaskCell
   self % bdyMaskEdge       = other % bdyMaskEdge
   self % bdyMaskVertex     = other % bdyMaskVertex

   self%afunctionspace = atlas_functionspace(other%afunctionspace%c_ptr())

   call fckit_log%debug('====> copy of geom corelist and domain')

   self % corelist => other % corelist
   self % domain   => other % domain
   write(message,*) 'inside geo_clone: other % domain % domainID = ', other % domain % domainID
   call fckit_log%debug(message)

   do ii = 1, size(geom_count)
      if (geom_count(ii)%id == self%domain%domainID) then
         geom_count(ii)%counter = geom_count(ii)%counter + 1
      end if
   end do

   call fckit_log%debug('====> copy of geom done')

end subroutine geo_clone

! ------------------------------------------------------------------------------

subroutine geo_delete(self)

   implicit none

   type(mpas_geom), intent(inout) :: self
   integer :: ii

   if (allocated(self%templated_fields)) deallocate(self%templated_fields)
   if (allocated(self%latCell)) deallocate(self%latCell)
   if (allocated(self%lonCell)) deallocate(self%lonCell)
   if (allocated(self%latEdge)) deallocate(self%latEdge)
   if (allocated(self%lonEdge)) deallocate(self%lonEdge)
   if (allocated(self%areaCell)) deallocate(self%areaCell)
   if (allocated(self%edgeNormalVectors)) deallocate(self%edgeNormalVectors)
   if (allocated(self%zgrid)) deallocate(self%zgrid)
   if (allocated(self%height)) deallocate(self%height)
   if (allocated(self%scaleheight)) deallocate(self%scaleheight)
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
   if (allocated(self%bdyMaskCell)) deallocate(self%bdyMaskCell)
   if (allocated(self%bdyMaskEdge)) deallocate(self%bdyMaskEdge)
   if (allocated(self%bdyMaskVertex)) deallocate(self%bdyMaskVertex)

   do ii = 1, size(geom_count)
      if (geom_count(ii)%id == self%domain%domainID) then
         geom_count(ii)%counter = geom_count(ii)%counter - 1
         if ((associated(self % corelist)).and.(associated(self % domain))) then
            if (geom_count(ii)%counter == 0) then
               ! Completely destroy corelist and domain
               !  + can only do this if they are not used by another copy of self
               !  + otherwise the persistance of the underlying domain
               !    and corelist is a known memory leak
               write(message,'(A,I3)') &
                  '==> destruct MPAS corelist and domain:', self%domain%domainID
               call fckit_log%info(message)
               call mpas_timer_set_context( self % domain )
               call mpas_finalize(self % corelist, self % domain)
            else
               nullify(self % corelist)
               nullify(self % domain)
            end if
            exit
         end if
      end if
   end do
   
   call self%afunctionspace%final()

#ifdef MPAS_EXTERNAL_ESMF_LIB
   ! Finalize ESMF on behalf of MPAS, keep MPI alive for JEDI to clean up
   call ESMF_Finalize(endflag=ESMF_END_KEEPMPI)
#endif

end subroutine geo_delete

! ------------------------------------------------------------------------------

subroutine geo_is_equal(is_equal, self, other)

   implicit none

   logical*1,       intent(inout) :: is_equal
   type(mpas_geom), intent(in)    :: self
   type(mpas_geom), intent(in)    :: other

   ! Only compares two of many attributes of the two geometries, but
   ! sufficient for current needs.
   is_equal = (self%nCellsGlobal .eq. other%nCellsGlobal) .and. &
              (self%nVertLevels  .eq. other%nVertLevels)

end subroutine geo_is_equal

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

subroutine geo_vert_coord(self, nlevel, CoordName, VertCoord)

   implicit none

   type(mpas_geom), intent(in) :: self
   integer, intent(in) :: nlevel
   character(len=*), intent(in) :: CoordName
   real(kind=kind_real), intent(inout) :: VertCoord(nlevel)

   ! local variables
   integer :: iz, nx
   real(kind=RKIND) :: var_local, var_global
   real(kind=RKIND), pointer :: pressure_base(:,:)

   nx = self%nCellsSolve 

   if (trim(CoordName) .eq. "modellevel") then
      do iz = 1, nlevel
         VertCoord(iz) = real(iz, kind_real)
      end do
   else if (trim(CoordName) .eq. "height") then
      do iz = 1, nlevel
         var_local = MPAS_JEDI_HALF_kr * ( sum(self%zgrid(iz,  1:nx)) &
                                         + sum(self%zgrid(iz+1,1:nx)) )
         call mpas_dmpar_sum_real(self%domain%dminfo,var_local,var_global)
         VertCoord(iz) = real(var_global / self%nCellsGlobal , kind_real)
      end do
   else if (trim(CoordName) .eq. "scaleheight" .or. trim(CoordName) .eq. "logp") then
      call mpas_pool_get_array(self % domain % blocklist % allFields,'pressure_base', pressure_base)
      do iz = 1, nlevel
         var_local = sum(log(pressure_base(iz,1:nx)))
         call mpas_dmpar_sum_real(self%domain%dminfo,var_local,var_global)
         VertCoord(iz) = real(var_global / self%nCellsGlobal , kind_real)
      end do
   else
      write(message,*)'--> geo_vert_coord: ',trim(CoordName), &
                      ' not supported'
      call abor1_ftn(message)
   end if

end subroutine geo_vert_coord

! ------------------------------------------------------------------------------

function field_is_templated(self, fieldname) result(is_templated)
   implicit none
   class(mpas_geom), intent(in) :: self
   character(len=*), intent(in) :: fieldname
   integer :: ii
   logical :: is_templated
   is_templated = .false.
   if (allocated(self%templated_fields)) then
      do ii = 1, size(self%templated_fields)
         if (trim(self%templated_fields(ii)%name) == trim(fieldname)) then
            is_templated = .true.
            return
         end if
      end do
   end if
end function field_is_templated

! ------------------------------------------------------------------------------

function template_fieldname(self, fieldname) result(template)
   implicit none
   class(mpas_geom), intent(in) :: self
   character(len=*), intent(in) :: fieldname
   integer :: ii
   character(len=MAXVARLEN) :: template
   do ii = 1, size(self%templated_fields)
      if (trim(self%templated_fields(ii)%name) == trim(fieldname)) then
         template = self%templated_fields(ii)%template
         return
      end if
   end do
   write(message,*) '--> mpas_geom % template_fieldname: fieldname is not templated, ',fieldname
   call abor1_ftn(message)
end function template_fieldname

! ------------------------------------------------------------------------------

function field_has_identity(self, fieldname) result(has_identity)
   implicit none
   class(mpas_geom), intent(in) :: self
   character(len=*), intent(in) :: fieldname
   integer :: ii
   logical :: has_identity
   has_identity = .false.
   if (allocated(self%templated_fields)) then
      do ii = 1, size(self%templated_fields)
         if (trim(self%templated_fields(ii)%name) == trim(fieldname)) then
            if (trim(self%templated_fields(ii)%identity) /= 'none') then
               has_identity = .true.
            end if
            return
         end if
      end do
   end if
end function field_has_identity

! ------------------------------------------------------------------------------

function identity_fieldname(self, fieldname) result(identity)
   implicit none
   class(mpas_geom), intent(in) :: self
   character(len=*), intent(in) :: fieldname
   integer :: ii
   character(len=MAXVARLEN) :: identity
   do ii = 1, size(self%templated_fields)
      if (trim(self%templated_fields(ii)%name) == trim(fieldname)) then
         identity = self%templated_fields(ii)%identity
         return
      end if
   end do
   write(message,*) '--> mpas_geom % identity_fieldname: fieldname does not have identity, ',fieldname
   call abor1_ftn(message)
end function identity_fieldname

! ------------------------------------------------------------------------------

function field_has_io_name(self, fieldname) result(has_io_name)
   implicit none
   class(mpas_geom), intent(in) :: self
   character(len=*), intent(in) :: fieldname
   integer :: ii
   logical :: has_io_name
   has_io_name = .false.
   if (allocated(self%templated_fields)) then
      do ii = 1, size(self%templated_fields)
         if (trim(self%templated_fields(ii)%name) == trim(fieldname)) then
            if (trim(self%templated_fields(ii)%io_name) /= 'none') then
               has_io_name = .true.
            end if
            return
         end if
      end do
   end if
end function field_has_io_name

! ------------------------------------------------------------------------------

function io_name_fieldname(self, fieldname) result(io_name)
   implicit none
   class(mpas_geom), intent(in) :: self
   character(len=*), intent(in) :: fieldname
   integer :: ii
   character(len=MAXVARLEN) :: io_name
   do ii = 1, size(self%templated_fields)
      if (trim(self%templated_fields(ii)%name) == trim(fieldname)) then
         io_name = self%templated_fields(ii)%io_name
         return
      end if
   end do
   write(message,*) '--> mpas_geom % io_name_fieldname: fieldname does not have io_name, ',fieldname
   call abor1_ftn(message)
end function io_name_fieldname

! ------------------------------------------------------------------------------
function field_has_io_scaling_factor(self, fieldname) result(has_io_scaling_factor)
   implicit none
   class(mpas_geom), intent(in) :: self
   character(len=*), intent(in) :: fieldname
   integer :: ii
   logical :: has_io_scaling_factor
   has_io_scaling_factor = .false.
   if (allocated(self%templated_fields)) then
      do ii = 1, size(self%templated_fields)
         if (trim(self%templated_fields(ii)%name) == trim(fieldname)) then
            if (self%templated_fields(ii)%io_scaling_factor .gt. 0.0_kind_real) then
               has_io_scaling_factor = .true.
            end if
            return
         end if
      end do
   end if
end function field_has_io_scaling_factor

! ------------------------------------------------------------------------------

function io_scaling_factor_value(self, fieldname) result(io_scaling_factor)
   implicit none
   class(mpas_geom), intent(in) :: self
   character(len=*), intent(in) :: fieldname
   integer :: ii
   real(kind=kind_real) :: io_scaling_factor
   do ii = 1, size(self%templated_fields)
      if (trim(self%templated_fields(ii)%name) == trim(fieldname)) then
         io_scaling_factor = self%templated_fields(ii)%io_scaling_factor
         return
      end if
   end do
   write(message,*) '--> mpas_geom % io_scaling_factor_value: fieldname does not have io_scaling_factor, ',fieldname
   call abor1_ftn(message)
end function io_scaling_factor_value

! ------------------------------------------------------------------------------

!***********************************************************************
!
!  subroutine pool_has_field
!
!> \brief Check for presence of field in pool
!
!-----------------------------------------------------------------------
function pool_has_field(pool, field) result(has)

   implicit none

   type(mpas_pool_type), pointer, intent(in) :: pool
   character(len=*), intent(in) :: field
   type(mpas_pool_iterator_type) :: poolItr
   logical :: has

   has = .false.

   call mpas_pool_begin_iteration(pool)

   do while ( mpas_pool_get_next_member(pool, poolItr) )
      if (poolItr % memberType == MPAS_POOL_FIELD .and. &
          trim(field) == trim(poolItr % memberName) ) then
         has = .true.
         exit
      endif
   end do

end function pool_has_field


!*******************************************************************************
!
! function isDecomposed
!
!> \brief whether horizontal dimension name is decomposed across MPI tasks
!
!*******************************************************************************

logical function isDecomposed(dimName)
   implicit none
   character(len=*), intent(in) :: dimName
   if (trim(dimName) == 'nCells') then
      isDecomposed = .true.
      return
   else if (trim(dimName) == 'nEdges') then
      isDecomposed = .true.
      return
   else if (trim(dimName) == 'nVertices') then
      isDecomposed = .true.
      return
   else
      isDecomposed = .false.
      return
   end if
end function isDecomposed


!***********************************************************************
!
!  function getSolveDimNames
!
!> \brief   Returns an array with the dimension names of a pool field member
!
!-----------------------------------------------------------------------

function getSolveDimNames(pool, key) result(solveDimNames)

   implicit none

   ! Arguments
   type(mpas_pool_type), pointer, intent(in) :: pool
   character(len=*), intent(in) :: key

   character (len=StrKIND), allocatable :: solveDimNames(:)

   type(mpas_pool_data_type), pointer :: dptr
   character (len=StrKIND) :: dimName
   integer :: id

   dptr => pool_get_member(pool, key, MPAS_POOL_FIELD)

   if (.not. associated(dptr)) then
      write(message,*)'--> getSolveDimNames: problem finding ',trim(key), &
                  ' MPAS_POOL_FIELD in pool'
      call abor1_ftn(message)
   end if

   allocate(solveDimNames(dptr%contentsDims))

   do id = 1, dptr%contentsDims

      ! fields with only one time level
      !real
      if (associated(dptr%r1)) then
         dimName = trim(dptr%r1%dimNames(id))
      else if (associated(dptr%r2)) then
         dimName = trim(dptr%r2%dimNames(id))
      else if (associated(dptr%r3)) then
         dimName = trim(dptr%r3%dimNames(id))
      !integer
      else if (associated(dptr%i1)) then
         dimName = trim(dptr%i1%dimNames(id))
      else if (associated(dptr%i2)) then
         dimName = trim(dptr%i2%dimNames(id))
      else if (associated(dptr%i3)) then
         dimName = trim(dptr%i3%dimNames(id))

      ! fields with multiple time levels (use first time level)
      !real
      else if (associated(dptr%r1a)) then
         dimName = trim(dptr%r1a(1)%dimNames(id))
      else if (associated(dptr%r2a)) then
         dimName = trim(dptr%r2a(1)%dimNames(id))
      else if (associated(dptr%r3a)) then
         dimName = trim(dptr%r3a(1)%dimNames(id))
      !integer
      else if (associated(dptr%i1a)) then
         dimName = trim(dptr%i1a(1)%dimNames(id))
      else if (associated(dptr%i2a)) then
         dimName = trim(dptr%i2a(1)%dimNames(id))
      else if (associated(dptr%i3a)) then
         dimName = trim(dptr%i3a(1)%dimNames(id))

      else
         write(message,*) '--> getSolveDimNames: unsupported field type for key,&
            & contentsDims, contentsType = ',trim(key),',',dptr%contentsDims, &
            ',',dptr%contentsType
         call abor1_ftn(message)
      end if

      if (isDecomposed(dimName)) then
         solveDimNames(id) = trim(dimName)//'Solve'
      else
         solveDimNames(id) = trim(dimName)
      end if
   end do
end function getSolveDimNames


!***********************************************************************
!
!  function getSolveDimSizes
!
!> \brief   Returns an array with the dimension sizes of a pool field member
!
!> \details
!    In some cases, the pool containing the field will not contain any
!    dimensions.  Then the optional dimPool argument can be used.
!
!-----------------------------------------------------------------------

function getSolveDimSizes(pool, key, dimPool) result(solveDimSizes)

   implicit none

   ! Arguments
   type (mpas_pool_type), pointer, intent(in) :: pool
   character(len=*), intent(in) :: key
   type (mpas_pool_type), pointer, optional, intent(in) :: dimPool

   integer, allocatable :: solveDimSizes(:)

   character (len=StrKIND), allocatable :: dimNames(:)
   integer :: id
   integer, pointer :: dimSize

   dimNames = getSolveDimNames(pool, key)

   allocate(solveDimSizes(size(dimNames)))

   do id = 1, size(dimNames)
      if (present(dimPool)) then
         call mpas_pool_get_dimension(dimPool, dimNames(id), dimSize)
      else
         call mpas_pool_get_dimension(pool, dimNames(id), dimSize)
      end if

      if (associated(dimSize)) then
         solveDimSizes(id) = dimSize
      else
         write(message,*)'--> getSolveDimSizes: ',trim(dimNames(id)), &
                     ' not available'
         call abor1_ftn(message)
      end if
   end do

   deallocate(dimNames)

end function getSolveDimSizes


!***********************************************************************
!
!  function getVertLevels
!
!> \brief   Returns the number of vertical levels for a pool field member
!
!> \details
!    In some cases, the pool containing the field will not contain any
!    dimensions.  Then the optional dimPool argument should be used.
!    See var_nlevels for an example.
!
!-----------------------------------------------------------------------

function getVertLevels(pool, key, dimPool) result(nlevels)

   ! Arguments
   type (mpas_pool_type), pointer, intent(in) :: pool
   character(len=*), intent(in) :: key
   type (mpas_pool_type), pointer, optional, intent(in) :: dimPool

   integer :: nlevels, id
   integer, allocatable :: dimSizes(:)
   character (len=StrKIND), allocatable :: dimNames(:)

   ! assume default value of 1, unless stated otherwise in field attributes
   nlevels = 1

   dimNames = getSolveDimNames(pool, key)
   if (present(dimPool)) then
      dimSizes = getSolveDimSizes(pool, key, dimPool)
   else
      dimSizes = getSolveDimSizes(pool, key)
   end if
   do id = 1, size(dimSizes)
      ! check if dimNames(id) is one of the available vertical coordinates
      if (ufo_vars_getindex(MPASVerticalCoordinates, dimNames(id)) > 0) then
         nlevels = dimSizes(id)
      end if
   end do

   deallocate(dimNames, dimSizes)

end function getVertLevels

! ------------------------------------------------------------------------------

subroutine var_nlevels(self, var, nlevels)

   implicit none

   class(mpas_geom), intent(in) :: self
   character(len=*), intent(in) :: var
   integer, intent(out) :: nlevels

   character(len=StrKIND) :: poolVar
   type(mpas_pool_type), pointer :: allFields

   allFields => self % domain % blocklist % allFields

   if (pool_has_field(allFields, var)) then
      poolVar = var
   else if (self % is_templated(var)) then
      poolVar = self%template(var)
      if (.not.pool_has_field(allFields, poolVar)) then
         write(message,*)'--> vars_nlevels: ',trim(poolVar), &
                     ' not available in MPAS domain'
         call abor1_ftn(message)
      end if
   else
      write(message,*)'--> vars_nlevels: ',trim(var), &
                  ' not available in MPAS domain or self % templated_fields'
      call abor1_ftn(message)
   end if

   nlevels = getVertLevels(allFields, poolVar, self % domain % blocklist % dimensions)

end subroutine var_nlevels

! ------------------------------------------------------------------------------

subroutine variables_nlevels(self, vars, nv, nlevels)

   implicit none

   class(mpas_geom), intent(in) :: self
   type(oops_variables), intent(in) :: vars
   integer(c_size_t), intent(in) :: nv
   integer(c_size_t), intent(inout) :: nlevels(nv)

   integer :: ivar, nn

   do ivar = 1, nv
      call self%nlevels(vars%variable(ivar), nn)
      nlevels(ivar) = int(nn, c_size_t)
   end do

end subroutine variables_nlevels

! ------------------------------------------------------------------------------

subroutine get_num_nodes_and_elements(self, num_nodes, num_tris)

   implicit none

   class(mpas_geom), intent(in) :: self
   integer, intent(out) :: num_nodes
   integer, intent(out) :: num_tris

   integer :: nVerticesBdy7

   num_nodes = self % nCells         ! Local + Halo
   num_tris = self % nVerticesSolve  ! Local

   ! for regional MPAS mesh
   if ( self % is_regional ) then
      ! count the number of outer-most Vertices ( bdyMaskVertex == 7 )
      nVerticesBdy7 = count ( self%bdyMaskVertex(1:self%nVerticesSolve) ==7 )
      ! number of "interior" own vertices for a given process
      num_tris = self % nVerticesSolve - nVerticesBdy7
      write(message,*) self%f_comm%rank(), 'this is regional, nVerticesSolve, nVerticesBdy7=', &
              self % nVerticesSolve, nVerticesBdy7
      call fckit_log%debug(message)
   end if

end subroutine get_num_nodes_and_elements

! ------------------------------------------------------------------------------

subroutine get_coords_and_connectivities(self, &
   num_nodes, num_tri_boundary_nodes, &
   lons, lats, ghosts, global_indices, remote_indices, partition, &
   raw_tri_boundary_nodes)

   implicit none

   class(mpas_geom), intent(in) :: self
   integer, intent(in) :: num_nodes
   integer, intent(in) :: num_tri_boundary_nodes
   real(kind_real), intent(out) :: lons(num_nodes)
   real(kind_real), intent(out) :: lats(num_nodes)
   integer, intent(out) :: ghosts(num_nodes)
   integer, intent(out) :: global_indices(num_nodes)
   integer, intent(out) :: remote_indices(num_nodes)
   integer, intent(out) :: partition(num_nodes)
   integer, intent(out) :: raw_tri_boundary_nodes(num_tri_boundary_nodes)

   integer :: i, iVertValid
   type (field1DInteger), pointer :: indexToCellID, indexToVertexID, iTmp

   lons = self % lonCell * MPAS_JEDI_RAD2DEG_kr
   lats = self % latCell * MPAS_JEDI_RAD2DEG_kr

   ghosts = 1
   ghosts(1:self%nCellsSolve) = 0 ! Own nodes (or Cells in MPAS)

   !indexToCellID & indexToVertexID from MPAS domain
   call mpas_pool_get_field(self%domain%blocklist%allFields, 'indexToCellID', indexToCellID)
   call mpas_pool_get_field(self%domain%blocklist%allFields, 'indexToVertexID', indexToVertexID)
   call mpas_duplicate_field(indexToCellID, iTmp)   ! intermediate for halo exchange

   global_indices(1:num_nodes) = indexToCellID % array(1:num_nodes)
   !No need halo exchange. b/c using self%nCells here.

   iTmp % array(:) = -1
   do i=1,self%nCellsSolve ! Own definition first
      iTmp % array(i) = i
   end do
   call mpas_dmpar_exch_halo_field(iTmp) ! halo exchange
   do i=1,num_nodes !=self%nCells
      remote_indices(i) = iTmp % array(i)
   end do

   iTmp % array(:) = -1
   do i=1,self%nCellsSolve ! Own definition first
      iTmp % array(i) = self%f_comm%rank()
   end do
   call mpas_dmpar_exch_halo_field(iTmp) ! halo exchange
   do i=1,num_nodes !=self%nCells
      partition(i) = iTmp % array(i)
   end do

   !clean up 
   call mpas_deallocate_field(iTmp)

   iVertValid=1 ! used for global indices of vertices.
   do i=1,self%nVerticesSolve ! Note self%nVerticesSolve = num_tris (for global mesh), but != num_tris (for regional mesh)
      if( self % bdyMaskVertex(i) .eq. 7 ) cycle ! regional MPAS mesh, Skip the outer-most Vertices
      raw_tri_boundary_nodes(3*(iVertValid-1)+1) = indexToCellID % array( self%cellsOnVertex(1,i) )
      raw_tri_boundary_nodes(3*(iVertValid-1)+2) = indexToCellID % array( self%cellsOnVertex(2,i) )
      raw_tri_boundary_nodes(3*(iVertValid-1)+3) = indexToCellID % array( self%cellsOnVertex(3,i) )
      iVertValid=iVertValid+1
   end do

end subroutine get_coords_and_connectivities

! ------------------------------------------------------------------------------

end module mpas_geom_mod
