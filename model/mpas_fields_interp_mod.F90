! Copyright (c) 2018, National Atmospheric for Atmospheric Research (NCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html

module mpas_fields_interp_mod

!***********************************************************************
!
!  Module mpas_fields_interp_mod contains all the utilities for mpas
!  interpolation at observation locations
!  NICAS is needed 
!  or from /somewhere/mpas-bundle/mpas/model (OOPS) 
!> \author  Bjung J. / Gael Descombes NCAR/MMMM
!> \date    January 2018
!
!-----------------------------------------------------------------------

use config_mod
use mpas_geom_mod
use mpas_vars_mod
!use mpas_locs_mod --> UFO not ready yet
use kinds
use tools_nc
use netcdf
use mpas_constants, only : pii
! call to ufo module
use ufo_locs_mod
use ufo_geovals_mod
use mpas_kind_types

use mpas_derived_types
use mpas_framework
use mpas_kind_types
!use init_atm_core_interface
use atm_core

implicit none

#define LISTED_TYPE mpas_field

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_field_registry
real(kind=kind_real), parameter :: rad2deg = 180./pii
real(kind=kind_real), parameter :: deg2rad = pii/180.

!> Fortran derived type to hold MPAS fields
#include "mpas_fields_oops_type.inc"


!-- DART Interp. test begin
type location_type
   private
   real(kind=kind_real) :: lon, lat        ! lon, lat are stored in radians
   real(kind=kind_real) :: vloc            ! units vary based on value of which_vert
   integer  :: which_vert      ! determines if vert is level, height, pressure, ...
end type location_type

real(kind=kind_real), parameter :: MISSING_R8 = -9999.99

type xyz_location_type
   private
   real(kind=kind_real) :: x, y, z
end type xyz_location_type

type box_type
   private
   integer, allocatable  :: obs_box(:)           ! (nobs); List of obs indices in boxes
   integer, allocatable  :: count(:, :, :)       ! (nx, ny, nz); # of obs in each box
   integer, allocatable  :: start(:, :, :)       ! (nx, ny, nz); Start of list of obs in this box
   real(kind=kind_real)          :: bot_x, top_x         ! extents in x, y, z
   real(kind=kind_real)          :: bot_y, top_y
   real(kind=kind_real)          :: bot_z, top_z
   real(kind=kind_real)          :: x_width, y_width, z_width    ! widths of boxes in x,y,z
   real(kind=kind_real)          :: nboxes_x, nboxes_y, nboxes_z ! based on maxdist how far to search
end type box_type

type xyz_get_close_type
   private
   integer           :: num
   real(kind=kind_real)          :: maxdist
   type(box_type)    :: box
end type xyz_get_close_type

real(kind=kind_real), parameter :: radius = 6.371e6 !1.0 !6.371e6 !6371000.0 !6371229.0 ! meters
integer, parameter:: nCells = 40962, nEdges = 122880, nVertices = 81920, maxEdges = 10, maxEdges2 = 20
real(kind=kind_real), parameter :: roundoff = 1.0e-12_kind_real
integer, DIMENSION(:), ALLOCATABLE :: nEdgesOnCell     !(nCells)
integer, DIMENSION(:,:), ALLOCATABLE :: verticesOnCell   !(nCells, maxEdges)
real(kind=kind_real), DIMENSION(:),   ALLOCATABLE :: xVertex, yVertex, zVertex  !(nVertices)
integer, DIMENSION(:,:), ALLOCATABLE :: edgesOnCell !(nCells, maxEdges)
integer, DIMENSION(:,:), ALLOCATABLE :: cellsOnEdge !(nEdges, TWO)
integer :: nx               = 20
integer :: ny               = 20
integer :: nz               = 20
interface xyz_set_location
!   module procedure set_location_single
!   module procedure set_location_array
   module procedure set_location_lonlat
end interface xyz_set_location

!-- DART Interp. test end

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

#define TEST_INTERP
#ifdef TEST_INTERP
! ------------------------------------------------------------------------------
! Routine to test interpolation inside MPAS, based on read_file
! ------------------------------------------------------------------------------

subroutine mpas_test_interp(fld, c_conf, vdate)

  use iso_c_binding
  use datetime_mod
  use fckit_log_module, only : log

  implicit none
  type(mpas_field), intent(inout) :: fld      !< Fields
  type(c_ptr),      intent(in)    :: c_conf   !< Configuration
  type(datetime),   intent(inout) :: vdate    !< DateTime

  character(len=StrKIND+50) :: record
  character(len=StrKIND) :: filename, string1
  character(len=StrKIND)  :: varname
  character(len=StrKIND)  :: sdate  
  character(len=1024)  :: buf
  !real(kind=kind_real), allocatable :: zz(:)
  integer :: i, k, iread, nf
  integer :: itime, istrlen
  integer :: ncid, dimid, varid
  integer :: ncells, nlevels, ndims !, nedges
  !BJJ NICAS-H
  type(ufo_locs)   :: locs
  type(ufo_geoval) :: gom
  !BJJ DART Interp. test
  real(kind=RKIND), allocatable :: x_dart(:), dval(:), x_nicas(:,:), x_func(:,:,:)
  integer :: nfunc = 2 
  integer :: iobs, ier
  integer, allocatable :: ival(:)
  type(location_type) :: locs_single
  real(kind=RKIND) :: weights(3)

  iread = 1
  if (config_element_exists(c_conf,"read_from_file")) then
     iread = config_get_int(c_conf,"read_from_file")
  end if
  filename = config_get_string(c_conf,len(filename),"filename")
  WRITE(buf,*) 'mpas_field:read_file: opening '//filename
  call log%info(buf)
  string1 = filename
  call ncerr(string1, nf90_open(trim(filename),nf90_nowrite,ncid))

  !> Grid dimensions
  call ncerr(string1, nf90_inq_dimid(ncid,'nCells',dimid))
  call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=ncells)) !self%nCells))
  call ncerr(string1, nf90_inq_dimid(ncid,'nVertLevels',dimid))
  call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=nlevels)) !self%nVertLevels))

  !> Check dimensions between mesh and fields
  if (ncells /= fld%geom%nCells .or.  nlevels /= fld%geom%nVertLevels) then
      write (record,*) "mpas_fields:read_file: wrong dimensions: ",fld%geom%nCells,fld%geom%nVertLevels
      call log%error(record)
      write (record,*) "mpas_fields:read_file: expected dimensions: ",ncells,nlevels
      call log%error(record)
      call abor1_ftn("mpas_fields:read_file: input fields have wrong dimensions")
  end if

  !> Time info
  call ncerr(string1, nf90_inq_varid(ncid,'xtime',varid))
  call ncerr(string1, nf90_inq_dimid(ncid, "Time", dimid))
  call ncerr(string1, nf90_inquire_dimension(ncid, dimid, len=itime))
  call ncerr(string1, nf90_inq_dimid(ncid, "StrLen", dimid))
  call ncerr(string1, nf90_inquire_dimension(ncid, dimid, len=istrlen))
  if (itime /= 1) then
      write(record,*) 'multiple timesteps (',itime,') in file ', trim(filename)
      write(record,*) 'We are using the LAST one, presumably, the LATEST timestep.'
      call log%error(record)
  end if

  call ncerr(string1, nf90_get_var(ncid, varid, sdate(1:istrlen), start = (/1, itime/), count = (/istrlen, 1 /)))
  sdate(11:11) = 'T'
  sdate(20:20) = 'Z'
  WRITE(buf,*) 'validity date is: '//trim(sdate)
  call log%info(buf)
  call datetime_set(sdate(1:20), vdate)

  !> Read variables fld (nCells, nVertLevels, nVariables)

  !do i=1,fld % nf
  !   call ncerr(string1, nf90_inq_varid(ncid, trim(fld%fldnames(i)), varid))
  !   do k=1,fld%geom%nVertLevels
  !      call ncerr(string1, nf90_get_var  (ncid, varid, fld%fld(:,k,i), start = (/k,1,1/), count = (/1,fld%geom%nCells,1/) ))
  !   end do
  !end do

  !--- BJJ TEST Function : f = lat + lon [degrees], which is totally linear
  allocate(x_func(ncells,nfunc))
  x_func(:,1) = fld%geom%latCell(:) + fld%geom%lonCell(:)
  x_func(:,2) = (cos(deg2rad*fld%geom%lonCell(:)) + sin(deg2rad*fld%geom%lonCell(:)) ) * cos(deg2rad*fld%geom%latCell(:))
  write(*,*) " BJJ TEST function(2), f = ",(cos(deg2rad*10.0_kind_real)+sin(deg2rad*10.0_kind_real))*cos(deg2rad*80.0_kind_real)

  !------------------------------------------------------------
  ! BJJ NICAS-H begin
  !------------------------------------------------------------

  write(*,*) "=== call interp_tl ===="
  locs%nloc = 2
  allocate( locs%lat(locs%nlocs) )
  allocate( locs%lon(locs%nlocs) )

  locs%lon(1) = fld%geom%lonCell(1)
  locs%lat(1) = fld%geom%latCell(1) !26.5650501506626 !30.0
  locs%lon(2) = 10.0_kind_real
  locs%lat(2) = 80.0_kind_real

  
  Nc = fld%geom%nCell
  ! NCAT = 5 !<=== NO GOOD !!!!! GET FROM GEOM
  if (.not.allocated(fld_src)) allocate(fld_src(Nc)) ! <--- Hack job, need to replace with pointers? ...
     do ii=1, nfunc
        fld_src = x_func(:,ii)
        call apply_linop(geov_hinterp_op, fld_src, fld_dst)
        write(*,*) " interpolated valued for level ",ii," is ", fld_dst
     end do
  end do
  deallocate(fld_src)

  !call interp_tl(fld, locs, gom)
  write(*,*) "=== exit interp_tl ===="
  
  !-------------------------------------------------------------
  ! BJJ NICAS-H end
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  !--- BJJ DART Interp. test begin
  !--- BJJ Requires additional GEOMETRY information
  !------------------------------------------------------------
  !> Grid dimensions
  call ncerr(string1, nf90_inq_dimid        (ncid,'nCells',dimid))
  call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=fld%geom%nCells))
  call ncerr(string1, nf90_inq_dimid        (ncid,'nEdges',dimid))
  call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=fld%geom%nEdges))
  call ncerr(string1, nf90_inq_dimid        (ncid,'nVertices',dimid))
  call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=fld%geom%nVertices))
  call ncerr(string1, nf90_inq_dimid        (ncid,'nVertLevels',dimid))
  call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=fld%geom%nVertLevels))
  call ncerr(string1, nf90_inq_dimid        (ncid,'nVertLevelsP1',dimid))
  call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=fld%geom%nVertLevelsP1))
  call ncerr(string1, nf90_inq_dimid        (ncid,'nSoilLevels',dimid))
  call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=fld%geom%nSoilLevels))
  call ncerr(string1, nf90_inq_dimid        (ncid,'vertexDegree',dimid))
  call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=fld%geom%vertexDegree))
  call ncerr(string1, nf90_inq_dimid        (ncid,'maxEdges',dimid))
  call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=fld%geom%maxEdges))
  allocate(nEdgesOnCell(fld%geom%nCells))
  allocate(xVertex(fld%geom%nVertices))
  allocate(yVertex(fld%geom%nVertices))
  allocate(zVertex(fld%geom%nVertices))
  allocate(edgesOnCell(fld%geom%maxEdges,fld%geom%nCells))
  allocate(cellsOnEdge(2,fld%geom%nEdges))
  allocate(verticesOnCell(fld%geom%maxEdges,fld%geom%nCells))
  call ncerr(string1, nf90_inq_varid(ncid,'nEdgesOnCell',varid))
  call ncerr(string1, nf90_get_var  (ncid,varid,nEdgesOnCell))
  call ncerr(string1, nf90_inq_varid(ncid,'xVertex',varid))
  call ncerr(string1, nf90_get_var  (ncid,varid,xVertex))
  call ncerr(string1, nf90_inq_varid(ncid,'yVertex',varid))
  call ncerr(string1, nf90_get_var  (ncid,varid,yVertex))
  call ncerr(string1, nf90_inq_varid(ncid,'zVertex',varid))
  call ncerr(string1, nf90_get_var  (ncid,varid,zVertex))
  call ncerr(string1, nf90_inq_varid(ncid,'edgesOnCell',varid))
  call ncerr(string1, nf90_get_var  (ncid,varid,edgesOnCell,(/1,1/),(/fld%geom%maxEdges,fld%geom%nCells/)))
  call ncerr(string1, nf90_inq_varid(ncid,'cellsOnEdge',varid))
  call ncerr(string1, nf90_get_var  (ncid,varid,cellsOnEdge,(/1,1/),(/2,fld%geom%nEdges/)))
  call ncerr(string1, nf90_inq_varid(ncid,'verticesOnCell',varid))
  call ncerr(string1, nf90_get_var  (ncid,varid,verticesOnCell,(/1,1/),(/fld%geom%maxEdges,fld%geom%nCells/)))

  allocate( x_dart(nfunc*fld%geom%nCells) )
  allocate( dval(nfunc) )
  allocate( ival(1) )

  do iobs=1, locs%nloc
     locs_single%lon = locs%lon(iobs)
     locs_single%lat = locs%lat(iobs)
     !locs_single%vloc = locs%xyz(3,iobs)
     do i=1, nfunc
        x_dart(:) = reshape(x_func(:,:), (/nfunc*fld%geom%nCells/))   ! fld(:,:,:) ~ [nCells,nVertLevels,nf]
        write(*,*) "BJJ------- call compute_scalar_with_barycentric"
        call compute_scalar_with_barycentric(fld%geom%lonCell, fld%geom%latCell, x_dart, locs_single, 1, ival, dval, ier, weights)
        write(*,*) "weights=",weights
        do k=1, nfunc
           write(*,*) " DART interp. valued for function ",k," is ", dval(k)
        end do
     end do
  end do

  deallocate (x_dart)
  deallocate (dval)
  deallocate (ival)
  deallocate(nEdgesOnCell)
  deallocate(xVertex)
  deallocate(yVertex)
  deallocate(zVertex)
  deallocate(edgesOnCell)
  deallocate(cellsOnEdge)
  deallocate(verticesOnCell)
  !--- BJJ DART Interp. test end

  deallocate(nlevels)

  return

end subroutine mpas_test_interp
#endif

! ------------------------------------------------------------------------------

subroutine mpas_nicas_interph_weight(lat1d, lon1d, locs, geov, geov_hinterp_op)

implicit none

    use type_linop
    use tools_interp, only: compute_interp
    use type_randgen, only: rng,initialize_sampling,create_randgen !randgentype
    use module_namelist, only: namtype

    type(ufo_locs), intent(in)       :: locs
    type(ufo_geoval), intent(inout)  :: geov
    real(kind=kind_real), intent(in) :: lon1d(:), lat1d(:)
    type(linoptype), intent(out) :: geov_hinterp_op 
    real(kind=kind_real), allocatable :: lon(:), lat(:), mask(:), masko(:), lono(:), lato(:)
    integer :: Nc, No

    Nc = size(lat1d)
    No = locs%nloc

    if (No>0) then
       print *,'nobs=',No, geov%nobs, Nc
       allocate(lon(Nc), lat(Nc), mask(Nc))
       allocate(masko(No), lono(No), lato(No)) 

       masko = .true. ! Figured out what's the use for masko????
       mask  = .true.
       if (.not.(geov_hinterp_initialized)) then
          print *,'INITIALIZE INTERP'
          rng  = create_randgen(nam)
          lono = deg2rad*locs%xyz(1,:)
          lato = deg2rad*locs%xyz(2,:)
          lon  = deg2rad*lon1d
          lat  = deg2rad*lat1d
          call compute_interp(rng, Nc, lon,  lat,  mask, No, lono, lato, masko, geov_hinterp_op)
          geov_hinterp_initialized = .true.
          write(*,*) " NICAS Weights    =", geov_hinterp_op%S(:)  
          write(*,*) " NICAS ColumnIndex=", geov_hinterp_op%col(:)
       end if
    end if

    deallocate(lon, lat, mask, masko, lono, lato)

end subroutine mpas_nicas_interph_weight

! ------------------------------------------------------------------------------

subroutine nicas_interph(fld, locs, geov, op_type)

    use type_linop
    use tools_interp, only: compute_interp
    use type_randgen, only: rng,initialize_sampling,create_randgen !randgentype
    use module_namelist, only: namtype    
    !use tools_const, only: deg2rad

    type(mpas_field), intent(in)    :: fld
    type(ufo_locs), intent(in)      :: locs
    type(ufo_geoval), intent(inout) :: geov
    character(2), intent(in)        :: op_type !('TL' or 'AD')

    integer :: Nc, No
    real(kind=kind_real), allocatable :: fld_dst(:)
    type(namtype) :: nam !< Namelist variables
    integer :: ii
    logical :: geov_hinterp_initialized = .false.
    type(linoptype) :: geov_hinterp_op

    Nc = fld%geom%nCells
    No = locs%nloc

    if (No>0) then

       allocate(fld_dst(No))

       if (.not.(geov_hinterp_initialized)) then
          call mpas_nicas_interph_weight(fld%geom%latCell, fld%geom%lonCell, locs, geov, geov_hinterp_op)
       end if       

       print *,'&&&&&&&&&&&&&&&&&&& APPLY NICAS_INTERPH TL OPERATOR ',No
       !
       ! Iterate over all fields in pool_b, adding them to fields of the same
       ! name in pool_a
       !
       if (.not.(geov_hinterp_initialized)) then
       call mpas_pool_begin_iteration(fld % subFields)

       do while ( mpas_pool_get_next_member(fld % subFields, poolItr) )

            ! Pools may in general contain dimensions, namelist options, fields, or other pools,
            ! so we select only those members of the pool that are fields
            if (poolItr % memberType == MPAS_POOL_FIELD) then
            ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
            if (poolItr % dataType == MPAS_POOL_REAL) then

               ! Depending on the dimensionality of the field, we need to set pointers of
               ! the correct type
               if (poolItr % nDims == 1) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r1d_ptr_a)
                  select case (op_type)
                  case ('TL') 
                     call apply_linop(geov_hinterp_op, r1d_ptr_a(:), fld_dst)
                  case ('AD')
                     !call apply_linop_ad(hinterp_op,fld_dst,r1d_ptr_a(:))
                  end select
                  
               else if (poolItr % nDims == 2) then
                  call mpas_pool_get_array(pool_a, trim(poolItr % memberName), r2d_ptr_a)
                  do ii=1,fld%geom%nVertLevels
                     select case (op_type)
                     case ('TL') 
                        call apply_linop(geov_hinterp_op, r2d_ptr_a(:,ii), fld_dst)
                     case ('AD')
                        !call apply_linop_ad(hinterp_op,fld_dst,r2d_ptr_a(:,ii))
                     end select
                     write(*,*) " interpolated valued for level ",ii," is ", fld_dst
                  end do

               else if (poolItr % nDims == 3) then
                  write(*,*)'Not implemented yet'
                  !call abort
               end if

            end if
            end if

       end do

       deallocate(fld_dst(No))

    end if

end subroutine nicas_interph


!===================================================================================================
! IMPORT SECTION FORM DART
!====================================================================================================
! BJJ DART Interp. Routines begin

subroutine compute_scalar_with_barycentric(lonCell, latCell, nfunc, x, loc, n, ival, dval, ier, weights)
               !BJJ Special config.
               !x= column vector for given (single) state variable.
               !n= 1         Always.
               !ival(1) = 1  Always.

real(kind=kind_real),            intent(in)  :: x(:), lonCell(:), latCell(:)
type(location_type), intent(in)  :: loc
integer,             intent(in)  :: n, nfunc
integer,             intent(in)  :: ival(:)
real(kind=kind_real),            intent(out) :: dval(:)
integer,             intent(out) :: ier

real(kind=kind_real) :: fract(3), lowval(3), uppval(3), fdata(3), weights(3)
integer  :: lower(3), upper(3), c(3), nvert, index1, k, i, nc, nCell

dval = -9999.99 !MISSING_R8

call find_triangle_vert_indices (lonCell, latCell, x, loc, nc, c, lower, upper, fract, weights, ier)
if(ier /= 0) return

dval = 0.0_kind_real
!do k=1, n
   ! get the starting index in the state vector
   index1 = 1 !BJJ progvar(ival(k))%index1
   nCell = size(lonCell)

   ! go around triangle and interpolate in the vertical
   ! t1, t2, t3 are the xyz of the cell centers
   ! c(3) are the cell ids
   !   do i = 1, nc
   !      lowval(i) = x(index1 + (c(i)-1) * nvert + lower(i)-1)
   !      uppval(i) = x(index1 + (c(i)-1) * nvert + upper(i)-1)
   !      fdata(i) = lowval(i)*(1.0_kind_real - fract(i)) + uppval(i)*fract(i)
   !!      if((debug > 9) .and. do_output()) &
   !      print '(A,I2,A,I2,5f12.5)','compute_scalar_with_barycentric: nv=',k,' ic =',i, &
!                                  lowval(i),uppval(i),fdata(i),fract(i),weights(i)
!
!   enddo

   ! now have vertically interpolated values at cell centers.
   ! use weights to compute value at interp point.
!   dval(k) = sum(weights(1:nc) * fdata(1:nc))
!print *, 'k, dval(k) = ', k, dval(k)
!enddo
do k=1, nfunc
   do i = 1, nc
!    x(:) = reshape(fld%fld(:,:,i), (/fld%geom%nVertLevels*fld%geom%nCells/))   ! fld(:,:,:) ~ [nCells,nVertLevels,nf]
     !wrong: dval(k) = dval(k) + weights(i) * x(index1 + (c(i)-1)*nvert + k-1 )     !~ x: vert. inner-column
     !wrong: if (i.eq.1) write(*,*) " testttttttttt x(index1 + (c(1)-1)*nvert + k-1) = ",x(index1 + (c(i)-1)*nvert + k-1 )
     dval(k) = dval(k) + weights(i) * x(index1 + c(i)-1 + (k-1)*nCell )
     !if (i.eq.1) write(*,*) " testttttttttt x(index1 + c(1)-1 + (k-1)*nCell ) = ",x(index1 + c(i)-1 + (k-1)*nCell )
   enddo
enddo

ier = 0

end subroutine compute_scalar_with_barycentric

!call find_triangle_vert_indices (x, loc, nc, c, lower, upper, fract, weights, ier)
subroutine find_triangle_vert_indices (lonCell, latCell, x, loc, nc, c, lower, upper, fract, weights, ier)

implicit none

real(kind=kind_real),            intent(in)  :: x(:), lonCell(:), latCell(:)
type(location_type), intent(in)  :: loc
integer,             intent(out) :: nc
integer,             intent(out) :: c(:)
integer,             intent(out) :: lower(:), upper(:)
real(kind=kind_real),            intent(out) :: fract(:)
real(kind=kind_real),            intent(out) :: weights(:)
integer,             intent(out) :: ier

! compute the values at the correct vertical level for each
! of the 3 cell centers defining a triangle that encloses the
! the interpolation point, then interpolate once in the horizontal
! using barycentric weights to get the value at the interpolation point.

integer, parameter :: listsize = 30
integer  :: nedges, i, neighborcells(maxEdges), edgeid
real(kind=kind_real) :: xdata(listsize), ydata(listsize), zdata(listsize)
real(kind=kind_real) :: t1(3), t2(3), t3(3), r(3)
integer  :: cellid, verts(listsize), closest_vert
real(kind=kind_real) :: lat, lon, vert, llv(3)
integer  :: verttype, vindex, v, vp1
logical  :: inside, foundit


! initialization
      c = MISSING_R8
  lower = MISSING_R8
  upper = MISSING_R8
  fract = 0.0_kind_real
weights = 0.0_kind_real
    ier = 0
     nc = 1

! unpack the location into local vars
!BJJ llv = get_location(loc)
lon  = loc%lon  !BJJ llv(1)
lat  = loc%lat  !BJJ llv(2)
vert = loc%vloc !BJJ llv(3)
!BJJ verttype = nint(query_location(loc))

cellid = find_closest_cell_center(lonCell, latCell, lat, lon)
!if ((xyzdebug > 5) .and. do_output()) &
   print *, 'closest cell center for lon/lat: ', lon, lat, cellid
if (cellid < 1) then
   !if(xyzdebug > 0)
   print *, 'closest cell center for lon/lat: ', lon, lat, cellid
   ier = 11
   return
endif

c(1) = cellid

!if (on_boundary(cellid)) then
!   ier = 12
!   return
!endif
!
!if (.not. inside_cell(cellid, lat, lon)) then
!   ier = 13
!   return
!endif

! closest vertex to given point.
closest_vert = closest_vertex_ll(cellid, lat, lon)
!if ((xyzdebug > 5) .and. do_output()) &
   print *, 'closest vertex for lon/lat: ', lon, lat, closest_vert

! collect the neighboring cell ids and vertex numbers
! this 2-step process avoids us having to read in the
! cellsOnCells() array which i think we only need here.
! if it comes up in more places, we can give up the space
! and read it in and then this is a direct lookup.
! also note which index is the closest vert and later on
! we can start the triangle search there.
vindex = 1
nedges = nEdgesOnCell(cellid)
do i=1, nedges
   edgeid = edgesOnCell(i, cellid)
   if (cellsOnEdge(1, edgeid) /= cellid) then
      neighborcells(i) = cellsOnEdge(1, edgeid)
   else
      neighborcells(i) = cellsOnEdge(2, edgeid)
   endif
   verts(i) = verticesOnCell(i, cellid)
   if (verts(i) == closest_vert) vindex = i
   call latlon_to_xyz(latCell(neighborcells(i)), lonCell(neighborcells(i)), &
      xdata(i), ydata(i), zdata(i))
enddo


! get the cartesian coordinates in the cell plane for the closest center
call latlon_to_xyz(latCell(cellid), lonCell(cellid), t1(1), t1(2), t1(3))

! and the observation point
call latlon_to_xyz(lat, lon, r(1), r(2), r(3))

if (all(abs(t1-r) < roundoff)) then   ! Located at a grid point (counting roundoff errors)

   ! need vert index for the vertical level
   nc = 1

   ! This is a grid point - horiz interpolation NOT needed
   weights(1) = 1.0_kind_real

else                       ! an arbitrary point

! find the cell-center-tri that encloses the obs point
! figure out which way vertices go around cell?
foundit = .false.
findtri: do i=vindex, vindex+nedges
   v = mod(i-1, nedges) + 1
   vp1 = mod(i, nedges) + 1
   t2(1) = xdata(v)
   t2(2) = ydata(v)
   t2(3) = zdata(v)
   t3(1) = xdata(vp1)
   t3(2) = ydata(vp1)
   t3(3) = zdata(vp1)
   call inside_triangle(t1, t2, t3, r, lat, lon, inside, weights)
   if (inside) then
      ! weights are the barycentric weights for the point r
      ! in the triangle formed by t1, t2, t3.
      ! v and vp1 are vert indices which are same indices
      ! for cell centers
! FIXME: i want to remove this code.  does it affect the answers?
      if(any(weights == 1.0_kind_real)) then
         nc = 1
      else
!end FIXME section
         nc = 3
         c(2) = neighborcells(v)
         c(3) = neighborcells(vp1)
! FIXME: i want to remove this code.  does it affect the answers?
      endif
!end FIXME section
      foundit = .true.
      exit findtri
   endif
enddo findtri
if (.not. foundit) then
   ier = 14     ! 11
   return
endif

endif     ! horizontal index search is done now.

!! need vert index for the vertical level
!call find_vert_level(x, loc, nc, c, .true., lower, upper, fract, ier)
!if(ier /= 0) return
!
!if ((debug > 9) .and. do_output()) then
!   write(string3,*) 'ier = ',ier, ' triangle = ',c(1:nc), ' vert_index = ',lower(1:nc)+fract(1:nc)
!   call error_handler(E_MSG, 'find_triangle_vert_indices', string3, source, revision, revdate)
!endif

!if ((debug > 8) .and. do_output()) then
   print *, 'nc = ', nc
   print *, 'c = ', c
   print *, 'lonCell(c) = ', lonCell(c(1:nc)) !, lonCell(c(2)), lonCell(c(3))
   print *, 'latCell(c) = ', latCell(c(1:nc)) !, latCell(c(2)), latCell(c(3))
   print *, 'lower = ', lower(1:nc)
   print *, 'upper = ', upper(1:nc)
!   print *, 'fract = ', fract(1:nc)
   print *, 'weights = ', weights
!   print *, 'ier = ', ier
!endif

end subroutine find_triangle_vert_indices

!-------------------------------------------------------------------------------------------------------------------------

subroutine latlon_to_xyz(lat, lon, x, y, z)

! Given a lat, lon in degrees, return the cartesian x,y,z coordinate
! on the surface of a specified radius relative to the origin
! at the center of the earth.  (this radius matches the one
! used at MPAS grid generation time and must agree in order
! to be consistent with the cartisian coordinate arrays in
! the MPAS data files.)

real(kind=kind_real), intent(in)  :: lat, lon
real(kind=kind_real), intent(out) :: x, y, z

real(kind=kind_real) :: rlat, rlon

rlat = lat * deg2rad
rlon = lon * deg2rad

x = radius * cos(rlon) * cos(rlat)
y = radius * sin(rlon) * cos(rlat)
z = radius * sin(rlat)

end subroutine latlon_to_xyz


!-------------------------------------------------------------------------------------------------------------------------


function closest_vertex_ll(cellid, lat, lon)

! Return the vertex id of the closest one to the given point
! this version uses lat/lon.  see closest_vertex_xyz for the
! cartesian version.

integer,  intent(in)  :: cellid
real(kind=kind_real), intent(in)  :: lat, lon
integer               :: closest_vertex_ll

real(kind=kind_real) :: px, py, pz

! use the same radius as MPAS for computing this
call latlon_to_xyz(lat, lon, px, py, pz)

closest_vertex_ll = closest_vertex_xyz(cellid, px, py, pz)
!if ((closest_vertex_ll < 0)  .and. &
!    (debug > 8) .and. do_output()) &
   print *, 'cannot find nearest vertex to lon, lat: ', lon, lat

end function closest_vertex_ll


!-------------------------------------------------------------------------------------------------------------------------


function closest_vertex_xyz(cellid, px, py, pz)

! Return the vertex id of the closest one to the given point
! see closest_vertex_ll for the lat/lon version (which calls this)

integer,  intent(in)  :: cellid
real(kind=kind_real), intent(in)  :: px, py, pz
integer               :: closest_vertex_xyz

integer :: nverts, i, vertexid
real(kind=kind_real) :: distsq, closest_dist, dx, dy, dz

! nedges and nverts is same in a closed figure
nverts = nEdgesOnCell(cellid)

closest_dist = 1.0e38_kind_real   ! something really big; these are meters not radians
closest_vertex_xyz = -1

do i=1, nverts
   vertexid = verticesOnCell(i, cellid)
   dx = xVertex(vertexid) - px
   dy = yVertex(vertexid) - py
   dz = zVertex(vertexid) - pz
   distsq = (dx * dx) + (dy * dy) + (dz * dz)
   if (distsq < closest_dist) then
      closest_dist = distsq
      closest_vertex_xyz = vertexid
   endif
enddo

end function closest_vertex_xyz

!-------------------------------------------------------------------------------------------------------------------------

function get_location(loc)

! Given a location type (in radians), 
! return the longitude, latitude (in degrees) and vertical value

type(location_type), intent(in) :: loc
real(kind=kind_real), dimension(3) :: get_location

get_location(1) = loc%lon * RAD2DEG
get_location(2) = loc%lat * RAD2DEG
get_location(3) = loc%vloc

end function get_location

!-------------------------------------------------------------------------------------------------------------------------

function find_closest_cell_center(lonCell, latCell, lat, lon)

! Determine the cell index for the closest center to the given point
! 2D calculation only.

real(kind=kind_real), intent(in)  :: lat, lon, lonCell(:), latCell(:)
integer               :: find_closest_cell_center

type(xyz_location_type) :: pointloc
integer :: closest_cell, rc
integer :: i
type(xyz_location_type), allocatable :: cell_locs(:)
type(xyz_get_close_type)             :: cc_gc

allocate(cell_locs(nCells))
do i=1, nCells
   cell_locs(i) = xyz_set_location(lonCell(i), latCell(i), 0.0_kind_real, radius)
enddo
call xyz_get_close_maxdist_init(cc_gc, 1.0_kind_real)
call xyz_get_close_obs_init(cc_gc, nCells, cell_locs)

pointloc = xyz_set_location(lon, lat, 0.0_kind_real, radius)

call xyz_find_nearest(cc_gc, pointloc, cell_locs, closest_cell, rc)

! decide what to do if we don't find anything.
if (rc /= 0 .or. closest_cell < 0) then
!   if ((debug > 8) .and. do_output()) &
       print *, 'cannot find nearest cell to lon, lat: ', lon, lat
   find_closest_cell_center = -1
   return
endif

! this is the cell index for the closest center
find_closest_cell_center = closest_cell

end function find_closest_cell_center

!-------------------------------------------------------------------------------------------------------------------------

function set_location_lonlat(lon, lat, height, radius)

! location semi-independent interface routine
! given a lon, lat, and radius, compute X,Y,Z and set location

real(kind_real), intent(in) :: lon, lat, height, radius
type (xyz_location_type) :: set_location_lonlat

real(kind_real) :: x, y, z
real(kind_real) :: rlat, rlon, rr

!if ( .not. module_initialized ) call initialize_module

rlat = lat * deg2rad
rlon = lon * deg2rad

rr = radius + height

x = rr * cos(rlon) * cos(rlat)
y = rr * sin(rlon) * cos(rlat)
z = rr * sin(rlat)

set_location_lonlat%x = x
set_location_lonlat%y = y
set_location_lonlat%z = z

end function set_location_lonlat


!-------------------------------------------------------------------------------------------------------------------------

subroutine xyz_find_nearest(gc, base_loc, loc_list, nearest, rc)
 type(xyz_get_close_type), intent(in), target  :: gc
 type(xyz_location_type),  intent(in)  :: base_loc
 type(xyz_location_type),  intent(in)  :: loc_list(:)
 integer,              intent(out) :: nearest
 integer,              intent(out) :: rc

! find the nearest point in the get close list to the specified point

   call find_nearest_boxes(gc, base_loc, loc_list, nearest, rc)

end subroutine xyz_find_nearest

!-------------------------------------------------------------------------------------------------------------------------

subroutine xyz_get_close_obs_init(gc, num, obs)

! Initializes part of get_close accelerator that depends on the particular obs

type(xyz_get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
type(xyz_location_type),  intent(in)    :: obs(num)

   call get_close_init_boxes(gc, num, obs)

end subroutine xyz_get_close_obs_init


!-------------------------------------------------------------------------------------------------------------------------

subroutine find_nearest_boxes(gc, base_loc, loc_list, nearest, rc)
 type(xyz_get_close_type), intent(in), target  :: gc
 type(xyz_location_type),  intent(in)  :: base_loc
 type(xyz_location_type),  intent(in)  :: loc_list(:)
 integer,              intent(out) :: nearest
 integer,              intent(out) :: rc

! find the nearest point in the get close list to the specified point

integer :: x_box, y_box, z_box, i, j, k, l
integer :: start_x, end_x, start_y, end_y, start_z, end_z
integer ::  n_in_box, st, t_ind, ghost
real(kind_real) :: this_dist, dist

! First, set the intent out arguments to a missing value
nearest = -99
rc = -1
dist = 1e38_kind_real                ! something big and positive.

! the list of locations in the obs() argument must be the same
! as the list of locations passed into get_close_obs_init(), so
! gc%num and size(obs) better be the same.   if the list changes,
! you have to destroy the old gc and init a new one.
if (size(loc_list) /= gc%num) then
   !write(errstring,*)'obs() array must match one passed to get_close_obs_init()'
   !call error_handler(E_ERR, 'get_close_obs', errstring, source, revision, revdate)
endif

! If num == 0, no point in going any further. 
if (gc%num == 0) return

!--------------------------------------------------------------

! Begin by figuring out which box the base loc is in
x_box = floor((base_loc%x - gc%box%bot_x) / gc%box%x_width) + 1
if(x_box > nx) x_box = nx
if(x_box < 1)  x_box = 1
y_box = floor((base_loc%y - gc%box%bot_y) / gc%box%y_width) + 1
if(y_box > ny) y_box = ny
if(y_box < 1)  y_box = 1
z_box = floor((base_loc%z - gc%box%bot_z) / gc%box%z_width) + 1
if(z_box > nz) z_box = nz
if(z_box < 1)  z_box = 1

!print *, 'base_loc box ', x_box, y_box, z_box
!print *, 'nx, ny, nz = ', nx, ny, nz

! If it is not in any box, then it is more than the maxdist away from everybody
if(x_box > nx .or. x_box < 1) return
if(y_box > ny .or. y_box < 1) return
if(z_box > nz .or. z_box < 1) return

!print *, 'good box'


! First, search all points in this box.

! Box to search is x_box,y_box,z_box
n_in_box = gc%box%count(x_box,y_box,z_box)
st = gc%box%start(x_box,y_box,z_box)

! find the closest one in this box
do l = 1, n_in_box

   t_ind = gc%box%obs_box(st - 1 + l)
!print *, 'l, t_ind = ', l, t_ind

   this_dist = xyz_get_dist(base_loc, loc_list(t_ind))
!print *, 'this_dist = ', this_dist
   ! If this obs' distance is less than current nearest, it's new nearest
   if(this_dist <= dist) then
      nearest = t_ind
      dist = this_dist
      if (rc < 0) rc = 0
   endif
end do

! if box small enough that no points match, expand search
ghost = 0

10 continue
if (nearest < 0 .or. ghost == 0)  then
   ghost = ghost + 1

   start_x = x_box - ghost
   if (start_x < 1) start_x = 1
   end_x = x_box + ghost
   if (end_x > nx) end_x = nx

   start_y = y_box - ghost
   if (start_y < 1) start_y = 1
   end_y = y_box + ghost
   if (end_y > ny) end_y = ny

   start_z = z_box - ghost
   if (start_z < 1) start_z = 1
   end_z = z_box + ghost
   if (end_z > nz) end_z = nz

   !print *, 'looping from '
   !print *, 'x: ', start_x, end_x
   !print *, 'y: ', start_y, end_y
   !print *, 'z: ', start_z, end_z

   ! Next, loop through each box that is close to this box
   do i = start_x, end_x
      do j = start_y, end_y
         do k = start_z, end_z

            ! Box to search is i,j,k
            n_in_box = gc%box%count(i, j, k)
            st = gc%box%start(i,j,k)


            ! Loop to check how close all obs in the box are; add those that are close
            do l = 1, n_in_box

               t_ind = gc%box%obs_box(st - 1 + l)
  ! print *, 'l, t_ind = ', l, t_ind

               this_dist = xyz_get_dist(base_loc, loc_list(t_ind))
  ! print *, 'this_dist = ', this_dist
               ! If this obs' distance is less than current nearest, it's new nearest
               if(this_dist <= dist) then
                  nearest = t_ind
                  dist = this_dist
                  if (rc < 0) rc = 0
               endif
            end do
         end do
      end do
   end do

   if (nearest < 0) then
      ! if we have searched the entire space, punt.
      if (start_x == 1 .and. end_x == nx .and. &
          start_y == 1 .and. end_y == ny .and. &
          start_z == 1 .and. end_z == nz) return

      ! repeat search with larger radius of boxes
      goto 10
   endif
endif

end subroutine find_nearest_boxes

!-------------------------------------------------------------------------------------------------------------------------

function xyz_get_dist(loc1, loc2, type1, kind2)

! returns the distance between 2 locations 

! the 3rd name is a specific type, the 4th is a generic kind.
! These types/kinds are part of
! the interface in case user-code wants to do a more sophisticated distance
! calculation based on the base type or target kind.  In the usual case this
! code still doesn't use the types, but there's an undocumented feature that
! allows you to maintain the original vertical normalization even when
! changing the cutoff distance in the horizontal.  For that to work we
! do need to know the type, and we use the type of loc1 to control it.
! 

type(xyz_location_type), intent(in) :: loc1, loc2
integer, optional,       intent(in) :: type1, kind2
real(kind=kind_real)                            :: xyz_get_dist

real(kind=kind_real) :: x_dif, y_dif, z_dif

!if ( .not. module_initialized ) call initialize_module

x_dif = loc1%x - loc2%x
y_dif = loc1%y - loc2%y
z_dif = loc1%z - loc2%z

xyz_get_dist = sqrt(x_dif * x_dif + y_dif * y_dif + z_dif * z_dif)

end function xyz_get_dist

!-------------------------------------------------------------------------------------------------------------------------

subroutine xyz_get_close_maxdist_init(gc, maxdist)

type(xyz_get_close_type), intent(inout) :: gc
real(kind_real),             intent(in)    :: maxdist

character(len=129) :: str1
integer :: i

! set the default value.
gc%maxdist = maxdist
!print *, 'setting maxdist to ', maxdist

!if (.not. use_octree) then
if (.not. allocated(gc%box%count)) then
   ! Allocate the storage for the grid dependent boxes
   allocate(gc%box%count(nx,ny,nz), gc%box%start(nx,ny,nz))
   gc%box%count  = -1
   gc%box%start  = -1
endif

end subroutine xyz_get_close_maxdist_init


!-------------------------------------------------------------------------------------------------------------------------

subroutine get_close_init_boxes(gc, num, obs)

! Initializes part of get_close accelerator that depends on the particular obs

type(xyz_get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
type(xyz_location_type),  intent(in)    :: obs(num)

integer :: i, j, k, cum_start, l
integer :: x_box(num), y_box(num), z_box(num)
integer :: tstart(nx, ny, nz)

!if ( .not. module_initialized ) call initialize_module

! Allocate storage for obs number dependent part
allocate(gc%box%obs_box(num))
gc%box%obs_box(:) = -1

! Set the value of num_obs in the structure
gc%num = num

! If num == 0, no point in going any further.
if (num == 0) return

! Determine where the boxes should be for this set of obs and maxdist
call find_box_ranges(gc, obs, num)

! Begin by computing the number of observations in each box in x,y,z
gc%box%count = 0
do i = 1, num

!print *, i, obs(i)%x, obs(i)%y, obs(i)%z
   x_box(i) = floor((obs(i)%x - gc%box%bot_x) / gc%box%x_width) + 1
   if(x_box(i) > nx) x_box(i) = nx
   if(x_box(i) < 1)  x_box(i) = 1

   y_box(i) = floor((obs(i)%y - gc%box%bot_y) / gc%box%y_width) + 1
   if(y_box(i) > ny) y_box(i) = ny
   if(y_box(i) < 1)  y_box(i) = 1

   z_box(i) = floor((obs(i)%z - gc%box%bot_z) / gc%box%z_width) + 1
   if(z_box(i) > nz) z_box(i) = nz
   if(z_box(i) < 1)  z_box(i) = 1

   gc%box%count(x_box(i), y_box(i), z_box(i)) = gc%box%count(x_box(i), y_box(i), z_box(i)) + 1
!print *, 'adding count to box ', x_box(i), y_box(i), z_box(i), &
!                                 gc%box%count(x_box(i), y_box(i), z_box(i))
end do

! Figure out where storage for each boxes members should begin
cum_start = 1
do i = 1, nx
   do j = 1, ny
      do k = 1, nz
         gc%box%start(i, j, k) = cum_start
         cum_start = cum_start + gc%box%count(i, j, k)
      end do
   end do
end do

! Now we know how many are in each box, get a list of which are in each box
tstart = gc%box%start
do i = 1, num
   gc%box%obs_box(tstart(x_box(i), y_box(i), z_box(i))) = i
   tstart(x_box(i), y_box(i), z_box(i)) = tstart(x_box(i), y_box(i), z_box(i)) + 1
end do

do i = 1, nx
   do j = 1, ny
      do k = 1, nz
!if (gc%box%count(i,j,k) > 0) print *, i,j,k, gc%box%count(i,j,k), gc%box%start(i,j,k)
         do l=1, gc%box%count(i,j,k)
!print *, l, gc%box%obs_box(l)
         enddo
      end do
   end do
end do


end subroutine get_close_init_boxes

!-------------------------------------------------------------------------------------------------------------------------

subroutine find_box_ranges(gc, locs, num)

! Finds boundaries for x,y,z boxes.  
! FIXME: ways boxes could be divided:
!  - evenly along each axis
!  - octree-like, divide each axis so roughly half the points are
!     on each side of the dividing plane.
!  - about 100 other schemes

type(xyz_get_close_type), intent(inout) :: gc
integer,              intent(in)    :: num
type(xyz_location_type),  intent(in)    :: locs(num)

logical :: old_out
integer :: i


! FIXME: this space could be very sparse

gc%box%bot_x = minval(locs(:)%x)
gc%box%bot_y = minval(locs(:)%y)
gc%box%bot_z = minval(locs(:)%z)

gc%box%top_x = maxval(locs(:)%x)
gc%box%top_y = maxval(locs(:)%y)
gc%box%top_z = maxval(locs(:)%z)

!gc%box%bot_x = locs(1)%x
!gc%box%bot_y = locs(1)%y
!gc%box%bot_z = locs(1)%z
!  
!gc%box%top_x = locs(1)%x
!gc%box%top_y = locs(1)%y
!gc%box%top_z = locs(1)%z
!
!do i=2, num
!   gc%box%bot_x = min(gc%box%bot_x, locs(i)%x)
!   gc%box%bot_y = min(gc%box%bot_y, locs(i)%y)
!   gc%box%bot_z = min(gc%box%bot_z, locs(i)%z)
!   
!   gc%box%top_x = max(gc%box%top_x, locs(i)%x)
!   gc%box%top_y = max(gc%box%top_y, locs(i)%y)
!   gc%box%top_z = max(gc%box%top_z, locs(i)%z)
!enddo

gc%box%x_width = (gc%box%top_x - gc%box%bot_x) / nx
gc%box%y_width = (gc%box%top_y - gc%box%bot_y) / ny
gc%box%z_width = (gc%box%top_z - gc%box%bot_z) / nz

! FIXME:  compute a sphere of radius maxdist and see how
! many boxes in x, y, z that would include.
gc%box%nboxes_x = aint((gc%maxdist + (gc%box%x_width-1)) / gc%box%x_width)
gc%box%nboxes_y = aint((gc%maxdist + (gc%box%y_width-1)) / gc%box%y_width)
gc%box%nboxes_z = aint((gc%maxdist + (gc%box%z_width-1)) / gc%box%z_width)

!print *, 'min xyz = ', gc%box%bot_x, gc%box%bot_y, gc%box%bot_z
!print *, 'max xyz = ', gc%box%top_x, gc%box%top_y, gc%box%top_z
!print *, 'wid xyz = ', gc%box%x_width, gc%box%y_width, gc%box%z_width
!print *, 'nbx xyz = ', nx, ny, nz
!print *, 'nbx xyz = ', gc%box%nboxes_x, gc%box%nboxes_y, gc%box%nboxes_z

end subroutine find_box_ranges

!-------------------------------------------------------------------------------------------------------------------------

subroutine inside_triangle(t1, t2, t3, r, lat, lon, inside, weights)

! given 3 corners of a triangle and an xyz point, compute whether
! the point is inside the triangle.  this assumes r is coplanar
! with the triangle - the caller must have done the lat/lon to
! xyz conversion with a constant radius and then this will be
! true (enough).  sets inside to true/false, and returns the
! weights if true.  weights are set to 0 if false.

real(kind=kind_real), intent(in)  :: t1(3), t2(3), t3(3)
real(kind=kind_real), intent(in)  :: r(3), lat, lon
logical,  intent(out) :: inside
real(kind=kind_real), intent(out) :: weights(3)

! check for degenerate cases first - is the test point located
! directly on one of the vertices?  (this case may be common
! if we're computing on grid point locations.)
if (all(abs(r - t1) < roundoff)) then
   inside = .true.
   weights = (/ 1.0_kind_real, 0.0_kind_real, 0.0_kind_real /)
   return
else if (all(abs(r - t2) < roundoff)) then
   inside = .true.
   weights = (/ 0.0_kind_real, 1.0_kind_real, 0.0_kind_real /)
   return
else if (all(abs(r - t3) < roundoff)) then
   inside = .true.
   weights = (/ 0.0_kind_real, 0.0_kind_real, 1.0_kind_real /)
   return
endif

! not a vertex. compute the weights.  if any are
! negative, the point is outside.  since these are
! real valued computations define a lower bound for
! numerical roundoff error and be sure that the
! weights are not just *slightly* negative.
call get_3d_weights(r, t1, t2, t3, lat, lon, weights)

if (any(weights < -roundoff)) then
   inside = .false.
   weights = 0.0_kind_real
   return
endif

! truncate barely negative values to 0
inside = .true.
where (weights < 0.0_kind_real) weights = 0.0_kind_real
return

end subroutine inside_triangle

!-------------------------------------------------------------------------------------------------------------------------

subroutine get_3d_weights(p, v1, v2, v3, lat, lon, weights)

! Given a point p (x,y,z) inside a triangle, and the (x,y,z)
! coordinates of the triangle corner points (v1, v2, v3),
! find the weights for a barycentric interpolation.  this
! computation only needs two of the three coordinates, so figure
! out which quadrant of the sphere the triangle is in and pick
! the 2 axes which are the least planar:
!  (x,y) near the poles,
!  (y,z) near 0 and 180 longitudes near the equator,
!  (x,z) near 90 and 270 longitude near the equator.
! (lat/lon are the coords of p. we could compute them here
! but since in all cases we already have them, pass them
! down for efficiency)

real(kind_real), intent(in)  :: p(3)
real(kind_real), intent(in)  :: v1(3), v2(3), v3(3)
real(kind_real), intent(in)  :: lat, lon
real(kind_real), intent(out) :: weights(3)

real(kind_real) :: cxs(3), cys(3)

! above or below 45 in latitude, where -90 < lat < 90:
if (lat >= 45.0_kind_real .or. lat <= -45.0_kind_real) then
   cxs(1) = v1(1)
   cxs(2) = v2(1)
   cxs(3) = v3(1)
   cys(1) = v1(2)
   cys(2) = v2(2)
   cys(3) = v3(2)
   call get_barycentric_weights(p(1), p(2), cxs, cys, weights)
   return
endif

! nearest 0 or 180 in longitude, where 0 < lon < 360:
if ( lon <= 45.0_kind_real .or. lon >= 315.0_kind_real .or. &
    (lon >= 135.0_kind_real .and. lon <= 225.0_kind_real)) then
   cxs(1) = v1(2)
   cxs(2) = v2(2)
   cxs(3) = v3(2)
   cys(1) = v1(3)
   cys(2) = v2(3)
   cys(3) = v3(3)
   call get_barycentric_weights(p(2), p(3), cxs, cys, weights)
   return
endif

! last option, nearest 90 or 270 in lon:
cxs(1) = v1(1)
cxs(2) = v2(1)
cxs(3) = v3(1)
cys(1) = v1(3)
cys(2) = v2(3)
cys(3) = v3(3)
call get_barycentric_weights(p(1), p(3), cxs, cys, weights)

end subroutine get_3d_weights

!-------------------------------------------------------------------------------------------------------------------------

subroutine get_barycentric_weights(x, y, cxs, cys, weights)

! Computes the barycentric weights for a 2d interpolation point
! (x,y) in a 2d triangle with the given (cxs,cys) corners.

real(kind_real), intent(in)  :: x, y, cxs(3), cys(3)
real(kind_real), intent(out) :: weights(3)

real(kind_real) :: denom

! Get denominator
denom = (cys(2) - cys(3)) * (cxs(1) - cxs(3)) + &
   (cxs(3) - cxs(2)) * (cys(1) - cys(3))

weights(1) = ((cys(2) - cys(3)) * (x - cxs(3)) + &
   (cxs(3) - cxs(2)) * (y - cys(3))) / denom

weights(2) = ((cys(3) - cys(1)) * (x - cxs(3)) + &
   (cxs(1) - cxs(3)) * (y - cys(3))) / denom

weights(3) = 1.0_kind_real - weights(1) - weights(2)

! FIXME: i want to remove this code.  does it affect the answers?
if (any(abs(weights) < roundoff)) then
   !print *, 'get_barycentric_weights due to roundoff errors: ', weights
   where (abs(weights) < roundoff) weights = 0.0_kind_real
   where (abs(1.0_kind_real - abs(weights)) < roundoff) weights = 1.0_kind_real
endif
if(abs(sum(weights)-1.0_kind_real) > roundoff) &
   print *, 'fail in get_barycentric_weights: sum(weights) = ',sum(weights)
!end FIXME section

end subroutine get_barycentric_weights

!
! BJJ DART Interp. Routines end
!=============================================================================================
! END OF DART SECTION
!=============================================================================================

end module mpas_fields_interp_mod
