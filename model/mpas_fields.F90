
!> Handle fields for the MPAS-A model

module mpas_fields

use config_mod
use mpas_geom_mod
use mpas_vars_mod
use kinds
use tools_nc
use netcdf

!MPAS
!use io_input
!use io_output
!use grid_types
use mpas_framework
!use mpas_subdriver
!use mpas_configs
!!use configure
!use cell_search


implicit none
private

public :: mpas_field, &
        & create, delete, zeros, random, copy, &
        & self_add, self_schur, self_sub, self_mul, axpy, &
        & dot_prod, add_incr, diff_incr, change_resol, &
        & read_file, write_file, gpnorm, fldrms, &
        & convert_to_ug, convert_from_ug, dirac
public :: mpas_field_registry

! ------------------------------------------------------------------------------
 integer, parameter :: max_string_length = 128
 integer, parameter :: max_variable_dims = 20
! ------------------------------------------------------------------------------

!> Fortran derived type to hold MPAS-A fields
type :: mpas_field
  type(mpas_geom), pointer :: geom                                 !< Number of unstructured grid cells
  integer :: nf                           !< Number of variables in fld
  integer :: ns                           !< Number of surface fields (x1d [nCells])
  real(kind=kind_real), pointer :: fld(:,:,:)       !< all 3D fields in each column ([nCells,nVertLevels,nf])
  real(kind=kind_real), pointer :: th(:,:)          !< Potential temperature
  real(kind=kind_real), pointer ::  r(:,:)          !< Dry density
  real(kind=kind_real), pointer ::  u(:,:)          !< Zonal wind reconstructed at cell centers
  real(kind=kind_real), pointer ::  v(:,:)          !< Meridional wind reconstructed at cell centers
  real(kind=kind_real), pointer ::  q(:,:)          !< Water vapor mixing ratio at cell centers (Qv)
  real(kind=kind_real), pointer :: qc(:,:)          !< Cloud mixing ratio at cell centers
  real(kind=kind_real), pointer :: qr(:,:)          !< Rain mixing ratio at cell centers
  real(kind=kind_real), pointer :: qi(:,:)          !< Ice mixing ratio at cell centers
  real(kind=kind_real), pointer :: qs(:,:)          !< Snow mixing ratio at cell centers
  real(kind=kind_real), pointer :: qg(:,:)          !< Graupel mixing ratio at cell centers
  real(kind=kind_real), pointer ::  s(:,:)          !< Normal speed at cell edges
  real(kind=kind_real), pointer ::  w(:,:)          !< Vertical velocity - not used for now
  real(kind=kind_real), pointer :: x1d(:,:)         !< multiple 1D surface fields in each column [nCells,nf]
  character(len=20), allocatable :: fldnames(:)      !< Variable identifiers
end type mpas_field

#define LISTED_TYPE mpas_field

!> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_field_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine create(self, geom, vars)
implicit none
type(mpas_field), intent(inout)       :: self
type(mpas_geom),  intent(in), pointer :: geom
type(mpas_vars),  intent(in)          :: vars
integer :: ioff, nl, nC

associate(nl=>geom%nVertLevels, nC=>geom%nCells)
!associate(nl=>self%geom%nVertLevels, nC=>self%geom%nCells)

self%geom => geom
self%nf = vars%nv

! self%fldnames(:) = (/"theta","rho","qv","uReconstructZonal","uReconstructMeridional"/)
allocate(self%fld(nC,nl,self%nf))
self%fld(:,:,:)=0.0_kind_real   ! initialize 3-D variables

ioff=0
self%th => self%fld(:,:,1)      ! potential temperature (theta)
!ioff =ioff + nl
self%s  => null()               ! Not used for now.

if (self%nf>1) then
  self%r  => self%fld(:,:,2)    ! dry density (rho)
  self%q  => self%fld(:,:,3)    ! water vapor mixing ratio (Qv)
  self%u  => self%fld(:,:,4)    ! reconstructed zonal wind
  self%v  => self%fld(:,:,5)    ! reconstructed meridional wind
  !self%r  => self%fld(:,ioff+1:ioff+nl)
  !ioff =ioff + nl

  !if (self%nf == 4) then	! FIXME: In this case, I need to use ioff for nE in self%s.
  !    self%s  => self%fld(:,:,4)
  !    self%u  => null()
  !    self%v  => null()
  !endif
else
  self%r  => null()
  self%q  => null()
  self%u  => null()
  self%v  => null()
endif
!if (ioff/=self%nf*self%nl)  call abor1_ftn ("mpas_fields:create error number of fields")

allocate(self%fldnames(self%nf))
self%fldnames(:)=vars%fldnames(:)

call check(self)

end associate

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)
implicit none
type(mpas_field), intent(inout) :: self

call check(self)

if (associated(self%fld)) deallocate(self%fld)
!if (associated(self%x1d)) deallocate(self%x1d)
!if (associated(self%th)) deallocate(self%th)
!if (associated(self%r)) deallocate(self%r)
!if (associated(self%u)) deallocate(self%u)
!if (associated(self%v)) deallocate(self%v)
!if (associated(self%s)) deallocate(self%s)
!if (associated(self%q)) deallocate(self%q)
!if (associated(self%w)) deallocate(self%w)
!if (associated(self%qc)) deallocate(self%qc)
!if (associated(self%qr)) deallocate(self%qr)
!if (associated(self%qi)) deallocate(self%qi)
!if (associated(self%qs)) deallocate(self%qs)
!if (associated(self%qg)) deallocate(self%qg)
if (allocated(self%fldnames)) deallocate(self%fldnames)

end subroutine delete

! ------------------------------------------------------------------------------

subroutine zeros(self)
implicit none
type(mpas_field), intent(inout) :: self

call check(self)

self%fld(:,:,:) = 0.0_kind_real

end subroutine zeros

! ------------------------------------------------------------------------------

subroutine random(self)
use random_vectors_mod
implicit none
type(mpas_field), intent(inout) :: self

call check(self)

call random_vector(self%fld(:,:,:))

end subroutine random

! ------------------------------------------------------------------------------

subroutine copy(self,rhs)
implicit none
type(mpas_field), intent(inout) :: self
type(mpas_field), intent(in)    :: rhs
integer :: nf

call check_resolution(self, rhs)

nf = common_vars(self, rhs)
self%fld(:,:,1:nf) = rhs%fld(:,:,1:nf)
if (self%nf>nf) self%fld(:,:,nf+1:self%nf) = 0.0_kind_real

return
end subroutine copy

! ------------------------------------------------------------------------------

subroutine self_add(self,rhs)
implicit none
type(mpas_field), intent(inout) :: self
type(mpas_field), intent(in)    :: rhs
integer :: nf

call check_resolution(self, rhs)

nf = common_vars(self, rhs)
self%fld(:,:,1:nf) = self%fld(:,:,1:nf) + rhs%fld(:,:,1:nf)

return
end subroutine self_add

! ------------------------------------------------------------------------------

subroutine self_schur(self,rhs)
implicit none
type(mpas_field), intent(inout) :: self
type(mpas_field), intent(in)    :: rhs
integer :: nf

call check_resolution(self, rhs)

nf = common_vars(self, rhs)
self%fld(:,:,1:nf) = self%fld(:,:,1:nf) * rhs%fld(:,:,1:nf)

return
end subroutine self_schur

! ------------------------------------------------------------------------------

subroutine self_sub(self,rhs)
implicit none
type(mpas_field), intent(inout) :: self
type(mpas_field), intent(in)    :: rhs
integer :: nf

call check_resolution(self, rhs)

nf = common_vars(self, rhs)
self%fld(:,:,1:nf) = self%fld(:,:,1:nf) - rhs%fld(:,:,1:nf)

return
end subroutine self_sub

! ------------------------------------------------------------------------------

subroutine self_mul(self,zz)
implicit none
type(mpas_field), intent(inout) :: self
real(kind=kind_real), intent(in) :: zz

call check(self)

self%fld(:,:,:) = zz * self%fld(:,:,:)

return
end subroutine self_mul

! ------------------------------------------------------------------------------

subroutine axpy(self,zz,rhs)
implicit none
type(mpas_field), intent(inout) :: self
real(kind=kind_real), intent(in) :: zz
type(mpas_field), intent(in)    :: rhs
integer :: nf

call check_resolution(self, rhs)

nf = common_vars(self, rhs)
self%fld(:,:,1:nf) = self%fld(:,:,1:nf) + zz * rhs%fld(:,:,1:nf)

return
end subroutine axpy

! ------------------------------------------------------------------------------

subroutine dot_prod(fld1,fld2,zprod)
implicit none
type(mpas_field),     intent(in)    :: fld1, fld2
real(kind=kind_real), intent(inout) :: zprod

integer :: jx,jy,jz,jj
integer :: nl, nC, nf
real(kind=kind_real), allocatable :: zz(:), zv(:,:)

associate(nl=>fld1%geom%nVertLevels, nC=>fld1%geom%nCells, nf=>fld1%nf)

call check_resolution(fld1, fld2)
if (fld1%nf /= fld2%nf .or. fld1%geom%nVertLevels /= fld2%geom%nVertLevels) then
  call abor1_ftn("mpas_fields:field_prod error number of fields or vertical levels")
endif

allocate(zz(nl*nf))
allocate(zv(nl,nf))
zv(:,:)=0.0_kind_real
do jx=1,nC
   zv(:,:) = zv(:,:) + fld1%fld(jx,:,:) * fld2%fld(jx,:,:)
enddo

jj=1
do jz=1,nf
do jy=1,nl
  zz(jj) = zv(jy,jz)
  jj=jj+1
enddo
enddo

zprod=0.0_kind_real
do jj=1,nl*nf
  zprod=zprod+zz(jj)
enddo
deallocate(zz)
deallocate(zv)

end associate

return

end subroutine dot_prod

! ------------------------------------------------------------------------------

subroutine add_incr(self,rhs)
implicit none
type(mpas_field), intent(inout) :: self
type(mpas_field), intent(in)    :: rhs

call check(self)
call check(rhs)

if (self%geom%nCells==rhs%geom%nCells .and. self%geom%nVertLevels==rhs%geom%nVertLevels) then
    self%fld(:,:,:) = self%fld(:,:,:) + rhs%fld(:,:,:)
else
  call abor1_ftn("mpas_fields:add_incr: dimension mismatch")
endif

! FIXME - not sure what to do with this. (SH)
!if (self%nf>1) self%fld(:,:,2:) = 0.0_kind_real

return
end subroutine add_incr

! ------------------------------------------------------------------------------

subroutine diff_incr(lhs,x1,x2)
implicit none
type(mpas_field), intent(inout) :: lhs
type(mpas_field), intent(in)    :: x1
type(mpas_field), intent(in)    :: x2

call check(lhs)
call check(x1)
call check(x2)

call zeros(lhs)
if (x1%geom%nCells==x2%geom%nCells .and. x1%geom%nVertLevels==x2%geom%nVertLevels) then
  if (lhs%geom%nCells==x1%geom%nCells .and. lhs%geom%nVertLevels==x1%geom%nVertLevels) then
      lhs%fld(:,:,:) = x1%fld(:,:,:) - x2%fld(:,:,:)
  else
    call abor1_ftn("mpas_fields:diff_incr: dimension mismatch between the two variables.")
  endif
else
  call abor1_ftn("mpas_fields:diff_incr: states not at same resolution")
endif

return
end subroutine diff_incr

! ------------------------------------------------------------------------------

subroutine change_resol(fld,rhs)
implicit none
type(mpas_field), intent(inout) :: fld
type(mpas_field), intent(in)    :: rhs
real(kind=kind_real), allocatable :: ztmp(:,:)
real(kind=kind_real) :: dy1, dy2, ya, yb, dx1, dx2, xa, xb
integer :: jx, jy, jf, iy, ia, ib

call check(fld)
call check(rhs)

! FIXME: We just copy rhs to fld for now. Need an actual interpolation routine later. (SH)
if (fld%geom%nCells == rhs%geom%nCells .and.  fld%geom%nVertLevels == rhs%geom%nVertLevels) then
  call copy(fld, rhs)
else
  write(0,*) fld%geom%nCells, rhs%geom%nCells, fld%geom%nVertLevels, rhs%geom%nVertLevels
  call abor1_ftn("mpas_fields:field_resol: dimension mismatch")
endif

return
end subroutine change_resol

! ------------------------------------------------------------------------------

subroutine read_file(fld, c_conf, vdate)

use iso_c_binding
use datetime_mod
use fckit_log_module, only : log

implicit none
type(mpas_field), intent(inout) :: fld      !< Fields
type(c_ptr),      intent(in)    :: c_conf   !< Configuration
type(datetime),   intent(inout) :: vdate    !< DateTime

character(len=max_string_length+50) :: record
character(len=max_string_length) :: filename, string1
character(len=max_string_length)  :: varname
character(len=max_string_length)  :: sdate  
character(len=1024)  :: buf
!real(kind=kind_real), allocatable :: zz(:)
integer :: i, k, iread, nf
integer :: itime, istrlen
integer :: ncid, dimid, varid
integer :: ncells, nlevels, ndims !, nedges

iread = 1
if (config_element_exists(c_conf,"read_from_file")) then
  iread = config_get_int(c_conf,"read_from_file")
endif
!if (iread==0) then
!  call log%warning("mpas_fields:read_file: Inventing State")
!  call invent_state(fld,c_conf)
!  sdate = config_get_string(c_conf,len(sdate),"date")
!  WRITE(buf,*) 'validity date is: '//sdate
!  call log%info(buf)
!  call datetime_set(sdate, vdate)
!else
  call zeros(fld)
  filename = config_get_string(c_conf,len(filename),"filename")
  WRITE(buf,*) 'mpas_field:read_file: opening '//filename
  call log%info(buf)
  string1 = filename
  call ncerr(string1, nf90_open(trim(filename),nf90_nowrite,ncid))

  !> Grid dimensions
  call ncerr(string1, nf90_inq_dimid(ncid,'nCells',dimid))
  call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=ncells)) !self%nCells))
  !call ncerr(string1, nf90_inq_dimid(ncid,'nEdges',dimid))
  !call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=nedges)) !self%nEdges))
  call ncerr(string1, nf90_inq_dimid(ncid,'nVertLevels',dimid))
  call ncerr(string1, nf90_inquire_dimension(ncid,dimid,len=nlevels)) !self%nVertLevels))

  !> Check dimensions between mesh and fields
  if (ncells /= fld%geom%nCells .or.  nlevels /= fld%geom%nVertLevels) then
      write (record,*) "mpas_fields:read_file: wrong dimensions: ",fld%geom%nCells,fld%geom%nVertLevels
      call log%error(record)
      write (record,*) "mpas_fields:read_file: expected dimensions: ",ncells,nlevels
      call log%error(record)
      call abor1_ftn("mpas_fields:read_file: input fields have wrong dimensions")
  endif

  !> Time info
  call ncerr(string1, nf90_inq_varid(ncid,'xtime',varid))
  !call ncerr(string1, nf90_inquire_variable(ncid, varid, dimids=dimIDs, ndims=ndims))
  !if (ndims /= 2) then
  !    write(record,*) 'xtime variable has unknown shape in ', trim(filename)
  !    call log%error(record)
  !    call abor1_ftn("mpas_fields:read_file: error in reading xtime")
  !endif


  call ncerr(string1, nf90_inq_dimid(ncid, "Time", dimid))
  call ncerr(string1, nf90_inquire_dimension(ncid, dimid, len=itime))

  call ncerr(string1, nf90_inq_dimid(ncid, "StrLen", dimid))
  call ncerr(string1, nf90_inquire_dimension(ncid, dimid, len=istrlen))
  if (itime /= 1) then
      write(record,*) 'multiple timesteps (',itime,') in file ', trim(filename)
      write(record,*) 'We are using the LAST one, presumably, the LATEST timestep.'
      call log%error(record)
  endif

  call ncerr(string1, nf90_get_var(ncid, varid, sdate(1:istrlen), start = (/1, itime/), count = (/istrlen, 1 /)))
  sdate(11:11) = 'T'
  sdate(20:20) = 'Z'
  WRITE(buf,*) 'validity date is: '//trim(sdate)
  call log%info(buf)
  call datetime_set(sdate(1:20), vdate)

  !> Read variables fld (nCells, nVertLevels, nVariables)
  do i=1,fld%nf
     call ncerr(string1, nf90_inq_varid(ncid, trim(fld%fldnames(i)), varid))
     do k=1,fld%geom%nVertLevels
        call ncerr(string1, nf90_get_var  (ncid, varid, fld%fld(:,k,i), start = (/k,1,1/), count = (/1,fld%geom%nCells,1/) ))
     enddo
  enddo
!endif

call check(fld)

return
end subroutine read_file

! ------------------------------------------------------------------------------

subroutine write_file(fld, c_conf, vdate)
use iso_c_binding
use datetime_mod
use fckit_log_module, only : log

implicit none
type(mpas_field), intent(in) :: fld    !< Fields
type(c_ptr),      intent(in) :: c_conf !< Configuration
type(datetime),   intent(in) :: vdate  !< DateTime

character(len=max_string_length+50) :: record
character(len=max_string_length) :: gridfile, filename, string1, string2
character(len=20) :: sdate
character(len=1024):: buf
integer :: i, iread, nf
integer, dimension(fld%nf)  :: varid
integer :: fid, ncid, cellid, levelid, xtype, dimlen, dimid, vid
integer :: numdims, nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: ntimes, n
integer :: latid, lonid

call check(fld)

!> Write to an output netcdf file
filename = genfilename(c_conf,max_string_length,vdate)
WRITE(buf,*) 'mpas_field:write_file: writing '//trim(filename)
call log%info(buf)
string2 = trim(filename)
call ncerr(string2, nf90_create(trim(filename), nf90_clobber, ncid))

call ncerr(string2, nf90_def_dim(  ncid, 'Time',        nf90_unlimited, unlimitedDimID))
call ncerr(string2, nf90_def_dim(  ncid, 'nCells',      fld%geom%nCells,      cellid))
call ncerr(string2, nf90_def_dim(  ncid, 'nVertLevels', fld%geom%nVertLevels, levelid))

!> Time info
call datetime_to_string(vdate, sdate)
call ncerr(string2, nf90_put_att(ncid, nf90_global, 'date', sdate))

do i=1,fld%nf

   call ncerr(string2, nf90_def_var(ncid, trim(fld%fldnames(i)), nf90_double, &
                            (/ levelid, cellid, unlimitedDimID /), varid(i)))
!                            (/ levelid, cellid /), varid(i)))
end do

!> lat/lon
call ncerr(string2, nf90_def_var(ncid, "latCell", nf90_double, (/cellid/), latid))
call ncerr(string2, nf90_def_var(ncid, "lonCell", nf90_double, (/cellid/), lonid))

call ncerr(string2, nf90_enddef(ncid))

call ncerr(string2, nf90_put_var(ncid, latid, fld%geom%latCell))
call ncerr(string2, nf90_put_var(ncid, lonid, fld%geom%lonCell))

!allocate(self%fld(nC,nl,self%nf))
do i=1,fld%nf
   call ncerr(string2, nf90_put_var(ncid, varid(i), fld%fld(:,:,i), &
   start=(/1,1,1/), count=(/1,fld%geom%nCells,fld%geom%nVertLevels/)))
!   start=(/1,1/), count=(/fld%geom%nVertLevels,fld%geom%nCells/)))
enddo

call ncerr(string2, nf90_close(ncid))
return
end subroutine write_file

! ------------------------------------------------------------------------------

subroutine gpnorm(fld, nf, pstat)
implicit none
type(mpas_field),     intent(in)    :: fld
integer,              intent(in)    :: nf
real(kind=kind_real), intent(inout) :: pstat(3, nf)
integer :: jj
integer :: nl, nC

associate(nl=>fld%geom%nVertLevels, nC=>fld%geom%nCells)

call check(fld)   !(fld(nC,nl,nf))
if (fld%nf /= nf) then
    write(*,*) 'mpas_fields.F90 in gpnorm: ',fld%nf, nf
    call abor1_ftn("mpas_fields_gpnorm: error number of fields")
endif

do jj=1,fld%nf
  pstat(1,jj)=minval(fld%fld(:,:,jj))
  pstat(2,jj)=maxval(fld%fld(:,:,jj))
  pstat(3,jj)=sqrt(sum(fld%fld(:,:,jj)**2) &
               & /real(nl*nC,kind_real))
enddo

end associate

return
end subroutine gpnorm

! ------------------------------------------------------------------------------

subroutine fldrms(fld, prms)
implicit none
type(mpas_field), intent(in) :: fld
real(kind=kind_real), intent(out) :: prms
real(kind=kind_real) :: zz
integer :: jf,jy,jx,ii,nl,nC

associate(nl=>fld%geom%nVertLevels, nC=>fld%geom%nCells)

call check(fld)

zz = 0.0

do jf=1,fld%nf
  do jy=1,nl
    do jx=1,nC
      zz = zz + fld%fld(jx,jy,jf)*fld%fld(jx,jy,jf)
    enddo
  enddo
enddo

ii = nl*nC*fld%nf

prms = sqrt(zz/real(ii,kind_real))

end associate

end subroutine fldrms

! ------------------------------------------------------------------------------

subroutine lin_weights(kk,delta1,delta2,k1,k2,w1,w2)
implicit none
integer, intent(in)  :: kk
real(kind=kind_real), intent(in)     :: delta1,delta2
integer, intent(out) :: k1,k2
real(kind=kind_real), intent(out)    :: w1,w2

integer :: ii
real(kind=kind_real) :: zz

zz=real(kk-1,kind_real)*delta1
zz=zz/delta2
ii=int(zz)
w1=zz-real(ii,kind_real)
w2=1.0_kind_real-w1
k1=ii+1
k2=ii+2

return
end subroutine lin_weights

! ------------------------------------------------------------------------------

function genfilename (c_conf,length,vdate)
use iso_c_binding
use datetime_mod
use duration_mod
type(c_ptr), intent(in)    :: c_conf  !< Configuration
integer, intent(in) :: length
character(len=length) :: genfilename
type(datetime), intent(in) :: vdate

character(len=length) :: fdbdir, expver, typ, validitydate, referencedate, sstep, &
                       & prefix, mmb
type(datetime) :: rdate
type(duration) :: step
integer lenfn

! here we should query the length and then allocate "string".
! But Fortran 90 does not allow variable-length allocatable strings.
! config_get_string checks the string length and aborts if too short.
fdbdir = config_get_string(c_conf,len(fdbdir),"datadir")
expver = config_get_string(c_conf,len(expver),"exp")
typ    = config_get_string(c_conf,len(typ)   ,"type")

if (typ=="ens") then
  mmb = config_get_string(c_conf, len(mmb), "member")
  lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ) + 1 + LEN_TRIM(mmb)
  prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ) // "." // TRIM(mmb)
else
  lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ)
  prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ)
endif

if (typ=="fc" .or. typ=="ens") then
  referencedate = config_get_string(c_conf,len(referencedate),"date")
  call datetime_to_string(vdate, validitydate)
  call datetime_create(TRIM(referencedate), rdate)
  call datetime_diff(vdate, rdate, step)
  call duration_to_string(step, sstep)
  lenfn = lenfn + 1 + LEN_TRIM(referencedate) + 1 + LEN_TRIM(sstep)
  genfilename = TRIM(prefix) // "." // TRIM(referencedate) // "." // TRIM(sstep)
endif

if (typ=="an") then
  call datetime_to_string(vdate, validitydate)
  lenfn = lenfn + 1 + LEN_TRIM(validitydate)
  genfilename = TRIM(prefix) // "." // TRIM(validitydate) // ".nc"
endif

! SH
!referencedate = config_get_string(c_conf,len(referencedate),"date")
!lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(referencedate) + 1 + LEN_TRIM(typ)
!prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(referencedate) 
!genfilename = TRIM(prefix) // "." // TRIM(typ)

if (lenfn>length) &
  & call abor1_ftn("mpas_fields:genfilename: filename too long")

end function genfilename

! ------------------------------------------------------------------------------

function common_vars(x1, x2)

implicit none
type(mpas_field), intent(in) :: x1, x2
integer :: common_vars
integer :: jf

! We assume here that one set of fields is a subset of the other,
! that fields are always in the same order starting with x,
! and that the common fields are the first ones.

common_vars = min(x1%nf, x2%nf)
do jf = 1, common_vars
  if (x1%fldnames(jf)/=x2%fldnames(jf)) &
    & call abor1_ftn("common_vars: fields do not match")
enddo
if (x1%geom%nVertLevels /= x2%geom%nVertLevels) call abor1_ftn("common_vars: error number of levels")
!common_vars = x1%geom%nVertLevels * common_vars

end function common_vars

! ------------------------------------------------------------------------------

subroutine check_resolution(x1, x2)

implicit none
type(mpas_field), intent(in) :: x1, x2

if (x1%geom%nCells /= x2%geom%nCells .or.  x1%geom%nVertLevels /= x2%geom%nVertLevels) then
  call abor1_ftn ("mpas_fields: resolution error")
endif
call check(x1)
call check(x2)

end subroutine check_resolution

! ------------------------------------------------------------------------------

subroutine check(self)
implicit none
type(mpas_field), intent(in) :: self
logical :: bad
integer :: nl,nC

associate(nl=>self%geom%nVertLevels, nC=>self%geom%nCells)

bad = .not.associated(self%fld)

bad = bad .or. (size(self%fld, 1) /= nC)
bad = bad .or. (size(self%fld, 2) /= nL)
bad = bad .or. (size(self%fld, 3) /= self%nf)

bad = bad .or. .not.associated(self%th)

if (self%nf>1) then
  bad = bad .or. .not.associated(self%r)
  bad = bad .or. .not.associated(self%q)
  bad = bad .or. .not.associated(self%u)
  bad = bad .or. .not.associated(self%v)
else
  bad = bad .or. associated(self%r)
  bad = bad .or. associated(self%q)
  bad = bad .or. associated(self%u)
  bad = bad .or. associated(self%v)
endif

!allocate(self%fldnames(self%nf))

if (bad) then
  write(0,*)'nCells, nVertLevels, nVariables = ',nC,nL,self%nf
  write(0,*) '3-D Variable names = ', self%fldnames(:)
  if (associated(self%fld)) write(0,*)'shape(fld) = ',shape(self%fld)
  call abor1_ftn ("mpas_fields: field not consistent")
endif

end associate

end subroutine check

! ------------------------------------------------------------------------------

subroutine convert_to_ug(self, ug)
use unstructured_grid_mod
implicit none
type(mpas_field),        intent(in) :: self
type(unstructured_grid), intent(inout) :: ug
real(kind=kind_real), allocatable :: zz(:)
integer,              allocatable :: cmask(:)
integer                           :: i,j,k,v,nl,nC

associate(nl=>self%geom%nVertLevels, nC=>self%geom%nCells)
allocate(zz(nl))
allocate(cmask(nl))

do k = 1, nl
   zz(k) = real(k,kind=kind_real)
enddo
cmask(:) = 1      ! SH: We do not need a binary option (0 or 1) for MPAS. Just set to 1 for all.

call create_unstructured_grid(ug, nl, zz)

do i=1,nC
  call add_column(ug, self%geom%latCell(i), self%geom%lonCell(i), self%geom%areaCell(i), nl, self%nf, self%ns, cmask, 1)
  j = 0
  do v=1,self%nf
  do k=1,nl
       j = j+1
       ug%last%column%fld3d(j) = self%fld(i,k,v)
  enddo
  enddo
enddo

deallocate(zz)
deallocate(cmask)
end associate

end subroutine convert_to_ug

! ------------------------------------------------------------------------------

subroutine convert_from_ug(self, ug)
use unstructured_grid_mod
implicit none
type(mpas_field), intent(inout) :: self
type(unstructured_grid), intent(in) :: ug
type(column_element), pointer :: current
integer :: i,j,k,v,nl,nC

associate(nl=>self%geom%nVertLevels, nC=>self%geom%nCells)

current => ug%head
do i=1,nC
   j = 0
   do v=1,self%nf
   do k=1,nl
      j = j+1
      self%fld(i,k,v) = current%column%fld3d(j)
   enddo
   enddo
   current => current%next
enddo

end associate

end subroutine convert_from_ug

! ------------------------------------------------------------------------------

subroutine dirac(self, c_conf)
use iso_c_binding
implicit none
type(mpas_field), intent(inout) :: self
type(c_ptr), intent(in)       :: c_conf   !< Configuration
integer :: ndir,idir,ildir,ifdir,ioff
integer,allocatable :: iCell(:)
character(len=3) :: idirchar

call check(self)

! Get Diracs positions
ndir = config_get_int(c_conf,"ndir")

allocate(iCell(ndir))
do idir=1,ndir
   write(idirchar,'(i3)') idir
   iCell(idir) = config_get_int(c_conf,"iCell("//trim(adjustl(idirchar))//")")
end do
ildir = config_get_int(c_conf,"ildir")
ifdir = config_get_int(c_conf,"ifdir")

! Check
if (ndir<1) call abor1_ftn("mpas_fields:dirac non-positive ndir")
if (any(iCell<1).or.any(iCell>self%geom%nCells)) call abor1_ftn("mpas_fields:dirac invalid iCell")
if ((ildir<1).or.(ildir>self%geom%nVertLevels))  call abor1_ftn("mpas_fields:dirac invalid ildir")
if ((ifdir<1).or.(ifdir>self%nf)) call abor1_ftn("mpas_fields:dirac invalid ifdir")

! Setup Diracs
call zeros(self)
!ioff = (ifdir-1)*self%geom%nVertLevels
do idir=1,ndir
   self%fld(iCell(idir),ildir,ifdir) = 1.0
end do

end subroutine dirac

! ------------------------------------------------------------------------------

end module mpas_fields
