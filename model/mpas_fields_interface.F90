
!> Interfaces to be called from C++ for Fortran handling of QG model fields

! ------------------------------------------------------------------------------

subroutine mpas_field_create_c(c_key_self, c_key_geom, c_key_vars) bind(c,name='mpas_field_create_f90')
use iso_c_binding
use mpas_fields
use mpas_geom_mod
use mpas_vars_mod
implicit none
integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(in) :: c_key_geom !< Geometry
integer(c_int), intent(in) :: c_key_vars !< List of variables

type(mpas_field), pointer :: self
type(mpas_geom),  pointer :: geom
type(mpas_vars),  pointer :: vars

call mpas_geom_registry%get(c_key_geom, geom)
call mpas_vars_registry%get(c_key_vars, vars)
call mpas_field_registry%init()
call mpas_field_registry%add(c_key_self)
call mpas_field_registry%get(c_key_self,self)

call create(self, geom, vars)

end subroutine mpas_field_create_c

! ------------------------------------------------------------------------------

subroutine mpas_field_delete_c(c_key_self) bind(c,name='mpas_field_delete_f90')
use iso_c_binding
use mpas_fields
implicit none
integer(c_int), intent(inout) :: c_key_self
type(mpas_field), pointer :: self

call mpas_field_registry%get(c_key_self,self)

call delete(self)

call mpas_field_registry%remove(c_key_self)

end subroutine mpas_field_delete_c

! ------------------------------------------------------------------------------

subroutine mpas_field_zero_c(c_key_self) bind(c,name='mpas_field_zero_f90')
use iso_c_binding
use mpas_fields
implicit none
integer(c_int), intent(in) :: c_key_self
type(mpas_field), pointer :: self

call mpas_field_registry%get(c_key_self,self)
call zeros(self)

end subroutine mpas_field_zero_c

! ------------------------------------------------------------------------------

subroutine mpas_field_dirac_file_c(c_key_fld, c_conf) bind(c,name='mpas_field_dirac_f90')
use iso_c_binding
use mpas_fields

implicit none
integer(c_int), intent(in) :: c_key_fld  !< Fields
type(c_ptr), intent(in)    :: c_conf !< Configuration

type(mpas_field), pointer :: fld

call mpas_field_registry%get(c_key_fld,fld)
call dirac(fld, c_conf)

end subroutine mpas_field_dirac_file_c

! ------------------------------------------------------------------------------

subroutine mpas_field_random_c(c_key_self) bind(c,name='mpas_field_random_f90')
use iso_c_binding
use mpas_fields
implicit none
integer(c_int), intent(in) :: c_key_self
type(mpas_field), pointer :: self

call mpas_field_registry%get(c_key_self,self)
call random(self)

end subroutine mpas_field_random_c

! ------------------------------------------------------------------------------

subroutine mpas_field_copy_c(c_key_self,c_key_rhs) bind(c,name='mpas_field_copy_f90')
use iso_c_binding
use mpas_fields
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(mpas_field), pointer :: self
type(mpas_field), pointer :: rhs
call mpas_field_registry%get(c_key_self,self)
call mpas_field_registry%get(c_key_rhs,rhs)

call copy(self, rhs)

end subroutine mpas_field_copy_c

! ------------------------------------------------------------------------------

subroutine mpas_field_self_add_c(c_key_self,c_key_rhs) bind(c,name='mpas_field_self_add_f90')
use iso_c_binding
use mpas_fields
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(mpas_field), pointer :: self
type(mpas_field), pointer :: rhs
call mpas_field_registry%get(c_key_self,self)
call mpas_field_registry%get(c_key_rhs,rhs)

call self_add(self,rhs)

end subroutine mpas_field_self_add_c

! ------------------------------------------------------------------------------

subroutine mpas_field_self_schur_c(c_key_self,c_key_rhs) bind(c,name='mpas_field_self_schur_f90')
use iso_c_binding
use mpas_fields
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(mpas_field), pointer :: self
type(mpas_field), pointer :: rhs
call mpas_field_registry%get(c_key_self,self)
call mpas_field_registry%get(c_key_rhs,rhs)

call self_schur(self,rhs)

end subroutine mpas_field_self_schur_c

! ------------------------------------------------------------------------------

subroutine mpas_field_self_sub_c(c_key_self,c_key_rhs) bind(c,name='mpas_field_self_sub_f90')
use iso_c_binding
use mpas_fields
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(mpas_field), pointer :: self
type(mpas_field), pointer :: rhs
call mpas_field_registry%get(c_key_self,self)
call mpas_field_registry%get(c_key_rhs,rhs)

call self_sub(self,rhs)

end subroutine mpas_field_self_sub_c

! ------------------------------------------------------------------------------

subroutine mpas_field_self_mul_c(c_key_self,c_zz) bind(c,name='mpas_field_self_mul_f90')
use iso_c_binding
use mpas_fields
use kinds
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
type(mpas_field), pointer :: self
real(kind=kind_real) :: zz

call mpas_field_registry%get(c_key_self,self)
zz = c_zz

call self_mul(self,zz)

end subroutine mpas_field_self_mul_c

! ------------------------------------------------------------------------------

subroutine mpas_field_axpy_c(c_key_self,c_zz,c_key_rhs) bind(c,name='mpas_field_axpy_f90')
use iso_c_binding
use mpas_fields
use kinds
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
integer(c_int), intent(in) :: c_key_rhs

type(mpas_field), pointer :: self
type(mpas_field), pointer :: rhs
real(kind=kind_real) :: zz

call mpas_field_registry%get(c_key_self,self)
call mpas_field_registry%get(c_key_rhs,rhs)
zz = c_zz

call axpy(self,zz,rhs)

end subroutine mpas_field_axpy_c

! ------------------------------------------------------------------------------

subroutine mpas_field_dot_prod_c(c_key_fld1,c_key_fld2,c_prod) bind(c,name='mpas_field_dot_prod_f90')
use iso_c_binding
use mpas_fields
use kinds
implicit none
integer(c_int), intent(in)    :: c_key_fld1, c_key_fld2
real(c_double), intent(inout) :: c_prod
real(kind=kind_real) :: zz
type(mpas_field), pointer :: fld1, fld2

call mpas_field_registry%get(c_key_fld1,fld1)
call mpas_field_registry%get(c_key_fld2,fld2)

call dot_prod(fld1,fld2,zz)

c_prod = zz

end subroutine mpas_field_dot_prod_c

! ------------------------------------------------------------------------------

subroutine mpas_field_add_incr_c(c_key_self,c_key_rhs) bind(c,name='mpas_field_add_incr_f90')
use iso_c_binding
use mpas_fields
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs
type(mpas_field), pointer :: self
type(mpas_field), pointer :: rhs

call mpas_field_registry%get(c_key_self,self)
call mpas_field_registry%get(c_key_rhs,rhs)

call add_incr(self,rhs)

end subroutine mpas_field_add_incr_c

! ------------------------------------------------------------------------------

subroutine mpas_field_diff_incr_c(c_key_lhs,c_key_x1,c_key_x2) bind(c,name='mpas_field_diff_incr_f90')
use iso_c_binding
use mpas_fields
implicit none
integer(c_int), intent(in) :: c_key_lhs
integer(c_int), intent(in) :: c_key_x1
integer(c_int), intent(in) :: c_key_x2
type(mpas_field), pointer :: lhs
type(mpas_field), pointer :: x1
type(mpas_field), pointer :: x2

call mpas_field_registry%get(c_key_lhs,lhs)
call mpas_field_registry%get(c_key_x1,x1)
call mpas_field_registry%get(c_key_x2,x2)

call diff_incr(lhs,x1,x2)

end subroutine mpas_field_diff_incr_c

! ------------------------------------------------------------------------------

subroutine mpas_field_change_resol_c(c_key_fld,c_key_rhs) bind(c,name='mpas_field_change_resol_f90')
use iso_c_binding
use mpas_fields
implicit none
integer(c_int), intent(in) :: c_key_fld
integer(c_int), intent(in) :: c_key_rhs
type(mpas_field), pointer :: fld, rhs

call mpas_field_registry%get(c_key_fld,fld)
call mpas_field_registry%get(c_key_rhs,rhs)

call change_resol(fld,rhs)

end subroutine mpas_field_change_resol_c

! ------------------------------------------------------------------------------

subroutine mpas_field_read_file_c(c_key_fld, c_conf, c_dt) bind(c,name='mpas_field_read_file_f90')
use iso_c_binding
use mpas_fields
use datetime_mod

implicit none
integer(c_int), intent(in) :: c_key_fld  !< Fields
type(c_ptr), intent(in)    :: c_conf !< Configuration
type(c_ptr), intent(inout) :: c_dt   !< DateTime

type(mpas_field), pointer :: fld
type(datetime) :: fdate

call mpas_field_registry%get(c_key_fld,fld)
call c_f_datetime(c_dt, fdate)
call read_file(fld, c_conf, fdate)

end subroutine mpas_field_read_file_c

! ------------------------------------------------------------------------------

subroutine mpas_field_write_file_c(c_key_fld, c_conf, c_dt) bind(c,name='mpas_field_write_file_f90')
use iso_c_binding
use mpas_fields
use datetime_mod

implicit none
integer(c_int), intent(in) :: c_key_fld  !< Fields
type(c_ptr), intent(in) :: c_conf !< Configuration
type(c_ptr), intent(in) :: c_dt   !< DateTime

type(mpas_field), pointer :: fld
type(datetime) :: fdate

call mpas_field_registry%get(c_key_fld,fld)
call c_f_datetime(c_dt, fdate)
call write_file(fld, c_conf, fdate)

end subroutine mpas_field_write_file_c

! ------------------------------------------------------------------------------

subroutine mpas_field_gpnorm_c(c_key_fld, kf, pstat) bind(c,name='mpas_field_gpnorm_f90')
use iso_c_binding
use mpas_fields
use kinds
implicit none
integer(c_int), intent(in) :: c_key_fld
integer(c_int), intent(in) :: kf
real(c_double), intent(inout) :: pstat(3*kf)

type(mpas_field), pointer :: fld
real(kind=kind_real) :: zstat(3, kf)
integer :: jj, js, jf

call mpas_field_registry%get(c_key_fld,fld)

call gpnorm(fld, kf, zstat)
jj=0
do jf = 1, kf
  do js = 1, 3
    jj=jj+1
    pstat(jj) = zstat(js,jf)
  enddo
enddo

end subroutine mpas_field_gpnorm_c

! ------------------------------------------------------------------------------

subroutine mpas_field_rms_c(c_key_fld, prms) bind(c,name='mpas_field_rms_f90')
use iso_c_binding
use mpas_fields
use kinds
implicit none
integer(c_int), intent(in) :: c_key_fld
real(c_double), intent(inout) :: prms

type(mpas_field), pointer :: fld
real(kind=kind_real) :: zz

call mpas_field_registry%get(c_key_fld,fld)

call fldrms(fld, zz)

prms = zz

end subroutine mpas_field_rms_c

! ------------------------------------------------------------------------------

subroutine mpas_fieldnum_c(c_key_fld, nx, ny, nf) bind(c,name='mpas_field_sizes_f90')
use iso_c_binding
use mpas_fields
implicit none
integer(c_int), intent(in) :: c_key_fld
integer(kind=c_int), intent(inout) :: nx, ny, nf
type(mpas_field), pointer :: fld

call mpas_field_registry%get(c_key_fld,fld)

nx = fld%geom%nCells
ny = fld%geom%nVertLevels
nf = fld%nf

end subroutine mpas_fieldnum_c

! ------------------------------------------------------------------------------

subroutine mpas_field_convert_to_c(c_key_fld, c_key_ug) bind (c,name='mpas_field_convert_to_f90')
use iso_c_binding
use mpas_fields
use unstructured_grid_mod
implicit none
integer(c_int), intent(in) :: c_key_fld
integer(c_int), intent(in) :: c_key_ug
type(mpas_field), pointer :: fld
type(unstructured_grid), pointer :: ug

call mpas_field_registry%get(c_key_fld,fld)
call unstructured_grid_registry%get(c_key_ug,ug)

call convert_to_ug(fld, ug)

end subroutine mpas_field_convert_to_c
! ------------------------------------------------------------------------------
subroutine mpas_field_convert_from_c(c_key_fld, c_key_ug) bind (c,name='mpas_field_convert_from_f90')
use iso_c_binding
use mpas_fields
use unstructured_grid_mod
implicit none
integer(c_int), intent(in) :: c_key_fld
integer(c_int), intent(in) :: c_key_ug
type(mpas_field), pointer :: fld
type(unstructured_grid), pointer :: ug

call mpas_field_registry%get(c_key_fld,fld)
call unstructured_grid_registry%get(c_key_ug,ug)

call convert_from_ug(fld, ug)

end subroutine mpas_field_convert_from_c
! ------------------------------------------------------------------------------
