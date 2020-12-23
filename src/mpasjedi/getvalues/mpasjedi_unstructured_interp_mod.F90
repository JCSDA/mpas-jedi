! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module mpasjedi_unstructured_interp_mod

use fckit_log_module,               only: fckit_log
use fckit_mpi_module,               only: fckit_mpi_comm, & 
                                       &  fckit_mpi_sum, &
                                       &  fckit_mpi_max, &
                                       &  fckit_mpi_min
! oops
use unstructured_interpolation_mod, only: unstrc_interp

!mpas-jedi
use mpas_constants_mod
use mpas_geom_mod

private 

public :: unsinterp_integer_apply, &
          unsinterp_nearest_apply

contains 

! ------------------------------------------------------------------------------

! This subroutine populates field_out(n) with the field_in variable value that
! has the greatest sum of weights contributed from all the neighbors of the
! output location. (Only works with a discrete-valued field_in)
! Note: assumes wtype='barycent'.
subroutine unsinterp_integer_apply(unsinterp, field_in, field_out)

  type(unstrc_interp), intent(in)    :: unsinterp
  real(kind=kind_real),    intent(in)    :: field_in(:) !Integer field in
  real(kind=kind_real),    intent(inout) :: field_out(:) !Integer field out

  integer :: maxtypel, mintypel, maxtype, mintype, ngrid_out
  integer :: i, j, k, n, index
  real(kind=kind_real), allocatable :: interp_w(:,:)
  real(kind=kind_real), allocatable :: field_out_tmp(:)
  real(kind=kind_real), allocatable :: field_neighbors(:,:)
  real(kind=kind_real), allocatable :: field_types(:)

  ! Inteprolation of integer fields

  ! Size of output
  ngrid_out = size(field_out)

  ! Get nearest neighbors
  allocate(field_neighbors(unsinterp%nn,ngrid_out))
  allocate(field_out_tmp(ngrid_out))
  call unsinterp%apply(field_in, field_out_tmp, field_neighbors)

  ! Global min and max integers in field
  maxtypel = int(maxval(field_in))
  mintypel = int(minval(field_in))
  call unsinterp%comm%allreduce(maxtypel,maxtype,fckit_mpi_max())
  call unsinterp%comm%allreduce(mintypel,mintype,fckit_mpi_min())

  ! Put weights into field type array and pick max for interpolated value
  allocate(field_types(mintype:maxtype))

  field_out = 0.0_kind_real
  do i = 1,ngrid_out
    field_types = 0.0
    do n = 1, unsinterp%nn
      index = int(field_neighbors(n,i))
      field_types(index) = field_types(index) + unsinterp%interp_w(n,i)
    enddo
    field_out(i) = real(maxloc(field_types,1)+(mintype-1),kind_real)
  enddo

end subroutine unsinterp_integer_apply

! ------------------------------------------------------------------------------

! This subroutine populates field_out(n) with the field_in variable value contained
! in the nearest neighbor (i.e. neighbor with greatest weight) of the output location.
! (Would work with either discrete- or continuous-valued field_in)
! Note: assumes wtype='barycent'.
subroutine unsinterp_nearest_apply(unsinterp, field_in, field_out)

  type(unstrc_interp), intent(in)    :: unsinterp
  real(kind=kind_real),    intent(in)    :: field_in(:) !Integer field in
  real(kind=kind_real),    intent(inout) :: field_out(:) !Integer field out

  integer :: n, ngrid_out
  real(kind=kind_real), allocatable :: field_out_tmp(:)
  real(kind=kind_real), allocatable :: field_neighbors(:,:)

  ! Inteprolation using the nearest neighbor

  ! Size of output
  ngrid_out = size(field_out)

  ! Get nearest neighbors
  allocate(field_neighbors(unsinterp%nn,ngrid_out))
  allocate(field_out_tmp(ngrid_out))
  call unsinterp%apply(field_in, field_out_tmp, field_neighbors)

  ! Find nearest neighbor
  do n = 1, ngrid_out
    ! The original code from lfric/fv3 had 'minloc' in the line below. Seems like a bug.
    field_out(n) = field_neighbors(maxloc(unsinterp%interp_w(:,n),1),n)
  enddo

end subroutine unsinterp_nearest_apply

! ------------------------------------------------------------------------------
end module mpasjedi_unstructured_interp_mod
