! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------
subroutine c_mpas_init(c_conf) bind(c,name='mpas_init_f90')
use iso_c_binding
!use mpi
!use gungho_driver, only: mpas_init
use config_mod
!use mpas_derived_types
use mpas_dmpar, only: mpas_mpi_init

implicit none
type(c_ptr), intent(in) :: c_conf

integer, parameter               :: max_string_length=800 ! Yuk!
character(len=max_string_length) :: filename

integer :: length, mpi_ierr, comm

write(0,*)'MPI MPAS initialized'
call mpas_mpi_init()

end subroutine c_mpas_init

! ------------------------------------------------------------------------------

subroutine c_mpi_finalize() bind(c,name='mpi_finalize_f90')

use mpas_dmpar, only: mpas_mpi_finalize

implicit none

write(0,*)'finalize MPI from run_interface'
!call mpas_mpi_finalize()

end subroutine c_mpi_finalize
