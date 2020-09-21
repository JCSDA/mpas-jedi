module mpas_run_mod

use fckit_configuration_module, only: fckit_configuration
use mpas_subdriver

implicit none

private

contains

! ------------------------------------------------------------------------------
!> MPAS-specific Run.cc initialization
!!
!! \details **mpas_run_init()** ensures the correct namelist.atmosphere is used
!!
subroutine mpas_run_init(c_conf) bind(c,name='mpas_run_init_f90')

   use fckit_mpi_module, only: fckit_mpi_comm
   use, intrinsic :: iso_c_binding, only: c_ptr

   implicit none
   type(c_ptr), value, intent(in) :: c_conf
   character(len=30) :: fn
   type(fckit_configuration) :: f_conf
   character(len=:), allocatable :: str

   f_conf = fckit_configuration(c_conf)
   if(f_conf%has("nml_file")) then
     call f_conf%get_or_die("nml_file",str)
     fn = str
     call system("rm namelist.atmosphere")
     write(*,*) " fn = "//trim(fn)
     call system("cp "//trim(fn)//" namelist.atmosphere")
   endif

end subroutine mpas_run_init

! ------------------------------------------------------------------------------
!> MPAS-specific Run.cc finalize
!!
!! \details **mpas_run_final()** is empty
!!
subroutine mpas_run_final() bind(c,name='mpas_run_final_f90')

implicit none

end subroutine mpas_run_final

! ------------------------------------------------------------------------------

end module mpas_run_mod
