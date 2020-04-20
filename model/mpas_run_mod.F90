module mpas_run_mod

use fckit_configuration_module, only: fckit_configuration
use mpas_subdriver

implicit none

private

!NOTE:  Can we potentially have the multiple domain/corelist for different resolution
!       by calling "mpas_init"/"mpas_finalize" multiple times ?
public :: domain, corelist

type (domain_type), pointer :: domain => null() 
type (core_type), pointer :: corelist => null()

contains

! ------------------------------------------------------------------------------
!> MPAS-specific Run.cc initialization
!!
!! \details **mpas_run_init()** initializes the whole MPAS framework.
!!          "domain" and "corelist" are public variables in mpas_run_mod.
!!
subroutine mpas_run_init(c_conf) bind(c,name='mpas_run_init_f90')

   use fckit_mpi_module, only: fckit_mpi_comm
   use, intrinsic :: iso_c_binding, only: c_ptr

   implicit none
   type(c_ptr), value, intent(in) :: c_conf
   character(len=30) :: fn
   type(fckit_mpi_comm) :: f_comm
   type(fckit_configuration) :: f_conf
   character(len=:), allocatable :: str

   f_comm = fckit_mpi_comm()

   f_conf = fckit_configuration(c_conf)
   if(f_conf%has("nml_file")) then
     call f_conf%get_or_die("nml_file",str)
     fn = str
     call system("rm namelist.atmosphere")
     write(*,*) " fn = "//trim(fn)
     call system("cp "//trim(fn)//" namelist.atmosphere")
   endif

   !> MPAS subdriver
   call mpas_init( corelist, domain, mpi_comm=f_comm%communicator() )

end subroutine mpas_run_init

! ------------------------------------------------------------------------------
!> MPAS-specific Run.cc finalize
!!
!! \details **mpas_run_final()** finalizes the whole MPAS framework.
!!
subroutine mpas_run_final() bind(c,name='mpas_run_final_f90')

implicit none

   !> MPAS subdriver
   call mpas_timer_set_context( domain )
   call mpas_finalize( corelist, domain)

end subroutine mpas_run_final

! ------------------------------------------------------------------------------

end module mpas_run_mod
