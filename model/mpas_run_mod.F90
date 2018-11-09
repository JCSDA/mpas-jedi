module mpas_run_mod

use iso_c_binding 
use mpas_subdriver
use config_mod

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

   implicit none
   type(c_ptr), intent(in) :: c_conf
   character(len=30) :: fn
   type(fckit_mpi_comm) :: f_comm
   f_comm = fckit_mpi_comm()

   if(config_element_exists(c_conf,"nml_file")) then
     fn = config_get_string(c_conf,len(fn),"nml_file")
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
