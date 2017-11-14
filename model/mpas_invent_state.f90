
!> Invent an initial state for the QG model.

!> This routine invent an initial state for the QG model. It is used to
!! initialise the "truth run". The initial state consists of a horizontally
!! uniform wind in each layer, with a vertical shear sufficient to produce
!! baroclinic instability. Povided the orography is non-zero and is not
!! symmetrically place in the domain, this is sufficient to generate a
!! non-trivial flow after a few days of integration.
!!
!! Two slightly different initial states may be created (according to whether
!! or not ctype is set to 'f').

subroutine invent_state(flds,config)

use mpas_fields
use iso_c_binding
use config_mod
use fckit_log_module, only : log
!use mpas_constants, only: u1,u2,bet,worog, domain_zonal, domain_meridional, &
!                      & dlogtheta,f0,g,horog,scale_length,rossby_number
use kinds

implicit none

type(mpas_field), intent(inout) :: flds    !< Model fields
type(c_ptr), intent(in)       :: config  !< Configuration structure


call log%info('Empty for now')
! ------------------------------------------------------------------------------

return
end subroutine invent_state
