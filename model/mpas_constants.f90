
!> Constants for the MPAS model 

module mpas_constants

use kinds
implicit none

!--- Dimensional parameters (copied from MPAS/src/framework/mpas_constants.F)

   real (kind=kind_real), parameter :: pii     = 3.141592653589793   !< Constant: Pi
   real (kind=kind_real), parameter :: raidus  = 6371229.0           !< Constant: Spherical Earth radius [m]
   real (kind=kind_real), parameter :: gravity = 9.80616             !< Constant: Acceleration due to gravity [m s-2]
   real (kind=kind_real), parameter :: rgas    = 287.0               !< Constant: Gas constant for dry air [J kg-1 K-1]
   real (kind=kind_real), parameter :: rv      = 461.6               !< Constant: Gas constant for water vapor [J kg-1 K-1]
   real (kind=kind_real), parameter :: rvord   = rv/rgas             ! => 1.608362
!  real (kind=kind_real), parameter :: cp      = 1003.0              !< Constant: Specific heat of dry air at constant pressure [J kg-1 K-1]
   real (kind=kind_real), parameter :: cp      = 7.*rgas/2.          !< Constant: Specific heat of dry air at constant pressure [J kg-1 K-1]
   real (kind=kind_real), parameter :: cv      = cp - rgas           !< Constant: Specific heat of dry air at constant volume [J kg-1 K-1]

!--- Dimensional parameters (copied from MPAS/src/core_atmosphere/physics/mpas_atmphys_constants.F)
   real(kind=kind_real), parameter :: deg2rad = pii/180.            !< conversion from degree to radiant  
   real(kind=kind_real), parameter :: p0      = 100000.             !< Constant: reference pressure 
!--- derived parameters
   real (kind=kind_real), parameter :: rcv    = rgas/(cp-rgas)      !< Constant

!--- Non-dimensional parameters

end module mpas_constants
