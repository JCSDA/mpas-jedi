! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_constants_mod

!oops
use kinds, only : kind_real

!MPAS-Model
use mpas_constants

implicit none

public

! ------------------------------------------------------------------------------

real(kind=kind_real), parameter :: deg2rad = pii/180.0_kind_real !-BJJ: TODO:  To-be-removed, when MPAS-release updated from Gael.

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

end module mpas_constants_mod
