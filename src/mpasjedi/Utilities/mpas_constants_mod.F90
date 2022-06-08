! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_constants_mod

!oops
use kinds, only : kind_real

!MPAS-Model
use mpas_constants, only: pii, rgas, rv

implicit none

public

! ------------------------------------------------------------------------------

!Strings
character(len=3), parameter :: MPAS_JEDI_OFF = 'off'

! ------------------------------------------------------------------------------

!Commonly used numbers
real(kind=kind_real), parameter :: MPAS_JEDI_ZERO_kr     = 0.0_kind_real
real(kind=kind_real), parameter :: MPAS_JEDI_HALF_kr     = 0.5_kind_real
real(kind=kind_real), parameter :: MPAS_JEDI_ONE_kr      = 1.0_kind_real
real(kind=kind_real), parameter :: MPAS_JEDI_TWO_kr      = 2.0_kind_real
real(kind=kind_real), parameter :: MPAS_JEDI_THREE_kr    = 3.0_kind_real
real(kind=kind_real), parameter :: MPAS_JEDI_HUNDRED_kr  = 100.0_kind_real
real(kind=kind_real), parameter :: MPAS_JEDI_THOUSAND_kr = 1000.0_kind_real
real(kind=kind_real), parameter :: MPAS_JEDI_MILLION_kr  = 1000000.0_kind_real

! ------------------------------------------------------------------------------

!Geometry
real(kind=kind_real), parameter :: MPAS_JEDI_PII_kr      = real(pii,kind_real)
real(kind=kind_real), parameter :: MPAS_JEDI_PIIo2_kr    = MPAS_JEDI_PII_kr/MPAS_JEDI_TWO_kr
real(kind=kind_real), parameter :: MPAS_JEDI_DEG2RAD_kr  = MPAS_JEDI_PII_kr/180.0_kind_real
real(kind=kind_real), parameter :: MPAS_JEDI_RAD2DEG_kr  = 180.0_kind_real/MPAS_JEDI_PII_kr

! ------------------------------------------------------------------------------

!Comparison
real(kind=kind_real), parameter :: MPAS_JEDI_GREATERZERO_kr = 0.0000000000000001_kind_real
real(kind=kind_real), parameter :: MPAS_JEDI_LESSONE_kr     = 0.9999999999999999_kind_real

! ------------------------------------------------------------------------------

! reference pressure p0 in Pa
real(kind=kind_real), parameter :: MPAS_JEDI_P0_kr  = 100000.0_kind_real

! ------------------------------------------------------------------------------

! For relative humidity conversion
real(kind=kind_real), parameter :: rd_over_rv  = rgas/rv
real(kind=kind_real), parameter :: rd_over_rv1 = MPAS_JEDI_ONE_kr - rd_over_rv
real(kind=kind_real), parameter :: t_kelvin    = 273.15_kind_real
!Saturation Vapour Pressure Constants(Rogers & Yau, 1989)
real(kind=kind_real), parameter :: es_alpha    = 611.2_kind_real
real(kind=kind_real), parameter :: es_beta     = 17.67_kind_real
real(kind=kind_real), parameter :: es_gamma    = 243.5_kind_real

contains

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

end module mpas_constants_mod
