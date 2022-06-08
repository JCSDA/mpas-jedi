! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_kinds
  use, intrinsic :: iso_c_binding
  implicit none

  public kind_double, c_real_type
  private
  integer, parameter :: kind_double=c_double   ! for special micro-physics accuracy
  integer, parameter :: c_real_type=c_double
  ! note in oops ./util/kinds.F90 :  kind_float = c_float
end module mpas_kinds
