! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_fields_mod

!use iso_c_binding , only: c_char, c_int
use iso_c_binding 
use config_mod
use datetime_mod
use mpas_geom_mod
use ufo_vars_mod
use mpas_kinds, only : kind_real
use ioda_locs_mod
use ufo_geovals_mod
use mpas_getvaltraj_mod, only: mpas_getvaltraj

use mpas_dmpar
use mpas_derived_types
use mpas_framework
use mpas_kind_types
!use init_atm_core_interface
use mpas_subdriver
use atm_core
use mpas4da_mod
use mpas2ufo_vars_mod
use mpi ! only MPI_COMM_WORLD
use mpas_stream_manager
use mpas_pool_routines
use mpas_field_routines
use mpas_constants

implicit none
private

public :: mpas_field, &
        & create, delete, zeros, random, copy, &
        & self_add, self_schur, self_sub, self_mul, axpy, &
        & dot_prod, add_incr, diff_incr, &
        & read_file, write_file, gpnorm, fldrms, &
        & change_resol, getvalues, getvalues_tl, getvalues_ad, &
        & ug_coord, field_to_ug, field_from_ug, &
        & dirac, analytic_IC
public :: mpas_field_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold MPAS fields
type :: mpas_field
  !type (domain_type), pointer :: domain          ! NOW in geom, For convenience not used 
  !type (core_type), pointer :: corelist          ! NOW in geom, For convenience not used
  type (mpas_geom), pointer :: geom              ! grid and MPI infos
  integer :: nf                                  ! Number of variables in fld
  character(len=MAXVARLEN), allocatable  :: fldnames(:) ! Variable identifiers
  type (mpas_pool_type), pointer  :: subFields   !---> state variables (to be analyzed)
  type (mpas_pool_type), pointer  :: auxFields   !---> auxiliary variables, such as pressure, t2m, u10, v10, Tsfc
  type (MPAS_streamManager_type), pointer :: manager
  type (MPAS_Clock_type), pointer :: clock
end type mpas_field

real(kind=kind_real), parameter :: deg2rad = pii/180.0_kind_real !-BJJ: TODO:  To-be-removed, when MPAS-release updated from Gael.
 
#define LISTED_TYPE mpas_field

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: mpas_field_registry

integer, parameter :: nf_aux = 21

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine create(self, geom, vars)

    use mpas_kind_types

    implicit none

    type(mpas_field), intent(inout)       :: self
    type(mpas_geom),  intent(in), pointer :: geom
    type(ufo_vars),   intent(in)          :: vars

    integer :: nsize, nfields
    integer :: ierr!, ii

    character(len=22), allocatable  :: fldnames_aux(:)

    !-- fortran level test (temporally sit here)
    real(kind=kind_real), allocatable :: pstat(:, :)
    real(kind=kind_real)              :: prms
    type (MPAS_Time_type) :: local_time, write_time, fld_time
    character (len=StrKIND) :: dateTimeString, dateTimeString2, streamID, time_string, filename
    character (len=StrKIND) :: dateTimeString_oops

!   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a
!   type (field2DReal), pointer :: field2d, field2d_src

    ! from the namelist
    self % nf =  vars % nv
    allocate(self % fldnames(self % nf))
    self % fldnames(:) = vars % fldnames(:)
    write(*,*)'self % nf =',self % nf
    write(*,*)'allocate ::',self % fldnames(:)
   
    ! link geom
    if (associated(geom)) then
      self % geom => geom
    else
      write(*,*)'mpas_fields: geom not associated'
      call abor1_ftn("mpas_fields: geom not associated")
    end if

    write(0,*)'-- Create a sub Pool from list of variable ',self % nf
    call da_make_subpool(self % geom % domain, self % subFields, self % nf, self % fldnames, nfields)

    if ( self % nf .ne. nfields  ) then
       call abor1_ftn("mpas_fields:create: dimension mismatch ",self % nf, nfields)
    end  if

    !--- TODO: aux test: BJJ  !- get this from json ???
    allocate(fldnames_aux(nf_aux))
    fldnames_aux = [ character(len=22) :: "theta", "rho", "u", &
                                          "landmask", "xice", "snowc", "skintemp", "ivgtyp", "isltyp", &
                                          "snowh", "vegfra", "u10", "v10", "lai", "smois", "tslb", "w", &
                                          "index_qc", "index_qi", "re_cloud", "re_ice" ]
                                          !BJJ- "w" is for dimension information of var_prsi
    write(0,*)'-- Create a sub Pool for auxFields'
    call da_make_subpool(self % geom % domain, self % auxFields, nf_aux, fldnames_aux, nfields)
    deallocate(fldnames_aux)
    if ( nf_aux .ne. nfields  ) then
       call abor1_ftn("mpas_fields:create: dimension mismatch ",nf_aux, nfields)
    end  if

    ! clock creation
    allocate(self % clock)
    call atm_simulation_clock_init(self % clock, self % geom % domain % blocklist % configs, ierr)
    if ( ierr .ne. 0 ) then
       call abor1_ftn("mpas_fields: atm_simulation_clock_init problem")
    end if

    !-------------------------------------------------------------
    ! Few temporary tests
    !-------------------------------------------------------------
!    call update_mpas_field(self % geom % domain, self % auxFields)
!!    call mpas_pool_get_subpool(self % geom % domain % blocklist % structs,'state',state)
!    call da_fldrms(self % subFields, self % geom % domain % dminfo, prms)
!    allocate(pstat(3, self % nf))
!    call da_gpnorm(self % subFields, self % geom % domain % dminfo, self % nf, pstat)
!    deallocate(pstat)
!!    call abor1_ftn("MPAS test")

     call zeros(self) !-- set zero for self % subFields

   ! mpas_get_time(curr_time, YYYY, MM, DD, DoY, H, M, S, S_n, S_d, dateTimeString, ierr)
   ! call mpas_get_time(curr_time, YYYY, MM, DD, H, M, S, dateTimeString, ierr)
   ! call mpas_set_timeInterval(runDuration, timeString=config_run_duration, ierr=local_err)
   ! call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=local_err)
   ! currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
   ! startTime = mpas_get_clock_time(clock, MPAS_START_TIME, ierr)
   ! currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
   ! xtimeTime = currTime - startTime
   ! call mpas_get_timeInterval(interval=xtimeTime, S=s, S_n=s_n, S_d=s_d, ierr=ierr)
   ! xtime_s = (s + s_n / s_d)


!--------------------------------------------------------------------------------------- 

!   write(0,*)''
!   write(0,*)'=========== MPAS WRITE A RESTART ===================='

   !call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
!   call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)
!   all mpas_pool_get_subpool(block_ptr % structs, 'diag_physics', diag_physics)
!   call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
!   call atm_compute_restart_diagnostics(state, 1, diag, mesh)
 
    return

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

   implicit none
   type(mpas_field), intent(inout) :: self
   integer :: ierr = 0 
  
   if (allocated(self % fldnames)) deallocate(self % fldnames)
   if (associated(self % subFields)) then
      write(*,*)'--> deallocate subFields Pool'
      call mpas_pool_empty_pool(self % subFields)
      call mpas_pool_destroy_pool(self % subFields)
   end if
   if (associated(self % auxFields)) then
      write(*,*)'--> deallocate auxFields Pool'
      call mpas_pool_empty_pool(self % auxFields)
      call mpas_pool_destroy_pool(self % auxFields)
   end if
   call mpas_destroy_clock(self % clock, ierr)
   if ( ierr .ne. 0  ) then
      write(*,*)'mpas_fields deallocate clock failed'
   end if
   write(*,*)'--> mpas_fields done deallocate'

   return

end subroutine delete

! ------------------------------------------------------------------------------

subroutine zeros(self)

   implicit none
   type(mpas_field), intent(inout) :: self

   call da_zeros(self % subFields)

end subroutine zeros

! ------------------------------------------------------------------------------

subroutine random(self)

   implicit none
   type(mpas_field), intent(inout) :: self
   
   call da_random(self % subFields)

end subroutine random

! ------------------------------------------------------------------------------

subroutine copy(self,rhs)

   implicit none
   type(mpas_field), intent(inout) :: self
   type(mpas_field), intent(in)    :: rhs
   
   write(*,*)'====> copy of mpas_field'

   !TODO: Do we need to empty/destroy subFields and re-create/clone it ?
   !    : Is "mpas_pool_clone_pool" enough?
   !    : Check this, considering the use of "copy" routine in OOPS.

   ! Duplicate the members of rhs into self and do a deep copy
   ! of the fields from self % subFields to rhs % subFields
   call mpas_pool_empty_pool(self % subFields)
   call mpas_pool_destroy_pool(self % subFields)
   self % nf = rhs % nf
   call mpas_pool_create_pool(self % subFields,self % nf)
   call mpas_pool_clone_pool(rhs % subFields, self % subFields)

   call mpas_pool_empty_pool(self % auxFields)
   call mpas_pool_destroy_pool(self % auxFields)
   call mpas_pool_create_pool(self % auxFields,nf_aux)
   call mpas_pool_clone_pool(rhs % auxFields, self % auxFields)

   ! We should consider adding a subroutine just updating the fields
   ! call mpas_pool_copy_fied() 
 
   write(*,*)'====> copy of mpas_field done'
  
end subroutine copy

! ------------------------------------------------------------------------------

subroutine self_add(self,rhs)

   implicit none
   type(mpas_field), intent(inout) :: self
   type(mpas_field), intent(in)    :: rhs
   character(len=StrKIND) :: kind_op

   kind_op = 'add'
   call da_operator(trim(kind_op), self % subFields, rhs % subFields)

end subroutine self_add

! ------------------------------------------------------------------------------

subroutine self_schur(self,rhs)

   implicit none
   type(mpas_field), intent(inout) :: self
   type(mpas_field), intent(in)    :: rhs
   character(len=StrKIND) :: kind_op

   kind_op = 'schur'
   call da_operator(trim(kind_op), self % subFields, rhs % subFields)

end subroutine self_schur

! ------------------------------------------------------------------------------

subroutine self_sub(self,rhs)

   implicit none
   type(mpas_field), intent(inout) :: self
   type(mpas_field), intent(in)    :: rhs
   character(len=StrKIND) :: kind_op

   kind_op = 'sub'
   call da_operator(trim(kind_op), self % subFields, rhs % subFields)

end subroutine self_sub

! ------------------------------------------------------------------------------

subroutine self_mul(self,zz)

   implicit none
   type(mpas_field),     intent(inout) :: self
   real(kind=kind_real), intent(in)    :: zz

   call da_self_mult(self % subFields, zz)

end subroutine self_mul

! ------------------------------------------------------------------------------

subroutine axpy(self,zz,rhs)

   implicit none
   type(mpas_field),     intent(inout) :: self
   real(kind=kind_real), intent(in)    :: zz
   type(mpas_field),     intent(in)    :: rhs

   call da_axpy(self % subFields, rhs % subFields, zz)

end subroutine axpy

! ------------------------------------------------------------------------------

subroutine dot_prod(fld1,fld2,zprod)

   implicit none
   type(mpas_field),     intent(in)    :: fld1, fld2
   real(kind=kind_real), intent(inout) :: zprod

   call da_dot_product(fld1 % subFields, fld2 % subFields, fld1 % geom % domain % dminfo, zprod)

end subroutine dot_prod

! ------------------------------------------------------------------------------
!> add increment to state
!!
!! \details **add_incr()** adds "increment" to "state", such as
!!          state (containing analysis) = state (containing guess) + increment
!!          Here, we also update "theta", "rho", and "u" (edge-normal wind), which are
!!          close to MPAS prognostic variable.
!!          While conversion to "theta" and "rho" uses full state variables,
!!          conversion to "u" from cell center winds uses their increment to reduce 
!!          the smoothing effect.
!!
subroutine add_incr(self,rhs)

   implicit none
   type(mpas_field), intent(inout) :: self !< state
   type(mpas_field), intent(in)    :: rhs  !< increment
   character(len=StrKIND) :: kind_op

   type (mpas_pool_type), pointer :: state, diag, mesh
   type (field2DReal), pointer :: field2d_t, field2d_p, field2d_qv, field2d_uRz, field2d_uRm, &
                                  field2d_th, field2d_rho, field2d_u, field2d_u_inc

   ! GD: I don''t see any difference than for self_add other than subFields can contain
   ! different variables than mpas_field and the resolution of incr can be different. 

   if (self%geom%nCells==rhs%geom%nCells .and. self%geom%nVertLevels==rhs%geom%nVertLevels) then
      !NOTE: first, get full state of "subFields" variables
      kind_op = 'add'
      call da_operator(trim(kind_op), self % subFields, rhs % subFields)

      !NOTE: second, also update variables which are closely related to MPAS prognostic vars.
      !  update theta from temperature and pressure
      !  update rho   from temperature, pressure, and index_qv
      call mpas_pool_get_field(self % subFields,            'temperature', field2d_t)
      call mpas_pool_get_field(self % subFields,               'pressure', field2d_p)
      call mpas_pool_get_field(self % subFields,               'index_qv', field2d_qv)
      call mpas_pool_get_field(self % subFields,      'uReconstructZonal', field2d_uRz)
      call mpas_pool_get_field(self % subFields, 'uReconstructMeridional', field2d_uRm)
      call mpas_pool_get_field(self % auxFields,                  'theta', field2d_th)
      call mpas_pool_get_field(self % auxFields,                    'rho', field2d_rho)

      field2d_th % array(:,:) = field2d_t % array(:,:) * &
                 ( 100000.0_kind_real / field2d_p % array(:,:) ) ** ( rgas / cp )
      write(*,*) 'add_inc: theta min/max = ', minval(field2d_th % array), maxval(field2d_th % array)
      field2d_rho % array(:,:) = field2d_p % array(:,:) /  ( rgas * field2d_t % array(:,:) * &
                 ( 1.0_kind_real + (rv/rgas - 1.0_kind_real) * field2d_qv % array(:,:) ) )
      write(*,*) 'add_inc: rho min/max = ', minval(field2d_rho % array), maxval(field2d_rho % array)

      !  update u     from uReconstructZonal and uReconstructMeridional "incrementally"
      call mpas_pool_get_field(self % auxFields,                      'u', field2d_u)
      call mpas_pool_get_field( rhs % subFields,      'uReconstructZonal', field2d_uRz)
      call mpas_pool_get_field( rhs % subFields, 'uReconstructMeridional', field2d_uRm)
      call mpas_pool_get_field( rhs % auxFields,                      'u', field2d_u_inc)
      write(*,*) 'add_inc: u_inc min/max = ', minval(field2d_uRz % array), maxval(field2d_uRz % array)
      write(*,*) 'add_inc: v_inc min/max = ', minval(field2d_uRm % array), maxval(field2d_uRm % array)

      call uv_cell_to_edges(self % geom % domain, field2d_uRz, field2d_uRm, field2d_u_inc, &
                 self%geom%latCell, self%geom%lonCell, self%geom%nCellsSolve, &
                 self%geom%edgeNormalVectors, self%geom%nEdgesOnCell, self%geom%edgesOnCell, &
                 self%geom%nVertLevels)
      write(*,*) 'add_inc: u_guess min/max = ', minval(field2d_u % array), maxval(field2d_u % array)
      write(*,*) 'add_inc: u_inc min/max = ', minval(field2d_u_inc % array), maxval(field2d_u_inc % array)
      field2d_u % array(:,:) = field2d_u % array(:,:) + field2d_u_inc % array(:,:)
      write(*,*) 'add_inc: u_analy min/max = ', minval(field2d_u % array), maxval(field2d_u % array)

      ! TODO: DO we need HALO exchange here or in ModelMPAS::initialize for model integration?

   else
      call abor1_ftn("mpas_fields:add_incr: dimension mismatch")
   endif

   return

end subroutine add_incr

! ------------------------------------------------------------------------------

subroutine diff_incr(lhs,x1,x2)

   implicit none
   type(mpas_field), intent(inout) :: lhs
   type(mpas_field), intent(in)    :: x1
   type(mpas_field), intent(in)    :: x2
   character(len=StrKIND) :: kind_op

   call zeros(lhs)
   if (x1%geom%nCells==x2%geom%nCells .and. x1%geom%nVertLevels==x2%geom%nVertLevels) then
     if (lhs%geom%nCells==x1%geom%nCells .and. lhs%geom%nVertLevels==x1%geom%nVertLevels) then
        kind_op = 'sub'
        call da_operator(trim(kind_op), lhs % subFields, x1 % subFields, x2 % subFields)
     else
       call abor1_ftn("mpas_fields:diff_incr: dimension mismatch between the two variables.")
     endif
   else
     call abor1_ftn("mpas_fields:diff_incr: states not at same resolution")
   endif

   return

end subroutine diff_incr

! ------------------------------------------------------------------------------

subroutine change_resol(fld,rhs)

   implicit none
   type(mpas_field), intent(inout) :: fld
   type(mpas_field), intent(in)    :: rhs

   ! FIXME: We just copy rhs to fld for now. Need an actual interpolation routine later. (SH)
   if (fld%geom%nCells == rhs%geom%nCells .and.  fld%geom%nVertLevels == rhs%geom%nVertLevels) then
     call copy(fld, rhs)
   else
     write(0,*) fld%geom%nCells, rhs%geom%nCells, fld%geom%nVertLevels, rhs%geom%nVertLevels
     call abor1_ftn("mpas_fields:field_resol: dimension mismatch")
   endif

end subroutine change_resol

! ------------------------------------------------------------------------------
!> Analytic Initialization for the MPAS Model
!!
!! \details **analytic_IC()** initializes the MPAS Field and State objects using one of
!! several alternative idealized analytic models.  This is intended to facilitate testing by
!! eliminating the need to read in the initial state from a file and by providing exact expressions
!! to test interpolations.  This function is activated by setting the "analytic_init" field in the
!! "initial" or "StateFile" section of the configuration file.
!!
!! Initialization options that begin with "dcmip" refer to tests defined by the multi-institutional
!! 2012 [Dynamical Core Intercomparison Project](https://earthsystealcmcog.org/projects/dcmip-2012)
!! and the associated Summer School, sponsored by NOAA, NSF, DOE, NCAR, and the University of Michigan.
!!
!! Currently implemented options for analytic_init include:
!! * invent-state: Backward compatibility with original analytic init option
!! * dcmip-test-1-1: 3D deformational flow
!! * dcmip-test-1-2: 3D Hadley-like meridional circulation
!! * dcmip-test-3-1: Non-hydrostatic gravity wave
!! * dcmip-test-4-0: Baroclinic instability
!!
!! \author J. Guerrette (adapted from fv3jedi code by M. Miesch)
!! \date July, 2018: Created
!!
subroutine analytic_IC(fld, geom, c_conf, vdate)

  use kinds
!  use iso_c_binding
!  use datetime_mod
!  use fckit_log_module, only : log
  use dcmip_initial_conditions_test_1_2_3, only : test1_advection_deformation, &
       test1_advection_hadley, test3_gravity_wave
  use dcmip_initial_conditions_test_4, only : test4_baroclinic_wave

!  !MPAS Test Cases
!  !JJG: This initialization requires the init_atmospher_core core_type 
!  !      in the MPAS library for OOPS, but currently it is not included
!  use init_atm_core, only: init_atm_core_run!, init_atm_core_finalize (could be used for cleanup...)

  implicit none

  type(mpas_field), intent(inout)     :: fld !< Fields
  type(mpas_geom), target, intent(in) :: geom    !< Geometry 
  type(c_ptr), intent(in)                :: c_conf   !< Configuration
  type(datetime), intent(inout)          :: vdate    !< DateTime

  character(len=30) :: IC
  character(len=20) :: sdate
  character(len=1024) :: buf
  Integer :: jlev,ii
  integer :: ierr = 0 
  real(kind=kind_real) :: rlat, rlon, z
  real(kind=kind_real) :: pk,pe1,pe2,ps
  real(kind=kind_real) :: u0,v0,w0,t0,phis0,ps0,rho0,hum0,q1,q2,q3,q4

  real(kind=kind_real)             :: DTdummy = 900.0
  logical, allocatable             :: grids_on_this_pe(:)
  integer                          :: p_split = 1
  real (kind=kind_real), dimension(:,:), pointer :: &
              u_ptr, v_ptr, temperature_ptr, p_ptr, &
              qv_ptr, qc_ptr, qr_ptr, qi_ptr, qs_ptr, &
              ln_p_ptr
  integer, pointer :: index_qv, index_qc, index_qr, index_qi, index_qs
  type (field3DReal), pointer :: field3d
  type (mpas_pool_type), pointer :: pool_a, pool_b, state
  real(kind=kind_real) :: zhalf

  ! Pointer to geometry component of field object
  fld%geom => geom

  If (config_element_exists(c_conf,"analytic_init")) Then
     IC = Trim(config_get_string(c_conf,len(IC),"analytic_init"))
  Else
     ! This default value is for backward compatibility
     IC = "invent-state"
  EndIf

  WRITE(*,*) "mpas_fields:analytic_init: "//IC

! Conflicts with natural log below
!  call log%warning("mpas_fields:analytic_init: "//IC)
  sdate = config_get_string(c_conf,len(sdate),"date")
  WRITE(buf,*) 'validity date is: '//sdate
  WRITE(*,*) buf
!  call log%info(buf)
  call datetime_set(sdate, vdate)

   ! Need to initialize variables that are used in interpolation/getVals
   ! In "create" and "read" subroutines, subFields and auxFields are 
   ! initialized from geom % domain % blocklist % allFields, zeroed,
   ! reread from file into allFields, then values copied to sub/auxFields
   ! -> must initialize allFields here and copy to sub/auxFields

   call mpas_pool_get_subpool(geom % domain % blocklist % structs, &
                              'state', state)

   pool_a => geom % domain % blocklist % allFields

   !Diagnostic vars (diag pool)
   call mpas_pool_get_array(pool_a, "pressure", p_ptr)
   call mpas_pool_get_array(pool_a, "uReconstructZonal", u_ptr)
   call mpas_pool_get_array(pool_a, "uReconstructMeridional", v_ptr)
   call mpas_pool_get_array(pool_a, "temperature", temperature_ptr)


   !Scalars (state pool)
   call mpas_pool_get_field(pool_a, "scalars", field3d)
   call mpas_pool_get_dimension(state, "index_qv", index_qv)
   if ( index_qv .gt. 0 ) &
      qv_ptr => field3d % array(index_qv,:,:)

   call mpas_pool_get_dimension(state, "index_qc", index_qc)
   if ( index_qc .gt. 0 ) &
      qc_ptr => field3d % array(index_qc,:,:)

   call mpas_pool_get_dimension(state, "index_qr", index_qr)
   if ( index_qr .gt. 0 ) &
      qr_ptr => field3d % array(index_qr,:,:)

   call mpas_pool_get_dimension(state, "index_qi", index_qi)
   if ( index_qi .gt. 0 ) &
      qi_ptr => field3d % array(index_qi,:,:)

   call mpas_pool_get_dimension(state, "index_qs", index_qs)
   if ( index_qs .gt. 0 ) &
      qs_ptr => field3d % array(index_qs,:,:)

  !===================================================================
  int_option: Select Case (IC)

     Case("invent-state")

        call invent_state(fld,c_conf)


!     !TODO: This case requires the init_atmospher_core core_type to be 
!     !      built as part of the MPAS library.
!     Case("mpas_init_case") 
!
!!Would use init_atm_setup_case in MPAS-Release/src/core_init_atmosphere/mpas_init_atm_cases.F
!!mpas_init has already been called at this point from geo_setup
!
!!init_atms_setup_case is normally called from the following set of subroutines:
!!mpas_run => core_run [init_atm_core_run] => init_atm_setup_case => [select from preset cases]
!!Can we bypass the first two somehow?  
!!Would use "config_init_case" in the json file, then check for matching with one of the ideal cases below... (not 7 or 8)
!
!!if ((config_init_case == 1) .or. (config_init_case == 2) .or. (config_init_case == 3)) then
!!   write(0,*) ' Jablonowski and Williamson baroclinic wave test case '
!!   if (config_init_case == 1) write(0,*) ' no initial perturbation '
!!   if (config_init_case == 2) write(0,*) ' initial perturbation included '
!!   if (config_init_case == 3) write(0,*) ' normal-mode perturbation included '
!!else if ((config_init_case == 4) .or. (config_init_case == 5)) then
!!   write(0,*) ' squall line - super cell test case '
!!   if (config_init_case == 4) write(0,*) ' squall line test case'
!!   if (config_init_case == 5) write(0,*) ' supercell test case'
!!      else if (config_init_case == 6 ) then
!!   write(0,*) ' mountain wave test case '
!!else if (config_init_case == 7 ) then
!!   write(0,*) ' real-data GFS test case '
!!else if (config_init_case == 8 ) then
!!   write(0,*) 'real-data surface (SST) update test case '
!
!       ierr = init_atm_core_run(geom % domain)
!       if ( ierr .ne. 0  ) then
!          call abor1_ftn("mpas_fields: init_atm_core_run failed")
!       end if

     Case ("dcmip-test-1-1")
        do ii = 1, geom%nCellsSolve
           rlat = geom%latCell(ii)
           rlon = geom%lonCell(ii)

           ! Now loop over all levels
           do jlev = 1, geom%nVertLevels

              zhalf = 0.5_kind_real * (geom%zgrid(jlev,ii) + geom%zgrid(jlev+1,ii))
              Call test1_advection_deformation(rlon,rlat,pk,zhalf,1,u0,v0,w0,t0,&
                                               phis0,ps0,rho0,hum0,q1,q2,q3,q4)
              p_ptr(jlev,ii) = pk

              u_ptr(jlev,ii) = u0 !MMiesch: ATTN Not going to necessary keep a-grid winds, u can be either a-
              v_ptr(jlev,ii) = v0 ! or staggered-grid so this needs to be generic. You cannot drive the model 
                                  ! with A grid winds
              if (index_qv.gt.0) qv_ptr(jlev,ii) = hum0 !set to zero for this test
              if (index_qc.gt.0) qc_ptr(jlev,ii) = q1
              if (index_qi.gt.0) qi_ptr(jlev,ii) = q2
              if (index_qr.gt.0) qr_ptr(jlev,ii) = q3
              if (index_qs.gt.0) qs_ptr(jlev,ii) = q4 

              temperature_ptr(jlev,ii) = t0
           enddo
        enddo

     Case ("dcmip-test-1-2")

        do ii = 1, geom%nCellsSolve
           rlat = geom%latCell(ii)
           rlon = geom%lonCell(ii)

           ! Now loop over all levels
           do jlev = 1, geom%nVertLevels

              zhalf = 0.5_kind_real * (geom%zgrid(jlev,ii) + geom%zgrid(jlev+1,ii))
              Call test1_advection_hadley(rlon,rlat,pk,zhalf,1,u0,v0,w0,&
                                          t0,phis0,ps0,rho0,hum0,q1)
              p_ptr(jlev,ii) = pk

              u_ptr(jlev,ii) = u0 !MMiesch: ATTN Not going to necessary keep a-grid winds, u can be either a-
              v_ptr(jlev,ii) = v0 ! or staggered-grid so this needs to be generic. You cannot drive the model 
                                  ! with A grid winds
              if (index_qv.gt.0) qv_ptr(jlev,ii) = hum0 !set to zero for this test
              if (index_qc.gt.0) qc_ptr(jlev,ii) = q1

              if (index_qi.gt.0) qi_ptr(jlev,ii) = 0._kind_real
              if (index_qr.gt.0) qr_ptr(jlev,ii) = 0._kind_real
              if (index_qs.gt.0) qs_ptr(jlev,ii) = 0._kind_real

              temperature_ptr(jlev,ii) = t0
           enddo
        enddo

     Case ("dcmip-test-3-1")

        do ii = 1, geom%nCellsSolve
           rlat = geom%latCell(ii)
           rlon = geom%lonCell(ii)

           ! Now loop over all levels
           do jlev = 1, geom%nVertLevels

              zhalf = 0.5_kind_real * (geom%zgrid(jlev,ii) + geom%zgrid(jlev+1,ii))
              Call test3_gravity_wave(rlon,rlat,pk,zhalf,1,u0,v0,w0,&
                                      t0,phis0,ps0,rho0,hum0)

              p_ptr(jlev,ii) = pk
              u_ptr(jlev,ii) = u0 !MMiesch: ATTN Not going to necessary keep a-grid winds, u can be either a-
              v_ptr(jlev,ii) = v0 ! or staggered-grid so this needs to be generic. You cannot drive the model 
                                  ! with A grid winds

              if (index_qv.gt.0) qv_ptr (jlev,ii) = hum0 !set to zero for this test
              if (index_qc.gt.0) qc_ptr(jlev,ii) = 0._kind_real
              if (index_qi.gt.0) qi_ptr(jlev,ii) = 0._kind_real
              if (index_qr.gt.0) qr_ptr(jlev,ii) = 0._kind_real
              if (index_qs.gt.0) qs_ptr(jlev,ii) = 0._kind_real

              temperature_ptr(jlev,ii) = t0
           enddo
        enddo

     Case ("dcmip-test-4-0")

        do ii = 1, geom%nCellsSolve
           rlat = geom%latCell(ii)
           rlon = geom%lonCell(ii)

           ! Now loop over all levels
           do jlev = 1, geom%nVertLevels

              zhalf = 0.5_kind_real * (geom%zgrid(jlev,ii) + geom%zgrid(jlev+1,ii))
              Call test4_baroclinic_wave(0,1.0_kind_real,rlon,rlat,pk,zhalf,1,u0,v0,w0,&
                                      t0,phis0,ps0,rho0,hum0,q1,q2)

              p_ptr(jlev,ii) = pk

              u_ptr(jlev,ii) = u0 !MMiesch: ATTN Not going to necessary keep a-grid winds, u can be either a-

              v_ptr(jlev,ii) = v0 ! or staggered-grid so this needs to be generic. You cannot drive the model 
                                  ! with A grid winds

              if (index_qv.gt.0) qv_ptr(jlev,ii) = hum0 !set to zero for this test
              if (index_qc.gt.0) qc_ptr(jlev,ii) = 0._kind_real
              if (index_qi.gt.0) qi_ptr(jlev,ii) = 0._kind_real
              if (index_qr.gt.0) qr_ptr(jlev,ii) = 0._kind_real
              if (index_qs.gt.0) qs_ptr(jlev,ii) = 0._kind_real

              temperature_ptr(jlev,ii) = t0
           enddo
        enddo

     Case Default

        call invent_state(fld,c_conf)

     End Select int_option

     call da_copy_all2sub_fields(fld % geom % domain, fld % subFields) 
     call da_copy_all2sub_fields(fld % geom % domain, fld % auxFields) 

   write(*,*)'==> end mpas_fields:analytic_init'

end subroutine analytic_IC



! ------------------------------------------------------------------------------
subroutine invent_state(fld,config)

   use kinds

   implicit none

   type(mpas_field), intent(inout) :: fld    !< Model fields
   type(c_ptr), intent(in)         :: config  !< Configuration structure
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   integer :: jlev,ii
   type (mpas_pool_type), pointer :: pool_a, state
   type (field3DReal), pointer :: field3d
   integer, pointer :: index_qv

   !- read/interp.

   call mpas_pool_get_subpool(fld % geom % domain % blocklist % structs, &
                              'state', state)
   pool_a => fld % geom % domain % blocklist % allFields

   !Diagnostic vars (diag pool)
   !u
   call mpas_pool_get_array(pool_a, "uReconstructZonal", r2d_ptr_a)
   do jlev = 1,fld % geom % nVertLevels
      do ii = 1, fld % geom % nCellsSolve
         r2d_ptr_a(jlev,ii) = cos(0.25*fld % geom % lonEdge(ii)) + cos(0.25*fld % geom % latEdge(ii))
      enddo
   enddo

   !v
   call mpas_pool_get_array(pool_a, "uReconstructMeridional", r2d_ptr_a)
   do jlev = 1,fld % geom % nVertLevels
      do ii = 1, fld % geom % nCellsSolve
         r2d_ptr_a(jlev,ii) = 1.0_kind_real
      enddo
   enddo

   !temperature
   call mpas_pool_get_array(pool_a, "temperature", r2d_ptr_a)
   do jlev = 1,fld % geom % nVertLevels
      do ii = 1, fld % geom % nCellsSolve
         r2d_ptr_a(jlev,ii) = cos(0.25*fld % geom % lonCell(ii)) + cos(0.25*fld % geom % latCell(ii))
      enddo
   enddo

   !pressure
   call mpas_pool_get_array(pool_a, "pressure", r2d_ptr_a)
   do jlev = 1,fld % geom % nVertLevels
      do ii = 1, fld % geom % nCellsSolve
         r2d_ptr_a(jlev,ii) = real(jlev,kind_real)
      enddo
   enddo

   !Scalars (state pool)
   !qv
   call mpas_pool_get_field(pool_a, 'scalars', field3d)
   call mpas_pool_get_dimension(state, 'index_qv', index_qv)
   if ( index_qv .gt. 0 ) then
      r2d_ptr_a => field3d % array(index_qv,:,:)
      do jlev = 1,fld % geom % nVertLevels
         do ii = 1, fld % geom % nCellsSolve
            r2d_ptr_a(jlev,ii) = 0.0_kind_real
         enddo
      enddo
   end if

return
end subroutine invent_state
! -----------------------------------------------------------------------------------------------------------

subroutine read_file(fld, c_conf, vdate)

   implicit none
   type(mpas_field), intent(inout) :: fld      !< Fields
   type(c_ptr),      intent(in)    :: c_conf   !< Configuration
   type(datetime),   intent(inout) :: vdate    !< DateTime
   character(len=20)       :: sdate
   type (MPAS_Time_type)   :: local_time
   character (len=StrKIND) :: dateTimeString, streamID, time_string, filename, temp_filename
   integer                 :: ierr = 0

   type (mpas_pool_type), pointer :: state, diag, mesh
   type (field2DReal), pointer :: field2d, field2d_b, field2d_c

   write(*,*)'==> read fields'
   sdate = config_get_string(c_conf,len(sdate),"date")
   call datetime_set(sdate, vdate)

   temp_filename = config_get_string(c_conf,len(temp_filename),&
                      "filename")
   write(*,*)'Reading ',trim(temp_filename)
   !temp_filename = 'restart.$Y-$M-$D_$h.$m.$s.nc'
   ! GD look at oops/src/util/datetime_mod.F90
   ! we probably need to extract from vdate a string to enforce the reading ..
   ! and then can be like this ....
   ! TODO: we can get streamID from json
   !streamID = 'restart'
   !streamID = 'input'
   streamID = 'output'
   !streamID = 'da'
   ierr = 0
   fld % manager => fld % geom % domain % streamManager
   dateTimeString = '$Y-$M-$D_$h:$m:$s'
   call cvt_oopsmpas_date(sdate,dateTimeString,1)
   write(*,*)'dateTimeString: ',trim(dateTimeString)
   call mpas_set_time(local_time, dateTimeString=dateTimeString, ierr=ierr)
   !call mpas_set_clock_time(fld % clock, local_time, MPAS_START_TIME)
   call mpas_set_clock_time(fld % geom % domain % clock, local_time, MPAS_START_TIME)
   call mpas_expand_string(dateTimeString, -1, temp_filename, filename)
   call MPAS_stream_mgr_set_property(fld % manager, streamID, MPAS_STREAM_PROPERTY_FILENAME, filename)
   write(*,*)'Reading ',trim(filename)
   call MPAS_stream_mgr_read(fld % manager, streamID=streamID, &
                           & when=dateTimeString, rightNow=.True., ierr=ierr)
   if ( ierr .ne. 0  ) then
      call abor1_ftn('MPAS_stream_mgr_read failed ierr=',ierr)
   end if
   !==TODO: Speific part when reading parameterEst. for BUMP.
   !      : They write/read a list of variables directly.
   If (config_element_exists(c_conf,"no_transf")) Then
      ierr = config_get_int(c_conf,"no_transf")
      if(ierr .eq. 1) then
        call da_copy_all2sub_fields(fld % geom % domain, fld % subFields) 
        call da_copy_all2sub_fields(fld % geom % domain, fld % auxFields) 
        return
      endif
   endif
   !--TODO: BJJ test. Do I need to "re-calculate"/"update" diagnostic variables ?
   !call update_mpas_field(fld % geom % domain, fld % auxFields) !--> this will construct "pressure" from pressure_base & pressure_p
   call mpas_pool_get_subpool(fld % geom % domain % blocklist % structs,'state',state)
   call mpas_pool_get_subpool(fld % geom % domain % blocklist % structs, 'diag', diag)
   call mpas_pool_get_subpool(fld % geom % domain % blocklist % structs, 'mesh', mesh)
   call atm_compute_output_diagnostics(state, 1, diag, mesh)

   call da_copy_all2sub_fields(fld % geom % domain, fld % subFields) 
   call da_copy_all2sub_fields(fld % geom % domain, fld % auxFields) 
   !TODO- special case: read theta, pressure and convert to temperature
   !NOTE: This formula is somewhat different with MPAS one's (in physics, they use "exner") 
   !    : If T diagnostic is added in, for example, subroutine atm_compute_output_diagnostics 6 lines above,
   !    : we need to include "exner" in stream_list.for.reading
   call mpas_pool_get_field(fld % auxFields, 'theta', field2d_b)
   call mpas_pool_get_field(fld % subFields, 'pressure', field2d_c)
   call mpas_pool_get_field(fld % subFields, 'temperature', field2d)
   field2d % array(:,:) = field2d_b % array(:,:) / &
             ( 100000.0_kind_real / field2d_c % array(:,:) ) ** ( rgas / cp )

end subroutine read_file

! ------------------------------------------------------------------------------

subroutine write_file(fld, c_conf, vdate)

use duration_mod
   implicit none
   type(mpas_field), intent(inout) :: fld    !< Fields
   type(c_ptr),      intent(in)    :: c_conf !< Configuration
   type(datetime),   intent(inout) :: vdate  !< DateTime
   character(len=20)       :: validitydate
   integer                 :: ierr
   type (MPAS_Time_type)   :: fld_time, write_time
   character (len=StrKIND) :: dateTimeString, dateTimeString2, streamID, time_string, filename, temp_filename

   call datetime_to_string(vdate, validitydate)
   write(*,*)'==> write fields at ',trim(validitydate)
   temp_filename = config_get_string(c_conf,len(temp_filename)&
                      ,"filename")
   write(*,*)'==> writing ',trim(temp_filename)
   !temp_filename = 'restart.$Y-$M-$D_$h.$m.$s.nc'
   ! GD look at oops/src/util/datetime_mod.F90
   ! we probably need to extract from vdate a string to enforce the reading ..
   ! and then can be like this ....
   dateTimeString = '$Y-$M-$D_$h:$m:$s'
   call cvt_oopsmpas_date(validitydate,dateTimeString,-1)
   ierr = 0
   call mpas_set_time(write_time, dateTimeString=dateTimeString, ierr=ierr)
   fld_time = mpas_get_clock_time(fld % clock, MPAS_NOW, ierr)
   call mpas_get_time(fld_time, dateTimeString=dateTimeString2, ierr=ierr)
   if ( fld_time .NE. write_time ) then
      write(*,*)'write_time,fld_time: ',trim(dateTimeString),trim(dateTimeString2)
      call abor1_ftn('Different times MPAS_stream_mgr_write failed ')
   end if
   call mpas_expand_string(dateTimeString, -1, trim(temp_filename), filename)
   fld % manager => fld % geom % domain % streamManager
   ! TODO: we can get streamID from json
   !streamID = 'restart'
   streamID = 'output'
   !streamID = 'da'
   call MPAS_stream_mgr_set_property(fld % manager, streamID, MPAS_STREAM_PROPERTY_FILENAME, filename)
   call da_copy_sub2all_fields(fld % geom % domain, fld % subFields)
   call da_copy_sub2all_fields(fld % geom % domain, fld % auxFields)
   write(*,*)'writing ',trim(filename)
   call mpas_stream_mgr_write(fld % geom % domain % streamManager, streamID=streamID, forceWriteNow=.true., ierr=ierr)
   if ( ierr .ne. 0  ) then
     call abor1_ftn('MPAS_stream_mgr_write failed ierr=',ierr)
   end if

end subroutine write_file

! ------------------------------------------------------------------------------

subroutine gpnorm(fld, nf, pstat)

   implicit none
   type(mpas_field),     intent(in)  :: fld
   integer,              intent(in)  :: nf
   real(kind=kind_real), intent(out) :: pstat(3, nf)

   call da_gpnorm(fld % subFields, fld % geom % domain % dminfo, fld%nf, pstat) 

end subroutine gpnorm

! ------------------------------------------------------------------------------

subroutine fldrms(fld, prms)

   implicit none
   type(mpas_field),     intent(in)  :: fld
   real(kind=kind_real), intent(out) :: prms

   call da_fldrms(fld % subFields, fld % geom % domain % dminfo, prms)

end subroutine fldrms

! ------------------------------------------------------------------------------

subroutine dirac(self, c_conf)

use iso_c_binding
use mpas_pool_routines

   implicit none
   type(mpas_field), intent(inout) :: self
   type(c_ptr),      intent(in)    :: c_conf   !< Configuration
   integer                :: ndir, idir, ildir, ndirlocal
   character(len=3)       :: idirchar
   character(len=StrKIND) :: dirvar
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a
   integer :: nearestCell
   integer, allocatable, dimension(:) :: dirOwned, dirOwnedGlobal
   real (kind=kind_real), allocatable, dimension(:) :: dirLats
   real (kind=kind_real), allocatable, dimension(:) :: dirLons
   integer, allocatable, dimension(:) :: dirCells

   ! Get number and positions of Diracs
   ndir = config_get_int(c_conf,"ndir")

   allocate( dirOwned(ndir) )
   allocate( dirLats(ndir) )
   allocate( dirLons(ndir) )
   allocate( dirCells(ndir) )

   do idir=1,ndir
      write(idirchar,'(i3)') idir
      dirLats(idir) = config_get_real(c_conf,"dirLats("//trim(adjustl(idirchar))//")")
      dirLons(idir) = config_get_real(c_conf,"dirLons("//trim(adjustl(idirchar))//")")
   end do
   ildir = config_get_int(c_conf,"ildir")
   dirvar = config_get_string(c_conf,len(dirvar),"dirvar")

   !Test if dir is owned and find the nearest local cell
   ! (repurposed from MPAS-Release/src/core_atmosphere/diagnostics/soundings.F)

   ndirlocal = 0
   do idir=1,ndir
      nearestCell = self % geom % nCellsSolve
      nearestCell = nearest_cell( (dirLats(idir) * deg2rad), &
                                  (dirLons(idir) * deg2rad), &
                                  nearestCell, self % geom % nCells, self % geom % maxEdges, &
                                  self % geom % nEdgesOnCell, self % geom % cellsOnCell, &
                                  self % geom % latCell, self % geom % lonCell )

      if (nearestCell <= self % geom % nCellsSolve) then
          dirOwned(idir) = 1
          dirCells(idir) = nearestCell
          ndirlocal = ndirlocal + 1
      else
          dirOwned(idir) = 0
          dirCells(idir) = self % geom % nCells + 1
      end if
   end do

   write(*,*) ' This processor owns ',ndirlocal, &
              ' dirac forcing locations'

   ! Check
   if (ndir<1) call abor1_ftn("mpas_fields:dirac non-positive ndir")

   allocate( dirOwnedGlobal(ndir) )
   call mpas_dmpar_max_int_array( self % geom % domain % dminfo, ndir, dirOwned, dirOwnedGlobal)
   if ( any(dirOwnedGlobal.lt.1) ) then
         call abor1_ftn("mpas_fields:dirac invalid Lat/Lon")
   end if

   call mpas_dmpar_sum_int_array( self % geom % domain % dminfo, ndir, dirOwned, dirOwnedGlobal)
   if ( any(dirOwnedGlobal.gt.1) ) then
         call abor1_ftn("mpas_fields:duplicated dirac on >1 processors")
   end if
   deallocate( dirOwnedGlobal )

  if ((ildir < 1) .or. (ildir > self % geom % nVertLevels)) then
      call abor1_ftn("mpas_fields:dirac invalid ildir")
   endif

   ! Setup Diracs
   call zeros(self)

   call mpas_pool_begin_iteration(self % subFields)

   do while ( mpas_pool_get_next_member(self % subFields, poolItr) )
        write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , poolItr % memberName
        ! Pools may in general contain dimensions, namelist options, fields, or other pools,
        ! so we select only those members of the pool that are fields
        if (poolItr % memberType == MPAS_POOL_FIELD) then
        ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
        if (poolItr % dataType == MPAS_POOL_REAL) then
           ! Depending on the dimensionality of the field, we need to set pointers of
           ! the correct type
           if (poolItr % nDims == 1) then
              write(*,*)'Not implemented yet'
          else if (poolItr % nDims == 2) then
              call mpas_pool_get_array(self % subFields, trim(poolItr % memberName), r2d_ptr_a)
              if( trim(dirvar) .eq. trim(poolItr % memberName) ) then
                ndirlocal = 0
                do idir=1, ndir
                   if ( dirOwned(idir).eq.1 ) then
                      r2d_ptr_a( ildir, dirCells(idir) ) = 1.0_kind_real
                      ndirlocal = ndirlocal + 1
                   end if
                end do
                write(*,*) ' Dirac is set in ',ndirlocal,'locations for',trim(poolItr % memberName)
              end if
           else if (poolItr % nDims == 3) then
              write(*,*)'Not implemented yet'
           end if
        end if
        end if
   end do

   deallocate( dirOwned )
   deallocate( dirLats )
   deallocate( dirLons )
   deallocate( dirCells )

end subroutine dirac

!!TODO: Alternatively could make the function nearest_cell public in MPAS
!!      and then define an interface to it within this module (cleaner)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finds the MPAS grid cell nearest to (target_lat, target_lon)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function nearest_cell(target_lat, target_lon, start_cell, nCells, maxEdges, &
                              nEdgesOnCell, cellsOnCell, latCell, lonCell)

   implicit none

   real (kind=kind_real), intent(in) :: target_lat, target_lon
   integer, intent(in) :: start_cell
   integer, intent(in) :: nCells, maxEdges
   integer, dimension(nCells), intent(in) :: nEdgesOnCell
   integer, dimension(maxEdges,nCells), intent(in) :: cellsOnCell
   real (kind=kind_real), dimension(nCells), intent(in) :: latCell, lonCell

   integer :: i
   integer :: iCell
   integer :: current_cell
   real (kind=kind_real) :: current_distance, d
   real (kind=kind_real) :: nearest_distance

   nearest_cell = start_cell
   current_cell = -1

   do while (nearest_cell /= current_cell)
      current_cell = nearest_cell
      current_distance = sphere_distance(latCell(current_cell), lonCell(current_cell), target_lat, &
                                         target_lon, 1.0_kind_real)
      nearest_cell = current_cell
      nearest_distance = current_distance
      do i = 1, nEdgesOnCell(current_cell)
         iCell = cellsOnCell(i,current_cell)
         if (iCell <= nCells) then
            d = sphere_distance(latCell(iCell), lonCell(iCell), target_lat, target_lon, 1.0_kind_real)
            if (d < nearest_distance) then
               nearest_cell = iCell
               nearest_distance = d
            end if
         end if
      end do
   end do
end function nearest_cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) 
!    on a sphere with given radius.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (kind=kind_real) function sphere_distance(lat1, lon1, lat2, lon2, radius)

   implicit none

   real (kind=kind_real), intent(in) :: lat1, lon1, lat2, lon2, radius
   real (kind=kind_real) :: arg1

   arg1 = sqrt( sin(0.5*(lat2-lat1))**2 + cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
   sphere_distance = 2.0 * radius * asin(arg1)

end function sphere_distance

! ------------------------------------------------------------------------------

subroutine getvalues(fld, locs, vars, gom, traj)

   use type_bump, only: bump_type
   use mpas2ufo_vars_mod !, only: usgs_to_crtm_mw, wrf_to_crtm_soil

   implicit none
   type(mpas_field),                        intent(in)    :: fld
   type(ioda_locs),                         intent(in)    :: locs
   type(ufo_vars),                          intent(in)    :: vars
   type(ufo_geovals),                       intent(inout) :: gom
   type(mpas_getvaltraj), optional, target, intent(inout) :: traj
   
   character(len=*), parameter :: myname = 'getvalues'

   type(bump_type), target  :: bump
   type(bump_type), pointer :: pbump
   logical,         target  :: bump_alloc
   logical,         pointer :: pbumpa
   
   integer :: ii, jj, ji, jvar, jlev, ngrid, nobs, ivar
   real(kind=kind_real), allocatable :: mod_field(:,:), mod_field_ext(:,:)
   real(kind=kind_real), allocatable :: obs_field(:,:)
   real(kind=kind_real), allocatable :: tmp_field(:,:)  !< for wspeed/wdir
   
   type (mpas_pool_type), pointer :: pool_ufo  !< pool with ufo variables
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   integer, dimension(:), pointer :: i1d_ptr_a, i1d_ptr_b
   integer, allocatable  :: index_nn(:)
   real (kind=kind_real), allocatable :: weight_nn(:)
   type (mpas_pool_type), pointer :: pool_tmp  !< temporary pool for setting trajectory
   type (field2DReal), pointer :: field2d => null()     !< for setting trajectory
   type (field2DReal), pointer :: field2d_src => null() !< for setting trajectory

   real(kind=kind_real) :: wdir           !< for wind direction
   integer :: ivarw, ivarl, ivari, ivars  !< for sfc fraction indices


   ! Get grid dimensions and checks
   ! ------------------------------
   ngrid = fld%geom%nCellsSolve
   nobs = locs%nlocs 
   write(*,*)'interp: ngrid, nobs = : ',ngrid, nobs
   call interp_checks("nl", fld, locs, vars, gom)

   ! Allocate and set trajectory for obsop.
   if (present(traj)) then

     pbump => traj % bump

     if (.not. traj%lalloc) then

       traj%ngrid = ngrid
       traj%nobs = nobs

       call mpas_pool_create_pool( pool_tmp, traj % nsize )
       call mpas_pool_get_field(fld % subFields, 'temperature', field2d_src)
       call mpas_pool_add_field(pool_tmp, 'temperature', field2d_src)
       call mpas_pool_get_field(fld % subFields, 'index_qv', field2d_src)
       call mpas_pool_add_field(pool_tmp, 'index_qv', field2d_src)

       call mpas_pool_clone_pool(pool_tmp, traj % pool_traj)

       pbumpa => traj%lalloc

    endif

  else

    pbump => bump
    bump_alloc = .false.
    pbumpa => bump_alloc

  endif
  
  if (.not. pbumpa) then 
    ! Calculate interpolation weight using BUMP
    ! ------------------------------------------
    write(*,*)'call initialize_interp(...)'
    call initialize_interp(fld%geom, locs, pbump)
    pbumpa = .true.
    write(*,*)'interp: after initialize_interp'
  endif
   
   !Make sure the return values are allocated and set
   !-------------------------------------------------
! BJJ: move allocation inside mpas_pool iteration
!   do jvar=1,vars%nv
!      if( .not. allocated(gom%geovals(jvar)%vals) )then
!         gom%geovals(jvar)%nval = fld%geom%nVertLevels
!         allocate( gom%geovals(jvar)%vals(fld%geom%nVertLevels,nobs) )
!         write(*,*) ' gom%geovals(n)%vals allocated'
!      endif
!   enddo
   gom%linit = .true.
   

   !Create Buffer for interpolated values
   !--------------------------------------
   allocate(mod_field(ngrid,1))
   allocate(obs_field(nobs,1))
   
   !Interpolate fields to obs locations using pre-calculated weights
   !----------------------------------------------------------------
   write(0,*)'interp: vars%nv       : ',vars%nv
   write(0,*)'interp: vars%fldnames : ',vars%fldnames
   

   !------- need some table matching UFO_Vars & related MPAS_Vars
   !------- for example, Tv @ UFO may require Theta, Pressure, Qv.
   !-------                               or  Theta_m, exner_base, Pressure_base, Scalar(&index_qv)
   call convert_mpas_field2ufo(fld % geom, fld % subFields, fld % auxFields, pool_ufo, vars % fldnames, vars % nv) !--pool_ufo is new pool with ufo_vars

   call mpas_pool_begin_iteration(pool_ufo)

   do while ( mpas_pool_get_next_member(pool_ufo, poolItr) )
        write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)
        if (poolItr % memberType == MPAS_POOL_FIELD) then

        if (poolItr % dataType == MPAS_POOL_INTEGER) then
           if (poolItr % nDims == 1) then
              call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), i1d_ptr_a)
              ivar = ufo_vars_getindex(vars, trim(poolItr % memberName) )
              write(*,*) "interp: var, ufo_var_index = ",trim(poolItr % memberName), ivar
              if( .not. allocated(gom%geovals(ivar)%vals) )then
                 gom%geovals(ivar)%nval = 1
                 allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
                 write(*,*) ' gom%geovals(n)%vals allocated'
              endif
              mod_field(:,1) = real( i1d_ptr_a(1:ngrid), kind_real)
              !write(*,*) 'MIN/MAX of ',trim(poolItr % memberName),minval(i1d_ptr_a),maxval(i1d_ptr_a)
              call pbump%apply_obsop(mod_field,obs_field)
              gom%geovals(ivar)%vals(1,:) = obs_field(:,1)
              !write(*,*) 'MIN/MAX of ',trim(poolItr % memberName),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
           end if
        end if

        if (poolItr % dataType == MPAS_POOL_REAL) then
           if (poolItr % nDims == 1) then
              call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r1d_ptr_a)
              ivar = ufo_vars_getindex(vars, trim(poolItr % memberName) )
              write(*,*) "interp: var, ufo_var_index = ",trim(poolItr % memberName), ivar
              if( .not. allocated(gom%geovals(ivar)%vals) )then
                 gom%geovals(ivar)%nval = 1
                 allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
                 write(*,*) ' gom%geovals(n)%vals allocated'
              endif
              mod_field(:,1) = r1d_ptr_a(1:ngrid)
              write(*,*) 'MIN/MAX of ',trim(poolItr % memberName),minval(r1d_ptr_a),maxval(r1d_ptr_a)
              call pbump%apply_obsop(mod_field,obs_field)
              gom%geovals(ivar)%vals(1,:) = obs_field(:,1)
              write(*,*) 'MIN/MAX of ',trim(poolItr % memberName),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)

           else if (poolItr % nDims == 2) then
              call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r2d_ptr_a)
              ivar = ufo_vars_getindex(vars, trim(poolItr % memberName) )
              write(*,*) "interp: var, ufo_var_index = ",trim(poolItr % memberName), ivar
              if( .not. allocated(gom%geovals(ivar)%vals) )then
                 gom%geovals(ivar)%nval = fld%geom%nVertLevels
                 if(trim(poolItr % memberName).eq.var_prsi) gom%geovals(ivar)%nval = fld%geom%nVertLevelsP1 !BJJ: Can we do this better ??
                 allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
                 write(*,*) ' gom%geovals(n)%vals allocated'
              endif
              !write(*,*) 'MIN/MAX of ',trim(poolItr % memberName),minval(r2d_ptr_a),maxval(r2d_ptr_a)
              do jlev = 1, gom%geovals(ivar)%nval
                 mod_field(:,1) = r2d_ptr_a(jlev,1:ngrid)
                 call pbump%apply_obsop(mod_field,obs_field)
                 !ORG- gom%geovals(ivar)%vals(jlev,:) = obs_field(:,1)
                 gom%geovals(ivar)%vals(gom%geovals(ivar)%nval - jlev + 1,:) = obs_field(:,1) !BJJ-tmp vertical flip, top-to-bottom for CRTM geoval
              end do
              !write(*,*) 'MIN/MAX of ',trim(poolItr % memberName),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)

           else if (poolItr % nDims == 3) then
           end if
        end if

        end if
   end do !- end of pool iteration


   !---add special cases: var_sfc_wspeed and/or var_sfc_wdir
   if ( (ufo_vars_getindex(vars,var_sfc_wspeed)    .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_wdir) .ne. -1) ) then

     write(*,*) ' BJJ: special cases: var_sfc_wspeed and/or var_sfc_wdir'

     !- allocate
     allocate(tmp_field(nobs,2))

     !- read/interp.
     call mpas_pool_get_array(fld % auxFields, "u10", r1d_ptr_a)
     mod_field(:,1) = r1d_ptr_a(1:ngrid)
     write(*,*) 'MIN/MAX of u10=',minval(mod_field(:,1)),maxval(mod_field(:,1))
     call pbump%apply_obsop(mod_field,obs_field)
     tmp_field(:,1)=obs_field(:,1)
     call mpas_pool_get_array(fld % auxFields, "v10", r1d_ptr_a)
     mod_field(:,1) = r1d_ptr_a(1:ngrid)
     write(*,*) 'MIN/MAX of v10=',minval(mod_field(:,1)),maxval(mod_field(:,1))
     call pbump%apply_obsop(mod_field,obs_field)
     tmp_field(:,2)=obs_field(:,1)

     !- allocate geoval & put values for var_sfc_wspeed
     ivar = ufo_vars_getindex(vars, var_sfc_wspeed)
     if(ivar .ne. -1) then
       gom%geovals(ivar)%nval = 1
       if(allocated(gom%geovals(ivar)%vals)) deallocate(gom%geovals(ivar)%vals) !BJJ: for 4D H
       allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
       gom%geovals(ivar)%vals(1,:) = sqrt( tmp_field(:,1)**2 + tmp_field(:,2)**2 ) ! ws = sqrt(u**2+v**2) [m/s]
       write(*,*) 'MIN/MAX of ',trim(var_sfc_wspeed),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif

     !- allocate geoval & put values for var_sfc_wdir
     ivar = ufo_vars_getindex(vars, var_sfc_wdir)
     if(ivar .ne. -1) then
       gom%geovals(ivar)%nval = 1
       if(allocated(gom%geovals(ivar)%vals)) deallocate(gom%geovals(ivar)%vals) !BJJ: for 4D H
       allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
       do ii=1,nobs
         call uv_to_wdir(tmp_field(ii,1), tmp_field(ii,2), wdir) ! uu, vv, wind10_direction in radian
         gom%geovals(ivar)%vals(1,ii) = wdir / deg2rad           ! radian -> degree
       enddo
       write(*,*) 'MIN/MAX of ',trim(var_sfc_wdir),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif

     !- deallocate
     deallocate(tmp_field)
    
   endif  !---end special cases




   !---add special cases: var_sfc_landtyp, var_sfc_vegtyp, var_sfc_soiltyp
   if ( (ufo_vars_getindex(vars,var_sfc_landtyp)      .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_vegtyp)  .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_soiltyp) .ne. -1) ) then

     write(*,*) ' BJJ: special cases: var_sfc_landtyp, var_sfc_vegtyp, or var_sfc_soiltyp'

     call mpas_pool_get_array(fld % auxFields, "ivgtyp", i1d_ptr_a)
     write(*,*) 'MIN/MAX of ivgtyp=',minval(i1d_ptr_a),maxval(i1d_ptr_a)

     call mpas_pool_get_array(fld % auxFields, "isltyp", i1d_ptr_b)
     write(*,*) 'MIN/MAX of isltyp=',minval(i1d_ptr_b),maxval(i1d_ptr_b)


     !initialize vector of nearest neighbor indices
     allocate( index_nn(nobs) )
     allocate( weight_nn(pbump%obsop%h%n_s) )

     do ii=1,nobs
       !Picks index of pbump%obsop%h%S containing maxium weight for obs ii
       !Generic method for any interpolation scheme
       weight_nn = 0.0_kind_real
       where ( pbump%obsop%h%row .eq. ii ) 
          weight_nn = pbump%obsop%h%S
       end where
       jj = maxloc(weight_nn,1)

!       !Cheaper method that works for BUMP unstructured "triangular mesh" ( 3 vertices per obs ) with Bilinear interp.
!       jj=3*(ii-1) + maxloc(pbump%obsop%h%S( 3*(ii-1)+1:3*(ii-1)+3 ),1) !nearest-interp. / maximum-weight specified.

       !Store index of BUMP extended vector
       index_nn(ii) = pbump%obsop%h%col(jj)
     enddo

     deallocate(weight_nn)

     !- allocate geoval & put values for var_sfc_landtyp
     ivar = ufo_vars_getindex(vars, var_sfc_landtyp)
     if(ivar .ne. -1) then
       gom%geovals(ivar)%nval = 1
       if(allocated(gom%geovals(ivar)%vals)) deallocate(gom%geovals(ivar)%vals) !BJJ: for 4D H
       allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
       mod_field(:,1) = real( i1d_ptr_a(1:ngrid), kind_real)
       allocate( mod_field_ext(pbump%obsop%nc0b,1) )
       call pbump%obsop%com%ext(pbump%mpl,1,mod_field,mod_field_ext)
       do ii=1,nobs
         gom%geovals(ivar)%vals(1,ii) = mod_field_ext( index_nn(ii), 1 )
       enddo
       deallocate( mod_field_ext )
       write(*,*) 'MIN/MAX of ',trim(var_sfc_landtyp),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif

     !- allocate geoval & put values for var_sfc_vegtyp
     ivar = ufo_vars_getindex(vars, var_sfc_vegtyp)
     if(ivar .ne. -1) then
       gom%geovals(ivar)%nval = 1
       if(allocated(gom%geovals(ivar)%vals)) deallocate(gom%geovals(ivar)%vals) !BJJ: for 4D H
       allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
       mod_field(:,1) = real( i1d_ptr_a(1:ngrid), kind_real)
       allocate( mod_field_ext(pbump%obsop%nc0b,1) )
       call pbump%obsop%com%ext(pbump%mpl,1,mod_field,mod_field_ext)
       do ii=1,nobs
         gom%geovals(ivar)%vals(1,ii) = real( convert_type_veg( int(mod_field_ext( index_nn(ii), 1 )) ) , kind_real)
       enddo
       deallocate( mod_field_ext )
       write(*,*) 'MIN/MAX of ',trim(var_sfc_vegtyp),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif

     !- allocate geoval & put values for var_sfc_soiltyp
     ivar = ufo_vars_getindex(vars, var_sfc_soiltyp)
     if(ivar .ne. -1) then
       gom%geovals(ivar)%nval = 1
       if(allocated(gom%geovals(ivar)%vals)) deallocate(gom%geovals(ivar)%vals) !BJJ: for 4D H
       allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
       mod_field(:,1) = real( i1d_ptr_b(1:ngrid), kind_real)
       allocate( mod_field_ext(pbump%obsop%nc0b,1) )
       call pbump%obsop%com%ext(pbump%mpl,1,mod_field,mod_field_ext)
       do ii=1,nobs
         gom%geovals(ivar)%vals(1,ii) = real( convert_type_soil( int(mod_field_ext( index_nn(ii), 1 )) ), kind_real)
       enddo
       deallocate( mod_field_ext )
       write(*,*) 'MIN/MAX of ',trim(var_sfc_soiltyp),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif

     deallocate(index_nn)

   endif  !---end special cases


   !---add special cases: var_sfc_wfrac, var_sfc_lfrac, var_sfc_ifrac, var_sfc_sfrac
   !---    simple interpolation now, but can be more complex: Consider FOV ??
   if ( (ufo_vars_getindex(vars,var_sfc_wfrac)      .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_lfrac) .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_ifrac) .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_sfrac) .ne. -1) ) then

     write(*,*) ' BJJ: special cases: var_sfc_wfrac, var_sfc_lfrac, var_sfc_ifrac, or var_sfc_sfrac'

     call mpas_pool_get_array(fld % auxFields, "landmask", i1d_ptr_a) !"land-ocean mask (1=land ; 0=ocean)"
     write(*,*) 'MIN/MAX of landmask=',minval(i1d_ptr_a),maxval(i1d_ptr_a)
     call mpas_pool_get_array(fld % auxFields, "xice", r1d_ptr_a)     !"fractional area coverage of sea-ice"
     write(*,*) 'MIN/MAX of xice=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
     call mpas_pool_get_array(fld % auxFields, "snowc", r1d_ptr_b)    !"flag for snow on ground (=0 no snow; =1,otherwise"
     write(*,*) 'MIN/MAX of snowc=',minval(r1d_ptr_b),maxval(r1d_ptr_b)
     write(*,*) 'MIN/MAX of lnad+xice+snowc=', &
                minval(real(i1d_ptr_a)+r1d_ptr_a+r1d_ptr_b),maxval(real(i1d_ptr_a)+r1d_ptr_a+r1d_ptr_b)

     ivarw = ufo_vars_getindex(vars, var_sfc_wfrac)
     ivarl = ufo_vars_getindex(vars, var_sfc_lfrac)
     ivari = ufo_vars_getindex(vars, var_sfc_ifrac)
     ivars = ufo_vars_getindex(vars, var_sfc_sfrac)

     !--- Land first. will be adjusted later
     ivar = ivarl
     if(ivar .ne. -1) then
       gom%geovals(ivar)%nval = 1
       if(allocated(gom%geovals(ivar)%vals)) deallocate(gom%geovals(ivar)%vals) !BJJ: for 4D H
       allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
       mod_field(:,1) = real(i1d_ptr_a(1:ngrid))
       call pbump%apply_obsop(mod_field,obs_field)
       gom%geovals(ivar)%vals(1,:) = obs_field(:,1)
       write(*,*) 'MIN/MAX of ',trim(var_sfc_lfrac),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif
     !--- determine ICE
     ivar = ivari
     if(ivar .ne. -1) then
       gom%geovals(ivar)%nval = 1
       if(allocated(gom%geovals(ivar)%vals)) deallocate(gom%geovals(ivar)%vals) !BJJ: for 4D H
       allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
       mod_field(:,1) = r1d_ptr_a(1:ngrid)
       call pbump%apply_obsop(mod_field,obs_field)
       gom%geovals(ivar)%vals(1,:) = obs_field(:,1)
       write(*,*) 'MIN/MAX of ',trim(var_sfc_ifrac),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif
     !--- detemine/adjust SNOW & SEA
     ivar = ivars
     if(ivar .ne. -1) then
       gom%geovals(ivar)%nval = 1
       gom%geovals(ivarw)%nval = 1
       if(allocated(gom%geovals(ivar)%vals)) deallocate(gom%geovals(ivar)%vals) !BJJ: for 4D H
       allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,nobs) )
       if(allocated(gom%geovals(ivarw)%vals)) deallocate(gom%geovals(ivarw)%vals) !BJJ: for 4D H
       allocate( gom%geovals(ivarw)%vals(gom%geovals(ivarw)%nval,nobs) )
       mod_field(:,1) = r1d_ptr_b(1:ngrid)
       call pbump%apply_obsop(mod_field,obs_field)
       gom%geovals(ivar)%vals(1,:) = obs_field(:,1)
       do ii=1, nobs
         if(gom%geovals(ivari)%vals(1,ii).gt.0.0_kind_real) then
           gom%geovals(ivar)%vals(1,ii) = min( gom%geovals(ivar)%vals(1,ii), 1.0_kind_real - gom%geovals(ivari)%vals(1,ii) )
           gom%geovals(ivarw)%vals(1,ii)= 1.0_kind_real - gom%geovals(ivari)%vals(1,ii) - gom%geovals(ivar)%vals(1,ii)
         else
           gom%geovals(ivarw)%vals(1,ii)= 1.0_kind_real - gom%geovals(ivarl)%vals(1,ii)
         endif
       enddo
       write(*,*) 'MIN/MAX of ',trim(var_sfc_sfrac),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
       write(*,*) 'MIN/MAX of ',trim(var_sfc_wfrac),minval(gom%geovals(ivarw)%vals),maxval(gom%geovals(ivarw)%vals)
     endif
     !--- Final adjust LAND
     ivar = ivarl
     if(ivar .ne. -1) then
       do ii=1, nobs
         gom%geovals(ivar)%vals(1,ii) = max( 1.0_kind_real - gom%geovals(ivarw)%vals(1,ii) - gom%geovals(ivari)%vals(1,ii) &
                                                          - gom%geovals(ivars)%vals(1,ii), 0.0_kind_real)
       enddo
       write(*,*) 'MIN/MAX of ',trim(var_sfc_lfrac),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif
     !do ii=17,19  !1,nobs
     !write(*,*) gom%geovals(ivarl)%vals(1,ii), gom%geovals(ivarw)%vals(1,ii), &
     !           gom%geovals(ivari)%vals(1,ii), gom%geovals(ivars)%vals(1,ii)
     !enddo
   endif  !---end special cases

   !--- OMG: adjust between sfc coverage & sfc type
   !         see wrfda/da_get_innov_vector_crtm.inc#L521
   if ( (ufo_vars_getindex(vars,var_sfc_wfrac)      .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_lfrac) .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_ifrac) .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_sfrac) .ne. -1) ) then
   if ( (ufo_vars_getindex(vars,var_sfc_landtyp)      .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_vegtyp)  .ne. -1) &
        .or. (ufo_vars_getindex(vars,var_sfc_soiltyp) .ne. -1) ) then
     do ii=1, nobs
       if(gom%geovals(ivarl)%vals(1,ii) .gt. 0.0_kind_real) then
         if(nint(gom%geovals(ufo_vars_getindex(vars,var_sfc_soiltyp))%vals(1,ii)) .eq. 9 .or. &
            nint(gom%geovals(ufo_vars_getindex(vars,var_sfc_vegtyp))%vals(1,ii)) .eq. 13 ) then
           gom%geovals(ivari)%vals(1,ii) = min( gom%geovals(ivari)%vals(1,ii) + gom%geovals(ivarl)%vals(1,ii), 1.0_kind_real )
           gom%geovals(ivarl)%vals(1,ii) = 0.0_kind_real
         endif
       endif
     enddo
   endif
   endif  !--- OMG: end


   nullify(pbump)

   deallocate(mod_field)
   deallocate(obs_field)

   call mpas_pool_empty_pool(pool_ufo)
   call mpas_pool_destroy_pool(pool_ufo)

   write(*,*) '---- Leaving getvalues ---'
end subroutine getvalues

! ------------------------------------------------------------------------------

subroutine getvalues_tl(fld, locs, vars, gom, traj)

   implicit none
   type(mpas_field),      intent(inout) :: fld
   type(ioda_locs),       intent(in)    :: locs
   type(ufo_vars),        intent(in)    :: vars
   type(ufo_geovals),     intent(inout) :: gom
   type(mpas_getvaltraj), intent(in)    :: traj

   character(len=*), parameter :: myname = 'getvalues_tl'
   
   integer :: ii, jj, ji, jvar, jlev, ngrid, nobs, ivar
   real(kind=kind_real), allocatable :: mod_field(:,:)
   real(kind=kind_real), allocatable :: obs_field(:,:)
   
   type (mpas_pool_type), pointer :: pool_ufo
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   
   ! Check traj is implemented
   ! -------------------------
   if (.not.traj%lalloc) &
   call abor1_ftn(trim(myname)//" trajectory for this obs op not found")

   ! Get grid dimensions and checks
   ! ------------------------------
   ngrid = fld%geom%nCellsSolve !or traj%ngrid
   nobs = locs%nlocs            !or traj%nobs
   write(*,*)'getvalues_tl: ngrid, nobs = : ',ngrid, nobs
   call interp_checks("tl", fld, locs, vars, gom)
   
   !Make sure the return values are allocated and set
   !-------------------------------------------------
   do jvar=1,vars%nv
      if( .not. allocated(gom%geovals(jvar)%vals) )then
         gom%geovals(jvar)%nval = fld%geom%nVertLevels
         allocate( gom%geovals(jvar)%vals(fld%geom%nVertLevels,nobs) )
         write(*,*) ' gom%geovals(n)%vals allocated'
      endif
   enddo
   gom%linit = .true.
   

   !Create Buffer for interpolated values
   !--------------------------------------
   allocate(mod_field(ngrid,1))
   allocate(obs_field(nobs,1))
   
   !Interpolate fields to obs locations using pre-calculated weights
   !----------------------------------------------------------------
   write(0,*)'getvalues_tl: vars%nv       : ',vars%nv
   write(0,*)'getvalues_tl: vars%fldnames : ',vars%fldnames
   
   !------- need some table matching UFO_Vars & related MPAS_Vars
   !------- for example, Tv @ UFO may require Theta, Pressure, Qv.
   !-------                               or  Theta_m, exner_base, Pressure_base, Scalar(&index_qv)
   call convert_mpas_field2ufoTL(traj % pool_traj, fld % subFields, fld % auxFields, pool_ufo, vars % fldnames, vars % nv) !--pool_ufo is new pool with ufo_vars

   call mpas_pool_begin_iteration(pool_ufo)
   do while ( mpas_pool_get_next_member(pool_ufo, poolItr) )
        write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)
        if (poolItr % memberType == MPAS_POOL_FIELD) then
        if (poolItr % dataType == MPAS_POOL_REAL) then
           if (poolItr % nDims == 1) then

           else if (poolItr % nDims == 2) then
              call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r2d_ptr_a)
              write(*,*) "getvalues_tl: var=",trim(poolItr % memberName)
              ivar = ufo_vars_getindex(vars, trim(poolItr % memberName) )
              write(*,*) "ufo_vars_getindex: ivar=",ivar
              do jlev = 1, gom%geovals(ivar)%nval
                 mod_field(:,1) = r2d_ptr_a(jlev,1:ngrid)
                 call traj%bump%apply_obsop(mod_field,obs_field)
                 !ORG- gom%geovals(ivar)%vals(jlev,:) = obs_field(:,1)
                 gom%geovals(ivar)%vals(gom%geovals(ivar)%nval - jlev + 1,:) = obs_field(:,1) !BJJ-tmp vertical flip, top-to-bottom for CRTM geoval
              end do

           else if (poolItr % nDims == 3) then

           end if
        end if
        end if
   end do


   deallocate(mod_field)
   deallocate(obs_field)

   call mpas_pool_empty_pool(pool_ufo)
   call mpas_pool_destroy_pool(pool_ufo)

   write(*,*) '---- Leaving getvalues_tl ---'
end subroutine getvalues_tl

! ------------------------------------------------------------------------------

subroutine getvalues_ad(fld, locs, vars, gom, traj)

   implicit none
   type(mpas_field),      intent(inout) :: fld
   type(ioda_locs),       intent(in)    :: locs
   type(ufo_vars),        intent(in)    :: vars
   type(ufo_geovals),     intent(inout) :: gom
   type(mpas_getvaltraj), intent(in)    :: traj

   character(len=*), parameter :: myname = 'getvalues_ad'

   integer :: ii, jj, ji, jvar, jlev, ngrid, nobs, ivar
   real(kind=kind_real), allocatable :: mod_field(:,:)
   real(kind=kind_real), allocatable :: obs_field(:,:)

   type (mpas_pool_type), pointer :: pool_ufo
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b

   ! Check traj is implemented
   ! -------------------------
   if (.not.traj%lalloc) &
   call abor1_ftn(trim(myname)//" trajectory for this obs op not found")

   ! Get grid dimensions and checks
   ! ------------------------------
   ngrid = fld%geom%nCellsSolve !or traj%ngrid
   nobs = locs%nlocs            !or traj%nobs

   call interp_checks("ad", fld, locs, vars, gom)

   !Create Buffer for interpolated values
   !--------------------------------------
   allocate(mod_field(ngrid,1))
   allocate(obs_field(nobs,1))

   !Interpolate fields to obs locations using pre-calculated weights
   !----------------------------------------------------------------
   write(0,*)'getvalues_ad: vars%nv       : ',vars%nv
   write(0,*)'getvalues_ad: vars%fldnames : ',vars%fldnames

   !NOTE: This TL routine is called JUST to create "pool_ufo". Their values from TL routine doesn't matter.
   !    : Actually their values are initialized as "zero" in following "do while" loop.
   call convert_mpas_field2ufoTL(traj % pool_traj, fld % subFields, fld % auxFields, pool_ufo, vars % fldnames, vars % nv) !--pool_ufo is new pool with ufo_vars

   call mpas_pool_begin_iteration(pool_ufo)
   do while ( mpas_pool_get_next_member(pool_ufo, poolItr) )
        write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)
        if (poolItr % memberType == MPAS_POOL_FIELD) then
        if (poolItr % dataType == MPAS_POOL_REAL) then
           if (poolItr % nDims == 1) then

           else if (poolItr % nDims == 2) then
              call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r2d_ptr_a)
              r2d_ptr_a=0.0
              write(*,*) "Interp. var=",trim(poolItr % memberName)
              ivar = ufo_vars_getindex(vars, trim(poolItr % memberName) )
              write(*,*) "ufo_vars_getindex, ivar=",ivar
              do jlev = 1, gom%geovals(ivar)%nval
                 !ORG- obs_field(:,1) = gom%geovals(ivar)%vals(jlev,:)
                 obs_field(:,1) = gom%geovals(ivar)%vals(gom%geovals(ivar)%nval - jlev + 1,:) !BJJ-tmp vertical flip, top-to-bottom for CRTM geoval
                 call traj%bump%apply_obsop_ad(obs_field,mod_field)
                 r2d_ptr_a(jlev,1:ngrid) = 0.0_kind_real
                 r2d_ptr_a(jlev,1:ngrid) = r2d_ptr_a(jlev,1:ngrid) + mod_field(:,1)
              end do

           else if (poolItr % nDims == 3) then

           end if
        end if
        end if
   end do

   call convert_mpas_field2ufoAD(traj % pool_traj, fld % subFields, fld % auxFields, pool_ufo, vars % fldnames, vars % nv) !--pool_ufo is new pool with ufo_vars

   deallocate(mod_field)
   deallocate(obs_field)

   call mpas_pool_empty_pool(pool_ufo)
   call mpas_pool_destroy_pool(pool_ufo)

   write(*,*) '---- Leaving getvalues_ad ---' 
end subroutine getvalues_ad

! ------------------------------------------------------------------------------

subroutine initialize_interp(grid, locs, bump)

   use mpas_geom_mod, only: mpas_geom
   use type_bump, only: bump_type
   
   implicit none
   type(mpas_geom),          intent(in)  :: grid
   type(ioda_locs),          intent(in)  :: locs
   type(bump_type), pointer, intent(out) :: bump
   
   logical, save :: interp_initialized = .FALSE.
   
   integer :: mod_nz,mod_num
   real(kind=kind_real), allocatable :: mod_lat(:), mod_lon(:) 
   
   real(kind=kind_real), allocatable :: area(:),vunit(:,:)
   logical, allocatable :: lmask(:,:)

   integer, save :: bumpcount = 0
   character(len=5) :: cbumpcount
   character(len=17) :: bump_nam_prefix
   
   integer :: ii, jj, ji, jvar, jlev
 

   ! Each bump%nam%prefix must be distinct
   ! -------------------------------------
   bumpcount = bumpcount + 1
   write(cbumpcount,"(I0.2)") bumpcount
   bump_nam_prefix = 'mpas_bump_data_'//cbumpcount

   !Get the Solution dimensions
   !---------------------------
   mod_nz  = grid%nVertLevels
   mod_num = grid%nCellsSolve
   write(*,*)'initialize_interp mod_num,obs_num = ', mod_num, locs%nlocs
   
   !Calculate interpolation weight using BUMP
   !------------------------------------------
    allocate( mod_lat(mod_num), mod_lon(mod_num) )
    mod_lat(:) = grid%latCell( 1:mod_num ) / deg2rad !- to Degrees
    mod_lon(:) = grid%lonCell( 1:mod_num ) / deg2rad !- to Degrees

    !Important namelist options
    call bump%nam%init

    bump%nam%prefix       = bump_nam_prefix  ! Prefix for files output
    bump%nam%nobs         = locs%nlocs       ! Number of observations
    bump%nam%obsop_interp = 'bilin'          ! Interpolation type (bilinear)
    bump%nam%obsdis       = 'local'          ! Observation distribution parameter ('random', 'local', or 'adjusted')
    bump%nam%diag_interp  = 'bilin'
    bump%nam%local_diag   = .false.

    !Less important namelist options (should not be changed)
    bump%nam%default_seed        = .true.
    bump%nam%new_hdiag           = .false.
    bump%nam%new_nicas           = .false.
    bump%nam%check_adjoints      = .false.
    bump%nam%check_pos_def       = .false.
    bump%nam%check_sqrt          = .false.
    bump%nam%check_dirac         = .false.
    bump%nam%check_randomization = .false.
    bump%nam%check_consistency   = .false.
    bump%nam%check_optimality    = .false.
    bump%nam%new_lct             = .false.
    bump%nam%new_obsop           = .true.

    !Initialize geometry
    allocate(area(mod_num))
    allocate(vunit(mod_num,1))
    allocate(lmask(mod_num,1))
    area  = 1.0          ! Dummy area, unit [m^2]
    vunit = 1.0          ! Dummy vertical unit
    lmask = .true.       ! Mask

    !Initialize BUMP
    write(*,*) 'call bump%setup_online, nobs=locs%nlocs=',locs%nlocs
    write(*,*) 'call bump%setup_online, lonobs=locs%lon(1:3)=',locs%lon(1:3)
    write(*,*) 'call bump%setup_online, lonobs=locs%lat(1:3)=',locs%lat(1:3)
    call bump%setup_online(mpi_comm_world,mod_num,1,1,1,mod_lon,mod_lat,area,vunit,lmask, &
                           nobs=locs%nlocs,lonobs=locs%lon(:),latobs=locs%lat(:))

    !Release memory
    deallocate(area)
    deallocate(vunit)
    deallocate(lmask)
    deallocate( mod_lat, mod_lon )

end subroutine initialize_interp

! ------------------------------------------------------------------------------

subroutine interp_checks(cop, fld, locs, vars, gom)

   implicit none
   character(len=2),  intent(in) :: cop
   type(mpas_field),  intent(in) :: fld
   type(ioda_locs),   intent(in) :: locs
   type(ufo_vars),    intent(in) :: vars
   type(ufo_geovals), intent(in) :: gom
   
   integer :: jvar
   character(len=26) :: cinfo
   
   cinfo="mpas_fields:checks "//cop//" : "
   
   !Check things are the sizes we expect
   !------------------------------------
   if (gom%nobs /= locs%nlocs ) then
      call abor1_ftn(cinfo//"geovals wrong size")
   endif
   if( gom%nvar .ne. vars%nv )then
      call abor1_ftn(cinfo//"nvar wrong size")
   endif
   if( .not. allocated(gom%geovals) )then
      call abor1_ftn(cinfo//"geovals unallocated")
   endif
   if( size(gom%geovals) .ne. vars%nv )then
      call abor1_ftn(cinfo//"geovals wrong size")
   endif
   if (cop/="tl" .and. cop/='nl') then !BJJ or cop='ad'.   why only check ad??? because ad is set/allocated in UFO side ??
                                                           ! also, the dimension is not general for all possible geovals.
   if (.not.gom%linit) then
      call abor1_ftn(cinfo//"geovals not initialized")
   endif
   do jvar=1,vars%nv
      if (allocated(gom%geovals(jvar)%vals)) then  
         if( gom%geovals(jvar)%nval .ne. fld%geom%nVertLevels )then
            write(*,*) jvar, gom%geovals(jvar)%nval, fld%geom%nVertLevels
            call abor1_ftn(cinfo//"nval wrong size")
         endif
         if( gom%geovals(jvar)%nobs .ne. locs%nlocs )then
            call abor1_ftn(cinfo//"nobs wrong size")
         endif
         if( size(gom%geovals(jvar)%vals, 1) .ne. fld%geom%nVertLevels )then
            call abor1_ftn(cinfo//"vals wrong size 1")
         endif
         if( size(gom%geovals(jvar)%vals, 2) .ne. locs%nlocs )then
            call abor1_ftn(cinfo//"vals wrong size 2")
         endif       
      else
        call abor1_ftn(cinfo//"vals not allocated")
      endif 
   enddo
   endif
   
   write(*,*)'interp_checks ',cinfo,' done'
   
end subroutine interp_checks

! ------------------------------------------------------------------------------

subroutine ug_size(self, ug)

   use unstructured_grid_mod
   
   implicit none
   type(mpas_field),        intent(in)    :: self
   type(unstructured_grid), intent(inout) :: ug
   integer :: igrid
   
   ! Set number of grids
   if (ug%colocated==1) then
      ! Colocatd
      ug%ngrid = 1
   else
      ! Not colocatedd
      ug%ngrid = 1
   end if

   ! Allocate grid instances
   if (.not.allocated(ug%grid)) allocate(ug%grid(ug%ngrid))

   if (ug%colocated==1) then ! colocated

      ! Set local number of points
      ug%grid(1)%nmga = self%geom%nCellsSolve

      ! Set number of levels
      ug%grid(1)%nl0 = self%geom%nVertLevels

      ! Set number of variables
      ug%grid(1)%nv = self%nf

      ! Set number of timeslots
      ug%grid(1)%nts = 1

   else ! Not colocated

      do igrid=1,ug%ngrid
         ! Set local number of points
         ug%grid(igrid)%nmga = self%geom%nCellsSolve

         ! Set number of levels
         ug%grid(igrid)%nl0 = self%geom%nVertLevels

         ! Set number of variables
         ug%grid(igrid)%nv = self%nf

         ! Set number of timeslots
         ug%grid(igrid)%nts = 1
      enddo
   end if
end subroutine ug_size

! ------------------------------------------------------------------------------

subroutine ug_coord(self, ug, colocated)

   use unstructured_grid_mod
   
   implicit none
   type(mpas_field),        intent(in)    :: self
   type(unstructured_grid), intent(inout) :: ug
   integer,                 intent(in)    :: colocated
   
   integer :: jl, igrid
   
   ! Define size
   call ug_size(self, ug)

   ! Alocate unstructured grid coordinates
   call allocate_unstructured_grid_coord(ug)

   ! Copy coordinates
   if (ug%colocated==1) then ! colocated
     ug%grid(1)%lon = self%geom%lonCell(1:ug%grid(1)%nmga) / deg2rad !- to Degrees
     ug%grid(1)%lat = self%geom%latCell(1:ug%grid(1)%nmga) / deg2rad !- to Degrees
     ug%grid(1)%area = self%geom%areaCell(1:ug%grid(1)%nmga)
     do jl=1,self%geom%nVertLevels
       ug%grid(1)%vunit(:,jl) = real(jl,kind=kind_real)
       ug%grid(1)%lmask(:,jl) = .true.
     enddo
   else ! Not colocated
     do igrid=1,ug%ngrid
       ug%grid(igrid)%lon = self%geom%lonCell(1:ug%grid(igrid)%nmga) / deg2rad !- to Degrees
       ug%grid(igrid)%lat = self%geom%latCell(1:ug%grid(igrid)%nmga) / deg2rad !- to Degrees
       ug%grid(igrid)%area = self%geom%areaCell(1:ug%grid(igrid)%nmga)
       do jl=1,self%geom%nVertLevels
         ug%grid(igrid)%vunit(:,jl) = real(jl,kind=kind_real)
         ug%grid(igrid)%lmask(:,jl) = .true.
       enddo
     enddo
   endif

end subroutine ug_coord

! ------------------------------------------------------------------------------

subroutine field_to_ug(self, ug, colocated)

   use mpas_pool_routines
   use unstructured_grid_mod
   
   implicit none
   type(mpas_field),        intent(in)    :: self
   type(unstructured_grid), intent(inout) :: ug
   integer,                 intent(in)    :: colocated
   
   integer :: idx_var,jC,jl  
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   type(ufo_vars) :: vars ! temporary to access variable "index" easily

   ! Set list of variables
   vars % nv = self % nf
   allocate(vars % fldnames(vars % nv))
   vars % fldnames(:) = self % fldnames(:)

   ! Define size
   call ug_size(self, ug)

   ! Allocate unstructured grid field
   call allocate_unstructured_grid_field(ug)

   ! Copy field
   call mpas_pool_begin_iteration(self % subFields)
   
   do while ( mpas_pool_get_next_member(self % subFields, poolItr) )
        write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)
        ! Pools may in general contain dimensions, namelist options, fields, or other pools,
        ! so we select only those members of the pool that are fields
        if (poolItr % memberType == MPAS_POOL_FIELD) then
        ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
        if (poolItr % dataType == MPAS_POOL_REAL) then
           ! Depending on the dimensionality of the field, we need to set pointers of
           ! the correct type
           if (poolItr % nDims == 1) then
              write(*,*)'Not implemented yet'
           else if (poolItr % nDims == 2) then
              call mpas_pool_get_array(self % subFields, trim(poolItr % memberName), r2d_ptr_a)
              idx_var = -999
              idx_var = ufo_vars_getindex(vars, trim(poolItr % memberName))
              if(idx_var.gt.0) then
                 write(*,*) '  sub. field_to_ug, poolItr % memberName=',trim(poolItr % memberName)
                 write(*,*) '  sub. field_to_ug, idx_var=',idx_var
                 do jC=1,self%geom%nCellsSolve
                   do jl=1,self%geom%nVertLevels
                     ug%grid(1)%fld(jC,jl,idx_var,1) = r2d_ptr_a(jl,jC)
                   enddo
                 enddo
              endif
           else if (poolItr % nDims == 3) then
              write(*,*)'Not implemented yet'
           end if
        end if
        end if
   end do

   ! Cleanup
   call ufo_vars_delete(vars)

end subroutine field_to_ug

! -----------------------------------------------------------------------------

subroutine field_from_ug(self, ug)

   use mpas_pool_routines
   use unstructured_grid_mod

   implicit none
   type(mpas_field),        intent(inout) :: self
   type(unstructured_grid), intent(in)    :: ug
   
   integer :: idx_var,jC,jl
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   type(ufo_vars) :: vars ! temporary to access variable "index" easily

   ! Set list of variables
   vars % nv = self % nf
   allocate(vars % fldnames(vars % nv))
   vars % fldnames(:) = self % fldnames(:)

   ! Copy field
   call mpas_pool_begin_iteration(self % subFields)
   
   do while ( mpas_pool_get_next_member(self % subFields, poolItr) )
        write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)
        ! Pools may in general contain dimensions, namelist options, fields, or other pools,
        ! so we select only those members of the pool that are fields
        if (poolItr % memberType == MPAS_POOL_FIELD) then
        ! Fields can be integer, logical, or real. Here, we operate only on real-valued fields
        if (poolItr % dataType == MPAS_POOL_REAL) then
           ! Depending on the dimensionality of the field, we need to set pointers of
           ! the correct type
           if (poolItr % nDims == 1) then
              write(*,*)'Not implemented yet'
           else if (poolItr % nDims == 2) then
              call mpas_pool_get_array(self % subFields, trim(poolItr % memberName), r2d_ptr_a)

              idx_var = -999
              idx_var = ufo_vars_getindex(vars, trim(poolItr % memberName))
              if(idx_var.gt.0) then
                 write(*,*) '  sub. convert_from_ug, poolItr % memberName=',trim(poolItr % memberName)
                 do jC=1,self%geom%nCellsSolve
                   do jl=1,self%geom%nVertLevels
                     r2d_ptr_a(jl,jC) = ug%grid(1)%fld(jC,jl,idx_var,1)
                   enddo
                 enddo
              end if
           else if (poolItr % nDims == 3) then
              write(*,*)'Not implemented yet'
           end if
        end if
        end if
   end do

   ! Cleanup
   call ufo_vars_delete(vars)

   ! TODO: Since only local locations are updated/transferred from ug, 
   !       need MPAS HALO comms before using these fields in MPAS

end subroutine field_from_ug

! ------------------------------------------------------------------------------

end module mpas_fields_mod
