! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_state_mod

use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum

!oops
use datetime_mod
use kinds, only: kind_real
use oops_variables_mod, only: oops_variables

!dcmip initialization
use dcmip_initial_conditions_test_1_2_3, only : test1_advection_deformation, &
       test1_advection_hadley, test3_gravity_wave
use dcmip_initial_conditions_test_4, only : test4_baroclinic_wave

!ufo
use ufo_locs_mod
use ufo_geovals_mod
use ufo_vars_mod

!MPAS-Model
use mpas_constants
use mpas_derived_types
use mpas_field_routines
use mpas_kind_types, only: StrKIND
use mpas_pool_routines
use mpas_dmpar, only: mpas_dmpar_exch_halo_field


!mpas-jedi
use mpas_constants_mod
use mpas_geom_mod
use mpas_field_utils_mod
use mpas2ufo_vars_mod
use mpas4da_mod

implicit none

private

public :: add_incr, &
        & analytic_IC, invent_state

! ------------------------------------------------------------------------------

contains


! ------------------------------------------------------------------------------
!> add increment to state
!!
!! \details **add_incr()** adds "increment" to "state", such as
!!          state (containing analysis) = state (containing guess) + increment
!!          Here, we also update "theta", "rho", and "u" (edge-normal wind), which are
!!          close to MPAS prognostic variable. 
!!          Intermediate 3D pressure is diagnosed with hydrostatic balance.
!!          While conversion to "theta" and "rho" uses full state variables,
!!          conversion to "u" from cell center winds uses their increment to reduce 
!!          the smoothing effect.
!!
subroutine add_incr(self,rhs)

   implicit none
   class(mpas_field), intent(inout) :: self !< state
   class(mpas_field), intent(in)    :: rhs  !< increment
   character(len=StrKIND) :: kind_op

   integer :: ngrid
   type (mpas_pool_type), pointer :: state, diag, mesh
   type (field2DReal), pointer :: fld2d_t, fld2d_sh, fld2d_p, fld2d_pp, fld2d_pb, &
                                  fld2d_uRz, fld2d_uRm, fld2d_u, fld2d_u_inc &
                                  fld2d_th, fld2d_qv, fld2d_rho
   type (field1DReal), pointer :: fld1d_ps

   ! GD: I don''t see any difference than for self_add other than subFields can contain
   ! different variables than mpas_field and the resolution of incr can be different. 

   if (self%geom%nCells==rhs%geom%nCells .and. self%geom%nVertLevels==rhs%geom%nVertLevels) then
      !NOTE: first, get full state of "subFields" variables
      kind_op = 'add'
      call da_operator(trim(kind_op), self % subFields, rhs % subFields, fld_select = rhs % fldnames_ci)

      ! Impose positive-definite limits on hydrometeors
      call da_posdef( self % subFields, mpas_hydrometeor_fields)

      !NOTE: second, also update variables which are closely related to MPAS prognostic vars.
      call mpas_pool_get_field(self % subFields,      'temperature', fld2d_t)
      call mpas_pool_get_field(self % subFields,          'spechum', fld2d_sh)
      call mpas_pool_get_field(self % subFields, 'surface_pressure', fld1d_ps)
      call mpas_pool_get_field(self % subFields,         'index_qv', fld2d_qv)
      call mpas_pool_get_field(self % subFields,         'pressure', fld2d_p)
      call mpas_pool_get_field(self % subFields,              'rho', fld2d_rho)
      call mpas_pool_get_field(self % subFields,            'theta', fld2d_th)

      ! Ensure positive sh
      call da_posdef( self % subFields, (/'spechum'/))

      ngrid = self%geom%nCellsSolve
      ! Update index_qv (water vapor mixing ratio) from spechum (specific humidity) [ w = q / (1 - q) ]
      call q_to_w( fld2d_sh % array(:,1:ngrid), fld2d_qv % array(:,1:ngrid) )
!      write(*,*) 'add_inc: index_qv min/max = ', minval(fld2d_sh % array), maxval(fld2d_sh % array)

      ! Diagnose 3D pressure and update full state theta and rho
      call hydrostatic_balance( ngrid, self%geom%nVertLevels, self%geom%zgrid(:,1:ngrid), &
                fld2d_t%array(:,1:ngrid), fld2d_qv%array(:,1:ngrid), &
                fld1d_ps%array(1:ngrid), fld2d_p%array(:,1:ngrid), &
                fld2d_rho%array(:,1:ngrid), fld2d_th%array(:,1:ngrid) )

      ! Update edge normal wind u from uReconstructZonal and uReconstructMeridional "incrementally"
      call mpas_pool_get_field(self % subFields,                      'u', fld2d_u)
      call mpas_pool_get_field( rhs % subFields,      'uReconstructZonal', fld2d_uRz)
      call mpas_pool_get_field( rhs % subFields, 'uReconstructMeridional', fld2d_uRm)

      call mpas_duplicate_field(fld2d_u, fld2d_u_inc)

!      write(*,*) 'add_inc: u_inc min/max = ', minval(fld2d_uRz % array), maxval(fld2d_uRz % array)
!      write(*,*) 'add_inc: v_inc min/max = ', minval(fld2d_uRm % array), maxval(fld2d_uRm % array)

      call mpas_dmpar_exch_halo_field(fld2d_uRz)
      call mpas_dmpar_exch_halo_field(fld2d_uRm)
      call uv_cell_to_edges(self % geom % domain, fld2d_uRz, fld2d_uRm, fld2d_u_inc, &
                 self%geom%lonCell, self%geom%latCell, self%geom%nCells, &
                 self%geom%edgeNormalVectors, self%geom%nEdgesOnCell, self%geom%edgesOnCell, &
                 self%geom%nVertLevels)
!      write(*,*) 'add_inc: u_guess min/max = ', minval(fld2d_u % array), maxval(fld2d_u % array)
!      write(*,*) 'add_inc: u_inc min/max = ', minval(fld2d_u_inc % array), maxval(fld2d_u_inc % array)
      ngrid = self%geom%nEdgesSolve
      fld2d_u % array(:,1:ngrid) = fld2d_u % array(:,1:ngrid) + fld2d_u_inc % array(:,1:ngrid)
!      write(*,*) 'add_inc: u_analy min/max = ', minval(fld2d_u % array), maxval(fld2d_u % array)

      ! TODO: DO we need HALO exchange here or in ModelMPAS::initialize for model integration?

      call mpas_deallocate_field( fld2d_u_inc )

      ! Update pressure_p (pressure perturbation) , which is a diagnostic variable
      call mpas_pool_get_field(self % subFields, 'pressure_p', fld2d_pp)
      call mpas_pool_get_field(self % geom % domain % blocklist % allFields, 'pressure_base', fld2d_pb)
      fld2d_pp % array(:,1:ngrid) = fld2d_p % array(:,1:ngrid) - fld2d_pb % array(:,1:ngrid)
   else
      call abor1_ftn("mpas_state:add_incr: dimension mismatch")
   endif

   return

end subroutine add_incr

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
subroutine analytic_IC(self, geom, f_conf, vdate)

!  !MPAS Test Cases
!  !JJG: This initialization requires the init_atmospher_core core_type 
!  !      in the MPAS library for OOPS, but currently it is not included
!  use init_atm_core, only: init_atm_core_run!, init_atm_core_finalize (could be used for cleanup...)

  implicit none

  class(mpas_field),         intent(inout) :: self   !< State
  type(mpas_geom), target,   intent(in)    :: geom   !< Geometry 
  type(fckit_configuration), intent(in)    :: f_conf !< Configuration
  type(datetime),            intent(inout) :: vdate  !< DateTime

  character(len=:), allocatable :: str
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
  real(kind=kind_real), dimension(:),pointer :: ps_ptr
  integer, pointer :: index_qv, index_qc, index_qr, index_qi, index_qs
  type (field3DReal), pointer :: field3d
  type (mpas_pool_type), pointer :: pool_a, pool_b, state
  real(kind=kind_real) :: zhalf

  ! Pointer to geometry component of field object
  self%geom => geom

  If (f_conf%has("analytic_init")) Then
     call f_conf%get_or_die("analytic_init",str)
     IC = str
  Else
     ! This default value is for backward compatibility
     IC = "invent-state"
  EndIf

  WRITE(*,*) "mpas_state:analytic_init: "//IC

! Conflicts with natural log below
!  call log%warning("mpas_state:analytic_init: "//IC)
  call f_conf%get_or_die("date",str)
  sdate = str
  WRITE(buf,*) 'validity date is: '//sdate
  WRITE(*,*) buf
!  call log%info(buf)
  call datetime_set(sdate, vdate)

   ! Need to initialize variables that are used in interpolation/getVals
   ! In "create" and "read" subroutines, subFields are 
   ! initialized from geom % domain % blocklist % allFields, zeroed,
   ! reread from file into allFields, then values copied to subFields
   ! -> must initialize allFields here and copy to subFields

   call mpas_pool_get_subpool(geom % domain % blocklist % structs, &
                              'state', state)

   pool_a => geom % domain % blocklist % allFields

   !Diagnostic vars (diag pool)
   call mpas_pool_get_array(pool_a, "pressure", p_ptr)
   call mpas_pool_get_array(pool_a, "uReconstructZonal", u_ptr)
   call mpas_pool_get_array(pool_a, "uReconstructMeridional", v_ptr)
   call mpas_pool_get_array(pool_a, "temperature", temperature_ptr)
   call mpas_pool_get_array(pool_a, "surface_pressure", ps_ptr)


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

        call invent_state(self,f_conf)


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
!!Would use "config_init_case" in the yaml file, then check for matching with one of the ideal cases below... (not 7 or 8)
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
!          call abor1_ftn("mpas_state: init_atm_core_run failed")
!       end if

     Case ("dcmip-test-1-1")
        do ii = 1, geom%nCellsSolve
           rlat = geom%latCell(ii)
           rlon = geom%lonCell(ii)

           ! Now loop over all levels
           do jlev = 1, geom%nVertLevels

              zhalf = MPAS_JEDI_HALF_kr * (geom%zgrid(jlev,ii) + geom%zgrid(jlev+1,ii))
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
           ps_ptr(ii) = ps0
        enddo

     Case ("dcmip-test-1-2")

        do ii = 1, geom%nCellsSolve
           rlat = geom%latCell(ii)
           rlon = geom%lonCell(ii)

           ! Now loop over all levels
           do jlev = 1, geom%nVertLevels

              zhalf = MPAS_JEDI_HALF_kr * (geom%zgrid(jlev,ii) + geom%zgrid(jlev+1,ii))
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
           ps_ptr(ii) = ps0
        enddo

     Case ("dcmip-test-3-1")

        do ii = 1, geom%nCellsSolve
           rlat = geom%latCell(ii)
           rlon = geom%lonCell(ii)

           ! Now loop over all levels
           do jlev = 1, geom%nVertLevels

              zhalf = MPAS_JEDI_HALF_kr * (geom%zgrid(jlev,ii) + geom%zgrid(jlev+1,ii))
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
           ps_ptr(ii) = ps0
        enddo

     Case ("dcmip-test-4-0")

        do ii = 1, geom%nCellsSolve
           rlat = geom%latCell(ii)
           rlon = geom%lonCell(ii)

           ! Now loop over all levels
           do jlev = 1, geom%nVertLevels

              zhalf = MPAS_JEDI_HALF_kr * (geom%zgrid(jlev,ii) + geom%zgrid(jlev+1,ii))
              Call test4_baroclinic_wave(0,MPAS_JEDI_ONE_kr,rlon,rlat,pk,zhalf,1,u0,v0,w0,&
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
           ps_ptr(ii) = ps0
        enddo

     Case Default

        call invent_state(self,f_conf)

     End Select int_option

     call da_copy_all2sub_fields(self % geom % domain, self % subFields) 

!   write(*,*)'==> end mpas_state:analytic_init'

end subroutine analytic_IC

! ------------------------------------------------------------------------------

subroutine invent_state(self,f_conf)

   implicit none

   class(mpas_field),         intent(inout) :: self    !< Model fields
   type(fckit_configuration), intent(in)    :: f_conf  !< Configuration structure
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   integer :: jlev,ii
   type (mpas_pool_type), pointer :: pool_a, state
   type (field3DReal), pointer :: field3d
   integer, pointer :: index_qv

   !- read/interp.

   call mpas_pool_get_subpool(self % geom % domain % blocklist % structs, &
                              'state', state)
   pool_a => self % geom % domain % blocklist % allFields

   !Diagnostic vars (diag pool)
   !u
   call mpas_pool_get_array(pool_a, "uReconstructZonal", r2d_ptr_a)
   do jlev = 1,self % geom % nVertLevels
      do ii = 1, self % geom % nCellsSolve
         r2d_ptr_a(jlev,ii) = cos(0.25*self % geom % lonEdge(ii)) + cos(0.25*self % geom % latEdge(ii))
      enddo
   enddo

   !v
   call mpas_pool_get_array(pool_a, "uReconstructMeridional", r2d_ptr_a)
   do jlev = 1,self % geom % nVertLevels
      do ii = 1, self % geom % nCellsSolve
         r2d_ptr_a(jlev,ii) = MPAS_JEDI_ONE_kr
      enddo
   enddo

   !temperature
   call mpas_pool_get_array(pool_a, "temperature", r2d_ptr_a)
   do jlev = 1,self % geom % nVertLevels
      do ii = 1, self % geom % nCellsSolve
         r2d_ptr_a(jlev,ii) = cos(0.25*self % geom % lonCell(ii)) + cos(0.25*self % geom % latCell(ii))
      enddo
   enddo

   !pressure
   call mpas_pool_get_array(pool_a, "pressure", r2d_ptr_a)
   do jlev = 1,self % geom % nVertLevels
      do ii = 1, self % geom % nCellsSolve
         r2d_ptr_a(jlev,ii) = real(jlev,kind_real)
      enddo
   enddo

   !surface_pressure
   call mpas_pool_get_array(pool_a, "surface_pressure", r1d_ptr_a)
   do ii = 1, self % geom % nCellsSolve
      r1d_ptr_a(ii) = MPAS_JEDI_ONE_kr
   enddo
   !Scalars (state pool)
   !qv
   call mpas_pool_get_field(pool_a, 'scalars', field3d)
   call mpas_pool_get_dimension(state, 'index_qv', index_qv)
   if ( index_qv .gt. 0 ) then
      r2d_ptr_a => field3d % array(index_qv,:,:)
      do jlev = 1,self % geom % nVertLevels
         do ii = 1, self % geom % nCellsSolve
            r2d_ptr_a(jlev,ii) = MPAS_JEDI_ZERO_kr
         enddo
      enddo
   end if

   return

end subroutine invent_state

! ------------------------------------------------------------------------------

end module mpas_state_mod
