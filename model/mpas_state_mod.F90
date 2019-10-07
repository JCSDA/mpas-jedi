! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module mpas_state_mod

use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm

!oops
use datetime_mod
use kinds, only: kind_real
use variables_mod

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

!mpas-jedi
use mpas_constants_mod
use mpas_geom_mod
use mpas_getvaltraj_mod, only: mpas_getvaltraj
use mpas_field_utils_mod
use mpas_increment_utils_mod, only: mpas_increment
use mpas_state_utils_mod
use mpas2ufo_vars_mod
use mpas4da_mod

implicit none

private

public :: add_incr, &
        & analytic_IC, invent_state, &
        & getvalues

! ------------------------------------------------------------------------------

contains


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

   use mpas2ufo_vars_mod, only: q_to_w, temp_to_theta, twp_to_rho

   implicit none
   class(mpas_state),     intent(inout) :: self !< state
   class(mpas_increment), intent(in)    :: rhs  !< increment
   character(len=StrKIND) :: kind_op

   type (mpas_pool_type), pointer :: state, diag, mesh
   type (field2DReal), pointer :: field2d_t, field2d_p, field2d_sh, field2d_uRz, field2d_uRm, &
                                  field2d_th, field2d_qv, field2d_rho, field2d_u, field2d_u_inc

   ! GD: I don''t see any difference than for self_add other than subFields can contain
   ! different variables than mpas_state and the resolution of incr can be different. 

   if (self%geom%nCells==rhs%geom%nCells .and. self%geom%nVertLevels==rhs%geom%nVertLevels) then
      !NOTE: first, get full state of "subFields" variables
      kind_op = 'add'
      call da_operator(trim(kind_op), self % subFields, rhs % subFields, fld_select = rhs % fldnames_ci)

      ! Impose positive-definite limits on hydrometeors
      call da_posdef( self % subFields, mpas_hydrometeor_fields)

      !NOTE: second, also update variables which are closely related to MPAS prognostic vars.
      !  update index_qv (water vapor mixing ratio) from spechum (specific humidity) [ w = q / (1 - q) ]
      !  update theta from temperature and pressure
      !  update rho   from temperature, pressure, and index_qv
      call mpas_pool_get_field(self % subFields,            'temperature', field2d_t)
      call mpas_pool_get_field(self % subFields,               'pressure', field2d_p)
      call mpas_pool_get_field(self % subFields,                'spechum', field2d_sh)
      call mpas_pool_get_field(self % subFields,      'uReconstructZonal', field2d_uRz)
      call mpas_pool_get_field(self % subFields, 'uReconstructMeridional', field2d_uRm)
      call mpas_pool_get_field(self % subFields,                  'theta', field2d_th)
      call mpas_pool_get_field(self % subFields,               'index_qv', field2d_qv)
      call mpas_pool_get_field(self % subFields,                    'rho', field2d_rho)

      ! Ensure positive sh
      call da_posdef( self % subFields, (/'spechum'/))

      call temp_to_theta( field2d_t % array(:,:), field2d_p % array(:,:), field2d_th % array(:,:))
!      write(*,*) 'add_inc: theta min/max = ', minval(field2d_th % array), maxval(field2d_th % array)
      call q_to_w( field2d_sh % array(:,:), field2d_qv % array(:,:) )
!      write(*,*) 'add_inc: index_qv min/max = ', minval(field2d_sh % array), maxval(field2d_sh % array)
      ! Ensure positive qv : BJJ Do we need this? just in case ? or positive sh would be enough ?
      call da_posdef( self % subFields, (/'index_qv'/))

      call twp_to_rho( field2d_t % array(:,:), field2d_qv % array(:,:), field2d_p % array(:,:), &
                       field2d_rho % array(:,:) )
!      write(*,*) 'add_inc: rho min/max = ', minval(field2d_rho % array), maxval(field2d_rho % array)

      !  update u     from uReconstructZonal and uReconstructMeridional "incrementally"
      call mpas_pool_get_field(self % subFields,                      'u', field2d_u)
      call mpas_pool_get_field( rhs % subFields,      'uReconstructZonal', field2d_uRz)
      call mpas_pool_get_field( rhs % subFields, 'uReconstructMeridional', field2d_uRm)

      call mpas_duplicate_field(field2d_u, field2d_u_inc)

!      write(*,*) 'add_inc: u_inc min/max = ', minval(field2d_uRz % array), maxval(field2d_uRz % array)
!      write(*,*) 'add_inc: v_inc min/max = ', minval(field2d_uRm % array), maxval(field2d_uRm % array)

      call uv_cell_to_edges(self % geom % domain, field2d_uRz, field2d_uRm, field2d_u_inc, &
                 self%geom%latCell, self%geom%lonCell, self%geom%nCellsSolve, &
                 self%geom%edgeNormalVectors, self%geom%nEdgesOnCell, self%geom%edgesOnCell, &
                 self%geom%nVertLevels)
!      write(*,*) 'add_inc: u_guess min/max = ', minval(field2d_u % array), maxval(field2d_u % array)
!      write(*,*) 'add_inc: u_inc min/max = ', minval(field2d_u_inc % array), maxval(field2d_u_inc % array)
      field2d_u % array(:,:) = field2d_u % array(:,:) + field2d_u_inc % array(:,:)
!      write(*,*) 'add_inc: u_analy min/max = ', minval(field2d_u % array), maxval(field2d_u % array)

      ! TODO: DO we need HALO exchange here or in ModelMPAS::initialize for model integration?

      call mpas_deallocate_field( field2d_u_inc )
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

  use dcmip_initial_conditions_test_1_2_3, only : test1_advection_deformation, &
       test1_advection_hadley, test3_gravity_wave
  use dcmip_initial_conditions_test_4, only : test4_baroclinic_wave

!  !MPAS Test Cases
!  !JJG: This initialization requires the init_atmospher_core core_type 
!  !      in the MPAS library for OOPS, but currently it is not included
!  use init_atm_core, only: init_atm_core_run!, init_atm_core_finalize (could be used for cleanup...)

  implicit none

  class(mpas_state),         intent(inout) :: self   !< State
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
           ps_ptr(ii) = ps0
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
           ps_ptr(ii) = ps0
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
           ps_ptr(ii) = ps0
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

   class(mpas_state),         intent(inout) :: self    !< Model fields
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
         r2d_ptr_a(jlev,ii) = 1.0_kind_real
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
      r1d_ptr_a(ii) = 1.0_kind_real
   enddo
   !Scalars (state pool)
   !qv
   call mpas_pool_get_field(pool_a, 'scalars', field3d)
   call mpas_pool_get_dimension(state, 'index_qv', index_qv)
   if ( index_qv .gt. 0 ) then
      r2d_ptr_a => field3d % array(index_qv,:,:)
      do jlev = 1,self % geom % nVertLevels
         do ii = 1, self % geom % nCellsSolve
            r2d_ptr_a(jlev,ii) = 0.0_kind_real
         enddo
      enddo
   end if

   return

end subroutine invent_state

! ------------------------------------------------------------------------------

subroutine getvalues(self, locs, vars, gom, traj)

   use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum
   use type_bump, only: bump_type
   use mpas2ufo_vars_mod !, only: usgs_to_crtm_mw, wrf_to_crtm_soil

   implicit none
   class(mpas_state),                       intent(in)    :: self
   type(ufo_locs),                          intent(in)    :: locs
   type(oops_vars),                         intent(in)    :: vars
   type(ufo_geovals),                       intent(inout) :: gom
   type(mpas_getvaltraj), optional, target, intent(inout) :: traj
   
   character(len=*), parameter :: myname = 'getvalues'

   type(fckit_mpi_comm)     :: f_comm
   type(bump_type), target  :: bump
   type(bump_type), pointer :: pbump => null()
   logical,         target  :: bump_alloc = .false.
   logical,         pointer :: pbump_alloc => null()
   integer, save, target    :: bumpid = 1000
   integer ,        pointer :: pbumpid => null()
   
   integer :: jj, jvar, jlev, ilev, jloc, ngrid, maxlevels, nlevels, ivar, nlocs, nlocsg
   real(kind=kind_real), allocatable :: mod_field(:,:), mod_field_ext(:,:)
   real(kind=kind_real), allocatable :: obs_field(:,:)
   real(kind=kind_real), allocatable :: tmp_field(:,:)  !< for wspeed/wdir
   
   type (mpas_pool_type), pointer :: pool_ufo  !< pool with ufo variables
   type (mpas_pool_iterator_type) :: poolItr
   real (kind=kind_real), pointer :: r0d_ptr_a, r0d_ptr_b
   real (kind=kind_real), dimension(:), pointer :: r1d_ptr_a, r1d_ptr_b
   real (kind=kind_real), dimension(:,:), pointer :: r2d_ptr_a, r2d_ptr_b
   real (kind=kind_real), dimension(:,:,:), pointer :: r3d_ptr_a, r3d_ptr_b
   integer, dimension(:), pointer :: i1d_ptr_a, i1d_ptr_b!, i2d_ptr_a
   integer, allocatable  :: index_nn(:)
   real (kind=kind_real), allocatable :: weight_nn(:)
   type (mpas_pool_type), pointer :: pool_tmp  !< temporary pool for setting trajectory
   type (field2DReal), pointer :: field2d_src => null() !< for setting trajectory

   real(kind=kind_real) :: wdir           !< for wind direction
   integer :: ivarw, ivarl, ivari, ivars  !< for sfc fraction indices
   character(len=MAXVARLEN) :: ufo_var_name
   character(len=1024)      :: buf

   ! Get grid dimensions and checks
   ! ------------------------------
   ngrid = self%geom%nCellsSolve
   nlocs = locs%nlocs ! # of location for given time slot, given processor

   !If no observations can early exit
   !---------------------------------
   f_comm = fckit_mpi_comm()
   call f_comm%allreduce(nlocs,nlocsg,fckit_mpi_sum())
   if (nlocsg == 0) then
     if (present(traj)) then
        traj%lalloc = .true.
        traj%noobs = .true.
     endif
     return
   endif

   call interp_checks("nl", self, locs, vars, gom)

   ! Initialize the interpolation trajectory
   ! ---------------------------------------
   if (present(traj)) then

     pbump => traj % bump

     if (.not. traj%lalloc) then

       traj%ngrid = ngrid

       call mpas_pool_create_pool( pool_tmp )

       call mpas_pool_get_field(self % subFields, 'temperature', field2d_src)
       call mpas_pool_add_field(pool_tmp, 'temperature', field2d_src)
       call mpas_pool_get_field(self % subFields, 'spechum', field2d_src)
       call mpas_pool_add_field(pool_tmp, 'spechum', field2d_src)
       call mpas_pool_get_field(self % subFields, 'pressure', field2d_src)
       call mpas_pool_add_field(pool_tmp, 'pressure', field2d_src)

       call mpas_pool_clone_pool(pool_tmp, traj % pool_traj)

       call mpas_pool_empty_pool(pool_tmp)
       call mpas_pool_destroy_pool(pool_tmp)

       pbump_alloc => traj%lalloc
       pbumpid => traj%bumpid

    endif

  else

    pbump => bump
    bump_alloc = .false.
    pbump_alloc => bump_alloc
    bumpid = bumpid + 1
    pbumpid => bumpid

  endif

  if (.not. pbump_alloc) then 
    ! Calculate interpolation weight using BUMP
    ! ------------------------------------------
!    write(*,*)'call initialize_bump(...)'
    call initialize_bump(self%geom, locs, pbump, pbumpid)
    pbump_alloc = .true.
!    write(*,*)'interp: after initialize_bump'
  endif

   !Make sure the return values are allocated and set below
   !-------------------------------------------------------
   gom%linit = .true.
      
   !Interpolate fields to obs locations using pre-calculated weights
   !----------------------------------------------------------------
   write(0,*)'getvalues   : vars%nv       : ',vars%nv
   write(0,*)'getvalues   : vars%fldnames : ',vars%fldnames
   write(0,*)'getvalues   : nlocs, nlocsg : ',nlocs, ',',  nlocsg

   !------- need some table matching UFO_Vars & related MPAS_Vars
   !------- for example, Tv @ UFO may require Theta, Pressure, Qv.
   !-------                               or  Theta_m, exner_base, Pressure_base, Scalar(&index_qv)
   call convert_mpas_field2ufo(self % geom, self % subFields, pool_ufo, vars % fldnames, vars % nv, ngrid) !--pool_ufo is new pool with ufo_vars

   maxlevels = self%geom%nVertLevelsP1
   allocate(mod_field(ngrid,maxlevels))
   allocate(obs_field(nlocs,maxlevels))

   call mpas_pool_begin_iteration(pool_ufo)
   do while ( mpas_pool_get_next_member(pool_ufo, poolItr) )
      if (poolItr % memberType == MPAS_POOL_FIELD) then
         ufo_var_name = trim(poolItr % memberName)
         ivar = str_match(ufo_var_name, vars%fldnames)
         if ( ivar == -1 ) cycle

!         write(*,*) 'poolItr % nDims , poolItr % memberName =', poolItr % nDims , trim(poolItr % memberName)
         if (poolItr % nDims == 1) then
            if( .not. allocated(gom%geovals(ivar)%vals) )then
               gom%geovals(ivar)%nval = 1
               allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
!               write(*,*) ' gom%geovals(n)%vals allocated'
            endif
            nlevels = gom%geovals(ivar)%nval

            if (poolItr % dataType == MPAS_POOL_INTEGER) then
               call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), i1d_ptr_a)
!               write(*,*) "interp: var, ufo_var_index = ",trim(poolItr % memberName), ivar
!               write(*,*) 'MIN/MAX of ',trim(poolItr % memberName),minval(i1d_ptr_a),maxval(i1d_ptr_a)
               mod_field(:,1) = real( i1d_ptr_a(1:ngrid), kind_real)
            else if (poolItr % dataType == MPAS_POOL_REAL) then
               call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r1d_ptr_a)
!               write(*,*) "interp: var, ufo_var_index = ",trim(poolItr % memberName), ivar
!               write(*,*) 'MIN/MAX of ',trim(poolItr % memberName),minval(r1d_ptr_a),maxval(r1d_ptr_a)
               mod_field(:,1) = r1d_ptr_a(1:ngrid)
            end if

         else if (poolItr % nDims == 2) then
            if( .not. allocated(gom%geovals(ivar)%vals) )then
               gom%geovals(ivar)%nval = self%geom%nVertLevels
               if(trim(poolItr % memberName).eq.var_prsi) gom%geovals(ivar)%nval = self%geom%nVertLevelsP1 !BJJ: Can we do this better ??
               allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
!               write(*,*) ' gom%geovals(n)%vals allocated'
            endif
            nlevels = gom%geovals(ivar)%nval

            if (poolItr % dataType == MPAS_POOL_INTEGER) then
!               call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), i1d_ptr_a)
!!               write(*,*) "interp: var, ufo_var_index = ",trim(poolItr % memberName), ivar
!!               write(*,*) 'MIN/MAX of ',trim(poolItr % memberName),minval(i2d_ptr_a),maxval(i2d_ptr_a)
!               mod_field(:,1:nlevels) = real( transpose(r2d_ptr_a(:,1:ngrid)), kind_real )

            else if (poolItr % dataType == MPAS_POOL_REAL) then
               call mpas_pool_get_array(pool_ufo, trim(poolItr % memberName), r2d_ptr_a)
!               write(*,*) "interp: var, ufo_var_index = ",trim(poolItr % memberName), ivar
!               write(*,*) 'MIN/MAX of ',trim(poolItr % memberName),minval(r2d_ptr_a),maxval(r2d_ptr_a)
               mod_field(:,1:nlevels) = transpose(r2d_ptr_a(:,1:ngrid))
            end if

!         else if (poolItr % nDims == 3) then

         end if

         !TODO - JJG: Reduce wall-time of getvalues/_tl/_ad
         ! + apply_obsop takes ~50% of wall-time of getvalues on cheyenne and
         !   scales with node count. Seems to have MPI-related issue.
         !
         ! + initialize_bump takes other ~50% on cheyennne
         !
         pbump%geom%nl0 = nlevels
         call pbump%apply_obsop(mod_field(:,1:nlevels),obs_field(:,1:nlevels))

         do jlev = 1, nlevels
            !BJJ-tmp vertical flip, top-to-bottom for CRTM geoval
            ! only selected obs (using locs%indx()) are filling "geovals"
            ilev = nlevels - jlev + 1
            gom%geovals(ivar)%vals(ilev, locs%indx(1:nlocs)) = obs_field(:,jlev)
         end do
!         write(*,*) 'MIN/MAX of ',trim(poolItr % memberName),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)

      end if
   end do !- end of pool iteration
   deallocate(mod_field)
   deallocate(obs_field)

   pbump%geom%nl0 = 1
   allocate(mod_field(ngrid,1))
   allocate(obs_field(nlocs,1))

   !---add special cases: var_sfc_wspeed and/or var_sfc_wdir
   if ( (str_match(var_sfc_wspeed,vars%fldnames)    .ne. -1) &
        .or. (str_match(var_sfc_wdir,vars%fldnames) .ne. -1) ) then

!     write(*,*) ' BJJ: special cases: var_sfc_wspeed and/or var_sfc_wdir'

     !- allocate
     allocate(tmp_field(nlocs,2))

     !- read/interp.
     call mpas_pool_get_array(self % subFields, "u10", r1d_ptr_a)
     mod_field(:,1) = r1d_ptr_a(1:ngrid)
!     write(*,*) 'MIN/MAX of u10=',minval(mod_field(:,1)),maxval(mod_field(:,1))
     call pbump%apply_obsop(mod_field,obs_field)
     tmp_field(:,1)=obs_field(:,1)
     call mpas_pool_get_array(self % subFields, "v10", r1d_ptr_a)
     mod_field(:,1) = r1d_ptr_a(1:ngrid)
!     write(*,*) 'MIN/MAX of v10=',minval(mod_field(:,1)),maxval(mod_field(:,1))
     call pbump%apply_obsop(mod_field,obs_field)
     tmp_field(:,2)=obs_field(:,1)

     !- allocate geoval & put values for var_sfc_wspeed
     ivar = str_match(var_sfc_wspeed,vars%fldnames)
     if(ivar .ne. -1) then
       if( .not. allocated(gom%geovals(ivar)%vals) )then
          gom%geovals(ivar)%nval = 1
          allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
       end if
       do jloc = 1, nlocs
         gom%geovals(ivar)%vals(1,locs%indx(jloc)) = sqrt( tmp_field(jloc,1)**2 + tmp_field(jloc,2)**2 ) ! ws = sqrt(u**2+v**2) [m/s]
       end do
!       write(*,*) 'MIN/MAX of ',trim(var_sfc_wspeed),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif

     !- allocate geoval & put values for var_sfc_wdir
     ivar = str_match(var_sfc_wdir,vars%fldnames)
     if(ivar .ne. -1) then
       if( .not. allocated(gom%geovals(ivar)%vals) )then
          gom%geovals(ivar)%nval = 1
          allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
       end if
       do jloc = 1, nlocs
         call uv_to_wdir(tmp_field(jloc,1), tmp_field(jloc,2), wdir) ! uu, vv, wind10_direction in radian
         gom%geovals(ivar)%vals(1,locs%indx(jloc)) = wdir / deg2rad           ! radian -> degree
       enddo
!       write(*,*) 'MIN/MAX of ',trim(var_sfc_wdir),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif

     !- deallocate
     deallocate(tmp_field)
    
   endif  !---end special cases




   !---add special cases: var_sfc_landtyp, var_sfc_vegtyp, var_sfc_soiltyp
   if ( (str_match(var_sfc_landtyp,vars%fldnames)      .ne. -1) &
        .or. (str_match(var_sfc_vegtyp,vars%fldnames)  .ne. -1) &
        .or. (str_match(var_sfc_soiltyp,vars%fldnames) .ne. -1) ) then

!     write(*,*) ' BJJ: special cases: var_sfc_landtyp, var_sfc_vegtyp, or var_sfc_soiltyp'

     call mpas_pool_get_array(self % subFields, "ivgtyp", i1d_ptr_a)
!     write(*,*) 'MIN/MAX of ivgtyp=',minval(i1d_ptr_a),maxval(i1d_ptr_a)

     call mpas_pool_get_array(self % subFields, "isltyp", i1d_ptr_b)
!     write(*,*) 'MIN/MAX of isltyp=',minval(i1d_ptr_b),maxval(i1d_ptr_b)


     !initialize vector of nearest neighbor indices
     allocate( index_nn(nlocs) )
     allocate( weight_nn(pbump%obsop%h%n_s) )

     do jloc = 1, nlocs
       !Picks index of pbump%obsop%h%S containing maxium weight for obs jloc
       !Generic method for any interpolation scheme
       weight_nn = 0.0_kind_real
       where ( pbump%obsop%h%row .eq. jloc ) 
          weight_nn = pbump%obsop%h%S
       end where
       jj = maxloc(weight_nn,1)

!       !Cheaper method that works for BUMP unstructured "triangular mesh" ( 3 vertices per obs ) with Bilinear interp.
!       jj=3*(jloc-1) + maxloc(pbump%obsop%h%S( 3*(jloc-1)+1:3*(jloc-1)+3 ),1) !nearest-interp. / maximum-weight specified.

       !Store index of BUMP extended vector
       index_nn(jloc) = pbump%obsop%h%col(jj)
     enddo

     deallocate(weight_nn)

     !- allocate geoval & put values for var_sfc_landtyp
     ivar = str_match(var_sfc_landtyp,vars%fldnames)
     if(ivar .ne. -1) then
       if( .not. allocated(gom%geovals(ivar)%vals) )then
          gom%geovals(ivar)%nval = 1
          allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
       end if
       mod_field(:,1) = real( i1d_ptr_a(1:ngrid), kind_real)
       allocate( mod_field_ext(pbump%obsop%nc0b,1) )
       call pbump%obsop%com%ext(pbump%mpl,1,mod_field,mod_field_ext)
       do jloc = 1, nlocs
         gom%geovals(ivar)%vals(1,locs%indx(jloc)) = mod_field_ext( index_nn(jloc), 1 )
       enddo
       deallocate( mod_field_ext )
!       write(*,*) 'MIN/MAX of ',trim(var_sfc_landtyp),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif

     !- allocate geoval & put values for var_sfc_vegtyp
     ivar = str_match(var_sfc_vegtyp,vars%fldnames)
     if(ivar .ne. -1) then
       if( .not. allocated(gom%geovals(ivar)%vals) )then
          gom%geovals(ivar)%nval = 1
          allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
       end if
       mod_field(:,1) = real( i1d_ptr_a(1:ngrid), kind_real)
       allocate( mod_field_ext(pbump%obsop%nc0b,1) )
       call pbump%obsop%com%ext(pbump%mpl,1,mod_field,mod_field_ext)
       do jloc = 1, nlocs
         gom%geovals(ivar)%vals(1,locs%indx(jloc)) = real( convert_type_veg( int(mod_field_ext( index_nn(jloc), 1 )) ) , kind_real)
       enddo
       deallocate( mod_field_ext )
!       write(*,*) 'MIN/MAX of ',trim(var_sfc_vegtyp),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif

     !- allocate geoval & put values for var_sfc_soiltyp
     ivar = str_match(var_sfc_soiltyp,vars%fldnames)
     if(ivar .ne. -1) then
       if( .not. allocated(gom%geovals(ivar)%vals) )then
          gom%geovals(ivar)%nval = 1
          allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
       end if
       mod_field(:,1) = real( i1d_ptr_b(1:ngrid), kind_real)
       allocate( mod_field_ext(pbump%obsop%nc0b,1) )
       call pbump%obsop%com%ext(pbump%mpl,1,mod_field,mod_field_ext)
       do jloc = 1, nlocs
         gom%geovals(ivar)%vals(1,locs%indx(jloc)) = real( convert_type_soil( int(mod_field_ext( index_nn(jloc), 1 )) ), kind_real)
       enddo
       deallocate( mod_field_ext )
!       write(*,*) 'MIN/MAX of ',trim(var_sfc_soiltyp),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif

     deallocate(index_nn)

   endif  !---end special cases


   !---add special cases: var_sfc_wfrac, var_sfc_lfrac, var_sfc_ifrac, var_sfc_sfrac
   !---    simple interpolation now, but can be more complex: Consider FOV ??
   if ( (str_match(var_sfc_wfrac,vars%fldnames)      .ne. -1) &
        .or. (str_match(var_sfc_lfrac,vars%fldnames) .ne. -1) &
        .or. (str_match(var_sfc_ifrac,vars%fldnames) .ne. -1) &
        .or. (str_match(var_sfc_sfrac,vars%fldnames) .ne. -1) ) then

!     write(*,*) ' BJJ: special cases: var_sfc_wfrac, var_sfc_lfrac, var_sfc_ifrac, or var_sfc_sfrac'

     call mpas_pool_get_array(self % subFields, "landmask", i1d_ptr_a) !"land-ocean mask (1=land ; 0=ocean)"
!     write(*,*) 'MIN/MAX of landmask=',minval(i1d_ptr_a),maxval(i1d_ptr_a)
     call mpas_pool_get_array(self % subFields, "xice", r1d_ptr_a)     !"fractional area coverage of sea-ice"
!     write(*,*) 'MIN/MAX of xice=',minval(r1d_ptr_a),maxval(r1d_ptr_a)
     call mpas_pool_get_array(self % subFields, "snowc", r1d_ptr_b)    !"flag for snow on ground (=0 no snow; =1,otherwise"
!     write(*,*) 'MIN/MAX of snowc=',minval(r1d_ptr_b),maxval(r1d_ptr_b)
!     write(*,*) 'MIN/MAX of lnad+xice+snowc=', &
!                minval(real(i1d_ptr_a)+r1d_ptr_a+r1d_ptr_b),maxval(real(i1d_ptr_a)+r1d_ptr_a+r1d_ptr_b)

     ivarw = str_match(var_sfc_wfrac,vars%fldnames)
     ivarl = str_match(var_sfc_lfrac,vars%fldnames)
     ivari = str_match(var_sfc_ifrac,vars%fldnames)
     ivars = str_match(var_sfc_sfrac,vars%fldnames)

     !--- Land first. will be adjusted later
     ivar = ivarl
     if(ivar .ne. -1) then
       if( .not. allocated(gom%geovals(ivar)%vals) )then
          gom%geovals(ivar)%nval = 1
          allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
       end if
       mod_field(:,1) = real(i1d_ptr_a(1:ngrid))
       call pbump%apply_obsop(mod_field,obs_field)
       do jloc = 1, nlocs
          gom%geovals(ivar)%vals(1,locs%indx(jloc)) = obs_field(jloc,1)
       enddo
!       write(*,*) 'MIN/MAX of ',trim(var_sfc_lfrac),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif
     !--- determine ICE
     ivar = ivari
     if(ivar .ne. -1) then
       if( .not. allocated(gom%geovals(ivar)%vals) )then
          gom%geovals(ivar)%nval = 1
          allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
       end if
       mod_field(:,1) = r1d_ptr_a(1:ngrid)
       call pbump%apply_obsop(mod_field,obs_field)
       do jloc = 1, nlocs
          gom%geovals(ivar)%vals(1,locs%indx(jloc)) = obs_field(jloc,1)
       enddo
!       write(*,*) 'MIN/MAX of ',trim(var_sfc_ifrac),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif
     !--- detemine/adjust SNOW & SEA
     ivar = ivars
     if(ivar .ne. -1) then
       if( .not. allocated(gom%geovals(ivar)%vals) )then
          gom%geovals(ivar)%nval = 1
          allocate( gom%geovals(ivar)%vals(gom%geovals(ivar)%nval,gom%geovals(ivar)%nlocs) )
       end if
       if( .not. allocated(gom%geovals(ivarw)%vals) )then
          gom%geovals(ivarw)%nval = 1
          allocate( gom%geovals(ivarw)%vals(gom%geovals(ivarw)%nval,gom%geovals(ivarw)%nlocs) )
       end if
       mod_field(:,1) = r1d_ptr_b(1:ngrid)
       call pbump%apply_obsop(mod_field,obs_field)
       do jloc = 1, nlocs
          gom%geovals(ivar)%vals(1,locs%indx(jloc)) = obs_field(jloc,1)
       enddo
       do jloc = 1, nlocs
         jj = locs%indx(jloc)
         if(gom%geovals(ivari)%vals(1,jj).gt.0.0_kind_real) then
           gom%geovals(ivar)%vals(1,jj) = min( gom%geovals(ivar)%vals(1,jj), 1.0_kind_real - gom%geovals(ivari)%vals(1,jj) )
           gom%geovals(ivarw)%vals(1,jj)= 1.0_kind_real - gom%geovals(ivari)%vals(1,jj) - gom%geovals(ivar)%vals(1,jj)
         else
           gom%geovals(ivarw)%vals(1,jj)= 1.0_kind_real - gom%geovals(ivarl)%vals(1,jj)
         endif
       enddo
!       write(*,*) 'MIN/MAX of ',trim(var_sfc_sfrac),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
!       write(*,*) 'MIN/MAX of ',trim(var_sfc_wfrac),minval(gom%geovals(ivarw)%vals),maxval(gom%geovals(ivarw)%vals)
     endif
     !--- Final adjust LAND
     ivar = ivarl
     if(ivar .ne. -1) then
       do jloc = 1, nlocs
         jj = locs%indx(jloc)
         gom%geovals(ivar)%vals(1,jj) = max( 1.0_kind_real - gom%geovals(ivarw)%vals(1,jj) - gom%geovals(ivari)%vals(1,jj) &
                                                          - gom%geovals(ivars)%vals(1,jj), 0.0_kind_real)
       enddo
!       write(*,*) 'MIN/MAX of ',trim(var_sfc_lfrac),minval(gom%geovals(ivar)%vals),maxval(gom%geovals(ivar)%vals)
     endif
     !do ii=17,19  !1,nlocs
     !write(*,*) gom%geovals(ivarl)%vals(1,ii), gom%geovals(ivarw)%vals(1,ii), &
     !           gom%geovals(ivari)%vals(1,ii), gom%geovals(ivars)%vals(1,ii)
     !enddo
   endif  !---end special cases

   !--- OMG: adjust between sfc coverage & sfc type
   !         see wrfda/da_get_innov_vector_crtm.inc#L521
   if ( (str_match(var_sfc_wfrac,vars%fldnames)      .ne. -1) &
        .or. (str_match(var_sfc_lfrac,vars%fldnames) .ne. -1) &
        .or. (str_match(var_sfc_ifrac,vars%fldnames) .ne. -1) &
        .or. (str_match(var_sfc_sfrac,vars%fldnames) .ne. -1) ) then
   if ( (str_match(var_sfc_landtyp,vars%fldnames)      .ne. -1) &
        .or. (str_match(var_sfc_vegtyp,vars%fldnames)  .ne. -1) &
        .or. (str_match(var_sfc_soiltyp,vars%fldnames) .ne. -1) ) then
     do jloc = 1, nlocs
       jj = locs%indx(jloc)
       if(gom%geovals(ivarl)%vals(1,jj) .gt. 0.0_kind_real) then
         if(nint(gom%geovals(str_match(var_sfc_soiltyp,vars%fldnames))%vals(1,jj)) .eq. 9 .or. &
            nint(gom%geovals(str_match(var_sfc_vegtyp,vars%fldnames))%vals(1,jj)) .eq. 13 ) then
           gom%geovals(ivari)%vals(1,jj) = min( gom%geovals(ivari)%vals(1,jj) + gom%geovals(ivarl)%vals(1,jj), 1.0_kind_real )
           gom%geovals(ivarl)%vals(1,jj) = 0.0_kind_real
         endif
       endif
     enddo
   endif
   endif  !--- OMG: end


   if (.not. present(traj)) then
      call pbump%dealloc()
   endif

   nullify(pbump)
   nullify(pbump_alloc)
   nullify(pbumpid)

   ! Deallocate local memory
   ! -----------------------
   deallocate(mod_field)
   deallocate(obs_field)

   call mpas_pool_destroy_pool(pool_ufo)

!   write(*,*) '---- Leaving getvalues ---'
end subroutine getvalues

! ------------------------------------------------------------------------------

subroutine initialize_bump(grid, locs, bump, bumpid)

   use mpas_geom_mod, only: mpas_geom
   use type_bump, only: bump_type
   
   implicit none
   type(mpas_geom),          intent(in)  :: grid
   type(ufo_locs),           intent(in)  :: locs
   type(bump_type), pointer, intent(out) :: bump
   integer,                  intent(in)  :: bumpid
   
   integer :: mod_nz,mod_num
   real(kind=kind_real), allocatable :: mod_lat(:), mod_lon(:) 
   real(kind=kind_real), allocatable :: area(:),vunit(:,:)
   logical, allocatable :: lmask(:,:)

   character(len=5)   :: cbumpcount
   character(len=255) :: bump_nam_prefix
   
   type(fckit_mpi_comm) :: f_comm

   ! Each bump%nam%prefix must be distinct
   ! -------------------------------------
   write(cbumpcount,"(I0.5)") bumpid
   bump_nam_prefix = 'mpas_bump_data_'//cbumpcount

   ! Get the Solution dimensions
   ! ---------------------------
   mod_num = grid%nCellsSolve
!   write(*,*)'initialize_bump mod_num,obs_num = ', mod_num, locs%nlocs
   
   !Calculate interpolation weight using BUMP
   !------------------------------------------
   allocate( mod_lat(mod_num) )
   allocate( mod_lon(mod_num) )
   mod_lat(:) = grid%latCell( 1:mod_num ) / deg2rad !- to Degrees
   mod_lon(:) = grid%lonCell( 1:mod_num ) / deg2rad !- to Degrees

   ! Namelist options
   ! ----------------

   !Important namelist options
   call bump%nam%init

   !Less important namelist options (should not be changed)
   bump%nam%prefix       = trim(bump_nam_prefix)  ! Prefix for files output
   bump%nam%default_seed = .true.
   bump%nam%new_obsop    = .true.
   bump%nam%write_obsop  = .false.
   bump%nam%verbosity    = "none"

   ! Initialize geometry
   ! -------------------
   allocate(area(mod_num))
   allocate(vunit(mod_num,1))
   allocate(lmask(mod_num,1))

   area  = 1.0          ! Dummy area, unit [m^2]
   vunit = 1.0          ! Dummy vertical unit
   lmask = .true.       ! Mask

   ! Initialize BUMP
   ! ---------------
   f_comm = fckit_mpi_comm()
   call bump%setup_online(f_comm,mod_num,1,1,1,mod_lon,mod_lat,area,vunit,lmask, &
                          nobs=locs%nlocs,lonobs=locs%lon(:),latobs=locs%lat(:))

   ! Run BUMP drivers
   call bump%run_drivers

   ! Partial deallocate option
   call bump%partial_dealloc

   ! Release memory
   ! --------------
   deallocate(area)
   deallocate(vunit)
   deallocate(lmask)
   deallocate(mod_lat)
   deallocate(mod_lon)

end subroutine initialize_bump

! ------------------------------------------------------------------------------

end module mpas_state_mod
